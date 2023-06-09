library(tidyverse)
library(fuzzyjoin)
library(pcaMethods)
library(biomaRt)
library(imputeLCMD)

# define pipeline
# pg_df: input dataframe
# removeContaminant, removeReverse: whether rows with possible contaminants or reverse sequences, respectively, are kept (if set to false) or removed (if set to true).
# type: type of intensity to extract (possible values: "MS2", "LFQ", "Raw")
# log2t: whether intensities are transformed using log2.
# oneProteinPerRow: if set to true, N proteins identified in the same row are split into N different rows.
runMQ <- function(pg_df, removeContaminant = TRUE, removeReverse = TRUE, type = "Raw", log2t = FALSE, oneProteinPerRow = FALSE, contaminantList = NULL) {
  print("Running Part 1 of 6...")
  # create copy of input dataframe
  print("Attempting to create copy of input tibble...")
  out_df <- tibble(pg_df)
  print("Successfully created copy of input tibble.")
  
  if (removeContaminant & is.null(contaminantList)) {
    stop("removeContaminant selected, contaminantList must contain an argument!")
  }
  print("*************************************")
  print("contaminantList")
  print(contaminantList)


  
  # create id column
  print("Attempting to detect if id column exists...")
  if ("id" %in% colnames(out_df)) {
    print("id column exists.")
  }
  else {
    print("id column does not exist, attempting to create id column...")
    out_df <- out_df %>%
      mutate(id = row_number())
  }
  
  # convert "." in columns to " "
  print("Attemtping to convert periods in column names to spaces...")
  names(out_df) <- sub("Protein.IDs", "Protein IDs", names(out_df))
  names(out_df) <- sub("Intensity.", "Intensity ", names(out_df))
  names(out_df) <- sub("LFQ.intensity.", "LFQ intensity ", names(out_df))
  names(out_df) <- sub("MS.MS.count.", "MS/MS count ", names(out_df))
  names(out_df) <- sub("Potential.contaminant", "Potential contaminant", names(out_df))
  names(out_df) <- sub("Gene.names", "Gene names", names(out_df))
  print("Successfully changed column names.")
  
  # remove unnecessary columns
  print("Attempting to remove unnecessary columns from input tibble...")
  print(paste("Number of columns:", ncol(out_df)))
  if ("MS/MS count" %in% colnames(out_df)) {
    out_df <- out_df %>%
      dplyr::select(-`MS/MS count`)
  }
  out_df <- out_df %>%
    dplyr::select(`Protein IDs`, starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"), Reverse, `Potential contaminant`, id, -Intensity, -`Gene names`)
  print(paste("Number of columns after removal:", ncol(out_df)))
  print("Successfully removed unnecessary columns from input tibble.")
  
  # perform log2 transformation, if specified
  # 0 is transformed to -Inf
  if (log2t) {
    print("Attempting to perform log2 transformation on input tibble...")
    out_df <- out_df %>%
      mutate(across(matches("Intensity") | matches("LFQ intensity") | matches("MS/MS count"), log2))
    print("Successfully performed log2 transformation on input tibble.")
    
    # change -Inf to NA
    print("Attempting to change -Inf values to NA values...")
    is.na(out_df) <- sapply(out_df, is.infinite)
    print("Successfully changed -Inf values to NA values.")
  }
  else {
    # change 0 to NA
    print("Attempting to change 0 values to NA values...")
    #is.na(out_df) <- sapply(out_df, is.zero) # MEW: error - is.zero is unknown
    out_df[out_df == 0] <- NA
    print("Successfully changed 0 values to NA values.")
  }
  
  # remove columns of other types of intensities
  if (type == "Raw") {
    print("Raw intensity selected, attempting to remove LFQ intensity and MS/MS count columns...")
    print(paste("Number of columns:", ncol(out_df)))
    out_df <- out_df %>%
      dplyr::select(-starts_with("LFQ intensity"), -starts_with("MS/MS count"))
    print(paste("Number of columns after removal:", ncol(out_df)))
    print("Successfully removed LFQ intensity and MS/MS count columns.")
  } else if (type == "LFQ") {
    print("LFQ intensity selected, attempting to remove raw intensity and MS/MS count columns...")
    print(paste("Number of columns:", ncol(out_df)))
    out_df <- out_df %>%
      dplyr::select(-starts_with("Intensity"), -starts_with("MS/MS count"))
    print(paste("Number of columns after removal:", ncol(out_df)))
    print("Successfully removed raw intensity and MS/MS count columns.")
  } else if (type == "MS2") {
    print("MS/MS count selected, attempting to remove raw intensity and LFQ intensity columns...")
    print(paste("Number of columns:", ncol(out_df)))
    out_df <- out_df %>%
      dplyr::select(-starts_with("Intensity"), -starts_with("LFQ intensity"))
    print(paste("Number of columns after removal:", ncol(out_df)))
    print("Successfully removed raw intensity and LFQ intensity columns.")
  } else {
    warning("type should be specified as 'Raw', 'LFQ', or 'MS2'.")
  }
  
  # remove rows with contaminants and reverses, if specified
  if (removeContaminant) {
    #out_df <- out_df %>%
    #  mutate(temp = c(unlist(str_split(`Protein IDs`, ";")) %in% contaminantList$`Entry (UniProtAC)`))
    contaminant_col = contaminantList$`Entry (UniProtAC)`
    print("contaminant_col")
    print(contaminant_col)
    print("head(out_df)")
    print(out_df)
    out_df_2 <- out_df %>%
      group_by(id) %>%
      separate_rows(`Protein IDs`, sep=";") %>%
      filter(!any(`Protein IDs` %in% contaminantList$`Entry (UniProtAC)`))
    print("head(out_df_2)")
    print(out_df_2)
    out_df <- out_df_2 %>%
      summarize(id=unique(id)) %>%
      left_join(out_df)
    print("head(out_df)")
    print(out_df)
  }
  
  if (removeReverse) {
    print("Remove reverses selected, attempting to rows that were potentially reversed...")
    print(paste("Number of rows:", nrow(out_df)))
    out_df <- out_df %>%
      filter(is.na(Reverse))
    print(paste("Number of rows after removal:", nrow(out_df)))
    print("Successfully removed rows that were potentially reversed.")
  }
  
  out_df <- out_df %>%
    dplyr::select(-`Potential contaminant`, -Reverse)
  
  # split each row into more rows depending on amount of identified proteins, if oneProteinPerRow = TRUE
  if (oneProteinPerRow) {
    print("One protein per row selected, attempting to separate each row corresponding to individual proteins...")
    print("---------------------------------------> Showing first 6 rows of out_df:")
    print(head(out_df))
    print(paste("Number of rows:", nrow(out_df)))
    print(head(dplyr::select(out_df, id, `Protein IDs`)))
    # NOTE: here we used to have each protein on a different row, and then later we call distinct(id)
    # on the data, which selected a random protein id for each row 
    # now here we make sure we are taking the first protein from each row
    # See: https://github.com/mskcc/chiosis-chaperomics-shiny/issues/17
    #out_df <- out_df %>%
    #  separate_rows(`Protein IDs`, sep = ';')
    out_df <- out_df %>% tidyr::separate(`Protein IDs`, into=c("Protein IDs", NA), sep="\\s?;\\s?", extra="drop")
    print(paste("Number of rows after separations:", nrow(out_df)))
    print(head(dplyr::select(out_df, id, `Protein IDs`)))
    print("---------------------------------------> Showing first 6 rows of out_df:")
    print(head(out_df))
    print("Successfully separated each row corresponding to individual proteins.")
  }
  
  # output a dataframe
  print("Part 1 of 6 complete, returning tibble...")
  return(out_df)
}

# pg_df: tibble containing protein groups
# mf_df: tibble representing sample manifest
# sample.id.col: name of the column containing sample IDs
# intensity_type: "Raw", "LFQ", or "MS2" (REMOVED)
combine_matrix <- function (pg_df, mf_df, sample.id.col='Sample.ID2') {
  print("Running part 2 of 6...")
  print("Attempting to create copy of input tibble...")
  proteingroups <- tibble(pg_df)
  print("Succesfully created copy of input tible.")
  
  print("Attempting to convert protein groups tibble to long format...")
  print(paste("Number of rows:", nrow(pg_df)))
  print(paste("Number of columns:", ncol(pg_df)))
  # convert protein groups tibble to long format
  proteingroups <- proteingroups %>%
    gather("Sample", "Intensity", starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"))
  # if (intensity_type == "Raw") {
  #   proteingroups <- proteingroups %>%
  #     gather("Sample", "Intensity", starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"))
  # } else if (intensity_type == "LFQ") {
  #   proteingroups <- proteingroups %>%
  #     gather("Sample", "LFQ intensity", starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"))
  # } else if (intensity_type == "MS2") {
  #   proteingroups <- proteingroups %>%
  #     gather("Sample", "MS/MS count", starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"))
  # } else { 
  #   warning("type should be specified as 'Raw', 'LFQ', or 'MS2'.")
  # }
  print(paste("Number of rows after conversion:", nrow(pg_df)))
  print(paste("Number of columns after conversion:", ncol(pg_df)))
  print("Successfully converted protein groups tibble to long format.")
  
  # add regex column to manifest
  print("Attempting to add regex column to sample manifest tibble...")
  manifest_local <- mf_df %>%
    mutate(Sample.ID.regex = paste("(Intensity|LFQ intensity|MS/MS count)\\s", !!sym(sample.id.col), "\\b", sep=""))
  
  gsub(".", "[.]", manifest_local$Sample.ID.regex)
  
  print("Successfully added regex column to the sample manifest tibble.")
  
  # use left join to add patient diagnosis to sample ID
  print("Attempting to join protein groups and sample manifest tibbles...")
  print(paste("Number of columns:", ncol(pg_df)))
  combined_df <- regex_left_join(proteingroups, manifest_local, by=c(Sample = "Sample.ID.regex"))[-ncol(proteingroups) - ncol(manifest_local)] %>%
    dplyr::select(-!!sym(sample.id.col)) %>%
    distinct()
  print(paste("Number of columns after join:", ncol(pg_df)))
  print("Successfully joined protein groups and sample manifest tibbles.")
  
  print("Part 2 of 6 complete, returning combined tibble...")
  return(combined_df)
}

# pg_df: tibble containing protein groups
# method: used to determine whether to filter NA values per group by percentage ("percentage") or raw amount ("raw")
# filterGroup: column used to create groups for filtering
# TODO 100% percentageNA doesn't work
# percentageNA: minimum percentage of non-NA entries needed within a filter group in order to keep the protein group
# maxNA: maximum amount of NA values allowed within a filter group
# log2t: Boolean that specifies whether the data was log2 transformed or not
filter_intensities <- function(pg_df, method='percentage', filterGroup='DIAGNOSIS',percentageNA = 0.5, maxNA=0) {
  print("Running part 3 of 6...")
  
  # throw error if percentageNA < 0 or > 1
  if (percentageNA < 0 || percentageNA > 1) {
    stop("percentageNA must be between 0 and 1")
  }
  
  # create NARate column names for every unique diagnosis
  print("Attempting to retrieve unique filter group names...")
  diagnoses <- unique(pull(pg_df, !!sym(filterGroup)))
  print(paste("Number of unique filter groups detected:", length(diagnoses)))
  print("Successfully retrieved unique filter group names.")
  
  # filter values that do not meet conditions
  # first, create columns of NA counts and total diagnoses for each protein group and diagnosis
  print("Attempting to add MinNARate and MaxNAcount columns...")
  print(paste("Number of columns:", ncol(pg_df)))
  filtered_df <- pg_df %>%
    mutate(MinNARate = 1, MaxNAcount = 0)
  print("Successfully added MinNARate and MaxNAcount columns.")
  print(paste("Number of columns:", ncol(filtered_df)))
  
  print("Attempting to create NARate/NACount columns for every unique diagnosis...")
  print(paste("Number of columns:", ncol(filtered_df)))
  
  print("---------------------------------------> Total missing value percentage of intensity")
  total_missing_na_intensity <- (sum(is.na(filtered_df$Intensity)) / (nrow(filtered_df)))
  print(total_missing_na_intensity)
  
  for (i in diagnoses) {
    filtered_df <- filtered_df %>%
      group_by(`Protein IDs`) %>%
      mutate("NACount.{{i}}" := sum((is.na(Intensity) | Intensity == 0) & !!sym(filterGroup) == i), "NARate.{{i}}" := sum((is.na(Intensity) | Intensity == 0) & !!sym(filterGroup) == i)/sum(!!sym(filterGroup) == i))
    filtered_df <- filtered_df %>%
      group_by(`Protein IDs`) %>%
      mutate(MinNARate = min(MinNARate, !!sym(paste0('NARate."', i, '"'))), MaxNAcount = max(MaxNAcount, !!sym(paste0('NACount."', i, '"'))))
  }
  print("---------------------------------------> Showing first 6 rows of filtered_df:")
  print(head(filtered_df))
  
  # if (log2t) {
  #   for (i in diagnoses) {
  #     pg_df <- pg_df %>%
  #       group_by(`Protein IDs`) %>%
  #       mutate("NACount.{{i}}" := sum(is.na(Intensity) & !!sym(filterGroup) == i), "NARate.{{i}}" := sum(is.na(Intensity) & !!sym(filterGroup) == i)/sum(!!sym(filterGroup) == i), MinNARate = min(MinNARate, sum(is.na(Intensity) & !!sym(filterGroup) == i)/sum(!!sym(filterGroup) == i)), MaxNAcount = max(MaxNAcount, sum(is.na(Intensity) & !!sym(filterGroup) == i)))
  #   }
  # } else {
  #   for (i in diagnoses) {
  #     pg_df <- pg_df %>%
  #       group_by(`Protein IDs`) %>%
  #       mutate("NACount.{{i}}" := sum(Intensity == 0 & !!sym(filterGroup) == i), "NARate.{{i}}" := sum(is.na(Intensity) & !!sym(filterGroup) == i)/sum(!!sym(filterGroup) == i), MinNARate = min(MinNARate, sum(is.na(Intensity) & !!sym(filterGroup) == i)/sum(!!sym(filterGroup) == i)), MaxNAcount = max(MaxNAcount, sum(is.na(Intensity) & !!sym(filterGroup) == i)))
  #   }
  # }
  
  print(paste("Number of columns after adding NARate/NAcount columns:", ncol(filtered_df)))
  print("Successfully created NARate/NACount columns for every unique diagnosis.")
  
  # then, filter based on percentage or raw amount
  if (method == "percentage") {
    print("Filter by percentage selected, only keep rows whose NA rates are lower than or equal to specified...")
    print(paste("Number of rows:", nrow(filtered_df)))
    # filter retains all rows that satisfy your condition
    filtered_df <- filtered_df %>%
      filter(MinNARate <= percentageNA)
    print(paste("Number of rows after filtering:", nrow(filtered_df)))
    print(paste("Successfully removed rows with all filter groups whose NA rates are lower than specified:", percentageNA))
  } else if (method == "raw") {
    print("Filter by raw amount selected, only keep rows with whose NA counts are lower than specified...")
    print(paste("Number of rows:", nrow(filtered_df)))
    print(paste("maxNA:", maxNA))
    # filter retains all rows that satisfy your condition
    filtered_df <- filtered_df %>%
      group_by(`Protein IDs`, !!sym(filterGroup)) %>%
      filter(MaxNAcount <= maxNA)
    print(paste("Number of rows after filtering:", nrow(filtered_df)))
    print(paste("Successfully removed rows with all filter groups whose NA counts are lower than specified:", maxNA))
  } else {
    warning("Matrix was not filtered; method was not specified as 'percentage' or 'raw'")
  }
  
  print("---------------------------------------> Showing first 6 rows of filtered_df after filtering:")
  print(head(filtered_df))
  
  print("Part 3 of 6 complete, returning filtered tibble...")
  return(list("pgdata_3" = filtered_df, "total_missing_na_intensity" = total_missing_na_intensity))
}

# pg_df: tibble containing protein groups
# method: "ppca", "bpca", or "mindet"
# nPcs: number of principal components
# imputeGroup: variable that is used to separate groups for imputation
# intensity_type: raw intensity, LFQ intensity, or MS/MS count
impute <- function(pg_df, method = "ppca", nPcs=3, imputeGroup = "DIAGNOSIS") {
  print("Running part 4 of 6...")
  
  print("------------------------------------------> head(pg_df)")
  print(head(pg_df))
  # split into groups
  print("Attempting to retrieve unique filter group names...")
  diagnoses <- unique(pull(pg_df, !!sym(imputeGroup)))
  print("Successfully retrieved unique filter group names.")
  
  print("Attempting to create tibbles for imputation...")
  groupdata_df <- pg_df %>%
    dplyr::select(`Protein IDs`, id, Sample, Intensity, !!sym(imputeGroup))
  
  print("------------------------------------------> head(groupdata_df)")
  print(head(groupdata_df))
  print("groupdata_df %>% filter(id == 1)")
  print(groupdata_df %>% filter(id == 1))
  
  impute_df <- tibble(groupdata_df) %>%
    ungroup() %>%
    dplyr::select(-!!sym(imputeGroup)) %>%
    spread(Sample, Intensity) %>%
    dplyr::select(-starts_with("Intensity"), -starts_with("LFQ intensity"), -starts_with("MS/MS count"))
  print("Successfully created tibbles for imputation.")
  # impute each group, combine with impute_df
  print("Attempting to perform imputations...")
  if (method == "mindet") {
    groupdata_df <- groupdata_df %>%
      ungroup() %>%
      dplyr::select(-!!sym(imputeGroup)) %>%
      spread(Sample, Intensity) %>%
      dplyr::select(id, starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"))
    groupdata_df <- groupdata_df %>%
      impute.MinDet()
    impute_df <- impute_df %>%
      inner_join(groupdata_df)
  }
  else {
    for (i in diagnoses) {
      cur_grp <- groupdata_df %>%
        ungroup() %>%
        filter(!!sym(imputeGroup) == i) %>%
        spread(Sample, Intensity) %>%
        dplyr::select(`Protein IDs`, starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"), id)
      pca_col_names = colnames(cur_grp)[grepl("^Intensity|^LFQ intensity|^MS/MS count",colnames(cur_grp))]
      print("******************************************")
      print("Including the following columns in pca method")
      print(pca_col_names)
      print("cur_grp")
      print(cur_grp)
      print("cur_grp %>% filter(id == 1)")
      print(cur_grp %>% filter(id == 1))
      #print("head(cur_grp)")
      #print(head(cur_grp)) # terminal says "Error: Data contains rows in which all elements are 'NA'. Remove them first"
      #print("checkData")
      #checkData(cur_grp, verbose=TRUE) # when 100% get "Error: Data is not numeric" in terminal 
      # and "Error in is.infinite(data): default method not implemented for type 'list'" in site
      print("****************************************** HERE")
      #save(cur_grp, file = "/Users/wilson/Work/Chiosis - Chaperomics/cur_grp.RData")
      # problem is that all data, including id row has NA in it for one row
      # but also this makes me think that the id row is being used in the pca method
      # see https://rdrr.io/bioc/pcaMethods/man/pca.html
      # where it says, "Can also be a data frame in which case all numberic variables are used to fit the PCA."
      tryCatch({
          #pca <- pca(cur_grp, "bpca", nPcs=3, subset=pca_col_names) 
          # TODO add back, because "id" column is being used but for now let it use it so I don't have to clean data
          pc <- pca(cur_grp, "bpca", nPcs=3) 
        }, error = function(e) {
          # TODO do we want to capture the error message?
          num <- vapply(cur_grp, is.numeric, logical(1))
          Matrix <- as.matrix(cur_grp[,num])
          errorMessage <- capture.output(checkData(Matrix, verbose=TRUE), type="message")
          stop(paste(errorMessage, ". At least one of the following columns must not be NA: ", paste(pca_col_names, collapse=", "), sep=""))
        })

      #pc <- pca(cur_grp, method=method, nPcs=nPcs, subset=pca_col_names)
      impute_df <- impute_df %>%
        inner_join(completeObs(pc), copy=TRUE)
    }
  }
  
  impute_df <- impute_df %>%
    distinct(id, .keep_all = TRUE)
  
  print("impute_df %>% filter(id == 1)")
  print(impute_df %>% filter(id == 1))
  
  print("Successfully performed imputations.")
  
  # convert to long format
  # print("Attempting to convert tibble with imputations to long format...")
  # print(paste("Number of rows:", nrow(pg_df)))
  # print(paste("Number of columns:", ncol(pg_df)))
  # if (intensity_type == "Raw") {
  #   impute_df <- impute_df %>%
  #     gather("Sample", "Intensity", starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"))
  # } else if (intensity_type == "LFQ") {
  #   impute_df <- impute_df %>%
  #     gather("Sample", "LFQ intensity", starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"))
  # } else if (intensity_type == "MS2") {
  #   impute_df <- impute_df %>%
  #     gather("Sample", "MS/MS count", starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"))
  # } else {
  #   warning("type should be specified as 'Raw', 'LFQ', or 'MS2'.")
  # }
  # print(paste("Number of rows after conversion:", nrow(impute_df)))
  # print(paste("Number of columns after conversion:", ncol(impute_df)))
  # print("Successfully converted tibble with imputations to long format.")
  
  print("Part 4 of 6 complete, returning imputed tibble...")
  return (impute_df)
}

# pg_df: tibble containing protein groups
# pg_col: column containing name of protein IDs
# server: biomaRt server
# dataset: ensembl dataset
add_bm <- function (pg_df, pg_col = "Protein IDs", server = "useast.ensembl.org", dataset = "hsapiens_gene_ensembl") {
  print("Running part 5 of 6...")
  
  print("Attempting to initialize biomaRt...")
  mart <- useMart("ensembl")
  mart = useDataset(dataset, mart=mart)
  print("Successfully initialized biomaRt.")
  
  print("Attempting to retrieve information from biomaRt...")
  
  # getBM_df <- getBM(attributes=c("hgnc_symbol","description","uniprotswissprot"),filters = "uniprotswissprot",values=pull(pg_df, !!sym(pg_col)),mart=mart)
  
  pg_df_2 <- pg_df %>%
    separate_rows(`Protein IDs`, sep = ';')
  
  getBM_df <- getBM(attributes=c("hgnc_symbol","description","uniprotswissprot"),filters = "uniprotswissprot",values=pull(pg_df_2, !!sym(pg_col)),mart=mart)
  
  getBM_df <- getBM_df %>%
    group_by(uniprotswissprot) %>%
    summarize(hgnc_symbol = paste0(hgnc_symbol, collapse=';'), description = paste0(description, collapse=';'), uniprotswissprot = unique(uniprotswissprot))
  
  new_df <- pg_df_2 %>%
    left_join(getBM_df, by=c("Protein IDs" = "uniprotswissprot"))

  new_df <- new_df %>%
    dplyr::select(hgnc_symbol, description, `Protein IDs`, id) %>%
    group_by(id) %>%
    summarize(hgnc_symbol = paste0(hgnc_symbol, collapse=';'), description = paste0(description, collapse=';'), `Protein IDs` = paste0(`Protein IDs`, collapse=';'), id = unique(id)) %>%
    right_join(pg_df)
  
  print("Successfully retrieved information from biomaRt.")
  
  print("Part 5 of 6 complete, returning updated tibble...")
  
  return (new_df)
}

# pg_df: tibble with protein groups
# df_type: specified as "metadata" or "matrixx"
output_dfs <- function (pg_df, df_type) {
  print("Running part 6 of 6...")
  
  out_df <- data.frame(pg_df)
  
  if (df_type == "matrixx") {
    print("matrixx selected, attempting to create data matrix...")
    row.names(out_df) <- out_df$id
    out_df <- out_df %>%
      dplyr::select(starts_with("Intensity"), starts_with("LFQ.intensity"), starts_with("MS.MS.count"))
    colnames(out_df) <- gsub("\\bIntensity.", "", colnames(out_df))
    colnames(out_df) <- gsub("\\bLFQ.intensity.", "", colnames(out_df))
    colnames(out_df) <- gsub("\\bMS.MS.count.", "", colnames(out_df))
    print("Data matrix successfully created.")
  } else if (df_type == "metadata") {
    print("metadata selected, attempting to create metadata...")
    row.names(out_df) <- out_df$id
    out_df <- out_df %>%
      dplyr::select(!starts_with("Intensity") & !starts_with("LFQ.intensity") & !starts_with("MS.MS.count"))
    print("Metadata successfully created.")
  } else {
    warning("Parameter 'df_type' should be specified as 'metadata' or 'matrixx'!")
  }
  
  print("Part 6 of 6 complete, returning requested tibble...")
  return(out_df)
}