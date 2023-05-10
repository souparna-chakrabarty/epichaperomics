# TODO do we need all of these libraries?
library(tidyverse)
library(limma)
library(plotly)
source("GSEAenricher.R")

get_top2 <- function(matrixx, metadata, doe, 
                     contrast="Gender", 
                     covariates=c("AGE", "Braak.stage", "Gel.Batch"), 
                     comparisonX="Female", 
                     comparisonY="Male")
  {
  doe_data <- tibble(doe)
  
  #Gender <- factor(make.names(doe$Gender))
  #Age <- factor(make.names(doe$AGE))
  #Age <- doe$AGE
  #Braak <- factor(make.names(doe$Braak.stage))
  #Batch <- factor(make.names(doe$Gel.Batch))
  #print(typeof(Age))
  
  #print("MEW TODO figure this out: ")
  #print(Gender)
  #print(Age)
  #print(Braak)
  #print(Batch)
  #print(GenderFemale)
  #print(GenderMale)
  
  formulaString <- "~0"
  
  assign(contrast, factor( make.names( doe_data[[contrast]] )))
  
  #factors = list()
  for (covariate in covariates) {
    
    #print(class( doe_data[[col]] ))
    #print(typeof(doe_data[[col]]))
    # TODO Gel.Batch should be factor(make.names()) -- how do we know?
    # TODO what about other columns?
    # TODO will input files definitely have valid column names or else we need make.names(col) and need to keep track of that (list or named vector)
    #print(col)
    if( class( doe_data[[covariate]] ) == "character" || covariate == "Gel.Batch") {
      #print("*********** COLUMN IS CHAR OR gel.batch *******************")
      assign(covariate, make.names(doe_data[[covariate]] ))
    } else {
      #print("*********** COLUMN IS SOMETHING ELSE *******************")
      assign(covariate, doe_data[[covariate]])
    }
  }
  formulaString <- paste("~0+", contrast, "+", paste(covariates, collapse="+"), sep = "")
  formulaString <- gsub('[+]{1}$', '', formulaString) # if formula string ends in '+', remove '+'
  print(formulaString)
  comparisonString <- paste(contrast, make.names(comparisonX), " - ", contrast, make.names(comparisonY), sep = "")
  
  formula <- as.formula(formulaString)
  
  #formula <- create.formula(outcome.name = NULL, input.names = cols)
  #print(typeof(~0+Gender+Age+Braak+Batch))
  #print(typeof(formula))
  #print(typeof(formula$formula))
  #print(~0+Gender+Age+Braak+Batch)
  #print(formula$formula)
  
  #Gender <- doe$Gender
  #Age <- doe$AGE
  #Braak <- doe$Braak.stage
  #Batch <- doe$Gel.Batch
  #print("***************** FORMULA STRING *************************")
  #print(formulaString)
  #print("***************** COMPARISON STRING *************************")
  #print(comparisonString)
  
  design <- model.matrix(formula)
  #print(design)
  #contrast <- makeContrasts(SelectedComparison = "GenderFemale - GenderMale", levels = design)
  # following https://support.bioconductor.org/p/27900/
  myargs = list(SelectedComparison = comparisonString, levels = design)
  #myargs = list(SelectedComparison = "GenderFemale - GenderMale", levels = design)
  mycontrast <- do.call(makeContrasts, myargs)
  
  fit <- lmFit(matrixx, design) 
  # fit <- eBayes(fit)
  fit2 <- contrasts.fit(fit, mycontrast) 
  fit2 <- eBayes (fit2)
  
  # top <- topTable(fit, coef = "GenderComparison", adjust.method = "BH",n=Inf, sort.by = "P")
  # BH = Benjamini & Hochberg/FDR adjust (adjust the p-values for multiple testing)
  top2 <- topTable(fit2, coef = "SelectedComparison", adjust.method = "BH", n=Inf, sort.by = "P")
  # top <- merge(top, metadata, by=0, all=TRUE)
  top2 <- merge(top2, metadata, by=0, all=TRUE)
  top2$Shown.ID <- top2$hgnc_symbol
  top2$Shown.ID[top2$Shown.ID==""] <- top2$Protein.IDs[top2$Shown.ID==""]
  #top2$setting <- index
  
  return(top2)
}

run_gsea <- function(top2, gmt.sel, group, minimumGSSize, maximumGSSize) {
  # "group_name" is displayed in the plot  
  load("gmt-go.Rdata")
  print("head(gmt.go):")
  print(head(gmt.go))
  print("head(gmt.sel)")
  print(head(gmt.sel))
  print("calling do_GSEA")
  g1 <- do_GSEA(top2, 1, gmt.sel, group_name = group, Protein.ID.Column ="hgnc_symbol", minimumGSSize=minimumGSSize, maximumGSSize=maximumGSSize)
  
  #g1$df$setting <- index
  return(g1)
}

plot_volcano <- function(top2) {
  print("********************** VOLCANO! ******************")
  # see https://www.biostars.org/p/214100/
  #table(top2$heading)
  print(colnames(top2))
  top2 <- top2[c("hgnc_symbol", "logFC", "adj.P.Val", "P.Value")]
  #print(head(top2))
  # change the grouping for the entries with significance but not a large enough Fold change
  top2["group"] <- "NotSignificant"
  #top2[which(top2['adj.P.Val'] >= 0.05 & top2['adj.P.Val'] < 0.5),"group"] <- "AlmostSignficant"
  top2[which(-log10(top2$adj.P.Val) < 0.05),"group"] <- "Significant"
  print(head(top2))
  print("*********************** ABOUT TO PLOT VOLCANO")
  p <- plot_ly(data = top2, x = top2$logFC, y = -log10(top2$P.Value), text = top2$hgnc_symbol, mode = "markers", color = top2$group) %>% 
    layout(title ="Volcano Plot",
           xaxis = list(title = "log fold change"),
           yaxis = list(title = "-log10(p-value)"),
           showlegend=TRUE)
  p
}

plot_gsea <- function(g1, gene_sets_to_display="20", useAdjustedPValue=FALSE) {
  
  #print("**************************************************")
  #print(colnames(g1$df))
  #print("**************************************************")
  
  g1_filtered <- g1$df %>%
   filter(pvalue < 0.05)
  
  # TODO should we use p.adjust?  for me that leaves no data
  if (useAdjustedPValue) {
    print("********************** using adjusted p-value ****************************")
    g1_filtered <- g1$df %>%
      filter(p.adjust < 0.05)
  }
 
  print(length(unique (g1_filtered$ID)))
  #print(gene_sets_to_display)
  if (gene_sets_to_display != "all") {
    #print("all not selected so filtering gene sets")
    # get lowest p-values and display those gene sets only
    g1_filtered <- g1_filtered %>% slice_min(pvalue, n = as.numeric(gene_sets_to_display))
  }
  #print(length(unique (g1_filtered$ID)))
  
  #print(unique (g1$df$ID))
  #hcal_pea <- max ((length (unique (g1$df$ID)) * 16.5 + 25 + 10 + 100), 500)
  hcal_pea <- max ((length (unique (g1_filtered$ID)) * 16.5 + 25 + 10 + 100), 500)
  
  ## TAI: Added these three line to resort the pathway (must keep the 2nd and the 3rd in a row to make sure the order in levels are correct!!!)
  g1_filtered <- g1_filtered [order (g1_filtered$NES, decreasing = F), ] 
  g1_filtered$ID <- factor (g1_filtered$ID)
  g1_filtered$ID <- factor (g1_filtered$ID, levels = unique (g1_filtered$ID))
  
  # MEW removed &nbsp;<b>Core enrichment:</b>', core_enrichment, '<br>' from pop up
  p2 <- ggplot (data=g1_filtered, aes (x=NES, y=ID)) +
    geom_point(aes (size=setSize, fill=NES, color=-log10(pvalue), stroke=0.5, text=paste('&nbsp;<b>ID:</b>', ID, '<br>',
                                                                                         '<b>Comparison:</b>', Group,'<br>',
                                                                                         '<b>Set size:</b>', setSize,'<br>',
                                                                                         '<b>NES:</b>',NES,'<br>',
                                                                                         '<b>p-value:</b>',pvalue,'<br>',
                                                                                         '<b>adjusted p-value:</b>',p.adjust,'<br>',
                                                                                         '<b>Rank:</b>',rank,'<br>',
                                                                                         '<b>Leading edge:</b>',leading_edge,'<br>'))) +
    scale_fill_gradient2(limits=c(-max(abs(g1_filtered$NES)), max(abs(g1_filtered$NES))), low = "#0199CC", mid = "white", high = "#FF3705", midpoint = 0) +
    scale_color_gradient2(limits=c(0, max(abs(log10(g1_filtered$pvalue)))), low = "white", mid = "white", high = "black", midpoint = 1) +
    theme_minimal () + 
    ylab("") + 
    theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_text(size = 8, angle = 0, hjust = 1, face = "plain"))# +
    # scale_y_discrete(labels = function(y) {
    #   y = gsub('^GO:[0-9]{2,}-', '', y) # remove the GO id (e.g. GO:0008152)    ### Tai: GMT libraries are not always of Gene Ontology, suggest to remove
    #   stringr::str_trunc(y, 30)
    # })

  gp.pt.2 <- ggplotly(p2, tooltip = "text", height = hcal_pea) %>% 
    layout (yaxis=(list(automargin = F)), margin=list (l=200))
}
