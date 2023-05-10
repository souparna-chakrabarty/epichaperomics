library(shiny)
library(plotly)
library(tidyverse)
library(fuzzyjoin)
library(scales)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
options(repos = BiocManager::repositories())
source("chaperomics.R")
source("gsea.R")
source ("./lib/func_readgmt.R") ## Based the same function from clusterProfiler package

options(shiny.maxRequestSize=35*1024^2)

contaminants <- read_tsv("Contaminants.tsv")

# You have to factor vector in order to do a comparision or it won't be in the design matrix correctly -- an individual column is created for each level in the factor

# TODO get list of authors for code
# TODO what is project id for?
DEFAULT_MAX_RAW_FILTER = 10

# anything defined here has a scope of the session
function(input, output, session) {


  disableTab <- function(tabId) {
    tabSelector = paste("#navBar li a[data-value=", tabId, "]")
    shinyjs::disable(selector = tabSelector)
    shinyjs::addClass(selector = tabSelector, class = "disabledTab")
  }

  enableTab <- function(tabId) {
    tabSelector = paste("#navBar li a[data-value=", tabId, "]")
    shinyjs::enable(selector = tabSelector)
    shinyjs::removeClass(selector = tabSelector, class = "disabledTab")
  }

  disableTab("quantitativeAnalysesPanel")
  disableTab("differentialAnalysisPanel")
  disableTab("geneSetEnrichementAnalysisPanel")

  pipelineResults <- NULL
  pgdata <- NULL
  manifest <- NULL
  top2Result <- reactiveValues(top2=NULL)
  gseaResult <- reactiveValues(g1=NULL)


  observeEvent(input$goToDiffAnalysis, {
    updateTabsetPanel(session, 'navBar', selected = 'differentialAnalysisPanel')
    enableTab("differentialAnalysisPanel")
  })

  observeEvent(input$goToGSEA, {
    updateTabsetPanel(session, 'navBar', selected = 'geneSetEnrichementAnalysisPanel')
    gmt.list <- list.files(path="./gmt/")
    updateSelectizeInput(inputId = 'gmtSets', choices = c('',gmt.list), selected = 'h.all.v2022.1.Hs.symbols.gmt')

    enableTab("geneSetEnrichementAnalysisPanel")
  })

  observeEvent(input$uploadData, { # wait until next button clicked
    # make sure other tabs are disabled again
    disableTab("quantitativeAnalysesPanel")
    disableTab("differentialAnalysisPanel")
    disableTab("geneSetEnrichementAnalysisPanel")
    req(input$sampleManifest)
    req(input$pgData)

    sampleManifestFile <- input$sampleManifest
    if (is.null(sampleManifestFile))
      return(NULL)

    pgDataFile <- input$pgData
    if (is.null(pgDataFile))
      return(NULL)

    sampleManifestFileExt <- tools::file_ext(input$sampleManifest$name)
    pgDataFileExt <- tools::file_ext(input$pgData$name)
    print("**************************************** ext = ?")
    print(sampleManifestFileExt)
    print(sampleManifestFileExt == "tsv")
    # TODO we probably can't have another validateTextArea but we do need validate to be contained inside a rendertext or something
    output$validateTextArea2 <- renderText({
      validate(need(sampleManifestFileExt == "tsv" || sampleManifestFileExt == "TSV", "Invalid sample manifest file. Please upload a .tsv file"),
               need(pgDataFileExt == "tsv" || pgDataFileExt == "TSV", "Invalid protein groups file. Please upload a .tsv file"))

      # TODO make sure this is still local to the session
      manifest <<- tibble(read_tsv(sampleManifestFile$datapath))
      pgdata <<- tibble(read_tsv(pgDataFile$datapath))
      required_columns <- c('Sample.ID2', 'DIAGNOSIS')
      column_names <- colnames(manifest)

      output$validateTextArea <- renderText({
        validate(
          need(all(required_columns %in% column_names),
               paste("Required columns are missing from sample manifest file:",
                     paste(required_columns[!(required_columns %in% column_names)],
                           collapse = ", "))),
          need(nrow(manifest) > 0, "No data found in uploaded sample manifest file."),
          need(nrow(pgdata) > 0, "No data found in uploaded protein groups file."),
        )
        updateTabsetPanel(session, "navBar", selected="quantitativeAnalysesPanel")

        enableTab("quantitativeAnalysesPanel")
        output$processedDataTableHeader <- renderUI({
          h3("Uploaded Data Summary")
        })

        output$numSamples <- renderText({
          paste("<b>Number of records in samples file:</b>", nrow(manifest))
        })

        output$numProteinGroups <- renderText({
          paste("<b>Number of rows in protein groups file:</b>", nrow(pgdata))
        })

        output$totalMissingNAIntensity <- renderText({
          paste("")
        })

        "" # don't display anything in validateTextArea if valid in validate area - this needs to be at the end
      })
      "" # don't display anything in validateTextArea2 if valid in validate area - this needs to be at the end
    })
  })

  observeEvent(input$runAnalysis, { # wait until next button clicked
    # make sure other tabs are disabled again
    disableTab("differentialAnalysisPanel")
    disableTab("geneSetEnrichementAnalysisPanel")
    tryCatch({
      pgdata_1 <- runMQ(pgdata,
                        removeContaminant = input$removeContaminant,
                        removeReverse = input$removeReverse,
                        type = input$type,
                        log2t = input$log2t,
                        oneProteinPerRow = input$oneProteinPerRow,
                        contaminantList = contaminants)

      # TODO: pass intensity_type to combine_matrix, once it is ready to use it
      #pgdata_2 <- combine_matrix(pgdata_1, manifest, sample.id.col='Sample.ID2', intensity_type=input$type)
      pgdata_2 <- combine_matrix(pgdata_1, manifest, sample.id.col='Sample.ID2')

      # default is to filter NA values per group by percentage ("percentage")
      # if percentage given, input$filteringMethod will give us minimum percentage of non-NA entries needed
      # within a filter group in order to keep the protein group
      percentageNA = input$filteringMethod
      print("****************************************")
      print("percentageNA")
      print(percentageNA)
      maxNA = DEFAULT_MAX_RAW_FILTER
      method = "percentage" # default method
      if (input$filteringMethod == -1) {
        # user has selected "other" in "Filtering" drop down
        # filter NA values per group by raw amount ("raw")
        method = "raw" # TODO get this value
        percentageNA = 0 # this is not used, use maxNA
        # user has given "raw" value to filter NA
        if (!is.null(input$filteringRawMaxValue))
          maxNA = input$filteringRawMaxValue
      }

      # TODO get filterGroup from sample file?
      filter_intensities_return_values_list <- filter_intensities(pgdata_2, method = method, percentageNA = percentageNA, maxNA = maxNA, filterGroup = 'DIAGNOSIS')
      pgdata_3 <- filter_intensities_return_values_list$pgdata_3
      total_missing_na_intensity <- filter_intensities_return_values_list$total_missing_na_intensity

      # TODO get imputeGroup from sample file?
      pgdata_4 <- impute(pgdata_3, nPcs = input$imputeNPC, method = input$imputeMethod, imputeGroup = "DIAGNOSIS")
      #pgdata_4 <- impute2(pgdata_3, nPcs = input$imputeNPC, method = input$imputeMethod)


      # TODO what are inputs for this?
      pgdata_5 <- add_bm(pgdata_4, pg_col = "Protein IDs", server = "useast.ensembl.org", dataset = "hsapiens_gene_ensembl")

      matrixx <- output_dfs(pgdata_5, "matrixx")
      metadata <- output_dfs(pgdata_5, "metadata")

      # only allow users to select non-numeric columns with more than one unique value for contrast
      non_numeric_choices = Filter(function(x) length(unique(x))>1, Filter(negate(is.numeric), manifest))
      contrast_choices = grep("Sample.ID",
                              colnames(non_numeric_choices),
                              ignore.case = TRUE,
                              value = TRUE,
                              fixed = FALSE, # TRUE invalidates ingore.case = TRUE
                              invert = TRUE)
      updateSelectInput(inputId = "contrast", choices = contrast_choices)

      # TODO check scope with Caitlin: https://shiny.rstudio.com/articles/scoping.html
      pipelineResults <<- list(processedPGData=pgdata_5, matrixx=matrixx, metadata=metadata, manifest=manifest, contrast_choices=contrast_choices)

      output$processedDataTableHeader <- renderUI({
        h3("Processed Data Summary")
      })

      #output$processedDataTable <- renderTable({
      #  pipelineResults[['processedPGData']]
      #})

      output$numSamples <- renderText({
        paste("<b>Number of samples:</b>", length(unique(manifest$Sample.ID2)))
      })

      output$numProteinGroups <- renderText({
        paste("<b>Number of protein groups:</b>", length(pipelineResults[['processedPGData']]$'Protein IDs'))
      })

      output$totalMissingNAIntensity <- renderText({
        paste("<b>Percent missing intensity values:</b>", label_percent(accuracy=.01)(total_missing_na_intensity))
      })

      output$goToDiffAnalysesButton <- renderUI({
        actionButton("goToDiffAnalysis", "Next")
      })

      output$downloadProcessedDataButton <- renderUI({
        downloadButton("downloadProcessedData", "Export processed data table")
      })

      output$downloadProcessedData <- downloadHandler(
        filename = function() {
          "pgdata_processed.tsv"
        },
        content = function(file) {
          write_tsv(pipelineResults[['processedPGData']], file)
        }
      )
    },
    #warning = function(warn){
    #  showNotification(paste0(warn), type = 'warning', duration = NULL, closeButton = TRUE)
    #},
    error = function(err){
      showNotification(paste0(err), type = 'err', duration = NULL, closeButton = TRUE)
    })
  })

  has_significant_data <- function(g1, useAdjustedPValue) {
    if (useAdjustedPValue) {
      print("Using adjusted p-value")
      g1_filtered <- g1$df %>% filter(p.adjust < 0.05)
      print(length(g1_filtered$ID))
      length(g1_filtered$ID) > 0
    } else {
      print("Using p-value")
      g1_filtered <- g1$df %>% filter(pvalue < 0.05)
      print(length(g1_filtered$ID))
      length(g1_filtered$ID) > 0
    }
  }

  observeEvent(input$runDiffAnalysis, { # wait until run button clicked on diff analysis page
    # make sure last tab is disabled again
    disableTab("geneSetEnrichementAnalysisPanel")

    req(pipelineResults)
    req(input$groupASelection)
    req(input$groupBSelection)
    req(input$runDiffAnalysis)

    tryCatch({
      top2 <- get_top2(pipelineResults[['matrixx']],
                pipelineResults[['metadata']],
                pipelineResults[['manifest']],
                contrast=input$contrast,
                covariates=input$covariates,
                comparisonX = input$groupASelection,
                comparisonY = input$groupBSelection)

      top2Result$top2 <- top2

      # make plot
      output$volcanoPlot <- renderPlotly({
        plot_volcano(top2)
      })

      # add download button
      output$downloadDiffAnalysisDataButton <- renderUI({
        downloadButton("downloadDiffAnalysisData", "Export differential analysis data table")
      })

      output$downloadDiffAnalysisData <- downloadHandler(
        filename = function() {
          "diff_analysis.tsv"
        },
        content = function(file) {
          write_tsv(top2, file)
        }
      )

      output$goToGSEAButton <- renderUI({
        actionButton("goToGSEA", "Next")
      })

    },
    #warning = function(warn){
    #  showNotification(paste0(warn), type = 'warning', duration = NULL, closeButton = TRUE)
    #},
    error = function(err){
      showNotification(paste0(err), type = 'err', duration = NULL, closeButton = TRUE)
    })
  })

  observeEvent (input$runGSEAButton, {
    req(input$minGeneSetSize)
    req(input$maxGeneSetSize)
    print("************************ input$minGeneSetSize")
    print(input$minGeneSetSize)
    print("************************ input$maxGeneSetSize")
    print(input$maxGeneSetSize)
    # now we can go to GSEA

    tryCatch({

          gmt.sel <- read.gmt(paste0 ("./gmt/", input$gmtSets))
          print("*********** gmt.sel ***************")
          if (is.null(top2Result$top2)) {
            showNotification("Must run DEA first")
            return (NULL)}

          g1 <- run_gsea(top2Result$top2, gmt.sel, group=paste(input$groupASelection, "vs", input$groupBSelection), input$minGeneSetSize, input$maxGeneSetSize)

          gseaResult$g1 <- g1

          # add download button
          output$downloadGSEADataButton <- renderUI({
            downloadButton("downloadGSEAData", "Export GSEA data table")
          })

          output$downloadGSEAData <- downloadHandler(
            filename = function() {
              "gsea.tsv"
            },
            content = function(file) {
              write_tsv(g1$df, file)
            }
          )

          # make plot
          output$gseaPlot <- renderPlotly({

            print("********************* pvalue; adjusted p-value ************************")
            g1_filtered <- g1$df %>% filter(pvalue < 0.05)
            print(length(g1_filtered$ID))
            g1_filtered_adjusted <- g1$df %>% filter(p.adjust < 0.05)
            print(length(g1_filtered_adjusted$ID))
            print(if (input$useAdjustedPValue) "use adjusted p-value" else "use p-value")

            validate(
              need(has_significant_data(g1, input$useAdjustedPValue), "No significant results to plot.")
            )

            plot_gsea(g1, gene_sets_to_display=input$numGeneSets, useAdjustedPValue=input$useAdjustedPValue)

          })
    },

      error = function(err){
        showNotification(paste0(err), type = 'err', duration = NULL, closeButton = TRUE)
      })

    })




  output$filteringRawMaxValueInput <- renderUI({
    if (input$filteringMethod != -1)
      return(NULL)
    numericInput("filteringRawMaxValue", "Filtering raw value", DEFAULT_MAX_RAW_FILTER, 0)
   })

  output$uploadDataButton <- renderUI({
    req(input$sampleManifest)
    req(input$pgData)
    actionButton("uploadData", "Upload")
  })

  output$runAnalysesButton <- renderUI({
    req(manifest)
    actionButton("runAnalysis", "Run")
  })

  output$runDiffAnalysesButton <- renderUI({
    actionButton("runDiffAnalysis", "Run")
  })

  observeEvent(input$contrast, {
    req(pipelineResults)
    updateSelectInput(inputId = "groupASelection",
                      choices = unique(pipelineResults[['manifest']][input$contrast[[1]]]),
                      selected = unique(pipelineResults[['manifest']][input$contrast[[1]]])[[1]][[1]])
    updateSelectInput(inputId = "groupBSelection",
                      choices = unique(pipelineResults[['manifest']][input$contrast[[1]]]),
                      selected = unique(pipelineResults[['manifest']][input$contrast[[1]]])[[1]][[2]]) # TODO make sure we have at least 2 items

    contrast_choices = pipelineResults[['contrast_choices']]
    covariate_choices = contrast_choices[ contrast_choices != input$contrast[[1]] ]

    updateSelectInput(inputId = "covariates",
                      choices = covariate_choices)

  })

  output$modelFormula <- renderText({
    req(input$contrast)
    model <- paste("Model: ", paste("~ ", input$contrast), " + ", paste(input$covariates, collapse=" + "))
    gsub('\\s*[+]\\s*$', '', model) # if formula string ends in ' + ', remove ' + '
  })

}
