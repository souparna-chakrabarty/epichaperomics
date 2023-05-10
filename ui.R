library(shiny)
library(shinyjs)
library(plotly)
library(shinybusy)

fluidPage(
  useShinyjs(),  # Include shinyjs
  inlineCSS(list(.disabledTab = "color: currentColor;
                 cursor: not-allowed;
                 opacity: 0;
                 text-decoration: none;")),
  navbarPage('Epichaperomics - Gabriela Chiosis Lab',
    id = 'navBar',
    tabPanel('Input',
      value = 'inputPanel',       
      sidebarLayout(
        sidebarPanel(
          h3('Input settings'),
          textInput('projectName', 'Project name'),               
          selectInput('type', 'Input type', c('MS2', 'LFQ', 'Raw'), 'Raw'),
          fileInput('sampleManifest', 'Sample manifest file (*.tsv)', accept = '.tsv'),
          fileInput('pgData', 'Protein groups file (*.tsv)', accept = '.tsv'),
          #uiOutput('runAnalysesButton')
          uiOutput('uploadDataButton')
        ), # end sidebarPanel
        mainPanel(
          #uiOutput('sampleCount'),
          verbatimTextOutput("validateTextArea", placeholder = FALSE),
          verbatimTextOutput("validateTextArea2", placeholder = FALSE),
          #h3("Introduction:"),
          p("‘Epichaperomics’ is an emerging ‘omics platform to study context-dependent protein-protein interaction (PPI) network (i.e., interactome network) perturbations in complex diseases. It enables proteome-wide identification, analysis, and modulation of context-specific PPI dysfunctions in native (unengineered) cells and tissues. To facilitate analysis of epichaperomics-derived large-scale datasets at single cell and regional levels, it is essential to organize system-level functions into an accessible “easy-to-use” data presentation graphics format for the scientific and biomedical communities to mine for discovery, validation, and hypothesis testing, among others."),
          p("Herein, we make available an open resource with a bioinformatics pipeline, ‘epichaperomics shiny app’ that members of the broad scientific community – with or without computational background – can access freely and without restriction."),
          p("The outcome is a dynamic site with interactive network plots and complex analysis input screens. This allows end users to change analysis parameters and interact with plots, individual proteins, and pathways, to analyze epichaperomics datasets to inform disease biology, as well as help in the development of precision medicine for treatment and prevention. The app aims to provide a platform for researchers with limited programming skills to facilitate rapid data processing and visualization."),
          h3("About the tool"),
          p(paste("Epichaperomics app is developed by R shiny (Version ", packageVersion("shiny"), ") and is free and open to all users with no login requirement. It can be readily accessed by all popular web browsers, including Google Chrome, Mozilla Firefox, and Safari. In the app, the majority of analysis and plotting functions are written in R programming language and popular packages, to produce high-quality images, similar to running bona fide R codes minus the coding process. We expect feedback from active users to refine the app moving forward.", sep = "")),
          h3("Citations"),
          p("Inda, M.C., Joshi, S., Wang, T. et al. The epichaperome is a mediator of toxic hippocampal stress and leads to protein connectivity-based dysfunction. Nat Commun 11, 319 (2020).", a("https://doi.org/10.1038/s41467-019-14082-5.", href="https://doi.org/10.1038/s41467-019-14082-5"), " PMCID: PMC6965647."),
          p("Joshi, S., Gomes, E.D., Wang, T. et al. Pharmacologically controlling protein-protein interactions through epichaperomes for therapeutic vulnerability in cancer. Commun Biol 4, 1333 (2021). ", a("https://doi.org/10.1038/s42003-021-02842-3.", href="https://doi.org/10.1038/s42003-021-02842-3"), " PMCID: PMC8617294."),
          p("Ginsberg SD, Sharma S, Norton L, Chiosis G. Targeting stressor-induced dysfunctions in protein-protein interaction networks via epichaperomes. Trends Pharmacol Sci. 2023 Jan;44(1):20-33. doi: 10.1016/j.tips.2022.10.006. PMCID: PMC9789192."),
          p("Ginsberg SD, Neubert TA, Sharma S, Digwal CS, Yan P, Timbus C, Wang T, Chiosis G. Disease-specific interactome alterations via epichaperomics: the case for Alzheimer's disease. FEBS J. 2022 Apr;289(8):2047-2066. doi: 10.1111/febs.16031. PMCID: PMC8611103.")
        ) # end mainPanel
      ), # end sidebarLayout
    ), # end input tab panel
    tabPanel('Pre-processing',
      value = 'quantitativeAnalysesPanel',
      sidebarLayout(
        sidebarPanel(
          h3('Settings'),
          checkboxInput('log2t', 'Log-transformation', T),
          checkboxInput('removeContaminant', 'Remove contaminant', F),
          checkboxInput('removeReverse', 'Remove reverse', T),
          checkboxInput('oneProteinPerRow', 'One protein per row', F),
          selectInput('filteringMethod', 'Filtering', c('10%' = 0.1, '20%' = 0.2, '30%' = 0.3, '40%' = 0.4, '50%' = 0.5, '60%' = 0.6, '70%' = 0.7, '80%' = 0.8, '90%' = 0.9, '100%' = 1.0, 'other' = -1), 0.5),
          uiOutput('filteringRawMaxValueInput'),
          selectInput('imputeMethod', 'Impute method', c('ppca', 'bpca', 'mindet'), 'bpca'), # TODO make mindet default when it is working
          # TODO what is min and max for imputeNPC?
          numericInput('imputeNPC', 'Number of principal components to use', 3, 0),
          
          uiOutput('downloadProcessedDataButton'),
          uiOutput('runAnalysesButton'),
          uiOutput('goToDiffAnalysesButton')
        ), # end sidebarPanel
        mainPanel(
          add_busy_spinner(spin = "fading-circle"),
          uiOutput('processedDataTableHeader'),
          uiOutput("numSamples", placeholder = FALSE),
          uiOutput("numProteinGroups", placeholder = FALSE),
          uiOutput("totalMissingNAIntensity", placeholder = FALSE),
          #div(tableOutput('processedDataTable'))
        ) # end mainPanel
      ), # end sidebarLayout
    ), # end pre-processing tab panel
    tabPanel('Differential Analysis',
      value = 'differentialAnalysisPanel',
      sidebarLayout(
        sidebarPanel(
          h3('Settings'),
          selectInput("contrast", "Contrast", choices = NULL),
          selectInput("groupASelection", "Group A", choices = NULL),
          selectInput("groupBSelection", "Group B", choices = NULL),
          selectInput("covariates", "Covariates (select all that apply)", choices = NULL, multiple = TRUE),
          uiOutput('downloadDiffAnalysisDataButton'),
          uiOutput('runDiffAnalysesButton'),
          uiOutput('goToGSEAButton')
        ), # end sidebarPanel
        mainPanel(
          add_busy_spinner(spin = "fading-circle"),
          textOutput("modelFormula"), # TODO highlight this and give it more space
          plotlyOutput("volcanoPlot"),
        ) # end mainPanel
      ), # end sidebarLayout
    ), # end quantitative analysis tab panel
    tabPanel('Gene Set Enrichment Analysis',
      value = 'geneSetEnrichementAnalysisPanel',
      sidebarLayout(
        sidebarPanel(
          h3('Settings'),
          selectizeInput("gmtSets", "Select GMT library", choices = NULL, selected=NULL),
          
          selectInput("numGeneSets", "Number of gene sets to display", choices = c("20", "50", "100", "all"), selected="20"),
          selectInput('minGeneSetSize', 'Minimum gene set size', c(2, 3, 4, 5, 6, 7, 8, 9, 10), 10),
          selectInput('maxGeneSetSize', 'Maximum gene set size', c(50, 100, 250, 500, 750, 1000, 999999), 500),
          checkboxInput('useAdjustedPValue', 'Use adjusted p-value', F),
          
          actionButton('runGSEAButton',"Run GSEA"),
          
          uiOutput('downloadGSEADataButton'),
        ), # end sidebarPanel
        mainPanel(
          add_busy_spinner(spin = "fading-circle"),
          plotlyOutput("gseaPlot"),
        ) # end mainPanel
      ), # end sidebarLayout
    ), # end gene set enrichment analsys panel
  ) # end navbarPage
)