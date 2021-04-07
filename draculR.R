library(dplyr)
library(plyr)
library(ggplot2)
library(shiny)
library(ggrepel)
library(scales) # For percent_format()
library(tidyr)
library(magrittr)
library(reshape)
library(edgeR)
library(readr)
library(DT)
library(psych)
library(shinyhelper)

##### User defined functions #####

# negate %in%
`%notin%` <- Negate(`%in%`)

# function (x) to count the number of non-zero records in each column (ie per sample)
nonzero <- function(x) sum(x != 0)

##################################

# import our distribution difference dataframe
plotData_distDiff_dCq <- read_csv("www/plotData_distDiff_dCq.csv")

# import the rank and distribtuion difference data for GSE153813
rank_GSE153813 <- read_csv(file = "www/rank_GSE153813.csv")
unlistDistributionDifference_GSE153813 <- read_csv(file = "www/unlist_distributionDifference_GSE153813.csv")
GSE153813_info <- paste("The data containined in GSE153813 were obtained from",
                        " NCBI GEO. There is no associated publication" , sep = "")

countsRaw_GSE153813 <- read_csv(file = "www/counts_raw_GSE153813.csv")

rankDist_GSE153813 <- dplyr::full_join(rank_GSE153813, unlistDistributionDifference_GSE153813, by = "samplename") %>% 
  dplyr::mutate(., project = rep("GSE153813", nrow(.)))

# import the rank and distribtuion difference data for GSE118038
rank_GSE118038 <- read_csv(file = "www/rank_GSE118038.csv")
unlistDistributionDifference_GSE118038 <- read_csv(file = "www/unlist_distributionDifference_GSE118038.csv")
GSE118038_info <- paste("The data containined in GSE118038 were obtained from",
                        " NCBI GEO and originate from the article \"A preliminary study of micro-RNAs",
                        " as minimally invasive biomarkers for the diagnosis of prostate cancer patients.\"",
                        " J Exp Clin Cancer Res 2021 Feb 23;40(1):79.", " Notably the data was ",
                        " obtained in a normalised table hence the \"Mature miRNA as a function of Read Counts\"",
                        " plot is not as expected.", sep = "")

countsRaw_GSE118038 <- read_csv(file = "www/counts_raw_GSE118038.csv")

rankDist_GSE118038 <- dplyr::full_join(rank_GSE118038, unlistDistributionDifference_GSE118038, by = "samplename") %>% 
  dplyr::mutate(., project = rep("GSE118038", nrow(.)))

# import the rank and distribtuion difference data for GSE105052
rank_GSE105052 <- read_csv(file = "www/rank_GSE105052.csv")
unlistDistributionDifference_GSE105052 <- read_csv(file = "www/unlist_distributionDifference_GSE105052.csv")
GSE105052_info <- paste("The data containined in GSE105052 were obtained from",
                        " NCBI GEO and originate from the article \"Small ",
                        " RNA-seq analysis of circulating miRNAs to identify phenotypic variability in",
                        " Friedreich's ataxia patients.\" Sci Data 2018 Mar 6;5:180021.", sep = "")

countsRaw_GSE105052 <- read_csv(file = "www/counts_raw_GSE105052.csv")

rankDist_GSE105052 <- dplyr::full_join(rank_GSE105052, unlistDistributionDifference_GSE105052, by = "samplename") %>% 
  dplyr::mutate(., project = rep("GSE105052", nrow(.)))

# import the rank and distribtuion difference data for GSE151341
rank_GSE151341 <- read_csv(file = "www/rank_GSE151341.csv")
unlistDistributionDifference_GSE151341 <- read_csv(file = "www/unlist_distributionDifference_GSE151341.csv")
GSE151341_info <- paste("The data containined in c were obtained from",
                        " NCBI GEO and originate from the article \"Sequencing",
                        " identifies a distinct signature of circulating microRNAs",
                        " in early radiographic knee osteoarthritis.\"",
                        " Osteoarthritis Cartilage 2020 Nov;28(11):1471-1481.", sep = "")

countsRaw_GSE151341 <- read_csv(file = "www/counts_raw_GSE151341.csv")

rankDist_GSE151341 <- dplyr::full_join(rank_GSE151341, unlistDistributionDifference_GSE151341, by = "samplename") %>% 
  dplyr::mutate(., project = rep("GSE151341", nrow(.)))

# global objects for imported data calculations

classifier_miRs <- data.frame(
  SYMBOL = c(
    "hsa-miR-106b-3p",
    "hsa-miR-140-3p",
    "hsa-miR-142-5p",
    "hsa-miR-532-5p",
    "hsa-miR-17-5p",
    "hsa-miR-19b-3p",
    "hsa-miR-30c-5p",
    "hsa-miR-324-5p",
    "hsa-miR-192-5p",
    "hsa-miR-660-5p",
    "hsa-miR-186-5p",
    "hsa-miR-425-5p",
    "hsa-miR-25-3p",
    "hsa-miR-363-3p",
    "hsa-miR-183-5p",
    "hsa-miR-451a",
    "hsa-miR-182-5p",
    "hsa-miR-191-5p",
    "hsa-miR-194-5p",
    "hsa-miR-20b-5p"
  ))

ui <- fluidPage(navbarPage(title = "draculR",
                           
                           tabPanel("Instructions",
                                    tags$h5("Welcome to draculR, a Shiny App designed to help you detect red blood cell content contamination in miR-Seq data from human plasma."),
                                    tags$h5("This App uses a new method to allocate individual samples into risk groups for haemolysis."),
                                    tags$h5(HTML(paste(
                                      "All code used to calculate data shown here is available at the following",
                                      tags$a(href="https://github.com/mxhp75/haemolysis_maternaPlasma.git", "git repository")
                                    ))),
                                    tags$h5(HTML(paste(
                                      "In case you have questions",
                                      tags$a(href="mailto::melanie.smith@adelaide.edu.au", "email me")
                                    ))),
                                    tags$br(),
                                    tags$h4("Getting started"),
                                    tags$br(),
                                    fluidRow(
                                      
                                      column(4,
                                             tags$img(src = "Picture_1.png", height = 300, width = 600)),
                                      column(4,
                                             tags$div(HTML(paste("The", tags$b("Public Data Example"), "tab allows you to click through the plot output and raw data of four datasets available on NCBI GEO that we have run through our method. Samples we consider", tags$b("Clear"), "are seen in light blue, those we consider should be used with", tags$b("Caution"), "are seen in scarlet.",
                                                           sep = " "
                                                     )))),
                                      column(4,
                                             tags$img(src = "drac.png", height = 300, width = 300))
                                      
                                    ),
                                    br(),
                                    fluidRow(
                                      
                                      column(8, 
                                      tags$div(HTML(paste("To import new data, move to the", tags$b("Import new data"), "tab. From here click", tags$b("Browse"), "to access the required raw counts table from your computer. Importing a new file will populate the main page with a set of tabs that allow you to navigate through the new information. Your count data needs to be either a comma or tab delimited file with samples in the columns and miRNA count observations in the rows. Please ensure your column of miRNA names is titled", tags$b(paste("miRNA", "name", sep = "_")), "and that the miRNA names are in a", tags$b(paste("hsa", "miR", "123", "3p", sep = "-")), "format. Samplenames need to be in a", tags$b(paste("sample", "condition", sep = "_")), "format and should not include white space or special characters.", sep = " ")))
                           )),
                                    fluidRow(
                             
                                    tags$img(src = "Picture_2.png")
                           )
                           ),
                           
                           tabPanel("Public Data Example",
                                    
                                    br(),
                                    
                                    tags$h5("Click one of these buttons to see the results"),
                                    
                                    fluidRow(
                                      column(12,
                                             actionButton(inputId = "GSE153813", label = "GSE153813"),
                                             actionButton(inputId = "GSE118038", label = "GSE118038"),
                                             actionButton(inputId = "GSE105052", label = "GSE105052"),
                                             actionButton(inputId = "GSE151341", label = "GSE151341"))
                                    ),
                                    
                                    br(),
                                    
                                    fluidRow(
                                      column = 6, offset = 0,
                                      textOutput("projectInfo")
                                    ),
                                    
                                    br(),
                                    
                                    fluidRow(column(8, offset = 0,
                                                    plotOutput(outputId = "distributionDifference") %>% 
                                                      helper(type = "inline",
                                                             icon = "exclamation",
                                                             title = "Legend Help",
                                                             colour = "red",
                                                             size = "l",
                                                             content = c("The colours used here are consistent across all tabs",
                                                                         "<b>Caution GSE_____<b> = samples from the public data example with a distribution difference > 1.9 and should be used with caution",
                                                                         "<b>Clear GSE_____<b> = samples from the public data example with a distribution difference < 1.9 and are considered clear for use",
                                                                         "<b>Haemolysed (dCq)<b> = samples used in this experiment and represent a background of samples with a dCq > 7",
                                                                         "<b>Clear (dCq)<b> = samples used in this experiment and represent a background of samples with a dCq < 7")))
                                             ),
                                    br(),
                                    
                                    fluidRow(
                                      column(8, offset = 0,
                                             DT::dataTableOutput("rawCounts"))
                                    )),
                           
                           tabPanel("Import new data",
                                    sidebarLayout(
                                      sidebarPanel(
                                        fileInput("rawDataFile","Upload the file"), # fileinput() function is used to get the file upload control option
                                        helpText("Max. file size is 5MB"),
                                        tags$hr(),
                                        h5(helpText("Select the input file parameters below")),
                                        checkboxInput(inputId = 'header',
                                                      label = 'Header?',
                                                      value = TRUE),
                                        fluidRow(
                                          column = 6,
                                          h5(helpText("Add a project title")),
                                          textInput(inputId = "project",
                                                    label = "Project",
                                                    "myProjectName"),
                                          verbatimTextOutput("value"),
                                          column = 6, offset = 6,
                                          h5(helpText("Apply your filtering value")),
                                          numericInput("filterNum", label = h5("Number in smallest group"), value = 1)) %>% 
                                          helper(icon = "exclamation",
                                                 colour = "green",
                                                 type = "markdown",
                                                 content = "filtering"),
                                        # h5(helpText("Add a project title")),
                                        # textInput(inputId = "project",
                                        #           label = "Project",
                                        #                   "myProjectName"),
                                        # verbatimTextOutput("value"),
                                        radioButtons(inputId = 'sep',
                                                     label = 'File separator',
                                                     choices = c(Comma = ',',
                                                                 Tab = '\t'),
                                                     selected = ','),
                                        h5(helpText("Select any miRNA that are differentially expressed between your groups ")) %>% 
                                          helper(icon = "exclamation",
                                                 colour = "red",
                                                 type = "markdown",
                                                 content = "drop"),
                                        checkboxGroupInput(inputId = 'drop_miRs',
                                                     label = "Drop?",
                                                     choices = c("miR-106b-3p" = 'hsa-miR-106b-3p',
                                                                 "miR-140-3p" = 'hsa-mir-140-3p',
                                                                 "miR-186-5p" = 'hsa-miR-186-5p',
                                                                 "miR-425-5p" = 'hsa-miR-425-5p',
                                                                 "miR-142-5p" = 'hsa-miR-142-5p',
                                                                 "miR-532-5p" = 'hsa-miR-532-5p',
                                                                 "miR-17-5p" = 'hsa-miR-17-5p',
                                                                 "miR-25-3p" = 'hsa-miR-25-3p',
                                                                 "miR-363-3p" = 'hsa-miR-363-3p',
                                                                 "miR-183-5p" = 'hsa-miR-183-5p',
                                                                 "miR-660-5p" = 'hsa-miR-660-5p',
                                                                 "miR-451a" = 'hsa-miR-451a',
                                                                 "miR-19b-3p" = 'hsa-miR-19b-3p',
                                                                 "miR-182-5p" = 'hsa-miR-182-5p',
                                                                 "miR-30c-5p" = 'hsa-miR-30c-5p',
                                                                 "miR-324-5p" = 'hsa-miR-324-5p',
                                                                 "miR-191-5p" = 'hsa-miR-191-5p',
                                                                 "miR-192-5p" = 'hsa-miR-192-5p',
                                                                 "miR-194-5p" = "hsa-miR-194-5p",
                                                                 "miR-20b-5p" = 'hsa-miR-20b-5p'),
                                                     textOutput("txt"))),
                                      
                                      mainPanel(
                                        uiOutput("tb")
                                      )
                                    )
                           # ),
                            
                           # tabPanel("Machine Learning",
                           #          h3("This tab is currently ",
                           #             tags$img(src = "underConstruction.png",
                           #                          heigth = 300,
                           #                          width = 300)))
                           # 
)))
  

server <- function(input, output) {
  
  # make sure help functions are active
  observe_helpers(help_dir = "helpfiles",
                  withMathJax = TRUE)
  
  # This reactive function will take the inputs from UI.R and use them for read.table() to read the data from the file. It returns the dataset in the form of a dataframe.
  # file$datapath -> gives the path of the file
  uploadData <- reactive({
    file1 <- input$rawDataFile
    if(is.null(file1)){return()} 
    read.table(file = file1$datapath,
               sep = input$sep,
               header = input$header,
               # stringsAsFactors = input$stringAsFactors)
               stringsAsFactors = FALSE)
    
  })
  
  # this reactive output will plot according to the public data reactive buttons
  plotDataPublic_miRNA <- reactiveValues(data = rankDist_GSE153813)
  
  observeEvent(input$GSE153813, { plotDataPublic_miRNA$data <- rankDist_GSE153813 })
  observeEvent(input$GSE118038, { plotDataPublic_miRNA$data <- rankDist_GSE118038 })
  observeEvent(input$GSE105052, { plotDataPublic_miRNA$data <- rankDist_GSE105052 })
  observeEvent(input$GSE151341, { plotDataPublic_miRNA$data <- rankDist_GSE151341 })
  
  output$projectInfo <- renderText({
    paste("You have selected", plotDataPublic_miRNA$data$project[1], ".",get(paste(plotDataPublic_miRNA$data$project[1],"_info", sep = "")))
  })
  
  # this reactive output will display the raw data according to the public data reactive buttons
  rawDataPublic_miRNA <- reactiveValues(data = countsRaw_GSE153813)
  
  observeEvent(input$GSE153813, { rawDataPublic_miRNA$data <- countsRaw_GSE153813 })
  observeEvent(input$GSE118038, { rawDataPublic_miRNA$data <- countsRaw_GSE118038 })
  observeEvent(input$GSE105052, { rawDataPublic_miRNA$data <- countsRaw_GSE105052 })
  observeEvent(input$GSE151341, { rawDataPublic_miRNA$data <- countsRaw_GSE151341 })
  

  # This reactive function will calculate the distribution difference and be 
  # available in all output tabs
  
  distDiff <- reactive({
    
    counts <- uploadData() %>% 
      dplyr::mutate_if(is.integer, as.numeric) %>%
      as.data.frame() %>% 
      tibble::column_to_rownames("mir_name") %>% 
      replace(is.na(.), 0)
    
    # identify samples with < 1 million reads
    lowCounts <- names(counts[, base::colSums(counts) < 1000000])
    
    # remove columns/samples with readcouns less than 1 million
    counts <- counts[, base::colSums(counts) > 1000000]
    
    # reduce any individual count less than five to zero
    counts[counts < 5] <- 0
    
    # remove miRNAs with zero counts in all samples
    counts <- counts[ base::rowSums(counts)!=0, ]
    
    # create the (super) minimal metadata table
    meta <- dplyr::data_frame(samplename = base::colnames(counts)) %>% 
      dplyr::mutate(., copy = samplename) %>% 
      tidyr::separate(., col = copy, into = c("ID", "condition"), sep = "_")
    
    # rank the samples by read counts and by unique miRs
    # this table will be joined downstream with the distribution difference table
    rank <- base::as.data.frame(base::colSums(counts)) %>%
      magrittr::set_colnames(., "readCounts") %>% 
      dplyr::arrange(., -(readCounts)) %>% 
      tibble::rownames_to_column("samplename") %>% 
      dplyr::left_join(., meta, by = "samplename") %>% 
      dplyr::select(., samplename, readCounts, condition) %>% 
      dplyr::mutate(., rank_readCounts = 1:nrow(.)) %>% 
      dplyr::full_join(.,
                       as.data.frame(t(numcolwise(nonzero)(as.data.frame(counts)))) %>%
                         tibble::rownames_to_column() %>%
                         magrittr::set_colnames(., c("samplename", "unique_miRs")) %>%
                         arrange(., desc(unique_miRs)) %>%
                         mutate(., rank_unique = 1:nrow(.)),
                       by = "samplename")
    
    
    # establish a DGEList object
    DGEList_public <- DGEList(counts = counts,
                              samples = meta)
    
    is.na(DGEList_public$counts) %>% table()
    
    # calculate normalisation factors and apply to the DGEList object
    DGEList_public <- calcNormFactors(DGEList_public, method = "TMM")
    
    # calculate the CPMs
    rawCPM <- cpm(DGEList_public, log = FALSE)
    # remove low expressed genes
    keep.exprs <- rowSums(rawCPM > 40) >= input$filterNum
    DGEList_public <- DGEList_public[keep.exprs,, keep.lib.sizes = FALSE]
    
    ## calculate the difference between the geometric mean of the distributions
    ### here we calculate the geometric mean of the classifier distribution and the
    ### "other" ensuring those taken from the classifier list are not included in
    ### other.
    # create a vector of sample names for use in the lapply
    varc <- dplyr::select(DGEList_public$samples, samplename) %>%
      tibble::remove_rownames() %>% 
      dplyr::pull(., samplename)
    
    # define the dropped classifiers as input from the groupCheckboxInput
    dropped <- subset(classifier_miRs, SYMBOL %in% input$drop_miRs)
    
    # define the final set of classifiers
    final_classifiers <- subset(classifier_miRs, SYMBOL %notin% input$drop_miRs)
    
    
    distributionDifference <- lapply(varc,function(x){
      # calculate the geometric mean of the two distribut ions (1 = classifier, 0 = other, 2 = dropped)
      dtmp <- dplyr::select(as.data.frame(edgeR::cpm(DGEList_public$counts, log = TRUE)), x) %>%
        tibble::rownames_to_column("mirna") %>% 
        mutate(., classifier = as.factor(ifelse(mirna %in% final_classifiers$SYMBOL, 1,
                                                ifelse(mirna %in% dropped$SYMBOL, 2,
                                                       ifelse(mirna %notin% classifier_miRs$SYMBOL, 0, NA)))))
      cdat_tmp <- with(dtmp,tapply(get(x),classifier,geometric.mean,na.rm=T))
      cdat <- data.frame("classifier"=rownames(cdat_tmp),"geometric.mean"=cdat_tmp)
      # calculate the difference between the two geometric means (classifier-other)  
      cdat_out <- dplyr::filter(cdat, classifier == 1)$geometric.mean - dplyr::filter(cdat, classifier == 0)$geometric.mean
      return(cdat_out)
    })
    
    names(distributionDifference) <- varc
    
    unlist_distributionDifference <- do.call(cbind.data.frame, distributionDifference) %>% 
      t() %>%
      magrittr::set_colnames("distributionDifference") %>%
      base::as.data.frame() %>% 
      tibble::rownames_to_column("samplename") %>% 
      dplyr::mutate(., haemoResult = ifelse(distributionDifference < 1.9, "Clear",
                                            ifelse(distributionDifference >= 1.9, "Caution", NA)))
    
    unlist_distributionDifference$haemoResult <- as.factor(unlist_distributionDifference$haemoResult)
    
    caution <- dim(filter(unlist_distributionDifference, haemoResult == "Caution"))
    
    dplyr::full_join(rank, unlist_distributionDifference, by = "samplename") %>%
      dplyr::mutate(., project = rep(input$project, nrow(.)))

  })
  
  
  output$rawCounts <- DT::renderDataTable({
    
    rawDataPublic_miRNA$data
    
  }, options = list(autoWidth = TRUE,
                    columnDefs = list(list(list(targets='_all',
                                                visible=TRUE,
                                                width='90') ))))
  
  output$distributionDifference <- renderPlot({
    
    caution <- dim(filter(plotDataPublic_miRNA$data, haemoResult == "Caution"))
  
    # plot as histogram side by side
    p <- ggplot() +
      geom_histogram(data = plotData_distDiff_dCq,
                     aes(x = distributionDifference,
                         fill = haemolysis,
                         colour = haemolysis,
                         y = 2*(..density..)/sum(..density..)),
                     breaks = seq(0,5,0.1),
                     alpha = 0.6, 
                     position = "identity",
                     lwd = 0.8) +
      geom_histogram(
        data = plotDataPublic_miRNA$data,
        aes(x = distributionDifference,
            fill = haemoResult,
            colour = haemoResult,
            y = 2*(..density..)/sum(..density..)),
        breaks = seq(0,5,0.1),
        alpha = 0.6, 
        position = "identity",
        lwd = 0.8) +
      scale_fill_manual(values = c("#8B0000", "#96CDCD",
                                   "#A9A9A9", "#CDCDCD"),
                        name = "Haemolysis",
                        labels = c(paste("Caution ", plotDataPublic_miRNA$data$project[1], sep = ""), paste("Clear ", plotDataPublic_miRNA$data$project[1], sep = ""),
                                   "Haemolysed (dCq)", "Clear (dCq)")) +
      scale_colour_manual(values = c("#8B0000", "#96CDCD",
                                     "#8B0000", "#96CDCD"),
                          name = "Haemolysis",
                          labels = c(paste("Caution ", plotDataPublic_miRNA$data$project[1], sep = ""), paste("Clear ", plotDataPublic_miRNA$data$project[1], sep = ""),
                                     "Haemolysed (dCq)", "Clear (dCq)")) +
      geom_vline(show.legend = FALSE,
                 xintercept = 1.9,
                 col = 2,
                 lty = 2) +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(
        x = "Distribution Difference",
        y = "% samples",
        subtitle = paste("draculR identified", caution[1], "samples to use with caution", sep = " ")
      ) +
      ggtitle(paste0(plotDataPublic_miRNA$data$project[1])) +
      theme_bw(base_size = 16) +
      theme(plot.subtitle=element_text(color="#8B0000"))
    
    p
    
    # p + annotate(
    #   geom = "text",
    #   size = 10,
    #   x = 3,
    #   y = .23,
    #   label = paste("we have identified", caution[1], "samples to use with caution", sep = " "),
    #   colour = "red"
    # )
    
  })
  
  ## Info for the Import New Data tab
  # this reactive output contains the dataset and display the file information
  output$about <- renderTable({
    if(is.null(uploadData())){return ()}
    
    # full table as uploaded by user
    input$rawDataFile
  })
  
  # this reactive output contains the summary of the dataset and display the summary in table format
  output$sum <- renderTable({
    if(is.null(uploadData())){return ()}
    
    # ordered summary table of results
    # cut the columns we don't need
    distDiff() %>% 
      dplyr::arrange(., samplename) %>% 
      dplyr::select(., Samplename = samplename,
                    `Distribution Difference` = distributionDifference,
                    `Haemolysis Result` = haemoResult,
                    Project = project)
    
  })
  
  DGEList_public <- reactive({
    
    counts <- uploadData() %>% 
      dplyr::mutate_if(is.integer, as.numeric) %>%
      as.data.frame() %>% 
      tibble::column_to_rownames("mir_name") %>% 
      replace(is.na(.), 0)
    
    # identify samples with < 1 million reads
    lowCounts <- names(counts[, base::colSums(counts) < 1000000])
    
    # remove columns/samples with readcouns less than 1 million
    counts <- counts[, base::colSums(counts) > 1000000]
    
    # reduce any individual count less than five to zero
    counts[counts < 5] <- 0
    
    # remove miRNAs with zero counts in all samples
    counts <- counts[ base::rowSums(counts)!=0, ]
    
    # create the (super) minimal metadata table
    meta <- dplyr::data_frame(samplename = base::colnames(counts)) %>% 
      dplyr::mutate(., copy = samplename) %>% 
      tidyr::separate(., col = copy, into = c("ID", "condition"), sep = "_")
    
    # rank the samples by read counts and by unique miRs
    # this table will be joined downstream with the distribution difference table
    rank <- base::as.data.frame(base::colSums(counts)) %>%
      magrittr::set_colnames(., "readCounts") %>% 
      dplyr::arrange(., -(readCounts)) %>% 
      tibble::rownames_to_column("samplename") %>% 
      dplyr::left_join(., meta, by = "samplename") %>% 
      dplyr::select(., samplename, readCounts, condition) %>% 
      dplyr::mutate(., rank_readCounts = 1:nrow(.)) %>% 
      dplyr::full_join(.,
                       as.data.frame(t(numcolwise(nonzero)(as.data.frame(counts)))) %>%
                         tibble::rownames_to_column() %>%
                         magrittr::set_colnames(., c("samplename", "unique_miRs")) %>%
                         arrange(., desc(unique_miRs)) %>%
                         mutate(., rank_unique = 1:nrow(.)),
                       by = "samplename")
    
    
    # establish a DGEList object
    DGEList_public <- DGEList(counts = counts,
                              samples = meta)
    
    is.na(DGEList_public$counts) %>% table()
    
    # calculate normalisation factors and apply to the DGEList object
    DGEList_public <- calcNormFactors(DGEList_public, method = "TMM")
    
    # calculate the CPMs
    rawCPM <- cpm(DGEList_public, log = FALSE)
    # remove low expressed genes
    keep.exprs <- rowSums(rawCPM > 40) >= input$filterNum
    DGEList_public <- DGEList_public[keep.exprs,, keep.lib.sizes = FALSE]
    
  })
  
  # This reactive output contains the dataset and display the dataset in table format
  output$table <- renderTable({
    if(is.null(uploadData())){return ()}
    uploadData() %>% 
      dplyr::mutate_if(is.integer, as.numeric) %>%
      replace(is.na(.), 0)
  })
  
  # This reactive output contains the dataset and will render the
  # mature miRNA as a function of read depth scatter plot
  output$distributions <- renderPlot({
    if(is.null(uploadData())){return ()}
    
    # define the final set of classifiers
    final_classifiers <- subset(classifier_miRs, SYMBOL %notin% input$drop_miRs)
    
    # subset the DGEList for one sample
    # temp <- as.data.frame(cpm(DGEList_public()$counts, log = TRUE)) %>%
    #   dplyr::select(., NPC0031_NPC)
    
    selectData <- reactive({
      
      as.data.frame(cpm(DGEList_public()$counts, log = TRUE)) %>%
        dplyr::select(., input$select)
      
    })
    

    # plot side by side densities
    ggplot() +
      geom_density(data = subset(selectData(), rownames(selectData()) %in% final_classifiers$SYMBOL) %>% 
                     dplyr::mutate(., colour = rep("Classifier")),
                   alpha = 0.5,
                   aes(x = !!sym(input$select),
                       fill = colour,
                       colour = colour)) +
      geom_density(data = subset(selectData(), rownames(selectData()) %notin% classifier_miRs) %>% 
                     dplyr::mutate(., colour = rep("Background")),
                   alpha = 0.5,
                   aes(x = !!sym(input$select),
                       fill = colour,
                       colour = colour)) +
      scale_fill_manual(values = c("#96CDCD", "#8B0000"),
                        name = "miRNA Set",
                        labels = c("Background", "Classifier")) +
      scale_colour_manual(values = c("#96CDCD", "#8B0000"),
                          name = "miRNA Set",
                          labels = c("Background", "Classifier")) +
      labs(title = paste(colnames(selectData()), " - ",
                         distDiff() %>%
                           dplyr::filter(., samplename == input$select) %>%
                           dplyr::select(., haemoResult) %>%
                           .[[1]],
                         sep = " "),
           x = "log2 CPM",
           y = "Density") +
      theme_bw(base_size = 16)
    

    
  })
  
  # This reactive output contains new objects and will display the histogram of 
  # distribution difference
  output$dist_diff <- renderPlot({
    if(is.null(uploadData())){return ()}
    
    # plot as histogram side by side
    q <- ggplot() +
      geom_histogram(data = plotData_distDiff_dCq,
                     aes(x = distributionDifference,
                         fill = haemolysis,
                         colour = haemolysis,
                         y = 2*(..density..)/sum(..density..)),
                     breaks = seq(0,5,0.1),
                     alpha = 0.4, 
                     position = "identity",
                     lwd = 0.8) +
      geom_histogram(
        data = distDiff(),
        aes(x = distributionDifference,
            fill = haemoResult,
            colour = haemoResult,
            y = 2*(..density..)/sum(..density..)),
        breaks = seq(0,5,0.1),
        alpha = 0.6,
        position = "identity",
        lwd = 0.8) +
      scale_fill_manual(values = c("#8B0000", "#96CDCD",
                                     "#A9A9A9", "#CDCDCD"),
                        name = "Haemolysis",
                        labels = c("Caution", "Clear",
                                   "Haemolysed (dCq)", "Clear (dCq)")) +
      scale_colour_manual(values = c("#8B0000", "#96CDCD",
                                   "#8B0000", "#96CDCD"),
                          name = "Haemolysis",
                          labels = c("Caution", "Clear",
                                     "Haemolysed (dCq)", "Clear (dCq)")) +
      geom_vline(show.legend = FALSE,
                 xintercept = 1.9,
                 col = 2,
                 lty = 2) +
      scale_y_continuous(labels = percent_format()) +
      labs(
        title = paste0("Distribution difference: ", input$project),
        subtitle = paste("draculR identified",
                         dim(filter(distDiff(),
                                    haemoResult == "Caution"))[1],
                         "samples to use with caution", sep = " "),
        x = "Distribution Difference",
        y = "% samples"
        # y = "Number of Samples"
      ) +
      theme_bw(base_size = 16) +
      theme(plot.subtitle=element_text(color="#8B0000"))
    
    q

    # print the plot to the screen
    # q + annotate(
    #   geom = "text",
    #   size = 10,
    #   x = 3,
    #   y = .23,
    #   label = paste("we have identified",
    #                 dim(filter(distDiff(),
    #                            haemoResult == "Caution"))[1],
    #                 "samples to use with caution", sep = " "),
    #   colour = "red"
    # )
    # 


  })
  
  # the following renderUI is used to dynamically generate the tabsets when the
  # file is loaded. Until the file is loaded, app will not show the tabset.
  output$tb <- renderUI({
    if(is.null(uploadData()))
      tags$h3("Test your own plasma miR-Seq data using ",
              tags$img(src = "drac.png",
                       heigth = 200,
                       width = 200))
    else
      tabsetPanel(tabPanel("About file", tableOutput("about")),
                  tabPanel("Data", tableOutput("table")),
                  tabPanel("Results Summary", tableOutput("sum")),
                  tabPanel("Distributions", 
                           br(),
                           

                           fluidRow(
                             selectInput(inputId = "select",
                                         label = "Samplename",
                                         choices = colnames(DGEList_public()$counts),
                                         selected = NULL)),
                             fluidRow(
                               plotOutput("distributions")
                             )),
                  tabPanel("Distribution Difference", plotOutput("dist_diff")))
  })
  
  
}

shinyApp(ui = ui, server = server)

