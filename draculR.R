library(dplyr)
library(plyr)
library(ggplot2)
library(shiny)
library(ggrepel)
library(scales) # For percent_format()
library(tidyr)
library(magrittr)
library(edgeR)
library(readr)
library(psych)

# User defined functions
# negate %in%
`%notin%` <- Negate(`%in%`)
# function (x) to count the number of non-zero records in each column (ie per sample)
nonzero <- function(x) sum(x != 0)

# import our distribution difference dataframe
plotData_distDiff_dCq <- read_csv("www/plotData_distDiff_dCq.csv")

# import the rank and distribtuion difference data for GSE153813
rank_GSE153813 <- read_csv(file = "www/rank_GSE153813.csv")
unlistDistributionDifference_GSE153813 <- read_csv(file = "www/unlist_distributionDifference_GSE153813.csv")
GSE153813_info <- paste("The data containined in GSE153813 were obtained from",
                        " NCBI GEO. There is no associated publication" , sep = "")

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

rankDist_GSE118038 <- dplyr::full_join(rank_GSE118038, unlistDistributionDifference_GSE118038, by = "samplename") %>% 
  dplyr::mutate(., project = rep("GSE118038", nrow(.)))

# import the rank and distribtuion difference data for GSE105052
rank_GSE105052 <- read_csv(file = "www/rank_GSE105052.csv")
unlistDistributionDifference_GSE105052 <- read_csv(file = "www/unlist_distributionDifference_GSE105052.csv")
GSE105052_info <- paste("The data containined in GSE105052 were obtained from",
                        " NCBI GEO and originate from the article \"Small ",
                        " RNA-seq analysis of circulating miRNAs to identify phenotypic variability in",
                        " Friedreich's ataxia patients.\" Sci Data 2018 Mar 6;5:180021.", sep = "")

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
                                    tags$h5("This app uses two methods to classify miR-Seq sequencing data from human plasma into risk groups for haemolysis contamination"),
                                    tags$h5(HTML(paste(
                                      "All code used to calculate data shown here is available at the following",
                                      tags$a(href="https://github.com/mxhp75/haemolysis_maternaPlasma.git", "git repository")
                                    ))),
                                    tags$h5(HTML(paste(
                                      "In case you have questions",
                                      tags$a(href="mailto::melanie.smith@adelaide.edu.au", "email me")
                                    )))
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
                                      column = 8, offset = 0,
                                      textOutput("projectInfo")
                                    ),
                                    br(),
                                    fluidRow(
                                      column(6, offset = 0,
                                             plotOutput("mature_miRNA")),
                                      column(6, offset = 6,
                                             plotOutput("distributionDifference"))
                                    )
                                    
                           ),
                           
                           tabPanel("Import new data",
                                    sidebarLayout(
                                      sidebarPanel(
                                        fileInput("rawDataFile","Upload the file"), # fileinput() function is used to get the file upload control option
                                        helpText("Default max. file size is 5MB"),
                                        tags$hr(),
                                        h5(helpText("Select the read.table parameters below")),
                                        checkboxInput(inputId = 'header',
                                                      label = 'Header?',
                                                      value = TRUE),
                                        checkboxInput(inputId = "stringAsFactors",
                                                      label = "stringAsFactors",
                                                      value = FALSE),
                                        h5(helpText("Add a project title")),
                                        textInput(inputId = "project",
                                                  label = "Project",
                                                          "myProjectName"),
                                        verbatimTextOutput("value"),
                                        radioButtons(inputId = 'sep',
                                                     label = 'Separator',
                                                     choices = c(Comma = ',',
                                                                 Tab = '\t'),
                                                     selected = ','),
                                        h5(helpText("Select any miRNA that are DE between your groups")),
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
                                                     textOutput("txt"))
                                      ),
                                      mainPanel(
                                        uiOutput("tb")
                                      )
                                    )
                           ),
                            
                           tabPanel("Machine Learning",
                                    h3("This tab is currently ",
                                       tags$img(src = "underConstruction.png",
                                                    heigth = 300,
                                                    width = 300)))
                           
))
  

server <- function(input, output) {
  
  # This reactive function will take the inputs from UI.R and use them for read.table() to read the data from the file. It returns the dataset in the form of a dataframe.
  # file$datapath -> gives the path of the file
  uploadData <- reactive({
    file1 <- input$rawDataFile
    if(is.null(file1)){return()} 
    read.table(file = file1$datapath,
               sep = input$sep,
               header = input$header,
               stringsAsFactors = input$stringAsFactors)
    
  })
  
  # this reactive output will plot according to the public data reactive buttons
  plotDataPublic_miRNA <- reactiveValues(data = rankDist_GSE153813)
  
  observeEvent(input$GSE153813, { plotDataPublic_miRNA$data <- rankDist_GSE153813 })
  observeEvent(input$GSE118038, { plotDataPublic_miRNA$data <- rankDist_GSE118038 })
  observeEvent(input$GSE105052, { plotDataPublic_miRNA$data <- rankDist_GSE105052 })
  observeEvent(input$GSE151341, { plotDataPublic_miRNA$data <- rankDist_GSE151341 })
  
  output$projectInfo <- renderText({
    paste("You have selected", plotDataPublic_miRNA$data$project[1], ".",get(paste(plotDataPublic_miRNA$data$project[1],"_info", sep = "")))
  }
  )
  
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
    keep.exprs <- rowSums(rawCPM > 40) >= 12
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
      # calculate the geometric mean of the two distributions (1 = classifier, 0 = other, 2 = dropped)
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
    
    # rankDist <- dplyr::full_join(rank, unlist_distributionDifference, by = "samplename") %>%
    #   dplyr::mutate(., project = rep(input$project, nrow(.)))
    
  })
  
  output$mature_miRNA <- renderPlot({
    
    ggplot(data = plotDataPublic_miRNA$data,
           aes(x = readCounts,
               y = unique_miRs)) +
      geom_point(aes(colour = haemoResult,
                     size = 4)) +
      scale_size(guide = "none") +
      scale_x_continuous(name = "Filtered Read Counts",
                         breaks = seq(round_any(min(plotDataPublic_miRNA$data$readCounts), 10, f = floor), max(plotDataPublic_miRNA$data$readCounts), 100000)) +
      scale_y_continuous(name = "Mature miRNA Identified",
                         breaks = seq(round_any(min(plotDataPublic_miRNA$data$unique_miRs), 10, f = floor), max(plotDataPublic_miRNA$data$unique_miRs), 50)) +
      stat_smooth(method = 'loess',
                  se = FALSE,
                  size = 2) +
      geom_text_repel(data = filter(plotDataPublic_miRNA$data, haemoResult == "Caution"),
                      box.padding = 1,
                      aes(label = samplename),
                      show.legend = FALSE) +
      ggtitle(paste0(plotDataPublic_miRNA$data$project[1])) +
      theme_bw(base_size = 16) +
      theme(
        axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm"))
      ) +
      theme(
        axis.text.x = element_text(angle = 270, hjust = 1)
        ) 
    
  })
  
  output$distributionDifference <- renderPlot({
  
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
            fill = project,
            colour = project,
            y = 2*(..density..)/sum(..density..)),
        breaks = seq(0,5,0.1),
        alpha = 0.6, 
        position = "identity",
        lwd = 0.8) +
      geom_vline(show.legend = FALSE,
                 xintercept = 1.9,
                 col = 2,
                 lty = 2) +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(
        x = "Distribution Difference",
        y = "% samples"
      ) +
      ggtitle(paste0(plotDataPublic_miRNA$data$project[1])) +
      theme_bw(base_size = 16)
    
    caution <- dim(filter(plotDataPublic_miRNA$data, haemoResult == "Caution"))
    
    p + annotate(
      geom = "text",
      x = 3,
      y = .23,
      label = paste("we have identified", caution[1], "samples to use with caution", sep = " "),
      colour = "red"
    )
    
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
    distDiff() %>% 
      dplyr::arrange(., samplename)
    
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
      scale_fill_manual(values = c("#999999", "#E69F00",
                                     "#56B4E9", "#009E73")) +
      scale_colour_manual(values = c("#999999", "#E69F00",
                                   "#56B4E9", "#009E73")) +
      geom_vline(show.legend = FALSE,
                 xintercept = 1.9,
                 col = 2,
                 lty = 2) +
      scale_y_continuous(labels = percent_format()) +
      labs(
        title = paste0("Distribution difference: ", input$project),
        x = "Distribution Difference",
        y = "% samples"
      ) +
      theme_bw(base_size = 16)

    # print the plot to the screen
    q + annotate(
      geom = "text",
      x = 3,
      y = .23,
      label = paste("we have identified", dim(filter(distDiff(), haemoResult == "Caution"))[1], "samples to use with caution", sep = " "),
      colour = "red"
    )
    


  })
  
  # the following renderUI is used to dynamically generate the tabsets when the
  # file is loaded. Until the file is loaded, app will not show the tabset.
  output$tb <- renderUI({
    if(is.null(uploadData()))
      tags$h3("Test your own plasma miR-Seq data using ", tags$img(src = "drac.png",
                                                      heigth = 200,
                                                      width = 200))
    else
      tabsetPanel(tabPanel("About file", tableOutput("about")),
                  tabPanel("Data", tableOutput("table")),
                  tabPanel("Results Summary", tableOutput("sum")),
                  tabPanel("Distributions", plotOutput("distributions")),
                  tabPanel("Distribution Difference", plotOutput("dist_diff")))
  })
  
  
}

shinyApp(ui = ui, server = server)

