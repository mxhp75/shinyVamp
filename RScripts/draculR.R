library(dplyr)
library(plyr)
library(ggplot2)
library(shiny)
library(ggrepel)
library(scales) # For percent_format()
library(tidyr)
library(readr)

# import our distribution difference dataframe
plotData_distDiff_dCq <- read_csv("/Users/a1627211/Bioinformatics/github/PhD/shinyVamp/www/plotData_distDiff_dCq.csv")

# import the rank and distribtuion difference data for GSE118038
rank_GSE118038 <- read_csv(file = "/Users/a1627211/Bioinformatics/github/PhD/shinyVamp/www/rank_GSE118038.csv")
unlistDistributionDifference_GSE118038 <- read_csv("/Users/a1627211/Bioinformatics/github/PhD/shinyVamp/www/unlist_distributionDifference_GSE118038.csv")

rankDist_GSE118038 <- dplyr::full_join(rank_GSE118038, unlistDistributionDifference_GSE118038, by = "samplename") %>% 
  dplyr::mutate(., project = rep("GSE118038", nrow(.)))

# import the rank and distribtuion difference data for GSE105052
rank_GSE105052 <- read_csv("/Users/a1627211/Bioinformatics/github/PhD/shinyVamp/www/rank_GSE105052.csv")
unlistDistributionDifference_GSE105052 <- read_csv("/Users/a1627211/Bioinformatics/github/PhD/shinyVamp/www/unlist_distributionDifference_GSE105052.csv")

rankDist_GSE105052 <- dplyr::full_join(rank_GSE105052, unlistDistributionDifference_GSE105052, by = "samplename") %>% 
  dplyr::mutate(., project = rep("GSE105052", nrow(.)))

# import the rank and distribtuion difference data for GSE151341
rank_GSE151341 <- read_csv("/Users/a1627211/Bioinformatics/github/PhD/shinyVamp/www/rank_GSE151341.csv")
unlistDistributionDifference_GSE151341 <- read_csv("/Users/a1627211/Bioinformatics/github/PhD/shinyVamp/www/unlist_distributionDifference_GSE151341.csv")

rankDist_GSE151341 <- dplyr::full_join(rank_GSE151341, unlistDistributionDifference_GSE151341, by = "samplename") %>% 
  dplyr::mutate(., project = rep("GSE151341", nrow(.)))

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
                                    fluidRow(
                                      column(12,
                                             actionButton(inputId = "GSE118038", label = "GSE118038"),
                                             actionButton(inputId = "GSE105052", label = "GSE105052"),
                                             actionButton(inputId = "GSE151341", label = "GSE151341"))
                                    ),
                                    br(),
                                    fluidRow(
                                      column(6, offset = 0,
                                             plotOutput("mature_miRNA")),
                                      column(6, offset = 6,
                                             plotOutput("distributionDifference"))
                                    )
                                    
                           ),
                           
                           tabPanel("Import your data",
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
                                        br(),
                                        radioButtons(inputId = 'sep',
                                                     label = 'Separator',
                                                     choices = c(Comma = ',',
                                                                 Semicolon = ';',
                                                                 Tab = '\t',
                                                                 Space = ''),
                                                     selected = ',')
                                      ),
                                      mainPanel(
                                        uiOutput("tb")
                                      )
                                    )
                           ),
                            
                           tabPanel("Machine Learning",
                                    verbatimTextOutput("prediction"))
))
  

server <- function(input, output) {
  # This reactive function will take the inputs from UI.R and use them for read.table() to read the data from the file. It returns the dataset in the form of a dataframe.
  # file$datapath -> gives the path of the file
  data <- reactive({
    file1 <- input$rawDataFile
    if(is.null(file1)){return()} 
    read.table(file = file1$datapath,
               sep = input$sep,
               header = input$header,
               stringsAsFactors = input$stringAsFactors)
    
  })
  
  # this reactive output will plot according to the public data reactive buttons
  plotDataPublic_miRNA <- reactiveValues(data = NULL)
  
  observeEvent(input$GSE118038, { plotDataPublic_miRNA$data <- rankDist_GSE118038 })
  observeEvent(input$GSE105052, { plotDataPublic_miRNA$data <- rankDist_GSE105052 })
  observeEvent(input$GSE151341, { plotDataPublic_miRNA$data <- rankDist_GSE151341 })
  
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
      ggtitle(paste(plotDataPublic_miRNA$data$project[1], ": Mature miRNA as a \n function of Read Counts", sep = "")) +
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
      scale_y_continuous(labels = percent_format()) +
      labs(
        title = "Distribution difference using final classifiers",
        subtitle = "based on three classification groups",
        x = "Distribution Difference",
        y = "% samples"
      ) +
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
  
  # this reactive output contains the summary of the dataset and display the summary in table format
  output$filedf <- renderTable({
    if(is.null(data())){return ()}
    input$rawDataFile
  })
  
  # this reactive output contains the summary of the dataset and display the summary in table format
  output$sum <- renderTable({
    if(is.null(data())){return ()}
    summary(data())
    
  })
  
  # This reactive output contains the dataset and display the dataset in table format
  output$table <- renderTable({
    if(is.null(data())){return ()}
    data()
  })
  
  # the following renderUI is used to dynamically generate the tabsets when the file is loaded. Until the file is loaded, app will not show the tabset.
  output$tb <- renderUI({
    if(is.null(data()))
      h5("You are using", tags$img(src = "www/teeth.png",
                                   heigth = 200,
                                   width = 200))
    else
      tabsetPanel(tabPanel("About file", tableOutput("filedf")),
                  tabPanel("Data", tableOutput("table")),
                  tabPanel("Summary", tableOutput("sum")))
  })
  
   
  
}

shinyApp(ui = ui, server = server)