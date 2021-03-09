library(dplyr)
library(ggplot2)
library(shiny)
library(ggrepel)
library(tidyr) 

ui <- fluidPage(navbarPage(title = "draculR",
                           
                           tabPanel("Instructions",
                                    h5("This app uses two methods to classify miR-Seq sequencing data from human plasma into risk groups for haemolysis contamination"),
                                    h5(HTML(paste(
                                      "All code used to calculate data shown here is available at the following",
                                      tags$a(href="https://github.com/mxhp75/haemolysis_maternaPlasma.git", "git repository")
                                    ))),
                                    h5(HTML(paste(
                                      "In case you have questions",
                                      tags$a(href="mailto::melanie.smith@adelaide.edu.au", "email me")
                                    )))
                           ),
                           
                           tabPanel("Public Data Example",
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
      h5("You are using", tags$img(src='shinyVamp/inputData/teeth.png', heigth=200, width=200))
    else
      tabsetPanel(tabPanel("About file", tableOutput("filedf")),
                  tabPanel("Data", tableOutput("table")),
                  tabPanel("Summary", tableOutput("sum")))
  })
  
   
  
}

shinyApp(ui = ui, server = server)