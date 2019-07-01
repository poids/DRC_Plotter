#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
# br()=linebreak

#For Pipeline
library(emdbook)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(utils)
library(Hmisc)
library(tidyverse)
library(dplyr)
library(lubridate)
library(readxl)

#For WebApp
library(shiny)
library(shinythemes)
library(DT)
# library(shinyWidgets)

# Define UI for application that draws a histogram
shinyUI(fluidPage(theme = shinytheme("cerulean"),
  
  # Application title
  # titlePanel("Dose Response Curve"),
  
  #Tab Architecture
  navbarPage("Dose Response Curve", id="tabs",

    #DRC Plot Tab
    tabPanel("Summary Statistics", fluid = TRUE,
    
     sidebarLayout(
       sidebarPanel(
         img(src='codabioOrig.png',
             align = "left",
             width='100%'),
         
         tags$hr(),
         
         fileInput("unprocesssed_sumstat_df", "Choose CSV File",
                   accept = c(
                     "text/csv",
                     "text/comma-separated-values,text/plain",
                     ".csv")
         ),
         
         # Horizontal line ----
         tags$hr(),
         h4("Subset by exp date and then check the box to calculate stats:", align='left'),
         
         
         dateRangeInput("exp_range",
                        label = "Exp Date Range",
                        start=NULL,
                        end=NULL,
                        format='yyyymmdd'
         ),
         
         #Checkbox that will run script after dataframe has been subset by exp date range
         checkboxInput('run', 
                       label='Subset and Calculate Stats',
                       value=F),
         
         tags$hr(),
         
         selectInput("ligand",
                     "Subset by Ligand",
                     ""),
         
         selectizeInput("key",
                        "Subset by Ligand-Key",
                        "",
                        multiple=T,
                        options = list(
                          placeholder = 'Select or Type'
                        )
         ),
         
          selectizeInput("construct",
                        "Select Constructs to Show",
                        multiple=T,
                        choices="",
                        options = list(
                          placeholder = 'Select or Type'
                        )
         ),
         
         checkboxInput('normalize',
                       label = 'Normalize to Top',
                       value = TRUE),
         
         
         
         # Horizontal line ----
         tags$hr(),
         
         
         #All the following only are selectable if df_type='auto' [useuiOutput instead?]
         h4("Only for Sumstat_master_df generated in pipeline:", align='left'),
         
         
         radioButtons('plate_type',
                      label = 'Select Plate Type',
                      choices = list(
                        'Both'=1,
                        'Ensemble'=2,
                        'Single'=3
                      ),
                      selected=1),
         
         
         sliderInput("success_rate",
                     "Percent Success Rate Threshold:",
                     min = 0,
                     max = 100,
                     value = 0)
       ),
       
       
       # Show a plot of the generated distribution
       mainPanel(
         
         DTOutput('final_df'),
         
         downloadButton("save_df", "Save Dataframe", icon("floppy-o"))
         
       )
     )       
     
    ),
    
    
    tabPanel("Dose Response Plot", fluid=TRUE,
        sidebarLayout(
          sidebarPanel(
            img(src='codabioOrig.png',
                align = "left",
                width='100%'),
            
            tags$hr(),
            
            #Show Error Bars
            checkboxInput('geombars',
                          label = 'Show Error Bars',
                          value = TRUE),

            #Save Plot
            downloadButton("save_plot_png", "Save Plot Img", icon("floppy-o"))
            
          ),
          mainPanel(
            plotOutput('drc'),
            
            tags$hr(),
            br(),
            br(),
            tags$hr(),
            
            tableOutput('pop_n')
            # DTOutput('pop_n')
            
          )
        )
    ),
    
    tabPanel("Log10 FC of EC50", fluid=TRUE,
             sidebarLayout(
               sidebarPanel(
                 img(src='codabioOrig.png',
                     align = "left",
                     width='100%'),
                 
                 tags$hr(),
                 
                 selectInput("constructFC",
                             "Select Construct to Normalize to:",
                             ""),
                 
                 #Save Plot
                 downloadButton("save_FCplot", "Save Plot", icon("floppy-o"))
                 
               ),
               mainPanel(
                 plotOutput('fc_plot')
               )
             )
    )

  )
))

#Written by Vasco Morais (March 2019) Cell: 415-845-2118