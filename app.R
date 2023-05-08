## Author: Jou-Hsuan Lee
## jhlee18@bu.edu
## BU BF591
## Final

library(shiny)
library(tidyverse)
library(DT)
library(colourpicker)
library(gplots)
library("scales")
source('functions.R')


# Define UI for application
ui <- fluidPage(
  titlePanel('BF591 Final Project'),
  'Created by Jou-Hsuan Lee', br(),
  "This is an R Shiny application designed to explore the RNA sequence analysis results of the human brain affected by Huntington's disease. The data is from ",
  a(href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4670106/', 'Labadorf et al., 2015'),
  tabsetPanel(
    # sample tab
    tabPanel('Samples',
             sidebarPanel(
               p('Explore the distinct values and distributions of sample information.'),
               fileInput('upload_sample', label='Load sample information matrix (csv or tsv file)' ,accept = c(".csv",".tsv"), placeholder = 'SraRunTable.csv'),
               actionButton('sample_button','Sumbit', width = '98%',
                            style='padding:5px; font-size:100%; color:#fff;  background-color:#326DA8'), width = 3
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel('Summary',tableOutput('sample_summary')),
                 tabPanel('Table',DTOutput('sample_table')),
                 tabPanel('Plot',
                          sidebarPanel(
                            radioButtons('x_axis','Choose a column for the histogram',
                                         choices = c('Age of Death','Age of Onset',"CAG","Duration",
                                                     "H-V Cortical Score","H-V Striatal Score","Vonsattel Grade"),
                                         selected = 'Age of Death'), width=3),
                          mainPanel(plotOutput('sample_plot')))
               )
            )
    ),
    # count tab
    tabPanel('Counts',
             sidebarPanel(
               p('Explore counts matrices to select gene count filtering strategies and evaluate their impact.'),
               fileInput('upload_counts', label='Load count data (csv or tsv file)' ,accept = c(".csv",".tsv"), placeholder = 'DESeq2_norm_counts_adjust.csv'),
               p('It might take 10-20 second to load the data and the plots.'),
               sliderInput('percent_slider','Gene filtering threshold based on the percentile(%) of variance', min = 0, max = 100, value = 20),
               sliderInput('Nsamples_slider','Genes that have at least N number of non-zero samples', min = 0, max = 69, value = 10),
               actionButton('count_button','Sumbit', width = '98%',
                            style='padding:5px; font-size:100%; color:#fff;  background-color:#326DA8'),
               width=3
             ),
             mainPanel(
                tabsetPanel(
                  tabPanel('Filtered Counts Summary',tableOutput('count_summary')),
                  tabPanel('Gene Filtering Diagnostics',plotOutput('count_diag')),
                  tabPanel('Filtered Counts Heatmap',plotOutput('count_heat')),
                  tabPanel('PCA', 
                           sidebarPanel(
                             sliderInput('x_slider','PCA x axis',min = 1, max = 30, value = 1),
                             sliderInput('y_slider','PCA y axis',min = 1, max = 30, value = 2),
                             width = 3
                           ),
                           mainPanel(plotOutput('count_pca_plot')))
               ),
             )
    ),
    tabPanel('DE',
             sidebarPanel(
               p('Explore a differential expression dataset.'),
               fileInput('upload_DEG', label='Load deg data (csv or tsv file)' ,accept = c(".csv",".tsv"), placeholder = 'DESeq2_diffexp_outlier_trimmed_adjust.csv'),
               actionButton('deg_button','Submit', width = '98%',
                            style='padding:5px; font-size:100%; color:#fff;  background-color:#326DA8'), width = 3
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel('Summary',DTOutput('deg_table')),
                 tabPanel('Volcano Plot',
                          sidebarPanel(
                            radioButtons('de_x_axis','Choose the column for the x-axis',
                                         choices = c('baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'),
                                         selected = 'log2FoldChange'),
                            radioButtons('de_y_axis','Choose the column for the y-axis',
                                         choices = c('baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'),
                                         selected = 'padj'),
                            colourInput('base','Base point color',value='darkblue'),
                            colourInput('highlight','Highlight point color',value='red'),
                            # chanege
                            sliderInput('de_slider','Select the magnitude of the -log10(p adjusted) coloring',
                                        min = -35, max = 0, value = -1.3, step = 0.1), width=4
                          ),
                          mainPanel(plotOutput('volcano')))
               )
             )
    ),
    tabPanel('GSEA',
             sidebarPanel(
               p('Explore gene set enrichment analysis results generated by fgsea.'), 
               fileInput('upload_fGSEA', label='Load fGSEA data (csv or tsv file)' ,accept = c(".csv",".tsv"), placeholder = 'fgsea_result.csv'),
               actionButton('fGSEA_button','Submit', width = '98%',
                            style='padding:5px; font-size:100%; color:#fff;  background-color:#326DA8'), width=3
            ),
             mainPanel(
               tabsetPanel(
                 tabPanel('Top Results',
                          sidebarPanel(
                            sliderInput('fgsea_slider','Numbers of Top Results by adjusted p-value',
                                        min = 1, max = 50, value = 10), width = 3
                          ),
                          mainPanel(plotOutput('fgsea_bar'))),
                 tabPanel('Table',
                          sidebarPanel(
                            sliderInput('fgsea_slider_padj','Select the magnitude of the -log10(p-adjusted value)',
                                        min = 0, max = 20, value = 1.3,step = 0.1),
                            radioButtons('fgsea_radio','Choose All/Positive/Negative NES pathways',
                                         choices = c('All','Positive','Negative'),
                                         selected = 'All'),
                            downloadButton("downloadData", "Download", style = "width:98%")
                          ),
                          mainPanel(DTOutput('fgsea_table'))),
                 tabPanel('Scatter Plots',
                          sidebarPanel(
                            sliderInput('fgsea_slider_padj_2','Select the magnitude of the -log10(p-adjusted value)',
                                        min = 0, max = 20, value = 1.3,step = 0.1)),
                          mainPanel(plotOutput('fgsea_scatter')))
              )
            )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  # maximum file size
  options(shiny.maxRequestSize=30*1024^2)
  
  # Sample tab
  # click sample_button to upload the new file.
  observeEvent(input$sample_button, {
    sample_data <- load_data(input$upload_sample) %>% summary_process()
    # summary tab
    output$sample_summary <- renderTable(summary_table(sample_data))
    # table tab
    output$sample_table <- renderDT(sample_data)
    # plot tab
    output$sample_plot <- renderPlot(summary_hist(sample_data,input$x_axis))
    })
  
  # Count tab
  # click count_button to upload the new file and get the filtered data.
  observeEvent(input$count_button,{
    count_df <- load_data(input$upload_counts) %>% count_process() %>% count_filter(input$percent_slider, input$Nsamples_slider)
    # summary
    output$count_summary <- renderTable(count_summary(count_df))
    # diagnostic scatter plots
    output$count_diag <- renderPlot(count_diagnostic(count_df))
    # heatmap
    output$count_heat <- renderPlot(count_heatmap(count_df))
    # pca
    output$count_pca_plot <- renderPlot(count_pca(count_df, input$x_slider, input$y_slider))
  })
  
  # DEG tab
  # click deg_button to upload the new file.
  observeEvent(input$deg_button, {
    deg_df <- load_data(input$upload_DEG) %>% deg_process()
    # summary
    output$deg_table <- renderDT(deg_df)
    # volcano
    output$volcano <- renderPlot(volcano_plot(deg_df, input$de_x_axis, input$de_y_axis, input$de_slider, input$base, input$highlight))
  })
  
  #fGSEA tab
  # click fGSEA_button to upload the new file.
  observeEvent(input$fGSEA_button, {
    fgsea_df <- load_data(input$upload_fGSEA)
    # top results
    output$fgsea_bar <- renderPlot(bar_pathways(fgsea_df,input$fgsea_slider))
    # table
    output$fgsea_table <- renderDT(download_fgsea(fgsea_df, input$fgsea_slider_padj, input$fgsea_radio))
    output$downloadData <- downloadHandler(
      filename = function() {
        paste0('fgsea',input$fgsea_radio, ".csv")
      },
      content = function(file) {
        write.csv(download_fgsea_table(), file)
      }
    )
    # scatter
    output$fgsea_scatter <- renderPlot(fgsea_scatter(fgsea_df, input$fgsea_slider_padj_2))
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)