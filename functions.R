library(shiny)
library(tidyverse)
library(gplots)
library("scales")
# functions for BF591 final project

#' load data
#' @param file The upload file. Should be input$file.
#' This function should be used to load every CSV and TSV file to 
#' the Shiny app as a first step 
#' before using any other processing functions on the file.
load_data <- function(file){
  req(file)
  ext <- tools::file_ext(file$name)
  data <-switch(ext,
                csv = vroom::vroom(file$datapath, delim = ","),
                tsv = vroom::vroom(file$datapath, delim = "\t"),
                validate("Invalid file; Please upload a .csv or .tsv file")
        )
  return(data)
}

#' summary_process
#' @param load_df files that have been loaded by load_data function
summary_process <- function(load_df){
  data <- load_df %>% 
    dplyr::select('Sample Name', Run, age_of_death, Diagnosis, Organism, source_name, Tissue, 
                  age_of_onset,cag,Duration,'h-v_cortical_score','h-v_striatal_score',vonsattel_grade) %>%
    dplyr::rename('Age of Death'=age_of_death, 'Source Name'=source_name, 'Age of Onset'=age_of_onset, 
                  'CAG'=cag, 'H-V Cortical Score'='h-v_cortical_score',
                  'H-V Striatal Score'='h-v_striatal_score','Vonsattel Grade'=vonsattel_grade )
  # col name should convert to factor
  names <- colnames(data)[sapply(data, is.character)][-1:-2]
  data[names] <- lapply(data[names], as_factor)
  return(data)
}


#' summary_table
#' @param process_df The data frame that has been processed.
#'  (by summary_process())
summary_table <- function(process_df){
  names <- colnames(process_df)[4:7]
  process_df[names] <- lapply(process_df[names], as_factor)
  H_start <- which(process_df$Diagnosis=="Huntington's Disease")[1] # row index of first Huntington's Disease
  df <- tibble(
    'Column Name' = colnames(process_df)[-1:-2],
    'Type' = sapply(process_df[-1:-2], class),
    'Mean (sd) or Distinct Values' = sapply(process_df[-1:-2], function(x){
      if (!is.na(x) && is.numeric(x)){paste0(round(mean(x),2),'(+/-',round(sd(x),2),')')}
      else if (is.factor(x)){paste0(levels(x),collapse=', ')}
      else{paste0(round(mean(x[H_start:length(x)]),2),' (+/-',round(sd(x[H_start:length(x)]),2),')'," in Huntington's Disease Samples only")}
    }) 
  )
  return(df)
}

#' summary_hist
#' @param process_df The dataframe that has been processed.
#' (by summary_process())
#' @param column The column to draw the histogram
summary_hist <- function(process_df, column){
  p <- ggplot(process_df, mapping = aes(x=!!sym(column))) + 
    theme_classic() + ggtitle(paste('Distrubution of',column))
  if(!!sym(column)!='Age of Death'){
    p <- p + labs(caption = "Huntington's Disease Samples only")
  }
  
  if (column == 'Vonsattel Grade'){
      p <- p + scale_x_continuous(limits = c(0, 5)) + 
      geom_bar(fill="lightblue",color="#e9ecef", alpha=0.9, stat = 'count')
  } else if (column == 'H-V Cortical Score' || column == 'H-V Striatal Score'){
    p <- p + xlim(0, 5) + 
      geom_bar(fill="lightblue",color="#e9ecef", alpha=0.9, stat = 'bin')
  } else if (column != 'Age of Death' && column != 'Age of Onset'){
    p <- p + scale_x_continuous(breaks = pretty_breaks()) + 
      geom_bar(fill="lightblue",color="#e9ecef", alpha=0.9, stat = 'count')
  } else {
    p <- p + geom_bar(fill="lightblue",color="#e9ecef", alpha=0.9, stat = 'bin')
  }
  return(p)
}

#' count_process
#' @param load_df files that have been loaded by load_data function
count_process <- function(load_df){
  df <- load_df %>% rename(Gene='...1')
  return(df)
}

#' count_filter
#' @param process_df The data frame that has been processed.
#'  (by count_process())
#' @param percent_slider filtering threshold based on the percentile(%) of variance. 
#' should be input$percent_slider
#' @param Nsample_slider Genes that have at least N number of non-zero samples.  
#' should be input$Nsamples_slider
#' Calculates the variance and number of non-zero samples for each gene and 
#' assigns a pass or no pass label based on whether the gene meets the filtering criteria
count_filter <- function(process_df,percent_slider,Nsample_slider){
  filter_df <- process_df %>% mutate(
    variance = apply(.[, -1], 1, var),
    non_zero = apply(.[, -1], 1, function(x) sum(x != 0)),
    pass = as.factor(case_when(variance >= quantile(variance, percent_slider/100, na.rm = TRUE) &
                      non_zero >= Nsample_slider ~ 'pass',
                     TRUE ~ 'no')))
  return(filter_df)
}

#' count_summary
#' @param filter_df data frame that has been processed by count_filter() 
#' summarizes the number of samples and genes and the number and percentage 
#' of genes that pass a certain filtering threshold
count_summary <- function(filter_df){
  dim_df <- filter_df %>% dplyr::select(!c('Gene','variance','non_zero','pass')) %>% dim()
  dim_cf <- filter_df[filter_df$pass=='pass',] %>% dplyr::select(!c('Gene','variance','non_zero','pass')) %>% dim()
  summary_df <- tibble(
      'Number of Samples' = as.integer(dim_df[2]),
      'Number of Genes' = dim_df[1],
      'Passed Genes' = paste0(dim_cf[1], '(', round(dim_cf[1]/dim_df[1]*100,2),'%)'),
      'Not Passed Genes' = paste0(dim_df[1]-dim_cf[1],'(',round((dim_df[1]-dim_cf[1])/dim_df[1]*100,2),'%)')
    )
  return(summary_df)
}

#' count_diagnostic
#' @param filter_df data frame that has been processed by count_filter() 
#' Creates two diagnostic plots. 
#' The first plot shows the relationship between the median count and variance for each gene, 
#' colored by whether it passed the filtering criteria or not. 
#' The second plot shows the relationship between the median count and the number of samples 
#' with zero counts for each gene, also colored by whether it passed the filtering criteria or not.
count_diagnostic <-function(filter_df){
  # calculate the median
  c_plot <- filter_df %>% mutate(Median=apply(.[-c(1,ncol(filter_df)-2,ncol(filter_df)-1, ncol(filter_df))], 1, median)) 
  # variance plot
  variance <- ggplot(c_plot) + geom_point(mapping = aes(x=log2(Median), y=log2(variance), color=pass)) +
    scale_color_manual(values = c("no" = "lightblue", "pass" = "#0B346E")) +
    theme_classic() + theme(legend.position = 'none') +
    ggtitle("Median count vs variance")
  # number of samples of with zero counts plot
  zero <- ggplot(c_plot) + geom_point(mapping = aes(x=log2(Median), y=69-non_zero, color=pass)) +
    scale_color_manual(values = c("no" = "lightblue", "pass" = "#0B346E")) +
    theme_classic() + ylab('Number of zeros')
    ggtitle("Median count vs zero")
  
  p <- gridExtra::grid.arrange(variance, zero, ncol=2)
  return(p)
}

#' to_matrix
#' @param filter_df data frame that has been processed by count_filter() 
#' Turn the filtered dataframe into matrix for ploting heatmap
#' It was used in count_heatmap and count_pca functions
to_matrix <- function(filter_df){
  count_matrix <- filter_df[filter_df$pass=='pass',] %>% select(-c('Gene','variance','non_zero','pass')) %>% as.matrix()
  rownames(count_matrix) <- filter_df$Gene[filter_df$pass=='pass']
  return(count_matrix)
}

#' count_heatmap
#' @param filter_df data frame that has been processed by count_filter() 
#' Turn the filter_df into heatmap
count_heatmap <- function(filter_df){
  count_matrix <- to_matrix(filter_df)
  log_matrix <- log2(count_matrix+1)
  p <- gplots::heatmap.2(log_matrix, labRow=FALSE, trace="none", main='Heatmap' ,ylab="Genes", xlab="Samples")
  return(p)
}

#' count_pca
#' @param filter_df data frame that has been processed by count_filter() 
#' @param x_axis The pc on x axis. should be input$x_slider
#' @param y_axis The pc on y axis.should be input$y_slider
#' Performs principal component analysis (PCA) on the count matrix, 
#' and creates a scatter plot of the PCA results, 
#' colored by diagnosis (Huntington's Disease or neurologically normal). 
#' The user can specify which principal components to plot on the x and y axes. 
count_pca <- function(filter_df, x_axis, y_axis){
  count_matrix <- to_matrix(filter_df)
  pca_results <- prcomp(scale(t(count_matrix)), center=FALSE, scale=FALSE)
  pca_df <- cbind(pca_results$x[,c(x_axis,y_axis)],rownames(pca_results$x)) %>% as_tibble() %>%
    mutate(Diagnosis=case_when(
      str_detect(V3, 'H_') ~ "Huntington's Disease", TRUE ~ 'Neurologically normal'))
  # % of each PC
  all_var <- pca_results$sdev^2
  pca_var <- all_var/sum(all_var)
  
  p <- ggplot(pca_df) + geom_point(mapping = aes(x=pca_results$x[,x_axis], y=pca_results$x[,y_axis], color=Diagnosis)) +
    theme_classic() + xlab(paste0('PC',x_axis,' (', round(pca_var[x_axis]*100,2),'%)'))+
    ylab(paste0('PC',y_axis,' (', round(pca_var[y_axis]*100,2),'%)')) + ggtitle('Principal Component Analysis') + 
    scale_color_manual(values=c("red", "darkblue"))
  return(p)
}

#' deg_process
#' @param load_df files that have been loaded by load_data function
deg_process <- function(load_df){
  df <- load_df %>% rename( 'Ensemble ID' ='...1') %>% select(-c('HD.mean','Control.mean'))
  return(df)
}


#' @param dataf The data frame that has been processed by deg_process()
#' @param x_name The column name to plot on the x-axis
#' @param y_name The column name to plot on the y-axis
#' @param slider A negative integer value representing the magnitude of
#' p-adjusted values to color. Most of our data will be between -1 and -300.
#' @param color1 One of the colors for the points.
#' @param color2 The other colors for the points.
#'
#' @return A ggplot object of a volcano plot
#' Generate a volcano plot
#' !!sym() may be required to access column names in ggplot aes().
volcano_plot <- function(deg_df, x_name, y_name, slider, color1, color2) {
  g <- ggplot(deg_df) +
    geom_point(mapping = aes(x=!!sym(x_name),y=-log10(!!sym(y_name)), color=deg_df$padj<10^(slider))) +
    scale_color_manual(paste0('padj<10^(',slider,')'),values = c('FALSE'=color1,'TRUE'=color2)) +
    theme_bw() + theme(legend.position = 'bottom') + ggtitle('Volcano Plot of choosen columns')
  return(g)
}

#' bar_pathways
#' @param fgsea_df the fgsea results in tibble format
#' @param num_paths the number of pathways. should be input$fgsea_slider
bar_pathways <- function(fgsea_df, num_paths){
  
  top <- fgsea_df %>% arrange(padj) %>% head(num_paths) %>%
    mutate(pathway = str_replace_all(pathway, '_', ' ')) %>%
    mutate(pathway = str_replace(pathway, ' ',': ')) %>%
    mutate(pathway = forcats::fct_reorder(pathway, NES)) 
  # change to pathway name's size accord to the number of pathways
  word_size <- base::ifelse(num_paths>15,5,8)
  
  p <- ggplot(top) + 
    geom_bar(aes(x=pathway, y=NES, fill = NES > 0), stat='identity') +
    scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
    theme_minimal() +
    ggtitle(paste0('Top ', num_paths,' pathways of fGSEA results')) +
    ylab('NES') +
    xlab('') + scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
    theme(axis.text.y = element_text(size=word_size)) +
    coord_flip() 
  return(p)
}

#' download_fgsea
#' @param fgsea_df the fgsea results in tibble format
#' @param padj_num The Slider to filter table by adjusted p-value
#' should be input$fgsea_slider_padj
#' @param pos_neg  select all, positive or negative NES pathways
#' should be  input$fgsea_radio
download_fgsea <- function(fgsea_df, padj_num, pos_neg){
  padj_num <- 10^(-padj_num)
  df <- fgsea_df %>% dplyr::filter(padj<padj_num)
  df$pval <- formatC(df$pval, format = 'e')
  df$padj <- formatC(df$padj, format = 'e')
  df$log2err <- formatC(df$log2err, digits = 3)
  df$NES <- formatC(df$NES, digits = 3)
  df$ES <- formatC(df$ES, digits = 3)
  if(pos_neg == 'All'){
    return(df)
  }else if(pos_neg=='Positive'){
    return(df[df$NES > 0,])
  }else{
    return(df[df$NES < 0,])
  }
}

#' fgsea_scatter
#' @param load_df the fgsea results in tibble format
#' @param padj_num The Slider to filter table by adjusted p-value
#' should be fgsea_slider_padj_2
#' Scatter plot of NES on x-axis and -log10 adjusted p-value on
#' y-axis, with gene sets below threshold in grey color
fgsea_scatter <- function(load_df, padj_num){
  p <- ggplot(load_df) + 
    geom_point(aes(x=NES, y=-log10(padj), color=-log10(padj)>padj_num)) +
    scale_color_manual(values = c('TRUE' = 'darkblue','FALSE' = 'grey')) + 
    theme_minimal() + labs(color=paste0('-log10(padj) > ',padj_num)) +
    ggtitle('NES vs. -log10 adjusted p-value')
  return(p)
}