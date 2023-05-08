library(fgsea)
library(tidyverse)

deg <- read_csv('BF591-Final/data/DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.csv') %>% dplyr::rename(Gene='...1')


deg_df <- deg %>% arrange(desc(log2FoldChange)) %>% dplyr::filter(padj<0.05)
rnk_list <- deg_df$log2FoldChange
names(rnk_list) <- deg_df$symbol

# fgsea
pathways_fgsea <- fgsea::gmtPathways('/usr4/bf527/jhlee18/R/final/BF591-Final/data/c2.cp.v2023.1.Hs.symbols.gmt')
fgsea_results <- fgsea(pathways_fgsea, rnk_list, 15, 500)
fgsea_results <- fgsea_results %>% as_tibble() %>% select(-leadingEdge)

write_csv(fgsea_results,'fgsea_result.csv')