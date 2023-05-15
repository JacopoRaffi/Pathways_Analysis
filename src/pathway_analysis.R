source("graph_lib.R")
library(biogridr)
library(dplyr)

bg_get_key("Jacopo", "Raffi", "j.raffi@studenti.unipi.it", "metabolic_disease")

drugs = read.csv("../drugs/pharmacologically_active.csv")

gene_names = drugs$Gene.Name
gene_names = gene_names[gene_names != "" & !is.na(gene_names)]

biogrid_df = as.data.frame(bg("interactions") %>% 
                              bg_constrain(geneList = gene_names[1]) %>% 
                              bg_get_results())

biogrid_df = as.data.frame(biogrid_df[gene_names[1] %in% biogrid_df$official_symbol_for_interactor_a || gene_names[1] %in% biogrid_df$official_symbol_for_interactor_b])

for(i in 2:length(gene_names)){
  tmp_df = as.data.frame(bg("interactions") %>% 
                               bg_constrain(geneList = gene_names[i]) %>% 
                               bg_get_results())
  
  tmp_df = as.data.frame(tmp_df[gene_names[i] %in% tmp_df$official_symbol_for_interactor_a || gene_names[i] %in% tmp_df$official_symbol_for_interactor_b])
  
  biogrid_df = rbind(biogrid_df, tmp_df)
}

ig_drugs = csv2igraph("../drugs/pharmacologically_active.csv", "../drugs/drugbank vocabulary.csv")

ig_diseases = folder2igraph("../data")
ig_diseases = splitFolderNodes(ig_diseases, "../data")