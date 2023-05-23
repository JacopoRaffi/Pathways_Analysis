source("graph_lib.R")
library(biogridr)
library(dplyr)

bg_get_key("Jacopo", "Raffi", "j.raffi@studenti.unipi.it", "metabolic_disease")

drugs = read.csv("../drugs/pharmacologically_active.csv")

#bgrd_data = read.csv("../BIOGRID-ORGANISM-Homo_sapiens-4.4.221.tab3.txt", sep = "\t", header = TRUE)

gene_names = drugs$Gene.Name
gene_names = gene_names[gene_names != "" & !is.na(gene_names)]
gene_names = append(gene_names, c("PAH", "FAH", "HPD", "TAT", "CBS", "MMADHC", "MTHFR", "MTR", "MTRR", "AMT", "GLDC", "OTC", "SLC3A1", "SLC7A9", "SLC7A7", "PCCA", "PCCB", "MCEE", "MMAA", "MMAB", "MMADHC", "MMUT","BCKDHA", "BCKDHB", "DBT", "PPM1K"))
gene_names = unique(gene_names)

biogrid_df = as.data.frame(bg("interactions") %>% 
                              bg_constrain(geneList = gene_names[1]) %>% 
                              bg_get_results())

biogrid_df = as.data.frame(biogrid_df[gene_names[1] %in% biogrid_df$official_symbol_for_interactor_a || gene_names[1] %in% biogrid_df$official_symbol_for_interactor_b])

for(i in 2:length(gene_names)){
  print(paste0("BIOGRID query: ", i))
  tmp_df = as.data.frame(bg("interactions") %>% 
                               bg_constrain(geneList = gene_names[i]) %>% 
                               bg_get_results())
  
  tmp_df = as.data.frame(tmp_df[gene_names[i] %in% tmp_df$official_symbol_for_interactor_a || gene_names[i] %in% tmp_df$official_symbol_for_interactor_b])
  
  biogrid_df = rbind(biogrid_df, tmp_df)
}

write.csv(biogrid_df, "../biogrid.csv")

int_ig = interaction2igraph("../biogrid.csv")
int_ig2 = interaction2igraph2("../BIOGRID-ORGANISM-Homo_sapiens-4.4.221.tab3.txt", gene_names)

ig_drugs = csv2igraph("../drugs/pharmacologically_active.csv", "../drugs/drugbank vocabulary.csv")

ig_diseases = folder2igraph("../data")
ig_diseases = splitFolderNodes(ig_diseases, "../data")

full_ig = ig_diseases %u% ig_drugs %u% int_ig %u% int_ig2


#Diseases: Name - Gene1 - Gene2 - ... - GeneN
phenylketonuria = c("Phenylketonuria", "PAH")
tyrosinemia = c("Tyrosinemia", "FAH", "HPD", "TAT")
homocystinuria = c("HomoCystinuria", "CBS", "MMADHC", "MTHFR", "MTR", "MTRR")
hyperglycinemia = c("Hyperglycinemia", "AMT", "GLDC")
otc_deficiency = c("OTC Deficiency", "OTC")
cystinuria = c("Cystinuria", "SLC3A1", "SLC7A9")
lysinuric_protein_intolerance = c("Lysinuric Protein Intolerance", "SLC7A7")
propionic_acidemia = c("Propionic Acidemia", "PCCA", "PCCB")
methylmalonic_acidemia = c("Methylmalonic Acidemia", "MCEE", "MMAA", "MMAB", "MMADHC", "MMUT")
maple_syrup_disease = c("Maple Syrup Disease", "BCKDHA", "BCKDHB", "DBT", "PPM1K")

#All Drugs Vector
tmp = drugs$Drug.IDs
adv = c()
for(i in 1:length(tmp)){
  adv = append(adv, strsplit(tmp[i], split = "; ")[[1]])
}
adv = unique(adv)

#Building Ranking Matrix - Undirected Graphs
rnks_und_phenylketonuria = ranking(full_ig, adv, phenylketonuria[2:length(phenylketonuria)], FALSE)
rnks_und_tyrosinemia = ranking(full_ig, adv, tyrosinemia[2:length(tyrosinemia)], FALSE)
rnks_und_homocystinuria = ranking(full_ig, adv, homocystinuria[2:length(homocystinuria)], FALSE)
rnks_und_hyperglycinemia = ranking(full_ig, adv, hyperglycinemia[2:length(hyperglycinemia)], FALSE)
rnks_und_otc_deficiency = ranking(full_ig, adv, otc_deficiency[2:length(otc_deficiency)], FALSE)
rnks_und_cystinuria = ranking(full_ig, adv, cystinuria[2:length(cystinuria)], FALSE)
rnks_und_lysinuric_protein_intolerance = ranking(full_ig, adv, lysinuric_protein_intolerance[2:length(lysinuric_protein_intolerance)], FALSE)
rnks_und_propionic_acidemia = ranking(full_ig, adv, propionic_acidemia[2:length(propionic_acidemia)], FALSE)
rnks_und_methylmalonic_acidemia = ranking(full_ig, adv, methylmalonic_acidemia[2:length(methylmalonic_acidemia)], FALSE)
rnks_und_maple_syrup_disease = ranking(full_ig, adv, maple_syrup_disease[2:length(maple_syrup_disease)], FALSE)

#Building Ranking Matrix - Directed Graphs
rnks_d_phenylketonuria = ranking(full_ig, adv, phenylketonuria[2:length(phenylketonuria)], TRUE)
rnks_d_tyrosinemia = ranking(full_ig, adv, tyrosinemia[2:length(tyrosinemia)], TRUE)
rnks_d_homocystinuria = ranking(full_ig, adv, homocystinuria[2:length(homocystinuria)], TRUE)
rnks_d_hyperglycinemia = ranking(full_ig, adv, hyperglycinemia[2:length(hyperglycinemia)], TRUE)
rnks_d_otc_deficiency = ranking(full_ig, adv, otc_deficiency[2:length(otc_deficiency)], TRUE)
rnks_d_cystinuria = ranking(full_ig, adv, cystinuria[2:length(cystinuria)], TRUE)
rnks_d_lysinuric_protein_intolerance = ranking(full_ig, adv, lysinuric_protein_intolerance[2:length(lysinuric_protein_intolerance)], TRUE)
rnks_d_propionic_acidemia = ranking(full_ig, adv, propionic_acidemia[2:length(propionic_acidemia)], TRUE)
rnks_d_methylmalonic_acidemia = ranking(full_ig, adv, methylmalonic_acidemia[2:length(methylmalonic_acidemia)], TRUE)
rnks_d_maple_syrup_disease = ranking(full_ig, adv, maple_syrup_disease[2:length(maple_syrup_disease)], TRUE)

DB00147

hist(as.numeric(rnks_d_phenylketonuria[,2]), col = "green4", breaks = 40, xlab = "Mean distance from disease components", main = "Tyrosinemia", ylim = c(0,1200))
hist(as.numeric(rnks_und_phenylketonuria[,2]) + 0.2, col = "lawngreen", add = TRUE, position = "right", breaks = 40)


hist(as.numeric(rnks_d_tyrosinemia[,2]), col = "green4", breaks = 50, xlab = "Mean distance from disease components", main = "Tyrosinemia", ylim = c(0,700))
hist(as.numeric(rnks_und_tyrosinemia[,2])+0.1, col = "lawngreen", add = TRUE, position = "right", breaks = 60)


hist(as.numeric(rnks_d_homocystinuria[,2]), col = "green4", breaks = 90, xlab = "Mean distance from disease components", main = "Homocystinuria", ylim = c(0,600))
hist(as.numeric(rnks_und_homocystinuria[,2]) + 0.05, col = "lawngreen", add = TRUE, position = "right", breaks = 100)


hist(as.numeric(rnks_d_hyperglycinemia[,2]), col = "green4", breaks = 50, xlab = "Mean distance from disease components", main = "Hyperglycinemia", ylim = c(0,700))
hist(as.numeric(rnks_und_hyperglycinemia[,2])+0.1, col = "lawngreen", add = TRUE, position = "right", breaks = 60)

hist(as.numeric(rnks_d_cystinuria[,2]), col = "green4", breaks = 60, xlab = "Mean distance from disease components", main = "Cystinuria", ylim = c(0,700))
hist(as.numeric(rnks_und_cystinuria[,2])+0.1, col = "lawngreen", add = TRUE, position = "right", breaks = 60)

hist(as.numeric(rnks_d_lysinuric_protein_intolerance[,2]), col = "green4", breaks = 60, xlab = "Mean distance from disease components", main = "Lysinuric Protein Intolerance", ylim = c(0,1200))
hist(as.numeric(rnks_und_lysinuric_protein_intolerance[,2])+0.1, col = "lawngreen", add = TRUE, position = "right", breaks = 60)

hist(as.numeric(rnks_d_propionic_acidemia[,2]), col = "green4", breaks = 60, xlab = "Mean distance from disease components", main = "Propionic Acidemia", ylim = c(0,1500))
hist(as.numeric(rnks_und_propionic_acidemia[,2])+0.1, col = "lawngreen", add = TRUE, position = "right", breaks = 60)

hist(as.numeric(rnks_und_methylmalonic_acidemia[,2]), col = "lawngreen", breaks = 60, xlab = "Mean distance from disease components", main = "Methylmalonic Acidemia", ylim = c(0,400))
#hist(as.numeric(rnks_und_methylmalonic_acidemia[,2])+0.1, col = "lawngreen", add = TRUE, position = "right", breaks = 60)

hist(as.numeric(rnks_d_maple_syrup_disease[,2]), col = "green4", breaks = 100, xlab = "Mean distance from disease components", main = "Maple Syrup Disease", ylim = c(0,750))
hist(as.numeric(rnks_und_maple_syrup_disease[,2])-0.2, col = "lawngreen", add = TRUE, position = "right", breaks = 100)

hist(as.numeric(rnks_und_otc_deficiency[,2]), col = "lawngreen", breaks = 30, xlab = "Mean distance from disease components", main = "OTC Deficiency", ylim = c(0,1450), xlim = c(2,9))
hist(as.numeric(rnks_d_otc_deficiency[,2])+0.01, col = "green4", add = TRUE, position = "right", breaks = 50)














