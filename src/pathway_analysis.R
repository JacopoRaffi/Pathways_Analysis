source("graph_lib.R")
library(biogridr)
library(dplyr)

# *** Download and Integration of diseases' and drugs' pathways ***

#Generation of graphs starting form Reactome's Pathways
ig_diseases = folder2igraph("../data")
#Data Processing: splitting of complex nodes in components
ig_diseases = splitFolderNodes(ig_diseases, "../data")
#Dataframe containing all the pharmacologically active drugs retrieved from DrugBank
drugs = read.csv("../drugs/pharmacologically_active.csv")
ig_drugs = csv2igraph("../drugs/pharmacologically_active.csv", "../drugs/drugbank vocabulary.csv")
#Extraction of drugs'tragets from DrugBank
gene_names = drugs$Gene.Name
#Filtering
gene_names = gene_names[gene_names != "" & !is.na(gene_names)]
#Redundant Integration of the genes identifying the studied diseases
gene_names = append(gene_names, c("PAH", "FAH", "HPD", "TAT", "CBS", "MMADHC", "MTHFR", "MTR", "MTRR", "AMT", "GLDC", "OTC", "SLC3A1", "SLC7A9", "SLC7A7", "PCCA", "PCCB", "MCEE", "MMAA", "MMAB", "MMADHC", "MMUT","BCKDHA", "BCKDHB", "DBT", "PPM1K"))
#Filtering
gene_names = unique(gene_names)


# *** Query to Biogrid in order to extract the (distance +1) interactions of our pathways'genes ***

#Key for accessing Biogrid from R
bg_get_key("Jacopo", "Raffi", "j.raffi@studenti.unipi.it", "metabolic_disease")

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

#Generation of the graph starting from BioGrid's Interaction obtained by query 
int_ig = interaction2igraph("../biogrid.csv")
#Generation of the graph starting from BioGrid's Interaction for human organism
int_ig2 = interaction2igraph2("../BIOGRID-ORGANISM-Homo_sapiens-4.4.221.tab3.txt", gene_names)


# *** Generation of the complete graph ***
full_ig = ig_diseases %u% ig_drugs %u% int_ig %u% int_ig2


#Diseases: Name - Gene1 - Gene2 - ... - GeneN
phenylketonuria = c("Phenylketonuria", "PAH")
tyrosinemia1 = c("Tyrosinemia T1", "FAH")
tyrosinemia2 = c("Tyrosinemia T2", "TAT")
tyrosinemia3 = c("Tyrosinemia T3", "HPD")
homocystinuria = c("HomoCystinuria", "CBS", "MMADHC", "MTHFR", "MTR", "MTRR")
hyperglycinemia = c("Hyperglycinemia", "AMT", "GLDC")
otc_deficiency = c("OTC Deficiency", "OTC")
cystinuria = c("Cystinuria", "SLC3A1", "SLC7A9")
lysinuric_protein_intolerance = c("Lysinuric Protein Intolerance", "SLC7A7")
propionic_acidemia = c("Propionic Acidemia", "PCCA", "PCCB")
methylmalonic_acidemia = c("Methylmalonic Acidemia", "MCEE", "MMAA", "MMAB", "MMADHC", "MMUT")
maple_syrup_disease = c("Maple Syrup Disease", "BCKDHA", "BCKDHB", "DBT", "PPM1K")

#Data Processing: unique vector  for all drugs'IDs
tmp = drugs$Drug.IDs
adv = c()
for(i in 1:length(tmp)){
  adv = append(adv, strsplit(tmp[i], split = "; ")[[1]])
}
adv = unique(adv)

full_ig = read_graph("../graph/full_ig.graphml", format = "graphml")
#Building Ranking Matrix - Undirected Graphs
rnks_und_phenylketonuria = ranking(full_ig, adv, phenylketonuria[2:length(phenylketonuria)], FALSE)
rnks_und_tyrosinemia1 = ranking(full_ig, adv, tyrosinemia1[2:length(tyrosinemia1)], FALSE)
rnks_und_tyrosinemia2 = ranking(full_ig, adv, tyrosinemia2[2:length(tyrosinemia2)], FALSE)
rnks_und_tyrosinemia3 = ranking(full_ig, adv, tyrosinemia3[2:length(tyrosinemia3)], FALSE)
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
rnks_d_tyrosinemia1 = ranking(full_ig, adv, tyrosinemia1[2:length(tyrosinemia1)], TRUE)
rnks_d_tyrosinemia2 = ranking(full_ig, adv, tyrosinemia2[2:length(tyrosinemia2)], TRUE)
rnks_d_tyrosinemia3 = ranking(full_ig, adv, tyrosinemia3[2:length(tyrosinemia3)], TRUE)
rnks_d_homocystinuria = ranking(full_ig, adv, homocystinuria[2:length(homocystinuria)], TRUE)
rnks_d_hyperglycinemia = ranking(full_ig, adv, hyperglycinemia[2:length(hyperglycinemia)], TRUE)
rnks_d_otc_deficiency = ranking(full_ig, adv, otc_deficiency[2:length(otc_deficiency)], TRUE)
rnks_d_cystinuria = ranking(full_ig, adv, cystinuria[2:length(cystinuria)], TRUE)
rnks_d_lysinuric_protein_intolerance = ranking(full_ig, adv, lysinuric_protein_intolerance[2:length(lysinuric_protein_intolerance)], TRUE)
rnks_d_propionic_acidemia = ranking(full_ig, adv, propionic_acidemia[2:length(propionic_acidemia)], TRUE)
rnks_d_methylmalonic_acidemia = ranking(full_ig, adv, methylmalonic_acidemia[2:length(methylmalonic_acidemia)], TRUE)
rnks_d_maple_syrup_disease = ranking(full_ig, adv, maple_syrup_disease[2:length(maple_syrup_disease)], TRUE)



hist(as.numeric(rnks_und_phenylketonuria[,2]), col = "orchid", breaks = 20, xlab = "Mean distance from disease components", main = "Phenylketonuria", ylim = c(0,1200))

hist(as.numeric(rnks_und_tyrosinemia1[,2]), col = "orchid", breaks = 20, xlab = "Mean distance from disease components", main = "Tyrosinemia Type 1", ylim = c(0,1200))

hist(as.numeric(rnks_und_tyrosinemia2[,2]), col = "orchid", breaks = 20, xlab = "Mean distance from disease components", main = "Tyrosinemia Type 2", ylim = c(0,1200))

hist(as.numeric(rnks_und_tyrosinemia3[,2]), col = "orchid", breaks = 20, xlab = "Mean distance from disease components", main = "Tyrosinemia Type 3", ylim = c(0,1200))

hist(as.numeric(rnks_und_homocystinuria[,2]), col = "orchid", breaks = 20, xlab = "Mean distance from disease components", main = "Homocystinuria", ylim = c(0,800))

hist(as.numeric(rnks_und_hyperglycinemia[,2]), col = "orchid", breaks = 20, xlab = "Mean distance from disease components", main = "Hyperglycinemia", ylim = c(0,800))

hist(as.numeric(rnks_und_otc_deficiency[,2]), col = "orchid", breaks = 20, xlab = "Mean distance from disease components", main = "OTC Deficiency", ylim = c(0,1200))

hist(as.numeric(rnks_und_cystinuria[,2]), col = "orchid", breaks = 20, xlab = "Mean distance from disease components", main = "Cystinuria", ylim = c(0,800))

hist(as.numeric(rnks_und_lysinuric_protein_intolerance[,2]), col = "orchid", breaks = 20, xlab = "Mean distance from disease components", main = "Lysinuric Protein Intolerance", ylim = c(0,1200))

hist(as.numeric(rnks_und_propionic_acidemia[,2]), col = "orchid", breaks = 20, xlab = "Mean distance from disease components", main = "Propionic Acidemia", ylim = c(0,1500))

hist(as.numeric(rnks_und_methylmalonic_acidemia[,2]), col = "orchid", breaks = 20, xlab = "Mean distance from disease components", main = "Methylmalonic Acidemia", ylim = c(0,1000))

hist(as.numeric(rnks_und_maple_syrup_disease[,2]), col = "orchid", breaks = 20, xlab = "Mean distance from disease components", main = "Maple Syrup Disease", ylim = c(0,800))


# *** STATISTICAL VALIDATION *** #

#Vector containing all rankings
arv = list(rnks_und_phenylketonuria, rnks_und_tyrosinemia1, rnks_und_tyrosinemia2, rnks_und_tyrosinemia3, rnks_und_homocystinuria, 
        rnks_und_hyperglycinemia, rnks_und_otc_deficiency, rnks_und_cystinuria, rnks_und_lysinuric_protein_intolerance, 
        rnks_und_propionic_acidemia, rnks_und_methylmalonic_acidemia, rnks_und_maple_syrup_disease)

#Drugs: DrugID - Disease1 - Disease2 - ... - DiseaseN
d1 = list("DB00131", 1, 3, 4, 7)
d2 = list("DB01373", 1, 2, 3, 4, 7)
d3 = list ("DB00360", 1, 2)
d4 = list ("DB00368", 2)
d5 = list ("DB00668", 2)
d6 = list ("DB00988", 2)
d7 = list ("DB01235", 2)
d8 = list ("DB00115", 5)
d9 = list ("DB00200", 5)
d10 = list ("DB03614", 5)
d11 = list ("DB03904", 7, 8)
d12 = list ("DB13146", 9)
d13 = list ("DB00147", 10, 12)
d14 = list ("DB00348", 3)

#Vector containing all known drugs for analysed diseases
akdv = list(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14)

# TESTS #
# - test description is found commented inside the test function (./graph_lib.R)

test1 = test1Ranking(akdv, arv)
if(test1$p.value < 0.05) print("TEST 1 PASSED!")

test2 = test2Ranking(akdv, arv)
print(test2$p.value)
if(test2$p.value < 0.05) print("TEST 2 PASSED!")








