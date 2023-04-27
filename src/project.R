library(rBiopaxParser)
library(igraph)

file_name = "C:\\Users\\raffi\\Downloads\\R-HSA-2160456_level3.owl"

#Creazione di un biopax object da trasformare in un igraph
bp = readBiopax(file_name)

df = bp$dt #dataframe of BioPax model (dataframe completo)
df_dn = df[df$property == "displayName"] #df con solo i displayName

g = pathway2Graph(bp, "Pathway1", useIDasNodenames = FALSE, verbose = TRUE)

ig = graph_from_graphnel(g, name = TRUE) #grafo estratto dal pathway (igraph)

vrtx = V(ig)

#sostituiamo i nomi dei nodi del grafo con i corrispettivi displayName del pathway 
for(i in 1:length(vrtx)){
  V(ig)[i]$name = df_dn[df_dn$id == V(ig)[i]$name]$property_value 
}

plot.igraph(ig, edge.arrow.width = 1, edge.arrow.size = 0.2)
