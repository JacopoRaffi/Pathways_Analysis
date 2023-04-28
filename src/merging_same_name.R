library(rBiopaxParser)
library(igraph)

biopax2igraph <- function(file_name){
  #Creazione di un biopax object da trasformare in un igraph
  bp = readBiopax(file_name)
  
  df = bp$dt #dataframe of BioPax model (dataframe completo)
  df_dn = df[df$property == "displayName"] #df con solo i displayName
  
  g = pathway2Graph(bp, "Pathway1", useIDasNodenames = FALSE, verbose = FALSE)
  
  ig = graph_from_graphnel(g, name = TRUE) #grafo estratto dal pathway (igraph)
  
  vrtx = V(ig)
  
  #sostituiamo i nomi dei nodi del grafo con i corrispettivi displayName del pathway 
  for(i in 1:length(vrtx)){
    V(ig)[i]$name = df_dn[df_dn$id == V(ig)[i]$name]$property_value 
  }
  
  names = V(ig)$name #vettore dei nomi dei nodi del grafo
  ids = c(1) #vettore dei futuri id del grafo
  sub_num = 0 #conteggio numero sostituzioni
  
  for(i in 2:length(names)){
    if(names[i] %in% names[1:i-1]){ #nome del node già incontrato
      sub_num = sub_num + 1 
      ids = append(ids, match(names[i], names[1:i-1])) #inseriamo id del nodo già incontrato
    }
    else{
      ids = append(ids, i - sub_num) #inseriamo id corretto secondo la progressione
    }
  }
  
  #fondiamo i nodi con lo stesso nome mantenendo una sola occorrenza per quel nome
  ig = contract(ig, ids, vertex.attr.comb = list(name="first"))
  
  return(ig)   
}

#plot.igraph(ig, edge.arrow.width = 0.1, edge.arrow.size = 0.1, vertex.size = 10, vertex.label.cex = 0.6)
