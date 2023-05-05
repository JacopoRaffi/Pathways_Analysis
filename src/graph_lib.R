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
  last_i_used = 1 #conteggio numero sostituzioni
  
  for(i in 2:length(names)){
    if(names[i] %in% names[1:i-1]){ #nome del node già incontrato
      ids = append(ids, ids[match(names[i], names[1:i-1])]) #inseriamo id del nodo già incontrato
    }
    else{
      last_i_used = last_i_used + 1
      ids = append(ids, last_i_used) #inseriamo id corretto secondo la progressione
    }
  }
  
  #fondiamo i nodi con lo stesso nome mantenendo una sola occorrenza per quel nome
  ig = contract(ig, ids, vertex.attr.comb = list(name="first"))
  
  return(ig)   
}

folder2igraph <- function(dir_name){
  files = list.files(dir_name, full.names = TRUE) #ottengo la lista di file biopax (path relativi)
  
  fullgraph = biopax2igraph(files[1])
  
  for(i in 2:length(files)){
    biopax2igraph(files[i]) 
      fullgraph = fullgraph %u% biopax2igraph(files[i]) 
  }
  
  return(fullgraph)
}

csv2igraph <- function(file_name){
  data = read.csv2(file_name)#leggo file csv
  
  targets = data$Targets #prendo i target della medicina
  graph_vect = c()
  
  #creo il vettore degli archi
  for(i in 1:length(targets)){
    graph_vect = append(graph_vect, c(data$Name[1], targets[i]))
  }
  
  ig = make_graph(graph_vect)
  
  return(ig)
}

foldercsv2igraph <- function(dir_name){
  files = list.files(dir_name, full.names = TRUE) #ottengo la lista di file csv (path relativi)
  
  fullgraph = csv2igraph(files[1])
  
  for(i in 2:length(files)){
    csv2igraph(files[i]) 
    fullgraph = fullgraph %u% csv2igraph(files[i]) 
  }
  
  return(fullgraph)
}