library(rBiopaxParser)
library(igraph)

biopax2igraph <- function(file_name){
  #Creazione di un biopax object da trasformare in un igraph
  bp = readBiopax(file_name)
  
  df = bp$dt #dataframe of BioPax model (dataframe completo)
  df_dn = df[df$property == "displayName"] #df con solo i displayName
  
  g = pathway2Graph(bp, "Pathway1", useIDasNodenames = FALSE, verbose = FALSE, splitComplexMolecules = TRUE)
  
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

csv2igraph <- function(file_name, file_dict){
  data = read.csv(file_name)
  genes = data$GenBank.Gene.ID
  prots = data$UniProt.ID
  drugs_ids = data$Drug.IDs
  edges = c()
  
  for(i in 1:length(genes)){
    drugs = strsplit(drugs_ids, "; ")[[1]]
    uniprot = prots[i]
    
    if(!(genes[i] == "")){
      for(j in length(drugs)){
        edges = append(edges, drugs[j])
        edges = append(edges, genes[i])
        #edges = append(edges, c(drugs[j], genes[i]))
        if (uniprot != ""){
          edges = append(edges, genes[i])
          edges = append(edges, uniprot)
        } 
          #edges = append(edges, c(genes[i], uniprot))
      }
    }
    else{
      for(j in length(drugs)){
        if (uniprot != ""){
          edges = append(edges, drugs[j])
          edges = append(edges, uniprot)
        } 
          #edges = append(edges, c(drugs[j], uniprot))
      }
    }
  }
  ig = make_graph(edges)
  
  return(ig)
}

splitComplexes <- function(biopax){
  complexes = biopax$dt[biopax$dt$class == "Complex"]
  complex_ids = unique(complexes$id)
  
  all_complex_components = list()
  
  for(i in 1:length(complex_ids)){
    components = c()
    df = splitComplex(biopax, complex_ids[i])
    ids = df$id
    for(j in 1:length(ids)){
      if(startsWith(ids[j], "Protein")){
        prot_df = biopax$dt[biopax$dt$id == ids[j]]
        prot_ref = sub("#", "", prot_df[prot_df$property == "entityReference"]$property_attr_value)
        
        prot_ref_df = biopax$dt[biopax$dt$id == prot_ref]
        uniprot_gene = prot_ref_df$property_value[startsWith(prot_ref_df$property_value, "UniProt")]
        
        uniprot_gene = sub("UniProt:", "", uniprot_gene)
        parts = strsplit(uniprot_gene, " ")[[1]]
        
        uniprotID = parts[1]
        gene = parts[2]
        components = append(components, c(gene, uniprotID))
      }
      else{
        components = append(components, df[df$id == ids[j]]$name)
      }
    }
    
    all_complex_components[[i]] = components
  }
  
  return(all_complex_components)
}

splitNode <- function(file_name){
  biopax = readBiopax(file_name)
  
  complexes_components = splitComplexes(biopax)
  complexes = biopax$dt[biopax$dt$class == "Complex"]
  complex_names = complexes[complexes$property == "displayName"]$property_value 
}
