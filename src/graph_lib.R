library(rBiopaxParser)
library(igraph)

biopax2igraph <- function(file_name){
  #used for converting BioPax3 files in igraph
  
  bp = readBiopax(file_name)
  df = bp$dt # complete dataframe extracted from the BioPax3
  df_dn = df[df$property == "displayName"] # dataframe filtering for mantaining only the rows with "displayName" as property
  
  # Pathway1 -> identifier of the pathway to convert from the dataframe 
  g = pathway2Graph(bp, "Pathway1", useIDasNodenames = FALSE, verbose = FALSE, splitComplexMolecules = TRUE)
  
  ig = graph_from_graphnel(g, name = TRUE) #graph extracted from the pathway (igraph)
  
  # substitution of vertex names with "displayName" -> uniformation of the data
  vrtx = V(ig)
  for(i in 1:length(vrtx)){
    V(ig)[i]$name = df_dn[df_dn$id == V(ig)[i]$name]$property_value 
  }
  
  # merging of nodes which are modeling the same component
  ig = contractNodes(ig)
  
  return(ig)   
}

folder2igraph <- function(dir_name){
  #used for converting all the BioPax3 in a folder to igraph
  
  files = list.files(dir_name, full.names = TRUE)
  fullgraph = biopax2igraph(files[1])
  
  for(i in 2:length(files)){
    biopax2igraph(files[i]) 
      fullgraph = fullgraph %u% biopax2igraph(files[i]) 
  }
  
  return(fullgraph)
}

csv2igraph <- function(file_name, file_dict){
  #used for converting CSVs of drugs'interactions in igraph format
  
  data = read.csv(file_name)
  genes = data$Gene.Name
  prots = data$UniProt.ID
  drugs_ids = data$Drug.IDs
  edges = c()
  
  for(i in 1:length(genes)){
    drugs = strsplit(drugs_ids[i], "; ")[[1]]
    print(drugs)
    uniprot = prots[i]
    
    if(!(is.na(genes[i]))){
      if(!(genes[i] == "")){
        for(j in 1:length(drugs)){
          edges = append(edges, drugs[j])
          edges = append(edges, genes[i])
        }
      }
    }
    for(j in 1:length(drugs)){
      if (uniprot != ""){
        edges = append(edges, drugs[j])
        edges = append(edges, uniprot)
      }
    }
  }
  ig = make_graph(edges)
  
  return(ig)
}

splitComplexes <- function(biopax){
  # used for identifying complexes and their components
  # The information retrieved from this processing are used to split the actual nodes in the already built graph
  
  complexes = biopax$dt[biopax$dt$class == "Complex"]
  complex_ids = unique(complexes$id)
  
  all_complex_components = list()
  
  for(i in 1:length(complex_ids)){
    components = c()
    df = splitComplex(biopax, complex_ids[i])
    ids = df$id
    
    if(length(ids) > 0){
      for(j in 1:length(ids)){
        if(startsWith(ids[j], "Protein")){
          prot_df = biopax$dt[biopax$dt$id == ids[j]]
          prot_ref = sub("#", "", prot_df[prot_df$property == "entityReference"]$property_attr_value)
          
          prot_ref_df = biopax$dt[biopax$dt$id == prot_ref]
          uniprot_gene = prot_ref_df$property_value[startsWith(prot_ref_df$property_value, "UniProt")]
          
          if(identical(uniprot_gene,character(0))){
            next
          }
          
          uniprot_gene = sub("UniProt:", "", uniprot_gene)
          parts = strsplit(uniprot_gene, " ")[[1]]
          
          uniprotID = parts[1]
          gene = parts[2]
          components = append(components, c(gene, uniprotID))
        }
        else{
          components = append(components, df[df$id == ids[j],]$name)
        }
      }
    }
    if(length(components) == 0){
      next
    }
    all_complex_components[[i]] = components
  }
  
  return(all_complex_components)
}

splitNodes <- function(graph, file_name){
  # used to split graph's nodes which are modeling complexes
  
  biopax = readBiopax(file_name)
  complexes = biopax$dt[biopax$dt$class == "Complex"]
  
  if(nrow(complexes) == 0){
    return(graph)
  }
  
  # complexes_components contains the information useful to split the nodes
  complexes_components = splitComplexes(biopax)
  complex_names = complexes[complexes$property == "displayName"]$property_value 
  
  
  for(i in 1:length(complex_names)){
    components = complexes_components[[i]]
    
    if(!(complex_names[i] %in% V(graph)$name)){
      next
    }
    
    i_edges = incident(graph, complex_names[i], mode = "in")
    i_nodes = ends(graph, E(graph)[i_edges])[,1]
    
    o_edges = incident(graph, complex_names[i], mode = "out")
    o_nodes = ends(graph, E(graph)[o_edges])[,2]
    
    # thanks to this iteration each component of a complex models all the interactions of the complex
    if(length(components) > 0){
      for(j in 1:length(components)){
        new_edges = c()
        
        graph = add_vertices(graph, 1, name = components[j]) 
        
        if(length(i_nodes) > 0){
          for(k in 1:length(i_nodes)){
            
            new_edges = append(new_edges, c(i_nodes[k], components[j]))
          }
          
          graph = add_edges(graph, new_edges)
        }
        
        if(length(o_nodes) > 0){
          new_edges = c()
          for(k in 1:length(o_nodes)){
            
            new_edges = append(new_edges, c(components[j], o_nodes[k])) 
          }
          
          graph = add_edges(graph, new_edges)
        }
      }  
    }
    
    graph = delete_vertices(graph, complex_names[i])
  }
  
  return(graph)
}

splitFolderNodes <- function(graph, dir_name){
  # used to split complexes of every pathway found in a folder
  files = list.files(dir_name, full.names = TRUE)
  
  for(i in 1:length(files)){
    graph = splitNodes(graph, files[i])
  }
  
  graph = contractNodes(graph)
  
  return(graph)
}

interaction2igraph <- function(file_name){
  # used for translating BioGrid interactions in actual edges of the graph
  df = read.csv(file_name)
  
  int_a = df$official_symbol_for_interactor_a
  int_b = df$official_symbol_for_interactor_b
  
  edges = c()
  
  for(i in 1:length(int_a)){
    edges = append(edges, c(int_a[i], int_b[i]))
  }
  
  ig = make_graph(edges)
  
  return(ig)
}

interaction2igraph2 <- function(file_name, gene_names){
  # used for translating BioGrid interactions in actual edges of the graph
  # the framework contains two version useful for the different formatting of those files
  
  df = read.csv(file_name, sep = "\t", header = TRUE)
  
  int_a = df$Official.Symbol.Interactor.A
  int_b = df$Official.Symbol.Interactor.B
  
  edges = c()
  
  for(i in 1:length(int_a)){
    print(paste0("INT2: ", i))
    if(int_a[i] %in% gene_names | int_b[i] %in% gene_names){
      edges = append(edges, c(int_a[i], int_b[i]))
    }
  }
  
  ig = make_graph(edges)
  
  return(ig)
}

contractNodes <- function(ig){
  # used to merge the nodes which are modeling the same "biological entity"
  
  vrtx = V(ig)
  
  names = V(ig)$name # nodes'names
  ids = c(1) # vector of IDs necessary for igraph.contract
  last_i_used = 1 # count of substitutions
  
  for(i in 2:length(names)){
    if(names[i] %in% names[1:i-1]){ # if the node has already been found (same name)
      ids = append(ids, ids[match(names[i], names[1:i-1])]) # insert the IDs referring to the already encountered node
    }
    else{
      last_i_used = last_i_used + 1
      ids = append(ids, last_i_used) # elsewise insert the ID according to the progression
    }
  }
  
  # this call result in the marging of all the nodes with the same displayName, preserving only one cumulative instance
  ig = contract(ig, ids, vertex.attr.comb = list(name="first"))
  
  return(ig)
}

ranking <- function(ig, drugs, genes, directed){
  
  # ig -> igraph representation of the global network
  # drugs -> vector containing all drugs
  # genes -> vector containing all the genes of the studied disease
  # directed -> boolean: TRUE: consider edges directions, FALSE: elsewise
  
  m = "all"
  if(directed)
    m = "out"
  
  # igraph library function that calculate the shortest path lenghts between two sets of nodes
  dist = distances(ig, v = drugs, to = genes, mode = m)
  
  # generation of a matrix containing association between drugs and mean distances from the genes of the disease
  # columns: 1 -> drug_id ; 2 -> mean_distance; 3 -> rank; 4 -> percentage of drugs falling under this rank
  means = matrix(nrow = length(drugs), ncol = 4)
  for(i in 1:nrow(dist)){
    means[i,1] = drugs[i]
    means[i,2] = mean(dist[i,])
  }
  
  # sorting of the matrix based on the mean distances: -> the drugs column represent the ranking
  means = means[order(as.numeric(means[,2]), decreasing = FALSE),]

  # populating the rank column
  rank = 1
  last_mean = means[1,2]
  rows = nrow(means)
  
  for(i in 1:rows){
    if(as.numeric(means[i,2]) > last_mean){
      rank = rank + 1
      last_mean = as.numeric(means[i,2])
    }
    means[i,3] = rank
  }
  
  # populating the "percentage of drugs falling under this rank" column
  rank_percentage = c()
  for(i in 1:rank){
    curr_rank_rows = nrow(means[as.numeric(means[,3]) == i,])
    lower_rank_rows = nrow(means[as.numeric(means[,3]) > i,])
    rank_percentage = append(rank_percentage, round( ( (lower_rank_rows + (curr_rank_rows/2)) / rows) * 100, 2))
  }
  
  if(genes[1] == "PAH" | genes[1] == "HPD" | genes[1] == "SLC7A7")
    rank_percentage = append(rank_percentage, round((2219 / rows) * 100, 2), after = 0)

  for(i in 1:rows){
    means[i,4] = rank_percentage[as.numeric(means[i,3])]
  }
  
  return (means);
}

test1Ranking <- function(akdv, arv){
  # after performing ttest between rank1 and rank2, if pval is reasonably low ( < 0,05), the test assess that the
  # already known drugs achieve a better ranking on their target diseases respect that on NON target diseases
  
  rank1 = c() # ranks achieved by known drugs in each of their target diseases
  rank2 = c() # ranks achieved by known drugs in each of their NON target diseases
  
  for(i in 1:length(akdv)){
    drugID = akdv[[i]][[1]]
    
    for(j in 1:length(arv)){
      tmp = arv[[j]]
      maxRank = as.numeric(tmp[nrow(tmp), 3])
      
      if(j %in% akdv[[i]]){ # drug i is known for disease j
        rank1 = append(rank1, as.numeric(tmp[tmp[,1] == drugID, 3]) / maxRank)
      }
      else{ 
        rank2 = append(rank2, as.numeric(tmp[tmp[,1] == drugID, 3]) / maxRank)
      }
    }
  }
  
  test = t.test(rank1, rank2)
  return (test)
}

test2Ranking <- function(akdv, arv){
  # after performing ttest between rank1 and rank3, if pval is reasonably low ( < 0,05), the test assess that the
  # already known drugs achieve a better ranking respect that random drugs
  
  rank1 = c() # ranks achieved by known drugs in each of their target diseases
  rank3 = c() # 10 iterations of the rankings for each disease of 14 random drugs (from the total 2220 drugs list)
  
  for(i in 1:length(akdv)){
    drugID = akdv[[i]][[1]]
    
    for(j in 1:length(arv)){
      tmp = arv[[j]]
      maxRank = as.numeric(tmp[nrow(tmp), 3])
      
      if(j %in% akdv[[i]]){ # drug i is known for disease j
        rank1 = append(rank1, as.numeric(tmp[tmp[,1] == drugID, 3]) / maxRank)
      }
    }
  }
  
  for(i in 1:10){
    randomDrugs = c()
    alreadyTaken = c()
    j = 1
    while (j <= 14){
      tmp = round(runif(1, 1, 2220))
      if(! (tmp %in% alreadyTaken)){
        randomDrugs = append(randomDrugs, arv[[1]][tmp, 1]) 
        j = j + 1
      }
    }
    
    for(k in 1:length(randomDrugs)){
      drugID = randomDrugs[k]
      
      for(m in 1:length(arv)){
        tmp = arv[[m]]
        maxRank = as.numeric(tmp[nrow(tmp), 3])
        rank3 = append(rank3, as.numeric(tmp[tmp[,1] == drugID, 3]) / maxRank)
      }
    }
  }
  
  print(rank3)
  
  test = t.test(rank1, rank3)
  return (test)
}
  

