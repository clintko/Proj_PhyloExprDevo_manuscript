# -----
# Integrate expression matrix and adjacency matrix
# Desc: helper function to integrate the expression and network
# -----

##############################################
# giant.component
# return: 
#     extract largest connected component
#     the largest cc of the input graph
##############################################
giant.component <- function(graph, ...) {
    cl <- clusters(graph, ...)
    induced_subgraph (graph, which(cl$membership == which.max(cl$csize)))
} # end function giant.component


####################################################
# DoIntegPIN
# return: 
#     integrate the network and expression data
####################################################
DoIntegPIN <- function(edgeList, exprDat) {
    # separate the matrix from expression matrix
    gene   <- exprDat$Gene
    expMat <- exprDat %>% dplyr::select(-Gene) %>% as.matrix
    
    # filter based on expression data
    edge <- edgeList %>% 
        filter(Gene_A %in% gene) %>%
        filter(Gene_B %in% gene) %>%
        distinct
    
    # get largest connected component
    g <- giant.component(
        simplify(
            graph.data.frame(
                edge, 
                directed=FALSE)))
    
    # combine the data
    #node   <- as.numeric(names(V(g)))
    node   <- as.character(names(V(g)))
    idx    <- match(node, gene)
    expMat <- expMat[idx,]
    
    return(list(a=get.adjacency(g),
                e=expMat,
                g=g,
                gene=gene[idx]))
} # end function DoIntegPIN