# -----
# Network Entropy Analysis
# Desc: functions that compute
# -----

### CompS: 
###   Computes local entropy of a node with stochastic vector p.v
CompS <- function(p.v){
    tmp.idx <- which(p.v>0); # get the value > 0; ignore zero
    S <-  - sum( p.v[tmp.idx]*log(p.v[tmp.idx]) ) # -sum(p*log(p))
    return(S);
} # end func CompS

### CompNS: 
###   Computes the normalised local entropy of a node 
###   with stochastic vector p.v
CompNS <- function(p.v){
    tmp.idx <- which(p.v>0); # get the value > 0; ignore zero
    if(length(tmp.idx)>1){
        # - 1/log(k) * sum(p*log(p)) ; k = degree(node)
        S <-  - sum( p.v[tmp.idx] * log(p.v[tmp.idx]) ) / log(length(tmp.idx));
    } else { ### one degree nodes have zero entropy, avoid singularity.
        S <- 0;
    } # end if-else
    return(S);
} # end func CompNS

### betaFn: 
###     Linear function defining gene expression values of a given sample (exp.v) 
### into low, intermediate and high expression levels according to two thresholds 
### specified by vector alpha.v. 
### Gene expression is assumed to come from Illumina or Affy arrays. 
### Returns expression values normalised between 
###     0 (low expression) and 1 (high expression) 
### with intermediate values taking a linear function between 0 and 1.
betaFn <- function(exp.v, alpha.v){
    high.idx <- which(exp.v >= alpha.v[2]);
    low.idx  <- which(exp.v <= alpha.v[1]);
    mid.idx  <- setdiff(1:length(exp.v),c(high.idx,low.idx));
    out.v <- exp.v;
    out.v[high.idx] <- 1;
    out.v[low.idx]  <- 0;
    out.v[mid.idx]  <- (exp.v[mid.idx]-alpha.v[1])/(alpha.v[2]-alpha.v[1]);
    return(out.v);
} # end func betaFn

### CompMaxSR: 
###   Computes the maximum entropy rate of a network 
###   with adjacency matrix adj.m (assumed connected).
CompMaxSR <- function(adj.m){
    require(rARPACK);
    #require(igraph)
    # find right eigenvector of adjacency matrix
    #fa <- function(x,extra=NULL) {
    #    as.vector(adj.m %*% x)
    #}
    #ap.o <- arpack(fa, options=list(n=nrow(adj.m),nev=1,which="LM"),sym=TRUE);
    #v <- ap.o$vectors;
    ap.o <- eigs(adj.m, k=1, which="LR", sym=TRUE)
    lambda <- ap.o$values;
    maxSR <- log(lambda); ### maximum entropy
    return(maxSR);
} # end func CompMaxSR

### Timer
### Desc:
###   convert time into readable form
### Example:
###   ptm <- proc.time()
###   res <- myFunc(parameters)
###   tt  <- proc.time() - ptm
###   Timer(tt)
Timer <- function(tt) {
    ttSec <- tt[1] %%  60
    ttMin <- tt[1] %/% 60 %%  60
    ttHr  <- tt[1] %/% 60 %/% 60
    return(c(ttSec, ttMin, ttHr))
} # end func Timer


# ===== Main Function =====
### CompSR: 
### Desc:
###     This function computes the entropy rate (SR) of 
###     a sample specified by a gene expression profile exp.v
###     in a network specified by an adjacency matrix adj.m. 
###     The latter is assumed connected.
### INPUT:
###     adj.m: 
###         an adjacency matrix representing a connected PPI network, 
###         with rownames/colnames of the matrix annotated to Entrez gene IDs.
###     exp.v: 
###         a genome-wide expression vector of a sample with 
###         names labeling the genes (also annotated to entrez gene IDs) 
###         of same length as rows of adj.m.
###     k.v: 
###         the degree distribution of adj.m, i.e. k.v=rowSums(adj.m)
###     local: 
###         logical, if TRUE also compute the local Shannon (normalised) entropies.
###     method: 
###         scalar specifying method to compute stochastic matrix: 
###         1 = stochastic vector of a node is independent of the 
###             node's gene expression value, 
###         2 = stochastic vector of a node is dependent on the 
###             node's gene expression value through use of the betaFn. 
###             Default is 1.
### quants.v: optionally, a vector of length 2 specifying the quantiles for defining low and high expression (only needed if method=2 is used).
### maxSR: optionally, the maximum entropy rate of the network.

### OUTPUT: a list of four objects
### PI           
### NE_loc_unnor  
### NE_loc_nor   
### NE_global    
### NE_SR (NE_global/maxSR)

CompSR <- function(expVec, adjMat, k.v, maxSR, method, FILE_TXT){
    # require package
    require(rARPACK);require(Matrix);
    
    # Method 2: the author's method
    if        (method == 2) {
        cat("  =====Apply Method 02---Author's method 2=====\n",
            file=FILE_TXT, append=TRUE)
        
        cat("    01 Calculate Network Matrix...\n",
            file=FILE_TXT, append=TRUE)
        netMat <- apply(adjMat, 2, function(x){x*expVec})
        
        cat("    02 Calculate Markov Transition Matrix...\n",
            file=FILE_TXT, append=TRUE)
        x <- colSums(netMat)
        mkvMat <- netMat %*% bandSparse(length(x), k=0, diag=list(1/x))
        mkvMat <- as.matrix(mkvMat)
        
        cat("    03 Calculate Markov Transition Matrix...\n",
            file=FILE_TXT, append=TRUE)
        quants.v <- c(0.1, 0.9)
        alpha.v  <- quantile(as.vector(expVec),probs=quants.v)
        beta.v   <- betaFn(expVec, alpha.v)
        
        mkvMat <- apply(mkvMat, 1, function(x){
            ( (1-beta.v) / k.v + beta.v * x ) * (x != 0)
        })
        mkvMat <- t(mkvMat)
        
    } else if (method == 3) {
        
        cat("  =====Apply Method 03 --- low/mediate/high expression=====",
            file=FILE_TXT, append=TRUE, sep="\n")
        
        cat("    01 Scaling Range Expression to (0,1)...",
            file=FILE_TXT, append=TRUE, sep="\n")
        quants.v <- c(0.0, 1.0)
        alpha.v  <- quantile(as.vector(expVec),probs=quants.v);
        expVec2  <- betaFn(expVec, alpha.v);
        expVec2  <- sapply(expVec, function(x){x-1})
        expVec2  <- sapply(expVec, function(x){x/2})
        
        cat("    02 Calculate Network Matrix...",
            file=FILE_TXT, append=TRUE, sep="\n")
        netMat <- apply(adjMat, 2, function(x){x*expVec2})
        
        cat("    03 Calculate Markov Transition Matrix...",
            file=FILE_TXT, append=TRUE, sep="\n")
        x <- colSums(netMat)
        mkvMat <- netMat %*% bandSparse(length(x), k=0, diag=list(1/x))
        
    } else { # method == 1
        cat("  =====Apply Method 01---Author's method 1=====\n",
            file=FILE_TXT, append=TRUE)
        
        cat("    01 Calculate Network Matrix...\n",
            file=FILE_TXT, append=TRUE)
        netMat <- apply(adjMat, 2, function(x){x*expVec})
        
        cat("    02 Calculate Markov Transition Matrix...\n",
            file=FILE_TXT, append=TRUE)
        x <- colSums(netMat)
        mkvMat <- netMat %*% bandSparse(length(x), k=0, diag=list(1/x))
        mkvMat <- as.matrix(mkvMat)
    } # end ifelse
    
    cat("=====Checking Matrix Dimension=====\n",                                   file=FILE_TXT, append=TRUE)
    cat(paste("  Exp Vector", "Len:  ", length(expVec), "\n"),                     file=FILE_TXT, append=TRUE)
    cat(paste("  Adj Matrix", "nrow: ", nrow(adjMat), "ncol", ncol(adjMat), "\n"), file=FILE_TXT, append=TRUE)
    cat(paste("  Net Matrix", "nrow: ", nrow(netMat), "ncol", ncol(netMat), "\n"), file=FILE_TXT, append=TRUE)
    cat(paste("  Mkv Matrix", "nrow: ", nrow(mkvMat), "ncol", ncol(mkvMat), "\n"), file=FILE_TXT, append=TRUE)
    
    print("=====Checking everything is alright 01 Matrix=====")
    print("***expVec***")
    print(head(expVec))
    print("***adjMat:***")
    print(adjMat[1:6, 1:6])
    print("***netMat:***")
    print(netMat[1:6, 1:6])
    print("***mkvMat:***")
    print(mkvMat[1:6, 1:6])
    
    cat("  =====Calculate Eigen Decomposition=====\n", file=FILE_TXT, append=TRUE)
    eigDec <- eigs(mkvMat, k=1, which="LR", sym=FALSE)
    
    cat("  =====Assign and Return Result=====\n", file=FILE_TXT, append=TRUE, sep="\n")
    res <- list()
    
    ### Stationary Distribution
    piVec <- abs(as.numeric(eigDec$vectors))
    res$PI<- piVec / sum(piVec)
    
    ### local network entropy
    res$NE_loc_unor <- apply(mkvMat, 2, CompS)
    res$NE_loc_nor  <- apply(mkvMat, 2, CompNS)
    
    ### global network entropy
    res$NE_global <- as.numeric(
        crossprod(
            res$PI,
            res$NE_loc_unor))
    
    ### entropy rate
    res$NE_SR <- res$NE_global/maxSR
    
    print("=====Checking everything is alright 02 Result=====")
    print("***PI***")
    print(head(res$PI))
    print("***loc_unor***")
    print(head(res$NE_loc_unor))
    print("***loc_nor***")
    print(head(res$NE_loc_nor))
    print("***global***")
    print(res$NE_global)
    return(res)
} # end func CompSR

