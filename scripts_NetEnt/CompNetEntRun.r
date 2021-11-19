# -----
# Network Entropy Analysis
# Desc: functions that compute network entropy
# Paramters
#     Perm:      the number of permutation
#     Method:    the method of network entropy calculation
#     Transform: transformation of gene expression matrix, can be "raw", "log" or "sqrt"
#     fileIn:    the name of input .RData file that will be load, the workspace should contain an list object `pin` where
#                pin$a is the network adjacency matrix and
#                pin$e is the gene expression matrix
#     fileOut:   the name of output .RData file
#     fileTxt:   text file that record each printed output during the analysis
# Example
#     Rscript CompNetEntRun.r Perm=1 Method=1 Transform=raw fileIn=... fileOut=...
# -----

### set up environment
library(tidyr)
library(dplyr)
library(Matrix)
library(rARPACK)
source("CompNetEnt.r")
source("getArgs.r")

### parameters: get arguments from command line
###             set global variable
myArgs <- getArgs()
FILE_IN  <- myArgs$fileIn
FILE_OUT <- myArgs$fileOut
FILE_TXT <- myArgs$fileTxt
TRANS    <- myArgs$Transform
PERM     <- as.numeric(myArgs$Perm)
METHOD   <- as.numeric(myArgs$Method)

cat(paste0("File In:  ", FILE_IN),  file=FILE_TXT, append=TRUE, sep="\n")
cat(paste0("File Out: ", FILE_OUT), file=FILE_TXT, append=TRUE, sep="\n")
cat(paste0("File Log: ", FILE_TXT), file=FILE_TXT, append=TRUE, sep="\n")
cat(paste0("# Perm:   ", PERM),     file=FILE_TXT, append=TRUE, sep="\n")
cat(paste0("# Method: ", METHOD),   file=FILE_TXT, append=TRUE, sep="\n")

### load workspace & get data
load(myArgs$fileIn)

# get network adjacency matrix
ADJMAT <- as.matrix(pin$a)

# get gene express and perform transformation if needed
if(TRANS == "raw"){
    EXPMAT <- pin$e
} else if (TRANSFORM == "log") {
    EXPMAT <- log(pin$e)
} else {
    EXPMAT <- sqrt(pin$e)
}

cat(paste("ADJMAT: ", "nrow", nrow(ADJMAT), "ncol", ncol(ADJMAT), "\n"), file=FILE_TXT, append=TRUE)
cat(paste("EXPMAT: ", "nrow", nrow(EXPMAT), "ncol", ncol(EXPMAT), "\n"), file=FILE_TXT, append=TRUE)

### set up random permutation
if (PERM > 1){
    numRand <- as.numeric(PERM)
    idx <- 1:nrow(EXPMAT)
    idxMat <- replicate(numRand, sample(idx))
    colnames(idxMat) <- paste0("Perm_", 1:ncol(idxMat))
    idxMat <- cbind(idx, idxMat)
    colnames(idxMat)[1] <- "Original"
} else {
    idxMat <- matrix(1:nrow(EXPMAT), ncol=1)
    colnames(idxMat) <- "Original"
} # end ifelse

### maximum entropy rate & degree
MAXSR <- CompMaxSR(ADJMAT)
KV <- apply(ADJMAT, 2, sum)

print("***KV***")
print(head(KV))
print("***MAXSR***")
print(MAXSR)


### run network entropy analysis with parallel programming
# start timer
#ptm <- proc.time()     
tt1 = Sys.time()

### helper function: check if near zero
isNearZero <- function(x, degree=10) {
    threshold <- 10^(-degree)
    return(x <= threshold & x >= -threshold)
} # end func CountZero

### prepare function to run network entropy
RunNetEnt <- function(idxSample, expMat){
    # expression value of a sample
    expVec = expMat[,idxSample]
    idx    = isNearZero(expVec, degree=6)
    expVec[idx] <- 10^(-6)
    return(CompSR(expVec,ADJMAT,KV,MAXSR,METHOD,FILE_TXT))
} # end function RunNetEnt

### start analysis
cat("Start Analysis\n",file=FILE_TXT, append=TRUE)
NetEnt <- list()
for(i in 1:ncol(idxMat)) {
    # initialization
    cat("-----\n", file=FILE_TXT, append=TRUE)
    idxRand <- idxMat[,i] # Current Permute Function
    expMatTempt <- EXPMAT[idxRand,]
    
    tempt <- list()
    for(j in 1:ncol(EXPMAT)){
        cat(psssssaste("Permutation: #", i, "| Sample:", j, "\n"), file=FILE_TXT, append=TRUE)
        tempt[[j]] <- RunNetEnt(j, expMatTempt)
    }
    NetEnt[[i]] <- tempt
    names(NetEnt[[i]]) <- colnames(EXPMAT)
    
    rm(expMatTempt)
} # end for loop
names(NetEnt) <- colnames(idxMat)

### count running time
tt2 = Sys.time()
ttd = tt2 - tt1
cat("\n", 
    "===== Running Time =====", "\n", 
    ttd, 
    file=FILE_TXT, append=TRUE)

# save the result
save.image(myArgs$fileOut)
