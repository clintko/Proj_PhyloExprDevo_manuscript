# -----
# GSVA Analysis
# Desc: The is a script that wrap up the GSVA analysis
# Paramters
# Paramters
#     Perm:      the number of permutation
#     fileIn:    the name of input .RData file that will be load, the workspace should contain a gene expression matrix `EXPMAT`
#     fileOut:   the name of output .RData file
#     fileTxt:   text file that record each printed output during the analysis
# Example
#     Rscript CompNetEntRun.r Perm=1 fileIn=... fileOut=... fileTxt
# -----

library(dplyr)
library(tidyr)
library(GSVA)
library(stringr)
library(foreach)
library(doParallel)
library(parallel)
source("getArgs.r")

### parameters: get arguments from command line
###             set global variable
myArgs   <- getArgs()
FILE_IN  <- myArgs$fileIn
FILE_OUT <- myArgs$fileOut
FILE_TXT <- myArgs$fileTxt
PERM     <- as.numeric(myArgs$Perm)
print(paste0("File In:  ", FILE_IN))
print(paste0("File Out: ", FILE_OUT))
print(paste0("File Log: ", FILE_TXT))
print(paste0("# Perm:   ", PERM))

### load and prepare data
# Load variable: EXPMAT(matrix), GS(list)
load(myArgs$fileIn)

### set up random permutation
numRand <- as.numeric(PERM)
idx <- 1:nrow(EXPMAT)
idxMat <- replicate(numRand, sample(idx))
colnames(idxMat) <- paste0("Perm_", 1:ncol(idxMat))
idxMat <- cbind(idx, idxMat)
colnames(idxMat)[1] <- "Original"

### function for running GSVA
runGSVA <- function(i, expMat, method){
    # record
    cat(paste("Permutation: #", i), file=FILE_TXT, append=TRUE, sep = "\n")
    # index
    idx <- idxMat[,i]
    rownames(expMat) <- GENEID[idx]
    
    # Perform GSVA analysis with method = 1 or 2
    if (method == 1){
        res <- gsva(expMat, GS, min.sz=10, mx.diff=FALSE, verbose=FALSE)$es.obs
    } else { # method ==2
        res <- gsva(expMat, GS, min.sz=10, mx.diff=TRUE,  verbose=FALSE)$es.obs
    } # end if-else
    
    # arrange data
    res <- data.frame(res)
    res$GeneSet <- rownames(res)
    res$Type1 <- rep(colnames(idxMat)[i], nrow(res))
    return(res)
}

### set up parallel computation
numCores <- detectCores()-2
cl <- makeCluster(numCores)
registerDoParallel(cl)

### run GSVA analysis
# analysis
resGSVAmax <- foreach(i = 1:ncol(idxMat)) %do% runGSVA(i, EXPMAT, method=1)
resGSVAdif <- foreach(i = 1:ncol(idxMat)) %do% runGSVA(i, EXPMAT, method=2)
# data name
names(resGSVAmax) <- colnames(idxMat)
names(resGSVAdif) <- colnames(idxMat)
# row bind the data
resGSVAmax <- do.call(rbind, resGSVAmax)
resGSVAdif <- do.call(rbind, resGSVAdif)

### close parallel computation
stopCluster(cl)

### save the result
save.image(myArgs$fileOut)