#input is a case phenotype and a OMIM/ORPHA disease.
##################################################################################
#INPUT
rm(list=ls())
#.libPaths( c( .libPaths(), "/home/dgrc/R/x86_64-pc-linux-gnu-library/3.6") )
args <- commandArgs(TRUE)
#indicate the path to the rda file
load(file="hpodata.rda")
library(parallel)
#casephen<-strsplit(args[1], ",")
#casephen<-c("HP:0011344", "HP:0000750", "HP:0001344", "HP:0002136", "HP:0001290", "HP:0011968", "HP:0001250","HP:0012443", "HP:0000708","HP:0000729","HP:0040082","HP:0040196","HP:0000280","HP:0430028","HP:0005280","HP:0000463","HP:0000453","HP:0000154","HP:0010800","HP:0010804","HP:0000215","HP:0000687","HP:0008897","HP:0000486","HP:0030680","HP:0001629","HP:0001636","HP:0012020","HP:0000032","HP:0001741","HP:0000034","HP:0001537","HP:0002558","HP:0002943","HP:0001776")
#entrez geneID to be compared with the case phenotype
geneID<-args[length(args)]
#geneID<-"OMIM:617616"
#number of random observations for calculating the p-value. The more the better, but it will increase the time span
iterations<-1500
#p-value cutoff
pvalue<-0.05
load("new_ancestors.Rdata")
ancestor<-new_list
##################################################################################
casephen<-vector()
aa=1
while (aa<length(args)){
  casephen <-append(casephen,as.character(args[aa]))
  aa=aa+1
}
###################################################################
 getTermListSim_edited<-function (anno1, anno2, combinemethod = "funSimMax", method = "Resnik",
    IC, verbose = FALSE)
{
    if (length(anno1) * length(anno2) == 0) {
        warning(paste("No HPO information for", anno1, anno2,
            ". Similarity set to NaN."))
        return(NaN)
    }
    if (method == "Resnik") {
        if (combinemethod == "funSimMax") {
            if (length(anno1) == length(anno2) && all(anno1 ==
                anno2)) {
                b <- as.double(IC[IC[, 1] %in% anno1, ][, 3])
                return(mean(b))
            }
            dim(anno1) <- length(anno1)
            anc1 <- unlist(apply(anno1, 1, function(x) ancestor[x]$HP))
            anc1 <- unique(union(anc1, anno1))
            anc1 <- anc1[!is.na(anc1)]
            dim(anno2) <- length(anno2)
            anc2 <- unlist(apply(anno2, 1, function(x) ancestor[x]$HP))
            anc2 <- unique(union(anc2, anno2))
            anc2 <- anc2[!is.na(anc2)]
            colMax <- fmax(anno1, anc2, IC)
            rowMax <- fmax(anno2, anc1, IC)
            return(max(colMax, rowMax))
        }
        if (combinemethod == "funSimAvg") {
            if (length(anno1) == length(anno2) && all(anno1 ==
                anno2)) {
                b <- as.double(IC[IC[, 1] %in% anno1, ][, 3])
                return(mean(b))
            }
            dim(anno1) <- length(anno1)
            anc1 <- unlist(apply(anno1, 1, function(x) ancestor[x]$HP))
            anc1 <- unique(union(anc1, anno1))
            anc1 <- anc1[!is.na(anc1)]
            dim(anno2) <- length(anno2)
            anc2 <- unlist(apply(anno2, 1, function(x) ancestor[x]$HP))
            anc2 <- unique(union(anc2, anno2))
            anc2 <- anc2[!is.na(anc2)]
            colMax <- fmax(anno1, anc2, IC)
            rowMax <- fmax(anno2, anc1, IC)
            return(0.5 * (rowMax + colMax))
        }
    }
    ker <- matrix(0, nrow = length(anno1), ncol = length(anno2),
        dimnames = list(anno1, anno2))
    kerzero <- matrix(0, nrow = length(anno1), ncol = length(anno2),
        dimnames = list(anno1, anno2))
    for (i in 1:length(anno1)) {
        for (j in 1:length(anno2)) {
            if (kerzero[i, j] == 0) {
                ker[i, j] <- calcTermSim(anno1[i], anno2[j],
                  method, IC)
                kerzero[i, j] <- 1
                if ((anno1[i] %in% anno2) && (anno2[j] %in% anno1)) {
                  ker[anno2[j], anno1[i]] <- ker[i, j]
                  kerzero[anno2[j], anno1[i]] <- 1
                }
            }
        }
    }
    if (combinemethod == "max") {
        return(max(ker))
    }
    else if (combinemethod == "mean") {
        return(mean(ker))
    }
    else if (combinemethod == "funSimAvg") {
        rowMax = mean(apply(ker, 1, max))
        colMax = mean(apply(ker, 2, max))
        return(0.5 * (rowMax + colMax))
    }
    else if (combinemethod == "funSimMax") {
        rowMax = mean(apply(ker, 1, max))
        colMax = mean(apply(ker, 2, max))
        return(max(rowMax, colMax))
    }
    else if (combinemethod == "BMA") {
        m = nrow(ker)
        n = ncol(ker)
        return((sum(apply(ker, 1, max)) + sum(apply(ker, 2, max)))/(m +
            n))
    }
    else stop(paste("getTermListSim: Unknown gene combinemethod",
        combinemethod, "!"))
}

#####################################################################
print(casephen)
print(geneID)
#Run HPOSim
library(HPO.db)
library(HPOSim)
.initialize()
combinemethod = "funSimMax"
method = "Resnik"
#IC <- get("termIC", envir = HPOSimEnv)
IC<-read.table("IC_inhouse", header = TRUE, stringsAsFactors = FALSE)
#get the terms associated to disease
Terms2<-as.character(unlist(disease[geneID]))
#if there is terms associated with disease, makes everything, otherwise output empty
if (length(Terms2)!=0){
  #input the case phenotype
  Terms1<-casephen
  Terms1 <- RemoveTermsWithoutIC(Terms1, "PA", IC)
  #calculate input max score
  Max_Ker <- getTermListSim_edited(Terms1, Terms1, combinemethod, method, IC, verbose)
  Terms2 <- RemoveTermsWithoutIC(Terms2, "PA", IC)
  #calculate Disease max score 
  disker <- getTermListSim_edited(Terms2, Terms2, combinemethod, method, IC, verbose)

#make the similarity between the gene of interest and the case phenotype
#Terms2<-as.character(unlist(disease[geneID]))
#print(geneID)
#print(Terms2)
  #Terms1 <- RemoveTermsWithoutIC(Terms1, "PA", IC)
  #Terms2 <- RemoveTermsWithoutIC(Terms2, "PA", IC)
  Case_Ker <- getTermListSim_edited(Terms1, Terms2, combinemethod, method, IC, verbose)
  
  #to modulate the p-value, we have to overlap the case phenotype with all possible phenotypes available

get_pval<-function (Terms1, geneID)
	{
	Terms2<-as.character(unlist(disease[geneID]))
    Terms1 <- RemoveTermsWithoutIC(Terms1, "PA", IC)
    Terms2 <- RemoveTermsWithoutIC(Terms2, "PA", IC)
    Ker <- getTermListSim_edited(Terms1, Terms2, combinemethod, method, IC, verbose)
	return (Ker)
	}
 
samples<-replicate(1000,sample(hpo_terms,length(casephen),replace = FALSE))
#transforma a matrix numa lista de listas
bb<-lapply(seq_len(ncol(samples)), function(i) samples[,i])


cl <- makeCluster(10)
clusterExport(cl=cl, varlist=c("disease", "RemoveTermsWithoutIC", "IC", "getTermListSim_edited", "method", "combinemethod", "ancestor", "fmax"), envir=environment())
results<-lapply(bb, get_pval, geneID)
stopCluster(cl)

  qq<-length(results[results >= Case_Ker])
  if (qq==0){
    qq=qq+1
  }
  result<-(qq+1)/iterations
  print(paste("PhenSSc", as.character(signif(Case_Ker, digits=3))," (P=", as.character(signif(result, digits=3)), "; MaxSSc",as.character(signif(Max_Ker, digits=3)),"; MaxDiseaseSSc",as.character(signif(disker, digits=3)), ")"))
}else{
  print("Sc ND ; P = ND")
}

