# -----------------------------
# FUNCTION: MEGA
# -----------------------------
# Description: MEGA was developed to identify predefined gene 
# sets (e.g. genes involved in the same pathway, or predisposing 
# to specific diseases) that show a significantly higher number 
# of mutations in a group of samples as compared to another 
# group of samples.
#
# Inputs:
# A and B: Boolean matrices of mutations. Coloums are samples, while rows are
# mutations. The first coloumn must always contain the name of thegene in which
# the mutation fall.
#
# Example:
#
# +-------------------------+
# Symbol  S1    S2    S3
# MAST2   TRUE  FALSE FALSE
# MAST2   TRUE  FALSE FALSE
# ABCA10  FALSE FALSE FALSE
# SLC4A9  TRUE  TRUE  FALSE
# ZNF572  FALSE FALSE FALSE
# ASPM    FALSE FALSE TRUE
# CSMD2   FALSE FALSE TRUE
# +-------------------------+
#
#
# gene.sets: List of gene sets. Each element of the list is set of genes and
# the name of each element of the list must be the name of the gene set.
# 
# Example:
# gene.sets = list(g1=c("MAST2","ABCA10","ASPM"),g2=c("TP53","ZNF572","MYC"))
#
# fdr_th: false discovery rate threshold (default is 0.1)
#
# bootstrapping: If equal to true the bootstrapping strategy is performed 
# to assess the effect of the sample size
# 
# nsim: number of iterations used in the bootsrapping strategy (default is 1000)
#
MEGA = function(A,B=NULL,gene.sets,fdr_th=0.1,bootstrapping=F,montecarlo=F,MC.genome="HG19",nsim=1000,test="W") {
  # print info
  # -----------
  print.logo()
  cat("\n========================================\n")
  cat("Input parameters:\n")
  cat("FDR threshold:",fdr_th,"\n")
  cat(paste("Number of Gene Sets:",length(gene.sets),"\n"))
  if (bootstrapping) {
    cat("Bootstrapping: YES\n")
    cat("Number of iterations:",nsim,"\n")
  } else {
    cat("Bootstrapping: NO\n")
    cat("Number of iterations: 0\n")
  }
  if (montecarlo) {
    cat("Monte Carlo: YES\n")
    cat("Number of iterations:",nsim,"\n")
  } else {
    cat("Monte Carlo: NO\n")
    cat("Number of iterations: 0\n")
  }
  cat("========================================\n")
  
  if (!montecarlo)
  {
  
  cat("\nStep 1: Mutation Gene Set Enrichement Analysis\n")
  
  # Execute the enrichement
  # ------------------------
  pb = txtProgressBar(min = 0, max = length(gene.sets), initial = 0,style=3)
  p = sapply(1:length(gene.sets), function(x,z=pb,a=A,b=B,gs=gene.sets,t=test) {
    setTxtProgressBar(z,x)
    r = MEGA.core(A,B,gs[[x]],t)
    names(r) = names(gs)[x]
    return(r)
  })
  
  res = data.frame("gene.set"=names(p),p.value=p,fdr= p.adjust(p,method = "fdr"),stringsAsFactors = F)
  res = res[order(res$fdr),]
  rownames(res) = NULL
  
  # Execute Bootstrapping if required
  # ----------------------------------
  if (bootstrapping) {
    if (sum(res$fdr<fdr_th)) {
      cat("\n\nStep 2: Bootstrapping for",sum(res$fdr<fdr_th),"significant gene sets\n")
      bs = MEGA.bootstrapping(A,B,gene.sets[res$gene.set[res$fdr<fdr_th]],nsim,test)
      res$success_percentage = NA
      if (!is.null(dim(bs))) {
        res$success_percentage[res$fdr<fdr_th] = apply(bs, 1, function(x,n=nsim) sum(x<0.05)*100/n)
      } else {
        res$success_percentage[res$fdr<fdr_th] = sum(bs<0.05)*100/nsim
      }
    } else {
      cat("\n\n Step 2: Bootstrapping can not be performed. No significant gene sets\n")
    }
  } 
  }
  
  if (montecarlo)
  {
    require(parallel)
    load(paste("./RData/",MC.genome,".gene.cds.length.RData",sep = ""))
    cpus = detectCores()
    if (!is.null(cpus))
    {
      cpus = cpus - 1
    } else {
      cpus = 2
    }
    res = MEGA.MC(A,gene.sets,gene.cds.length,nsim,cpus)
    res$fdr = p.adjust(res$p.value,method = "fdr")
  }
  
  cat("\n\nResults:\n")
  cat("Significant Gene sets before FDR:",sum(res$p.value<0.05),"\n")
  cat("Significant Gene sets after FDR:",sum(res$fdr<fdr_th),"\n")
  
  rownames(res) = NULL
  res = res[order(res$fdr),]
  
  return(res)
}

# -----------------------------
# FUNCTION: MEGA.core
# -----------------------------
# Description: MEGA was developed to identify predefined gene 
# sets (e.g. genes involved in the same pathway, or predisposing 
# to specific diseases) that show a significantly higher number 
# of mutations in a group of samples as compared to another 
# group of samples.
#
# Input:
# A and B: Boolean matrices of mutations. Coloums are samples, while rows are
# mutations. The first coloumn must always contain the name of thegene in which
# the mutation fall.
#
# Example:
# +-------------------------+
# Symbol  S1    S2    S3
# MAST2   TRUE  FALSE FALSE
# MAST2   TRUE  FALSE FALSE
# ABCA10  FALSE FALSE FALSE
# SLC4A9  TRUE  TRUE  FALSE
# ZNF572  FALSE FALSE FALSE
# ASPM    FALSE FALSE TRUE
# CSMD2   FALSE FALSE TRUE
# +-------------------------+
#
# X --> Gene set: a list of genes
# Example: X = c("ABCA10","ZNF572","ASPM","CSMD2")
# 
MEGA.core = function(A,B,X,test="W") {
  ix = A[,1] %in% X
  if (sum(ix) > 0) {
    A = A[ix,2:ncol(A)]
    Da = apply(A,2,sum)
  } else {
    Da = matrix(data=0, nrow = 1, ncol = ncol(A))
  }
  
  ix = B[,1] %in% X
  if (sum(ix)) {
    B = B[ix,2:ncol(B)]
    Db = apply(B,2,sum)
  } else {
    Db = matrix(data=0, nrow = 1, ncol = ncol(B))
  }
  
  if (test == "W")
  {
    p = wilcox.test(Da,Db,alternative = "greater",exact = F)$p.value
  } else {
    p = ks.test(Da,Db, alternative = "less", exact = F)$p.value
  }
  return(p)
}

# -----------------------------
# FUNCTION: MEGA.bootstrapping
# -----------------------------
# Description: At each iteration, the larger cohort is randomly downsampled to reach the 
# sample size of the smaller cohort. The procedure is repeted for nsim time times and 
# MEGA is executed.
MEGA.bootstrapping = function(A,B,gene.sets,nsim,test="W") {
  gA = A[,1]
  gB = B[,1]
  A = A[,-1]
  B = B[,-1]
   
  pb = txtProgressBar(min = 0, max = nsim, initial = 0,style=3)
  out = sapply(1:nsim, function(x,a=A,b=B,gs=gene.sets,ga=gA,gb=gB,z=pb,t=test) {
                o = NULL
                N = min(ncol(a),ncol(b))
                setTxtProgressBar(z,x)
                if (ncol(a)<ncol(b)) {
                  ix = sample(1:ncol(b),N)
                  o = sapply(gs, function(xx,aa=cbind(ga,a),bb=cbind(gb,b[,ix]),tt=t) MEGA.core(aa,bb,xx,tt))
                } else {
                  ix = sample(1:ncol(a),N)
                  o = sapply(gs, function(xx,aa=cbind(ga,a[,ix]),bb=cbind(gb,b),tt=t) MEGA.core(aa,bb,xx,tt))
                }
                return(o)
              }
        )
}

# -----------------------------
# FUNCTION: print.logo
# -----------------------------
# Description: print MEGA logo
print.logo = function() {
  cat(".___  ___.  _______   _______      ___      \n")
  cat("|   \\/   | |   ____| /  _____|    /   \\     \n")
  cat("|  \\  /  | |  |__   |  |  __     /  ^  \\    \n")
  cat("|  |\\/|  | |   __|  |  | |_ |   /  /_\\  \\   \n")
  cat("|  |  |  | |  |____ |  |__| |  /  _____  \\  \n")
  cat("|__|  |__| |_______| \\______| /__/     \\__\\ \n")
}


read.gmt.file = function(pathMsigDbFile) {		
  inputFile <- pathMsigDbFile		
  con  <- file(inputFile, open = "r")		
  
  c = 1		
  pathway.list <- vector(mode="list",length=0)		
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0)
  {		
    myVector <- do.call("rbind",strsplit(oneLine, "\t"))		
    t = vector(mode="list",length=1)		
    t[[1]] = myVector[3:length(myVector)]		
    names(t) = myVector[1]		
    pathway.list = c(pathway.list,t)		
    c = c+1		
  }		
  
  close(con)		
  return(pathway.list)		
}

MEGA.MC.core = function(A,X)
{
  ix = A[,1] %in% X
  if (sum(ix) > 0) {
    A = A[ix,2:ncol(A)]
    Da = apply(A,2,sum)
  } else {
    Da = matrix(data=0, nrow = 1, ncol = ncol(A))
  }
  return(Da)
}

shuffle.patient.mutations = function(A,gene.cds.length)
{
  csum <- colSums(A[,2:ncol(A)]*1)
  random.mut.genes <- sapply(csum, function(x,y=gene.cds.length) sample(y$SYMBOL,size = x,replace = T,prob = y$prob))
  A.random <- data.frame(symbol=unique(unlist(random.mut.genes)),stringsAsFactors = F)
  A.random <- cbind(A.random,sapply(random.mut.genes, function(x,y=A.random$symbol) y %in% x)*1)
  return(A.random)
}

MEGA.MC = function(A,gene.sets,gene.cds.length,th=0.05,nsim=1000,cores=2)
{
  print("Performing Monte Carlo Permutations...")
  print(paste("Number of Cores using:",cores))
  
  if (cores > 1)
  {
    cl = makeCluster(cores)
    clusterExport(cl, c("A", "gene.cds.length","MEGA.MC.core","shuffle.patient.mutations"),envir=environment())
  }
  
  
  Da.empirical.matrix <- matrix(data = 0,nrow = length(gene.sets),ncol = nsim)
  rownames(Da.empirical.matrix) <- names(gene.sets)
  for (i in 1:nsim)
  {
    A_rand <- shuffle.patient.mutations(A,gene.cds.length)
    if (cores > 1) 
    {
      clusterExport(cl,"A_rand",envir=environment())
      Da.empirical.matrix[,i] <- parSapply(cl,gene.sets,function(x) sum(MEGA.MC.core(A_rand,x)))
    } else {
      Da.empirical.matrix[,i] <- sapply(gene.sets, function(x,y=A_rand) sum(MEGA.MC.core(y,x)))
    }
  }
  p.values <- sapply(1:length(gene.sets), function(x,y=Da.empirical.matrix,z=A,g=gene.sets) sum(y[x,] > sum(MEGA.MC.core(z,g[[x]]))))/nsim
  res <- data.frame(pathway=names(gene.sets),empirical.pvalue=p.values,stringsAsFactors = F)
  res <- res[order(res$empirical.pvalue),]

  if (cores > 1)
      stopCluster(cl)
  
  colnames(res) <- c("pathway","p.value")
  
  return(res)
}