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
MEGA = function(A,B,gene.sets,fdr_th=0.1,bootstrapping=F,nsim=1000) {
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
  cat("========================================\n")
  
  cat("\nStep 1: Enrichement Gene Set Enrichement Anlysis\n")
  
  # Execute the enrichement
  # ------------------------
  pb = txtProgressBar(min = 0, max = length(gene.sets), initial = 0,style=3)
  p = sapply(1:length(gene.sets), function(x,z=pb,a=A,b=B,gs=gene.sets) {
    setTxtProgressBar(z,x)
    r = MEGA.core(A,B,gs[[x]])
    names(r) = names(gs)[x]
    return(r)
  })
  
  res = data.frame("gene.set"=names(p),p.value=p,fdr= p.adjust(p,method = "fdr"),stringsAsFactors = F)
  rownames(res) = NULL
  res = res[order(res$fdr),]
  
  # Execute Bootstrapping if required
  # ----------------------------------
  if (bootstrapping) {
    if (sum(res$fdr<fdr_th)) {
      cat("\n\nStep 2: Bootstrapping for",sum(res$fdr<fdr_th),"significant gene sets\n")
      bs = MEGA.bootstrapping(A,B,gene.sets[res$gene.set[res$fdr<fdr_th]],nsim)
      res$success_percentage = NA
      res$success_percentage[res$fdr<fdr_th] = apply(bs, 1, function(x,n=nsim) sum(x<0.05)*100/n)
    } else {
      cat("\n\n Step 2: Bootstrapping can not be performed. No significant gene sets\n")
    }
  } 
  
  cat("\n\nResults:\n")
  cat("Significant Gene sets before FDR:",sum(res$p.value<0.05),"\n")
  cat("Significant Gene sets after FDR:",sum(res$fdr<fdr_th),"\n")
  
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
MEGA.core = function(A,B,X) {
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
  
  p = wilcox.test(Da,Db,alternative = "greater",exact = F)$p.value
  return(p)
}

# -----------------------------
# FUNCTION: MEGA.bootstrapping
# -----------------------------
# Description: At each iteration, the larger cohort is randomly downsampled to reach the 
# sample size of the smaller cohort. The procedure is repeted for nsim time times and 
# MEGA is executed.
MEGA.bootstrapping = function(A,B,gene.sets,nsim) {
  gA = A[,1]
  gB = B[,1]
  A = A[,-1]
  B = B[,-1]
   
  pb = txtProgressBar(min = 0, max = nsim, initial = 0,style=3)
  out = sapply(1:nsim, function(x,a=A,b=B,gs=gene.sets,ga=gA,gb=gB,z=pb) {
                o = NULL
                N = min(ncol(a),ncol(b))
                setTxtProgressBar(z,x)
                if (ncol(a)<ncol(b)) {
                  ix = sample(1:ncol(b),N)
                  o = sapply(gs, function(xx,aa=cbind(ga,a),bb=cbind(gb,b[,ix])) MEGA.core(aa,bb,xx))
                } else {
                  ix = sample(1:ncol(a),N)
                  o = sapply(gs, function(xx,aa=cbind(ga,a[,ix]),bb=cbind(gb,b)) MEGA.core(aa,bb,xx))
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

