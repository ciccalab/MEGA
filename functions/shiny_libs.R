require(MASS)
require(parallel)

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
MEGA = function(A,B,gene.sets,fdr_th=0.1,bootstrapping=F,nsim=1000, test="W",montecarlo=F,gene.cds.length=NULL,cpus=2) {
  
  # Execute the enrichement
  # ------------------------
  withProgress(message = 'Running MEGA', value = 0, {
    #inc <- 1/(length(gene.sets)+1)
    p = sapply(1:length(gene.sets), function(x,z=inc,a=A,b=B,gs=gene.sets,t=test) {
      incProgress(1/(length(gs)+1), detail = paste(x,"out of", length(gs)))
      r = MEGA.core(A,B,gs[[x]],t)
      names(r) = names(gs)[x]
      return(r)
    })
    
    res = data.frame("gene.set"=names(p),p.value=p,fdr= p.adjust(p,method = "fdr"),stringsAsFactors = F)
    setProgress(1)
  })
  
  
  # Execute Bootstrapping if required
  # ----------------------------------
  if (bootstrapping) {
    if (sum(res$fdr<fdr_th)) {
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
  
  # Execute Monte Carlo if required
  # ----------------------------------
  if (montecarlo) {
    res.mc = MEGA.MC(A,gene.sets,gene.cds.length,nsim,cpus)
    res$MC.pvalue = NA
    id1 = match(res$gene.set,res.mc$pathway)
    id2 <- id1[!is.na(id1)]
    res$MC.pvalue[!is.na(id1)] = res.mc$empirical.pvalue[id2]
  }
  
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
MEGA.core = function(A,B,X,test=2,montecarlo=F,cores=2,nsim=1000) {
  p = 1
  ix = A[,1] %in% X
  
  if (sum(ix) > 0)
  {
    A = A[ix,2:ncol(A)]
    Da = apply(A,2,sum)
  } else {
    Da = matrix(data=0, nrow = 1, ncol = ncol(A))
  }
  
  ix = B[,1] %in% X
  if (sum(ix))
  {
    B = B[ix,2:ncol(B)]
    Db = apply(B,2,sum)
  } else {
    Db = matrix(data=0, nrow = 1, ncol = ncol(B))
  }
  
  p = do.test(Da,Db,test)
  
  if (montecarlo)
  {
    if (cores > 1)
    {
      cl = makeCluster(cores)
      clusterExport(cl, c("Da", "Db","test","do.test","montecarlo.core"),envir=environment())
      p.empirical = parSapply(cl,1:nsim,function(x) montecarlo.core(Da,Db,test))
      stopCluster(cl)
    } else {
      p.empirical = sapply(1:nsim, function(x,a=Da,b=Db,t=test) montecarlo.core(a,b,t))
    }
    pe = sum(p.empirical<=p)/nsim
    return(pe)
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
  
  inc <- 1/nsim
  withProgress(message = 'Starting the Bootstrapping process...', value = 0, {
    out = sapply(1:nsim, function(x,a=A,b=B,gs=gene.sets,ga=gA,gb=gB,z=inc,tot=nsim,t=test) {
                  o = NULL
                  N = min(ncol(a),ncol(b))
                  incProgress(z,message = paste("random sampling",x,"of",nsim))
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
  })
  
  return(out)
}

get.gene.set = function(ix) {
  c("./gene_sets/c2.cp.kegg.v5.1.symbols.gmt","./gene_sets/c5.bp.v5.1.symbols.gmt","./gene_sets/c5.cc.v5.1.symbols.gmt","./gene_sets/c5.mf.v5.1.symbols.gmt","./gene_sets/c2.cp.reactome.v5.1.symbols.gmt","./gene_sets/c2.cp.biocarta.v5.1.symbols.gmt","./gene_sets/c1.positional.v5.1.symbols.gmt","./gene_sets/h.all.v5.1.symbols.gmt","./gene_sets/c3.tft.v5.1.symbols.gmt","./gene_sets/c4.cgn.v5.1.symbols.gmt","./gene_sets/c4.cm.v5.1.symbols.gmt","./gene_sets/c6.all.v5.1.symbols.gmt","./gene_sets/ncomm.cereda.gwas.symbols.gmt")[ix]
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


do.test = function(Da,Db,ix)
{
  p = 1
  if (sum(Da)>0 | sum(Db) > 0)
  {
    if (ix == 2)
    {
      p = wilcox.test(Da,Db,alternative = "greater",exact = F)$p.value
    }
    if (ix ==1){
      options(warn=-1)
      p = anova(glm.nb(c(Da,Db) ~ c(rep("t",length(Da)),rep("c",length(Db)))))$`Pr(>Chi)`[2]
      options(warn=0)
    }
    if (ix == 3)
    {
      p = ks.test(Da,Db, alternative = 'greater', exact = F)$p.value
    }
  }
  return(p)
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

MEGA.MC = function(A,gene.sets,gene.cds.length,nsim=1000,cores=2)
{
  if (cores > 1)
  {
    cl = makeCluster(cores)
    clusterExport(cl, c("A", "gene.cds.length","MEGA.MC.core","shuffle.patient.mutations"),envir=environment())
  }
  
  withProgress(
    message = 'MC Simulations', value = 0,
    {
      Da.empirical.matrix <- matrix(data = 0,nrow = length(gene.sets),ncol = nsim)
      rownames(Da.empirical.matrix) <- names(gene.sets)
      for (i in 1:nsim)
      {
        incProgress(1/(nsim+2), detail = paste(i, "out of",nsim))
        
        A_rand <- shuffle.patient.mutations(A,gene.cds.length)
        if (cores > 1) 
        {
          clusterExport(cl,"A_rand",envir=environment())
          Da.empirical.matrix[,i] <- parSapply(cl,gene.sets,function(x) sum(MEGA.MC.core(A_rand,x)))
        } else {
          Da.empirical.matrix[,i] <- sapply(gene.sets, function(x,y=A_rand) sum(MEGA.MC.core(y,x)))
        }
      }
      
      incProgress(1/(nsim+1), detail = "merge results")
      p.values <- sapply(1:length(gene.sets), function(x,y=Da.empirical.matrix,z=A,g=gene.sets) sum(y[x,] > sum(MEGA.MC.core(z,g[[x]]))))/nsim
      res <- data.frame(pathway=names(gene.sets),empirical.pvalue=p.values,stringsAsFactors = F)
      res <- res[order(res$empirical.pvalue),]
      if (cores > 1)
        stopCluster(cl)
      
      #res$MC.FDR <- p.adjust(res$empirical.pvalue,method = "fdr")
      
    })
  
  return(res)
}

