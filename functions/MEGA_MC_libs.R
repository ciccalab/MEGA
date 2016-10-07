require(parallel)

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
  if (cores > 1)
  {
    cl = makeCluster(cores)
    clusterExport(cl, c("A", "gene.cds.length","MEGA.MC.core","shuffle.patient.mutations"),envir=environment())
  }
  
  withProgress(
  message = 'Monte Carlo Simulations', value = 0,
  {
  Da.empirical.matrix <- matrix(data = 0,nrow = length(gene.sets),ncol = nsim)
  rownames(Da.empirical.matrix) <- names(gene.sets)
  for (i in 1:nsim)
  {
    #incProgress(1/(nsim+2), detail = paste(i, "out of",nsim))
    
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
  
  res$MC.FDR <- p.adjust(res$empirical.pvalue,method = "fdr")
  
  })
  
  return(res)
}

