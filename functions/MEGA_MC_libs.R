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
  
  res <- NULL
  for (i in 1:length(gene.sets))
  {
    print(i)
    
    X <- gene.sets[[i]]
    Da <- sum(MEGA.MC.core(A,X))
    if (Da > 0)
    {
      if (cores > 1) 
      {
        clusterExport(cl,"X",envir=environment())
        Da.empirical.distributions <- parSapply(cl,1:nsim,function(x) sum(MEGA.MC.core(shuffle.patient.mutations(A,gene.cds.length),X)))
      } else {
        Da.empirical.distributions <- sapply(1:nsim, function(x,y=A,z=X,k=gene.cds.length) sum(MEGA.MC.core(shuffle.patient.mutations(y,k),z)))
      }
      res <- rbind(res,data.frame(pathway=names(gene.sets)[i],empirical.pvalue=sum(Da.empirical.distributions>=Da)/nsim,stringsAsFactors = F))
    } else {
      res <- rbind(res,data.frame(pathway=names(gene.sets)[i],empirical.pvalue=1,stringsAsFactors = F))
    }
    #print(res)
  }
  
  if (cores > 1)
    stopCluster(cl)
  
  res$FDR <- p.adjust(res$empirical.pvalue,method = "fdr")
  res <- res[order(res$FDR),]
  return(res)
}

