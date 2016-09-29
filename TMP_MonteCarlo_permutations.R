make_pathway_list = function(pathMsigDbFile) {
  inputFile <- pathMsigDbFile
  con  <- file(inputFile, open = "r")

  c = 1
  pathway.list <- vector(mode="list",length=0)
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    print(c)
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
    p0 = wilcox.test(Da,Db,alternative = "greater",exact = F)$p.value
    tot=c(Da,Db)
    set.seed(1)
    p = rep(NA, 1000)
    for(i in 1:1000){
      x = sample(tot, length(Da), replace=T )
      y = tot[!names(tot)%in%names(x)]
      p[i] = wilcox.test(x,y, alternative = 'g', exact = F)$p.value
    }
    pe = sum(p<=p0)/1000

  }else if (test == "K")
  {
    p0 = ks.test(Da,Db,alternative = "less",exact = F)$p.value
    tot=c(Da,Db)
    set.seed(1)
    p = rep(NA, 1000)
    for(i in 1:1000){
      x = sample(tot, length(Da), replace=T )
      y = tot[!names(tot)%in%names(x)]
      p[i] = ks.test(x,y, alternative = 'l', exact = F)$p.value
    }
    pe = sum(p<=p0)/1000


  } else {
  }
  return(c(p0,pe))
}

pl = make_pathway_list("mSigDB_gene_sets/c2.cp.kegg.v5.1.symbols.gmt")


syn = read.table("example_dataset/A.tsv.gz", h=T)
con = read.table("example_dataset/B.tsv.gz", h=T)

A=syn
B=con

# ESEMPIO WILCOXON ============================
test="W"
pb = txtProgressBar(min = 0, max = length(pl), initial = 0,style=3)
p = lapply(1:length(pl), function(x,z=pb,a=A,b=B,gs=pl,t=test) {
  setTxtProgressBar(z,x)
  r = MEGA.core(A,B,gs[[x]],t)
  names(r) = names(gs)[x]
  return(r)
})
df = do.call(rbind, p)
colnames(df) = c("PV","PE")
rownames(df) = names(pl)
df = as.data.frame(df)

df$FDR=p.adjust(df$PV, "fdr")
subset(df, PV<0.05 & PE<0.05)

# ESEMPIO kolmogorov ============================

test="K"
pb = txtProgressBar(min = 0, max = length(pl), initial = 0,style=3)
p = lapply(1:length(pl), function(x,z=pb,a=A,b=B,gs=pl,t=test) {
  setTxtProgressBar(z,x)
  r = MEGA.core(A,B,gs[[x]],t)
  names(r) = names(gs)[x]
  return(r)
})
dfk = do.call(rbind, p)
colnames(dfk) = c("PV","PE")
rownames(dfk) = names(pl)
dfk = as.data.frame(dfk)

dfk$FDR=p.adjust(dfk$PV, "fdr")
dfk = dfk[order(dfk[,1]),]
subset(dfk, PV<0.05 & PE<0.05)
