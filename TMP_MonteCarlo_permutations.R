
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

pl = make_pathway_list("mSigDB_gene_sets/c2.cp.kegg.v5.1.symbols.gmt")

tmp= pl["KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"]

syn = read.table("example_dataset/A.tsv.gz", h=T)
con = read.table("example_dataset/B.tsv.gz", h=T)


gene=syn[,1]
syn=syn[,2:ncol(syn)]; rownames(syn) = gene
gene=con[,1]
con=con[,2:ncol(con)]; rownames(con) = gene

syn = syn[which(rownames(syn)%in%tmp[[1]]),]
con = con[which(rownames(con)%in%tmp[[1]]),]
dim(syn)
dim(con)

a = apply(syn,2,sum)
b = apply(con,2,sum)

tot=c(a,b)
p0=wilcox.test(a,b,alternative = "greater",exact = F)$p.value

set.seed(1)
p = rep(NA, 1000)
for(i in 1:1000){
 x = sample(tot, length(a), replace=T )
 y = tot[!names(tot)%in%names(x)]
 p[i] = wilcox.test(x,y, alternative = 'g', exact = F)$p.value
}
pe = sum(p<p0)/1000


p0=ks.test(a,b,alternative = "less",exact = F)$p.value

set.seed(1)
p = rep(NA, 1000)
for(i in 1:1000){
  x = sample(tot, length(a), replace=T )
  y = tot[!names(tot)%in%names(x)]
  p[i] = ks.test(x,y, alternative = 'less', exact = F)$p.value
}
pe = sum(p<=p0)/1000


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


  } else {
  }
  return(c(p0,pe))
}

res = lapply(pl, MEGA.core, A=syn, B=con, test="W")

A=syn
B=con
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
