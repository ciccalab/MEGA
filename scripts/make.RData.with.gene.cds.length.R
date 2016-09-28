get.cds.length = function(TxDb)
{
  cds.loc.all <- cdsBy(TxDb, by="gene")
  entrezids <- names(cds.loc.all)
  # filter out entrez ids with no gene symbols
  map.entrez.symbol <- subset(select(Homo.sapiens, keys=entrezids, columns=c("SYMBOL"), keytype="ENTREZID"),!is.na(SYMBOL))
  cds.loc.all <- cds.loc.all[map.entrez.symbol$ENTREZID]
  
  # get cds length for the remaning genes (entrez ids)
  tmp <- do.call("rbind",lapply(cds.loc.all, function(x) data.frame(cds_id = x@elementMetadata$cds_id[which.max(x@ranges@width)], width = x@ranges@width[which.max(x@ranges@width)],stringsAsFactors = F)))
  
  # make final data frame and return it
  gene.cds.length <- cbind(map.entrez.symbol,tmp)
  gene.cds.length$prob <- gene.cds.length$width / sum(gene.cds.length$width)
  return(gene.cds.length)
}

# HG19
require(Homo.sapiens)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
wd <- "~/MEGA/"
setwd(wd)
gene.cds.length <- get.cds.length(TxDb.Hsapiens.UCSC.hg19.knownGene)
save(gene.cds.length,file = "./RData/HG19.gene.cds.length.RData")

# HG38
require(Homo.sapiens)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
wd <- "~/MEGA/"
setwd(wd)
gene.cds.length <- get.cds.length(TxDb.Hsapiens.UCSC.hg38.knownGene)
save(gene.cds.length,file = "./RData/HG38.gene.cds.length.RData")
