NumMCL <- na.omit(read.table("~/Nakano_RNAseq/network_analysis/CWRtoscape/NumMCL.txt", sep = "\t", header = T))
rownames(NumMCL) <- NumMCL$name
####CWR enrichment####
CWR <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/CWR.txt", sep = "\t", header = T)
#全クラスターのCWRのDEGs数
NumCWR <- sum(!is.na(match(CWR$Gene.ID, NumMCL$name)))
#全クラスターのCWR以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumCWR)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
CWR_enrichment_pvalue <- c()
check <- c()
CWR_check_MCL <- list()
subnodes <- list()
for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCWR_DEGs <- length(intersect(CWR$Gene.ID, rownames(test)))
  Numsubcluster <- nrow(test)
  CWR_enrichment_pvalue <- c(CWR_enrichment_pvalue, phyper(c(NumCWR_DEGs - 1), NumCWR, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CWR_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumCWR_DEGs)
  T_INT <- intersect(CWR$Gene.ID, rownames(test))
  names(T_INT) <- rep(n, times = length(T_INT))
  CWR_check_MCL <- c(CWR_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(CWR_check_MCL)){
  data <- unlist(CWR_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}
NumCWR_subcluster <- as.numeric(NumCWR_subcluster)

NumSubCluster <- c()
i <- 1
for(i in i:maxNumMCL){
  NumSubCluster <- c(NumSubCluster, paste0("SubCluster", i))
  i <- i+1
}
CWR_statistics <- data.frame(NumSubcluster = NumSubCluster,
                              NumCWR_DEGs = check, 
                              CWR_DEGs_AGI = unlist(AGI),
                              NumSubcluster = unlist(subnodes), 
                              p_value = CWR_enrichment_pvalue, 
                              enrichment_score = -log2(CWR_enrichment_pvalue)
)

####output####
write.table(CWR_statistics, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/genes_set/CWR_statistics.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
