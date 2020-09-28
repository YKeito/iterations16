NumMCL <- na.omit(read.table("~/Nakano_RNAseq/network_analysis/cytoscape/NumMCL.txt", sep = "\t", header = T))
rownames(NumMCL) <- NumMCL$name
####BTH enrichment####
BTH <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/BTH_affected_AGI.txt", sep = "\t", header = T)
#全クラスターのBTHのDEGs数
NumBTH <- sum(!is.na(match(BTH$AGI, NumMCL$name)))
#全クラスターのBTH以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumBTH)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
BTH_enrichment_pvalue <- c()
check <- c()
BTH_check_MCL <- list()
subnodes <- list()
for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumBTH_DEGs <- length(intersect(BTH$AGI, rownames(test)))
  Numsubcluster <- nrow(test)
  BTH_enrichment_pvalue <- c(BTH_enrichment_pvalue, phyper(c(NumBTH_DEGs - 1), NumBTH, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(BTH_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumBTH_DEGs)
  T_INT <- intersect(BTH$AGI, rownames(test))
  names(T_INT) <- rep(n, times = length(T_INT))
  BTH_check_MCL <- c(BTH_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(BTH_check_MCL)){
  data <- unlist(BTH_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}
NumBTH_subcluster <- as.numeric(NumBTH_subcluster)

NumSubCluster <- c()
i <- 1
for(i in i:maxNumMCL){
  NumSubCluster <- c(NumSubCluster, paste0("SubCluster", i))
  i <- i+1
}
BTH_statistics <- data.frame(NumSubcluster = NumSubCluster,
                             NumBTH_DEGs = check, 
                             BTH_DEGs_AGI = unlist(AGI),
                             NumSubcluster = unlist(subnodes), 
                             p_value = BTH_enrichment_pvalue, 
                             enrichment_score = -log2(BTH_enrichment_pvalue)
)

####output####
write.table(BTH_statistics, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/genes_set/BTH_statistics.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
