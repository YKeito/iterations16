iterations14 <- na.omit(read.table("~/Nakano_RNAseq/network_analysis/cytoscape/iteration14_node.csv", sep = ",", header = T))
rownames(iterations14) <- iterations14$name
####TF_family enrichment####
TF_family <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/arabidopsis_TF_family.txt", sep = "\t", header = T)
maxiterations14 <- max(iterations14$X__mclCluster)
maxTF <- length(unique(TF_family$TF))
i <- 1
temp <- c()
all_TF <- as.character(unique(TF_family$TF))
for(i in i:maxTF){
  temp <- TF_family[grep(unique(TF_family$TF)[i], TF_family$TF), ]$AGI
  #全クラスターのTFのDEGs数
  NumTF <- length(intersect(temp, iterations14$name))
  #全クラスターのTF以外のDEGs数
  Numothergenes <- c(length(iterations14$name)-NumTF)
  n <- 1
  TF_enrichment_pvalue <- c()
  check <- c()
  TF_check_MCL <- list()
  subnodes <- list()
  for(n in n:maxiterations14){
    test <- iterations14[iterations14$X__mclCluster == n, ]
    NumTF_DEGs <- length(intersect(TF_family$AGI, test$name))
    Numsubcluster <- nrow(test)
    TF_enrichment_pvalue <- c(TF_enrichment_pvalue, phyper(c(NumTF_DEGs - 1), NumTF, Numothergenes, Numsubcluster, lower.tail = FALSE))
    names(TF_enrichment_pvalue)[n] <- n
    subnodes <- c(subnodes, list(Numsubcluster))
    check <- c(check, NumTF_DEGs)
    T_INT <- intersect(rownames(test), iterations14$name)
    names(T_INT) <- rep(n, times = length(T_INT))
    TF_check_MCL <- c(TF_check_MCL, list(T_INT))
    n <- n + 1
  }
  print(maxTF - i)
  i <- i+1
}


#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(TF_check_MCL)){
  data <- unlist(TF_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}
NumTF_subcluster <- as.numeric(NumTF_subcluster)

NumSubCluster <- c()
i <- 1
for(i in i:maxiterations14){
  NumSubCluster <- c(NumSubCluster, paste0("SubCluster", i))
  i <- i+1
}
TF_statistics <- data.frame(NumSubcluster = NumSubCluster,
                             NumTF_DEGs = check, 
                             TF_DEGs_AGI = unlist(AGI),
                             NumSubcluster = unlist(subnodes), 
                             p_value = TF_enrichment_pvalue, 
                             enrichment_score = -log2(TF_enrichment_pvalue)
)

####output####
write.table(TF_statistics, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/genes_set/TF_statistics.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
