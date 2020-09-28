defensegrowth_list <- read.table("~/Nakano_RNAseq/network_analysis/base/defensegrowth_FDR005union.txt", sep = "\t", row.names = 1, header = T)
####defense####
defense_list <- defensegrowth_list[defensegrowth_list$Defense_gene == "D", ]
####defense enrichment####
#全クラスターのdefenseのDEGs数
Numdefense_list <- sum(!is.na(match(rownames(defense_list), NumMCL$name)))
#全クラスターのdefense以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-Numdefense_list)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
defense_enrichment_pvalue <- c()
check <- c()
defense_check_MCL <- list()
subnodes <- list()
for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  Numdefense_DEGs <- length(intersect(rownames(defense_list), rownames(test)))
  Numsubcluster <- nrow(test)
  defense_enrichment_pvalue <- c(defense_enrichment_pvalue, phyper(c(Numdefense_DEGs - 1), Numdefense_list, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(defense_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, Numdefense_DEGs)
  T_INT <- intersect(rownames(defense_list), rownames(test))
  names(T_INT) <- rep(n, times = length(T_INT))
  defense_check_MCL <- c(defense_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(defense_check_MCL)){
  data <- unlist(defense_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}
Numdefense_subcluster <- as.numeric(Numdefense_subcluster)

defense_statistics <- data.frame(NumSubcluster = NumSubCluster, 
                                 Numdefense_DEGs = check,
                                 defense_DEGs_AGI = unlist(AGI),
                                 NumSubcluster = unlist(subnodes),
                                 p_value = defense_enrichment_pvalue, 
                                 enrichment_score = -log2(defense_enrichment_pvalue)
                                 )
####growth####
growth_list <- defensegrowth_list[defensegrowth_list$Growth_gene == "G", ]
####growth enrichment####
#全クラスターのgrowthのDEGs数
Numgrowth_list <- sum(!is.na(match(rownames(growth_list), NumMCL$name)))
#全クラスターのgrowth以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-Numgrowth_list)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
growth_enrichment_pvalue <- c()
check <- c()
growth_check_MCL <- list()
subnodes <- list()
for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  Numgrowth_DEGs <- length(intersect(rownames(growth_list), rownames(test)))
  Numsubcluster <- nrow(test)
  growth_enrichment_pvalue <- c(growth_enrichment_pvalue, phyper(c(Numgrowth_DEGs - 1), Numgrowth_list, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(growth_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, Numgrowth_DEGs)
  T_INT <- intersect(rownames(growth_list), rownames(test))
  names(T_INT) <- rep(n, times = length(T_INT))
  growth_check_MCL <- c(growth_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(growth_check_MCL)){
  data <- unlist(growth_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}
Numgrowth_subcluster <- as.numeric(Numgrowth_subcluster)

growth_statistics <- data.frame(NumSubcluster = NumSubCluster, 
                                Numgrowth_DEGs = check,
                                growth_DEGs_AGI = unlist(AGI),
                                NumSubcluster = unlist(subnodes),
                                p_value = growth_enrichment_pvalue, 
                                enrichment_score = -log2(growth_enrichment_pvalue)
                                )