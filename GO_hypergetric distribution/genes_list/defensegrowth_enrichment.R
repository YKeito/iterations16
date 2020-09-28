NumMCL <- na.omit(read.table("~/Nakano_RNAseq/network_analysis/cytoscape/NumMCL.txt", sep = "\t", header = T))
rownames(NumMCL) <- NumMCL$name
defensegrowth_list <- read.table("~/Nakano_RNAseq/network_analysis/base/defensegrowth_FDR005union.txt", sep = "\t", row.names = 1, header = T)
####defense####
defense_list <- defensegrowth_list[defensegrowth_list$Defense_gene == "D", ]
#全クラスターのdefenseのDEGs数
Numdefense_list <- sum(!is.na(match(rownames(defense_list), NumMCL$name)))
#全クラスターのdefense以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-Numdefense_list)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
defense_enrichment_pvalue <- c()
check <- c()
defense_check_MCL <- list()
 for(n in n:maxNumMCL){
        test <- NumMCL[NumMCL$X__mclCluster == n, ]
        NumCY_DEGs <- length(intersect(rownames(defense_list), rownames(test)))
        Numsubcluster <- nrow(test)
        defense_enrichment_pvalue <- c(defense_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), Numdefense_list, Numothergenes, Numsubcluster, lower.tail = FALSE))
        names(defense_enrichment_pvalue)[n] <- n
        check <- c(check, NumCY_DEGs)
          if(NumCY_DEGs != 0){
              T_INT <- intersect(rownames(defense_list), rownames(test))
              names(T_INT) <- rep(n, times = length(T_INT))
              defense_check_MCL <- c(defense_check_MCL, list(T_INT))
             }
          print(maxNumMCL - n)
          n <- n + 1
        }
defense_enrichment_qvalue <- p.adjust(defense_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
defense_enriched_Node <- defense_check_MCL[defense_enrichment_qvalue < 0.05]
Numdefense_subcluster <- c()
i <- 1
for(i in i:length(defense_enriched_Node)){
    Numdefense_subcluster <- c(Numdefense_subcluster, names(defense_enriched_Node[[i]])[1])
    i <- i+1
    }
Numdefense_subcluster <- as.numeric(Numdefense_subcluster)
#names(CY20_enriched_Node[[11]])[1]でクラスター番号とってこれる。

defense_enriched_Node <- data.frame(unlist(defense_enriched_Node))
write.table(defense_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/defense_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
save(defense_enrichment_qvalue, file = "~/Nakano_RNAseq/network_analysis/.RData/defense_enrichment_qvalue.RData")


####growth####
growth_list <- defensegrowth_list[defensegrowth_list$Growth_gene == "G", ]
#全クラスターのgrowthのDEGs数
Numgrowth_list <- sum(!is.na(match(rownames(growth_list), NumMCL$name)))
#全クラスターのgrowth以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-Numgrowth_list)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
growth_enrichment_pvalue <- c()
check <- c()
growth_check_MCL <- list()
for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(growth_list), rownames(test)))
  Numsubcluster <- nrow(test)
  growth_enrichment_pvalue <- c(growth_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), Numgrowth_list, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(growth_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(growth_list), rownames(test))
    names(T_INT) <- rep(n, times = length(T_INT))
    growth_check_MCL <- c(growth_check_MCL, list(T_INT))
  }
  print(maxNumMCL - n)
  n <- n + 1
}
growth_enrichment_qvalue <- p.adjust(growth_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
growth_enriched_Node <- growth_check_MCL[growth_enrichment_qvalue < 0.05]
Numgrowth_subcluster <- c()
i <- 1
for(i in i:length(growth_enriched_Node)){
  Numgrowth_subcluster <- c(Numgrowth_subcluster, names(growth_enriched_Node[[i]])[1])
  i <- i+1
}
Numgrowth_subcluster <- as.numeric(Numgrowth_subcluster)
#names(CY20_enriched_Node[[11]])[1]でクラスター番号とってこれる。

growth_enriched_Node <- data.frame(unlist(growth_enriched_Node))
write.table(growth_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/growth_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
save(growth_enrichment_qvalue, file = "~/Nakano_RNAseq/network_analysis/.RData/growth_enrichment_qvalue.RData")