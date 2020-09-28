####nottimecource####
#unique(EnrichedGO_Summ_CY15$CluNum)
#subcluster1
subcluster1_top10 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 1, ]
subcluster1_top10 <- subcluster1_top10[order(subcluster1_top10$enrichment_score, decreasing = T), ]
subcluster1_top10 <- subcluster1_top10[1:10, ]
library(ggplot2)
g <- ggplot(
  subcluster1_top10,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(FDR)")
g <- g + xlab("")
g <- g + labs(title = "subcluster1 GO")
g <- g + theme(axis.text=element_text(size=60), axis.title=element_text(size=60,face="bold"))
g <- g + theme_bw(base_size = 60)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster/subcluster1_top10_GO.png", plot = g, dpi = 100, width = 21.6, height = 12)
#subcluster2
subcluster2_top10 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 2, ]
subcluster2_top10 <- subcluster2_top10[order(subcluster2_top10$enrichment_score, decreasing = T), ]
subcluster2_top10 <- subcluster2_top10[1:10, ]
library(ggplot2)
g <- ggplot(
  subcluster2_top10,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(FDR)")
g <- g + xlab("")
g <- g + labs(title = "subcluster2 GO")
g <- g + theme(axis.text=element_text(size=60), axis.title=element_text(size=60,face="bold"))
g <- g + theme_bw(base_size = 60)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster/subcluster2_top10_GO.png", plot = g, dpi = 100, width = 28.8, height = 20)
#subcluster7
subcluster7_top10 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 7, ]
subcluster7_top10 <- subcluster7_top10[order(subcluster7_top10$enrichment_score, decreasing = T), ]
subcluster7_top10 <- subcluster7_top10[1:10, ]
library(ggplot2)
g <- ggplot(
  subcluster7_top10,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(FDR)")
g <- g + xlab("")
g <- g + labs(title = "subcluster7 GO")
g <- g + theme(axis.text=element_text(size=60), axis.title=element_text(size=60,face="bold"))
g <- g + theme_bw(base_size = 60)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster/subcluster7_top10_GO.png", plot = g, dpi = 100, width = 25, height = 18)
#subcluster35
subcluster35_top10 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 35, ]
subcluster35_top10 <- subcluster35_top10[order(subcluster35_top10$enrichment_score, decreasing = T), ]
subcluster35_top10 <- subcluster35_top10[1:10, ]
library(ggplot2)
subcluster35_top10$GO_term <- c("nitrile biosynthetic process", 
                                "cytoplasmic stress granule", 
                                "response to molecule of fungal origin", 
                                "defense response signaling pathway", 
                                "glucosinolate catabolic process", 
                                "monovalent inorganic cation transport", 
                                "regulation of pH", 
                                "monovalent cation:proton antiporter activity", 
                                "solute:proton antiporter activity", 
                                "response to ozone")
g <- ggplot(
  subcluster35_top10,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(FDR)")
g <- g + xlab("")
g <- g + labs(title = "subcluster35 GO")
g <- g + theme(axis.text=element_text(size=60), axis.title=element_text(size=60,face="bold"))
g <- g + theme_bw(base_size = 60)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster/subcluster35_top10_GO.png", plot = g, dpi = 100, width = 25, height = 15)
#####GO termに含まれるAGIの発現値をまとめる。####
#subcluster1
#FDR0.005_CYunion
n <- 1
subcluster1_Node <- c()
subcluster1_symbol <- c()
subcluster1_expression <- c()
subcluster1_Nodes <- list()
test <- c()
temp <- c()
for(n in 1:10){
  subcluster1_Node <- unlist(strsplit(as.character(unlist(EnrichedGO_Summ_MCLFDR005$Node[n])), " "))
  subcluster1_Node <- subcluster1_Node[nchar(subcluster1_Node) == 9]
  temp <- as.character(unlist(allRNASeq$genesymbol[match(subcluster1_Node, rownames(allRNASeq))]))
  subcluster1_symbol <- c(subcluster1_symbol, temp)
  test <- allRNASeq[match(subcluster1_Node, rownames(allRNASeq)), ]
  test <- test[, 1:36]
  subcluster1_expression <- rbind(subcluster1_expression, test)
  subcluster1_Nodes <- c(subcluster1_Nodes, list(subcluster1_Node))
  subcluster1_Node <- c()
  n <- n+1
}
subcluster1_GOexp <- data.frame(subcluster1_AGI = as.character(unlist(subcluster1_Nodes)),
                                subcluster1_symbol = subcluster1_symbol,
                                subcluster1_expression)

#積み上げグラフ
#library(ggplot2)
#p <- ggplot(subcluster1_PieChart, aes(sample, NumCYNodes, fill = Category, order = desc(Category))) + geom_bar(stat = "identity", col="black") + theme_bw() + ggtitle("subcluster1 FDR0.005") 
#p <- p + geom_text(aes(label = NumCYNodes), size = 5, hjust = 0.5, vjust = 2, position = "stack") 

subcluster1 <- data.frame(sample = rep(c("CY15", "CY16", "CY20"), each = 5), 
                          NumCYNodes = c(sum(subcluster1_GOexp$CY15_1h_q_value < 0.005), sum(subcluster1_GOexp$CY15_3h_q_value < 0.005), sum(subcluster1_GOexp$CY15_12h_q_value < 0.005), sum(subcluster1_GOexp$CY15_24h_q_value < 0.005), sum(subcluster1_GOexp$CY15_48h_q_value < 0.005),
                                         sum(subcluster1_GOexp$CY16_1h_q_value < 0.005), sum(subcluster1_GOexp$CY16_3h_q_value < 0.005), sum(subcluster1_GOexp$CY16_12h_q_value < 0.005), sum(subcluster1_GOexp$CY16_24h_q_value < 0.005), sum(subcluster1_GOexp$CY16_48h_q_value < 0.005), 
                                         sum(subcluster1_GOexp$CY20_1h_q_value < 0.005), sum(subcluster1_GOexp$CY20_3h_q_value < 0.005), sum(subcluster1_GOexp$CY20_12h_q_value < 0.005), sum(subcluster1_GOexp$CY20_24h_q_value < 0.005), sum(subcluster1_GOexp$CY20_48h_q_value < 0.005)
                          ), 
                          Category = rep(c("1h", "3h", "12h", "24h", "48h"), times = 3)
                          )
#複数線グラフ#
library(ggplot2)
g <- ggplot(
  subcluster1_PieChart,
  aes(
    x = Category,
    y = NumCYNodes,
    group = sample, 
    color = sample
  )
)
g <- g + geom_line()
g <- g + ylab("Number of DEGs")
g <- g + xlab("time-course")
g <- g + labs(title = "subcluster1 Number of DEGs")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)

#subcluster2
n <- 1
subcluster2_Node <- c()
subcluster2_symbol <- c()
subcluster2_expression <- c()
subcluster2_Nodes <- list()
test <- c()
temp <- c()
for(n in 1:10){
  subcluster2 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 2, ]
  subcluster2_Node <- unlist(strsplit(as.character(unlist(subcluster2$Node[n])), " "))
  subcluster2_Node <- subcluster2_Node[nchar(subcluster2_Node) == 9]
  temp <- as.character(unlist(allRNASeq$genesymbol[match(subcluster2_Node, rownames(allRNASeq))]))
  subcluster2_symbol <- c(subcluster2_symbol, temp)
  test <- allRNASeq[match(subcluster2_Node, rownames(allRNASeq)), ]
  test <- test[, 1:36]
  subcluster2_expression <- rbind(subcluster2_expression, test)
  subcluster2_Nodes <- c(subcluster2_Nodes, list(subcluster2_Node))
  subcluster2_Node <- c()
  n <- n+1
}
subcluster2_GOexp <- data.frame(subcluster2_AGI = as.character(unlist(subcluster2_Nodes)),
                                subcluster2_symbol = subcluster2_symbol,
                                subcluster2_expression)

subcluster2 <- data.frame(sample = rep(c("CY15", "CY16", "CY20"), each = 5), 
                          NumCYNodes = c(sum(subcluster2_GOexp$CY15_1h_q_value < 0.005), sum(subcluster2_GOexp$CY15_3h_q_value < 0.005), sum(subcluster2_GOexp$CY15_12h_q_value < 0.005), sum(subcluster2_GOexp$CY15_24h_q_value < 0.005), sum(subcluster2_GOexp$CY15_48h_q_value < 0.005),
                                         sum(subcluster2_GOexp$CY16_1h_q_value < 0.005), sum(subcluster2_GOexp$CY16_3h_q_value < 0.005), sum(subcluster2_GOexp$CY16_12h_q_value < 0.005), sum(subcluster2_GOexp$CY16_24h_q_value < 0.005), sum(subcluster2_GOexp$CY16_48h_q_value < 0.005), 
                                         sum(subcluster2_GOexp$CY20_1h_q_value < 0.005), sum(subcluster2_GOexp$CY20_3h_q_value < 0.005), sum(subcluster2_GOexp$CY20_12h_q_value < 0.005), sum(subcluster2_GOexp$CY20_24h_q_value < 0.005), sum(subcluster2_GOexp$CY20_48h_q_value < 0.005)
                          ), 
                          Category = rep(c("1h", "3h", "12h", "24h", "48h"), times = 3)
                          )

library(ggplot2)
g <- ggplot(
  subcluster2,
  aes(
    x = Category,
    y = NumCYNodes,
    group = sample, 
    color = sample
  )
)
g <- g + geom_line()
g <- g + ylab("Number of DEGs")
g <- g + xlab("time-course")
g <- g + labs(title = "subcluster2 Number of DEGs")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)

#subcluster7
n <- 1
subcluster7_Node <- c()
subcluster7_symbol <- c()
subcluster7_expression <- c()
subcluster7_Nodes <- list()
test <- c()
temp <- c()
for(n in 1:10){
  subcluster7 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 7, ]
  subcluster7_Node <- unlist(strsplit(as.character(unlist(subcluster7$Node[n])), " "))
  subcluster7_Node <- subcluster7_Node[nchar(subcluster7_Node) == 9]
  temp <- as.character(unlist(allRNASeq$genesymbol[match(subcluster7_Node, rownames(allRNASeq))]))
  subcluster7_symbol <- c(subcluster7_symbol, temp)
  test <- allRNASeq[match(subcluster7_Node, rownames(allRNASeq)), ]
  test <- test[, 1:36]
  subcluster7_expression <- rbind(subcluster7_expression, test)
  subcluster7_Nodes <- c(subcluster7_Nodes, list(subcluster7_Node))
  subcluster7_Node <- c()
  n <- n+1
}
subcluster7_GOexp <- data.frame(subcluster7_AGI = as.character(unlist(subcluster7_Nodes)),
                                subcluster7_symbol = subcluster7_symbol,
                                subcluster7_expression)

subcluster7 <- data.frame(sample = rep(c("CY15", "CY16", "CY20"), each = 5), 
                          NumCYNodes = c(sum(subcluster7_GOexp$CY15_1h_q_value < 0.005), sum(subcluster7_GOexp$CY15_3h_q_value < 0.005), sum(subcluster7_GOexp$CY15_13h_q_value < 0.005), sum(subcluster7_GOexp$CY15_34h_q_value < 0.005), sum(subcluster7_GOexp$CY15_48h_q_value < 0.005),
                                         sum(subcluster7_GOexp$CY16_1h_q_value < 0.005), sum(subcluster7_GOexp$CY16_3h_q_value < 0.005), sum(subcluster7_GOexp$CY16_13h_q_value < 0.005), sum(subcluster7_GOexp$CY16_34h_q_value < 0.005), sum(subcluster7_GOexp$CY16_48h_q_value < 0.005), 
                                         sum(subcluster7_GOexp$CY20_1h_q_value < 0.005), sum(subcluster7_GOexp$CY20_3h_q_value < 0.005), sum(subcluster7_GOexp$CY20_12h_q_value < 0.005), sum(subcluster7_GOexp$CY20_24h_q_value < 0.005), sum(subcluster7_GOexp$CY20_48h_q_value < 0.005)
                          ), 
                          Category = rep(c("1h", "3h", "12h", "24h", "48h"), times = 3)
)

library(ggplot2)
g <- ggplot(
  subcluster7,
  aes(
    x = Category,
    y = NumCYNodes,
    group = sample, 
    color = sample
  )
)
g <- g + geom_line()
g <- g + ylab("Number of DEGs")
g <- g + xlab("time-course")
g <- g + labs(title = "subcluster7 Number of DEGs")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)


#subcluster35
n <- 1
subcluster35_Node <- c()
subcluster35_symbol <- c()
subcluster35_expression <- c()
subcluster35_Nodes <- list()
test <- c()
temp <- c()
for(n in 1:10){
  subcluster35 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 35, ]
  subcluster35_Node <- unlist(strsplit(as.character(unlist(subcluster35$Node[n])), " "))
  subcluster35_Node <- subcluster35_Node[nchar(subcluster35_Node) == 9]
  temp <- as.character(unlist(allRNASeq$genesymbol[match(subcluster35_Node, rownames(allRNASeq))]))
  subcluster35_symbol <- c(subcluster35_symbol, temp)
  test <- allRNASeq[match(subcluster35_Node, rownames(allRNASeq)), ]
  test <- test[, 1:36]
  subcluster35_expression <- rbind(subcluster35_expression, test)
  subcluster35_Nodes <- c(subcluster35_Nodes, list(subcluster35_Node))
  subcluster35_Node <- c()
  n <- n+1
}
subcluster35_GOexp <- data.frame(subcluster35_AGI = as.character(unlist(subcluster35_Nodes)),
                                 subcluster35_symbol = subcluster35_symbol,
                                 subcluster35_expression)

subcluster35 <- data.frame(sample = rep(c("CY15", "CY16", "CY20"), each = 5), 
                           NumCYNodes = c(sum(subcluster35_GOexp$CY15_1h_q_value < 0.005), sum(subcluster35_GOexp$CY15_3h_q_value < 0.005), sum(subcluster35_GOexp$CY15_13h_q_value < 0.005), sum(subcluster35_GOexp$CY15_34h_q_value < 0.005), sum(subcluster35_GOexp$CY15_48h_q_value < 0.005),
                                          sum(subcluster35_GOexp$CY16_1h_q_value < 0.005), sum(subcluster35_GOexp$CY16_3h_q_value < 0.005), sum(subcluster35_GOexp$CY16_13h_q_value < 0.005), sum(subcluster35_GOexp$CY16_34h_q_value < 0.005), sum(subcluster35_GOexp$CY16_48h_q_value < 0.005), 
                                          sum(subcluster35_GOexp$CY20_1h_q_value < 0.005), sum(subcluster35_GOexp$CY20_3h_q_value < 0.005), sum(subcluster35_GOexp$CY20_12h_q_value < 0.005), sum(subcluster35_GOexp$CY20_24h_q_value < 0.005), sum(subcluster35_GOexp$CY20_48h_q_value < 0.005)
                           ), 
                           Category = rep(c("1h", "3h", "12h", "24h", "48h"), times = 3)
)

library(ggplot2)
g <- ggplot(
  subcluster35,
  aes(
    x = Category,
    y = NumCYNodes,
    group = sample, 
    color = sample
  )
)
g <- g + geom_line()
g <- g + ylab("Number of DEGs")
g <- g + xlab("time-course")
g <- g + labs(title = "subcluster35 Number of DEGs")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
