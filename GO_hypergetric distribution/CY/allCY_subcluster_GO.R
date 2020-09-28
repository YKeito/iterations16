####CY15####
i <- 1
total <- length(CY15_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY15 <- list()
EnrichedGO_pvalue_CY15 <- list()

CY15_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY15_enriched_Node[[i]])
  T_GO <- unlist(GO_List[match(T_AGI, names(GO_List))])
  k <- length(T_GO)
  T_GO_Uni <- unique(T_GO)
  
  j <- 1
  total_j <- length(T_GO_Uni)
  for (j in j:total_j) {
    AGI_name <- c()
    q <- sum(T_GO == T_GO_Uni[j])
    m <- sum(T_GO_Full == T_GO_Uni[j])
    n <- length(T_GO_Full) - m
    
    temp <- phyper(q-1, m, n, k, lower.tail=FALSE)
    temp_all_pV <- c(temp_all_pV, temp)
    temp_all_GO <- c(temp_all_GO, T_GO_Uni[j])
    j <- j+1
  }
  temp_CY15 <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY15)))
  names(AGI_name) <- NumCY15_subcluster[i]
  CY15_Node <- c(CY15_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY15 <- c(EnrichedGO_pvalue_CY15, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY15)[i] <- NumCY15_subcluster[i]
  
  EnrichedGO_CY15 <- c(EnrichedGO_CY15, list(temp_all_GO))
  names(EnrichedGO_CY15)[i] <- NumCY15_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY15_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY15_Node.RData")
####CY16####
i <- 1
total <- length(CY16_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY16 <- list()
EnrichedGO_pvalue_CY16 <- list()

CY16_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY16_enriched_Node[[i]])
  T_GO <- unlist(GO_List[match(T_AGI, names(GO_List))])
  k <- length(T_GO)
  T_GO_Uni <- unique(T_GO)
  
  j <- 1
  total_j <- length(T_GO_Uni)
  for (j in j:total_j) {
    AGI_name <- c()
    q <- sum(T_GO == T_GO_Uni[j])
    m <- sum(T_GO_Full == T_GO_Uni[j])
    n <- length(T_GO_Full) - m
    
    temp <- phyper(q-1, m, n, k, lower.tail=FALSE)
    temp_all_pV <- c(temp_all_pV, temp)
    temp_all_GO <- c(temp_all_GO, T_GO_Uni[j])
    j <- j+1
  }
  temp_CY16 <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY16)))
  names(AGI_name) <- NumCY16_subcluster[i]
  CY16_Node <- c(CY16_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY16 <- c(EnrichedGO_pvalue_CY16, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY16)[i] <- NumCY16_subcluster[i]
  
  EnrichedGO_CY16 <- c(EnrichedGO_CY16, list(temp_all_GO))
  names(EnrichedGO_CY16)[i] <- NumCY16_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY16_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY16_Node.RData")
####CY20####
i <- 1
total <- length(CY20_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY20 <- list()
EnrichedGO_pvalue_CY20 <- list()

CY20_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY20_enriched_Node[[i]])
  T_GO <- unlist(GO_List[match(T_AGI, names(GO_List))])
  k <- length(T_GO)
  T_GO_Uni <- unique(T_GO)
  
  j <- 1
  total_j <- length(T_GO_Uni)
  for (j in j:total_j) {
    AGI_name <- c()
    q <- sum(T_GO == T_GO_Uni[j])
    m <- sum(T_GO_Full == T_GO_Uni[j])
    n <- length(T_GO_Full) - m
    
    temp <- phyper(q-1, m, n, k, lower.tail=FALSE)
    temp_all_pV <- c(temp_all_pV, temp)
    temp_all_GO <- c(temp_all_GO, T_GO_Uni[j])
    j <- j+1
  }
  temp_CY20 <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY20)))
  names(AGI_name) <- NumCY20_subcluster[i]
  CY20_Node <- c(CY20_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY20 <- c(EnrichedGO_pvalue_CY20, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY20)[i] <- NumCY20_subcluster[i]
  
  EnrichedGO_CY20 <- c(EnrichedGO_CY20, list(temp_all_GO))
  names(EnrichedGO_CY20)[i] <- NumCY20_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY20_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY20_Node.RData")