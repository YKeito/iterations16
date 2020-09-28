genesymbol <- read.table("~/Nakano_RNAseq/network_analysis/base/allRNASeq_genesymbol.txt", sep = "\t", header = F, quote = "", row.names = 1)
genesymbol <- genesymbol[order(rownames(genesymbol)), ]

genesymbol <- genesymbol$V2[match(rownames(genesymbol), rownames(allRNASeq))]
allRNASeq <- data.frame(allRNASeq[, 1:36], 
                        symbol = genesymbol)