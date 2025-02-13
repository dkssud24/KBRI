library(data.table)
library(dplyr)

significant_snps <- data.table()


for (i in 1:22){ 
var1 <- paste("discovery_chr",i,"_A1BG_P04217_OID30771_v1_Inflammation_II",sep="")
var2 <- paste("discovery_RS2_chunk1473_chr",i,"_A1BG_P04217_OID30771_v1_Inflammation_II.regenie",sep="")
var3 <- paste(var1,'/',var2,sep="")
var4 <- fread(var3)
var5 <- var4 [ var4$A1FREQ > 0.01,]
var5 <- data.frame(var5)
var5$P <- 10^(-var5$LOG10P)
var6 <- var5 [ var5$P < 3.40E-11,]
print(i)
print(dim(var6))
  if (nrow(var6) > 0) {
    significant_snps <- rbind(significant_snps, var6, fill = TRUE)
  }
}

significant_snps <- data.frame(significant_snps)
root_bim <- fread("C:/Users/KBRI/Desktop/정해운/0014_UKB_data/all_bim.txt.v2")

significant_snps$original <- sub("(:[^:]+){2}$", "", significant_snps$ID)


significant_snps <- significant_snps %>%
  mutate(ID_clean = sub("(:[^:]+){2}$", "", ID))




