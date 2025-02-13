
# 📌 1. 패키지 로드
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(MR.SPI)

# 📌 2. 작업 디렉토리 설정
setwd("C:/Users/KBRI/Downloads/A1BG_P04217_OID30771_v1_Inflammation_II/A1BG_P04217_OID30771_v1_Inflammation_II")

# 📌 3. 유의한 SNP 필터링 (크로모좀별 데이터 처리)
significant_snps <- data.table()
for (i in 1:22) { 
  var1 <- paste0("discovery_chr", i, "_A1BG_P04217_OID30771_v1_Inflammation_II")
  var2 <- paste0("discovery_RS2_chunk1473_chr", i, "_A1BG_P04217_OID30771_v1_Inflammation_II.regenie")
  var3 <- file.path(var1, var2)
  
  var4 <- fread(var3)
  var5 <- var4[var4$A1FREQ > 0.01, ]
  var5$P <- 10^(-var5$LOG10P)
  var6 <- var5[var5$P < 3.40E-11, ]
  
  if (nrow(var6) > 0) {
    significant_snps <- rbind(significant_snps, var6, fill = TRUE)
  }
}

# 📌 4. SNP ID 가공
significant_snps <- data.frame(significant_snps)
significant_snps$original <- sub("(:[^:]+){2}$", "", significant_snps$ID)

# 📌 5. BIM 파일 로드 및 SNP 정보 매칭
root_bim <- fread("C:/Users/KBRI/Desktop/정해운/0014_UKB_data/all_bim.txt.v2")
m <- left_join(significant_snps, root_bim, by = "original")

# 📌 6. 필요한 컬럼만 선택 후 정리
m2 <- m[, c(1, 23, 26, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15)]
colnames(m2) <- c("CHR", "SNP", "BP", "A1", "A2", "A2FREQ", "INFO", "N", "BETA", "SE", "CHISQ", "LOG10P", "P")
m3 <- na.omit(m2)
fwrite(m3,'C:/Users/KBRI/Downloads/tmp/A1BG_sig.txt',quote=FALSE,row.names=FALSE,sep="\t")

#Plink
C:\Users\KBRI\Desktop\정해운\0015_plink\plink_win64_20241022\plink --bfile C:\Users\KBRI\Desktop\정해운\0014_UKB_data\ref --clump C:\Users\KBRI\Downloads\tmp\A1BG_sig.txt --clump-kb 1000 --clump-r2 0.01 --clump-p1 1 --clump-p2 1 --out :\Users\KBRI\Downloads\tmp\result_A1BG_sig.txt

# 📌 7. PLINK Clumping 결과 파일 로드 및 매칭
m4 <- fread("C:/Users/KBRI/Downloads/tmp/result_A1BG_sig.txt.clumped")
m5 <- data.frame(m4$SNP)
colnames(m5) <- c("SNP")
m6 <- left_join(m5, m3, by = "SNP")
rm(root_bim)

# 📌 8. Outcome GWAS 데이터 로드
outcome <- fread("C:/Users/KBRI/Desktop/GCST90027158_buildGRCh38.tsv")

# 📌 9. pQTL 데이터 컬럼 정리
pqtl_data <- m6[, c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "A2FREQ")]
colnames(pqtl_data) <- c(
  "SNP", "CHR", "BP", "effect_allele.exposure", "other_allele.exposure", 
  "beta.exposure", "se.exposure", "pval.exposure", "eaf.exposure"
)
pqtl_data$exposure <- "Protein"
pqtl_data$id.exposure <- pqtl_data$SNP

# 📌 10. Outcome 데이터 컬럼 정리
colnames(outcome) <- c(
  "SNP", "pval.outcome", "CHR", "BP", "effect_allele.outcome", "other_allele.outcome", 
  "eaf.outcome", "odds_ratio", "ci_lower", "ci_upper", "beta.outcome", 
  "se.outcome", "n_cases", "n_controls", "het_isq", "het_pvalue", "variant_alternate_id"
)
ad_gwas_data <- outcome
rm(outcome)  # 메모리 절약

ad_gwas_data$outcome <- "Alzheimer's Disease"
ad_gwas_data$id.outcome <- ad_gwas_data$SNP

# 📌 11. 데이터 프레임 변환 (TwoSampleMR 사용을 위해)
pqtl_data <- data.frame(pqtl_data)
ad_gwas_data <- data.frame(ad_gwas_data)

# 📌 12. Harmonisation 수행
dat <- harmonise_data(pqtl_data, ad_gwas_data, action = 1)

# 📌 13. MR-SPI 입력 데이터 생성
mr_spi_input <- dat %>%
  dplyr::select(
    SNP, beta.exposure, se.exposure, beta.outcome, se.outcome, 
    eaf.exposure
  ) %>%
  na.omit()  # NA 값 제거 (없어야 오류 안 남)

# 📌 14. 샘플 사이즈 설정
n1 <- 54306  # exposure sample size
n2 <- 455258 # outcome sample size

# 📌 15. MR-SPI 실행
mr_spi_result <- MR.SPI(
  gamma = mr_spi_input$beta.exposure,
  Gamma = mr_spi_input$beta.outcome,
  se_gamma = mr_spi_input$se.exposure,
  se_Gamma = mr_spi_input$se.outcome,
  n1 = n1,
  n2 = n2,
  freq = mr_spi_input$eaf.exposure,
  max_clique = TRUE,
  verbose = TRUE
)

# 📌 16. 결과 출력
print(mr_spi_result)

