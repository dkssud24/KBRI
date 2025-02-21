# 📌 1. 패키지 로드 및 설치
if (!require(archive)) install.packages("archive")
if (!require(data.table)) install.packages("data.table")
if (!require(dplyr)) install.packages("dplyr")

library(archive)
library(data.table)
library(dplyr)

# 📌 2. .tar 파일 내부 파일 목록 확인
tar_file <- "APOE_P02649_OID30727_v1_Inflammation_II.tar"
archive_contents <- archive::archive(tar_file)

# 📌 3. Chromosome 관련 파일 경로 추출
chromosome_files <- archive_contents$path[grepl("chr[0-9]+|chr1[0-9]|chr2[0-2]", archive_contents$path)]

# 확인
print(chromosome_files)

# 📌 4. Chromosome 데이터 읽기 (에러 없는 버전)
chromosome_data <- lapply(chromosome_files, function(file_path) {
  
  # .tar 내부에서 파일 읽기
  archive_read <- archive::archive_read(tar_file, file = file_path)
  
  # gzcon으로 압축 해제 없이 메모리에서 읽기
  con <- gzcon(archive_read)
  
  # 텍스트 데이터로 읽어오기
  raw_text <- readLines(con, warn = FALSE)
  
  # fread가 문자열을 직접 읽을 수 있도록 변환
  text_input <- paste(raw_text, collapse = "\n")
  
  # fread로 데이터프레임 생성
  data <- fread(text_input, na.strings = c("", "NA"))
  
  # 연결 닫기
  close(archive_read)
  
  return(data)
})

# 📌 5. 리스트 확인
# 리스트 구조 확인
str(chromosome_data)

# 첫 번째 chromosome 데이터 확인
head(chromosome_data[[1]])

# 테이블 형태 확인
table(chromosome_data[[1]])

# 📌 6. 리스트에 이름 부여 (chr1, chr2, ..., chr22)
names(chromosome_data) <- paste0("chr", 1:length(chromosome_data))

# 예시: chr1 데이터 확인
head(chromosome_data$chr1)

# 📌 7. 필요한 경우 모든 데이터를 하나의 데이터프레임으로 병합
combined_data <- rbindlist(chromosome_data, fill = TRUE)

# 확인
head(combined_data)
dim(combined_data)  # 전체 행/열 개수
