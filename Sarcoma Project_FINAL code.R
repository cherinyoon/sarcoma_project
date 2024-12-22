#install.packages("Rcpp")  # Rcpp 패키지 설치
#install.packages("parallel")
#install.packages("survival")

library(survival) # survival 패키지 로드
library(Rcpp) # Rcpp 패키지 로드
library(parallel) # parallel 패키지 로드
library(dplyr) 




##STEP1. Data Preparation##
##STEP1. Data Preparation(drug_martix)##

# 작업 디렉토리 설정
setwd("C:/Users/유민성/Desktop/YMS/M4/M4_2/M4 AI Class/Project")

# CSV 파일 읽기
data <- read.csv("drug_to_gene_list.csv")

# geneA와 geneB 열만 선택
gene_matrix_raw <- data[c("Gene_A", "Gene_B")]

# 확인을 위해 gene_from_drug_matrix_raw.csv 저장
write.csv(gene_matrix_raw, "gene_from_drug_matrix_raw.csv", row.names = FALSE)

# NA가 하나라도 있는 행 제거
gene_matrix <- na.omit(gene_matrix_raw)

# geneA와 geneB 조합 중복 제거
gene_matrix <- gene_matrix[!duplicated(t(apply(gene_matrix, 1, sort))), ]

# gene_matrix.csv 저장
write.csv(gene_matrix, "gene_from_drug_matrix.csv", row.names = FALSE)



##STEP1. Data Preparation##
##STEP1. Data Preparation(sarcoma_martix)##

# 파일 경로 및 파일 이름 설정
file_directory <- "C:/Users/유민성/Desktop/YMS/M4/M4_2/M4 AI Class/Project/sarcoma_GEO_forPresentation"
file_names <- c("exprs_GSE21257.csv", "exprs_GSE159848.csv", "exprs_GSE159847.csv", 
                "exprs_GSE63156.csv", "exprs_GSE63155.csv")

# 전체 파일 경로
file_paths <- file.path(file_directory, file_names)

# 기준 파일 처리 (첫 번째 파일: GSE21257)
base_file <- read.csv(file_paths[1], header = TRUE, stringsAsFactors = FALSE)

# 중복된 gene symbol 제거 및 row name 설정
base_file <- base_file[!duplicated(base_file$Gene_Symbol), ]
rownames(base_file) <- base_file$Gene_Symbol
base_file$Gene_Symbol <- NULL  # Gene_Symbol 열 제거

# 기준 파일 저장
merged_data <- base_file
merged_data <- tibble::rownames_to_column(merged_data, var = "Gene_Symbol")

# 다른 파일 병합
for (i in 2:length(file_paths)) {
  cat("병합 중인 파일:", file_names[i], "\n")
  
  # 파일 불러오기
  new_file <- read.csv(file_paths[i], header = TRUE, stringsAsFactors = FALSE)
  
  # 중복된 gene symbol 제거 및 row name 설정
  new_file <- new_file[!duplicated(new_file$Gene_Symbol), ]
  rownames(new_file) <- new_file$Gene_Symbol
  new_file$Gene_Symbol <- NULL
  
  # 행 이름을 컬럼으로 변환
  new_file <- tibble::rownames_to_column(new_file, var = "Gene_Symbol")
  
  # 병합 (Gene_Symbol 기준)
  merged_data <- full_join(merged_data, new_file, by = "Gene_Symbol")
}

# 다시 rownames를 설정
rownames(merged_data) <- merged_data$Gene_Symbol
merged_data$Gene_Symbol <- NULL

# 결과 확인
head(merged_data)

# NA가 하나라도 있는 행 제거
merged_data_filtered <- merged_data[complete.cases(merged_data), ]

# 1열 1행에 "Gene_Symbol" 추가
merged_data_filtered <- tibble::rownames_to_column(merged_data_filtered, var = "Gene_Symbol")
colnames(merged_data_filtered)[1] <- "Gene_Symbol"

# 결과 저장
output_path <- "C:/Users/유민성/Desktop/YMS/M4/M4_2/M4 AI Class/Project/sarcoma_matrix_5GSE.csv"
write.csv(merged_data_filtered, output_path, row.names = FALSE)

cat("파일 병합이 완료되었습니다. 저장된 경로: ", output_path, "\n")




##STEP2. Data Loading##
##STEP2. Data Loading##

# 데이터 파일 읽기 (row.names 설정 없이)
sarcoma <- read.csv("C:/Users/유민성/Desktop/YMS/M4/M4_2/M4 AI Class/Project/sarcoma_matrix_5GSE.csv", header = TRUE)

# Gene_Symbol 열의 중복 제거(중복된 gene name 중 첫 번째만 유지)
duplicated_genes <- sarcoma$Gene_Symbol[duplicated(sarcoma$Gene_Symbol)]
sarcoma <- sarcoma[!duplicated(sarcoma$Gene_Symbol), ]

# 첫 번째 열의 내용을 sarcoma_row 변수에 복사 (후반부에서 사용될 예정)
sarcoma_row <- sarcoma[[1]]

# 첫 번째 열을 row.names로 설정하고 데이터프레임에서 제거
rownames(sarcoma) <- sarcoma[[1]]  # 첫 번째 열을 row.names로 설정
sarcoma <- sarcoma[, -1]           # 첫 번째 열 제거

# 파일 로드
drug_matrix <- read.csv("C:/Users/유민성/Desktop/YMS/M4/M4_2/M4 AI Class/Project/gene_from_drug_matrix.csv")
head(sarcoma)

# Convert sarcoma dataset row names to a vector
sarcoma_row_names <- row.names(sarcoma)
head(sarcoma_row_names)

# Find the numeric indices for each gene in the drug matrix
sl_pairs <- apply(drug_matrix, c(1, 2), function(gene) {
  which(sarcoma_row_names == gene)
})
head(sl_pairs)

# 데이터 프레임의 각 열을 숫자형으로 변환
sarcoma <- as.data.frame(lapply(sarcoma, function(x) as.numeric(as.character(x))))




##STEP3. Under-represented negative selection##
##STEP3. Under-represented negative selection ##
## STEP3. Under-represented negative selection ##

# meta 데이터 불러오기
meta <- read.csv("C:/Users/유민성/Desktop/YMS/M4/M4_2/M4 AI Class/Project/meta_sarcoma_types.csv", stringsAsFactors = FALSE)

# sarcoma 데이터 확인 (각 열은 sample 번호)
sarcoma_samples <- colnames(sarcoma)

# meta 데이터에서 sample 순서를 sarcoma와 맞추기
meta <- meta[match(sarcoma_samples, meta$sample), ]

# GEO 데이터셋별로 분리하여 가중치 부여
unique_GEO <- unique(meta$GEO)  # 고유한 GEO 번호 추출
mRNAq2_list <- list()           # 각 GEO별 mRNAq2 저장용 리스트

for (geo in unique_GEO) {
  # 각 GEO 데이터셋의 sample 인덱스 추출
  geo_samples <- meta$sample[meta$GEO == geo]
  sample_indices <- which(colnames(sarcoma) %in% geo_samples)
  
  # 해당 GEO 데이터셋만 추출
  sarcoma_geo <- sarcoma[, sample_indices]
  
  # 각 GEO별로 mRNA 발현에 따른 가중치 부여
  mRNA1_geo <- sarcoma_geo > apply(sarcoma_geo, 1, quantile, probs = 1/3, na.rm = TRUE)
  mRNA2_geo <- sarcoma_geo > apply(sarcoma_geo, 1, quantile, probs = 2/3, na.rm = TRUE)
  
  # 결과 저장
  mRNAq2_geo <- mRNA1_geo + mRNA2_geo
  mRNAq2_list[[geo]] <- mRNAq2_geo  # 리스트에 저장
}

# 모든 GEO 데이터셋별 결과를 병합
mRNAq2 <- do.call(cbind, mRNAq2_list)

setwd("C:/Users/유민성/ISLE")
sourceCpp("C:/Users/유민성/Desktop/YMS/M4/M4_2/M4 AI Class/Project/HyperGeometricTest.pair.cpp")

# 0. sl_pairs 중 공란이 포함된 행 제거
sl_pairs_clean <- sl_pairs[apply(sl_pairs, 1, function(row) all(row != "integer(0)")), ]

# 1. sl_pairs_clean을 행렬로 변환
sl_pairs_matrix <- as.matrix(sl_pairs_clean)  # 행렬로 변환

# 2. 행렬 내부 데이터를 integer 형식으로 강제 변환
storage.mode(sl_pairs_matrix) <- "integer"

# molecular omics screening (based on SCNA)
mol.mRNA = hypergeometricTestPair(scnaq=mRNAq2, pairs=sl_pairs_matrix)
i.mol=which(mol.mRNA[,1]<0.05) #adjust p.value threshold # mol.mRNA[,1]이 geneA,B =0 인 경우
sl_pair_ns =sl_pairs_matrix[i.mol,]




##STEP4. Survival Analysis##
##STEP4. Survival Analysis##

meta <- read.csv("C:/Users/유민성/Desktop/YMS/M4/M4_2/M4 AI Class/Project/meta_sarcoma_types.csv", stringsAsFactors = FALSE)

# sarcoma의 열 이름 중 meta에 존재하는 sample만 선택
valid_samples <- colnames(sarcoma) %in% meta$sample

# sarcoma 데이터 정리
sarcoma <- sarcoma[, valid_samples]

# meta 데이터 정리 (순서 맞춤)
meta <- meta[match(colnames(sarcoma), meta$sample), ]


### 1.sample/gene number
numGenes = nrow(sarcoma)
numSamples = ncol(sarcoma)


### 2.survival table(dead = 1)
survival_status <- ifelse(meta$vital_status.demographic == "Dead", 1, 0)
# survival time은 NA가 많아 최종 followup을 2024 기준으로 작성(updated date)
survival_time <- case_when(
  !is.na(meta$days_to_death.demographic) ~ meta$days_to_death.demographic,
  !is.na(meta$days_to_last_follow_up.diagnoses) ~ meta$days_to_last_follow_up.diagnoses,
  TRUE ~ 365.25 * (2024 - meta$year_of_diagnosis.diagnoses)
)
surv.dt <- data.frame(
  sample = meta$sample,
  status = survival_status,
  time = survival_time
)


### 3.age, race, sex, cancer types
# quantile normalization
qnorm.array <- function(mat){
  mat.back = mat 
  mat = mat[!is.na(mat)]
  mat = rank(mat, ties.method = "average");
  mat = qnorm(mat / (length(mat)+1));
  mat.back[!is.na(mat.back)] = mat 
  mat.back
}
age=qnorm.array(meta$age_at_index.demographic)
sex=as.character(meta$gender.demographic)
types=meta$disease_type

## cox survival analysis 



cox.pair.sl <- function(pair, use.mRNA = FALSE) {
  pair <- as.numeric(pair)
  res <- rep(NA, 8)  # 결과 저장 (4개 모델의 coef와 p-value)
  
  if (use.mRNA) {
    # mRNA 데이터만 사용
    f1 <- mRNAq2[pair[1], ]
    f2 <- mRNAq2[pair[2], ]
    
    #[임시] 기준 조건을 완화(100->10)
    if (sum(!is.na(f1)) > 10 & sum(!is.na(f2)) > 10) {
      cov <- ifelse(f1 == 0 & f2 == 0, 1, 0)
      
      # race와 types 제거한 데이터프레임
      dt1 <- data.frame(
        surv.dt, 
        cov = cov, 
        age = age, 
        sex = sex
      )
      
      # Cox 모델 실행 및 결과 저장
      cox.out <- coxph(Surv(time, status) ~ cov + strata(types), data = dt1)
      res[1:2] <- summary(cox.out)$coefficients["cov", c(1, 5)]
      
      cox.out <- coxph(Surv(time, status) ~ cov + strata(types), data = dt1)
      res[3:4] <- summary(cox.out)$coefficients["cov", c(1, 5)]
      
      cox.out <- coxph(Surv(time, status) ~ cov + age + strata(types,sex), data = dt1)
      res[5:6] <- summary(cox.out)$coefficients["cov", c(1, 5)]
      
      cox.out <- coxph(Surv(time, status) ~ cov + age + strata(types,sex), data = dt1)
      res[7:8] <- summary(cox.out)$coefficients["cov", c(1, 5)]
    }
    
  } else {
    # SCNA 데이터 사용 (만약 필요할 경우)
    stop("SCNA 데이터 분석은 현재 구현되지 않았습니다.")
  }
  
  return(res)
}


# cox.mRNA에 해당 function 적용
cox.mRNA = mclapply(
  1:nrow(sl_pair_ns),
  function(tt) cox.pair.sl(as.numeric(sl_pair_ns[tt, ]), use.mRNA = TRUE),
  mc.cores = 1
)
cox.sl.mRNA = t(do.call(cbind, cox.mRNA))
head(cox.sl.mRNA)

# FDR값 설정
FDR <- 0.05

# 4가지 cox model에 대한 결과임
i.cox1=which(cox.sl.mRNA[,1]<0 & cox.sl.mRNA[,2]<FDR)
i.cox2=which(cox.sl.mRNA[,3]<0 & cox.sl.mRNA[,4]<FDR)
i.cox3=which(cox.sl.mRNA[,5]<0 & cox.sl.mRNA[,6]<FDR)
i.cox4=which(cox.sl.mRNA[,7]<0 & cox.sl.mRNA[,8]<FDR)

# 이 중 1개라도 해당하는 pair는 모두 남김
# 빈 벡터 필터링 후 결합
all_indices <- list(i.cox1, i.cox2, i.cox3, i.cox4)  # 모든 벡터를 리스트에 저장

# Filter()를 사용하여 길이가 0이 아닌 벡터만 남김
all_indices <- Filter(function(x) length(x) > 0, all_indices)

# Reduce()와 union()을 사용해 모든 벡터를 하나로 합치기
i.cox <- Reduce(union, all_indices)

# 조건을 만족하는 sl_pair_ns 추출
if (length(i.cox) > 0) {
  sl_pair_surv <- sl_pair_ns[i.cox, , drop = FALSE]  # drop = FALSE로 항상 행렬 유지
  cSL.pairs <- cbind(sarcoma_row[sl_pair_surv[, 1]], sarcoma_row[sl_pair_surv[, 2]])
} else {
  message("조건을 만족하는 Cox 결과가 없습니다.")
  sl_pair_surv <- NULL
}
