# compare with original RDS
# 필요한 라이브러리
library(dplyr)

# 1. 메타 데이터 로딩
mut_info <- readRDS('/data/workbench/scRSEQ_AML/data/aml.test.bd.mut3.Rds')
meta_df <- mut_info@meta.data
meta_df$cell_id <- rownames(meta_df)
meta_df

# 2. pt01_bm01과 pt01_bm02만 필터링
meta_pt01 <- meta_df %>%
  filter(Sample_Tag %in% c("pt01_bm01", "pt01_bm02"))

# 3. 숫자형 cell_id 추출 (p51_123 → 123)
meta_pt01$numeric_id <- gsub("^p\\d+_", "", meta_pt01$cell_id)
head(meta_pt01)


# 4. 비교 대상 RDS 파일들 로딩
rds_dir <- "/mnt/S1/sdata/processed_data/scRSEQ_AML/marker_macrogen/Adjusted_bam_dir/mutation_cellbarcode/merged_all_celltype/RDS_DIR"

ref_file <- file.path(rds_dir, "Combined_320-2_Bioproduct.final.sorted.valid_chrom.chr2.25234373.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds")
var_file <- file.path(rds_dir, "Combined_320-2_Bioproduct.final.sorted.valid_chrom.chr2.25234373.T.DNMT3A.var.bam.cell_barcodes.umi.merged.rds")

ref_df <- readRDS(ref_file)
var_df <- readRDS(var_file)

# 5. cell_id (숫자만 있는) 기준으로 유효성 확인
ref_valid <- ref_df %>% filter(cell_id %in% meta_pt01$numeric_id)
var_valid <- var_df %>% filter(cell_id %in% meta_pt01$numeric_id)

# 각각 Allele 열 추가
ref_valid$Allele <- "REF"
var_valid$Allele <- "ALT"


# 6. 결과 출력
cat("✅ 유효한 CellBarcode 수 (ref):", nrow(ref_valid), "/", nrow(ref_df), "\n")
cat("✅ 유효한 CellBarcode 수 (var):", nrow(var_valid), "/", nrow(var_df), "\n")


# 7. singleR & final.clus 일치 여부 확인 ----------------------------------

# 먼저 비교용 meta 정보 준비 (CellBarcode, singleR, final.clus만)
# meta_check도 cell_id로 컬럼명 바꿔주기
meta_check <- meta_pt01 %>%
  select(numeric_id, singleR, final.clus) %>%
  rename(cell_id = numeric_id)

# 조인 전에 cell_id를 모두 character로 변환
ref_valid$cell_id <- as.character(ref_valid$cell_id)
var_valid$cell_id <- as.character(var_valid$cell_id)
meta_check$cell_id <- as.character(meta_check$cell_id)

head(meta_check)

# ref: join 후 일치 여부 체크
ref_check <- left_join(ref_valid, meta_check, by = "cell_id", suffix = c("_rds", "_meta"))
ref_check <- ref_check %>%
  mutate(match = (singleR_rds == singleR_meta & final.clus_rds == final.clus_meta))

var_check <- left_join(var_valid, meta_check, by = "cell_id", suffix = c("_rds", "_meta"))
var_check <- var_check %>%
  mutate(match = (singleR_rds == singleR_meta & final.clus_rds == final.clus_meta))

head(ref_check)
head(var_check)
# 결과 요약 출력
cat("\n [REF] singleR + final.clus 일치\n")
cat(" - 일치 :", sum(ref_check$match, na.rm = TRUE), "\n")
cat(" - 불일치 :", sum(ref_check$match == FALSE, na.rm = TRUE), "\n")

cat("\n [VAR] singleR + final.clus 일치\n")
cat(" - 일치 :", sum(var_check$match, na.rm = TRUE), "\n")
cat(" - 불일치 :", sum(var_check$match == FALSE, na.rm = TRUE), "\n")
