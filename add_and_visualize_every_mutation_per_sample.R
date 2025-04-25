# 필요한 라이브러리 불러오기
library(dplyr)
library(Seurat)
library(tibble)
library(glue)

# 1. Seurat 객체에서 메타데이터 추출
mut_info <- readRDS('/data/workbench/scRSEQ_AML/data/aml.test.bd.mut3.Rds')
meta_df <- mut_info@meta.data
meta_df$cell_id <- rownames(meta_df)  # cell_id 열 생성 (행 이름 → 열)

# 2. pt01, pt02, pt03, pt04, pt05, pt06, pt07 샘플 필터링
meta_pt01 <- meta_df %>% filter(Sample_Tag %in% c("pt01_bm01", "pt01_bm02"))
meta_pt02 <- meta_df %>% filter(Sample_Tag %in% c("pt02_bm01", "pt02_bm02"))
meta_pt03 <- meta_df %>% filter(Sample_Tag %in% c("pt03_bm01", "pt03_bm02"))
meta_pt04 <- meta_df %>% filter(Sample_Tag %in% c("pt04_bm01", "pt04_bm02"))
meta_pt05 <- meta_df %>% filter(Sample_Tag %in% c("pt05_bm01", "pt05_bm02"))
meta_pt06 <- meta_df %>% filter(Sample_Tag %in% c("pt06_bm01", "pt06_bm02"))
meta_pt07 <- meta_df %>% filter(Sample_Tag %in% c("pt07_bm01", "pt07_bm02"))

# 3. cell_id에서 숫자만 추출하여 numeric_id 생성 (ex: p11_123456 → 123456)
meta_pt01$numeric_id <- gsub("^p\\d+_", "", meta_pt01$cell_id)
meta_pt02$numeric_id <- gsub("^p\\d+_", "", meta_pt02$cell_id)
meta_pt03$numeric_id <- gsub("^p\\d+_", "", meta_pt03$cell_id)
meta_pt04$numeric_id <- gsub("^p\\d+_", "", meta_pt04$cell_id)
meta_pt05$numeric_id <- gsub("^p\\d+_", "", meta_pt05$cell_id)
meta_pt06$numeric_id <- gsub("^p\\d+_", "", meta_pt06$cell_id)
meta_pt07$numeric_id <- gsub("^p\\d+_", "", meta_pt07$cell_id)


# 4. ID 매핑 테이블 생성 (full ID <-> 숫자 ID)
pt01_id_map <- meta_pt01 %>% select(cell_id, numeric_id)
pt02_id_map <- meta_pt02 %>% select(cell_id, numeric_id)
pt03_id_map <- meta_pt03 %>% select(cell_id, numeric_id)
pt04_id_map <- meta_pt04 %>% select(cell_id, numeric_id)
pt05_id_map <- meta_pt05 %>% select(cell_id, numeric_id)
pt06_id_map <- meta_pt06 %>% select(cell_id, numeric_id)
pt07_id_map <- meta_pt07 %>% select(cell_id, numeric_id)

# 5. BAM 기반 변이 분석 결과 RDS 파일 로딩
rds_dir <- "/mnt/S1/sdata/processed_data/scRSEQ_AML/marker_macrogen/Adjusted_bam_dir/mutation_cellbarcode/all_cb/RDS_Dir"

# ---------- pt03 ----------
# (1) chr1.114713908.NRAS
pt03_chr1.114713908.NRAS.ref_df <- readRDS(file.path(rds_dir, "Combined_327-1_Bioproduct.final.sorted.valid_chrom.chr1.114713908.T.NRAS.ref.bam.cell_barcodes.umi.merged.rds"))
pt03_chr1.114713908.NRAS.var_df <- readRDS(file.path(rds_dir, "Combined_327-1_Bioproduct.final.sorted.valid_chrom.chr1.114713908.C.NRAS.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt03_chr1.114713908.NRAS.ref_df$cell_id <- as.character(pt03_chr1.114713908.NRAS.ref_df$cell_id)
pt03_chr1.114713908.NRAS.var_df$cell_id <- as.character(pt03_chr1.114713908.NRAS.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt03_chr1.114713908.NRAS.ref_valid <- pt03_chr1.114713908.NRAS.ref_df %>%
  filter(cell_id %in% pt03_id_map$numeric_id) %>%
  left_join(pt03_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt03_chr1.114713908.NRAS.var_valid <- pt03_chr1.114713908.NRAS.var_df %>%
  filter(cell_id %in% pt03_id_map$numeric_id) %>%
  left_join(pt03_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt03_chr1.114713908.NRAS.ref_valid$Mut_info <- "chr1.114713908.T>C.NRAS"
pt03_chr1.114713908.NRAS.var_valid$Mut_info <- "chr1.114713908.T>C.NRAS"
pt03_chr1.114713908.NRAS.ref_valid$Allele <- "REF"
pt03_chr1.114713908.NRAS.var_valid$Allele <- "ALT"

# (2) chr1.114716126.NRAS
pt03_chr1.114716126.NRAS.ref_df <- readRDS(file.path(rds_dir, "Combined_327-1_Bioproduct.final.sorted.valid_chrom.chr1.114716126.C.NRAS.ref.bam.cell_barcodes.umi.merged.rds"))
pt03_chr1.114716126.NRAS.var_df <- readRDS(file.path(rds_dir, "Combined_327-1_Bioproduct.final.sorted.valid_chrom.chr1.114716126.A.NRAS.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt03_chr1.114716126.NRAS.ref_df$cell_id <- as.character(pt03_chr1.114716126.NRAS.ref_df$cell_id)
pt03_chr1.114716126.NRAS.var_df$cell_id <- as.character(pt03_chr1.114716126.NRAS.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt03_chr1.114716126.NRAS.ref_valid <- pt03_chr1.114716126.NRAS.ref_df %>%
  filter(cell_id %in% pt03_id_map$numeric_id) %>%
  left_join(pt03_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt03_chr1.114716126.NRAS.var_valid <- pt03_chr1.114716126.NRAS.var_df %>%
  filter(cell_id %in% pt03_id_map$numeric_id) %>%
  left_join(pt03_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt03_chr1.114716126.NRAS.ref_valid$Mut_info <- "chr1.114716126.C>A.NRAS"
pt03_chr1.114716126.NRAS.var_valid$Mut_info <- "chr1.114716126.C>A.NRAS"
pt03_chr1.114716126.NRAS.ref_valid$Allele <- "REF"
pt03_chr1.114716126.NRAS.var_valid$Allele <- "ALT"

# (3) chr15.90088702.IDH2
pt03_chr15.90088702.IDH2.ref_df <- readRDS(file.path(rds_dir, "Combined_327-1_Bioproduct.final.sorted.valid_chrom.chr15.90088702.T.IDH2.ref.bam.cell_barcodes.umi.merged.rds"))
pt03_chr15.90088702.IDH2.var_df <- readRDS(file.path(rds_dir, "Combined_327-1_Bioproduct.final.sorted.valid_chrom.chr15.90088702.C.IDH2.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt03_chr15.90088702.IDH2.ref_df$cell_id <- as.character(pt03_chr15.90088702.IDH2.ref_df$cell_id)
pt03_chr15.90088702.IDH2.var_df$cell_id <- as.character(pt03_chr15.90088702.IDH2.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt03_chr15.90088702.IDH2.ref_valid <- pt03_chr15.90088702.IDH2.ref_df %>%
  filter(cell_id %in% pt03_id_map$numeric_id) %>%
  left_join(pt03_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt03_chr15.90088702.IDH2.var_valid <- pt03_chr15.90088702.IDH2.var_df %>%
  filter(cell_id %in% pt03_id_map$numeric_id) %>%
  left_join(pt03_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt03_chr15.90088702.IDH2.ref_valid$Mut_info <- "chr15.90088702.T>C.IDH2"
pt03_chr15.90088702.IDH2.var_valid$Mut_info <- "chr15.90088702.T>C.IDH2"
pt03_chr15.90088702.IDH2.ref_valid$Allele <- "REF"
pt03_chr15.90088702.IDH2.var_valid$Allele <- "ALT"

# (4) chr17.7674890.TP53
pt03_chr17.7674890.TP53.ref_df <- readRDS(file.path(rds_dir, "Combined_327-1_Bioproduct.final.sorted.valid_chrom.chr17.7674890.T.TP53.ref.bam.cell_barcodes.umi.merged.rds"))
pt03_chr17.7674890.TP53.var_df <- readRDS(file.path(rds_dir, "Combined_327-1_Bioproduct.final.sorted.valid_chrom.chr17.7674890.C.TP53.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt03_chr17.7674890.TP53.ref_df$cell_id <- as.character(pt03_chr17.7674890.TP53.ref_df$cell_id)
pt03_chr17.7674890.TP53.var_df$cell_id <- as.character(pt03_chr17.7674890.TP53.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt03_chr17.7674890.TP53.ref_valid <- pt03_chr17.7674890.TP53.ref_df %>%
  filter(cell_id %in% pt03_id_map$numeric_id) %>%
  left_join(pt03_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt03_chr17.7674890.TP53.var_valid <- pt03_chr17.7674890.TP53.var_df %>%
  filter(cell_id %in% pt03_id_map$numeric_id) %>%
  left_join(pt03_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

# 변이 정보 추가
pt03_chr17.7674890.TP53.ref_valid$Mut_info <- "chr17.7674890.T>C.TP53" 
pt03_chr17.7674890.TP53.var_valid$Mut_info <- "chr17.7674890.T>C.TP53"
pt03_chr17.7674890.TP53.ref_valid$Allele <- "REF"
pt03_chr17.7674890.TP53.var_valid$Allele <- "ALT"


# (5) chr17.7674918.TP53
pt03_chr17.7674918.TP53.ref_df <- readRDS(file.path(rds_dir, "Combined_327-1_Bioproduct.final.sorted.valid_chrom.chr17.7674918.C.TP53.ref.bam.cell_barcodes.umi.merged.rds"))
pt03_chr17.7674918.TP53.var_df <- readRDS(file.path(rds_dir, "Combined_327-1_Bioproduct.final.sorted.valid_chrom.chr17.7674918.A.TP53.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt03_chr17.7674918.TP53.ref_df$cell_id <- as.character(pt03_chr17.7674918.TP53.ref_df$cell_id)
pt03_chr17.7674918.TP53.var_df$cell_id <- as.character(pt03_chr17.7674918.TP53.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt03_chr17.7674918.TP53.ref_valid <- pt03_chr17.7674918.TP53.ref_df %>%
  filter(cell_id %in% pt03_id_map$numeric_id) %>%
  left_join(pt03_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt03_chr17.7674918.TP53.var_valid <- pt03_chr17.7674918.TP53.var_df %>%
  filter(cell_id %in% pt03_id_map$numeric_id) %>%
  left_join(pt03_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt03_chr17.7674918.TP53.ref_valid$Mut_info <- "chr17.7674918.C>A.TP53"
pt03_chr17.7674918.TP53.var_valid$Mut_info <- "chr17.7674918.C>A.TP53"
pt03_chr17.7674918.TP53.ref_valid$Allele <- "REF"
pt03_chr17.7674918.TP53.var_valid$Allele <- "ALT"

# (6) chr2.25234339.DNMT3A
pt03_chr2.25234339.DNMT3A.ref_df <- readRDS(file.path(rds_dir, "Combined_327-1_Bioproduct.final.sorted.valid_chrom.chr2.25234339.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds"))
pt03_chr2.25234339.DNMT3A.var_df <- readRDS(file.path(rds_dir, "Combined_327-1_Bioproduct.final.sorted.valid_chrom.chr2.25234339.A.DNMT3A.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt03_chr2.25234339.DNMT3A.ref_df$cell_id <- as.character(pt03_chr2.25234339.DNMT3A.ref_df$cell_id)
pt03_chr2.25234339.DNMT3A.var_df$cell_id <- as.character(pt03_chr2.25234339.DNMT3A.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt03_chr2.25234339.DNMT3A.ref_valid <- pt03_chr2.25234339.DNMT3A.ref_df %>%
  filter(cell_id %in% pt03_id_map$numeric_id) %>%
  left_join(pt03_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt03_chr2.25234339.DNMT3A.var_valid <- pt03_chr2.25234339.DNMT3A.var_df %>%
  filter(cell_id %in% pt03_id_map$numeric_id) %>%
  left_join(pt03_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt03_chr2.25234339.DNMT3A.ref_valid$Mut_info <- "chr2.25234339.C>A.DNMT3A"
pt03_chr2.25234339.DNMT3A.var_valid$Mut_info <- "chr2.25234339.C>A.DNMT3A"
pt03_chr2.25234339.DNMT3A.ref_valid$Allele <- "REF"
pt03_chr2.25234339.DNMT3A.var_valid$Allele <- "ALT"

# (7) chr2.25234373.DNMT3A
pt03_chr2.25234373.DNMT3A.ref_df <- readRDS(file.path(rds_dir, "Combined_327-1_Bioproduct.final.sorted.valid_chrom.chr2.25234373.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds"))
pt03_chr2.25234373.DNMT3A.var_df <- readRDS(file.path(rds_dir, "Combined_327-1_Bioproduct.final.sorted.valid_chrom.chr2.25234373.T.DNMT3A.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt03_chr2.25234373.DNMT3A.ref_df$cell_id <- as.character(pt03_chr2.25234373.DNMT3A.ref_df$cell_id)
pt03_chr2.25234373.DNMT3A.var_df$cell_id <- as.character(pt03_chr2.25234373.DNMT3A.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt03_chr2.25234373.DNMT3A.ref_valid <- pt03_chr2.25234373.DNMT3A.ref_df %>%
  filter(cell_id %in% pt03_id_map$numeric_id) %>%
  left_join(pt03_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt03_chr2.25234373.DNMT3A.var_valid <- pt03_chr2.25234373.DNMT3A.var_df %>%
  filter(cell_id %in% pt03_id_map$numeric_id) %>%
  left_join(pt03_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt03_chr2.25234373.DNMT3A.ref_valid$Mut_info <- "chr2.25234373.C>T.DNMT3A"
pt03_chr2.25234373.DNMT3A.var_valid$Mut_info <- "chr2.25234373.C>T.DNMT3A"
pt03_chr2.25234373.DNMT3A.ref_valid$Allele <- "REF"
pt03_chr2.25234373.DNMT3A.var_valid$Allele <- "ALT"


# ---------- pt04 ----------
# (1) chr1.114713908.NRAS
pt04_chr1.114713908.NRAS.ref_df <- readRDS(file.path(rds_dir, "Combined_320-1_Bioproduct.final.sorted.valid_chrom.chr1.114713908.T.NRAS.ref.bam.cell_barcodes.umi.merged.rds"))
pt04_chr1.114713908.NRAS.var_df <- readRDS(file.path(rds_dir, "Combined_320-1_Bioproduct.final.sorted.valid_chrom.chr1.114713908.C.NRAS.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt04_chr1.114713908.NRAS.ref_df$cell_id <- as.character(pt04_chr1.114713908.NRAS.ref_df$cell_id)
pt04_chr1.114713908.NRAS.var_df$cell_id <- as.character(pt04_chr1.114713908.NRAS.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt04_chr1.114713908.NRAS.ref_valid <- pt04_chr1.114713908.NRAS.ref_df %>%
  filter(cell_id %in% pt04_id_map$numeric_id) %>%
  left_join(pt04_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt04_chr1.114713908.NRAS.var_valid <- pt04_chr1.114713908.NRAS.var_df %>%
  filter(cell_id %in% pt04_id_map$numeric_id) %>%
  left_join(pt04_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt04_chr1.114713908.NRAS.ref_valid$Mut_info <- "chr1.114713908.T>C.NRAS"
pt04_chr1.114713908.NRAS.var_valid$Mut_info <- "chr1.114713908.T>C.NRAS"
pt04_chr1.114713908.NRAS.ref_valid$Allele <- "REF"
pt04_chr1.114713908.NRAS.var_valid$Allele <- "ALT"

# (2) chr1.114716126.NRAS
pt04_chr1.114716126.NRAS.ref_df <- readRDS(file.path(rds_dir, "Combined_320-1_Bioproduct.final.sorted.valid_chrom.chr1.114716126.C.NRAS.ref.bam.cell_barcodes.umi.merged.rds"))
pt04_chr1.114716126.NRAS.var_df <- readRDS(file.path(rds_dir, "Combined_320-1_Bioproduct.final.sorted.valid_chrom.chr1.114716126.A.NRAS.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt04_chr1.114716126.NRAS.ref_df$cell_id <- as.character(pt04_chr1.114716126.NRAS.ref_df$cell_id)
pt04_chr1.114716126.NRAS.var_df$cell_id <- as.character(pt04_chr1.114716126.NRAS.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt04_chr1.114716126.NRAS.ref_valid <- pt04_chr1.114716126.NRAS.ref_df %>%
  filter(cell_id %in% pt04_id_map$numeric_id) %>%
  left_join(pt04_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt04_chr1.114716126.NRAS.var_valid <- pt04_chr1.114716126.NRAS.var_df %>%
  filter(cell_id %in% pt04_id_map$numeric_id) %>%
  left_join(pt04_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt04_chr1.114716126.NRAS.ref_valid$Mut_info <- "chr1.114716126.C>A.NRAS"
pt04_chr1.114716126.NRAS.var_valid$Mut_info <- "chr1.114716126.C>A.NRAS"
pt04_chr1.114716126.NRAS.ref_valid$Allele <- "REF"
pt04_chr1.114716126.NRAS.var_valid$Allele <- "ALT"

# (3) chr15.90088702.IDH2
pt04_chr15.90088702.IDH2.ref_df <- readRDS(file.path(rds_dir, "Combined_320-1_Bioproduct.final.sorted.valid_chrom.chr15.90088702.T.IDH2.ref.bam.cell_barcodes.umi.merged.rds"))
pt04_chr15.90088702.IDH2.var_df <- readRDS(file.path(rds_dir, "Combined_320-1_Bioproduct.final.sorted.valid_chrom.chr15.90088702.C.IDH2.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt04_chr15.90088702.IDH2.ref_df$cell_id <- as.character(pt04_chr15.90088702.IDH2.ref_df$cell_id)
pt04_chr15.90088702.IDH2.var_df$cell_id <- as.character(pt04_chr15.90088702.IDH2.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt04_chr15.90088702.IDH2.ref_valid <- pt04_chr15.90088702.IDH2.ref_df %>%
  filter(cell_id %in% pt04_id_map$numeric_id) %>%
  left_join(pt04_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt04_chr15.90088702.IDH2.var_valid <- pt04_chr15.90088702.IDH2.var_df %>%
  filter(cell_id %in% pt04_id_map$numeric_id) %>%
  left_join(pt04_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt04_chr15.90088702.IDH2.ref_valid$Mut_info <- "chr15.90088702.T>C.IDH2"
pt04_chr15.90088702.IDH2.var_valid$Mut_info <- "chr15.90088702.T>C.IDH2"
pt04_chr15.90088702.IDH2.ref_valid$Allele <- "REF"
pt04_chr15.90088702.IDH2.var_valid$Allele <- "ALT"

# (4) chr17.7674890.TP53
pt04_chr17.7674890.TP53.ref_df <- readRDS(file.path(rds_dir, "Combined_320-1_Bioproduct.final.sorted.valid_chrom.chr17.7674890.T.TP53.ref.bam.cell_barcodes.umi.merged.rds"))
pt04_chr17.7674890.TP53.var_df <- readRDS(file.path(rds_dir, "Combined_320-1_Bioproduct.final.sorted.valid_chrom.chr17.7674890.C.TP53.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt04_chr17.7674890.TP53.ref_df$cell_id <- as.character(pt04_chr17.7674890.TP53.ref_df$cell_id)
pt04_chr17.7674890.TP53.var_df$cell_id <- as.character(pt04_chr17.7674890.TP53.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt04_chr17.7674890.TP53.ref_valid <- pt04_chr17.7674890.TP53.ref_df %>%
  filter(cell_id %in% pt04_id_map$numeric_id) %>%
  left_join(pt04_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt04_chr17.7674890.TP53.var_valid <- pt04_chr17.7674890.TP53.var_df %>%
  filter(cell_id %in% pt04_id_map$numeric_id) %>%
  left_join(pt04_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

# 변이 정보 추가
pt04_chr17.7674890.TP53.ref_valid$Mut_info <- "chr17.7674890.T>C.TP53" 
pt04_chr17.7674890.TP53.var_valid$Mut_info <- "chr17.7674890.T>C.TP53"
pt04_chr17.7674890.TP53.ref_valid$Allele <- "REF"
pt04_chr17.7674890.TP53.var_valid$Allele <- "ALT"


# (5) chr17.7674918.TP53
pt04_chr17.7674918.TP53.ref_df <- readRDS(file.path(rds_dir, "Combined_320-1_Bioproduct.final.sorted.valid_chrom.chr17.7674918.C.TP53.ref.bam.cell_barcodes.umi.merged.rds"))
pt04_chr17.7674918.TP53.var_df <- readRDS(file.path(rds_dir, "Combined_320-1_Bioproduct.final.sorted.valid_chrom.chr17.7674918.A.TP53.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt04_chr17.7674918.TP53.ref_df$cell_id <- as.character(pt04_chr17.7674918.TP53.ref_df$cell_id)
pt04_chr17.7674918.TP53.var_df$cell_id <- as.character(pt04_chr17.7674918.TP53.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt04_chr17.7674918.TP53.ref_valid <- pt04_chr17.7674918.TP53.ref_df %>%
  filter(cell_id %in% pt04_id_map$numeric_id) %>%
  left_join(pt04_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt04_chr17.7674918.TP53.var_valid <- pt04_chr17.7674918.TP53.var_df %>%
  filter(cell_id %in% pt04_id_map$numeric_id) %>%
  left_join(pt04_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt04_chr17.7674918.TP53.ref_valid$Mut_info <- "chr17.7674918.C>A.TP53"
pt04_chr17.7674918.TP53.var_valid$Mut_info <- "chr17.7674918.C>A.TP53"
pt04_chr17.7674918.TP53.ref_valid$Allele <- "REF"
pt04_chr17.7674918.TP53.var_valid$Allele <- "ALT"


# (6) chr2.25234339.DNMT3A
pt04_chr2.25234339.DNMT3A.ref_df <- readRDS(file.path(rds_dir, "Combined_320-1_Bioproduct.final.sorted.valid_chrom.chr2.25234339.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds"))
pt04_chr2.25234339.DNMT3A.var_df <- readRDS(file.path(rds_dir, "Combined_320-1_Bioproduct.final.sorted.valid_chrom.chr2.25234339.A.DNMT3A.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt04_chr2.25234339.DNMT3A.ref_df$cell_id <- as.character(pt04_chr2.25234339.DNMT3A.ref_df$cell_id)
pt04_chr2.25234339.DNMT3A.var_df$cell_id <- as.character(pt04_chr2.25234339.DNMT3A.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt04_chr2.25234339.DNMT3A.ref_valid <- pt04_chr2.25234339.DNMT3A.ref_df %>%
  filter(cell_id %in% pt04_id_map$numeric_id) %>%
  left_join(pt04_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt04_chr2.25234339.DNMT3A.var_valid <- pt04_chr2.25234339.DNMT3A.var_df %>%
  filter(cell_id %in% pt04_id_map$numeric_id) %>%
  left_join(pt04_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt04_chr2.25234339.DNMT3A.ref_valid$Mut_info <- "chr2.25234339.C>A.DNMT3A"
pt04_chr2.25234339.DNMT3A.var_valid$Mut_info <- "chr2.25234339.C>A.DNMT3A"
pt04_chr2.25234339.DNMT3A.ref_valid$Allele <- "REF"
pt04_chr2.25234339.DNMT3A.var_valid$Allele <- "ALT"


# (7) chr2.25234373.DNMT3A
pt04_chr2.25234373.DNMT3A.ref_df <- readRDS(file.path(rds_dir, "Combined_320-1_Bioproduct.final.sorted.valid_chrom.chr2.25234373.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds"))
pt04_chr2.25234373.DNMT3A.var_df <- readRDS(file.path(rds_dir, "Combined_320-1_Bioproduct.final.sorted.valid_chrom.chr2.25234373.T.DNMT3A.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt04_chr2.25234373.DNMT3A.ref_df$cell_id <- as.character(pt04_chr2.25234373.DNMT3A.ref_df$cell_id)
pt04_chr2.25234373.DNMT3A.var_df$cell_id <- as.character(pt04_chr2.25234373.DNMT3A.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt04_chr2.25234373.DNMT3A.ref_valid <- pt04_chr2.25234373.DNMT3A.ref_df %>%
  filter(cell_id %in% pt04_id_map$numeric_id) %>%
  left_join(pt04_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt04_chr2.25234373.DNMT3A.var_valid <- pt04_chr2.25234373.DNMT3A.var_df %>%
  filter(cell_id %in% pt04_id_map$numeric_id) %>%
  left_join(pt04_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt04_chr2.25234373.DNMT3A.ref_valid$Mut_info <- "chr2.25234373.C>T.DNMT3A"
pt04_chr2.25234373.DNMT3A.var_valid$Mut_info <- "chr2.25234373.C>T.DNMT3A"
pt04_chr2.25234373.DNMT3A.ref_valid$Allele <- "REF"
pt04_chr2.25234373.DNMT3A.var_valid$Allele <- "ALT"



# ---------- pt05 ----------
# (1) chr1.114713908.NRAS
pt05_chr1.114713908.NRAS.ref_df <- readRDS(file.path(rds_dir, "Combined_327-2_Bioproduct.final.sorted.valid_chrom.chr1.114713908.T.NRAS.ref.bam.cell_barcodes.umi.merged.rds"))
pt05_chr1.114713908.NRAS.var_df <- readRDS(file.path(rds_dir, "Combined_327-2_Bioproduct.final.sorted.valid_chrom.chr1.114713908.C.NRAS.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt05_chr1.114713908.NRAS.ref_df$cell_id <- as.character(pt05_chr1.114713908.NRAS.ref_df$cell_id)
pt05_chr1.114713908.NRAS.var_df$cell_id <- as.character(pt05_chr1.114713908.NRAS.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt05_chr1.114713908.NRAS.ref_valid <- pt05_chr1.114713908.NRAS.ref_df %>%
  filter(cell_id %in% pt05_id_map$numeric_id) %>%
  left_join(pt05_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt05_chr1.114713908.NRAS.var_valid <- pt05_chr1.114713908.NRAS.var_df %>%
  filter(cell_id %in% pt05_id_map$numeric_id) %>%
  left_join(pt05_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt05_chr1.114713908.NRAS.ref_valid$Mut_info <- "chr1.114713908.T>C.NRAS"
pt05_chr1.114713908.NRAS.var_valid$Mut_info <- "chr1.114713908.T>C.NRAS"
pt05_chr1.114713908.NRAS.ref_valid$Allele <- "REF"
pt05_chr1.114713908.NRAS.var_valid$Allele <- "ALT"

# (2) chr1.114716126.NRAS
pt05_chr1.114716126.NRAS.ref_df <- readRDS(file.path(rds_dir, "Combined_327-2_Bioproduct.final.sorted.valid_chrom.chr1.114716126.C.NRAS.ref.bam.cell_barcodes.umi.merged.rds"))
pt05_chr1.114716126.NRAS.var_df <- readRDS(file.path(rds_dir, "Combined_327-2_Bioproduct.final.sorted.valid_chrom.chr1.114716126.A.NRAS.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt05_chr1.114716126.NRAS.ref_df$cell_id <- as.character(pt05_chr1.114716126.NRAS.ref_df$cell_id)
pt05_chr1.114716126.NRAS.var_df$cell_id <- as.character(pt05_chr1.114716126.NRAS.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt05_chr1.114716126.NRAS.ref_valid <- pt05_chr1.114716126.NRAS.ref_df %>%
  filter(cell_id %in% pt05_id_map$numeric_id) %>%
  left_join(pt05_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt05_chr1.114716126.NRAS.var_valid <- pt05_chr1.114716126.NRAS.var_df %>%
  filter(cell_id %in% pt05_id_map$numeric_id) %>%
  left_join(pt05_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt05_chr1.114716126.NRAS.ref_valid$Mut_info <- "chr1.114716126.C>A.NRAS"
pt05_chr1.114716126.NRAS.var_valid$Mut_info <- "chr1.114716126.C>A.NRAS"
pt05_chr1.114716126.NRAS.ref_valid$Allele <- "REF"
pt05_chr1.114716126.NRAS.var_valid$Allele <- "ALT"

# (3) chr15.90088702.IDH2
pt05_chr15.90088702.IDH2.ref_df <- readRDS(file.path(rds_dir, "Combined_327-2_Bioproduct.final.sorted.valid_chrom.chr15.90088702.T.IDH2.ref.bam.cell_barcodes.umi.merged.rds"))
pt05_chr15.90088702.IDH2.var_df <- readRDS(file.path(rds_dir, "Combined_327-2_Bioproduct.final.sorted.valid_chrom.chr15.90088702.C.IDH2.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt05_chr15.90088702.IDH2.ref_df$cell_id <- as.character(pt05_chr15.90088702.IDH2.ref_df$cell_id)
pt05_chr15.90088702.IDH2.var_df$cell_id <- as.character(pt05_chr15.90088702.IDH2.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt05_chr15.90088702.IDH2.ref_valid <- pt05_chr15.90088702.IDH2.ref_df %>%
  filter(cell_id %in% pt05_id_map$numeric_id) %>%
  left_join(pt05_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt05_chr15.90088702.IDH2.var_valid <- pt05_chr15.90088702.IDH2.var_df %>%
  filter(cell_id %in% pt05_id_map$numeric_id) %>%
  left_join(pt05_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt05_chr15.90088702.IDH2.ref_valid$Mut_info <- "chr15.90088702.T>C.IDH2"
pt05_chr15.90088702.IDH2.var_valid$Mut_info <- "chr15.90088702.T>C.IDH2"
pt05_chr15.90088702.IDH2.ref_valid$Allele <- "REF"
pt05_chr15.90088702.IDH2.var_valid$Allele <- "ALT"

# (4) chr17.7674890.TP53
pt05_chr17.7674890.TP53.ref_df <- readRDS(file.path(rds_dir, "Combined_327-2_Bioproduct.final.sorted.valid_chrom.chr17.7674890.T.TP53.ref.bam.cell_barcodes.umi.merged.rds"))
pt05_chr17.7674890.TP53.var_df <- readRDS(file.path(rds_dir, "Combined_327-2_Bioproduct.final.sorted.valid_chrom.chr17.7674890.C.TP53.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt05_chr17.7674890.TP53.ref_df$cell_id <- as.character(pt05_chr17.7674890.TP53.ref_df$cell_id)
pt05_chr17.7674890.TP53.var_df$cell_id <- as.character(pt05_chr17.7674890.TP53.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt05_chr17.7674890.TP53.ref_valid <- pt05_chr17.7674890.TP53.ref_df %>%
  filter(cell_id %in% pt05_id_map$numeric_id) %>%
  left_join(pt05_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt05_chr17.7674890.TP53.var_valid <- pt05_chr17.7674890.TP53.var_df %>%
  filter(cell_id %in% pt05_id_map$numeric_id) %>%
  left_join(pt05_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

# 변이 정보 추가
pt05_chr17.7674890.TP53.ref_valid$Mut_info <- "chr17.7674890.T>C.TP53" 
pt05_chr17.7674890.TP53.var_valid$Mut_info <- "chr17.7674890.T>C.TP53"
pt05_chr17.7674890.TP53.ref_valid$Allele <- "REF"
pt05_chr17.7674890.TP53.var_valid$Allele <- "ALT"


# (5) chr17.7674918.TP53
pt05_chr17.7674918.TP53.ref_df <- readRDS(file.path(rds_dir, "Combined_327-2_Bioproduct.final.sorted.valid_chrom.chr17.7674918.C.TP53.ref.bam.cell_barcodes.umi.merged.rds"))
pt05_chr17.7674918.TP53.var_df <- readRDS(file.path(rds_dir, "Combined_327-2_Bioproduct.final.sorted.valid_chrom.chr17.7674918.A.TP53.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt05_chr17.7674918.TP53.ref_df$cell_id <- as.character(pt05_chr17.7674918.TP53.ref_df$cell_id)
pt05_chr17.7674918.TP53.var_df$cell_id <- as.character(pt05_chr17.7674918.TP53.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt05_chr17.7674918.TP53.ref_valid <- pt05_chr17.7674918.TP53.ref_df %>%
  filter(cell_id %in% pt05_id_map$numeric_id) %>%
  left_join(pt05_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt05_chr17.7674918.TP53.var_valid <- pt05_chr17.7674918.TP53.var_df %>%
  filter(cell_id %in% pt05_id_map$numeric_id) %>%
  left_join(pt05_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt05_chr17.7674918.TP53.ref_valid$Mut_info <- "chr17.7674918.C>A.TP53"
pt05_chr17.7674918.TP53.var_valid$Mut_info <- "chr17.7674918.C>A.TP53"
pt05_chr17.7674918.TP53.ref_valid$Allele <- "REF"
pt05_chr17.7674918.TP53.var_valid$Allele <- "ALT"

# (6) chr2.25234339.DNMT3A
pt05_chr2.25234339.DNMT3A.ref_df <- readRDS(file.path(rds_dir, "Combined_327-2_Bioproduct.final.sorted.valid_chrom.chr2.25234339.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds"))
pt05_chr2.25234339.DNMT3A.var_df <- readRDS(file.path(rds_dir, "Combined_327-2_Bioproduct.final.sorted.valid_chrom.chr2.25234339.A.DNMT3A.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt05_chr2.25234339.DNMT3A.ref_df$cell_id <- as.character(pt05_chr2.25234339.DNMT3A.ref_df$cell_id)
pt05_chr2.25234339.DNMT3A.var_df$cell_id <- as.character(pt05_chr2.25234339.DNMT3A.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt05_chr2.25234339.DNMT3A.ref_valid <- pt05_chr2.25234339.DNMT3A.ref_df %>%
  filter(cell_id %in% pt05_id_map$numeric_id) %>%
  left_join(pt05_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt05_chr2.25234339.DNMT3A.var_valid <- pt05_chr2.25234339.DNMT3A.var_df %>%
  filter(cell_id %in% pt05_id_map$numeric_id) %>%
  left_join(pt05_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt05_chr2.25234339.DNMT3A.ref_valid$Mut_info <- "chr2.25234339.C>A.DNMT3A"
pt05_chr2.25234339.DNMT3A.var_valid$Mut_info <- "chr2.25234339.C>A.DNMT3A"
pt05_chr2.25234339.DNMT3A.ref_valid$Allele <- "REF"
pt05_chr2.25234339.DNMT3A.var_valid$Allele <- "ALT"

# (7) chr2.25234373.DNMT3A
pt05_chr2.25234373.DNMT3A.ref_df <- readRDS(file.path(rds_dir, "Combined_327-2_Bioproduct.final.sorted.valid_chrom.chr2.25234373.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds"))
pt05_chr2.25234373.DNMT3A.var_df <- readRDS(file.path(rds_dir, "Combined_327-2_Bioproduct.final.sorted.valid_chrom.chr2.25234373.T.DNMT3A.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt05_chr2.25234373.DNMT3A.ref_df$cell_id <- as.character(pt05_chr2.25234373.DNMT3A.ref_df$cell_id)
pt05_chr2.25234373.DNMT3A.var_df$cell_id <- as.character(pt05_chr2.25234373.DNMT3A.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt05_chr2.25234373.DNMT3A.ref_valid <- pt05_chr2.25234373.DNMT3A.ref_df %>%
  filter(cell_id %in% pt05_id_map$numeric_id) %>%
  left_join(pt05_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt05_chr2.25234373.DNMT3A.var_valid <- pt05_chr2.25234373.DNMT3A.var_df %>%
  filter(cell_id %in% pt05_id_map$numeric_id) %>%
  left_join(pt05_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt05_chr2.25234373.DNMT3A.ref_valid$Mut_info <- "chr2.25234373.C>T.DNMT3A"
pt05_chr2.25234373.DNMT3A.var_valid$Mut_info <- "chr2.25234373.C>T.DNMT3A"
pt05_chr2.25234373.DNMT3A.ref_valid$Allele <- "REF"
pt05_chr2.25234373.DNMT3A.var_valid$Allele <- "ALT"


# ---------- pt06 ----------
# (1) chr1.114713908.NRAS
pt06_chr1.114713908.NRAS.ref_df <- readRDS(file.path(rds_dir, "Combined_9030_Bioproduct.final.sorted.valid_chrom.chr1.114713908.T.NRAS.ref.bam.cell_barcodes.umi.merged.rds"))
pt06_chr1.114713908.NRAS.var_df <- readRDS(file.path(rds_dir, "Combined_9030_Bioproduct.final.sorted.valid_chrom.chr1.114713908.C.NRAS.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt06_chr1.114713908.NRAS.ref_df$cell_id <- as.character(pt06_chr1.114713908.NRAS.ref_df$cell_id)
pt06_chr1.114713908.NRAS.var_df$cell_id <- as.character(pt06_chr1.114713908.NRAS.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt06_chr1.114713908.NRAS.ref_valid <- pt06_chr1.114713908.NRAS.ref_df %>%
  filter(cell_id %in% pt06_id_map$numeric_id) %>%
  left_join(pt06_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt06_chr1.114713908.NRAS.var_valid <- pt06_chr1.114713908.NRAS.var_df %>%
  filter(cell_id %in% pt06_id_map$numeric_id) %>%
  left_join(pt06_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt06_chr1.114713908.NRAS.ref_valid$Mut_info <- "chr1.114713908.T>C.NRAS"
pt06_chr1.114713908.NRAS.var_valid$Mut_info <- "chr1.114713908.T>C.NRAS"
pt06_chr1.114713908.NRAS.ref_valid$Allele <- "REF"
pt06_chr1.114713908.NRAS.var_valid$Allele <- "ALT"

# (2) chr1.114716126.NRAS
pt06_chr1.114716126.NRAS.ref_df <- readRDS(file.path(rds_dir, "Combined_9030_Bioproduct.final.sorted.valid_chrom.chr1.114716126.C.NRAS.ref.bam.cell_barcodes.umi.merged.rds"))
pt06_chr1.114716126.NRAS.var_df <- readRDS(file.path(rds_dir, "Combined_9030_Bioproduct.final.sorted.valid_chrom.chr1.114716126.A.NRAS.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt06_chr1.114716126.NRAS.ref_df$cell_id <- as.character(pt06_chr1.114716126.NRAS.ref_df$cell_id)
pt06_chr1.114716126.NRAS.var_df$cell_id <- as.character(pt06_chr1.114716126.NRAS.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt06_chr1.114716126.NRAS.ref_valid <- pt06_chr1.114716126.NRAS.ref_df %>%
  filter(cell_id %in% pt06_id_map$numeric_id) %>%
  left_join(pt06_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt06_chr1.114716126.NRAS.var_valid <- pt06_chr1.114716126.NRAS.var_df %>%
  filter(cell_id %in% pt06_id_map$numeric_id) %>%
  left_join(pt06_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt06_chr1.114716126.NRAS.ref_valid$Mut_info <- "chr1.114716126.C>A.NRAS"
pt06_chr1.114716126.NRAS.var_valid$Mut_info <- "chr1.114716126.C>A.NRAS"
pt06_chr1.114716126.NRAS.ref_valid$Allele <- "REF"
pt06_chr1.114716126.NRAS.var_valid$Allele <- "ALT"

# (3) chr15.90088702.IDH2
pt06_chr15.90088702.IDH2.ref_df <- readRDS(file.path(rds_dir, "Combined_9030_Bioproduct.final.sorted.valid_chrom.chr15.90088702.T.IDH2.ref.bam.cell_barcodes.umi.merged.rds"))
pt06_chr15.90088702.IDH2.var_df <- readRDS(file.path(rds_dir, "Combined_9030_Bioproduct.final.sorted.valid_chrom.chr15.90088702.C.IDH2.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt06_chr15.90088702.IDH2.ref_df$cell_id <- as.character(pt06_chr15.90088702.IDH2.ref_df$cell_id)
pt06_chr15.90088702.IDH2.var_df$cell_id <- as.character(pt06_chr15.90088702.IDH2.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt06_chr15.90088702.IDH2.ref_valid <- pt06_chr15.90088702.IDH2.ref_df %>%
  filter(cell_id %in% pt06_id_map$numeric_id) %>%
  left_join(pt06_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt06_chr15.90088702.IDH2.var_valid <- pt06_chr15.90088702.IDH2.var_df %>%
  filter(cell_id %in% pt06_id_map$numeric_id) %>%
  left_join(pt06_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt06_chr15.90088702.IDH2.ref_valid$Mut_info <- "chr15.90088702.T>C.IDH2"
pt06_chr15.90088702.IDH2.var_valid$Mut_info <- "chr15.90088702.T>C.IDH2"
pt06_chr15.90088702.IDH2.ref_valid$Allele <- "REF"
pt06_chr15.90088702.IDH2.var_valid$Allele <- "ALT"

# (4) chr17.7674890.TP53
pt06_chr17.7674890.TP53.ref_df <- readRDS(file.path(rds_dir, "Combined_9030_Bioproduct.final.sorted.valid_chrom.chr17.7674890.T.TP53.ref.bam.cell_barcodes.umi.merged.rds"))
pt06_chr17.7674890.TP53.var_df <- readRDS(file.path(rds_dir, "Combined_9030_Bioproduct.final.sorted.valid_chrom.chr17.7674890.C.TP53.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt06_chr17.7674890.TP53.ref_df$cell_id <- as.character(pt06_chr17.7674890.TP53.ref_df$cell_id)
pt06_chr17.7674890.TP53.var_df$cell_id <- as.character(pt06_chr17.7674890.TP53.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt06_chr17.7674890.TP53.ref_valid <- pt06_chr17.7674890.TP53.ref_df %>%
  filter(cell_id %in% pt06_id_map$numeric_id) %>%
  left_join(pt06_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt06_chr17.7674890.TP53.var_valid <- pt06_chr17.7674890.TP53.var_df %>%
  filter(cell_id %in% pt06_id_map$numeric_id) %>%
  left_join(pt06_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

# 변이 정보 추가
pt06_chr17.7674890.TP53.ref_valid$Mut_info <- "chr17.7674890.T>C.TP53" 
pt06_chr17.7674890.TP53.var_valid$Mut_info <- "chr17.7674890.T>C.TP53"
pt06_chr17.7674890.TP53.ref_valid$Allele <- "REF"
pt06_chr17.7674890.TP53.var_valid$Allele <- "ALT"


# (5) chr17.7674918.TP53
pt06_chr17.7674918.TP53.ref_df <- readRDS(file.path(rds_dir, "Combined_9030_Bioproduct.final.sorted.valid_chrom.chr17.7674918.C.TP53.ref.bam.cell_barcodes.umi.merged.rds"))
pt06_chr17.7674918.TP53.var_df <- readRDS(file.path(rds_dir, "Combined_9030_Bioproduct.final.sorted.valid_chrom.chr17.7674918.A.TP53.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt06_chr17.7674918.TP53.ref_df$cell_id <- as.character(pt06_chr17.7674918.TP53.ref_df$cell_id)
pt06_chr17.7674918.TP53.var_df$cell_id <- as.character(pt06_chr17.7674918.TP53.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt06_chr17.7674918.TP53.ref_valid <- pt06_chr17.7674918.TP53.ref_df %>%
  filter(cell_id %in% pt06_id_map$numeric_id) %>%
  left_join(pt06_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt06_chr17.7674918.TP53.var_valid <- pt06_chr17.7674918.TP53.var_df %>%
  filter(cell_id %in% pt06_id_map$numeric_id) %>%
  left_join(pt06_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt06_chr17.7674918.TP53.ref_valid$Mut_info <- "chr17.7674918.C>A.TP53"
pt06_chr17.7674918.TP53.var_valid$Mut_info <- "chr17.7674918.C>A.TP53"
pt06_chr17.7674918.TP53.ref_valid$Allele <- "REF"
pt06_chr17.7674918.TP53.var_valid$Allele <- "ALT"

# (6) chr2.25234339.DNMT3A
pt06_chr2.25234339.DNMT3A.ref_df <- readRDS(file.path(rds_dir, "Combined_9030_Bioproduct.final.sorted.valid_chrom.chr2.25234339.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds"))
pt06_chr2.25234339.DNMT3A.var_df <- readRDS(file.path(rds_dir, "Combined_9030_Bioproduct.final.sorted.valid_chrom.chr2.25234339.A.DNMT3A.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt06_chr2.25234339.DNMT3A.ref_df$cell_id <- as.character(pt06_chr2.25234339.DNMT3A.ref_df$cell_id)
pt06_chr2.25234339.DNMT3A.var_df$cell_id <- as.character(pt06_chr2.25234339.DNMT3A.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt06_chr2.25234339.DNMT3A.ref_valid <- pt06_chr2.25234339.DNMT3A.ref_df %>%
  filter(cell_id %in% pt06_id_map$numeric_id) %>%
  left_join(pt06_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt06_chr2.25234339.DNMT3A.var_valid <- pt06_chr2.25234339.DNMT3A.var_df %>%
  filter(cell_id %in% pt06_id_map$numeric_id) %>%
  left_join(pt06_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt06_chr2.25234339.DNMT3A.ref_valid$Mut_info <- "chr2.25234339.C>A.DNMT3A"
pt06_chr2.25234339.DNMT3A.var_valid$Mut_info <- "chr2.25234339.C>A.DNMT3A"
pt06_chr2.25234339.DNMT3A.ref_valid$Allele <- "REF"
pt06_chr2.25234339.DNMT3A.var_valid$Allele <- "ALT"

# (7) chr2.25234373.DNMT3A
pt06_chr2.25234373.DNMT3A.ref_df <- readRDS(file.path(rds_dir, "Combined_9030_Bioproduct.final.sorted.valid_chrom.chr2.25234373.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds"))
pt06_chr2.25234373.DNMT3A.var_df <- readRDS(file.path(rds_dir, "Combined_9030_Bioproduct.final.sorted.valid_chrom.chr2.25234373.T.DNMT3A.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt06_chr2.25234373.DNMT3A.ref_df$cell_id <- as.character(pt06_chr2.25234373.DNMT3A.ref_df$cell_id)
pt06_chr2.25234373.DNMT3A.var_df$cell_id <- as.character(pt06_chr2.25234373.DNMT3A.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt06_chr2.25234373.DNMT3A.ref_valid <- pt06_chr2.25234373.DNMT3A.ref_df %>%
  filter(cell_id %in% pt06_id_map$numeric_id) %>%
  left_join(pt06_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt06_chr2.25234373.DNMT3A.var_valid <- pt06_chr2.25234373.DNMT3A.var_df %>%
  filter(cell_id %in% pt06_id_map$numeric_id) %>%
  left_join(pt06_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt06_chr2.25234373.DNMT3A.ref_valid$Mut_info <- "chr2.25234373.C>T.DNMT3A"
pt06_chr2.25234373.DNMT3A.var_valid$Mut_info <- "chr2.25234373.C>T.DNMT3A"
pt06_chr2.25234373.DNMT3A.ref_valid$Allele <- "REF"
pt06_chr2.25234373.DNMT3A.var_valid$Allele <- "ALT"






# ---------- pt07 ----------
# (1) chr1.114713908.NRAS
pt07_chr1.114713908.NRAS.ref_df <- readRDS(file.path(rds_dir, "Combined_316_Bioproduct.final.sorted.valid_chrom.chr1.114713908.T.NRAS.ref.bam.cell_barcodes.umi.merged.rds"))
pt07_chr1.114713908.NRAS.var_df <- readRDS(file.path(rds_dir, "Combined_316_Bioproduct.final.sorted.valid_chrom.chr1.114713908.C.NRAS.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt07_chr1.114713908.NRAS.ref_df$cell_id <- as.character(pt07_chr1.114713908.NRAS.ref_df$cell_id)
pt07_chr1.114713908.NRAS.var_df$cell_id <- as.character(pt07_chr1.114713908.NRAS.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt07_chr1.114713908.NRAS.ref_valid <- pt07_chr1.114713908.NRAS.ref_df %>%
  filter(cell_id %in% pt07_id_map$numeric_id) %>%
  left_join(pt07_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt07_chr1.114713908.NRAS.var_valid <- pt07_chr1.114713908.NRAS.var_df %>%
  filter(cell_id %in% pt07_id_map$numeric_id) %>%
  left_join(pt07_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt07_chr1.114713908.NRAS.ref_valid$Mut_info <- "chr1.114713908.T>C.NRAS"
pt07_chr1.114713908.NRAS.var_valid$Mut_info <- "chr1.114713908.T>C.NRAS"
pt07_chr1.114713908.NRAS.ref_valid$Allele <- "REF"
pt07_chr1.114713908.NRAS.var_valid$Allele <- "ALT"

# (2) chr1.114716126.NRAS
pt07_chr1.114716126.NRAS.ref_df <- readRDS(file.path(rds_dir, "Combined_316_Bioproduct.final.sorted.valid_chrom.chr1.114716126.C.NRAS.ref.bam.cell_barcodes.umi.merged.rds"))
pt07_chr1.114716126.NRAS.var_df <- readRDS(file.path(rds_dir, "Combined_316_Bioproduct.final.sorted.valid_chrom.chr1.114716126.A.NRAS.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt07_chr1.114716126.NRAS.ref_df$cell_id <- as.character(pt07_chr1.114716126.NRAS.ref_df$cell_id)
pt07_chr1.114716126.NRAS.var_df$cell_id <- as.character(pt07_chr1.114716126.NRAS.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt07_chr1.114716126.NRAS.ref_valid <- pt07_chr1.114716126.NRAS.ref_df %>%
  filter(cell_id %in% pt07_id_map$numeric_id) %>%
  left_join(pt07_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt07_chr1.114716126.NRAS.var_valid <- pt07_chr1.114716126.NRAS.var_df %>%
  filter(cell_id %in% pt07_id_map$numeric_id) %>%
  left_join(pt07_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt07_chr1.114716126.NRAS.ref_valid$Mut_info <- "chr1.114716126.C>A.NRAS"
pt07_chr1.114716126.NRAS.var_valid$Mut_info <- "chr1.114716126.C>A.NRAS"
pt07_chr1.114716126.NRAS.ref_valid$Allele <- "REF"
pt07_chr1.114716126.NRAS.var_valid$Allele <- "ALT"

# (3) chr15.90088702.IDH2
pt07_chr15.90088702.IDH2.ref_df <- readRDS(file.path(rds_dir, "Combined_316_Bioproduct.final.sorted.valid_chrom.chr15.90088702.T.IDH2.ref.bam.cell_barcodes.umi.merged.rds"))
pt07_chr15.90088702.IDH2.var_df <- readRDS(file.path(rds_dir, "Combined_316_Bioproduct.final.sorted.valid_chrom.chr15.90088702.C.IDH2.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt07_chr15.90088702.IDH2.ref_df$cell_id <- as.character(pt07_chr15.90088702.IDH2.ref_df$cell_id)
pt07_chr15.90088702.IDH2.var_df$cell_id <- as.character(pt07_chr15.90088702.IDH2.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt07_chr15.90088702.IDH2.ref_valid <- pt07_chr15.90088702.IDH2.ref_df %>%
  filter(cell_id %in% pt07_id_map$numeric_id) %>%
  left_join(pt07_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt07_chr15.90088702.IDH2.var_valid <- pt07_chr15.90088702.IDH2.var_df %>%
  filter(cell_id %in% pt07_id_map$numeric_id) %>%
  left_join(pt07_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt07_chr15.90088702.IDH2.ref_valid$Mut_info <- "chr15.90088702.T>C.IDH2"
pt07_chr15.90088702.IDH2.var_valid$Mut_info <- "chr15.90088702.T>C.IDH2"
pt07_chr15.90088702.IDH2.ref_valid$Allele <- "REF"
pt07_chr15.90088702.IDH2.var_valid$Allele <- "ALT"

# (4) chr17.7674890.TP53
pt07_chr17.7674890.TP53.ref_df <- readRDS(file.path(rds_dir, "Combined_316_Bioproduct.final.sorted.valid_chrom.chr17.7674890.T.TP53.ref.bam.cell_barcodes.umi.merged.rds"))
pt07_chr17.7674890.TP53.var_df <- readRDS(file.path(rds_dir, "Combined_316_Bioproduct.final.sorted.valid_chrom.chr17.7674890.C.TP53.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt07_chr17.7674890.TP53.ref_df$cell_id <- as.character(pt07_chr17.7674890.TP53.ref_df$cell_id)
pt07_chr17.7674890.TP53.var_df$cell_id <- as.character(pt07_chr17.7674890.TP53.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt07_chr17.7674890.TP53.ref_valid <- pt07_chr17.7674890.TP53.ref_df %>%
  filter(cell_id %in% pt07_id_map$numeric_id) %>%
  left_join(pt07_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt07_chr17.7674890.TP53.var_valid <- pt07_chr17.7674890.TP53.var_df %>%
  filter(cell_id %in% pt07_id_map$numeric_id) %>%
  left_join(pt07_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

# 변이 정보 추가
pt07_chr17.7674890.TP53.ref_valid$Mut_info <- "chr17.7674890.T>C.TP53" 
pt07_chr17.7674890.TP53.var_valid$Mut_info <- "chr17.7674890.T>C.TP53"
pt07_chr17.7674890.TP53.ref_valid$Allele <- "REF"
pt07_chr17.7674890.TP53.var_valid$Allele <- "ALT"


# (5) chr17.7674918.TP53
pt07_chr17.7674918.TP53.ref_df <- readRDS(file.path(rds_dir, "Combined_316_Bioproduct.final.sorted.valid_chrom.chr17.7674918.C.TP53.ref.bam.cell_barcodes.umi.merged.rds"))
pt07_chr17.7674918.TP53.var_df <- readRDS(file.path(rds_dir, "Combined_316_Bioproduct.final.sorted.valid_chrom.chr17.7674918.A.TP53.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt07_chr17.7674918.TP53.ref_df$cell_id <- as.character(pt07_chr17.7674918.TP53.ref_df$cell_id)
pt07_chr17.7674918.TP53.var_df$cell_id <- as.character(pt07_chr17.7674918.TP53.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt07_chr17.7674918.TP53.ref_valid <- pt07_chr17.7674918.TP53.ref_df %>%
  filter(cell_id %in% pt07_id_map$numeric_id) %>%
  left_join(pt07_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt07_chr17.7674918.TP53.var_valid <- pt07_chr17.7674918.TP53.var_df %>%
  filter(cell_id %in% pt07_id_map$numeric_id) %>%
  left_join(pt07_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt07_chr17.7674918.TP53.ref_valid$Mut_info <- "chr17.7674918.C>A.TP53"
pt07_chr17.7674918.TP53.var_valid$Mut_info <- "chr17.7674918.C>A.TP53"
pt07_chr17.7674918.TP53.ref_valid$Allele <- "REF"
pt07_chr17.7674918.TP53.var_valid$Allele <- "ALT"

# (6) chr2.25234339.DNMT3A
pt07_chr2.25234339.DNMT3A.ref_df <- readRDS(file.path(rds_dir, "Combined_316_Bioproduct.final.sorted.valid_chrom.chr2.25234339.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds"))
pt07_chr2.25234339.DNMT3A.var_df <- readRDS(file.path(rds_dir, "Combined_316_Bioproduct.final.sorted.valid_chrom.chr2.25234339.A.DNMT3A.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt07_chr2.25234339.DNMT3A.ref_df$cell_id <- as.character(pt07_chr2.25234339.DNMT3A.ref_df$cell_id)
pt07_chr2.25234339.DNMT3A.var_df$cell_id <- as.character(pt07_chr2.25234339.DNMT3A.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt07_chr2.25234339.DNMT3A.ref_valid <- pt07_chr2.25234339.DNMT3A.ref_df %>%
  filter(cell_id %in% pt07_id_map$numeric_id) %>%
  left_join(pt07_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt07_chr2.25234339.DNMT3A.var_valid <- pt07_chr2.25234339.DNMT3A.var_df %>%
  filter(cell_id %in% pt07_id_map$numeric_id) %>%
  left_join(pt07_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt07_chr2.25234339.DNMT3A.ref_valid$Mut_info <- "chr2.25234339.C>A.DNMT3A"
pt07_chr2.25234339.DNMT3A.var_valid$Mut_info <- "chr2.25234339.C>A.DNMT3A"
pt07_chr2.25234339.DNMT3A.ref_valid$Allele <- "REF"
pt07_chr2.25234339.DNMT3A.var_valid$Allele <- "ALT"

# (7) chr2.25234373.DNMT3A
pt07_chr2.25234373.DNMT3A.ref_df <- readRDS(file.path(rds_dir, "Combined_316_Bioproduct.final.sorted.valid_chrom.chr2.25234373.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds"))
pt07_chr2.25234373.DNMT3A.var_df <- readRDS(file.path(rds_dir, "Combined_316_Bioproduct.final.sorted.valid_chrom.chr2.25234373.T.DNMT3A.var.bam.cell_barcodes.umi.merged.rds"))
# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt07_chr2.25234373.DNMT3A.ref_df$cell_id <- as.character(pt07_chr2.25234373.DNMT3A.ref_df$cell_id)
pt07_chr2.25234373.DNMT3A.var_df$cell_id <- as.character(pt07_chr2.25234373.DNMT3A.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt07_chr2.25234373.DNMT3A.ref_valid <- pt07_chr2.25234373.DNMT3A.ref_df %>%
  filter(cell_id %in% pt07_id_map$numeric_id) %>%
  left_join(pt07_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt07_chr2.25234373.DNMT3A.var_valid <- pt07_chr2.25234373.DNMT3A.var_df %>%
  filter(cell_id %in% pt07_id_map$numeric_id) %>%
  left_join(pt07_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
# 변이 정보 추가
pt07_chr2.25234373.DNMT3A.ref_valid$Mut_info <- "chr2.25234373.C>T.DNMT3A"
pt07_chr2.25234373.DNMT3A.var_valid$Mut_info <- "chr2.25234373.C>T.DNMT3A"
pt07_chr2.25234373.DNMT3A.ref_valid$Allele <- "REF"
pt07_chr2.25234373.DNMT3A.var_valid$Allele <- "ALT"



# 6. 변이 주석 데이터 통합 (full_cell_id, Mut_info, Allele)
mutation_annot_df <- bind_rows(
  pt03_chr1.114713908.NRAS.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt03_chr1.114716126.NRAS.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt03_chr1.114716126.NRAS.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt03_chr15.90088702.IDH2.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt03_chr17.7674890.TP53.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt03_chr17.7674890.TP53.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt03_chr17.7674918.TP53.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt03_chr2.25234339.DNMT3A.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt03_chr2.25234339.DNMT3A.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt03_chr2.25234373.DNMT3A.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt03_chr2.25234373.DNMT3A.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt04_chr1.114713908.NRAS.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt04_chr1.114716126.NRAS.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt04_chr1.114716126.NRAS.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt04_chr15.90088702.IDH2.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt04_chr17.7674890.TP53.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt04_chr17.7674918.TP53.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt04_chr2.25234339.DNMT3A.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt04_chr2.25234339.DNMT3A.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt04_chr2.25234373.DNMT3A.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt04_chr2.25234373.DNMT3A.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt05_chr1.114713908.NRAS.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt05_chr1.114713908.NRAS.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt05_chr1.114716126.NRAS.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt05_chr15.90088702.IDH2.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt05_chr17.7674890.TP53.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt05_chr17.7674918.TP53.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt05_chr2.25234339.DNMT3A.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt05_chr2.25234373.DNMT3A.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt05_chr2.25234373.DNMT3A.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt06_chr1.114713908.NRAS.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt06_chr1.114716126.NRAS.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt06_chr17.7674890.TP53.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt06_chr17.7674918.TP53.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt06_chr17.7674918.TP53.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt06_chr2.25234339.DNMT3A.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt06_chr2.25234373.DNMT3A.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt07_chr1.114713908.NRAS.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt07_chr1.114713908.NRAS.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt07_chr1.114716126.NRAS.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt07_chr1.114716126.NRAS.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt07_chr15.90088702.IDH2.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt07_chr17.7674890.TP53.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt07_chr17.7674890.TP53.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt07_chr17.7674918.TP53.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt07_chr2.25234339.DNMT3A.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt07_chr2.25234373.DNMT3A.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt07_chr2.25234373.DNMT3A.var_valid %>% select(full_cell_id, Mut_info, Allele)
)

# 7. REF와 ALT 모두 존재하는 Hetero 셀 식별
allele_check <- mutation_annot_df %>%
  group_by(full_cell_id) %>%
  summarise(n_alleles = n_distinct(Allele), .groups = "drop") %>%
  filter(n_alleles > 1)

if (nrow(allele_check) > 0) {
  warning_cells <- allele_check$full_cell_id
  message(glue("⚠️ Hetero 셀 {length(warning_cells)}개 발견됨. Hetero로 처리합니다."))
}

hetero_cells <- allele_check$full_cell_id
hetero_cells

# 8. Allele 하나만 있는 셀 필터링
mutation_annot_clean <- mutation_annot_df %>%
  filter(!full_cell_id %in% hetero_cells) %>%
  distinct(full_cell_id, .keep_all = TRUE)

# 9. Hetero 셀 따로 정리 (Mut_info는 대표값 하나만 사용)
mutation_hetero_df <- mutation_annot_df %>%
  filter(full_cell_id %in% hetero_cells) %>%
  group_by(full_cell_id) %>%
  summarise(
    Mut_info = Mut_info[1],
    Allele = "Hetero",
    .groups = "drop"
  )

# 10. REF/ALT + Hetero 병합
mutation_for_merge <- bind_rows(mutation_annot_clean, mutation_hetero_df)

# 11. meta_df에 변이 주석 붙이기
meta_df$cell_id <- rownames(meta_df)
meta_df_annotated <- meta_df %>%
  left_join(mutation_for_merge, by = c("cell_id" = "full_cell_id"))

# 12. Seurat 객체에 Mut_info + Allele 추가
mutation_metadata <- meta_df_annotated %>%
  select(cell_id, Mut_info, Allele) %>%
  column_to_rownames("cell_id")

mut_info <- AddMetaData(mut_info, metadata = mutation_metadata)
head(mut_info@meta.data)

data.frame(table(mut_info@meta.data$Mut_info))
# 13. 결과 저장
saveRDS(mut_info, file = "/mnt/S1/sdata/processed_data/scRSEQ_AML/marker_macrogen/Adjusted_bam_dir/mutation_cellbarcode/aml.test.bd.mut3_all_sample_full_mut_meta_added.Rds")
system(paste0("chmod 777 ", "/mnt/S1/sdata/processed_data/scRSEQ_AML/marker_macrogen/Adjusted_bam_dir/mutation_cellbarcode/aml.test.bd.mut3_all_sample_full_mut_meta_added.Rds"))





## visualize
mut_info <- readRDS("/mnt/S1/sdata/processed_data/scRSEQ_AML/marker_macrogen/Adjusted_bam_dir/mutation_cellbarcode/aml.test.bd.mut3_all_sample_full_mut_meta_added.Rds")

# plot_mutation_umap_plotly(mut_info, "pt03")

library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(Seurat)

plot_mutation_umap_by_patient <- function(seurat_obj, patient_tag, umap_reduction = "umap.rpca2", k = 10) {
  umap_df <- as.data.frame(Embeddings(seurat_obj, reduction = umap_reduction)) %>%
    rownames_to_column("cell_id") %>%
    left_join(seurat_obj@meta.data %>% rownames_to_column("cell_id") %>%
                select(cell_id, Sample_Tag, final.clus, Allele), by = "cell_id")
  
  # 좌표 열 이름 추출
  coord_cols <- setdiff(colnames(umap_df), c("cell_id", "Sample_Tag", "final.clus", "Allele"))
  x_col <- coord_cols[1]
  y_col <- coord_cols[2]
  
  # k-means로 클러스터링
  set.seed(42)
  umap_df$umap_cluster <- kmeans(umap_df[, coord_cols], centers = k)$cluster %>% as.character()
  
  # dominant cluster 매핑 (final.clus -> umap_cluster)
  mapping <- umap_df %>%
    count(final.clus, umap_cluster) %>%
    group_by(final.clus) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # 매핑 정보와 병합
  umap_df_mapped <- umap_df %>%
    left_join(mapping, by = "final.clus", suffix = c("", ".map"))
  
  # 유효 셀만 필터링: 해당 final.clus가 주로 속한 umap_cluster와 일치하는 셀
  umap_df_highlight <- umap_df_mapped %>%
    filter(startsWith(Sample_Tag, paste0(patient_tag, "_")),
           Allele %in% c("REF", "Hetero", "ALT"),
           umap_cluster == umap_cluster.map)
  
  # 전체 셀 배경용
  umap_df_all <- umap_df %>%
    mutate(bg_color = ifelse(final.clus %in% c("CD4_T_cell", "CD8_T_cell"), "highlight", "default"))
  
  # 클러스터 라벨 좌표 계산
  cluster_labels <- umap_df_all %>%
    filter(!is.na(final.clus)) %>%
    group_by(final.clus) %>%
    summarise(x = median(.data[[x_col]]), y = median(.data[[y_col]]), .groups = "drop")
  
  # UMAP plot 생성
  p <- ggplot() +
    geom_point(data = umap_df_all,
               aes(x = .data[[x_col]], y = .data[[y_col]], color = bg_color),
               size = 0.3, alpha = 0.4) +
    scale_color_manual(values = c("default" = "gray85", "highlight" = "lightblue"), guide = "none")
  
  # 강조 셀 추가
  p <- p +
    geom_point(data = umap_df_highlight,
               aes(x = .data[[x_col]], y = .data[[y_col]], color = Allele),
               size = 1.3, alpha = 0.9) +
    scale_color_manual(values = c("REF" = "green3", "Hetero" = "orange", "ALT" = "red",
                                  "default" = "gray85", "highlight" = "lightblue")) +
    geom_text_repel(data = cluster_labels,
                    aes(x = x, y = y, label = final.clus),
                    size = 3.5, fontface = "bold", color = "black",
                    segment.color = "gray70", segment.size = 0.3) +
    labs(title = paste("UMAP: Patient", patient_tag), color = "Allele") +
    theme_minimal()
  
  # 강조 셀 정보 반환
  highlighted_cells <- umap_df_highlight %>% select(cell_id, Sample_Tag, final.clus, Allele)
  
  return(list(plot = p, highlighted_cells = highlighted_cells))
}



# 실행: 각 환자별로 개별 plot 출력
#patients <- c("pt03", "pt04", "pt05", "pt06", "pt07")
# for (pt in patients) {
#   print(plot_mutation_umap_by_patient(mut_info, pt))
# }
print(plot_mutation_umap_by_patient(mut_info, "pt03"))
print(plot_mutation_umap_by_patient(mut_info, "pt04"))
print(plot_mutation_umap_by_patient(mut_info, "pt05"))
print(plot_mutation_umap_by_patient(mut_info, "pt06"))
print(plot_mutation_umap_by_patient(mut_info, "pt07"))


library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

count_cells_by_patient <- function(seurat_obj) {
  # 모든 셀 정보
  meta_all <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    mutate(pt_id = str_extract(Sample_Tag, "^pt[0-9]+"))  # pt03 등 추출
  
  # ALT+HET 셀 정보
  alt_het_meta <- meta_all %>%
    filter(Allele %in% c("ALT", "Hetero"))
  
  # cell type 목록
  all_cell_types <- levels(factor(seurat_obj@meta.data$final.clus))
  
  # 각 cell type별 ALT+HET count
  alt_het_count <- alt_het_meta %>%
    count(pt_id, final.clus) %>%
    complete(pt_id, final.clus = all_cell_types, fill = list(n = 0)) %>%
    rename(ALT_HET = n)
  
  # 각 cell type별 전체 셀 수 (NA 포함)
  total_count <- meta_all %>%
    count(pt_id, final.clus) %>%
    complete(pt_id, final.clus = all_cell_types, fill = list(n = 0)) %>%
    rename(TOTAL = n)
  
  # 병합 및 포맷: n / N 형태
  merged <- left_join(alt_het_count, total_count, by = c("pt_id", "final.clus")) %>%
    mutate(label = sprintf("%d / %d", ALT_HET, TOTAL)) %>%
    select(pt_id, final.clus, label)
  
  # 총합
  total_summary <- left_join(
    alt_het_count %>% group_by(pt_id) %>% summarise(ALL_ALT_Hetero_cell = sum(ALT_HET), .groups = "drop"),
    total_count %>% group_by(pt_id) %>% summarise(ALL_Total_cell = sum(TOTAL), .groups = "drop"),
    by = "pt_id"
  )
  
  # wide format 출력
  wide_table <- merged %>%
    pivot_wider(names_from = final.clus, values_from = label) %>%
    left_join(total_summary, by = "pt_id") %>%
    rename(Patient = pt_id) %>%
    select(Patient, everything())
  
  return(wide_table)
}

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

count_cells_by_patient <- function(seurat_obj) {
  # 모든 셀 정보
  meta_all <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    mutate(pt_id = str_extract(Sample_Tag, "^pt[0-9]+"))  # pt03 등 추출
  
  # ALT+HET 셀 정보
  alt_het_meta <- meta_all %>%
    filter(Allele %in% c("ALT", "Hetero"))
  
  # cell type 목록
  all_cell_types <- levels(factor(seurat_obj@meta.data$final.clus))
  
  # 각 cell type별 ALT+HET count
  alt_het_count <- alt_het_meta %>%
    count(pt_id, final.clus) %>%
    complete(pt_id, final.clus = all_cell_types, fill = list(n = 0)) %>%
    rename(ALT_HET = n)
  
  # 각 cell type별 전체 셀 수 (NA 포함)
  total_count <- meta_all %>%
    count(pt_id, final.clus) %>%
    complete(pt_id, final.clus = all_cell_types, fill = list(n = 0)) %>%
    rename(TOTAL = n)
  
  # 병합 및 포맷: n / N 형태
  merged <- left_join(alt_het_count, total_count, by = c("pt_id", "final.clus")) %>%
    mutate(label = sprintf("%d / %d", ALT_HET, TOTAL)) %>%
    select(pt_id, final.clus, label)
  
  # 총합
  total_summary <- left_join(
    alt_het_count %>% group_by(pt_id) %>% summarise(ALL_ALT_Hetero_cell = sum(ALT_HET), .groups = "drop"),
    total_count %>% group_by(pt_id) %>% summarise(ALL_Total_cell = sum(TOTAL), .groups = "drop"),
    by = "pt_id"
  )
  
  # wide format 출력
  wide_table <- merged %>%
    pivot_wider(names_from = final.clus, values_from = label) %>%
    left_join(total_summary, by = "pt_id") %>%
    rename(Patient = pt_id) %>%
    select(Patient, everything())
  
  return(wide_table)
}
patient_celltype_table <- count_cells_by_patient(mut_info)
view(patient_celltype_table)

# 실행
vaf_table_by_patient <- count_vaf_by_patient(mut_info)
vaf_table_by_patient
