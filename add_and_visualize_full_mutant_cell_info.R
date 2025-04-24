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
rds_dir <- "/mnt/S1/sdata/processed_data/scRSEQ_AML/marker_macrogen/Adjusted_bam_dir/mutation_cellbarcode/merged_all_celltype/RDS_DIR"

# ---------- pt01 ----------
# (1) chr2.25234373.DNMT3A
pt01_chr2.25234373.DNMT3A.ref_df <- readRDS(file.path(rds_dir, "Combined_320-2_Bioproduct.final.sorted.valid_chrom.chr2.25234373.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds"))
pt01_chr2.25234373.DNMT3A.var_df <- readRDS(file.path(rds_dir, "Combined_320-2_Bioproduct.final.sorted.valid_chrom.chr2.25234373.T.DNMT3A.var.bam.cell_barcodes.umi.merged.rds"))

# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt01_chr2.25234373.DNMT3A.ref_df$cell_id <- as.character(pt01_chr2.25234373.DNMT3A.ref_df$cell_id)
pt01_chr2.25234373.DNMT3A.var_df$cell_id <- as.character(pt01_chr2.25234373.DNMT3A.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt01_chr2.25234373.DNMT3A.ref_valid <- pt01_ref_df %>%
  filter(cell_id %in% pt01_id_map$numeric_id) %>%
  left_join(pt01_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt01_chr2.25234373.DNMT3A.var_valid <- pt01_var_df %>%
  filter(cell_id %in% pt01_id_map$numeric_id) %>%
  left_join(pt01_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

# 변이 정보 추가
pt01_chr2.25234373.DNMT3A.ref_valid$Mut_info <- "chr2.25234373.C>T.DNMT3A"
pt01_chr2.25234373.DNMT3A.var_valid$Mut_info <- "chr2.25234373.C>T.DNMT3A"
pt01_chr2.25234373.DNMT3A.ref_valid$Allele <- "REF"
pt01_chr2.25234373.DNMT3A.var_valid$Allele <- "ALT"


# (2) chr17.7674918.TP53
pt01_chr17.7674918.TP53.ref_df <- readRDS(file.path(rds_dir, "Combined_320-2_Bioproduct.final.sorted.valid_chrom.chr17.7674918.A.TP53.ref.bam.cell_barcodes.umi.merged.rds"))
pt01_chr17.7674918.TP53.var_df <- readRDS(file.path(rds_dir, "Combined_320-2_Bioproduct.final.sorted.valid_chrom.chr17.7674918.C.TP53.var.bam.cell_barcodes.umi.merged.rds"))

# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt01_chr17.7674918.TP53.ref_df$cell_id <- as.character(pt01_chr17.7674918.TP53.ref_df$cell_id)
pt01_chr17.7674918.TP53.var_df$cell_id <- as.character(pt01_chr17.7674918.TP53.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt01_chr17.7674918.TP53.ref_valid <- pt01_chr17.7674918.TP53.ref_df %>%
  filter(cell_id %in% pt01_id_map$numeric_id) %>%
  left_join(pt01_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt01_chr17.7674918.TP53.var_valid <- pt01_var_df %>%
  filter(cell_id %in% pt01_id_map$numeric_id) %>%
  left_join(pt01_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

# 변이 정보 추가
pt01_chr17.7674918.TP53.ref_valid$Mut_info <- "chr17.7674918.C>A.TP53"
pt01_chr17.7674918.TP53.var_valid$Mut_info <- "chr17.7674918.C>A.TP53"
pt01_chr17.7674918.TP53.ref_valid$Allele <- "REF"
pt01_chr17.7674918.TP53.var_valid$Allele <- "ALT"

# (3) chr15.90088702.IDH2
pt01_chr15.90088702.IDH2.ref_df <- readRDS(file.path(rds_dir, "Combined_320-2_Bioproduct.final.sorted.valid_chrom.chr15.90088702.C.IDH2.ref.bam.cell_barcodes.umi.merged.rds"))
pt01_chr15.90088702.IDH2.var_df <- readRDS(file.path(rds_dir, "Combined_320-2_Bioproduct.final.sorted.valid_chrom.chr15.90088702.T.IDH2.var.bam.cell_barcodes.umi.merged.rds"))

# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt01_chr15.90088702.IDH2.ref_df$cell_id <- as.character(pt01_chr15.90088702.IDH2.ref_df$cell_id)
pt01_chr15.90088702.IDH2.var_df$cell_id <- as.character(pt01_chr15.90088702.IDH2.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt01_chr15.90088702.IDH2.ref_valid <- pt01_chr15.90088702.IDH2.ref_df %>%
  filter(cell_id %in% pt01_id_map$numeric_id) %>%
  left_join(pt01_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt01_chr15.90088702.IDH2.var_valid <- pt01_chr15.90088702.IDH2.var_df %>%
  filter(cell_id %in% pt01_id_map$numeric_id) %>%
  left_join(pt01_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

# 변이 정보 추가
pt01_chr15.90088702.IDH2.ref_valid$Mut_info <- "chr15.90088702.T>C.IDH2"
pt01_chr15.90088702.IDH2.var_valid$Mut_info <- "chr15.90088702.T>C.IDH2"
pt01_chr15.90088702.IDH2.ref_valid$Allele <- "REF"
pt01_chr15.90088702.IDH2.var_valid$Allele <- "ALT"


# ---------- pt02 ----------
# (1) chr15.90088702.IDH2
pt02_chr15.90088702.IDH2.ref_df <- readRDS(file.path(rds_dir, "Combined_8966_Bioproduct.final.sorted.valid_chrom.chr15.90088702.C.IDH2.ref.bam.cell_barcodes.umi.merged.rds"))
pt02_chr15.90088702.IDH2.var_df <- readRDS(file.path(rds_dir, "Combined_8966_Bioproduct.final.sorted.valid_chrom.chr15.90088702.T.IDH2.var.bam.cell_barcodes.umi.merged.rds"))

# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt02_chr15.90088702.IDH2.ref_df$cell_id <- as.character(pt02_chr15.90088702.IDH2.ref_df$cell_id)
pt02_chr15.90088702.IDH2.var_df$cell_id <- as.character(pt02_chr15.90088702.IDH2.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt02_chr15.90088702.IDH2.ref_valid <- pt02_chr15.90088702.IDH2.ref_df %>%
  filter(cell_id %in% pt02_id_map$numeric_id) %>%
  left_join(pt02_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt02_chr15.90088702.IDH2.var_valid <- pt02_chr15.90088702.IDH2.var_df %>%
  filter(cell_id %in% pt02_id_map$numeric_id) %>%
  left_join(pt02_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

# 변이 정보 추가
pt02_chr15.90088702.IDH2.ref_valid$Mut_info <- "chr15.90088702.T>C.IDH2"
pt02_chr15.90088702.IDH2.var_valid$Mut_info <- "chr15.90088702.T>C.IDH2"
pt02_chr15.90088702.IDH2.ref_valid$Allele <- "REF"
pt02_chr15.90088702.IDH2.var_valid$Allele <- "ALT"

#[2] chr2.25234339.DNMT3A
pt02_chr2.25234339.DNMT3A.ref_df <- readRDS(file.path(rds_dir, "Combined_8966_Bioproduct.final.sorted.valid_chrom.chr2.25234339.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds"))
pt02_chr2.25234339.DNMT3A.var_df <- readRDS(file.path(rds_dir, "Combined_8966_Bioproduct.final.sorted.valid_chrom.chr2.25234339.A.DNMT3A.var.bam.cell_barcodes.umi.merged.rds"))

# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt02_chr2.25234339.DNMT3A.ref_df$cell_id <- as.character(pt02_chr2.25234339.DNMT3A.ref_df$cell_id)
pt02_chr2.25234339.DNMT3A.var_df$cell_id <- as.character(pt02_chr2.25234339.DNMT3A.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt02_chr2.25234339.DNMT3A.ref_valid <- pt02_chr2.25234339.DNMT3A.ref_df %>%
  filter(cell_id %in% pt02_id_map$numeric_id) %>%
  left_join(pt02_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt02_chr2.25234339.DNMT3A.var_valid <- pt02_chr2.25234339.DNMT3A.var_df %>%
  filter(cell_id %in% pt02_id_map$numeric_id) %>%
  left_join(pt02_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

# 변이 정보 추가
pt02_chr2.25234339.DNMT3A.ref_valid$Mut_info <- "chr2.25234339.C>A.DNMT3A"
pt02_chr2.25234339.DNMT3A.var_valid$Mut_info <- "chr2.25234339.C>A.DNMT3A"
pt02_chr2.25234339.DNMT3A.ref_valid$Allele <- "REF"
pt02_chr2.25234339.DNMT3A.var_valid$Allele <- "ALT"

# (3) chr1.114716126.NRAS
pt02_chr1.114716126.NRAS.ref_df <- readRDS(file.path(rds_dir, "Combined_8966_Bioproduct.final.sorted.valid_chrom.chr1.114716126.C.NRAS.ref.bam.cell_barcodes.umi.merged.rds"))
pt02_chr1.114716126.NRAS.var_df <- readRDS(file.path(rds_dir, "Combined_8966_Bioproduct.final.sorted.valid_chrom.chr1.114716126.A.NRAS.var.bam.cell_barcodes.umi.merged.rds"))

# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt02_chr1.114716126.NRAS.ref_df$cell_id <- as.character(pt02_chr1.114716126.NRAS.ref_df$cell_id)
pt02_chr1.114716126.NRAS.var_df$cell_id <- as.character(pt02_chr1.114716126.NRAS.var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt02_chr1.114716126.NRAS.ref_valid <- pt02_chr1.114716126.NRAS.ref_df %>%
  filter(cell_id %in% pt02_id_map$numeric_id) %>%
  left_join(pt02_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt02_chr1.114716126.NRAS.var_valid <- pt02_chr1.114716126.NRAS.var_df %>%
  filter(cell_id %in% pt02_id_map$numeric_id) %>%
  left_join(pt02_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

# 변이 정보 추가
pt02_chr1.114716126.NRAS.ref_valid$Mut_info <- "chr1.114716126.C>A.NRAS" 
pt02_chr1.114716126.NRAS.var_valid$Mut_info <- "chr1.114716126.C>A.NRAS"
pt02_chr1.114716126.NRAS.ref_valid$Allele <- "REF"
pt02_chr1.114716126.NRAS.var_valid$Allele <- "ALT"

# ---------- pt03 ----------
# (1) chr17.7674890.TP53
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


# ---------- pt07 ----------
#[1] 
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

# 6. 변이 주석 데이터 통합 (full_cell_id, Mut_info, Allele)
mutation_annot_df <- bind_rows(
  pt01_chr2.25234373.DNMT3A.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt01_chr17.7674918.TP53.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt01_chr15.90088702.IDH2.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt02_chr15.90088702.IDH2.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt02_chr2.25234339.DNMT3A.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt02_chr1.114716126.NRAS.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt03_chr17.7674890.TP53.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt07_chr1.114713908.NRAS.ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt01_chr2.25234373.DNMT3A.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt01_chr17.7674918.TP53.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt01_chr15.90088702.IDH2.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt02_chr15.90088702.IDH2.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt02_chr2.25234339.DNMT3A.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt02_chr1.114716126.NRAS.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt03_chr17.7674890.TP53.var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt07_chr1.114713908.NRAS.var_valid %>% select(full_cell_id, Mut_info, Allele)
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

table(mut_info@meta.data$Mut_info)
# 13. 결과 저장
saveRDS(mut_info, file = "/mnt/S1/sdata/processed_data/scRSEQ_AML/marker_macrogen/Adjusted_bam_dir/mutation_cellbarcode/Adjusted_bam_diraml.test.bd.mut3_full_mut_meta_added.Rds")
system(paste0("chmod 777 ", "/mnt/S1/sdata/processed_data/scRSEQ_AML/marker_macrogen/Adjusted_bam_dir/mutation_cellbarcode/Adjusted_bam_diraml.test.bd.mut3_full_mut_meta_added.Rds"))


mut_to_pt <- list(
  "chr2.25234373.C>T.DNMT3A" = c("pt01"),
  "chr17.7674918.C>A.TP53" = c("pt01"),
  "chr15.90088702.T>C.IDH2" = c("pt01", "pt02"),
  "chr2.25234339.C>A.DNMT3A" = c("pt02"),
  "chr1.114716126.C>A.NRAS" = c("pt02"),
  "chr17.7674890.T>C.TP53" = c("pt03"),
  "chr1.114713908.T>C.NRAS" = c("pt07")
)

library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library(glue)

plot_mutation_umap_highlight <- function(seurat_obj, mut_name, mut_to_pt, umap_reduction = "umap.rpca2") {
  # UMAP 좌표 + 메타데이터 병합
  umap_df <- as.data.frame(Embeddings(seurat_obj, reduction = umap_reduction)) %>%
    tibble::rownames_to_column("cell_id") %>%
    left_join(
      seurat_obj@meta.data %>%
        tibble::rownames_to_column("cell_id") %>%
        dplyr::select(cell_id, Mut_info, Allele, Sample_Tag, final.clus),
      by = "cell_id"
    )
  
  # 좌표 열 자동 탐지
  coord_cols <- setdiff(colnames(umap_df), c("cell_id", "Mut_info", "Allele", "Sample_Tag", "final.clus"))
  x_col <- coord_cols[1]
  y_col <- coord_cols[2]
  
  # 해당 변이에 속한 pt
  pts <- mut_to_pt[[mut_name]]
  pt_pattern <- paste0("^(", paste(pts, collapse = "|"), ")_bm0[1-2]$")
  
  # 배경: UMAP 전체
  umap_df_gray <- umap_df
  
  # 강조: 해당 mutation + 해당 pt의 셀만
  umap_df_highlight <- umap_df %>%
    filter(Mut_info == mut_name & grepl(pt_pattern, Sample_Tag)) %>%
    mutate(Allele = factor(Allele, levels = c("REF", "ALT", "Hetero")))
  
  # Cell type 라벨 좌표 계산
  cluster_labels <- umap_df %>%
    filter(!is.na(final.clus)) %>%
    group_by(final.clus) %>%
    summarise(
      x = median(.data[[x_col]]),
      y = median(.data[[y_col]]),
      .groups = "drop"
    )
  
  # Plot
  ggplot() +
    geom_point(data = umap_df_gray, aes(x = .data[[x_col]], y = .data[[y_col]]),
               color = "gray85", size = 0.3, alpha = 0.5) +
    geom_point(data = umap_df_highlight,
               aes(x = .data[[x_col]], y = .data[[y_col]], color = Allele),
               size = 1.2, alpha = 0.9) +
    scale_color_manual(values = c("REF" = "green3", "ALT" = "red", "Hetero" = "orange")) +
    geom_text_repel(data = cluster_labels,
                    aes(x = x, y = y, label = final.clus),
                    size = 3.5, fontface = "bold", color = "black",
                    segment.color = "gray70", segment.size = 0.3) +
    labs(
      title = glue::glue("UMAP: {mut_name} (Samples: {paste(pts, collapse = ', ')})"),
      color = "Allele"
    ) +
    theme_minimal()
}
# 시각화
plot_mutation_umap_highlight(mut_info, "chr15.90088702.T>C.IDH2", mut_to_pt)
plot_mutation_umap_highlight(mut_info, "chr17.7674918.C>A.TP53", mut_to_pt)
plot_mutation_umap_highlight(mut_info, "chr2.25234373.C>T.DNMT3A", mut_to_pt)
plot_mutation_umap_highlight(mut_info, "chr17.7674890.T>C.TP53", mut_to_pt)
plot_mutation_umap_highlight(mut_info, "chr2.25234339.C>A.DNMT3A", mut_to_pt)
plot_mutation_umap_highlight(mut_info, "chr1.114716126.C>A.NRAS", mut_to_pt)
plot_mutation_umap_highlight(mut_info, "chr1.114713908.T>C.NRAS", mut_to_pt)



count_alt_hetero_tcells_with_total <- function(seurat_obj, mut_name, mut_to_pt) {
  pts <- mut_to_pt[[mut_name]]
  pt_pattern <- paste0("^(", paste(pts, collapse = "|"), ")_bm0[1-2]$")
  
  meta_all <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(Mut_info == mut_name,
           Allele %in% c("ALT", "Hetero"),
           grepl(pt_pattern, Sample_Tag))
  
  meta_tcells <- meta_all %>%
    filter(final.clus %in% c("CD4_T_cell", "CD8_T_cell"))
  
  # 클러스터 수 세기 (고정 factor levels 설정)
  meta_tcells$final.clus <- factor(meta_tcells$final.clus, levels = c("CD4_T_cell", "CD8_T_cell"))
  
  count_tbl <- meta_tcells %>%
    count(final.clus, .drop = FALSE) %>%
    tidyr::pivot_wider(names_from = final.clus, values_from = n, values_fill = 0) %>%
    mutate(Mut_info = mut_name,
           ALL_ALT_Hetero_cell = nrow(meta_all)) %>%
    select(Mut_info, everything())
  
  return(count_tbl)
}

alt_hetero_tcell_counts_final <- lapply(names(mut_to_pt), function(mut) {
  count_alt_hetero_tcells_with_total(mut_info, mut, mut_to_pt)
}) %>%
  bind_rows()

alt_hetero_tcell_counts_final

data.frame(table(mut_info@meta.data$final.clus))


count_all_cells_with_vaf <- function(seurat_obj, mut_name, mut_to_pt) {
  pts <- mut_to_pt[[mut_name]]
  pt_pattern <- paste0("^(", paste(pts, collapse = "|"), ")_bm0[1-2]$")
  
  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(Mut_info == mut_name,
           grepl(pt_pattern, Sample_Tag),
           !is.na(Allele))
  
  # Allele 범주 고정
  meta$Allele <- factor(meta$Allele, levels = c("REF", "ALT", "Hetero"))
  
  # 전체 cell type 목록 고정
  all_cell_types <- levels(factor(seurat_obj@meta.data$final.clus))
  
  # 각 cell type별 ALT+Hetero count
  alt_het_count <- meta %>%
    filter(Allele %in% c("ALT", "Hetero")) %>%
    count(final.clus) %>%
    complete(final.clus = all_cell_types, fill = list(n = 0)) %>%
    rename(ALT_HET = n)
  
  # 각 cell type별 ALT+HET+REF 전체 count
  total_count <- meta %>%
    count(final.clus) %>%
    complete(final.clus = all_cell_types, fill = list(n = 0)) %>%
    rename(TOTAL = n)
  
  # 병합 후 VAF 계산
  merged <- left_join(alt_het_count, total_count, by = "final.clus") %>%
    mutate(
      VAF = ifelse(TOTAL > 0, ALT_HET / TOTAL, NA_real_),
      label = sprintf("%d (%.2f)", ALT_HET, VAF)
    ) %>%
    select(final.clus, label)
  
  # 전체 ALT+Hetero 셀 수
  all_alt_het_count <- sum(alt_het_count$ALT_HET)
  
  # Wide format + 추가 열
  wide_table <- merged %>%
    pivot_wider(names_from = final.clus, values_from = label) %>%
    mutate(
      Mut_info = mut_name,
      ALL_ALT_Hetero_cell = all_alt_het_count
    ) %>%
    select(Mut_info, everything())  # ensure column order
  
  return(wide_table)
}

all_celltype_vaf_table <- lapply(names(mut_to_pt), function(mut) {
  count_all_cells_with_vaf(mut_info, mut, mut_to_pt)
}) %>% bind_rows()

all_celltype_vaf_table





