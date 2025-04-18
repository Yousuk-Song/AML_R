# 필요한 라이브러리 불러오기
library(dplyr)
library(Seurat)
library(tibble)
library(glue)

# 1. Seurat 객체에서 메타데이터 추출
mut_info <- readRDS('/data/workbench/scRSEQ_AML/data/aml.test.bd.mut3.Rds')
meta_df <- mut_info@meta.data
meta_df$cell_id <- rownames(meta_df)  # cell_id 열 생성 (행 이름 → 열)

# 2. pt01, pt02 샘플만 필터링
meta_pt01 <- meta_df %>% filter(Sample_Tag %in% c("pt01_bm01", "pt01_bm02"))
meta_pt02 <- meta_df %>% filter(Sample_Tag %in% c("pt02_bm01", "pt02_bm02"))

# 3. cell_id에서 숫자만 추출하여 numeric_id 생성 (ex: p11_123456 → 123456)
meta_pt01$numeric_id <- gsub("^p\\d+_", "", meta_pt01$cell_id)
meta_pt02$numeric_id <- gsub("^p\\d+_", "", meta_pt02$cell_id)

# 4. ID 매핑 테이블 생성 (full ID <-> 숫자 ID)
pt01_id_map <- meta_pt01 %>% select(cell_id, numeric_id)
pt02_id_map <- meta_pt02 %>% select(cell_id, numeric_id)

# 5. BAM 기반 변이 분석 결과 RDS 파일 로딩
rds_dir <- "/mnt/S1/sdata/processed_data/scRSEQ_AML/marker_macrogen/Adjusted_bam_dir/mutation_cellbarcode/merged_all_celltype/RDS_DIR"

# ---------- pt01 ----------
pt01_ref_df <- readRDS(file.path(rds_dir, "Combined_320-2_Bioproduct.final.sorted.valid_chrom.chr2.25234373.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds"))
pt01_var_df <- readRDS(file.path(rds_dir, "Combined_320-2_Bioproduct.final.sorted.valid_chrom.chr2.25234373.T.DNMT3A.var.bam.cell_barcodes.umi.merged.rds"))

# cell_id를 문자형으로 변환 (ID 매칭 오류 방지)
pt01_ref_df$cell_id <- as.character(pt01_ref_df$cell_id)
pt01_var_df$cell_id <- as.character(pt01_var_df$cell_id)

# 유효 cell_id만 필터링 후 full ID 매핑
pt01_ref_valid <- pt01_ref_df %>%
  filter(cell_id %in% pt01_id_map$numeric_id) %>%
  left_join(pt01_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt01_var_valid <- pt01_var_df %>%
  filter(cell_id %in% pt01_id_map$numeric_id) %>%
  left_join(pt01_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

# 변이 정보 추가
pt01_ref_valid$Mut_info <- "chr2.25234373.C>T.DNMT3A"
pt01_var_valid$Mut_info <- "chr2.25234373.C>T.DNMT3A"
pt01_ref_valid$Allele <- "REF"
pt01_var_valid$Allele <- "ALT"

# ---------- pt02 ----------
pt02_ref_df <- readRDS(file.path(rds_dir, "Combined_8966_Bioproduct.final.sorted.valid_chrom.chr2.25234339.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds"))
pt02_var_df <- readRDS(file.path(rds_dir, "Combined_8966_Bioproduct.final.sorted.valid_chrom.chr2.25234339.A.DNMT3A.var.bam.cell_barcodes.umi.merged.rds"))

pt02_ref_df$cell_id <- as.character(pt02_ref_df$cell_id)
pt02_var_df$cell_id <- as.character(pt02_var_df$cell_id)

pt02_ref_valid <- pt02_ref_df %>%
  filter(cell_id %in% pt02_id_map$numeric_id) %>%
  left_join(pt02_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)
pt02_var_valid <- pt02_var_df %>%
  filter(cell_id %in% pt02_id_map$numeric_id) %>%
  left_join(pt02_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

pt02_ref_valid$Mut_info <- "chr2.25234339.C>A.DNMT3A"
pt02_var_valid$Mut_info <- "chr2.25234339.C>A.DNMT3A"
pt02_ref_valid$Allele <- "REF"
pt02_var_valid$Allele <- "ALT"

# 6. 변이 주석 데이터 통합 (full_cell_id, Mut_info, Allele)
mutation_annot_df <- bind_rows(
  pt01_ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt01_var_valid %>% select(full_cell_id, Mut_info, Allele),
  pt02_ref_valid %>% select(full_cell_id, Mut_info, Allele),
  pt02_var_valid %>% select(full_cell_id, Mut_info, Allele)
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

# 13. 결과 저장
saveRDS(mut_info, file = "/data/processed_data/scRSEQ_AML/sevenbridge_modified_bam_dir/aml.test.bd.mut3_DNMT3A_mut_meta_added.Rds")
system(paste0("chmod 777 ", "/data/processed_data/scRSEQ_AML/sevenbridge_modified_bam_dir/aml.test.bd.mut3_DNMT3A_mut_meta_added.Rds"))

