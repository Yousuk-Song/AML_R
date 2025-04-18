# 필요한 라이브러리
library(dplyr)

# 1. 메타 데이터 로딩
mut_info <- readRDS('/data/workbench/scRSEQ_AML/data/aml.test.bd.mut3.Rds')
meta_df <- mut_info@meta.data
meta_df$cell_id <- rownames(meta_df)

# 2. pt01과 pt02 샘플별 필터링
meta_pt01 <- meta_df %>%
  filter(Sample_Tag %in% c("pt01_bm01", "pt01_bm02"))
meta_pt02 <- meta_df %>%
  filter(Sample_Tag %in% c("pt02_bm01", "pt02_bm02"))

# 3. 숫자형 cell_id 추출 (p51_123 → 123)
meta_pt01$numeric_id <- gsub("^p\\d+_", "", meta_pt01$cell_id)
meta_pt02$numeric_id <- gsub("^p\\d+_", "", meta_pt02$cell_id)

# 4. ID 매핑 테이블 생성
pt01_id_map <- meta_pt01 %>% select(cell_id, numeric_id)
pt02_id_map <- meta_pt02 %>% select(cell_id, numeric_id)


# 5. RDS 파일 로딩
rds_dir <- "/mnt/S1/sdata/processed_data/scRSEQ_AML/marker_macrogen/Adjusted_bam_dir/mutation_cellbarcode/merged_all_celltype/RDS_DIR"

# ---------- pt01 ----------
pt01_ref_df <- readRDS(file.path(rds_dir, "Combined_320-2_Bioproduct.final.sorted.valid_chrom.chr2.25234373.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds"))
pt01_var_df <- readRDS(file.path(rds_dir, "Combined_320-2_Bioproduct.final.sorted.valid_chrom.chr2.25234373.T.DNMT3A.var.bam.cell_barcodes.umi.merged.rds"))

pt01_ref_df$cell_id <- as.character(pt01_ref_df$cell_id)
pt01_var_df$cell_id <- as.character(pt01_var_df$cell_id)


# 유효 cell_id 필터링 + full ID 추가
pt01_ref_valid <- pt01_ref_df %>%
  filter(cell_id %in% pt01_id_map$numeric_id) %>%
  left_join(pt01_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

pt01_var_valid <- pt01_var_df %>%
  filter(cell_id %in% pt01_id_map$numeric_id) %>%
  left_join(pt01_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

# 추가 정보
pt01_ref_valid$mut_info <- "chr2.25234373.C>T.DNMT3A"
pt01_var_valid$mut_info <- "chr2.25234373.C>T.DNMT3A"
pt01_ref_valid$Allele <- "REF"
pt01_var_valid$Allele <- "ALT"

# ---------- pt02 ----------
pt02_ref_df <- readRDS(file.path(rds_dir, "Combined_8966_Bioproduct.final.sorted.valid_chrom.chr2.25234339.C.DNMT3A.ref.bam.cell_barcodes.umi.merged.rds"))
pt02_var_df <- readRDS(file.path(rds_dir, "Combined_8966_Bioproduct.final.sorted.valid_chrom.chr2.25234339.A.DNMT3A.var.bam.cell_barcodes.umi.merged.rds"))

pt02_ref_df$cell_id <- as.character(pt02_ref_df$cell_id)
pt02_var_df$cell_id <- as.character(pt02_var_df$cell_id)


# 유효 cell_id 필터링 + full ID 추가
pt02_ref_valid <- pt02_ref_df %>%
  filter(cell_id %in% pt02_id_map$numeric_id) %>%
  left_join(pt02_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

pt02_var_valid <- pt02_var_df %>%
  filter(cell_id %in% pt02_id_map$numeric_id) %>%
  left_join(pt02_id_map, by = c("cell_id" = "numeric_id")) %>%
  rename(numeric_cell_id = cell_id, full_cell_id = cell_id.y)

# 추가 정보
pt02_ref_valid$mut_info <- "chr2.25234339.C>A.DNMT3A"
pt02_var_valid$mut_info <- "chr2.25234339.C>A.DNMT3A"
pt02_ref_valid$Allele <- "REF"
pt02_var_valid$Allele <- "ALT"


# 1. 네 개의 데이터프레임에서 필요한 열만 추출하여 통합
mutation_annot_df <- bind_rows(
  pt01_ref_valid %>% select(full_cell_id, mut_info, Allele),
  pt01_var_valid %>% select(full_cell_id, mut_info, Allele),
  pt02_ref_valid %>% select(full_cell_id, mut_info, Allele),
  pt02_var_valid %>% select(full_cell_id, mut_info, Allele)
)

head(mutation_annot_df)
str(mutation_annot_df)

# 2. 셀 단위로 중복 Allele 확인
allele_check <- mutation_annot_df %>%
  group_by(full_cell_id) %>%
  summarise(n_alleles = n_distinct(Allele), .groups = "drop") %>%
  filter(n_alleles > 1)
allele_check

# 3. 경고 후 중단
if (nrow(allele_check) > 0) {
  warning_cells <- allele_check$full_cell_id
  stop(glue::glue("⚠️ 경고: 다음 셀들은 REF와 ALT 모두 포함하고 있어 분석을 중단합니다:\n{paste(warning_cells, collapse = ', ')}"))
}

# 2. REF & ALT 모두 포함한 셀 ID 추출
cells_with_both <- mutation_annot_df %>%
  group_by(full_cell_id) %>%
  summarise(n_alleles = n_distinct(Allele), .groups = "drop") %>%
  filter(n_alleles > 1) %>%
  pull(full_cell_id)

# 3. 해당 셀들의 read 정보 리스트로 추출
conflict_reads_list <- mutation_annot_df %>%
  filter(full_cell_id %in% cells_with_both) %>%
  group_by(full_cell_id) %>%
  group_split() %>%
  setNames(cells_with_both)

# 4. 리스트 확인 예시
print(conflict_reads_list)




library(dplyr)

# 1. conflict_reads_list의 셀들만 따로 벡터로 추출
hetero_cells <- names(conflict_reads_list)
hetero_cells

# 2. mutation_annot_df에서 Hetero 셀 제외한 순수 ALT/REF 단일 셀들만 사용
mutation_annot_clean <- mutation_annot_df %>%
  filter(!full_cell_id %in% hetero_cells) %>%
  distinct(full_cell_id, .keep_all = TRUE)
mutation_annot_clean

# 3. Hetero 셀 정보 테이블 생성
mutation_hetero_df <- mutation_annot_df %>%
  filter(full_cell_id %in% hetero_cells) %>%
  group_by(full_cell_id) %>%
  summarise(
    mut_info = mut_info[1],
    Allele = "Hetero",
    .groups = "drop"
  )
mutation_hetero_df


# 4. 둘을 합치기
mutation_for_merge <- bind_rows(mutation_annot_clean, mutation_hetero_df)

# 5. meta_df에 병합
meta_df$cell_id <- rownames(meta_df)

meta_df_annotated <- meta_df %>%
  left_join(mutation_for_merge, by = c("cell_id" = "full_cell_id"))


head(meta_df_annotated)
meta_df_annotated %>% filter(meta_df_annotated$Allele == "REF")



# cell ID 이름 확인
head(colnames(mut_info))       # 예: "p11_1234567"
head(meta_df_annotated$cell_id)

# 일치하는 경우 바로 metadata에 붙이기
mut_info <- AddMetaData(
  object = mut_info,
  metadata = meta_df_annotated %>% select(cell_id, Allele) %>% tibble::column_to_rownames("cell_id")
)

# cell ID 이름 확인
head(colnames(mut_info))       # 예: "p11_1234567"
head(meta_df_annotated$cell_id)

# 일치하는 경우 바로 metadata에 붙이기
mut_info <- AddMetaData(
  object = mut_info,
  metadata = meta_df_annotated %>% select(cell_id, Allele) %>% tibble::column_to_rownames("cell_id")
)

head(mut_info@meta.data)

saveRDS(mut_info, file = "/data/processed_data/scRSEQ_AML/sevenbridge_modified_bam_dir/aml.test.bd.mut3_DNMT3A_mut_meta_added.Rds")
system(paste0("chmod 777", " ","/data/processed_data/scRSEQ_AML/sevenbridge_modified_bam_dir/aml.test.bd.mut3_DNMT3A_mut_meta_added.Rds"))
