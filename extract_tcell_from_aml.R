# 데이터 로딩
mut_info <- readRDS('/data/workbench/scRSEQ_AML/data/aml.test.bd.mut3.Rds')
head(mut_info@meta.data)
table(mut_info@meta.data$final.clus)

# 저장 디렉토리 생성
output_dir <- "/data/processed_data/scRSEQ_AML/MUTECT_MUTATION/normal_tumor_pair/by_sample/cell_type_annotation_all"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# cell_id 컬럼 추가
mut_info@meta.data$cell_id <- rownames(mut_info@meta.data)

# 필요한 컬럼 선택
meta_df <- mut_info@meta.data[, c("cell_id", "Sample_Tag", "orig.ident", "singleR", "final.clus")]
head(meta_df)

# Sample_Tag 기준으로 나누어 저장 (NA로만 이루어진 행 제거 포함)
unique_samples <- unique(meta_df$Sample_Tag)
for (sample in unique_samples) {
  subset_df <- meta_df[meta_df$Sample_Tag == sample, ]
  # NA로만 이루어진 행 제거
  subset_df <- subset_df[rowSums(is.na(subset_df)) < ncol(subset_df), ]

  if (nrow(subset_df) > 0) {
    out_path <- file.path(output_dir, paste0("cell_type_id_list_", sample, ".tsv"))
    write.table(subset_df, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}


# 대상 셀만 필터링
normal_meta <- mut_info@meta.data[mut_info@meta.data$final.clus %in% c("CD4_T_cell", "CD8_T_cell"), ]

head(normal_meta)

# 유일한 (singleR, final.clus) 조합 구하기
#unique_combos <- unique(normal_meta[, c("singleR", "final.clus")])
#unique_combos


# 저장할 기본 경로
base_path <- "/mnt/S1/sdata/processed_data/scRSEQ_AML/MUTECT_MUTATION/normal_tumor_pair/by_sample"
system(paste0("chmod 777", " ", base_path))

# 디렉토리 없으면 생성
dir.create(base_path, recursive = TRUE, showWarnings = FALSE)

# orig.ident별로 분리하여 저장
split_list <- split(normal_meta, normal_meta$Sample_Tag)

# 4. 그룹별로 파일 저장
for (ident in names(split_list)) {
  sub_meta <- split_list[[ident]]

  df <- data.frame(
    cell_id = rownames(sub_meta),
    Sample_Tag = sub_meta$Sample_Tag,
    orig.ident = sub_meta$orig.ident,
    singleR = sub_meta$singleR,
    final.clus = sub_meta$final.clus,
    stringsAsFactors = FALSE
  )

  file_name <- paste0("normal_cell_id_list_", ident, ".tsv")
  file_path <- file.path(base_path, file_name)

  write.table(df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  system(paste("chmod 777", shQuote(file_path)))

  message("✅ 저장 완료: ", file_path)
}
