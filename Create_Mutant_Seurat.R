library(Seurat)
library(dplyr)
library(tibble)
library(stringr)
library(purrr)

# 1. 🔍 Seurat 객체 로딩
aml_data <- readRDS("/data/workbench/scRSEQ_AML/data/aml.test.bd.mut3_DNMT3A_mut_meta_added.Rds")

# 2. 🧬 cell_id → full Seurat cell ID 매핑
meta_df <- aml_data@meta.data %>%
  tibble::rownames_to_column("full_cell_id") %>%
  mutate(cell_id = str_extract(full_cell_id, "[^_]+$"))

# 3. 📂 RDS 파일 목록 불러오기
rds_dir <- "/data/processed_data/scRSEQ_AML/marker_macrogen/Adjusted_bam_dir/mutation_cellbarcode/merged_all_celltype/RDS_DIR"
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

# 4. 🧬 각 RDS에서 변이 어노테이션 추출
mutation_annot_df <- purrr::map_dfr(rds_files, function(f) {
  file <- basename(f)
  fields <- unlist(strsplit(file, "\\."))
  if (length(fields) < 8) return(NULL)
  
  chrom <- fields[5]
  pos <- fields[6]
  var <- fields[7]
  gene <- fields[8]
  allele_type <- fields[9]
  allele_label <- ifelse(allele_type == "ref", "REF",
                         ifelse(allele_type == "var", "ALT", NA))
  mut_id <- paste0(gene, ":", chrom, ":", pos)
  
  cell_barcodes <- tryCatch({
    cb_rds <- readRDS(f)
    rownames(cb_rds)
  }, error = function(e) NULL)
  
  if (is.null(cell_barcodes)) return(NULL)
  
  tibble(cell_id = cell_barcodes,
         Allele = allele_label,
         Mut_info = mut_id)
})

# 5. 🔗 Seurat 셀 ID에 매핑
annot_merged <- mutation_annot_df %>%
  left_join(meta_df %>% select(cell_id, full_cell_id), by = "cell_id") %>%
  filter(!is.na(full_cell_id))

# 6. ⚠️ Hetero 처리 (같은 변이에 REF + ALT 모두 있으면 → "Hetero")
hetero_cells <- annot_merged %>%
  group_by(full_cell_id, Mut_info) %>%
  summarise(n_alleles = n_distinct(Allele), .groups = "drop") %>%
  filter(n_alleles > 1) %>%
  pull(full_cell_id)

annot_clean <- annot_merged %>%
  filter(!full_cell_id %in% hetero_cells) %>%
  distinct(full_cell_id, Mut_info, Allele)

annot_hetero <- annot_merged %>%
  filter(full_cell_id %in% hetero_cells) %>%
  group_by(full_cell_id, Mut_info) %>%
  summarise(Allele = "Hetero", .groups = "drop")

# 7. ✅ 최종 어노테이션 테이블
annot_final <- bind_rows(annot_clean, annot_hetero)

# 8. 🧹 메타데이터 초기화 (덮어쓰기 전 NA → "" 처리)
aml_data@meta.data$Mut_info <- ifelse(is.na(aml_data@meta.data$Mut_info), "", aml_data@meta.data$Mut_info)
aml_data@meta.data$Allele <- ifelse(is.na(aml_data@meta.data$Allele), "", aml_data@meta.data$Allele)

# 9. 🪄 메타데이터 덮어쓰기 (벡터화 → 빠름 & 완전 대체)
target_cells <- annot_final$full_cell_id
new_mut_info <- annot_final$Mut_info
new_allele <- annot_final$Allele
names(new_mut_info) <- target_cells
names(new_allele) <- target_cells

aml_data@meta.data[target_cells, "Mut_info"] <- new_mut_info
aml_data@meta.data[target_cells, "Allele"] <- new_allele

# 10. 💾 저장 + 권한 설정
saveRDS(aml_data, file = "/data/processed_data/scRSEQ_AML/marker_macrogen/Adjusted_bam_dir/mutation_cellbarcode/aml.test.bd.mut3_DNMT3A_mut_meta_with_all_mutinfo.Rds")
system("chmod 777 /data/processed_data/scRSEQ_AML/marker_macrogen/Adjusted_bam_dir/mutation_cellbarcode/aml.test.bd.mut3_DNMT3A_mut_meta_with_all_mutinfo.Rds")
