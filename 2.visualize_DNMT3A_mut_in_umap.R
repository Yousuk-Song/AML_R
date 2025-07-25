#!/usr/bin/env Rscript

library(Seurat)
library(SingleR)
library(SummarizedExperiment)
library(celldex)
library(ggplot2)
library(dplyr)
library(tibble)

# ▶ 인자에서 샘플명 받기
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("샘플명을 인자로 넣어주세요. 예: Rscript run_dnmt3a_umap.R 10396")
}
sample <- args[1]

# ▶ 경로 설정
output_dir <- "/data/processed_data/scRSEQ_AML/HN00247118/mutation/2020"
rds_path <- file.path(output_dir, paste0(sample, "_seurat_obj.qc_final.rds"))
seurat_obj <- readRDS(rds_path)

# ▶ SingleR로 cell type annotation
ref <- celldex::HumanPrimaryCellAtlasData()
sce <- as.SingleCellExperiment(seurat_obj)
singler <- SingleR(test = sce, ref = ref, labels = ref$label.main)
names(singler$labels) <- colnames(sce)

# SingleR annotation (full label)
full_labels <- singler$labels
names(full_labels) <- colnames(sce)

# 희귀 라벨은 "Other"로 치환 (UMAP용)
label_tbl <- table(full_labels)
keep_labels <- names(label_tbl[label_tbl >= 10])
filtered_labels <- ifelse(full_labels %in% keep_labels, full_labels, "Other")
names(filtered_labels) <- names(full_labels)

# UMAP에 사용할 필터된 라벨만 Seurat 메타데이터에 저장
seurat_obj$SingleR <- filtered_labels


# RDS 저장
saveRDS(full_labels, file.path(output_dir, paste0(sample, "_singler_result.rds")))
saveRDS(filtered_labels, file.path(output_dir, paste0(sample, "_filtered_singler_result.rds")))
system(paste("chmod 777", shQuote(file.path(output_dir, paste0(sample, "_singler_result.rds")))))
system(paste("chmod 777", shQuote(file.path(output_dir, paste0(sample, "_filtered_singler_result.rds")))))

# ▶ mutation 위치 정의
mutation_info <- list(
  "8825"  = list(chr = "chr2", pos = "25234373", ref = "C", alt = "T"),
  "10396" = list(chr = "chr2", pos = "25234373", ref = "C", alt = "T"),
  "10090" = list(chr = "chr2", pos = "25234374", ref = "G", alt = "A"),
  "9740"  = list(chr = "chr2", pos = "25240672", ref = "C", alt = "G")
)
info <- mutation_info[[sample]]

# ▶ 바코드 파일
ref_file <- file.path(output_dir, paste0(sample, "_sample_alignments.", info$chr, ".", info$pos, ".", info$ref, ".DNMT3A.ref.bam.cell_barcodes.umi.tsv"))
var_file <- file.path(output_dir, paste0(sample, "_sample_alignments.", info$chr, ".", info$pos, ".", info$alt, ".DNMT3A.var.bam.cell_barcodes.umi.tsv"))
ref_barcodes <- read.delim(ref_file, header = TRUE)
var_barcodes <- read.delim(var_file, header = TRUE)

# ▶ UMAP 및 메타데이터
umap_df <- Embeddings(seurat_obj, "umap") %>%
  as.data.frame() %>%
  rownames_to_column("Barcode")
colnames(umap_df)[2:3] <- c("UMAP_1", "UMAP_2")

meta_df <- seurat_obj@meta.data %>%
  rownames_to_column("Barcode")

plot_df <- left_join(umap_df, meta_df, by = "Barcode")

# ▶ Mutation status 지정
plot_df$MutationStatus <- "Other"
plot_df$MutationStatus[plot_df$Barcode %in% ref_barcodes$CellBarcode] <- "Ref"
plot_df$MutationStatus[plot_df$Barcode %in% var_barcodes$CellBarcode] <- "Variant"
plot_df$MutationStatus <- factor(plot_df$MutationStatus, levels = c("Other", "Ref", "Variant"))

# ▶ UMAP 범위 지정
x_range <- range(plot_df$UMAP_1, na.rm = TRUE)
y_range <- range(plot_df$UMAP_2, na.rm = TRUE)
x_margin <- 0.1 * diff(x_range)
y_margin <- 0.1 * diff(y_range)

# ▶ Plot
p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(data = subset(plot_df, MutationStatus == "Other"),
             aes(color = SingleR), size = 0.6, alpha = 0.3) +
  geom_point(data = subset(plot_df, MutationStatus == "Ref"), color = "green", size = 1.5) +
  geom_point(data = subset(plot_df, MutationStatus == "Variant"), color = "red", size = 1.5) +
  ggtitle(paste("UMAP -", sample, "| DNMT3A", info$chr, info$pos, paste0(info$ref, ">", info$alt))) +
  coord_fixed() +
  xlim(x_range[1] - x_margin, x_range[2] + x_margin) +
  ylim(y_range[1] - y_margin, y_range[2] + y_margin) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

# ▶ 저장
png_file <- file.path(output_dir, paste0(sample, "_DNMT3A_", info$chr, "_", info$pos, "_annotated_umap.png"))
ggsave(png_file, p, width = 9, height = 9, dpi = 300)
system(paste("chmod 777", shQuote(png_file)))
