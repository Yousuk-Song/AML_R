
#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript qc_seurat_doublet_final_professor_feedback.R <sample>")
sample <- args[1]

# 기본 설정
cutoff_percent_mt <- 10
cutoff_nFeature_lower <- 300
cutoff_nFeature_upper <- 6000
cutoff_complexity <- 0.8
expected_rate <- switch(sample, "9740" = 0.10, "10090" = 0.08, 0.06)

# 경로 설정
base_dir <- file.path("/data/processed_data/scRSEQ_AML/HN00247118/cellranger_work", paste0(sample, "_multi"), "outs", "per_sample_outs", paste0(sample, "_multi"))
h5_path <- file.path(base_dir, "count", "sample_filtered_feature_bc_matrix.h5")
output_dir <- "/data/processed_data/scRSEQ_AML/HN00247118/mutation/2020"
output_rds <- file.path(output_dir, paste0(sample, "_seurat_obj.qc_final.rds"))
output_doublet_umap <- file.path(output_dir, paste0(sample, "_doublet_umap.png"))
output_log <- file.path(output_dir, paste0(sample, "_qc_full_log.txt"))
dir.create(dirname(output_log), showWarnings = FALSE)

# ▶ QC 플롯 저장 함수
save_plot <- function(p, path) {
  ggsave(path, p, width = 18, height = 5, dpi = 300)
  system(paste("chmod 777", shQuote(path)))
}

# ▶ QC 플롯 생성 함수
#make_qc_plots <- function(meta, prefix) {
#  save_plot(
#    ggplot(meta, aes(x = nFeature_RNA)) +
#      geom_density(fill = "orange", alpha = 0.4) + scale_x_log10() +
#      geom_vline(xintercept = c(cutoff_nFeature_lower, cutoff_nFeature_upper), color = "red") +
#      theme_classic() + ggtitle(paste(sample, prefix, "nFeature_RNA Distribution")),
#    file.path(output_dir, paste0(sample, "_QC_", prefix, "_nFeature_density.png"))
#  )

#  save_plot(
#    ggplot(meta, aes(x = percent.mt)) +
#      geom_density(fill = "green", alpha = 0.4) +
#      geom_vline(xintercept = cutoff_percent_mt, color = "red") +
#      theme_classic() + ggtitle(paste(sample, prefix, "Mitochondrial %")),
#    file.path(output_dir, paste0(sample, "_QC_", prefix, "_mito_ratio.png"))
#  )

#  save_plot(
#    ggplot(meta, aes(x = log10GenesPerUMI)) +
#      geom_density(fill = "purple", alpha = 0.4) +
#      geom_vline(xintercept = cutoff_complexity, color = "red") +
#      theme_classic() + ggtitle(paste(sample, prefix, "Gene Complexity")),
#    file.path(output_dir, paste0(sample, "_QC_", prefix, "_gene_complexity.png"))
#  )

# save_plot(
#  ggplot(meta, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
#    geom_point(size = 0.5) +
#    scale_colour_gradient(low = "gray90", high = "black") +
#    stat_smooth(method = "lm") +
#    scale_x_log10() + scale_y_log10() +
#    geom_hline(yintercept = cutoff_nFeature_lower, color = "red") +
#    geom_hline(yintercept = cutoff_nFeature_upper, color = "red") +
#    theme_classic() + ggtitle(paste(sample, "nCount vs nFeature")),
#file.path(output_dir, paste0(sample, "_QC_", prefix, "_umi_vs_gene.png"))
#)
#}

make_qc_plots <- function(meta_list) {
  # ▶ 종류별 플롯 저장
  # nFeature_RNA
  p1 <- lapply(names(meta_list), function(name) {
    ggplot(meta_list[[name]], aes(x = nFeature_RNA)) +
      geom_density(fill = "orange", alpha = 0.4) + scale_x_log10() +
      geom_vline(xintercept = c(cutoff_nFeature_lower, cutoff_nFeature_upper), color = "red") +
      theme_classic() + ggtitle(name)
  }) %>% wrap_plots(nrow = 1) + plot_annotation(title = paste(sample, "nFeature_RNA Distribution"))
  save_plot(p1, file.path(output_dir, paste0(sample, "_QC_nFeature_all.png")))

  # percent.mt
  p2 <- lapply(names(meta_list), function(name) {
    ggplot(meta_list[[name]], aes(x = percent.mt)) +
      geom_density(fill = "green", alpha = 0.4) +
      geom_vline(xintercept = cutoff_percent_mt, color = "red") +
      theme_classic() + ggtitle(name)
  }) %>% wrap_plots(nrow = 1) + plot_annotation(title = paste(sample, "Mitochondrial %"))
  save_plot(p2, file.path(output_dir, paste0(sample, "_QC_mito_ratio_all.png")))

  # complexity
  p3 <- lapply(names(meta_list), function(name) {
    ggplot(meta_list[[name]], aes(x = log10GenesPerUMI)) +
      geom_density(fill = "purple", alpha = 0.4) +
      geom_vline(xintercept = cutoff_complexity, color = "red") +
      theme_classic() + ggtitle(name)
  }) %>% wrap_plots(nrow = 1) + plot_annotation(title = paste(sample, "Gene Complexity"))
  save_plot(p3, file.path(output_dir, paste0(sample, "_QC_gene_complexity_all.png")))

  # scatter
  p4 <- lapply(names(meta_list), function(name) {
    ggplot(meta_list[[name]], aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
      geom_point(size = 0.5) +
      scale_colour_gradient(low = "gray90", high = "black") +
      stat_smooth(method = "lm") +
      scale_x_log10() + scale_y_log10() +
      geom_hline(yintercept = cutoff_nFeature_lower, color = "red") +
      geom_hline(yintercept = cutoff_nFeature_upper, color = "red") +
      theme_classic() + ggtitle(name)
  }) %>% wrap_plots(nrow = 1) + plot_annotation(title = paste(sample, "nCount vs nFeature"))
  save_plot(p4, file.path(output_dir, paste0(sample, "_QC_umi_vs_gene_all.png")))
}


# ▶ 데이터 로딩 및 QC 계산
seurat_obj <- Read10X_h5(h5_path) %>%
  CreateSeuratObject(project = sample, min.cells = 3, min.features = 200)

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["log10GenesPerUMI"]] <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
initial_cell_n <- ncol(seurat_obj)

meta_list <- list(
  before_filter = seurat_obj@meta.data
)


# ▶ QC 플롯 (before filter)
#make_qc_plots(seurat_obj@meta.data, "before_filter")

# ▶ QC Filtering
seurat_obj <- subset(seurat_obj, subset = (
  nFeature_RNA >= cutoff_nFeature_lower &
  nFeature_RNA <= cutoff_nFeature_upper &
  percent.mt < cutoff_percent_mt &
  log10GenesPerUMI > cutoff_complexity
))
meta_list$after_filter <- seurat_obj@meta.data
filtered_cell_n <- ncol(seurat_obj)

# ▶ QC 플롯 (after filter)
#make_qc_plots(seurat_obj@meta.data, "after_filter")

# ▶ Seurat 전처리
seurat_obj <- NormalizeData(seurat_obj) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# ▶ DoubletFinder
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj@assays$RNA <- as(Class = "Assay", seurat_obj@assays$RNA)
sweep_res <- paramSweep_v3(seurat_obj, PCs=1:10, sct=FALSE)
sweep_stats <- summarizeSweep(sweep_res, GT=FALSE)
pk <- find.pK(sweep_stats)$pK[which.max(find.pK(sweep_stats)$BCmetric)]
nExp <- round(ncol(seurat_obj) * expected_rate)

seurat_obj <- doubletFinder_v3(seurat_obj, PCs=1:10, pN=0.25,
                                pK=as.numeric(as.character(pk)),
                                nExp=nExp, reuse.pANN=FALSE, sct=FALSE)

doublet_col <- grep("DF.classifications", colnames(seurat_obj@meta.data), value=TRUE)
seurat_obj$doublet_status <- seurat_obj[[doublet_col]][,1]

# ▶ QC 플롯 (after doublet)
#make_qc_plots(seurat_obj@meta.data, "after_filter_doublet")
meta_list$after_filter_doublet <- seurat_obj@meta.data
# ▶ 모든 단계 QC plot 한 번에
make_qc_plots(meta_list)

# ▶ UMAP 저장
umap_p <- DimPlot(seurat_obj, reduction = "umap", group.by = "doublet_status",
                  cols = c("gray70", "red")) +
  ggtitle(paste("DoubletFinder:", sample))

ggsave(output_doublet_umap, umap_p, width = 7, height = 6, dpi = 300)
system(paste("chmod 777", shQuote(output_doublet_umap)))

# ▶ RDS 저장
saveRDS(seurat_obj, output_rds)
system(paste("chmod 777", shQuote(output_rds)))

# ▶ 로그 작성
filter_ratio <- 1 - (filtered_cell_n / initial_cell_n)
flag <- ifelse(filter_ratio > 0.6, "⚠️ Potential over-filtering", "")
log_lines <- c(
  paste0("Sample: ", sample),
  paste0("Input cells: ", initial_cell_n),
  sprintf("Filtered cells: %d (%.1f%% removed) %s", filtered_cell_n, filter_ratio * 100, flag),
  "",
  "QC Filtering Parameters:",
  paste0("Manual nFeature_RNA cutoff: ", cutoff_nFeature_lower, " ~ ", cutoff_nFeature_upper),
  paste0("Mitochondrial % cutoff: ", cutoff_percent_mt),
  paste0("Gene Complexity cutoff: ", cutoff_complexity),
  "",
  "DoubletFinder Parameters:",
  paste0("Best pK: ", pk),
  paste0("Expected doublet rate: ", expected_rate),
  paste0("Expected doublets: ", nExp)
)
writeLines(log_lines, con = output_log)
system(paste("chmod 777", shQuote(output_log)))


