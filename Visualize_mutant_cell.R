


library(Seurat)
library(ggplot2)


aml_data <- readRDS("/data/workbench/scRSEQ_AML/data/aml.test.bd.mut3_DNMT3A_mut_meta_added.Rds")

head(aml_data@meta.data)
table(aml_data@meta.data$Mut_info)
table(aml_data@meta.data$Allele)

Reductions(aml_data)

# UMAP 좌표에 있는 셀 ID
head(rownames(Embeddings(aml_data, "umap.rpca2")))

cells_in_umap <- rownames(Embeddings(aml_data, reduction = "umap.rpca2"))
head(cells_in_umap)
length(cells_in_umap)
nrow(aml_data@meta.data)

cells_with_allele <- rownames(aml_data@meta.data)[!is.na(aml_data@meta.data$Allele)]
cells_allele_in_umap <- intersect(cells_in_umap, cells_with_allele)
length(cells_allele_in_umap)


# Allele 있는 셀 ID 중 몇 개가 UMAP에 포함되는가?
sum(!is.na(aml_data@meta.data$Allele))                             # 전체 변이 셀 수
sum(!is.na(aml_data@meta.data$Allele) & rownames(aml_data@meta.data) %in% rownames(Embeddings(aml_data, "umap.rpca2")))  # UMAP에 포함된 변이 셀 수


DefaultAssay(aml_data) <- "RNA"

# # UMAP plot (색상 지정)
# DimPlot(aml_data, group.by = "Allele", reduction = "umap.rpca2") +
#   scale_color_manual(
#     values = c(
#       "REF" = "green3",
#       "ALT" = "red",
#       "Hetero" = "orange"
#     )
#   ) +
#   ggtitle("UMAP by Allele (DNMT3A mutation)") +
#   theme_minimal()


# 1. 변이 있는 셀 목록
cells_with_allele <- colnames(aml_data)[!is.na(aml_data@meta.data$Allele)]
cells_with_allele

length(cells_with_allele)


cells_ALT <- WhichCells(aml_data, expression = Allele == "ALT")
cells_REF <- WhichCells(aml_data, expression = Allele == "REF")
cells_HET <- WhichCells(aml_data, expression = Allele == "Hetero")

length(cells_ALT)
length(cells_REF)
length(cells_HET)


DimPlot(aml_data, reduction = "umap.rpca2",
        cells.highlight = list(ALT = cells_ALT, Hetero = cells_HET, REF = cells_REF),
        cols.highlight = c("green3", "orange", "red"),
        sizes.highlight = 1.5,
        pt.size = 0.2) +
  ggtitle("UMAP: DNMT3A mutation highlight") +
  theme_minimal()


# UMAP 좌표 + 메타데이터 결합

umap_df <- as.data.frame(Embeddings(aml_data, reduction = "umap.rpca2")) %>%
  tibble::rownames_to_column("cell_id") %>%
  left_join(aml_data@meta.data %>% select(Allele, final.clus) %>% tibble::rownames_to_column("cell_id"),
            by = "cell_id")

# 1. NA(회색 점)과 컬러 점을 분리
umap_df_gray <- umap_df %>% filter(is.na(Allele))
umap_df_colored <- umap_df %>% filter(!is.na(Allele))

# 2. ggplot: 회색 먼저, 컬러 점 나중에!
ggplot() +
  geom_point(data = umap_df_gray, aes(x = umaprpca2_1, y = umaprpca2_2),
             color = "gray85", size = 0.3, alpha = 0.5) +
  geom_point(data = umap_df_colored,
             aes(x = umaprpca2_1, y = umaprpca2_2, color = Allele),
             size = 1.2, alpha = 0.9) +
  scale_color_manual(values = c("REF" = "green3", "ALT" = "red", "Hetero" = "orange")) +
  # Cell type 텍스트 추가
  geom_text(data = umap_df %>%
              filter(!is.na(final.clus)) %>%
              group_by(final.clus) %>%
              summarise(x = median(umaprpca2_1), y = median(umaprpca2_2)),
            aes(x = x, y = y, label = final.clus),
            size = 3.5, fontface = "bold") +
  theme_minimal() +
  ggtitle("UMAP with DNMT3A Allele and Cell Type (Colored points on top)") +
  theme(legend.position = "right")

table(aml_data@meta.data$Allele, useNA = "ifany")

allele_counts <- table(aml_data@meta.data$Allele)
allele_counts




# UMAP 좌표 + 메타데이터 병합
umap_df <- as.data.frame(Embeddings(aml_data, reduction = "umap.rpca2")) %>%
  tibble::rownames_to_column("cell_id") %>%
  left_join(aml_data@meta.data %>% select(Allele, final.clus, Mut_info) %>%
              tibble::rownames_to_column("cell_id"),
            by = "cell_id")

# T cell 클러스터만 추출
tcell_df <- aml_data@meta.data %>%
  dplyr::filter(final.clus %in% c("CD4_T_cell", "CD8_T_cell")) %>%
  tibble::rownames_to_column("cell_id")

# UMAP 좌표와 병합
umap_tcell <- Embeddings(aml_data, reduction = "umap.rpca2") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::inner_join(tcell_df, by = "cell_id")

umap_tcell_label <- umap_tcell %>%
  filter(!is.na(Mut_info))  # 변이 있는 셀만 라벨 표시

library(ggplot2)
library(ggrepel)

ggplot(umap_tcell, aes(x = umaprpca2_1, y = umaprpca2_2)) +
  geom_point(color = "gray70", size = 0.8) +
  geom_point(data = umap_tcell_label, aes(color = Allele), size = 2) +
  geom_text_repel(data = umap_tcell_label,
                  aes(label = Mut_info, color = Allele),
                  size = 3.5,
                  max.overlaps = 50,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.size = 0.2,
                  segment.color = "gray50",
                  show.legend = FALSE) +
  scale_color_manual(values = c("REF" = "green3", "ALT" = "red", "Hetero" = "orange")) +
  ggtitle("T cell 영역 내 DNMT3A 변이 셀 라벨링 (with ggrepel)") +
  theme_minimal()

# T cell 클러스터 중 Mut_info 있는 셀만 필터링
tcell_mut_info_count <- aml_data@meta.data %>%
  dplyr::filter(final.clus %in% c("CD4_T_cell", "CD8_T_cell"),
                !is.na(Mut_info)) %>%
  dplyr::count(Mut_info, sort = TRUE)

# 결과 출력
tcell_mut_info_count


