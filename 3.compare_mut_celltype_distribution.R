#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("샘플명을 인자로 제공해주세요. 예: Rscript compare_celltype_distribution.R 10396")
}
sample <- args[1]

mutation_info <- list(
  "8825"  = list(chr = "chr2", pos = "25234373", ref = "C", alt = "T"),
  "10396" = list(chr = "chr2", pos = "25234373", ref = "C", alt = "T"),
  "10090" = list(chr = "chr2", pos = "25234374", ref = "G", alt = "A"),
  "9740"  = list(chr = "chr2", pos = "25240672", ref = "C", alt = "G")
)

info <- mutation_info[[sample]]
base_dir <- "/data/processed_data/scRSEQ_AML/HN00247118/mutation/2020"

#filtered_labels_path <- file.path(base_dir, paste0(sample, "_filtered_singler_result.rds"))
filtered_labels_path <- file.path(base_dir, paste0(sample, "_singler_result.rds"))
ref_file <- file.path(base_dir, paste0(sample, "_sample_alignments.", info$chr, ".", info$pos, ".", info$ref, ".DNMT3A.ref.bam.cell_barcodes.umi.tsv"))
var_file <- file.path(base_dir, paste0(sample, "_sample_alignments.", info$chr, ".", info$pos, ".", info$alt, ".DNMT3A.var.bam.cell_barcodes.umi.tsv"))

filtered_labels <- readRDS(filtered_labels_path)
ref_barcodes <- read.delim(ref_file, header = TRUE)$CellBarcode
var_barcodes <- read.delim(var_file, header = TRUE)$CellBarcode

if (is.null(names(filtered_labels))) {
  stop("filtered_labels has no names (barcodes). Cannot proceed.")
}
singler_df <- data.frame(CellBarcode = names(filtered_labels), CellType = filtered_labels)

ref_df <- singler_df %>% filter(CellBarcode %in% ref_barcodes) %>% mutate(Group = "Ref")
var_df <- singler_df %>% filter(CellBarcode %in% var_barcodes) %>% mutate(Group = "Variant")

# ▶ 요약 테이블 생성
plot_df_wide <- bind_rows(ref_df, var_df) %>%
  group_by(Group, CellType) %>%
  summarise(Count = n(), .groups = "drop") %>%
  complete(Group, CellType, fill = list(Count = 0)) %>%
  pivot_wider(names_from = Group, values_from = Count, values_fill = 0) %>%
  mutate(
    Total = Ref + Variant,
    Ratio = ifelse(Total > 0, Variant / Total, NA),
    RatioPercent = round(Ratio * 100, 1),
    Label = ifelse(Total > 0, paste0(Variant, " / ", Total, " (", RatioPercent, "%)"), "")
  )

# ▶ 전체 합계 요약
grand_total_ref <- sum(plot_df_wide$Ref)
grand_total_var <- sum(plot_df_wide$Variant)
grand_total <- grand_total_ref + grand_total_var
grand_ratio <- ifelse(grand_total > 0, grand_total_var / grand_total, NA)
grand_ratio_percent <- round(grand_ratio * 100, 1)
grand_label <- paste0("Total: ", grand_total_var, " / ", grand_total, " (", grand_ratio_percent, "%)")

# ▶ 저장용 CSV
csv_out <- file.path(base_dir, paste0(sample, "_DNMT3A_", info$chr, "_", info$pos, "_celltype_distribution.csv"))
write_csv(plot_df_wide %>% select(CellType, Ref, Variant, Total, RatioPercent, Label), csv_out)
system(paste("chmod 777", shQuote(csv_out)))

# ▶ plot용 long format
plot_df_long <- plot_df_wide %>%
  pivot_longer(cols = c("Ref", "Variant"), names_to = "Group", values_to = "Count")

# ▶ label은 Variant만
text_df <- plot_df_long %>%
  filter(Group == "Variant") %>%
  mutate(label_y = Count + max(Count) * 0.03)

# ▶ plot 생성
plot_title <- paste0("Cell Type Distribution — Sample ", sample,
                     " (", info$chr, ":", info$pos, " ", info$ref, ">", info$alt, ")\n", grand_label)

output_file <- file.path(base_dir, paste0(sample, "_DNMT3A_", info$chr, "_", info$pos, "_celltype_distribution.png"))

p <- ggplot(plot_df_long, aes(x = CellType, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(data = text_df, aes(y = label_y, label = Label), vjust = 0, size = 3.2) +
  scale_fill_manual(values = c("Ref" = "green", "Variant" = "red")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(plot_title) +
  ylab("Cell Count") +
  xlab("Cell Type")

ggsave(output_file, plot = p, width = 10, height = 6, dpi = 300)
system(paste("chmod 777", shQuote(output_file)))
