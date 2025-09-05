#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(utils)
})

# ---------- CLI Arguments ----------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input directory with GFF3 files"),
  make_option(c("-b", "--blast"), type = "character", help = "Input directory with BLAST format6 files"),
  make_option(c("-o", "--output"), type = "character", help = "Output image file (e.g., out.svg)"),
  make_option(c("-l", "--list"), type = "character", help = "File with desired genome ID order (no extensions)"),
  make_option(c("--output_format"), type = "character", default = "svg", help = "Output format (svg|png|pdf) [default: %default]"),
  make_option(c("--distance"), type = "double", default = 0.7, help = "Vertical spacing between genomes [default: %default]"),
  make_option(c("--show_legends"), action = "store_true", default = FALSE, help = "Show gene legend"),
  make_option(c("--labels_on_top"), action = "store_true", default = FALSE, help = "Show gene names on top of arrows"),
  make_option("--use_product", action = "store_true", default = FALSE, help = "Use 'product' instead of 'gene'")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$output)) {
  stop("Please specify input GFF dir (-i) and output file (-o)")
}

# ---------- Read Order List ----------
file_order <- NULL
if (!is.null(opt$list)) {
  if (!file.exists(opt$list)) stop("Provided order list file does not exist.")
  file_order <- read_lines(opt$list) %>% str_trim() %>% discard(~ .x == "")
  file_order <- rev(file_order)
}

# ---------- Load GFF ----------
gff_files <- list.files(opt$input, pattern = "\\.gff$", full.names = TRUE)

# Modify the parse_gff function to extract 'product'
parse_gff <- function(file) {
  raw_lines <- readLines(file)
  fasta_start <- which(grepl("^##FASTA", raw_lines))[1]
  if (!is.na(fasta_start)) {
    raw_lines <- raw_lines[1:(fasta_start - 1)]
  }
  parsed <- read_tsv(paste(raw_lines, collapse = "\n"), comment = "#", col_names = FALSE, show_col_types = FALSE)
  parsed %>%
    filter(X3 == "CDS") %>%
    mutate(
      contig = X1,
      start = as.numeric(X4),
      end = as.numeric(X5),
      strand = X7,
      gene = str_extract(X9, "gene=([^;]+)") %>% str_replace("gene=", "") %>% replace_na("unknown") %>% URLdecode(),
      product = str_extract(X9, "product=([^;]+)") %>% str_replace("product=", "") %>% replace_na("unknown") %>% URLdecode(),
      id = str_extract(X9, "ID=([^;]+)") %>% str_replace("ID=", ""),
      file = basename(file),
      file_id = tools::file_path_sans_ext(file)
    ) %>%
    select(file, file_id, contig, start, end, strand, gene, product, id)
}

# Determine which column to use based on the flag
all_genes_raw <- map_dfr(gff_files, parse_gff) %>%
  mutate(label = if (opt$use_product) product else gene)

# If order list is provided, arrange files accordingly
if (!is.null(file_order)) {
  all_genes_raw <- all_genes_raw %>%
    mutate(file_id = factor(file_id, levels = file_order)) %>%
    arrange(file_id)
}

# Update unique_genes and other references to use 'label'
unique_genes <- sort(unique(all_genes_raw$label))
all_files <- unique(all_genes_raw$file_id)
full_grid <- expand_grid(file_id = all_files, label = unique_genes)

# Update gene_lengths_reference to use 'label'
gene_lengths_reference <- all_genes_raw %>%
  arrange(file_id) %>%
  group_by(label) %>%
  slice_head(n = 1) %>%
  mutate(gene_length = if_else(!is.na(end) & !is.na(start), (end - start) * 0.25, 200)) %>%
  select(label, gene_length)

# Merge lengths and build layout
all_genes <- full_grid %>%
  left_join(all_genes_raw, by = c("file_id", "label")) %>%
  group_by(file_id, label) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  left_join(gene_lengths_reference, by = "label") %>%
  mutate(y = as.numeric(factor(file_id))) %>%
  arrange(file_id, label) %>%
  group_by(file_id) %>%
  mutate(
    spacing = 50,
    pseudo_start = cumsum(lag(gene_length + spacing, default = 0)),
    pseudo_end = pseudo_start + gene_length,
    pseudo_mid = (pseudo_start + pseudo_end) / 2,
    gene_id = paste0(label, "_", row_number())
  ) %>%
  ungroup()

# Set horizontal shift
x_offset <- 15000  # or adjust as needed

# Shift gene positions
all_genes <- all_genes %>%
  mutate(
    pseudo_start = pseudo_start + x_offset,
    pseudo_end = pseudo_end + x_offset,
    pseudo_mid = pseudo_mid + x_offset
  )

# Shift genome labels
file_labels <- all_genes %>%
  group_by(file_id, y) %>%
  summarise(
    label = first(file_id),
    x = min(pseudo_start) - 100,
    .groups = "drop"
  )

# Shift center lines
center_lines <- all_genes %>%
  group_by(file_id, y) %>%
  summarise(
    x_start = min(pseudo_start) - 100,
    x_end = max(pseudo_end) + 100,
    .groups = "drop"
  )

# Gene label coordinates
gene_labels <- all_genes %>%
  select(file_id, gene_id, label, pseudo_mid, y)

# Gene arrow polygons
feature_height <- 0.3
gene_polys <- all_genes %>%
  rowwise() %>%
  mutate(
    strand = "+",
    length = pseudo_end - pseudo_start,
    body_length = length * 0.75,
    head_length = length * 0.25,
    body_start = ifelse(strand == "+", pseudo_start, pseudo_end - body_length),
    body_end = ifelse(strand == "+", pseudo_start + body_length, pseudo_end),
    tip_x = ifelse(strand == "+", pseudo_end, pseudo_start),
    y0 = y,
    poly = list(tibble(
      x_poly = if (strand == "+") {
        c(pseudo_start, body_end, tip_x, body_end, pseudo_start)
      } else {
        c(pseudo_end, body_start, tip_x, body_start, pseudo_end)
      },
      y_poly = c(y0 - feature_height, y0 - feature_height, y0, y0 + feature_height, y0 + feature_height)
    )),
    has_data = !is.na(start)
  ) %>%
  unnest(poly) %>%
  mutate(poly_id = paste(file_id, gene_id, sep = "_"))

# Calculate outside the ggplot
min_x_label <- min(file_labels$x, na.rm = TRUE)

# ---- Plot ----
p <- ggplot() +
  geom_text(data = file_labels,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5, fontface = "bold", size = 3) +
  geom_segment(data = center_lines,
               aes(x = x_start, xend = x_end, y = y, yend = y),
               color = "black", linewidth = 0.6) +
  geom_polygon(data = gene_polys %>% filter(has_data),
               aes(x = x_poly, y = y_poly, group = poly_id, fill = label),
               color = "black", show.legend = opt$show_legends) +
  geom_polygon(data = gene_polys %>% filter(!has_data),
               aes(x = x_poly, y = y_poly, group = poly_id),
               color = "black", fill = "white") +
  {
    if (opt$labels_on_top) {
      geom_text(data = gene_labels %>% filter(file_id == rev(file_order)[1]),
                aes(x = pseudo_mid, y = y + 0.5, label = label),
                size = 2.5, angle = 45, hjust = 0, vjust = 0.5)
    } else NULL
  } +
  scale_fill_discrete(name = "Genes") +
  # coord_cartesian(
  #   xlim = c(min(min_x_label, min(gene_polys$x_poly, na.rm = TRUE)) - 100,
  #            max(gene_polys$x_poly, na.rm = TRUE) + 500),
  #   ylim = c(0.5, max(gene_polys$y) + 1)
  # ) +
  theme_void() +
  coord_cartesian(clip = "off") + 
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(60, 20, 20, 40),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 7),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.5, "cm"),
    legend.key.width = unit(1.0, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(0.1, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0)
  ) +
  guides(fill = guide_legend(
    ncol = 2,
    byrow = TRUE,
    override.aes = list(color = "black", size = 0.2)              # Thinner borders
  ))

ggsave(opt$output, plot = p,
       width = 10,
       height = 1 + opt$distance * n_distinct(all_genes$file_id),
       dpi = 300,
       device = opt$output_format,
       limitsize = FALSE)

