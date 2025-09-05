#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(utils)
  library(cowplot) # For combining plots
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
  make_option(c("--char"), type = "integer", default = NULL, help = "Maximum number of characters for product names")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$output) || is.null(opt$blast)) {
  stop("Please specify input GFF dir (-i), blast folder (-b), and output file (-o)")
}

# ---------- Read Order List ----------
file_order <- NULL
if (!is.null(opt$list)) {
  if (!file.exists(opt$list)) stop("Provided order list file does not exist.")
  file_order <- readr::read_lines(opt$list) %>% str_trim() %>% discard(~ .x == "")
  file_order <- rev(file_order)  # Top-to-bottom
}

# ---------- Load GFF ----------
gff_files <- list.files(opt$input, pattern = "\\.gff$", full.names = TRUE)

parse_gff <- function(file) {
  raw_lines <- readLines(file)
  fasta_start <- which(grepl("^##FASTA", raw_lines))[1]
  if (!is.na(fasta_start)) {
    raw_lines <- raw_lines[1:(fasta_start - 1)]
  }

  parsed <- readr::read_tsv(paste(raw_lines, collapse = "\n"), comment = "#", col_names = FALSE, show_col_types = FALSE)

  parsed %>%
    dplyr::filter(X3 == "CDS") %>%
    dplyr::mutate(
      contig = X1,
      start = as.numeric(X4),
      end = as.numeric(X5),
      strand = X7,
      product = str_extract(X9, "product=([^;]+)") %>%
        str_replace("product=", "") %>%
        replace_na("unknown") %>%
        URLdecode() %>%
        str_replace_all(" ", "_"),  # Replace spaces with underscores
      id = str_extract(X9, "ID=([^;]+)") %>% str_replace("ID=", ""),
      file = basename(file),
      file_id = tools::file_path_sans_ext(file)
    ) %>%
    dplyr::select(file, file_id, contig, start, end, strand, product, id)
}

all_genes <- purrr::map_dfr(gff_files, parse_gff)

# ---------- Keep full product and make display version (honors --char) ----------
all_genes <- all_genes %>%
  mutate(
    product_full = product,                                               # full cleaned product for logic/legend
    product_disp = if (!is.null(opt$char)) str_sub(product, 1, opt$char) else product,
    is_hypo = str_detect(product_full, regex("^hypothetical(_|\\s)*protein$", ignore_case = TRUE))
  )

# ---------- Reorder Based on List ----------
if (!is.null(file_order)) {
  all_genes <- all_genes %>%
    mutate(file_id = factor(file_id, levels = file_order)) %>%
    arrange(file_id)
}

# Assign a unique vertical position (y) for each file_id
all_genes <- all_genes %>%
  mutate(y = as.numeric(factor(file_id)))

# Calculate max_y for tree scaling
max_y <- max(all_genes$y)

# ---------- Arrange Contigs by coordinate order ----------
all_genes <- all_genes %>%
  mutate(contig = factor(contig, levels = unique(contig)[order(as.numeric(str_extract(unique(contig), "\\d+")))]))

# Gene IDs (use full product so truncation doesn't affect IDs)
all_genes <- all_genes %>%
  group_by(file_id, contig, start, end, product_full) %>%
  mutate(gene_id = paste0(product_full, "_", cur_group_id())) %>%
  ungroup()

# Sort contigs naturally (again, safe)
all_genes <- all_genes %>%
  mutate(contig = factor(contig, levels = unique(contig)[order(as.numeric(str_extract(unique(contig), "\\d+")))]))

# Block detection
all_genes <- all_genes %>%
  group_by(file_id, contig) %>%
  arrange(start) %>%
  mutate(
    gap = start - lag(end, default = first(start)),
    block = cumsum(if_else(gap > 1000, 1, 0))
  ) %>%
  ungroup()

# Treat all blocks (regardless of contig) equally in offset spacing
spacing <- 600

block_sizes <- all_genes %>%
  group_by(file_id, contig, block) %>%
  summarise(block_width = sum(end - start), .groups = "drop") %>%
  arrange(file_id, contig, block) %>%
  group_by(file_id) %>%
  mutate(offset = cumsum(lag(block_width + spacing, default = 0)))

# Pseudo coords
all_genes <- all_genes %>%
  left_join(block_sizes, by = c("file_id", "contig", "block")) %>%
  group_by(file_id, contig, block) %>%
  arrange(start) %>%
  mutate(
    pseudo_start = cumsum(lag(end - start + 1, default = 0)) + 1 + offset,
    pseudo_end = pseudo_start + (end - start),
    pseudo_mid = (pseudo_start + pseudo_end) / 2
  ) %>%
  ungroup()

# Center lines: break around blocks
center_lines <- all_genes %>%
  group_by(file_id, y, contig, block) %>%
  summarise(
    x_start = min(pseudo_start) - 200,
    x_end = max(pseudo_end) + 200,
    .groups = "drop"
  )

# Identify inter-block breaks for labeling
block_transitions <- center_lines %>%
  group_by(file_id, y) %>%
  arrange(x_start) %>%
  mutate(next_contig = lead(contig),
         next_x_start = lead(x_start)) %>%
  filter(!is.na(next_x_start)) %>%
  mutate(
    symbol = if_else(contig != next_contig, "//", "."),
    x = (x_end + next_x_start) / 2
  ) %>%
  select(file_id, y, x, symbol)

block_transitions <- block_transitions %>%
  mutate(y_adj = case_when(
    symbol == "." ~ y + 0.09,
    symbol == "//" ~ y + 0.05
  ))

# File labels
file_labels <- all_genes %>%
  distinct(file_id, y) %>%
  mutate(label = file_id, x = -500)

# Gene labels (precompute label; NA hides text; display is truncated)
gene_labels <- all_genes %>%
  transmute(file_id, gene_id, pseudo_mid, y,
            label = if_else(is_hypo, NA_character_, product_disp))

# Arrows
feature_height <- 0.3

gene_polys <- all_genes %>%
  rowwise() %>%
  mutate(
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
    ))
  ) %>%
  unnest(poly) %>%
  mutate(poly_id = paste(file_id, gene_id, sep = "_")) %>%
  ungroup()

# Read BLAST
blast_files <- list.files(opt$blast, pattern = "\\.\\w+$", full.names = TRUE)

blast_df <- purrr::map_dfr(blast_files, readr::read_tsv, col_names = FALSE, show_col_types = FALSE)

colnames(blast_df) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                        "qstart", "qend", "sstart", "send", "evalue", "bitscore")

hits <- blast_df %>%
  left_join(select(all_genes, id, pseudo_start, pseudo_end, y, file_id), by = c("qseqid" = "id")) %>%
  rename(q_start = pseudo_start, q_end = pseudo_end, q_y = y, q_file = file_id) %>%
  left_join(select(all_genes, id, pseudo_start, pseudo_end, y, file_id), by = c("sseqid" = "id")) %>%
  rename(s_start = pseudo_start, s_end = pseudo_end, s_y = y, s_file = file_id) %>%
  filter(!is.na(q_start) & !is.na(s_start))

blast_polys <- hits %>%
  mutate(poly_id = paste(qseqid, sseqid, sep = "_")) %>%
  rowwise() %>%
  mutate(poly = list(tibble(
    x = c(q_start, q_end, s_end, s_start),
    y = c(q_y - feature_height, q_y - feature_height,
          s_y + feature_height, s_y + feature_height)
  ))) %>%
  unnest(poly) %>%
  ungroup()

# ---------- Plot ----------
p <- ggplot() +
  geom_text(data = file_labels,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5, fontface = "bold", size = 3) +
  geom_segment(data = center_lines,
               aes(x = x_start, xend = x_end, y = y, yend = y),
               color = "black", linewidth = 0.6) +
  geom_polygon(data = blast_polys,
               aes(x = x, y = y, group = poly_id),
               fill = "grey80", alpha = 0.4, color = NA) +
  geom_polygon(data = gene_polys,
               aes(x = x_poly, y = y_poly, group = poly_id, fill = product_full),
               color = "black", show.legend = opt$show_legends) +
  geom_text(data = block_transitions,
            aes(x = x, y = y_adj, label = symbol),
            size = 4.5, fontface = "bold", vjust = 0.5) +
  {
    if (opt$labels_on_top) {
      geom_text(
        data = gene_labels,
        aes(x = pseudo_mid, y = y + 0.4, label = label),
        size = 2.5, angle = 45, hjust = 0, vjust = 0.5,
        na.rm = TRUE
      )
    } else NULL
  } +
  scale_fill_discrete(name = "Products") +
  coord_cartesian(clip = "off") +                 # prevent clipping of labels outside panel
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(face = "bold"),
    legend.key.width = unit(1.5, "cm"),
    plot.margin = margin(10, 10, 50, 10)         # larger left margin to avoid truncation
  ) +
  guides(fill = guide_legend(
    ncol = 1,
    override.aes = list(color = "black", size = 0.3)
  ))

# Save the plot
ggsave(opt$output, plot = p,
       width = 18,
       height = 1 + opt$distance * n_distinct(all_genes$file_id),
       dpi = 300,
       device = opt$output_format)
