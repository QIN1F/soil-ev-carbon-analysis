#!/usr/bin/env Rscript

# =============================================================================
# Publication-quality rectangular MAG phylogeny with two aligned bar tracks
#
# Required files in one folder:
#   1) MAG-result.csv
#   2) unroot.tree
#   3) one FASTA assembly per target MAG: *.fna, *.fa, or *.fasta
#
# Figure layout:
#   Left:  rectangular MAG phylogeny and phylum-coloured MAG label band
#   Right: Genome size (Mb) and Carbon_Gene_Total (0-200) horizontal bars
#   Bars:  no per-MAG numeric labels; both quantitative axes start at zero
#
# Dependency: ape + base R only
# =============================================================================

options(stringsAsFactors = FALSE)
SCRIPT_VERSION <- "1.0-rectangular-bars-2026-09-02"
message("MAG rectangular phylogeny script: ", SCRIPT_VERSION)

# ---- 0. Load the only non-base dependency ----------------------------------
if (!requireNamespace("ape", quietly = TRUE)) {
  message("Package 'ape' is not installed. Installing it from CRAN ...")
  try(
    install.packages("ape", repos = "https://cloud.r-project.org", dependencies = TRUE),
    silent = TRUE
  )
}
if (!requireNamespace("ape", quietly = TRUE)) {
  stop("Package 'ape' is required. Run: install.packages('ape')", call. = FALSE)
}
suppressPackageStartupMessages(library(ape))

# ---- 1. Select the input folder --------------------------------------------
# Optional: enter an absolute folder path here. Leave as "" to select the
# folder interactively in RStudio.
INPUT_DIR_OVERRIDE <- ""

has_required_inputs <- function(path) {
  is.character(path) && length(path) == 1L && dir.exists(path) &&
    file.exists(file.path(path, "MAG-result.csv")) &&
    file.exists(file.path(path, "unroot.tree"))
}

script_directory <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", command, value = TRUE)
  if (length(hit) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE)))
  }
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (!is.null(frames[[i]]$ofile)) {
      return(dirname(normalizePath(frames[[i]]$ofile, mustWork = FALSE)))
    }
  }
  getwd()
}

locate_input_dir <- function() {
  if (nzchar(trimws(INPUT_DIR_OVERRIDE))) {
    candidate <- normalizePath(trimws(INPUT_DIR_OVERRIDE), mustWork = FALSE)
    if (!has_required_inputs(candidate)) {
      stop(
        "INPUT_DIR_OVERRIDE must contain MAG-result.csv and unroot.tree: ",
        candidate,
        call. = FALSE
      )
    }
    return(candidate)
  }

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1L) {
    candidate <- normalizePath(args[1], mustWork = FALSE)
    if (!has_required_inputs(candidate)) {
      stop(
        "The supplied input folder must contain MAG-result.csv and unroot.tree: ",
        candidate,
        call. = FALSE
      )
    }
    return(candidate)
  }

  if (interactive()) {
    if (.Platform$OS.type == "windows") {
      candidate <- utils::choose.dir(
        default = getwd(),
        caption = "Select the folder containing MAG-result.csv, unroot.tree, and FASTA files"
      )
      if (is.na(candidate) || !nzchar(candidate)) {
        stop("No input folder was selected.", call. = FALSE)
      }
    } else {
      message("Select MAG-result.csv; its containing folder will be used.")
      candidate <- dirname(file.choose())
    }
    candidate <- normalizePath(candidate, mustWork = TRUE)
    if (!has_required_inputs(candidate)) {
      stop("The selected folder lacks MAG-result.csv or unroot.tree.", call. = FALSE)
    }
    return(candidate)
  }

  script_dir <- script_directory()
  candidates <- unique(c(
    getwd(), script_dir, file.path(getwd(), "upload"),
    file.path(dirname(script_dir), "upload"),
    file.path(dirname(dirname(script_dir)), "upload")
  ))
  for (candidate in candidates) {
    if (has_required_inputs(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  stop(
    paste0(
      "Input folder not found. Run with:\n",
      "Rscript draw_MAG_phylogeny_rectangular_bars.R /path/to/input_folder"
    ),
    call. = FALSE
  )
}

input_dir <- locate_input_dir()
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2L) {
  output_dir <- normalizePath(args[2], mustWork = FALSE)
} else {
  output_dir <- file.path(input_dir, "MAG_phylogeny_rectangular_R_output")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Input folder:  ", input_dir)
message("Output folder: ", output_dir)

# ---- 2. Publication palette ------------------------------------------------
COL <- list(
  ink = "#252525",
  tree = "#303030",
  mid_grey = "#777777",
  light_grey = "#D9D9D9",
  grid = "#E4E5E5",
  bar_background = "#F0F1F1",
  genome = "#7896A9",
  genome_edge = "#405C6A",
  cgenes = "#829B72",
  cgenes_edge = "#4F6645"
)

PHYLUM_COLORS <- c(
  Patescibacteria = "#6F8FB3",
  Proteobacteria = "#D49A5B",
  Bacteroidota = "#88A875",
  Actinobacteriota = "#A982A8",
  Cyanobacteria = "#C9AE58",
  Other = "#AFAFAF"
)

phylum_color <- function(phylum) {
  phylum <- trimws(as.character(phylum))
  ifelse(phylum %in% names(PHYLUM_COLORS), PHYLUM_COLORS[phylum], PHYLUM_COLORS["Other"])
}

pale_color <- function(color, white_fraction = 0.78) {
  rgb_value <- grDevices::col2rgb(color) / 255
  mixed <- rgb_value * (1 - white_fraction) + white_fraction
  grDevices::rgb(mixed[1], mixed[2], mixed[3])
}

# ---- 3. Read and validate MAG-result.csv -----------------------------------
result_file <- file.path(input_dir, "MAG-result.csv")
dat <- read.csv(result_file, check.names = FALSE, stringsAsFactors = FALSE)

required_columns <- c("MAG", "Phylum", "Carbon_Gene_Total")
missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns) > 0L) {
  stop(
    "MAG-result.csv is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

dat$MAG <- trimws(as.character(dat$MAG))
dat$Phylum <- trimws(as.character(dat$Phylum))
dat$Carbon_Gene_Total <- suppressWarnings(as.numeric(dat$Carbon_Gene_Total))

if (nrow(dat) != 33L || length(unique(dat$MAG)) != 33L) {
  stop("MAG-result.csv must contain exactly 33 unique target MAGs.", call. = FALSE)
}
if (any(!nzchar(dat$MAG)) || anyDuplicated(dat$MAG)) {
  stop("MAG IDs must be non-empty and unique.", call. = FALSE)
}
if (any(!is.finite(dat$Carbon_Gene_Total)) || any(dat$Carbon_Gene_Total < 0)) {
  stop("Carbon_Gene_Total must contain finite, nonnegative values.", call. = FALSE)
}
if (max(dat$Carbon_Gene_Total) > 200) {
  stop("Carbon_Gene_Total exceeds the fixed 0-200 plotting scale.", call. = FALSE)
}

# ---- 4. Calculate genome sizes from FASTA sequence lengths -----------------
canonical_mag_id <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  trimws(sub("\\([0-9]+\\)$", "", stem))
}

fasta_statistics <- function(path) {
  lines <- readLines(path, warn = FALSE)
  is_header <- startsWith(lines, ">")
  sequence_lines <- gsub("[[:space:]]", "", lines[!is_header])
  c(
    bases = sum(nchar(sequence_lines, type = "bytes")),
    contigs = sum(is_header)
  )
}

fasta_files <- list.files(
  input_dir,
  pattern = "\\.(fna|fa|fasta)$",
  full.names = TRUE,
  ignore.case = TRUE
)
if (length(fasta_files) == 0L) {
  stop("No .fna, .fa, or .fasta files were found.", call. = FALSE)
}

fasta_groups <- split(fasta_files, vapply(fasta_files, canonical_mag_id, character(1)))
genome_rows <- lapply(names(fasta_groups), function(mag) {
  paths <- fasta_groups[[mag]]
  statistics <- lapply(paths, fasta_statistics)
  signatures <- unique(vapply(statistics, function(x) paste(x, collapse = ":"), character(1)))
  if (length(signatures) > 1L) {
    stop(
      "Conflicting duplicate FASTA files for ", mag, ": ",
      paste(basename(paths), collapse = ", "),
      call. = FALSE
    )
  }

  stems <- tools::file_path_sans_ext(basename(paths))
  has_upload_suffix <- grepl("\\([0-9]+\\)$", stems)
  chosen <- order(has_upload_suffix, basename(paths))[1]
  values <- statistics[[chosen]]
  data.frame(
    MAG = mag,
    Genome_size_Mb = unname(values["bases"]) / 1e6,
    Contig_count = as.integer(unname(values["contigs"])),
    Genome_file = basename(paths[chosen]),
    stringsAsFactors = FALSE
  )
})
genomes <- do.call(rbind, genome_rows)

genome_match <- match(dat$MAG, genomes$MAG)
if (anyNA(genome_match)) {
  stop(
    "Missing FASTA file(s) for: ",
    paste(dat$MAG[is.na(genome_match)], collapse = ", "),
    call. = FALSE
  )
}

dat$Genome_size_Mb <- genomes$Genome_size_Mb[genome_match]
dat$Contig_count <- genomes$Contig_count[genome_match]
dat$Genome_file <- genomes$Genome_file[genome_match]

if (any(!is.finite(dat$Genome_size_Mb)) || any(dat$Genome_size_Mb <= 0)) {
  stop("All calculated genome sizes must be finite and positive.", call. = FALSE)
}
if (max(dat$Genome_size_Mb) > 6.5) {
  stop("Genome_size_Mb exceeds the fixed 0-6.5 Mb plotting scale.", call. = FALSE)
}

# ---- 5. Read and repair the local Newick tree ------------------------------
tree_file <- file.path(input_dir, "unroot.tree")

read_local_newick <- function(path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    stop("unroot.tree is missing or empty: ", path, call. = FALSE)
  }

  lines <- suppressWarnings(readLines(path, warn = FALSE, encoding = "UTF-8"))
  text <- paste(lines, collapse = "")
  text <- gsub(intToUtf8(65279), "", text, fixed = TRUE)
  text <- trimws(text)

  first_open <- regexpr("\\(", text, perl = TRUE)[1]
  close_positions <- gregexpr("\\)", text, perl = TRUE)[[1]]
  close_positions <- close_positions[close_positions > 0]
  if (!nzchar(text) || first_open < 1L || length(close_positions) == 0L) {
    stop("No readable Newick tree was found in unroot.tree.", call. = FALSE)
  }

  newick_body <- substr(text, first_open, max(close_positions))
  open_count <- lengths(regmatches(newick_body, gregexpr("\\(", newick_body, perl = TRUE)))
  close_count <- lengths(regmatches(newick_body, gregexpr("\\)", newick_body, perl = TRUE)))
  if (open_count != close_count) {
    stop("Unbalanced parentheses in unroot.tree.", call. = FALSE)
  }

  tryCatch(
    ape::read.tree(text = paste0(newick_body, ";")),
    error = function(e) {
      stop("ape could not parse unroot.tree: ", conditionMessage(e), call. = FALSE)
    }
  )
}

normalize_tree_tip <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub(intToUtf8(34), "", x, fixed = TRUE)
  x <- gsub("^'+|'+$", "", x)
  x <- basename(gsub("\\\\", "/", x))
  x <- sub("\\.(fna|fa|fasta)(\\.gz)?$", "", x, ignore.case = TRUE)
  x <- sub("\\([0-9]+\\)$", "", x)
  x <- sub("^X(?=S[0-9]+C[0-9]+$)", "", x, perl = TRUE)
  trimws(x)
}

recover_target_mag <- function(raw_label, normalized_label, targets) {
  if (normalized_label %in% targets) return(normalized_label)
  hits <- targets[vapply(targets, function(id) grepl(id, raw_label, fixed = TRUE), logical(1))]
  if (length(hits) == 1L) hits else normalized_label
}

phy <- read_local_newick(tree_file)
raw_tips <- phy$tip.label
normalized_tips <- normalize_tree_tip(raw_tips)
phy$tip.label <- mapply(
  recover_target_mag,
  raw_label = raw_tips,
  normalized_label = normalized_tips,
  MoreArgs = list(targets = dat$MAG),
  USE.NAMES = FALSE
)

if (anyDuplicated(phy$tip.label)) {
  stop("Tree-tip normalization produced duplicated MAG IDs.", call. = FALSE)
}

missing_from_tree <- setdiff(dat$MAG, phy$tip.label)
if (length(missing_from_tree) > 0L) {
  stop(
    "MAG/tree matching failed. Missing from unroot.tree: ",
    paste(missing_from_tree, collapse = ", "),
    call. = FALSE
  )
}

pruned_tips <- setdiff(phy$tip.label, dat$MAG)
if (length(pruned_tips) > 0L) {
  phy <- ape::drop.tip(phy, pruned_tips)
}
phy <- ape::ladderize(phy, right = FALSE)
phy <- ape::reorder.phylo(phy, order = "cladewise")

if (ape::Ntip(phy) != 33L) {
  stop("The pruned tree must contain exactly 33 target tips.", call. = FALSE)
}

message("Local tree path:  ", normalizePath(tree_file, mustWork = TRUE))
message("Local tree MD5:   ", unname(tools::md5sum(tree_file)))
message("Original tips:    ", length(raw_tips))
message("MAG/tree matched: 33/33")
message("Pruned tips:      ", if (length(pruned_tips)) paste(pruned_tips, collapse = ", ") else "none")

# ---- 6. Drawing helpers ----------------------------------------------------
draw_bar_panel <- function(values, y_positions, scale_max, display_max, ticks,
                           title, fill, border) {
  n <- length(values)
  par(mar = c(2.2, 0.55, 3.3, 0.45), xpd = FALSE, family = "sans")
  plot.new()
  plot.window(
    xlim = c(0, display_max), ylim = c(0.5, n + 0.5),
    xaxs = "i", yaxs = "i"
  )

  abline(v = ticks, col = COL$grid, lwd = 0.70)
  segments(0, 0.5, 0, n + 0.5, col = COL$ink, lwd = 0.80)
  rect(
    0, y_positions - 0.31, scale_max, y_positions + 0.31,
    col = COL$bar_background, border = "white", lwd = 0.25
  )
  rect(
    0, y_positions - 0.31, values, y_positions + 0.31,
    col = fill, border = border, lwd = 0.45
  )

  axis(
    3, at = ticks, labels = ticks, las = 1,
    cex.axis = 0.76, tck = -0.012, padj = -0.45,
    col = COL$light_grey, col.axis = COL$ink, col.ticks = COL$ink
  )
  mtext(title, side = 3, line = 1.65, cex = 1.04, font = 2, col = COL$ink)
  box(bty = "n")
}

draw_legend_panel <- function(observed_phyla) {
  par(mar = c(1.3, 0.5, 2.2, 0.2), xpd = NA, family = "sans")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  y <- 0.965

  heading <- function(label, gap = 0.055) {
    text(0.02, y, labels = label, adj = c(0, 1), cex = 1.04, font = 2, col = COL$ink)
    y <<- y - gap
  }
  item <- function(fill, label, border = "white") {
    rect(0.025, y - 0.018, 0.087, y + 0.013, col = fill, border = border, lwd = 0.55)
    text(0.105, y - 0.003, labels = label, adj = c(0, 0.5), cex = 0.83, col = COL$ink)
    y <<- y - 0.044
  }

  heading("Phylum")
  phylum_order <- c(
    "Patescibacteria", "Proteobacteria", "Bacteroidota",
    "Actinobacteriota", "Cyanobacteria"
  )
  for (phylum in phylum_order[phylum_order %in% observed_phyla]) {
    item(PHYLUM_COLORS[phylum], phylum)
  }

  y <- y - 0.030
  heading("Data tracks")
  item(COL$genome, "Genome size (Mb)", border = COL$genome_edge)
  item(COL$cgenes, "Carbon_Gene_Total (0-200)", border = COL$cgenes_edge)

  y <- y - 0.030
  heading("Node support", gap = 0.045)
  points(0.055, y, pch = 16, cex = 0.72, col = COL$ink)
  text(0.105, y, labels = ">=0.95", adj = c(0, 0.5), cex = 0.83, col = COL$ink)
  y <- y - 0.046
  points(0.055, y, pch = 21, cex = 0.75, bg = "white", col = COL$mid_grey, lwd = 0.8)
  text(0.105, y, labels = "0.70-0.95", adj = c(0, 0.5), cex = 0.83, col = COL$ink)
}

prepare_rectangular_layout <- function(tree) {
  n_tip <- ape::Ntip(tree)
  n_all <- n_tip + tree$Nnode
  parents <- tree$edge[, 1]
  children <- tree$edge[, 2]
  root_candidates <- setdiff(unique(parents), unique(children))
  if (length(root_candidates) != 1L) {
    stop("The pruned tree does not have one identifiable root.", call. = FALSE)
  }
  root <- root_candidates[1]
  child_map <- split(children, parents)

  get_children <- function(node) {
    result <- child_map[[as.character(node)]]
    if (is.null(result)) integer(0) else as.integer(result)
  }

  node_height <- rep(NA_integer_, n_all)
  calculate_height <- function(node) {
    if (!is.na(node_height[node])) return(node_height[node])
    kids <- get_children(node)
    if (length(kids) == 0L) {
      node_height[node] <<- 0L
    } else {
      node_height[node] <<- 1L + max(vapply(kids, calculate_height, integer(1)))
    }
    node_height[node]
  }
  max_height <- calculate_height(root)

  x <- max_height - node_height
  y <- rep(NA_real_, n_all)
  ordered_tip_nodes <- integer(0)
  tip_counter <- 0L

  assign_y <- function(node) {
    kids <- get_children(node)
    if (length(kids) == 0L) {
      tip_counter <<- tip_counter + 1L
      # The first depth-first tip is placed at the top of the figure.
      y[node] <<- n_tip - tip_counter + 1L
      ordered_tip_nodes <<- c(ordered_tip_nodes, node)
    } else {
      child_y <- vapply(kids, assign_y, numeric(1))
      y[node] <<- mean(child_y)
    }
    y[node]
  }
  assign_y(root)

  if (length(ordered_tip_nodes) != n_tip || any(ordered_tip_nodes > n_tip)) {
    stop("Failed to derive a complete rectangular tip order.", call. = FALSE)
  }

  support <- rep(NA_real_, n_all)
  if (!is.null(tree$node.label)) {
    support[n_tip + seq_along(tree$node.label)] <- suppressWarnings(
      as.numeric(tree$node.label)
    )
  }

  list(
    n_tip = n_tip,
    n_all = n_all,
    root = root,
    children = get_children,
    x = x,
    y = y,
    support = support,
    tip_nodes = ordered_tip_nodes,
    tip_names = tree$tip.label[ordered_tip_nodes]
  )
}

rect_layout <- prepare_rectangular_layout(phy)

draw_figure <- function() {
  n <- rect_layout$n_tip
  tree_max <- max(rect_layout$x[rect_layout$tip_nodes])
  band_left <- tree_max + 0.16
  band_right <- band_left + 3.35
  stripe_width <- 0.19
  label_x <- band_left + stripe_width + 0.14

  layout(matrix(1:4, nrow = 1), widths = c(7.7, 2.15, 2.35, 2.35))

  # Tree and phylum label band.
  par(mar = c(2.2, 0.45, 3.3, 0.05), xpd = NA, family = "sans")
  plot.new()
  plot.window(
    xlim = c(-0.20, band_right + 0.08),
    ylim = c(0.5, n + 0.5),
    xaxs = "i", yaxs = "i"
  )

  tip_x <- rect_layout$x[rect_layout$tip_nodes]
  tip_y <- rect_layout$y[rect_layout$tip_nodes]
  tip_names <- rect_layout$tip_names
  dat_by_mag <- dat[match(tip_names, dat$MAG), , drop = FALSE]

  # Phylum-coloured label rows are drawn first, behind all tree elements.
  for (i in seq_len(n)) {
    phylum <- dat_by_mag$Phylum[i]
    color <- phylum_color(phylum)
    rect(
      band_left, tip_y[i] - 0.39, band_right, tip_y[i] + 0.39,
      col = pale_color(color), border = "white", lwd = 0.25
    )
    rect(
      band_left, tip_y[i] - 0.39, band_left + stripe_width, tip_y[i] + 0.39,
      col = color, border = "white", lwd = 0.25
    )
  }

  # Rectangular cladogram: vertical connector per internal node, followed by
  # one horizontal parent-child segment for every edge.
  for (node in seq.int(n + 1L, rect_layout$n_all)) {
    kids <- rect_layout$children(node)
    if (length(kids) == 0L) next
    segments(
      rect_layout$x[node], min(rect_layout$y[kids]),
      rect_layout$x[node], max(rect_layout$y[kids]),
      col = COL$tree, lwd = 0.82, lend = "butt"
    )
    for (child in kids) {
      segments(
        rect_layout$x[node], rect_layout$y[child],
        rect_layout$x[child], rect_layout$y[child],
        col = COL$tree, lwd = 0.82, lend = "butt"
      )
    }
  }

  for (i in seq_len(n)) {
    segments(tip_x[i], tip_y[i], band_left, tip_y[i], col = COL$tree, lwd = 0.72)
    text(
      label_x, tip_y[i], labels = tip_names[i],
      adj = c(0, 0.5), cex = 0.71, font = 2, col = COL$ink
    )
  }

  # Node-support markers from Newick node labels.
  for (node in seq.int(n + 1L, rect_layout$n_all)) {
    value <- rect_layout$support[node]
    if (!is.finite(value) || value < 0.70) next
    if (value >= 0.95) {
      points(rect_layout$x[node], rect_layout$y[node], pch = 16, cex = 0.42, col = COL$ink)
    } else {
      points(
        rect_layout$x[node], rect_layout$y[node],
        pch = 21, cex = 0.45, bg = "white", col = COL$mid_grey, lwd = 0.75
      )
    }
  }

  mtext("MAG phylogeny", side = 3, line = 1.65, adj = 0, cex = 1.08, font = 2, col = COL$ink)
  text(
    (band_left + band_right) / 2, n + 0.90, labels = "MAG",
    adj = c(0.5, 0), cex = 0.84, font = 2, col = COL$ink
  )
  mtext("(A)", side = 3, line = 2.55, adj = 0, cex = 2.1, font = 2, col = COL$ink)

  # Use the exact tree-tip y coordinates in both bar panels.
  draw_bar_panel(
    values = dat_by_mag$Genome_size_Mb,
    y_positions = tip_y,
    scale_max = 6.5,
    display_max = 7.25,
    ticks = c(0, 3, 6),
    title = "Genome size (Mb)",
    fill = COL$genome,
    border = COL$genome_edge
  )

  draw_bar_panel(
    values = dat_by_mag$Carbon_Gene_Total,
    y_positions = tip_y,
    scale_max = 200,
    display_max = 225,
    ticks = c(0, 100, 200),
    title = "Carbon_Gene_Total (0-200)",
    fill = COL$cgenes,
    border = COL$cgenes_edge
  )

  draw_legend_panel(unique(dat_by_mag$Phylum))

  # Return the exact displayed top-to-bottom order for the output data table.
  tip_names[order(tip_y, decreasing = TRUE)]
}

# ---- 7. Export PDF, SVG, 600-dpi PNG, and preview ---------------------------
pdf_file <- file.path(output_dir, "MAG_phylogeny_rectangular_bars_R_SCI.pdf")
svg_file <- file.path(output_dir, "MAG_phylogeny_rectangular_bars_R_SCI.svg")
png_file <- file.path(output_dir, "MAG_phylogeny_rectangular_bars_R_SCI.png")
preview_file <- file.path(output_dir, "MAG_phylogeny_rectangular_bars_R_preview.png")

grDevices::pdf(
  pdf_file, width = 15.6, height = 10.8,
  family = "Helvetica", useDingbats = FALSE, bg = "white"
)
tree_order <- draw_figure()
grDevices::dev.off()

grDevices::svg(
  svg_file, width = 15.6, height = 10.8,
  pointsize = 10, family = "sans", bg = "white"
)
draw_figure()
grDevices::dev.off()

grDevices::png(
  png_file, width = 15.6, height = 10.8,
  units = "in", res = 600, pointsize = 10, bg = "white"
)
draw_figure()
grDevices::dev.off()

grDevices::png(
  preview_file, width = 15.6, height = 10.8,
  units = "in", res = 180, pointsize = 10, bg = "white"
)
draw_figure()
grDevices::dev.off()

# ---- 8. Save exact plotted data and validation report ----------------------
plot_data <- dat[match(tree_order, dat$MAG), c(
  "MAG", "Phylum", "Genome_size_Mb", "Contig_count",
  "Genome_file", "Carbon_Gene_Total"
), drop = FALSE]
plot_data$Tree_order_top_to_bottom <- seq_len(nrow(plot_data))
plot_data <- plot_data[, c(
  "Tree_order_top_to_bottom", "MAG", "Phylum", "Genome_size_Mb",
  "Contig_count", "Genome_file", "Carbon_Gene_Total"
)]

write.csv(
  plot_data,
  file.path(output_dir, "MAG_phylogeny_rectangular_bars_R_plot_data.csv"),
  row.names = FALSE,
  quote = TRUE
)

validation <- paste0(
  "Validation status: ready to plot.\n",
  "Script version: ", SCRIPT_VERSION, "\n",
  "Tree file: ", normalizePath(tree_file, mustWork = TRUE), "\n",
  "Tree MD5: ", unname(tools::md5sum(tree_file)), "\n",
  "Target MAGs: 33/33 matched.\n",
  "Pruned tree tips: ", if (length(pruned_tips)) paste(pruned_tips, collapse = ", ") else "none", "\n",
  "FASTA coverage: 33/33.\n",
  sprintf("Genome-size range: %.6f-%.6f Mb.\n", min(dat$Genome_size_Mb), max(dat$Genome_size_Mb)),
  sprintf(
    "Carbon_Gene_Total range: %d-%d; source: MAG-result.csv.\n",
    min(dat$Carbon_Gene_Total), max(dat$Carbon_Gene_Total)
  ),
  "C gene num.csv is not used because it represents a different annotation scope.\n",
  "Per-MAG numeric labels are intentionally omitted; quantitative axes start at zero.\n"
)
writeLines(
  validation,
  file.path(output_dir, "MAG_phylogeny_rectangular_bars_R_validation.txt"),
  useBytes = TRUE
)

message(validation)
message("Finished. Files were written to: ", output_dir)
