# utility functions shared across phaseTime

#' Chromosome sizes of the GRCh38 (hg38) autosomes
#'
#' Named numeric vector with the lengths (bp) of chr1-chr22 in GRCh38, as
#' given by `BSgenome.Hsapiens.UCSC.hg38`. Used to lay out genome-wide plots
#' and simulations. Pass your own named vector to the `chr_sizes` argument of
#' the plotting functions to use a different reference.
#'
#' @format Named numeric vector of length 22.
#' @export
hg38_chr_sizes <- c(
  chr1 = 248956422, chr2 = 242193529, chr3 = 198295559, chr4 = 190214555,
  chr5 = 181538259, chr6 = 170805979, chr7 = 159345973, chr8 = 145138636,
  chr9 = 138394717, chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
  chr13 = 114364328, chr14 = 107043718, chr15 = 101991189, chr16 = 90338345,
  chr17 = 83257441, chr18 = 80373285, chr19 = 58617616, chr20 = 64444167,
  chr21 = 46709983, chr22 = 50818468
)

# colours used throughout the thesis figures
.cols <- list(
  hapA_cn   = "darkred",
  hapB_cn   = "#313695",
  hapA_mut  = "#d73027",
  hapB_mut  = "#4575b4",
  unphased  = "#bababa",
  grey      = "#878787",
  early     = "#e94e1b",
  late      = "#662483",
  state_one = "#d8b365",
  state_two = "#5ab4ac"
)

#' Make sure chromosome names carry the "chr" prefix
#' @param chr character or numeric vector of chromosome names.
#' @return character vector with "chr" prefix.
#' @export
#' @examples
#' add_chr_prefix(c(1, "2", "chr3"))
add_chr_prefix <- function(chr) {
  chr <- as.character(chr)
  ifelse(grepl("^chr", chr), chr, paste0("chr", chr))
}

#' Build segment identifiers
#'
#' Segment identifiers have the form `chr_start_end`, as used throughout
#' phaseTime output tables.
#' @param chr,seg_start,seg_end vectors of chromosome, segment start and end.
#' @return character vector of segment identifiers.
#' @export
#' @examples
#' seg_id("chr1", 1, 1000)
seg_id <- function(chr, seg_start, seg_end) {
  paste0(chr, "_", format(seg_start, scientific = FALSE, trim = TRUE), "_",
         format(seg_end, scientific = FALSE, trim = TRUE))
}

#' Split segment identifiers into chromosome, start and end
#' @param segID character vector of identifiers created by [seg_id()].
#' @return data.frame with columns `chr`, `seg_start`, `seg_end`.
#' @export
#' @examples
#' split_seg_id("chr14_47734531_83979663")
split_seg_id <- function(segID) {
  parts <- do.call(rbind, strsplit(segID, "_", fixed = TRUE))
  data.frame(chr = parts[, 1],
             seg_start = as.numeric(parts[, 2]),
             seg_end = as.numeric(parts[, 3]),
             stringsAsFactors = FALSE)
}

#' Convert chromosome positions to continuous genome coordinates
#'
#' Adds columns `adjust_pos`, `adjust_seg_start` and `adjust_seg_end` (where
#' the corresponding columns `pos`, `seg_start`, `seg_end` exist) by adding the
#' cumulative length of all preceding chromosomes. Rows on chromosomes not in
#' `chr_sizes` are dropped.
#'
#' @param data data.frame with a `chr` column and any of `pos`, `seg_start`,
#'   `seg_end`.
#' @param chr_sizes named vector of chromosome sizes (default GRCh38
#'   autosomes).
#' @return `data` restricted to chromosomes in `chr_sizes` with the adjusted
#'   coordinate columns added.
#' @export
#' @examples
#' genome_continuous(data.frame(chr = c("chr1", "chr2"), pos = c(10, 10)))
genome_continuous <- function(data, chr_sizes = hg38_chr_sizes) {
  data <- as.data.frame(data)
  data$chr <- add_chr_prefix(data$chr)
  data <- data[data$chr %in% names(chr_sizes), , drop = FALSE]
  offset <- c(0, cumsum(as.numeric(chr_sizes)))[seq_along(chr_sizes)]
  names(offset) <- names(chr_sizes)
  off <- offset[data$chr]
  if ("pos" %in% names(data)) data$adjust_pos <- data$pos + off
  if ("seg_start" %in% names(data)) data$adjust_seg_start <- data$seg_start + off
  if ("seg_end" %in% names(data)) data$adjust_seg_end <- data$seg_end + off
  data
}

# chromosome boundary lines and label positions for genome-wide plots
.genome_axis <- function(chr_sizes = hg38_chr_sizes) {
  lines <- cumsum(as.numeric(chr_sizes))
  list(lines = lines,
       midpoints = lines - as.numeric(chr_sizes) / 2,
       labels = names(chr_sizes))
}

# map point positions to segments (per chromosome, assumes non-overlapping
# segments). Returns the row index in `segments` or NA.
.match_segments <- function(chr, pos, seg_chr, seg_start, seg_end) {
  idx <- rep(NA_integer_, length(pos))
  for (ch in unique(chr)) {
    w <- which(chr == ch)
    s <- which(seg_chr == ch)
    if (length(s) == 0) next
    s <- s[order(seg_start[s])]
    k <- findInterval(pos[w], seg_start[s])
    ok <- k > 0
    hit <- rep(NA_integer_, length(w))
    hit[ok] <- s[k[ok]]
    ok2 <- ok & pos[w] <= seg_end[hit]
    hit[!ok2] <- NA_integer_
    idx[w] <- hit
  }
  idx
}

# cap a numeric vector for plotting (the thesis plots cap copy number at 10)
.cap <- function(x, max_value = 10) pmin(x, max_value)

# check that required columns exist
.require_cols <- function(df, cols, what = "data") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(what, " is missing required column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}

# is x (numerically) an integer
.is_whole <- function(x, tol = 1e-8) abs(x - round(x)) < tol

#' Colour palette for gain order
#'
#' Returns colours from the ColorBrewer "Paired" palette (hard coded so that
#' RColorBrewer is not needed), recycled with decreasing opacity if more than
#' 12 colours are requested. Used to colour consecutive copy number gains.
#' @param n number of colours.
#' @return character vector of hex colours.
#' @export
#' @examples
#' gain_colours(3)
gain_colours <- function(n) {
  paired <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
              "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
  rep_len(paired, n)
}

# the y axis used for all gain timing plots
.timing_y_scale <- function() {
  ggplot2::scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                              labels = c("0.00\n(early)", "0.25", "0.50", "0.75", "1.00\n(late)"),
                              limits = c(0, 1))
}

.genome_theme <- function(show_x = TRUE) {
  th <- ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "white"))
  if (show_x) {
    th + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5))
  } else {
    th + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                        axis.ticks.x = ggplot2::element_blank())
  }
}

.genome_x_scale <- function(chr_sizes) {
  ax <- .genome_axis(chr_sizes)
  ggplot2::scale_x_continuous(breaks = ax$midpoints, labels = ax$labels,
                              expand = c(0, 0),
                              limits = c(0, max(ax$lines)))
}
