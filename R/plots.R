# plotting functions reproducing the phaseTime figures of thesis chapter 6

# copy number segments for plotting from a cn table or from the mutation table
.segments_for_plot <- function(cn, mutations, chr_sizes, max_cn) {
  src <- if (!is.null(cn)) cn else mutations
  seg <- unique(as.data.frame(src)[, c("chr", "seg_start", "seg_end", "cnA", "cnB")])
  seg$chr <- add_chr_prefix(seg$chr)
  seg$allelic_imbalance <- factor(abs(seg$cnA - seg$cnB) > 0.1, levels = c(TRUE, FALSE))
  # offset equal copy numbers so that both haplotypes remain visible
  eq <- abs(seg$cnA - seg$cnB) < 1e-8
  seg$yA <- .cap(seg$cnA, max_cn) + ifelse(eq, 0.05, 0)
  seg$yB <- .cap(seg$cnB, max_cn) - ifelse(eq, 0.05, 0)
  genome_continuous(seg, chr_sizes)
}

.filter_chr <- function(d, chr) {
  if (is.null(chr)) return(d)
  d[d$chr %in% add_chr_prefix(chr), , drop = FALSE]
}

.base_cn_plot <- function(seg, xs, xe, haplotypes = c("A", "B"), alpha_ai = TRUE) {
  p <- ggplot2::ggplot()
  if ("A" %in% haplotypes) {
    p <- p + ggplot2::geom_segment(
      data = seg, ggplot2::aes(x = .data[[xs]], xend = .data[[xe]], y = .data$yA, yend = .data$yA,
                               colour = "haplotype A", alpha = .data$allelic_imbalance),
      linewidth = 2)
  }
  if ("B" %in% haplotypes) {
    p <- p + ggplot2::geom_segment(
      data = seg, ggplot2::aes(x = .data[[xs]], xend = .data[[xe]], y = .data$yB, yend = .data$yB,
                               colour = "haplotype B", alpha = .data$allelic_imbalance),
      linewidth = 2)
  }
  p + ggplot2::scale_colour_manual(name = "copy number\nhaplotype",
                                   values = c("haplotype A" = .cols$hapA_cn, "haplotype B" = .cols$hapB_cn)) +
    ggplot2::scale_alpha_manual(name = "allelic imbalance", values = c("TRUE" = 0.8, "FALSE" = 0.3),
                                drop = FALSE)
}

.x_layer <- function(p, chr, chr_sizes, show_x = TRUE) {
  if (is.null(chr)) {
    ax <- .genome_axis(chr_sizes)
    p + ggplot2::geom_vline(xintercept = ax$lines, linetype = 2, colour = "#bababa") +
      .genome_x_scale(chr_sizes) + ggplot2::xlab("") + .genome_theme(show_x)
  } else {
    p + ggplot2::xlab(paste(chr, collapse = ", ")) + .genome_theme(TRUE) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0))
  }
}

#' Plot SNP B-allele frequencies coloured by haplotype
#'
#' Reproduces thesis Figure 6-5: BAF of heterozygous SNPs across the genome,
#' coloured by the haplotype obtained from this tumour region alone (top) and
#' by the multi-region reference haplotype (bottom).
#'
#' @param snps output of [phase_snps()].
#' @param region tumour region to plot.
#' @param which `"both"`, `"single"` or `"multi"`.
#' @param chr optional chromosome(s) to plot; default whole genome.
#' @param chr_sizes chromosome sizes.
#' @param point_size point size.
#' @return a ggplot / patchwork object.
#' @export
plot_snp_haplotypes <- function(snps, region = NULL, which = c("both", "single", "multi"),
                                chr = NULL, chr_sizes = hg38_chr_sizes, point_size = 0.4) {
  which <- match.arg(which)
  if (is.null(region)) region <- snps$region[1]
  d <- snps[snps$region == region, , drop = FALSE]
  d <- .filter_chr(genome_continuous(d, chr_sizes), chr)
  xcol <- if (is.null(chr)) "adjust_pos" else "pos"
  mk <- function(col, title) {
    d$h <- factor(ifelse(is.na(d[[col]]), "unphased", d[[col]]), levels = c("A", "B", "unphased"))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[xcol]], y = .data$baf, colour = .data$h)) +
      ggplot2::geom_point(size = point_size, alpha = 0.6, stroke = 0) +
      ggplot2::scale_colour_manual(name = title, values = c(A = .cols$hapA_mut, B = .cols$hapB_mut,
                                                            unphased = .cols$unphased),
                                   labels = c("haplotype A", "haplotype B", "not phasable"), drop = FALSE) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2))) +
      ggplot2::ylab("B-allele frequency (BAF)") + ggplot2::ylim(0, 1)
    .x_layer(p, chr, chr_sizes)
  }
  p1 <- mk("single_region_haplotype", "SNP single region\nhaplotype")
  p2 <- mk("haplotype", "SNP multi region\nhaplotype")
  switch(which,
         single = p1 + ggplot2::ggtitle(region),
         multi = p2 + ggplot2::ggtitle(region),
         both = p1 / p2 + patchwork::plot_annotation(title = region))
}

#' Plot phased mutations and haplotype specific copy number
#'
#' Reproduces thesis Figure 6-8 (whole genome) or its per-chromosome version.
#' Mutations are plotted at their multiplicity (mutation copy number, capped
#' at `max_cn`) together with the haplotype specific copy number segments.
#'
#' @param mutations phased mutations with `chr`, `pos`, `mut_cpn`, `phase`,
#'   `cnA`, `cnB`, `seg_start`, `seg_end`.
#' @param cn optional haplotype copy number table ([haplotype_copy_number()])
#'   used to draw segments without mutations.
#' @param type `"copy_number"` (Figure 6-8A), `"all"` (6-8B, all mutations in
#'   grey), `"haplotype"` (6-8C, coloured by haplotype), `"A"` (6-8D) or
#'   `"B"` (6-8E).
#' @param chr optional chromosome(s); default whole genome.
#' @param max_cn copy numbers and multiplicities are capped at this value.
#' @param chr_sizes chromosome sizes.
#' @param point_size point size.
#' @return a ggplot object.
#' @export
plot_phased_mutations <- function(mutations, cn = NULL,
                                  type = c("haplotype", "all", "copy_number", "A", "B"),
                                  chr = NULL, max_cn = 10, chr_sizes = hg38_chr_sizes,
                                  point_size = 0.5) {
  type <- match.arg(type)
  seg <- .filter_chr(.segments_for_plot(cn, mutations, chr_sizes, max_cn), chr)
  m <- .filter_chr(genome_continuous(mutations, chr_sizes), chr)
  m$y <- .cap(m$mut_cpn, max_cn)
  m$ph <- factor(ifelse(is.na(m$phase), "unphased", m$phase), levels = c("A", "B", "unphased"))
  genome <- is.null(chr)
  xs <- if (genome) "adjust_seg_start" else "seg_start"
  xe <- if (genome) "adjust_seg_end" else "seg_end"
  xp <- if (genome) "adjust_pos" else "pos"
  haps <- switch(type, A = "A", B = "B", c("A", "B"))
  if (type == "A") m <- m[m$ph == "A", ]
  if (type == "B") m <- m[m$ph == "B", ]

  p <- .base_cn_plot(seg, xs, xe, haps)
  if (type != "copy_number") {
    vals <- c(A = .cols$hapA_mut, B = .cols$hapB_mut, unphased = .cols$unphased)
    if (type == "all") vals[] <- .cols$grey
    p <- p + ggplot2::geom_point(data = m, ggplot2::aes(x = .data[[xp]], y = .data$y, fill = .data$ph),
                                 shape = 21, stroke = 0, size = point_size, alpha = 0.6) +
      ggplot2::scale_fill_manual(name = "mutation haplotype", values = vals,
                                 labels = c(A = "haplotype A", B = "haplotype B", unphased = "unphased"),
                                 drop = FALSE) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 2.5)))
  }
  p <- p + ggplot2::scale_y_continuous(limits = c(-0.1, max_cn + 0.1), breaks = seq(0, max_cn, 2)) +
    ggplot2::ylab(if (type == "copy_number") "copy number" else "mutation copy number")
  .x_layer(p, chr, chr_sizes)
}

#' Plot early and late mutations per haplotype
#'
#' Reproduces thesis Figure 6-9: phased mutations on each haplotype coloured
#' by whether they occurred before (early) or after (late) the copy number
#' gain of the segment, with the haplotype specific copy number.
#'
#' @param mutations output of [classify_mutation_timing()].
#' @param chr optional chromosome(s); default whole genome.
#' @param haplotypes haplotypes to show.
#' @param colour_by `"timing"` (Figure 6-9D) or `"haplotype"` (Figure 6-9A-C).
#' @param max_cn cap for copy number and multiplicity.
#' @param chr_sizes chromosome sizes.
#' @param point_size point size.
#' @return a ggplot object faceted by haplotype.
#' @export
plot_mutation_timing <- function(mutations, chr = NULL, haplotypes = c("A", "B"),
                                 colour_by = c("timing", "haplotype"), max_cn = 10,
                                 chr_sizes = hg38_chr_sizes, point_size = 0.6) {
  colour_by <- match.arg(colour_by)
  if (!"timing" %in% names(mutations)) mutations <- classify_mutation_timing(mutations)
  m <- .filter_chr(genome_continuous(mutations[!is.na(mutations$phase), ], chr_sizes), chr)
  m <- m[m$phase %in% haplotypes, ]
  m$y <- .cap(m$mut_cpn, max_cn)
  m$haplotype <- paste("haplotype", m$phase)
  genome <- is.null(chr)
  xp <- if (genome) "adjust_pos" else "pos"
  xs <- if (genome) "adjust_seg_start" else "seg_start"
  xe <- if (genome) "adjust_seg_end" else "seg_end"
  seg <- unique(m[, c("chr", "seg_start", "seg_end", xs, xe, "cnA", "cnB")])
  seg <- rbind(data.frame(seg, haplotype = "haplotype A", y = .cap(seg$cnA, max_cn)),
               data.frame(seg, haplotype = "haplotype B", y = .cap(seg$cnB, max_cn)))
  seg <- seg[seg$haplotype %in% paste("haplotype", haplotypes), ]
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg, ggplot2::aes(x = .data[[xs]], xend = .data[[xe]], y = .data$y,
                                                   yend = .data$y, colour = .data$haplotype),
                          linewidth = 2, alpha = 0.7) +
    ggplot2::scale_colour_manual(name = "copy number", values = c("haplotype A" = .cols$hapA_cn,
                                                                  "haplotype B" = .cols$hapB_cn))
  if (colour_by == "timing") {
    p <- p + ggplot2::geom_point(data = m, ggplot2::aes(x = .data[[xp]], y = .data$y, fill = .data$timing),
                                 shape = 21, stroke = 0, size = point_size, alpha = 0.8) +
      ggplot2::scale_fill_manual(name = "mutation timing",
                                 values = c(early = .cols$early, late = .cols$late, not_possible = .cols$unphased),
                                 labels = c("early", "late", "not possible"), drop = FALSE)
  } else {
    p <- p + ggplot2::geom_point(data = m, ggplot2::aes(x = .data[[xp]], y = .data$y, fill = .data$phase),
                                 shape = 21, stroke = 0, size = point_size, alpha = 0.8) +
      ggplot2::scale_fill_manual(name = "mutation haplotype",
                                 values = c(A = .cols$hapA_mut, B = .cols$hapB_mut))
  }
  p <- p + ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 2.5))) +
    ggplot2::facet_wrap(~haplotype, ncol = 1, strip.position = "right") +
    ggplot2::ylab("mutation copy number")
  .x_layer(p, chr, chr_sizes)
}

#' Plot mutations coloured by (ancestral) copy number state
#'
#' Reproduces thesis Figure 6-10C: mutations of one haplotype coloured with
#' alternating colours by the integer copy number state they were assigned to
#' by [assign_cn_states()].
#'
#' @param mutations output of [assign_cn_states()].
#' @param haplotype haplotype to plot.
#' @param chr optional chromosome; `segID` optional segment to zoom into.
#' @param segID optional segment identifier.
#' @param max_cn cap for multiplicity.
#' @param chr_sizes chromosome sizes.
#' @param point_size point size.
#' @return a ggplot object.
#' @export
plot_cn_states <- function(mutations, haplotype = "A", chr = NULL, segID = NULL, max_cn = 10,
                           chr_sizes = hg38_chr_sizes, point_size = 1) {
  if (!"cn_state" %in% names(mutations)) mutations <- assign_cn_states(mutations)
  m <- mutations[!is.na(mutations$phase) & mutations$phase == haplotype, , drop = FALSE]
  if (!is.null(segID)) {
    m <- m[m$segID %in% segID, , drop = FALSE]
    chr <- unique(m$chr)
  }
  m <- .filter_chr(genome_continuous(m, chr_sizes), chr)
  m$y <- .cap(m$mut_cpn, max_cn)
  m$state_col <- ifelse(is.na(m$cn_state), "none", ifelse(m$cn_state %% 2 == 1, "odd", "even"))
  m$hcn <- .cap(if (haplotype == "A") m$cnA else m$cnB, max_cn)
  genome <- is.null(chr)
  xp <- if (genome) "adjust_pos" else "pos"
  xs <- if (genome) "adjust_seg_start" else "seg_start"
  xe <- if (genome) "adjust_seg_end" else "seg_end"
  seg <- unique(m[, c(xs, xe, "hcn")])
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg, ggplot2::aes(x = .data[[xs]], xend = .data[[xe]], y = .data$hcn,
                                                   yend = .data$hcn),
                          colour = if (haplotype == "A") .cols$hapA_cn else .cols$hapB_cn,
                          linewidth = 2, alpha = 0.7) +
    ggplot2::geom_hline(yintercept = seq(0.5, max_cn, 1), linetype = 2, colour = "#d9d9d9") +
    ggplot2::geom_point(data = m, ggplot2::aes(x = .data[[xp]], y = .data$y, fill = .data$state_col),
                        shape = 21, stroke = 0, size = point_size, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c(odd = .cols$state_one, even = .cols$state_two,
                                          none = .cols$unphased), guide = "none") +
    ggplot2::ylab("mutation copy number") +
    ggplot2::ggtitle(paste("haplotype", haplotype))
  .x_layer(p, chr, chr_sizes)
}

#' Plot genome-wide haplotype specific gain timing
#'
#' Reproduces thesis Figure 6-15. From top to bottom: haplotype specific copy
#' number; haplotype A copy number and mutations (yellow: acquired before the
#' first gain, blue: after the first gain, grey: not timed); haplotype A gain
#' timing estimates with 95% confidence intervals; and the same two panels for
#' haplotype B.
#'
#' @param timing output of [time_gains()].
#' @param mutations phased mutations used for timing (with `cn_state`, as
#'   returned by [assign_cn_states()]; computed if missing).
#' @param cn optional haplotype copy number table for the top panel.
#' @param gain_index which gain to show in the timing panels (default the
#'   first gain).
#' @param clusters optional cluster(s) whose timing estimates are shown;
#'   default all clusters in `timing`.
#' @param title plot title.
#' @param max_cn cap for copy number / multiplicity.
#' @param chr_sizes chromosome sizes.
#' @param point_size point size.
#' @return a patchwork object.
#' @export
plot_gain_timing <- function(timing, mutations, cn = NULL, gain_index = 1, clusters = NULL,
                             title = NULL, max_cn = 10, chr_sizes = hg38_chr_sizes,
                             point_size = 0.5) {
  ax <- .genome_axis(chr_sizes)
  if (!"cn_state" %in% names(mutations)) mutations <- assign_cn_states(mutations)
  t <- timing[!is.na(timing$gain_timing) & timing$gain_index %in% gain_index, , drop = FALSE]
  if (!is.null(clusters)) t <- t[t$cluster %in% clusters, , drop = FALSE]
  t <- genome_continuous(t, chr_sizes)
  seg <- .segments_for_plot(cn, mutations, chr_sizes, max_cn)

  m <- genome_continuous(mutations[!is.na(mutations$phase), ], chr_sizes)
  m$y <- .cap(m$mut_cpn, max_cn)
  timed_keys <- paste(t$segID, t$haplotype)
  top <- stats::ave(ifelse(is.na(m$cn_state), -Inf, m$cn_state), paste(m$segID, m$phase), FUN = max)
  m$col <- "not timed"
  is_timed <- paste(m$segID, m$phase) %in% timed_keys & !is.na(m$cn_state)
  m$col[is_timed & m$cn_state == top] <- "before first gain"
  m$col[is_timed & m$cn_state < top] <- "after first gain"
  m$col <- factor(m$col, levels = c("before first gain", "after first gain", "not timed"))

  vlines <- ggplot2::geom_vline(xintercept = ax$lines, linetype = 2, colour = "#bababa")
  p_cn <- .base_cn_plot(seg, "adjust_seg_start", "adjust_seg_end") + vlines +
    ggplot2::scale_y_continuous(breaks = 0:max_cn, limits = c(-0.2, max_cn + 0.2)) +
    ggplot2::ylab("copy number") + .genome_x_scale(chr_sizes) + .genome_theme(FALSE) +
    ggplot2::guides(alpha = "none") + ggplot2::xlab("")

  mut_panel <- function(h) {
    s <- seg
    s$y <- if (h == "A") s$yA else s$yB
    ggplot2::ggplot() + vlines +
      ggplot2::geom_segment(data = s, ggplot2::aes(x = .data$adjust_seg_start, xend = .data$adjust_seg_end,
                                                   y = .data$y, yend = .data$y),
                            colour = if (h == "A") .cols$hapA_cn else .cols$hapB_cn, linewidth = 1, alpha = 0.7) +
      ggplot2::geom_point(data = m[m$phase == h, ], ggplot2::aes(x = .data$adjust_pos, y = .data$y,
                                                                 fill = .data$col),
                          shape = 21, stroke = 0, size = point_size, alpha = 0.7) +
      ggplot2::scale_fill_manual(name = "mutations", values = c("before first gain" = .cols$state_one,
                                                                "after first gain" = .cols$state_two,
                                                                "not timed" = .cols$unphased),
                                 drop = FALSE) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 2.5))) +
      ggplot2::facet_wrap(~"x", labeller = ggplot2::as_labeller(c(x = paste("haplotype", h))),
                          strip.position = "right") +
      ggplot2::scale_y_continuous(limits = c(-0.2, max_cn + 0.2)) +
      ggplot2::ylab("mutation\ncopy number") + .genome_x_scale(chr_sizes) + .genome_theme(FALSE) +
      ggplot2::xlab("")
  }
  time_panel <- function(h, show_x) {
    th <- t[t$haplotype == h, , drop = FALSE]
    col <- if (h == "A") .cols$hapA_cn else .cols$hapB_cn
    ggplot2::ggplot() + vlines +
      ggplot2::geom_rect(data = th, ggplot2::aes(xmin = .data$adjust_seg_start, xmax = .data$adjust_seg_end,
                                                 ymin = .data$loCI, ymax = .data$hiCI),
                         fill = col, alpha = 0.15) +
      ggplot2::geom_segment(data = th, ggplot2::aes(x = .data$adjust_seg_start, xend = .data$adjust_seg_end,
                                                    y = .data$gain_timing, yend = .data$gain_timing),
                            colour = col, linewidth = 2) +
      ggplot2::facet_wrap(~"x", labeller = ggplot2::as_labeller(c(x = paste("haplotype", h))),
                          strip.position = "right") +
      .timing_y_scale() + ggplot2::ylab("gain timing") + .genome_x_scale(chr_sizes) +
      .genome_theme(show_x) + ggplot2::xlab("")
  }
  p <- p_cn / mut_panel("A") / time_panel("A", FALSE) / mut_panel("B") / time_panel("B", TRUE) +
    patchwork::plot_layout(guides = "collect")
  if (!is.null(title)) p <- p + patchwork::plot_annotation(title = title)
  p
}

#' Heatmap of gain timing across tumour regions
#'
#' Reproduces thesis Figure 6-16A: timing of the first clonal gain of every
#' timed segment, one row per tumour region, separately for each haplotype.
#'
#' @param timing combined output of [time_gains()] for several samples.
#' @param gain_index gain to show (default 1, the first gain).
#' @param chr_sizes chromosome sizes.
#' @return a ggplot object.
#' @export
plot_timing_heatmap <- function(timing, gain_index = 1, chr_sizes = hg38_chr_sizes) {
  t <- timing[!is.na(timing$gain_timing) & timing$gain_index %in% gain_index, , drop = FALSE]
  t <- genome_continuous(t, chr_sizes)
  samples <- unique(timing$sample)
  t$row <- match(t$sample, samples)
  ax <- .genome_axis(chr_sizes)
  ggplot2::ggplot(t) +
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$adjust_seg_start, xmax = .data$adjust_seg_end,
                                    ymin = .data$row - 0.5, ymax = .data$row + 0.5,
                                    fill = .data$gain_timing)) +
    ggplot2::geom_vline(xintercept = ax$lines, linetype = 2, colour = "black", linewidth = 0.3) +
    ggplot2::scale_fill_gradientn(name = "gain timing",
                                  colours = c("#8e0152", "#c51b7d", "#de77ae", "#f1b6da", "#fde0ef",
                                              "#e6f5d0", "#b8e186", "#7fbc41", "#4d9221", "#276419"),
                                  limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                                  labels = c("0 (early)", "0.25", "0.5", "0.75", "1 (late)")) +
    ggplot2::scale_x_continuous(breaks = ax$midpoints, labels = ax$labels, expand = c(0, 0),
                                limits = c(0, max(ax$lines))) +
    ggplot2::scale_y_continuous(breaks = seq_along(samples), labels = samples, expand = c(0, 0),
                                limits = c(0.5, length(samples) + 0.5)) +
    ggplot2::facet_wrap(~haplotype, ncol = 1, strip.position = "right",
                        labeller = ggplot2::as_labeller(c(A = "haplotype A", B = "haplotype B"))) +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5),
                   panel.grid = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill = "#d9d9d9"),
                   strip.background = ggplot2::element_rect(fill = "white"))
}

#' Plot the concordance of gain timing between tumour regions
#'
#' Reproduces thesis Figure 6-16B and C: histograms of the absolute
#' difference of gain timing estimates of the same segments in two regions of
#' a tumour and of the z-test p-values.
#'
#' @param comparison output of [compare_region_timing()].
#' @param bins number of histogram bins.
#' @return a patchwork object.
#' @export
plot_region_comparison <- function(comparison, bins = 20) {
  p1 <- ggplot2::ggplot(comparison, ggplot2::aes(x = .data$difference)) +
    ggplot2::geom_histogram(bins = bins, colour = "darkgrey", fill = "darkblue") +
    ggplot2::xlab("gain timing difference") + ggplot2::theme_bw()
  p2 <- ggplot2::ggplot(comparison, ggplot2::aes(x = .data$p_value)) +
    ggplot2::geom_histogram(bins = bins, colour = "darkgrey", fill = "darkblue", boundary = 0) +
    ggplot2::geom_vline(xintercept = 0.05, colour = "darkred") +
    ggplot2::xlab("p-value") + ggplot2::theme_bw()
  p1 + p2
}

#' Plot mutation multiplicities of clusters on one segment
#'
#' Reproduces thesis Figure 6-17A: mutations of one haplotype on one segment,
#' coloured by mutation cluster, for all clusters together and for each
#' cluster separately.
#'
#' @param mutations phased mutations with a cluster column.
#' @param segID segment to plot.
#' @param haplotype haplotype to plot.
#' @param cluster_col name of the cluster column.
#' @param clusters clusters to show (default all with mutations on the
#'   segment).
#' @return a patchwork object.
#' @export
plot_cluster_segment <- function(mutations, segID, haplotype = "A", cluster_col = "cluster",
                                 clusters = NULL) {
  m <- mutations[mutations$segID == segID & !is.na(mutations$phase) & mutations$phase == haplotype, ]
  m$cl <- as.character(m[[cluster_col]])
  if (is.null(clusters)) clusters <- sort(unique(m$cl))
  m <- m[m$cl %in% clusters, ]
  m$cl <- factor(m$cl, levels = clusters)
  cols <- stats::setNames(gain_colours(max(as.numeric(clusters), length(clusters)))[
    if (all(!is.na(suppressWarnings(as.numeric(clusters))))) as.numeric(clusters) else seq_along(clusters)],
    clusters)
  ymax <- max(4, ceiling(max(m$mut_cpn)))
  panel <- function(d, show_x, legend = FALSE) {
    ggplot2::ggplot(d, ggplot2::aes(x = .data$pos, y = .data$mut_cpn, colour = .data$cl)) +
      ggplot2::geom_hline(yintercept = seq(0.5, ymax, 1), linetype = 2, colour = "#bababa") +
      ggplot2::geom_point(size = 1.2, show.legend = legend) +
      ggplot2::scale_colour_manual(name = "Cluster", values = cols, drop = FALSE) +
      ggplot2::ylim(0, ymax) + ggplot2::ylab("mutation copy number") +
      ggplot2::xlab(if (show_x) unique(m$chr) else "") + ggplot2::theme_bw() +
      if (!show_x) ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                  axis.ticks.x = ggplot2::element_blank()) else NULL
  }
  ps <- c(list(panel(m, FALSE, legend = TRUE)),
          lapply(seq_along(clusters), function(i) panel(m[m$cl == clusters[i], ], i == length(clusters))))
  patchwork::wrap_plots(ps, ncol = 1) + patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title = paste0(segID, " - haplotype ", haplotype))
}

# simple top-down tree layout: leaves spaced evenly, parents centred
.tree_layout <- function(edges, trunk) {
  kids <- function(n) edges[edges[, 1] == n, 2]
  pos <- list()
  next_x <- 0
  walk <- function(n, depth) {
    ch <- kids(n)
    if (length(ch) == 0) {
      x <- next_x
      next_x <<- next_x + 1
    } else {
      xs <- vapply(ch, function(k) walk(k, depth + 1), numeric(1))
      x <- mean(xs)
    }
    pos[[n]] <<- c(x = x, y = depth)
    x
  }
  walk(trunk, 1)
  data.frame(node = names(pos), x = vapply(pos, `[`, numeric(1), "x"),
             y = vapply(pos, `[`, numeric(1), "y"), stringsAsFactors = FALSE)
}

#' Plot copy number gains on a phylogenetic tree
#'
#' Reproduces thesis Figure 6-17B: the gains timed for one segment and
#' haplotype are drawn on the edges leading to the clusters in which they
#' occurred, at the position corresponding to their timing within the
#' subclonal expansion of the cluster.
#'
#' @param tree tree (see [as_phaseTime_tree()]).
#' @param timing output of [time_gains()] with a `cluster` column.
#' @param segID segment to show.
#' @param haplotype haplotype to show.
#' @param title plot title (default segment and haplotype).
#' @return a ggplot object.
#' @export
plot_timing_tree <- function(tree, timing, segID, haplotype = "A", title = NULL) {
  tr <- as_phaseTime_tree(tree)
  lay <- .tree_layout(tr$edges, tr$trunk)
  root <- lay[lay$node == tr$trunk, ]
  ed <- rbind(data.frame(parent = "0", child = tr$trunk),
              data.frame(parent = tr$edges[, 1], child = tr$edges[, 2]))
  ed$x0 <- ifelse(ed$parent == "0", root$x, lay$x[match(ed$parent, lay$node)])
  ed$y0 <- ifelse(ed$parent == "0", root$y - 0.6, lay$y[match(ed$parent, lay$node)])
  ed$x1 <- lay$x[match(ed$child, lay$node)]
  ed$y1 <- lay$y[match(ed$child, lay$node)]
  t <- timing[timing$segID == segID & timing$haplotype == haplotype & !is.na(timing$gain_timing), ]
  t <- t[t$cluster %in% ed$child, , drop = FALSE]
  if (nrow(t) > 0) {
    e <- ed[match(t$cluster, ed$child), ]
    # nodes are drawn with radius ~0.18 in y; place the mark on the visible part of the edge
    frac <- t$gain_timing
    t$mx <- e$x0 + (e$x1 - e$x0) * frac
    t$my <- e$y0 + (e$y1 - e$y0) * frac
    t$label <- paste0(round(100 * t$gain_timing, 1), "%")
  }
  lay$fill <- gain_colours(12)[((suppressWarnings(as.integer(lay$node)) - 1) %% 12) + 1]
  lay$fill[is.na(lay$fill)] <- "grey"
  if (is.null(title)) {
    ss <- split_seg_id(segID)
    title <- paste0(ss$chr, ":", ss$seg_start, "-", ss$seg_end, " - haplotype ", haplotype)
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = ed, ggplot2::aes(x = .data$x0, y = .data$y0, xend = .data$x1, yend = .data$y1),
                          linewidth = 0.8)
  if (nrow(t) > 0) {
    p <- p + ggplot2::geom_point(data = t, ggplot2::aes(x = .data$mx, y = .data$my), shape = 95, size = 8,
                                 colour = if (haplotype == "A") .cols$hapA_cn else .cols$hapB_cn) +
      ggplot2::geom_text(data = t, ggplot2::aes(x = .data$mx, y = .data$my, label = .data$label),
                         hjust = -0.25, size = 4)
  }
  p + ggplot2::geom_point(data = lay, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$fill),
                          shape = 21, size = 12, stroke = 0.8) +
    ggplot2::geom_text(data = lay, ggplot2::aes(x = .data$x, y = .data$y, label = .data$node), size = 5) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::ggtitle(title) + ggplot2::theme_void()
}

#' Plot phasing performance per sample
#'
#' Reproduces thesis Figure 6-7A-C: percentage of mutations phased to
#' haplotypes, percentage of all mutations phased incorrectly (estimated from
#' LOH regions) and percentage of phased mutations in LOH regions phased
#' incorrectly.
#'
#' @param summary output of [phasing_summary()] (several samples combined
#'   with `rbind`).
#' @return a patchwork object.
#' @export
plot_phasing_summary <- function(summary) {
  s <- summary[order(-summary$pct_phased), ]
  s$sample <- factor(s$sample, levels = s$sample)
  bar <- function(col, lab, show_x) {
    ggplot2::ggplot(s, ggplot2::aes(x = .data$sample, y = .data[[col]])) +
      ggplot2::geom_col(fill = "darkblue") + ggplot2::ylab(lab) + ggplot2::xlab("") +
      ggplot2::theme_bw() +
      if (show_x) ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) else
        ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }
  bar("pct_phased", "% mutations\nphased to haplotypes", FALSE) /
    bar("pct_incorrect", "% mutations\nincorrectly phased", FALSE) /
    bar("pct_incorrect_loh", "% mutations in LOH\nregions incorrectly phased", TRUE)
}
