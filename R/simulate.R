# simulations: phasing feasibility, gain timing validation and a synthetic
# multi-region tumour for examples

#' Simulate the fraction of mutations that can be phased to SNPs
#'
#' Reproduces the simulation of thesis Figure 6-4B: random genomic positions
#' of `n_mutations` mutations and `n_snps` heterozygous SNPs are drawn and the
#' genome is tiled into windows of the sequencing read length. A mutation can
#' be phased if a SNP lies in the same window, i.e. on the same read.
#'
#' @param read_lengths read lengths to test.
#' @param n_mutations number of mutations (default 50,000, the median number
#'   per tumour region in the TRACERx long-read cohort).
#' @param n_snps number of heterozygous SNPs (default 2,250,000).
#' @param genome_size genome size in bp.
#' @param seed optional random seed.
#' @return data.frame with `read_length`, `number_mutations`,
#'   `number_phased_muts` and `pct_phased`.
#' @export
#' @examples
#' simulate_snp_phasing(read_lengths = c(150, 5000), n_mutations = 5000, seed = 1)
simulate_snp_phasing <- function(read_lengths = c(100, 150, 200, 300, 500, 1000, 2000, 5000, 7500, 10000),
                                 n_mutations = 50000, n_snps = 2250000, genome_size = 3054815472,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  out <- lapply(read_lengths, function(rl) {
    snp <- sample.int(genome_size, n_snps)
    mut <- sample.int(genome_size, n_mutations)
    # mutations may not sit on a SNP position
    while (any(clash <- mut %in% snp)) mut[clash] <- sample.int(genome_size, sum(clash))
    phased <- sum((mut %/% rl) %in% (snp %/% rl))
    data.frame(read_length = rl, number_mutations = n_mutations, number_phased_muts = phased,
               pct_phased = 100 * phased / n_mutations)
  })
  do.call(rbind, out)
}

#' Plot the phasing simulation
#'
#' Bar chart of thesis Figure 6-4B.
#' @param sim output of [simulate_snp_phasing()].
#' @return a ggplot object.
#' @export
plot_snp_phasing_simulation <- function(sim) {
  sim$read_length <- factor(sim$read_length, levels = sim$read_length)
  ggplot2::ggplot(sim, ggplot2::aes(x = .data$read_length, y = .data$pct_phased)) +
    ggplot2::geom_col(fill = "darkblue") +
    ggplot2::ylab("% mutations phased to\nheterozygous SNPs") + ggplot2::xlab("read length") +
    ggplot2::theme_bw()
}

#' Simulate mutation accumulation on a segment with copy number gains
#'
#' Simulation framework of thesis section 6.3.5.3. A DNA segment of
#' `segment_size` bp is followed over cell divisions. At every division
#' mutations are drawn from a Poisson distribution with
#' \eqn{\lambda = \mu \cdot L \cdot CN}, where \eqn{\mu} is the mutation rate
#' per bp, `L` the segment size and `CN` the current copy number. After the
#' number of divisions given in `timepoints` the copy number changes to the
#' next value of `cn_steps`.
#'
#' @param timepoints cumulative numbers of cell divisions at which a gain
#'   happens; the last value is the time of sampling.
#' @param cn_steps copy number in each interval; must have length
#'   `length(timepoints) + 1` (the last value is the copy number at
#'   sampling).
#' @param segment_size segment size in bp.
#' @param mutation_rate mutation rate per bp per cell division.
#' @return data.frame with one row per interval: `interval_index`,
#'   `mutation_rate_bp`, `segment_size`, `copy_number_before_gain`,
#'   `copy_number_after_gain`, `number_cell_divisions_since_last_gain`,
#'   `number_muts_accumulated`, `adjustment_factor`, `number_muts_cn_adjusted`.
#' @export
#' @examples
#' set.seed(1)
#' simulate_gain_history(timepoints = c(5000, 10000), cn_steps = c(1, 2, 2))
simulate_gain_history <- function(timepoints = c(1000, 2000, 5000, 10000), cn_steps = c(1, 2, 3, 5, 5),
                                  segment_size = 5e6, mutation_rate = 1.33e-8) {
  if (length(cn_steps) != length(timepoints) + 1) {
    stop("cn_steps must have length(timepoints) + 1 elements", call. = FALSE)
  }
  divisions <- diff(c(0, timepoints))
  muts <- vapply(seq_along(timepoints), function(i) {
    sum(stats::rpois(divisions[i], lambda = segment_size * cn_steps[i] * mutation_rate))
  }, numeric(1))
  data.frame(interval_index = seq_along(timepoints),
             mutation_rate_bp = mutation_rate,
             segment_size = segment_size,
             copy_number_before_gain = cn_steps[seq_along(timepoints)],
             copy_number_after_gain = cn_steps[seq_along(timepoints) + 1],
             number_cell_divisions_since_last_gain = divisions,
             number_muts_accumulated = muts,
             adjustment_factor = cn_steps[seq_along(timepoints)],
             number_muts_cn_adjusted = muts / cn_steps[seq_along(timepoints)])
}

#' Gain timing of a simulated gain history
#'
#' Applies the phaseTime timing estimate to a simulated history and compares
#' it with the ground truth (the fraction of cell divisions after which each
#' gain occurred).
#'
#' @param sim output of [simulate_gain_history()].
#' @param conf confidence level.
#' @return data.frame with one row per gain: `gain_index`,
#'   `mutation_rate_bp`, `segment_size`, `copy_number_before_gain`,
#'   `copy_number_after_gain`, `number_early_muts_cn_adjusted`,
#'   `number_all_muts_cn_adjusted`, `gain_timing`, `loCI`, `hiCI`,
#'   `theoretical_gain_timing`.
#' @export
simulated_gain_timing <- function(sim, conf = 0.95) {
  n <- nrow(sim)
  if (n < 2) stop("need at least one gain (two intervals)", call. = FALSE)
  total <- sum(sim$number_muts_cn_adjusted)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  early <- cumsum(sim$number_muts_cn_adjusted)[-n]
  t <- early / total
  err <- z * sqrt(t * (1 - t) / total)
  data.frame(gain_index = seq_len(n - 1),
             mutation_rate_bp = sim$mutation_rate_bp[1],
             segment_size = sim$segment_size[1],
             copy_number_before_gain = sim$copy_number_before_gain[-n],
             copy_number_after_gain = sim$copy_number_after_gain[-n],
             number_early_muts_cn_adjusted = early,
             number_all_muts_cn_adjusted = total,
             gain_timing = t,
             loCI = pmax(0, t - err), hiCI = pmin(1, t + err),
             theoretical_gain_timing = cumsum(sim$number_cell_divisions_since_last_gain)[-n] /
               sum(sim$number_cell_divisions_since_last_gain))
}

#' Plot a simulated gain timing scenario
#'
#' Reproduces thesis Figure 6-14A-C/E-G. Left: simulated mutations on the
#' segment at their final multiplicity (mutations acquired in the first
#' interval have the highest multiplicity). Right: calculated gain timing
#' (coloured bars with transparent 95% confidence intervals) and ground truth
#' gain timing (black outlined lines).
#'
#' @param sim output of [simulate_gain_history()].
#' @param timing output of [simulated_gain_timing()]; computed if `NULL`.
#' @param max_points maximum number of mutations drawn per interval.
#' @return a patchwork object.
#' @export
plot_simulated_gain_timing <- function(sim, timing = NULL, max_points = 1500) {
  if (is.null(timing)) timing <- simulated_gain_timing(sim)
  n <- nrow(sim)
  seg_size <- sim$segment_size[1]
  # multiplicity of mutations from each interval, see gain_timing_from_counts()
  mult <- rev(sim$copy_number_before_gain)
  muts <- do.call(rbind, lapply(seq_len(n), function(i) {
    k <- min(sim$number_muts_accumulated[i], max_points)
    if (k == 0) return(NULL)
    data.frame(pos = stats::runif(k, 0, seg_size), multiplicity = jitter(rep(mult[i], k), amount = 0.08),
               state = mult[i])
  }))
  muts$col <- factor(ifelse(match(muts$state, sort(unique(muts$state))) %% 2 == 1, "one", "two"))
  p1 <- ggplot2::ggplot(muts, ggplot2::aes(x = .data$pos, y = .data$multiplicity, colour = .data$col)) +
    ggplot2::geom_point(size = 0.3, alpha = 0.6) +
    ggplot2::scale_colour_manual(values = c(one = .cols$state_one, two = .cols$state_two), guide = "none") +
    ggplot2::ylab("mutation copy number") + ggplot2::xlab("copy number segment") +
    ggplot2::ylim(0, max(mult) + 0.5) + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  timing$gain <- factor(timing$gain_index)
  cols <- gain_colours(nrow(timing))
  p2 <- ggplot2::ggplot(timing) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 1, ymin = .data$loCI, ymax = .data$hiCI, fill = .data$gain),
                       alpha = 0.2) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 1, y = .data$gain_timing, yend = .data$gain_timing,
                                       colour = .data$gain), linewidth = 3) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 1, y = .data$theoretical_gain_timing,
                                       yend = .data$theoretical_gain_timing), colour = "#525252", linewidth = 1.6) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 1, y = .data$theoretical_gain_timing,
                                       yend = .data$theoretical_gain_timing, colour = .data$gain), linewidth = 0.7) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::scale_colour_manual(name = "gain order\n(bar: calculated,\noutlined: ground truth)", values = cols) +
    .timing_y_scale() + ggplot2::ylab("gain timing") + ggplot2::xlab("copy number segment") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
  (p1 + p2) + patchwork::plot_annotation(title = paste0("mutation rate: ", signif(sim$mutation_rate_bp[1], 3)))
}

#' Simulate many random gain histories
#'
#' Reproduces the validation of thesis Figure 6-14D and H: random numbers of
#' gains (1 to `max_gains`), random copy numbers (up to `max_cn`) and random
#' gain times are simulated for each mutation rate and the calculated gain
#' timing is compared with the ground truth.
#'
#' @param n_sim number of simulations per mutation rate.
#' @param mutation_rates mutation rates per bp per division.
#' @param segment_size segment size in bp.
#' @param n_divisions total number of cell divisions.
#' @param max_gains maximum number of gains.
#' @param max_cn maximum copy number.
#' @param seed optional random seed.
#' @return data.frame of [simulated_gain_timing()] results with the columns
#'   `iteration`, `abs_difference` and `ci_width` added.
#' @export
simulate_gain_timing_accuracy <- function(n_sim = 1000, mutation_rates = c(1.33e-8, 1.33e-9, 1.33e-10),
                                          segment_size = 5e6, n_divisions = 10000, max_gains = 5,
                                          max_cn = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  out <- vector("list", n_sim * length(mutation_rates))
  k <- 0
  for (i in seq_len(n_sim)) {
    for (mu in mutation_rates) {
      n_gains <- sample.int(max_gains, 1)
      cn <- sort(sample(2:max_cn, n_gains))
      cn_steps <- c(1, cn, max(cn))
      tp <- c(sort(sample.int(n_divisions - 1, n_gains)), n_divisions)
      sim <- simulate_gain_history(tp, cn_steps, segment_size, mu)
      tm <- simulated_gain_timing(sim)
      tm$iteration <- i
      k <- k + 1
      out[[k]] <- tm
    }
  }
  out <- do.call(rbind, out)
  out$abs_difference <- abs(out$gain_timing - out$theoretical_gain_timing)
  out$ci_width <- out$hiCI - out$loCI
  out
}

#' Plot the accuracy of simulated gain timing
#'
#' Reproduces thesis Figure 6-14D (absolute difference between calculated and
#' ground truth gain timing at the highest mutation rate) and 6-14H (width of
#' the confidence intervals per mutation rate).
#'
#' @param acc output of [simulate_gain_timing_accuracy()].
#' @param bins number of histogram bins.
#' @return a patchwork object.
#' @export
plot_gain_timing_accuracy <- function(acc, bins = 30) {
  top <- acc[acc$mutation_rate_bp == max(acc$mutation_rate_bp), ]
  p1 <- ggplot2::ggplot(top, ggplot2::aes(x = .data$abs_difference)) +
    ggplot2::geom_histogram(bins = bins, colour = "darkgrey", fill = "darkblue") +
    ggplot2::xlab("absolute difference") + ggplot2::theme_bw() +
    ggplot2::ggtitle(paste0("mutation rate ", signif(max(acc$mutation_rate_bp), 3)))
  acc$rate <- factor(signif(acc$mutation_rate_bp, 3),
                     levels = signif(sort(unique(acc$mutation_rate_bp), decreasing = TRUE), 3))
  p2 <- ggplot2::ggplot(acc, ggplot2::aes(x = .data$ci_width)) +
    ggplot2::geom_histogram(bins = bins, colour = "darkgrey", fill = "darkblue") +
    ggplot2::facet_wrap(~rate, nrow = 1) + ggplot2::xlab("CI width") + ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"))
  p1 + p2 + patchwork::plot_layout(widths = c(1, 2.2))
}

#' Copy number above which single-copy mutations are hard to detect
#'
#' In highly amplified regions a mutation present on one copy has a low VAF
#' and may be missed by mutation callers, which leads to an underestimation
#' of late mutations (thesis section 6.4). This returns, for each purity, the
#' largest total copy number at which a mutation with multiplicity 1 still
#' has a VAF of at least `min_vaf`.
#'
#' @param purity tumour purity (vector).
#' @param min_vaf minimum detectable VAF.
#' @param normal_cn copy number of the normal cells.
#' @return data.frame with `purity` and `max_copy_number`.
#' @export
#' @examples
#' detection_limit_cn(purity = c(0.3, 0.5, 0.8))
detection_limit_cn <- function(purity, min_vaf = 0.05, normal_cn = 2) {
  data.frame(purity = purity,
             max_copy_number = floor((purity / min_vaf - normal_cn * (1 - purity)) / purity))
}

# ---------------------------------------------------------------------------
# synthetic multi-region tumour

# simulate the copies of one haplotype of one segment along the tumour tree.
# events: data.frame(cluster, time (fraction of branch), to_cn)
.sim_haplotype <- function(branches, events, len_mb, rate, id_offset) {
  end_copies <- list()
  muts <- list()
  next_id <- id_offset
  for (b in seq_len(nrow(branches))) {
    cl <- branches$cluster[b]
    par <- branches$parent[b]
    copies <- if (is.na(par)) list(integer(0)) else end_copies[[par]]
    ev <- events[events$cluster == cl, , drop = FALSE]
    ev <- ev[order(ev$time), , drop = FALSE]
    times <- c(0, ev$time, 1)
    for (j in seq_len(length(times) - 1)) {
      dt <- (times[j + 1] - times[j]) * branches$duration[b]
      ncp <- length(copies)
      if (ncp > 0 && dt > 0) {
        n <- stats::rpois(1, rate * len_mb * ncp * dt)
        if (n > 0) {
          ids <- next_id + seq_len(n)
          next_id <- next_id + n
          on <- sample.int(ncp, n, replace = TRUE)
          for (cpi in unique(on)) copies[[cpi]] <- c(copies[[cpi]], ids[on == cpi])
          muts[[length(muts) + 1]] <- data.frame(id = ids, cluster = cl)
        }
      }
      if (j <= nrow(ev)) {
        to <- ev$to_cn[j]
        cur <- length(copies)
        if (to > cur) {
          if (to == 2 * cur) {
            copies <- c(copies, copies)  # doubling of all copies
          } else {
            dup <- sample.int(cur, to - cur, replace = (to - cur) > cur)
            copies <- c(copies, copies[dup])
          }
        } else if (to < cur) {
          copies <- copies[sort(sample.int(cur, to))]
        }
      }
    }
    end_copies[[cl]] <- copies
  }
  list(muts = if (length(muts)) do.call(rbind, muts) else data.frame(id = integer(), cluster = character()),
       end = end_copies, next_id = next_id)
}

.final_cn <- function(events, path) {
  cn <- 1
  for (cl in path) {
    ev <- events[events$cluster == cl, , drop = FALSE]
    if (nrow(ev)) cn <- ev$to_cn[which.max(ev$time)]
  }
  cn
}

#' Simulate a multi-region tumour for testing and examples
#'
#' Generates a synthetic tumour with two tumour regions, sequenced by short
#' reads (mutation and SNP read counts, copy number) and long reads (reads
#' linking mutations to heterozygous SNPs). Every haplotype of every copy
#' number segment is followed along the tumour phylogeny (1 -> 2 -> 3 and
#' 1 -> 4 -> 5): mutations are acquired on individual DNA copies at a rate
#' proportional to the number of copies, and copy number gains duplicate DNA
#' copies including the mutations they carry. Cluster 2 is clonal in region 1
#' and cluster 4 in region 2; clusters 3 and 5 are subclonal. The returned
#' ground truth contains the simulated timing of every gain as the fraction
#' of the evolutionary time of the cluster in which it occurred.
#'
#' The tumour contains two hand-picked example segments: chr5 (haplotype A
#' gained 1 -> 2 -> 3 in the trunk) and chr14 (haplotype A gained in the
#' trunk and doubled again in cluster 2, haplotype B lost), resembling the
#' examples in Figures 6-10 and 6-17 of the thesis.
#'
#' @param tumour tumour identifier; regions are named `<tumour>_R1`,
#'   `<tumour>_R2`.
#' @param purity purity of the two regions.
#' @param mutation_rate mutations per Mb per DNA copy per unit of trunk
#'   evolutionary time.
#' @param sr_depth,lr_depth mean short-read and long-read depth.
#' @param snp_density heterozygous SNPs per Mb used for BAF.
#' @param linked_snps mean number of SNPs on the long reads of a mutation.
#' @param chr_sizes chromosome sizes.
#' @param seed random seed.
#' @return a list with elements `mutations`, `cn`, `snps`, `linkage`,
#'   `purity`, `tree` and `truth`.
#' @export
#' @examples
#' sim <- simulate_tumour(seed = 1)
#' names(sim)
simulate_tumour <- function(tumour = "SIM001", purity = c(0.7, 0.55), mutation_rate = 4,
                            sr_depth = 200, lr_depth = 17, snp_density = 15, linked_snps = 1.2,
                            chr_sizes = hg38_chr_sizes, seed = 1) {
  set.seed(seed)
  regions <- paste0(tumour, "_R", 1:2)
  names(purity) <- regions
  branches <- data.frame(cluster = c("1", "2", "3", "4", "5"),
                         parent = c(NA, "1", "2", "1", "4"),
                         duration = c(1, 0.35, 0.25, 0.35, 0.25), stringsAsFactors = FALSE)
  tree <- list(edges = cbind(parent = c("1", "2", "1", "4"), child = c("2", "3", "4", "5")), trunk = "1")
  # clonal lineage and subclonal clusters (with CCF) of each region
  lineage <- list(c("1", "2"), c("1", "4"))
  subclonal <- list(c("3" = 0.4), c("5" = 0.5))
  names(lineage) <- names(subclonal) <- regions

  # segments
  segs <- do.call(rbind, lapply(names(chr_sizes), function(ch) {
    n <- sample(1:3, 1, prob = c(0.4, 0.4, 0.2))
    br <- sort(round(stats::runif(n - 1, 0.15, 0.85) * chr_sizes[ch]))
    data.frame(chr = ch, seg_start = c(1, br + 1), seg_end = c(br, chr_sizes[ch]))
  }))
  ev <- function(cluster = character(), time = numeric(), to_cn = numeric()) {
    data.frame(cluster = cluster, time = time, to_cn = to_cn, stringsAsFactors = FALSE)
  }
  scen_A <- function() {
    s <- sample(c("none", "gain2", "gain3", "double", "sub_R2"), 1, prob = c(0.3, 0.38, 0.1, 0.1, 0.12))
    t <- sort(round(stats::runif(2, 0.05, 0.95), 3))
    switch(s, none = ev(), gain2 = ev("1", t[1], 2), gain3 = ev(c("1", "1"), t, c(2, 3)),
           double = ev(c("1", "1"), t, c(2, 4)), sub_R2 = ev("4", t[1], 2))
  }
  scen_B <- function(evA) {
    s <- sample(c("none", "loss", "gain2"), 1, prob = c(0.6, 0.22, 0.18))
    t <- round(stats::runif(1, 0.05, 0.95), 3)
    out <- switch(s, none = ev(), loss = ev("1", t, 0), gain2 = ev("1", t, 2))
    # keep haplotype A the major haplotype in every region
    for (r in regions) {
      if (.final_cn(out, lineage[[r]]) > .final_cn(evA, lineage[[r]])) out <- ev()
    }
    out
  }
  events <- vector("list", nrow(segs))
  for (i in seq_len(nrow(segs))) {
    a <- scen_A()
    events[[i]] <- list(A = a, B = scen_B(a))
  }
  # hand-picked example segments: one segment spanning chr5 and chr14
  segs$idx <- seq_len(nrow(segs))
  for (ch in c("chr5", "chr14")) {
    w <- which(segs$chr == ch)
    if (length(w) > 1) segs <- segs[-w[-1], ]
    segs$seg_end[segs$chr == ch] <- chr_sizes[ch]
  }
  events <- events[segs$idx]
  segs$idx <- NULL
  rownames(segs) <- NULL
  i5 <- which(segs$chr == "chr5")
  i14 <- which(segs$chr == "chr14")
  events[[i5]] <- list(A = ev(c("1", "1"), c(0.3, 0.65), c(2, 3)), B = ev())
  events[[i14]] <- list(A = ev(c("1", "2"), c(0.442, 0.3), c(2, 4)), B = ev("1", 0.2, 0))

  # simulate every haplotype of every segment
  muts <- list()
  truth <- list()
  id <- 0
  for (i in seq_len(nrow(segs))) {
    len_mb <- (segs$seg_end[i] - segs$seg_start[i] + 1) / 1e6
    for (h in c("A", "B")) {
      e <- events[[i]][[h]]
      sim <- .sim_haplotype(branches, e, len_mb, mutation_rate, id)
      id <- sim$next_id
      sid <- seg_id(segs$chr[i], segs$seg_start[i], segs$seg_end[i])
      if (nrow(e)) {
        e <- e[order(e$cluster, e$time), ]
        cn_before <- stats::ave(e$to_cn, e$cluster, FUN = function(x) c(NA, x[-length(x)]))
        for (j in seq_len(nrow(e))) {
          if (is.na(cn_before[j])) {
            par_path <- c()
            cl <- e$cluster[j]
            p <- branches$parent[branches$cluster == cl]
            while (!is.na(p)) { par_path <- c(p, par_path); p <- branches$parent[branches$cluster == p] }
            cn_before[j] <- .final_cn(e, par_path)
          }
        }
        tr <- data.frame(segID = sid, haplotype = h, cluster = e$cluster,
                         copy_number_before_gain = cn_before,
                         copy_number_after_gain = e$to_cn, true_timing = e$time,
                         stringsAsFactors = FALSE)
        truth[[length(truth) + 1]] <- tr[tr$copy_number_after_gain > tr$copy_number_before_gain, ]
      }
      if (nrow(sim$muts) == 0) next
      m <- sim$muts
      m$chr <- segs$chr[i]
      m$pos <- round(stats::runif(nrow(m), segs$seg_start[i], segs$seg_end[i]))
      m$true_phase <- h
      m$seg <- i
      for (r in regions) {
        leaf <- utils::tail(lineage[[r]], 1)
        cnt <- tabulate(as.integer(unlist(sim$end[[leaf]])), nbins = id)
        m[[paste0("mult_", r)]] <- cnt[m$id]
        for (sc in names(subclonal[[r]])) {
          cs <- tabulate(as.integer(unlist(sim$end[[sc]])), nbins = id)
          w <- m$cluster == sc
          m[[paste0("mult_", r)]][w] <- cs[m$id[w]]
        }
      }
      muts[[length(muts) + 1]] <- m
    }
  }
  muts <- do.call(rbind, muts)
  muts$ref <- sample(c("A", "C", "G", "T"), nrow(muts), replace = TRUE)
  muts$alt <- vapply(muts$ref, function(x) sample(setdiff(c("A", "C", "G", "T"), x), 1), character(1))

  # copy number per region
  cn <- do.call(rbind, lapply(regions, function(r) {
    a <- vapply(events, function(e) .final_cn(e$A, lineage[[r]]), numeric(1))
    b <- vapply(events, function(e) .final_cn(e$B, lineage[[r]]), numeric(1))
    data.frame(sample = r, chr = segs$chr, seg_start = segs$seg_start, seg_end = segs$seg_end,
               nMajor = pmax(a, b), nMinor = pmin(a, b), true_cnA = a, true_cnB = b,
               stringsAsFactors = FALSE)
  }))

  # observed mutations (short reads) and long-read linkage per region
  mutations <- list()
  linkage <- list()
  for (r in regions) {
    p <- purity[r]
    cr <- cn[cn$sample == r, ]
    present <- muts$cluster %in% c(lineage[[r]], names(subclonal[[r]]))
    m <- muts[present & muts[[paste0("mult_", r)]] > 0, ]
    ccf <- ifelse(m$cluster %in% names(subclonal[[r]]), subclonal[[r]][m$cluster], 1)
    tot <- cr$true_cnA[m$seg] + cr$true_cnB[m$seg]
    true_vaf <- p * m[[paste0("mult_", r)]] * ccf / (p * tot + 2 * (1 - p))
    depth <- stats::rpois(nrow(m), sr_depth)
    var <- stats::rbinom(nrow(m), depth, pmin(true_vaf, 1))
    keep <- var >= 3
    m <- m[keep, ]
    obs <- data.frame(sample = r, chr = m$chr, pos = m$pos, ref = m$ref, alt = m$alt,
                      var_count = var[keep], depth = depth[keep], vaf = var[keep] / depth[keep],
                      cluster = m$cluster, true_phase = m$true_phase,
                      true_multiplicity = m[[paste0("mult_", r)]], stringsAsFactors = FALSE)
    mutations[[r]] <- obs
    # long reads
    lr_var <- stats::rbinom(nrow(m), stats::rpois(nrow(m), lr_depth), pmin(true_vaf[keep], 1))
    n_snp <- pmin(stats::rpois(nrow(m), linked_snps), 16)
    n_snp[lr_var == 0] <- 0
    w <- rep(seq_len(nrow(m)), n_snp)
    if (length(w)) {
      lk <- data.frame(sample = r, chr = m$chr[w], pos = m$pos[w],
                       snp_pos = pmax(1, m$pos[w] + round(stats::runif(length(w), -4000, 4000))),
                       mut_hap = m$true_phase[w], n_reads = lr_var[w], stringsAsFactors = FALSE)
      linkage[[r]] <- lk
    }
  }
  mutations <- do.call(rbind, mutations)
  rownames(mutations) <- NULL
  linkage <- do.call(rbind, linkage)
  # each linked SNP gets a germline haplotype for its alternative allele, shared across regions
  snp_key <- unique(paste(linkage$chr, linkage$snp_pos))
  snp_hap <- stats::setNames(sample(c("A", "B"), length(snp_key), replace = TRUE), snp_key)
  linkage$snp_alt_hap <- snp_hap[paste(linkage$chr, linkage$snp_pos)]
  covering <- stats::rbinom(nrow(linkage), linkage$n_reads, 0.85)
  p_alt <- ifelse(linkage$snp_alt_hap == linkage$mut_hap, 0.97, 0.03)
  linkage$mt_snp_alt <- stats::rbinom(nrow(linkage), covering, p_alt)
  linkage$mt_snp_ref <- covering - linkage$mt_snp_alt
  linkage <- linkage[, c("sample", "chr", "pos", "snp_pos", "mt_snp_alt", "mt_snp_ref")]
  rownames(linkage) <- NULL

  # heterozygous SNPs with BAF in each region
  bg <- do.call(rbind, lapply(names(chr_sizes), function(ch) {
    n <- stats::rpois(1, snp_density * chr_sizes[ch] / 1e6)
    data.frame(chr = ch, pos = sort(round(stats::runif(n, 1, chr_sizes[ch]))))
  }))
  bg$true_haplotype <- sample(c("A", "B"), nrow(bg), replace = TRUE)
  lsnp <- unique(data.frame(chr = sub(" .*$", "", snp_key), pos = as.numeric(sub("^.* ", "", snp_key)),
                            true_haplotype = unname(snp_hap), stringsAsFactors = FALSE))
  allsnp <- rbind(bg, lsnp)
  allsnp <- allsnp[!duplicated(paste(allsnp$chr, allsnp$pos)), ]
  sidx <- .match_segments(allsnp$chr, allsnp$pos, segs$chr, segs$seg_start, segs$seg_end)
  allsnp <- allsnp[!is.na(sidx), ]
  sidx <- sidx[!is.na(sidx)]
  snps <- do.call(rbind, lapply(regions, function(r) {
    p <- purity[r]
    cr <- cn[cn$sample == r, ]
    cnh <- ifelse(allsnp$true_haplotype == "A", cr$true_cnA[sidx], cr$true_cnB[sidx])
    tot <- cr$true_cnA[sidx] + cr$true_cnB[sidx]
    frac <- (p * cnh + (1 - p)) / (p * tot + 2 * (1 - p))
    d <- pmax(stats::rpois(nrow(allsnp), sr_depth), 1)
    data.frame(region = r, chr = allsnp$chr, pos = allsnp$pos, baf = stats::rbinom(nrow(allsnp), d, frac) / d,
               true_haplotype = allsnp$true_haplotype, stringsAsFactors = FALSE)
  }))
  snps <- snps[order(snps$region, match(snps$chr, names(chr_sizes)), snps$pos), ]
  rownames(snps) <- NULL

  truth <- do.call(rbind, truth)
  list(mutations = mutations, cn = cn, snps = snps, linkage = linkage, purity = purity,
       tree = tree, truth = truth)
}
