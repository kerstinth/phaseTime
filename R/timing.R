# mutation timing, ancestral copy number states and copy number gain timing

# haplotype specific copy number of the haplotype a mutation was phased to
.hap_cn <- function(d) {
  ifelse(is.na(d$phase), NA_real_, ifelse(d$phase == "A", d$cnA, d$cnB))
}

#' Classify mutations as early or late relative to a copy number gain
#'
#' Thesis section 6.3.5. Only phased mutations on clonal (integer) gained
#' haplotype copy number segments can be timed. Mutation multiplicities are
#' rounded to the nearest integer: if the rounded multiplicity equals the
#' haplotype copy number the mutation occurred before the gain (`early`),
#' if it is smaller the mutation occurred after the gain (`late`).
#'
#' @param mutations phased mutations with `phase`, `cnA`, `cnB`, `mut_cpn`.
#' @param tolerance tolerance used to decide whether a copy number is an
#'   integer (clonal).
#' @return `mutations` with the columns `hap_cn` (copy number of the
#'   haplotype the mutation is on), `mut_cpn_round` and `timing` (factor with
#'   levels `early`, `late`, `not_possible`).
#' @export
classify_mutation_timing <- function(mutations, tolerance = 1e-6) {
  mutations <- as.data.frame(mutations)
  .require_cols(mutations, c("phase", "cnA", "cnB", "mut_cpn"), "mutations")
  mutations$hap_cn <- .hap_cn(mutations)
  mutations$mut_cpn_round <- round(mutations$mut_cpn)
  timeable <- !is.na(mutations$hap_cn) & mutations$hap_cn > 1 &
    .is_whole(mutations$hap_cn, tolerance)
  timing <- rep("not_possible", nrow(mutations))
  timing[timeable & mutations$mut_cpn_round >= round(mutations$hap_cn)] <- "early"
  timing[timeable & mutations$mut_cpn_round < round(mutations$hap_cn)] <- "late"
  mutations$timing <- factor(timing, levels = c("early", "late", "not_possible"))
  mutations
}

#' Assign mutations to (ancestral) copy number states
#'
#' Thesis section 6.3.5.1. On a gained haplotype, mutation multiplicities
#' cluster around integer values. Mutations at integer values below the
#' observed copy number define "ancestral" copy number states that existed
#' before a later gain. Every phased mutation on a clonal haplotype segment
#' with copy number > 1 is assigned to the closest integer multiplicity
#' (capped at the haplotype copy number). Mutations rounding to 0 and states
#' with fewer than `min_muts` mutations are removed (`cn_state = NA`).
#'
#' @param mutations phased mutations with `segID`, `phase`, `cnA`, `cnB`,
#'   `mut_cpn`.
#' @param min_muts minimum number of mutations for a copy number state to be
#'   kept (default 5).
#' @param by optional column name (e.g. `"cluster"`): the `min_muts` filter is
#'   applied per segment, haplotype and value of this column.
#' @param tolerance tolerance to decide whether a copy number is clonal.
#' @return `mutations` with columns `hap_cn` and `cn_state`.
#' @export
assign_cn_states <- function(mutations, min_muts = 5, by = NULL, tolerance = 1e-6) {
  mutations <- as.data.frame(mutations)
  .require_cols(mutations, c("segID", "phase", "cnA", "cnB", "mut_cpn"), "mutations")
  hap_cn <- .hap_cn(mutations)
  mutations$hap_cn <- hap_cn
  ok <- !is.na(hap_cn) & hap_cn > 1 & .is_whole(hap_cn, tolerance)
  state <- pmin(round(mutations$mut_cpn), round(hap_cn))
  state[!ok | state < 1] <- NA
  grp <- paste(mutations$segID, mutations$phase,
               if (is.null(by)) "" else mutations[[by]], state)
  n <- stats::ave(rep(1, length(state)), grp, FUN = length)
  state[n < min_muts] <- NA
  mutations$cn_state <- state
  mutations
}

#' Copy number gain timing from mutation counts per copy number state
#'
#' The core timing calculation of phaseTime (thesis sections 6.3.5.2 and
#' Figure 6-13). Given the mutation counts at each integer multiplicity
#' (copy number state) of a gained haplotype segment, the gain history is
#' reconstructed as a series of gains through the ancestral copy number
#' states `s_1 < s_2 < ... < s_k`. Mutations with the highest multiplicity
#' were acquired first (before the first gain), mutations with the lowest
#' multiplicity last. Because more copies of DNA accumulate more mutations,
#' the number of mutations acquired in each interval is divided by the number
#' of copies present while they were acquired. The timing of gain `g` is the
#' copy number adjusted number of mutations acquired before the gain divided
#' by all copy number adjusted mutations:
#'
#' \deqn{T_g = \sum_{early} (n_j / c_j) / \sum_{all} (n_j / c_j)}
#'
#' A value close to 0 means the gain happened early, close to 1 late in the
#' evolutionary time of the segment (or cluster). 95% confidence intervals use
#' the normal approximation \eqn{T \pm z \sqrt{T (1 - T) / N}} with `N` the
#' total adjusted number of mutations.
#'
#' @param states integer multiplicities (copy number states) observed.
#' @param counts number of mutations at each state.
#' @param copies number of DNA copies present when the mutations of each
#'   state were acquired. By default the ancestral copy number path is taken
#'   to be the sorted states, so that the mutations at the highest state were
#'   acquired when `s_1` copies were present, the next highest when `s_2`
#'   copies were present, etc. (Figure 6-13).
#' @param adjust `"copies_present"` (default, as described in the thesis and
#'   validated by simulation) divides counts by `copies`;
#'   `"multiplicity"` divides counts by the multiplicity of the state, which
#'   reproduces the original `segment.cluster.timing.fun()` script.
#' @param keep optional logical vector of length `length(states) - 1`
#'   indicating which boundaries between consecutive states (ordered from the
#'   highest multiplicity) are gains to be timed. Mutations between dropped
#'   boundaries are pooled into the same interval. Used to exclude gains that
#'   happened in a descendant clone.
#' @param conf confidence level of the interval.
#' @return data.frame with one row per gain: `gain_index` (1 = first gain),
#'   `copy_number_before_gain`, `copy_number_after_gain`,
#'   `number_early_muts_cn_adjusted`, `number_all_muts_cn_adjusted`,
#'   `gain_timing`, `se`, `loCI`, `hiCI`.
#' @export
#' @examples
#' # Figure 6-13: one mutation before the first gain (multiplicity 3), two
#' # mutations between the gains (multiplicity 2) and three after the second
#' # gain (multiplicity 1) give gain timings of 1/3 and 2/3.
#' gain_timing_from_counts(states = c(1, 2, 3), counts = c(3, 2, 1))
gain_timing_from_counts <- function(states, counts, copies = NULL,
                                    adjust = c("copies_present", "multiplicity"),
                                    keep = NULL, conf = 0.95) {
  adjust <- match.arg(adjust)
  o <- order(states, decreasing = TRUE)
  states <- states[o]
  counts <- counts[o]
  k <- length(states)
  empty <- data.frame(gain_index = integer(), copy_number_before_gain = numeric(),
                      copy_number_after_gain = numeric(),
                      number_early_muts_cn_adjusted = numeric(),
                      number_all_muts_cn_adjusted = numeric(), gain_timing = numeric(),
                      se = numeric(), loCI = numeric(), hiCI = numeric())
  if (k < 2) return(empty)
  if (is.null(copies)) {
    copies <- sort(states)  # s_1 ... s_k: copies present for states in decreasing order
  } else {
    copies <- copies[o]
  }
  if (is.null(keep)) keep <- rep(TRUE, k - 1)
  if (!any(keep)) return(empty)
  # interval membership of each state (ordered by decreasing multiplicity)
  interval <- cumsum(c(1, keep))
  if (adjust == "copies_present") {
    # copies present during an interval: those of its highest multiplicity state
    first_in_interval <- !duplicated(interval)
    factor <- copies[first_in_interval][interval]
  } else {
    factor <- states
  }
  factor[factor < 1] <- 1
  adj <- counts / factor
  total <- sum(adj)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  b <- which(keep)
  out <- lapply(seq_along(b), function(g) {
    early <- sum(adj[seq_len(b[g])])
    t <- early / total
    se <- sqrt(t * (1 - t) / total)
    data.frame(gain_index = g,
               copy_number_before_gain = copies[b[g]],
               copy_number_after_gain = copies[b[g] + 1],
               number_early_muts_cn_adjusted = early,
               number_all_muts_cn_adjusted = total,
               gain_timing = t, se = se,
               loCI = max(0, t - z * se), hiCI = min(1, t + z * se))
  })
  do.call(rbind, out)
}

#' Convert a tree to the edge list format used by phaseTime
#'
#' @param tree either a two column matrix / data.frame of edges (parent,
#'   child) or a CONIPHER / TRACERx tree object (a list with
#'   `graph_pyclone$default_tree` and `graph_pyclone$trunk`).
#' @param trunk the trunk (root) cluster; taken from the tree object if not
#'   given.
#' @return list with `edges` (character matrix with columns parent, child) and
#'   `trunk`.
#' @export
as_phaseTime_tree <- function(tree, trunk = NULL) {
  if (is.list(tree) && !is.data.frame(tree) && !is.null(tree$edges)) return(tree)
  if (is.list(tree) && !is.data.frame(tree) && !is.null(tree$graph_pyclone)) {
    trunk <- if (is.null(trunk)) tree$graph_pyclone$trunk else trunk
    tree <- tree$graph_pyclone$default_tree
  }
  edges <- as.matrix(tree)
  edges <- matrix(as.character(edges), ncol = 2, dimnames = list(NULL, c("parent", "child")))
  if (is.null(trunk)) {
    trunk <- setdiff(edges[, 1], edges[, 2])[1]
  }
  list(edges = edges, trunk = as.character(trunk))
}

.descendants <- function(edges, node) {
  kids <- edges[edges[, 1] == node, 2]
  if (length(kids) == 0) return(character())
  unique(c(kids, unlist(lapply(kids, function(k) .descendants(edges, k)))))
}

#' Estimate haplotype specific copy number gain timing
#'
#' Times every gain on every clonal, gained, haplotype specific copy number
#' segment (thesis sections 6.3.5.2, 6.3.6 and 6.3.7):
#'
#' 1. Phased mutations are assigned to integer copy number states
#'    ([assign_cn_states()]), separately for each haplotype and (if given)
#'    each mutation cluster.
#' 2. The copy number path of the segment is given by the states observed
#'    across all clusters, which defines the number of copies present when
#'    each state's mutations were acquired.
#' 3. For each cluster with at least two states the gains between
#'    consecutive states are timed with [gain_timing_from_counts()].
#' 4. If a tree is given (or a trunk without a tree, in which case every
#'    other cluster is treated as a descendant of the trunk), a gain in an
#'    ancestral cluster is only kept if the same gain (same copy number
#'    before and after) is not also detected in one of its descendant
#'    clusters; otherwise it occurred in the descendant (Figure 6-17).
#'
#' With `clonal_only = TRUE` only gains in the trunk are returned (section
#' 6.3.6, "timing clonal copy number gains").
#'
#' @param mutations phased mutations (output of [phase_mutations()] after
#'   [annotate_mutations()]) with `chr`, `segID`, `seg_start`, `seg_end`,
#'   `phase`, `cnA`, `cnB`, `mut_cpn` and optionally a cluster column.
#' @param cluster_col name of the column with mutation cluster assignments,
#'   or `NULL` to treat all mutations as one (clonal) cluster.
#' @param tree optional tree (see [as_phaseTime_tree()]).
#' @param trunk trunk cluster. Taken from `tree` if not given.
#' @param clonal_only only return gains timed in the trunk cluster.
#' @param min_muts minimum number of mutations per copy number state.
#' @param adjust adjustment of mutation counts, see
#'   [gain_timing_from_counts()].
#' @param include_untimed also return one row per segment and haplotype that
#'   could not be timed, with the reason in column `status`.
#' @param conf confidence level.
#' @param sample optional sample name added to the output.
#' @return data.frame with one row per timed gain: `sample`, `segID`, `chr`,
#'   `seg_start`, `seg_end`, `haplotype`, `hap_cn`, `cluster`, `gain_index`,
#'   `copy_number_before_gain`, `copy_number_after_gain`,
#'   `number_early_muts_cn_adjusted`, `number_all_muts_cn_adjusted`,
#'   `n_mutations`, `gain_timing`, `se`, `loCI`, `hiCI`, `status`.
#' @export
time_gains <- function(mutations, cluster_col = NULL, tree = NULL, trunk = NULL,
                       clonal_only = !is.null(trunk) || !is.null(tree),
                       min_muts = 5, adjust = c("copies_present", "multiplicity"),
                       include_untimed = FALSE, conf = 0.95, sample = NULL) {
  adjust <- match.arg(adjust)
  mutations <- as.data.frame(mutations)
  .require_cols(mutations, c("chr", "segID", "seg_start", "seg_end", "phase", "cnA", "cnB", "mut_cpn"),
                "mutations")
  if (is.null(sample)) {
    sample <- if ("sample" %in% names(mutations)) as.character(mutations$sample[1]) else NA_character_
  }
  edges <- NULL
  if (!is.null(tree)) {
    tr <- as_phaseTime_tree(tree, trunk)
    edges <- tr$edges
    trunk <- tr$trunk
  }
  if (is.null(cluster_col)) {
    mutations$.cluster <- if (is.null(trunk)) "all" else as.character(trunk)
    if (is.null(trunk)) trunk <- NULL
  } else {
    .require_cols(mutations, cluster_col, "mutations")
    mutations$.cluster <- as.character(mutations[[cluster_col]])
  }
  descendants_of <- function(cl) {
    if (!is.null(edges)) return(.descendants(edges, cl))
    if (!is.null(trunk) && cl == trunk) return(setdiff(unique(mutations$.cluster), trunk))
    character()
  }

  mutations <- assign_cn_states(mutations, min_muts = min_muts, by = ".cluster")
  segs <- unique(mutations[, c("segID", "chr", "seg_start", "seg_end", "cnA", "cnB")])
  segs <- segs[!duplicated(segs$segID), ]

  untimed <- function(s, h, hcn, status) {
    data.frame(sample = sample, segID = s$segID, chr = s$chr, seg_start = s$seg_start,
               seg_end = s$seg_end, haplotype = h, hap_cn = hcn, cluster = NA_character_,
               gain_index = NA_integer_, copy_number_before_gain = NA_real_,
               copy_number_after_gain = NA_real_, number_early_muts_cn_adjusted = NA_real_,
               number_all_muts_cn_adjusted = NA_real_, n_mutations = NA_integer_,
               gain_timing = NA_real_, se = NA_real_, loCI = NA_real_, hiCI = NA_real_,
               status = status, stringsAsFactors = FALSE)
  }

  res <- list()
  for (i in seq_len(nrow(segs))) {
    s <- segs[i, ]
    for (h in c("A", "B")) {
      hcn <- if (h == "A") s$cnA else s$cnB
      d <- mutations[mutations$segID == s$segID & !is.na(mutations$phase) &
                       mutations$phase == h, , drop = FALSE]
      status <- if (is.na(hcn) || hcn <= 1) "no_gain" else
        if (!.is_whole(hcn, 1e-6)) "subclonal_copy_number" else
          if (nrow(d) == 0) "no_phased_mutations" else NA
      d <- d[!is.na(d$cn_state), , drop = FALSE]
      if (is.na(status) && length(unique(d$cn_state)) < 2) status <- "fewer_than_two_states"
      if (!is.na(status)) {
        if (include_untimed) res[[length(res) + 1]] <- untimed(s, h, hcn, status)
        next
      }
      all_states <- sort(unique(d$cn_state))
      copies_map <- stats::setNames(rev(all_states), all_states)
      clusters <- unique(d$.cluster)
      # gains (as copy number transitions) detected within each cluster
      per_cluster <- lapply(clusters, function(cl) {
        st <- sort(unique(d$cn_state[d$.cluster == cl]), decreasing = TRUE)
        if (length(st) < 2) return(character())
        cp <- copies_map[as.character(st)]
        paste(cp[-length(cp)], cp[-1], sep = ">")
      })
      names(per_cluster) <- clusters
      out_h <- list()
      for (cl in clusters) {
        if (clonal_only && !is.null(trunk) && cl != trunk) next
        dc <- d[d$.cluster == cl, , drop = FALSE]
        tab <- table(dc$cn_state)
        st <- as.numeric(names(tab))
        if (length(st) < 2) next
        trans <- per_cluster[[cl]]
        desc <- intersect(descendants_of(cl), clusters)
        desc_trans <- unlist(per_cluster[desc])
        keep <- !(trans %in% desc_trans)
        # ordering of keep must follow decreasing states
        tm <- gain_timing_from_counts(st, as.numeric(tab), copies = copies_map[as.character(st)],
                                      adjust = adjust, keep = keep, conf = conf)
        if (nrow(tm) == 0) next
        out_h[[cl]] <- cbind(
          data.frame(sample = sample, segID = s$segID, chr = s$chr, seg_start = s$seg_start,
                     seg_end = s$seg_end, haplotype = h, hap_cn = hcn,
                     cluster = if (is.null(cluster_col) && is.null(trunk)) NA_character_ else cl,
                     stringsAsFactors = FALSE),
          tm[, c("gain_index", "copy_number_before_gain", "copy_number_after_gain",
                 "number_early_muts_cn_adjusted", "number_all_muts_cn_adjusted")],
          n_mutations = nrow(dc),
          tm[, c("gain_timing", "se", "loCI", "hiCI")],
          status = "timed", stringsAsFactors = FALSE)
      }
      if (length(out_h) == 0) {
        if (include_untimed) res[[length(res) + 1]] <- untimed(s, h, hcn, "no_gain_in_trunk")
      } else {
        res <- c(res, out_h)
      }
    }
  }
  if (length(res) == 0) return(untimed(segs[0, ], character(), numeric(), character()))
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}

#' Compare gain timing between two tumour regions
#'
#' Gains timed independently in two regions of the same tumour should have
#' similar timing estimates (thesis Figure 6-16B, C). Gains are matched by
#' haplotype, gain index and segment (identical segment, or segments with a
#' reciprocal overlap of at least `min_overlap`). The difference of the
#' timing estimates is tested with a two-sided z-test.
#'
#' @param timing output of [time_gains()] for several samples (column
#'   `sample`), combined with `rbind`.
#' @param tumour optional vector of tumour identifiers, one per row of
#'   `timing`. Defaults to the sample name up to the first underscore.
#' @param min_overlap minimum reciprocal overlap of segments.
#' @return data.frame with one row per matched gain: `tumour`, `haplotype`,
#'   `gain_index`, `sample1`, `segID1`, `timing1`, `sample2`, `segID2`,
#'   `timing2`, `difference`, `z`, `p_value`.
#' @export
compare_region_timing <- function(timing, tumour = NULL, min_overlap = 0.5) {
  timing <- as.data.frame(timing)
  timing <- timing[!is.na(timing$gain_timing), , drop = FALSE]
  if (is.null(tumour)) tumour <- sub("_.*$", "", timing$sample)
  timing$tumour <- tumour
  out <- list()
  for (tu in unique(timing$tumour)) {
    d <- timing[timing$tumour == tu, ]
    samples <- unique(d$sample)
    if (length(samples) < 2) next
    pairs <- utils::combn(samples, 2, simplify = FALSE)
    for (p in pairs) {
      a <- d[d$sample == p[1], ]
      b <- d[d$sample == p[2], ]
      for (i in seq_len(nrow(a))) {
        cand <- b[b$haplotype == a$haplotype[i] & b$gain_index == a$gain_index[i] &
                    b$chr == a$chr[i], , drop = FALSE]
        if (nrow(cand) == 0) next
        ov <- pmin(cand$seg_end, a$seg_end[i]) - pmax(cand$seg_start, a$seg_start[i])
        recip <- pmin(ov / (cand$seg_end - cand$seg_start), ov / (a$seg_end[i] - a$seg_start[i]))
        j <- which(recip >= min_overlap)
        if (length(j) == 0) next
        j <- j[which.max(recip[j])]
        diff <- abs(a$gain_timing[i] - cand$gain_timing[j])
        se <- sqrt(a$se[i]^2 + cand$se[j]^2)
        z <- if (se > 0) diff / se else NA_real_
        out[[length(out) + 1]] <- data.frame(
          tumour = tu, haplotype = a$haplotype[i], gain_index = a$gain_index[i],
          sample1 = p[1], segID1 = a$segID[i], timing1 = a$gain_timing[i],
          sample2 = p[2], segID2 = cand$segID[j], timing2 = cand$gain_timing[j],
          difference = diff, z = z, p_value = if (is.na(z)) NA_real_ else 2 * stats::pnorm(-abs(z)),
          stringsAsFactors = FALSE)
      }
    }
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}
