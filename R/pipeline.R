# convenience wrapper running the full phaseTime workflow

#' Run the phaseTime workflow on one or more tumour regions
#'
#' Runs all steps of phaseTime for every sample (tumour region):
#'
#' 1. multi-region phasing of SNPs ([phase_snps()]),
#' 2. assignment of copy number to haplotypes ([haplotype_copy_number()]),
#' 3. mutation multiplicity ([annotate_mutations()]),
#' 4. phasing of mutations ([phase_mutations()]),
#' 5. early / late mutation timing ([classify_mutation_timing()]) and copy
#'    number states ([assign_cn_states()]),
#' 6. copy number gain timing ([time_gains()]).
#'
#' @param mutations mutation table with `sample`, `chr`, `pos` and `vaf` (or
#'   `var_count` and `depth`), optionally a cluster column.
#' @param cn copy number per sample with `sample`, `chr`, `seg_start`,
#'   `seg_end` and major / minor copy number (see
#'   [haplotype_copy_number()]).
#' @param purity named numeric vector of purities (names = samples).
#' @param snps SNP BAF table with `region` (= sample), `chr`, `pos`, `baf`.
#' @param linkage long-read linkage table with `sample` and the columns
#'   required by [phase_mutations()].
#' @param cluster_col optional mutation cluster column.
#' @param tree optional tree ([as_phaseTime_tree()]).
#' @param clonal_only passed to [time_gains()]; by default only trunk gains
#'   are returned when a tree is given.
#' @param min_muts minimum mutations per copy number state.
#' @param adjust adjustment of mutation counts (see
#'   [gain_timing_from_counts()]).
#' @param output_dir optional directory to write result tables to.
#' @return list with `snps` (phased SNPs), `cn` (haplotype copy number),
#'   `mutations` (phased and timed mutations), `timing` (gain timing) and
#'   `phasing` (phasing summary), each combined over samples.
#' @export
run_phaseTime <- function(mutations, cn, purity, snps, linkage, cluster_col = NULL, tree = NULL,
                          clonal_only = !is.null(tree), min_muts = 5,
                          adjust = c("copies_present", "multiplicity"), output_dir = NULL) {
  adjust <- match.arg(adjust)
  .require_cols(mutations, "sample", "mutations")
  .require_cols(cn, "sample", "cn")
  .require_cols(linkage, "sample", "linkage")
  samples <- unique(mutations$sample)
  missing <- setdiff(samples, names(purity))
  if (length(missing)) stop("no purity given for: ", paste(missing, collapse = ", "), call. = FALSE)

  cn_long <- cn
  cn_long$region <- cn_long$sample
  snps_phased <- phase_snps(snps[snps$region %in% samples, ], cn = cn_long)

  res <- lapply(samples, function(s) {
    cn_s <- haplotype_copy_number(cn[cn$sample == s, ], snps = snps_phased[snps_phased$region == s, ])
    m <- annotate_mutations(mutations[mutations$sample == s, ], cn_s, purity[[s]])
    m <- phase_mutations(linkage[linkage$sample == s, ], snps_phased, m)
    m <- classify_mutation_timing(m)
    m <- assign_cn_states(m, min_muts = min_muts, by = cluster_col)
    tm <- time_gains(m, cluster_col = cluster_col, tree = tree, clonal_only = clonal_only,
                     min_muts = min_muts, adjust = adjust, sample = s)
    list(cn = cbind(sample = s, cn_s[, setdiff(names(cn_s), "sample")]), mutations = m, timing = tm,
         phasing = phasing_summary(m))
  })
  out <- list(snps = snps_phased,
              cn = do.call(rbind, lapply(res, `[[`, "cn")),
              mutations = do.call(rbind, lapply(res, `[[`, "mutations")),
              timing = do.call(rbind, lapply(res, `[[`, "timing")),
              phasing = do.call(rbind, lapply(res, `[[`, "phasing")))
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    for (nm in c("cn", "mutations", "timing", "phasing")) {
      utils::write.table(out[[nm]], file.path(output_dir, paste0("phaseTime_", nm, ".tsv")),
                         sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
  out
}
