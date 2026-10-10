#' phaseTime: haplotype-specific mutation phasing and copy number gain timing
#'
#' phaseTime integrates long-read and short-read whole genome sequencing of
#' (multi-region) tumour samples to
#'
#' * phase somatic mutations to haplotypes via heterozygous SNPs that are
#'   physically linked to the mutation on long reads
#'   ([select_phasing_snps()], [phase_snps()], [phase_mutations()]),
#' * time mutations relative to copy number gains and define ancestral copy
#'   number states ([classify_mutation_timing()], [assign_cn_states()]),
#' * estimate the timing of first and subsequent copy number gains for each
#'   haplotype, including gains placed on the edges of a phylogenetic tree
#'   ([time_gains()], [gain_timing_from_counts()]),
#' * validate the approach by simulation ([simulate_gain_history()],
#'   [simulate_gain_timing_accuracy()], [simulate_snp_phasing()]).
#'
#' [run_phaseTime()] runs the complete workflow; [simulate_tumour()] creates
#' an example data set. See `vignette("phaseTime")` for a tutorial.
#'
#' The methods are described in Chapter 6 ("Haplotype specific somatic copy
#' number gain timing using long-read sequencing of lung tumours") and
#' section 2.2.3 of Kerstin Thol's PhD thesis (University College London).
#'
#' @keywords internal
#' @importFrom ggplot2 .data
#' @importFrom patchwork wrap_plots plot_annotation plot_layout
"_PACKAGE"
