# phasing of SNPs and somatic mutations to haplotypes

#' Select heterozygous SNPs close to mutations
#'
#' Step 1 of mutation phasing (thesis section 6.3.1): for every mutation keep
#' the `max_snps` heterozygous SNPs closest to it within `max_distance` bp.
#' The resulting mutation-SNP pairs can then be passed to a pileup of the
#' long-read BAM (e.g. `samtools mpileup`) to count reads carrying both the
#' mutation and a SNP allele. This filter only reduces run time.
#'
#' @param mutations data.frame with `chr`, `pos`.
#' @param snps data.frame with `chr`, `pos` of heterozygous germline SNPs.
#' @param max_snps maximum number of SNPs per mutation (default 16).
#' @param max_distance maximum distance between mutation and SNP (default
#'   20,000 bp).
#' @return data.frame with one row per mutation-SNP pair: `chr`, `pos`,
#'   `snp_pos`, `distance`.
#' @export
#' @examples
#' muts <- data.frame(chr = "chr1", pos = c(1000, 50000))
#' snps <- data.frame(chr = "chr1", pos = c(500, 1500, 30000, 49000))
#' select_phasing_snps(muts, snps, max_snps = 2)
select_phasing_snps <- function(mutations, snps, max_snps = 16, max_distance = 20000) {
  .require_cols(mutations, c("chr", "pos"), "mutations")
  .require_cols(snps, c("chr", "pos"), "snps")
  mutations$chr <- add_chr_prefix(mutations$chr)
  snps$chr <- add_chr_prefix(snps$chr)
  out <- list()
  for (ch in unique(mutations$chr)) {
    mp <- unique(mutations$pos[mutations$chr == ch])
    sp <- sort(unique(snps$pos[snps$chr == ch]))
    if (length(sp) == 0) next
    lo <- findInterval(mp - max_distance - 1, sp) + 1
    hi <- findInterval(mp + max_distance, sp)
    res <- lapply(seq_along(mp), function(i) {
      if (hi[i] < lo[i]) return(NULL)
      cand <- sp[lo[i]:hi[i]]
      d <- abs(cand - mp[i])
      keep <- order(d)[seq_len(min(max_snps, length(d)))]
      data.frame(chr = ch, pos = mp[i], snp_pos = cand[keep], distance = d[keep])
    })
    out[[ch]] <- do.call(rbind, res)
  }
  out <- do.call(rbind, out)
  if (is.null(out)) {
    return(data.frame(chr = character(), pos = numeric(), snp_pos = numeric(), distance = numeric()))
  }
  rownames(out) <- NULL
  out
}

# create a common segmentation across regions (minimum consistent regions)
.disjoin_segments <- function(segments) {
  out <- lapply(split(segments, segments$chr), function(s) {
    br <- sort(unique(c(s$seg_start, s$seg_end + 1)))
    if (length(br) < 2) return(NULL)
    d <- data.frame(chr = s$chr[1], seg_start = br[-length(br)], seg_end = br[-1] - 1)
    # keep pieces covered by at least one input segment
    mid <- floor((d$seg_start + d$seg_end) / 2)
    covered <- vapply(mid, function(m) any(s$seg_start <= m & s$seg_end >= m), logical(1))
    d[covered, , drop = FALSE]
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' Phase heterozygous SNPs to haplotypes using B-allele frequencies
#'
#' Step 2 of mutation phasing (thesis section 6.3.2). In a segment with
#' allelic imbalance, the B-allele frequency (BAF) of a heterozygous SNP
#' deviates from 0.5: SNPs whose alternative allele lies on the more abundant
#' haplotype have BAF > 0.5. The haplotype carrying the alternative allele of
#' each SNP is therefore
#'
#' * `A` if BAF >= 0.5 and
#' * `B` if BAF < 0.5,
#'
#' where haplotype A is the more abundant haplotype of the reference region.
#'
#' With multi-region data (`cn` given for several regions), SNPs are phased
#' per segment using the tumour region with the strongest allelic imbalance
#' as the reference, and the resulting haplotypes are applied to all regions
#' ("multi-region reference phasing", Figure 6-5). SNPs in segments without
#' allelic imbalance in any region remain unphased (`haplotype = NA`).
#'
#' @param snps long data.frame with one row per SNP and region: `chr`, `pos`,
#'   `region`, `baf`.
#' @param cn optional long data.frame of copy number per region with
#'   `region`, `chr`, `seg_start`, `seg_end` and major / minor copy number
#'   (`nAraw`/`nBraw`, `nMajor`/`nMinor` or `maj_cn`/`min_cn`). Used to decide
#'   which segments are in allelic imbalance in each region. If `NULL`,
#'   allelic imbalance is called from the BAFs alone (mean |BAF - 0.5| >
#'   `min_imbalance`) on segments given by `segments`.
#' @param segments optional data.frame with `chr`, `seg_start`, `seg_end`
#'   used when `cn` is `NULL`. Defaults to one segment per chromosome.
#' @param min_imbalance minimum mean |BAF - 0.5| for a region to be
#'   considered imbalanced when `cn` is not given.
#' @param ai_tolerance minimum difference between major and minor copy number
#'   for allelic imbalance when `cn` is given.
#' @return `snps` with additional columns `segID` (common segment),
#'   `ref_region` (region used as reference), `haplotype` (multi-region
#'   haplotype of the SNP alternative allele) and `single_region_haplotype`
#'   (haplotype obtained from this region alone, `NA` if the region is not
#'   imbalanced at this SNP).
#' @export
phase_snps <- function(snps, cn = NULL, segments = NULL, min_imbalance = 0.1,
                       ai_tolerance = 0.1) {
  snps <- as.data.frame(snps)
  .require_cols(snps, c("chr", "pos", "baf"), "snps")
  if (!"region" %in% names(snps)) snps$region <- "region1"
  snps$chr <- add_chr_prefix(snps$chr)

  if (!is.null(cn)) {
    cn <- as.data.frame(cn)
    if (!"region" %in% names(cn)) cn$region <- snps$region[1]
    .require_cols(cn, c("chr", "seg_start", "seg_end"), "cn")
    cn$chr <- add_chr_prefix(cn$chr)
    maj <- intersect(c("nAraw", "nMajor", "maj_cn", "cnA"), names(cn))[1]
    min <- intersect(c("nBraw", "nMinor", "min_cn", "cnB"), names(cn))[1]
    if (is.na(maj) || is.na(min)) stop("cn needs major and minor copy number columns", call. = FALSE)
    common <- .disjoin_segments(unique(cn[, c("chr", "seg_start", "seg_end")]))
  } else {
    if (is.null(segments)) {
      rng <- stats::aggregate(pos ~ chr, data = snps, FUN = range)
      segments <- data.frame(chr = rng$chr, seg_start = rng$pos[, 1], seg_end = rng$pos[, 2])
    }
    common <- segments
    common$chr <- add_chr_prefix(common$chr)
  }
  common$segID <- seg_id(common$chr, common$seg_start, common$seg_end)
  sidx <- .match_segments(snps$chr, snps$pos, common$chr, common$seg_start, common$seg_end)
  snps$segID <- common$segID[sidx]

  # allelic imbalance strength per segment and region
  snps$dev <- abs(snps$baf - 0.5)
  ok <- !is.na(snps$segID) & !is.na(snps$baf)
  strength <- stats::aggregate(dev ~ segID + region, data = snps[ok, ], FUN = mean)
  if (!is.null(cn)) {
    mid <- floor((common$seg_start + common$seg_end) / 2)
    ai <- do.call(rbind, lapply(unique(cn$region), function(r) {
      cr <- cn[cn$region == r, ]
      k <- .match_segments(common$chr, mid, cr$chr, cr$seg_start, cr$seg_end)
      data.frame(segID = common$segID, region = r,
                 ai = !is.na(k) & abs(cr[[maj]][k] - cr[[min]][k]) > ai_tolerance)
    }))
    strength <- merge(strength, ai, by = c("segID", "region"), all.x = TRUE)
    strength$ai[is.na(strength$ai)] <- FALSE
  } else {
    strength$ai <- strength$dev > min_imbalance
  }
  imb <- strength[strength$ai, , drop = FALSE]
  imb <- imb[order(imb$segID, -imb$dev), , drop = FALSE]
  ref <- imb[!duplicated(imb$segID), c("segID", "region"), drop = FALSE]
  names(ref) <- c("segID", "ref_region")

  snps <- merge(snps, ref, by = "segID", all.x = TRUE, sort = FALSE)
  # haplotype from the reference region
  ref_baf <- snps[!is.na(snps$ref_region) & snps$region == snps$ref_region,
                  c("chr", "pos", "baf"), drop = FALSE]
  ref_baf$haplotype <- ifelse(ref_baf$baf >= 0.5, "A", "B")
  key <- paste(snps$chr, snps$pos)
  snps$haplotype <- ref_baf$haplotype[match(key, paste(ref_baf$chr, ref_baf$pos))]
  # haplotype from this region only (Figure 6-5, top)
  reg_ai <- paste(strength$segID, strength$region)[strength$ai]
  this_ai <- paste(snps$segID, snps$region) %in% reg_ai
  snps$single_region_haplotype <- ifelse(this_ai, ifelse(snps$baf >= 0.5, "A", "B"), NA)
  snps$dev <- NULL
  snps <- snps[order(snps$region, match(snps$chr, names(hg38_chr_sizes)), snps$pos), ]
  rownames(snps) <- NULL
  snps
}

#' Phase somatic mutations to haplotypes
#'
#' Step 3 of mutation phasing (thesis section 6.3.3, Figure 6-6). For each
#' SNP that is on the same long reads as a mutation, the fraction of mutated
#' reads carrying the SNP alternative allele is calculated. Using the
#' haplotype of the SNP alternative allele this is converted into the
#' fraction of mutated reads on haplotype A (one minus the fraction for SNPs
#' whose alternative allele is on haplotype B). The mean over all linked SNPs
#' decides the haplotype of the mutation: A if >= 0.5, B otherwise.
#'
#' Mutations whose mutated reads carry both alleles of a linked SNP are
#' flagged as possible parallel mutations (`parallel_mutation`).
#'
#' @param linkage data.frame with one row per mutation-SNP pair: `chr`,
#'   `pos` (mutation position), `snp_pos`, `mt_snp_alt` (number of reads with
#'   the mutation and the SNP alternative allele), `mt_snp_ref` (reads with
#'   the mutation and the SNP reference allele). Optional `wt_snp_alt`,
#'   `wt_snp_ref` (reads without the mutation) are kept.
#' @param snp_haplotypes phased SNPs, e.g. from [phase_snps()], with `chr`,
#'   `pos`, `haplotype`.
#' @param mutations optional mutation table (e.g. from
#'   [annotate_mutations()]). If given, phasing results are added to it and
#'   mutations that could not be phased get `phase = NA`.
#' @param min_mut_reads minimum number of mutated reads covering a SNP for the
#'   pair to be used.
#' @return If `mutations` is `NULL`, a data.frame with one row per phased
#'   mutation (`chr`, `pos`, `frac_hapA`, `n_linked_snps`,
#'   `n_linked_reads`, `parallel_mutation`, `phase`); otherwise `mutations`
#'   with these columns added.
#' @export
phase_mutations <- function(linkage, snp_haplotypes, mutations = NULL, min_mut_reads = 1) {
  linkage <- as.data.frame(linkage)
  .require_cols(linkage, c("chr", "pos", "snp_pos", "mt_snp_alt", "mt_snp_ref"), "linkage")
  .require_cols(snp_haplotypes, c("chr", "pos", "haplotype"), "snp_haplotypes")
  linkage$chr <- add_chr_prefix(linkage$chr)
  linkage$n_mut_reads <- linkage$mt_snp_alt + linkage$mt_snp_ref
  linkage <- linkage[linkage$n_mut_reads >= min_mut_reads, , drop = FALSE]

  # parallel mutations: mutated reads with both SNP alleles
  mut_key <- paste(linkage$chr, linkage$pos)
  par <- tapply(seq_len(nrow(linkage)), mut_key, function(i) {
    sum(linkage$mt_snp_alt[i]) > 0 && sum(linkage$mt_snp_ref[i]) > 0
  })

  hap <- unique(as.data.frame(snp_haplotypes)[, c("chr", "pos", "haplotype")])
  hap$chr <- add_chr_prefix(hap$chr)
  hap <- hap[!is.na(hap$haplotype), , drop = FALSE]
  hap <- hap[!duplicated(paste(hap$chr, hap$pos)), , drop = FALSE]
  linkage$snp_haplotype <- hap$haplotype[match(paste(linkage$chr, linkage$snp_pos),
                                               paste(hap$chr, hap$pos))]
  linkage <- linkage[!is.na(linkage$snp_haplotype), , drop = FALSE]

  linkage$snp_vaf_on_mut <- linkage$mt_snp_alt / linkage$n_mut_reads
  linkage$frac_hapA <- ifelse(linkage$snp_haplotype == "A",
                              linkage$snp_vaf_on_mut, 1 - linkage$snp_vaf_on_mut)
  key <- paste(linkage$chr, linkage$pos)
  res <- data.frame(
    key = unique(key),
    stringsAsFactors = FALSE
  )
  res$frac_hapA <- as.numeric(tapply(linkage$frac_hapA, key, mean)[res$key])
  res$n_linked_snps <- as.integer(tapply(linkage$frac_hapA, key, length)[res$key])
  res$n_linked_reads <- as.numeric(tapply(linkage$n_mut_reads, key, sum)[res$key])
  res$parallel_mutation <- as.logical(par[res$key])
  res$phase <- ifelse(res$frac_hapA >= 0.5, "A", "B")
  res$chr <- sub(" .*$", "", res$key)
  res$pos <- as.numeric(sub("^.* ", "", res$key))
  res <- res[, c("chr", "pos", "frac_hapA", "n_linked_snps", "n_linked_reads",
                 "parallel_mutation", "phase")]

  if (is.null(mutations)) return(res)
  mutations <- as.data.frame(mutations)
  mutations$chr <- add_chr_prefix(mutations$chr)
  mutations <- mutations[, setdiff(names(mutations), setdiff(names(res), c("chr", "pos"))), drop = FALSE]
  m <- match(paste(mutations$chr, mutations$pos), paste(res$chr, res$pos))
  for (col in setdiff(names(res), c("chr", "pos"))) mutations[[col]] <- res[[col]][m]
  mutations
}

#' Summarise phasing performance
#'
#' Computes, per sample, the percentage of mutations phased to haplotypes and
#' the percentage of incorrectly phased mutations estimated from regions of
#' loss of heterozygosity (LOH): a mutation phased to a haplotype with copy
#' number 0 must be phased incorrectly (thesis Figure 6-7).
#'
#' @param phased output of [phase_mutations()] with `mutations` given, i.e.
#'   containing `phase`, `cnA`, `cnB` and optionally `sample`.
#' @return data.frame with one row per sample: `n_mutations`, `n_phased`,
#'   `pct_phased`, `n_loh_phased`, `n_incorrect`, `pct_incorrect` (of all
#'   mutations) and `pct_incorrect_loh` (of phased mutations in LOH regions).
#' @export
phasing_summary <- function(phased) {
  phased <- as.data.frame(phased)
  if (!"sample" %in% names(phased)) phased$sample <- "sample"
  do.call(rbind, lapply(split(phased, phased$sample), function(d) {
    ph <- !is.na(d$phase)
    loh <- (d$cnA < 0.5 | d$cnB < 0.5) & ph
    wrong <- loh & ((d$phase == "A" & d$cnA < 0.5) | (d$phase == "B" & d$cnB < 0.5))
    data.frame(sample = d$sample[1], n_mutations = nrow(d), n_phased = sum(ph),
               pct_phased = 100 * mean(ph), n_loh_phased = sum(loh),
               n_incorrect = sum(wrong), pct_incorrect = 100 * sum(wrong) / nrow(d),
               pct_incorrect_loh = ifelse(sum(loh) > 0, 100 * sum(wrong) / sum(loh), NA))
  }))
}
