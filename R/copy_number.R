# copy number processing: Battenberg input, haplotype specific copy number,
# mutation multiplicity

#' Mutation multiplicity (mutation copy number)
#'
#' Computes the number of chromosomal copies carrying a mutation (n_mut) from
#' its variant allele frequency (thesis section 6.3.5; McGranahan et al. 2015):
#'
#' \deqn{n_{mut} = VAF \frac{1}{p} [p CN_t + CN_n (1 - p)]}
#'
#' @param vaf variant allele frequency.
#' @param purity tumour purity (p).
#' @param total_cn tumour locus-specific total copy number (CN_t).
#' @param normal_cn normal locus-specific copy number (CN_n), default 2.
#' @return numeric vector of mutation multiplicities.
#' @export
#' @examples
#' mutation_multiplicity(vaf = 0.35, purity = 0.7, total_cn = 2)
mutation_multiplicity <- function(vaf, purity, total_cn, normal_cn = 2) {
  (vaf / purity) * (purity * total_cn + normal_cn * (1 - purity))
}

#' Collapse Battenberg subclonal copy number to an average copy number
#'
#' Battenberg reports up to two copy number states per segment
#' (`nMaj1_A`, `nMin1_A`, `frac1_A`, `nMaj2_A`, `nMin2_A`, `frac2_A`). As in the
#' thesis, the copy number of a segment is the average of the clone states
#' weighted by their cellular fractions. Clonal segments therefore have
#' integer copy numbers and subclonal segments non-integer copy numbers.
#'
#' @param subclones data.frame read from a Battenberg `*_subclones.txt` file.
#' @return data.frame with `chr` (with "chr" prefix), `seg_start`, `seg_end`,
#'   `nAraw`, `nBraw` (average major / minor copy number), `nMajor`, `nMinor`
#'   (rounded) and `cnTotal`.
#' @export
#' @examples
#' bb <- data.frame(chr = 1, startpos = 1, endpos = 1e6,
#'                  nMaj1_A = 2, nMin1_A = 1, frac1_A = 0.6,
#'                  nMaj2_A = 1, nMin2_A = 1, frac2_A = 0.4)
#' battenberg_average_cn(bb)
battenberg_average_cn <- function(subclones) {
  .require_cols(subclones, c("chr", "startpos", "endpos", "nMaj1_A", "nMin1_A", "frac1_A"),
                "subclones")
  sc <- as.data.frame(subclones)
  frac1 <- sc$frac1_A
  maj2 <- if ("nMaj2_A" %in% names(sc)) sc$nMaj2_A else NA
  min2 <- if ("nMin2_A" %in% names(sc)) sc$nMin2_A else NA
  frac2 <- if ("frac2_A" %in% names(sc)) sc$frac2_A else NA
  clonal <- is.na(frac2) | frac1 == 1
  nAraw <- ifelse(clonal, sc$nMaj1_A, sc$nMaj1_A * frac1 + maj2 * frac2)
  nBraw <- ifelse(clonal, sc$nMin1_A, sc$nMin1_A * frac1 + min2 * frac2)
  out <- data.frame(chr = add_chr_prefix(sc$chr), seg_start = sc$startpos, seg_end = sc$endpos,
                    nAraw = nAraw, nBraw = nBraw, stringsAsFactors = FALSE)
  out$nMajor <- round(out$nAraw)
  out$nMinor <- round(out$nBraw)
  out$cnTotal <- out$nMajor + out$nMinor
  out
}

#' Read the Battenberg multi-sample allelic imbalance (MSAI) table
#'
#' Converts the Battenberg `multisample_MSAI.txt` output into one row per
#' segment and tumour region, with a `flip` column indicating whether the
#' major allele of that region is haplotype B of the multi-region reference
#' haplotype (`frac_consensus` < 0.5), and whether the region shows allelic
#' imbalance. This mirrors the processing in the original phaseTime scripts.
#'
#' @param msai data.frame with columns `seqnames`, `start`, `end`, `msai` and
#'   `frac_consensus` (semicolon separated `region=value` pairs).
#' @return long data.frame with columns `chr`, `seg_start`, `seg_end`, `msai`,
#'   `region`, `frac_consensus`, `flip`, `allelic_imbalance`.
#' @export
parse_battenberg_msai <- function(msai) {
  .require_cols(msai, c("seqnames", "start", "end", "msai", "frac_consensus"), "msai")
  rows <- lapply(seq_len(nrow(msai)), function(i) {
    kv <- strsplit(as.character(msai$frac_consensus[i]), ";", fixed = TRUE)[[1]]
    kv <- kv[nzchar(kv)]
    if (length(kv) == 0) return(NULL)
    key <- sub("=.*$", "", kv)
    val <- suppressWarnings(as.numeric(sub("^.*=", "", kv)))
    # Battenberg appends a 10 character suffix to region names; remove it
    key <- sub("--..........$", "", key)
    data.frame(chr = add_chr_prefix(msai$seqnames[i]), seg_start = msai$start[i],
               seg_end = msai$end[i], msai = msai$msai[i], region = key,
               frac_consensus = val, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$flip <- !is.na(out$frac_consensus) & out$frac_consensus < 0.5
  out$allelic_imbalance <- !is.na(out$frac_consensus)
  out
}

#' Assign copy number to haplotypes A and B
#'
#' Copy number callers report a major and a minor allele copy number per
#' segment. To compare copy number between haplotypes (and between tumour
#' regions) the major/minor copy number has to be assigned to the haplotypes
#' A and B defined by SNP phasing (thesis section 6.3.2). This can be done in
#' two ways:
#'
#' * `snps`: phased SNPs (output of [phase_snps()]) with their B-allele
#'   frequency in this tumour region. In each segment the mean fraction of
#'   reads supporting haplotype A is computed; if it is below 0.5, haplotype A
#'   is the minor allele in this region (`flip = TRUE`).
#' * `flip`: a data.frame with `chr`, `seg_start`, `seg_end` and a logical
#'   `flip` column (e.g. from [parse_battenberg_msai()]) giving the assignment
#'   directly.
#'
#' Segments without any information keep the major allele as haplotype A.
#'
#' @param cn data.frame with `chr`, `seg_start`, `seg_end` and major / minor
#'   copy number columns (see `major_col`, `minor_col`).
#' @param snps optional phased SNP table with `chr`, `pos`, `haplotype` (the
#'   haplotype carrying the SNP alternative allele, "A" or "B") and `baf`
#'   (B-allele frequency of the SNP in this region).
#' @param flip optional data.frame with columns `chr`, `seg_start`, `seg_end`,
#'   `flip`.
#' @param major_col,minor_col names of the major / minor copy number columns.
#'   By default the first of `nAraw`/`nMajor`/`maj_cn` (resp.
#'   `nBraw`/`nMinor`/`min_cn`) found is used.
#' @param ai_tolerance difference between `cnA` and `cnB` above which a
#'   segment is considered to be in allelic imbalance.
#' @return `cn` with additional columns `cnA`, `cnB`, `tot_cn`, `flip`,
#'   `allelic_imbalance` and `segID`.
#' @export
haplotype_copy_number <- function(cn, snps = NULL, flip = NULL,
                                  major_col = NULL, minor_col = NULL,
                                  ai_tolerance = 0.1) {
  cn <- as.data.frame(cn)
  .require_cols(cn, c("chr", "seg_start", "seg_end"), "cn")
  pick <- function(cands, given) {
    if (!is.null(given)) return(given)
    hit <- cands[cands %in% names(cn)]
    if (length(hit) == 0) stop("could not find copy number column, tried: ",
                               paste(cands, collapse = ", "), call. = FALSE)
    hit[1]
  }
  major_col <- pick(c("nAraw", "nMajor", "maj_cn"), major_col)
  minor_col <- pick(c("nBraw", "nMinor", "min_cn"), minor_col)
  cn$chr <- add_chr_prefix(cn$chr)
  cn$flip <- FALSE

  if (!is.null(flip)) {
    .require_cols(flip, c("chr", "seg_start", "seg_end", "flip"), "flip")
    flip$chr <- add_chr_prefix(flip$chr)
    key_cn <- paste(cn$chr, cn$seg_start, cn$seg_end)
    key_fl <- paste(flip$chr, flip$seg_start, flip$seg_end)
    m <- match(key_cn, key_fl)
    # segments that do not match exactly are matched by midpoint
    nm <- which(is.na(m))
    if (length(nm) > 0) {
      mid <- floor((cn$seg_start[nm] + cn$seg_end[nm]) / 2)
      m[nm] <- .match_segments(cn$chr[nm], mid, flip$chr, flip$seg_start, flip$seg_end)
    }
    cn$flip[!is.na(m)] <- as.logical(flip$flip[m[!is.na(m)]])
  } else if (!is.null(snps)) {
    .require_cols(snps, c("chr", "pos", "haplotype", "baf"), "snps")
    snps <- snps[!is.na(snps$haplotype) & !is.na(snps$baf), , drop = FALSE]
    snps$chr <- add_chr_prefix(snps$chr)
    frac_A <- ifelse(snps$haplotype == "A", snps$baf, 1 - snps$baf)
    seg_idx <- .match_segments(snps$chr, snps$pos, cn$chr, cn$seg_start, cn$seg_end)
    ok <- !is.na(seg_idx)
    if (any(ok)) {
      mean_frac <- tapply(frac_A[ok], seg_idx[ok], mean)
      idx <- as.integer(names(mean_frac))
      cn$frac_hapA <- NA_real_
      cn$frac_hapA[idx] <- as.numeric(mean_frac)
      cn$flip[idx] <- mean_frac < 0.5
    }
  }
  cn$flip[is.na(cn$flip)] <- FALSE
  cn$cnA <- ifelse(cn$flip, cn[[minor_col]], cn[[major_col]])
  cn$cnB <- ifelse(cn$flip, cn[[major_col]], cn[[minor_col]])
  cn$tot_cn <- cn$cnA + cn$cnB
  cn$allelic_imbalance <- abs(cn$cnA - cn$cnB) > ai_tolerance
  cn$segID <- seg_id(cn$chr, cn$seg_start, cn$seg_end)
  cn
}

#' Annotate mutations with copy number and multiplicity
#'
#' Overlaps mutations with haplotype specific copy number segments and
#' computes the mutation multiplicity ([mutation_multiplicity()]) of every
#' mutation from its VAF, the tumour purity and the total copy number of the
#' segment.
#'
#' @param mutations data.frame with `chr`, `pos` and either `vaf` or
#'   `var_count` and `depth`. Further columns (e.g. `cluster`, `ref`, `alt`)
#'   are kept.
#' @param cn haplotype specific copy number, output of
#'   [haplotype_copy_number()].
#' @param purity tumour purity of the region.
#' @param normal_cn normal copy number (2 for autosomes).
#' @return mutations on copy number segments with the added columns
#'   `seg_start`, `seg_end`, `segID`, `cnA`, `cnB`, `tot_cn`,
#'   `allelic_imbalance` and `mut_cpn`.
#' @export
annotate_mutations <- function(mutations, cn, purity, normal_cn = 2) {
  mutations <- as.data.frame(mutations)
  .require_cols(mutations, c("chr", "pos"), "mutations")
  .require_cols(cn, c("chr", "seg_start", "seg_end", "cnA", "cnB"), "cn")
  if (!"vaf" %in% names(mutations)) {
    .require_cols(mutations, c("var_count", "depth"), "mutations")
    mutations$vaf <- mutations$var_count / mutations$depth
  }
  mutations$chr <- add_chr_prefix(mutations$chr)
  idx <- .match_segments(mutations$chr, mutations$pos, cn$chr, cn$seg_start, cn$seg_end)
  mutations <- mutations[!is.na(idx), , drop = FALSE]
  idx <- idx[!is.na(idx)]
  drop <- intersect(names(mutations),
                    c("seg_start", "seg_end", "segID", "cnA", "cnB", "tot_cn", "allelic_imbalance", "mut_cpn"))
  mutations <- mutations[, setdiff(names(mutations), drop), drop = FALSE]
  mutations$seg_start <- cn$seg_start[idx]
  mutations$seg_end <- cn$seg_end[idx]
  mutations$segID <- seg_id(mutations$chr, mutations$seg_start, mutations$seg_end)
  mutations$cnA <- cn$cnA[idx]
  mutations$cnB <- cn$cnB[idx]
  mutations$tot_cn <- mutations$cnA + mutations$cnB
  mutations$allelic_imbalance <- if ("allelic_imbalance" %in% names(cn)) cn$allelic_imbalance[idx] else
    abs(mutations$cnA - mutations$cnB) > 0.1
  mutations$mut_cpn <- mutation_multiplicity(mutations$vaf, purity, mutations$tot_cn, normal_cn)
  rownames(mutations) <- NULL
  mutations
}

#' Detect whole genome doubling
#'
#' Determines whether one or two whole genome doubling events occurred, based
#' on the fraction of the autosomal genome with a major copy number of at
#' least 2 (resp. 3). Written by Alex Frankell and adapted for phaseTime.
#'
#' @param cn data.frame with `chr`, `seg_start`, `seg_end` and haplotype copy
#'   numbers `cnA`, `cnB` (or `maj_cn`, `min_cn`).
#' @return one-row data.frame with `first_gd`, `genome_frac_nMaj2`,
#'   `second_gd`, `genome_frac_nMaj3`, `genome_frac_nMaj4`,
#'   `genome_frac_balance`, `mean_allelic_difference`.
#' @export
detect_wgd <- function(cn) {
  cn <- as.data.frame(cn)
  a <- if ("cnA" %in% names(cn)) cn$cnA else cn$maj_cn
  b <- if ("cnB" %in% names(cn)) cn$cnB else cn$min_cn
  seg <- unique(data.frame(chr = add_chr_prefix(cn$chr), seg_start = cn$seg_start,
                           seg_end = cn$seg_end, maj = pmax(a, b), min = pmin(a, b)))
  autosomal <- seg$chr %in% paste0("chr", 1:22)
  w <- seg$seg_end - seg$seg_start
  frac <- w / sum(w)
  data.frame(
    first_gd = sum(frac[autosomal & seg$maj >= 2]) > 0.5,
    genome_frac_nMaj2 = sum(frac[autosomal & seg$maj >= 2]),
    second_gd = sum(frac[autosomal & seg$maj >= 3]) > 0.5,
    genome_frac_nMaj3 = sum(frac[autosomal & seg$maj >= 3]),
    genome_frac_nMaj4 = sum(frac[autosomal & seg$maj >= 4]),
    genome_frac_balance = sum(frac[seg$maj == seg$min]),
    mean_allelic_difference = sum((seg$maj - seg$min) * frac)
  )
}
