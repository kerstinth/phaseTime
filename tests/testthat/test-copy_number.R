test_that("mutation multiplicity follows the thesis formula", {
  # clonal heterozygous mutation on a diploid segment in a pure tumour
  expect_equal(mutation_multiplicity(0.5, purity = 1, total_cn = 2), 1)
  # n_mut = VAF / p * (p * CNt + CNn * (1 - p))
  expect_equal(mutation_multiplicity(0.3, 0.6, 3), 0.3 / 0.6 * (0.6 * 3 + 2 * 0.4))
})

test_that("Battenberg subclonal copy number is averaged", {
  bb <- data.frame(chr = c(1, 2), startpos = c(1, 1), endpos = c(100, 100),
                   nMaj1_A = c(2, 3), nMin1_A = c(1, 1), frac1_A = c(1, 0.6),
                   nMaj2_A = c(NA, 1), nMin2_A = c(NA, 1), frac2_A = c(NA, 0.4))
  out <- battenberg_average_cn(bb)
  expect_equal(out$chr, c("chr1", "chr2"))
  expect_equal(out$nAraw, c(2, 3 * 0.6 + 1 * 0.4))
  expect_equal(out$nBraw, c(1, 1))
})

test_that("haplotype copy number uses SNP haplotypes to flip major/minor", {
  cn <- data.frame(chr = "chr1", seg_start = c(1, 1001), seg_end = c(1000, 2000),
                   nMajor = c(2, 3), nMinor = c(1, 1))
  snps <- data.frame(chr = "chr1", pos = c(10, 20, 1100, 1200), haplotype = c("A", "B", "A", "B"),
                     baf = c(0.7, 0.3, 0.25, 0.75))
  out <- haplotype_copy_number(cn, snps = snps)
  # segment 1: haplotype A is the major allele
  expect_equal(c(out$cnA[1], out$cnB[1]), c(2, 1))
  # segment 2: haplotype A SNP alleles have low BAF, so A is the minor allele
  expect_equal(c(out$cnA[2], out$cnB[2]), c(1, 3))
  expect_true(all(out$allelic_imbalance))
})

test_that("haplotype copy number accepts an explicit flip table", {
  cn <- data.frame(chr = 1, seg_start = 1, seg_end = 100, nAraw = 3, nBraw = 1)
  out <- haplotype_copy_number(cn, flip = data.frame(chr = "chr1", seg_start = 1, seg_end = 100, flip = TRUE))
  expect_equal(c(out$cnA, out$cnB), c(1, 3))
})

test_that("MSAI table is parsed into region specific flips", {
  msai <- data.frame(seqnames = 1, start = 1, end = 100, msai = TRUE,
                     frac_consensus = "R1--0123456789=0.2;R2--0123456789=0.9")
  out <- parse_battenberg_msai(msai)
  expect_equal(out$region, c("R1", "R2"))
  expect_equal(out$flip, c(TRUE, FALSE))
})

test_that("mutations are annotated with segment copy number and multiplicity", {
  cn <- haplotype_copy_number(data.frame(chr = "chr1", seg_start = c(1, 1001), seg_end = c(1000, 2000),
                                         nMajor = c(2, 1), nMinor = c(1, 1)))
  muts <- data.frame(chr = c(1, 1, 1), pos = c(500, 1500, 5000), vaf = c(0.4, 0.25, 0.3))
  out <- annotate_mutations(muts, cn, purity = 0.8)
  expect_equal(nrow(out), 2)  # the third mutation is outside all segments
  expect_equal(out$tot_cn, c(3, 2))
  expect_equal(out$mut_cpn, mutation_multiplicity(c(0.4, 0.25), 0.8, c(3, 2)))
})

test_that("whole genome doubling is detected", {
  cn <- data.frame(chr = paste0("chr", 1:3), seg_start = 1, seg_end = 100, cnA = c(2, 2, 1), cnB = c(2, 1, 1))
  expect_true(detect_wgd(cn)$first_gd)
  expect_false(detect_wgd(cn)$second_gd)
})
