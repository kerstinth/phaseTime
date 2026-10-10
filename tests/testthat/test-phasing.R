test_that("closest SNPs within the maximum distance are selected", {
  muts <- data.frame(chr = "chr1", pos = c(1000, 50000))
  snps <- data.frame(chr = 1, pos = c(500, 1500, 1600, 30000, 49000, 90000))
  out <- select_phasing_snps(muts, snps, max_snps = 2, max_distance = 20000)
  expect_equal(out$snp_pos[out$pos == 1000], c(500, 1500))
  expect_equal(out$snp_pos[out$pos == 50000], c(49000, 30000))
})

test_that("multi-region reference phasing uses the region with the strongest imbalance", {
  # region 1 balanced, region 2 imbalanced on the same segment
  snps <- data.frame(region = rep(c("R1", "R2"), each = 4), chr = "chr1",
                     pos = rep(c(10, 20, 30, 40), 2),
                     baf = c(0.5, 0.52, 0.48, 0.5, 0.75, 0.25, 0.8, 0.2))
  cn <- data.frame(region = c("R1", "R2"), chr = "chr1", seg_start = 1, seg_end = 100,
                   nMajor = c(1, 2), nMinor = c(1, 1))
  out <- phase_snps(snps, cn = cn)
  r1 <- out[out$region == "R1", ]
  expect_true(all(r1$ref_region == "R2"))
  expect_equal(r1$haplotype, c("A", "B", "A", "B"))
  expect_true(all(is.na(r1$single_region_haplotype)))
  expect_equal(out$single_region_haplotype[out$region == "R2"], c("A", "B", "A", "B"))
})

test_that("mutations are phased by the mean fraction of reads on haplotype A", {
  snp_hap <- data.frame(chr = "chr1", pos = c(100, 200, 300), haplotype = c("A", "B", "A"))
  linkage <- data.frame(chr = "chr1",
                        pos = c(150, 150, 250, 250),
                        snp_pos = c(100, 200, 200, 300),
                        mt_snp_alt = c(5, 0, 4, 1),
                        mt_snp_ref = c(0, 5, 0, 3))
  out <- phase_mutations(linkage, snp_hap)
  # mutation 150: alt of A-SNP (1) and ref of B-SNP (1 - 0) -> A
  expect_equal(out$phase[out$pos == 150], "A")
  expect_equal(out$frac_hapA[out$pos == 150], 1)
  # mutation 250: alt of B-SNP (frac A = 0) and 1/4 alt of A-SNP -> mean 0.125 -> B
  expect_equal(out$phase[out$pos == 250], "B")
  expect_equal(out$frac_hapA[out$pos == 250], 0.125)
  expect_true(out$parallel_mutation[out$pos == 250])
})

test_that("phasing summary estimates errors from LOH regions", {
  d <- data.frame(phase = c("A", "B", NA, "B"), cnA = c(2, 2, 1, 1), cnB = c(0, 0, 1, 1))
  s <- phasing_summary(d)
  expect_equal(s$n_phased, 3)
  expect_equal(s$n_incorrect, 1)
  expect_equal(s$pct_incorrect_loh, 50)
})
