test_that("simulated gain timing recovers the ground truth at high mutation rates", {
  set.seed(1)
  sim <- simulate_gain_history(c(3000, 10000), c(1, 2, 2), segment_size = 5e7)
  tm <- simulated_gain_timing(sim)
  expect_equal(tm$theoretical_gain_timing, 0.3)
  expect_lt(abs(tm$gain_timing - 0.3), 0.02)
})

test_that("phasing simulation increases with read length", {
  s <- simulate_snp_phasing(c(150, 5000), n_mutations = 5000, seed = 1)
  expect_lt(s$pct_phased[1], 20)
  expect_gt(s$pct_phased[2], 90)
})

test_that("detection limit is computed analytically", {
  d <- detection_limit_cn(0.5, min_vaf = 0.05, normal_cn = 2)
  # VAF of a single copy mutation at the limit is >= 5%, one copy more is below
  vaf <- function(cn) 0.5 / (0.5 * cn + 2 * 0.5)
  expect_gte(vaf(d$max_copy_number), 0.05)
  expect_lt(vaf(d$max_copy_number + 1), 0.05)
})

test_that("the full workflow recovers simulated gain timing", {
  sim <- simulate_tumour(seed = 2)
  res <- run_phaseTime(sim$mutations, sim$cn, sim$purity, sim$snps, sim$linkage,
                       cluster_col = "cluster", tree = sim$tree)
  expect_true(all(c("snps", "cn", "mutations", "timing", "phasing") %in% names(res)))
  # phasing is accurate
  m <- res$mutations[!is.na(res$mutations$phase), ]
  expect_gt(mean(m$phase == m$true_phase), 0.97)
  # first clonal gains agree with the truth
  tm <- res$timing[res$timing$gain_index == 1, ]
  tr <- sim$truth[sim$truth$cluster == "1" & sim$truth$copy_number_before_gain == 1, ]
  tm$truth <- tr$true_timing[match(paste(tm$segID, tm$haplotype), paste(tr$segID, tr$haplotype))]
  tm <- tm[!is.na(tm$truth) & tm$n_mutations >= 100, ]
  expect_gt(nrow(tm), 10)
  expect_gt(cor(tm$gain_timing, tm$truth), 0.9)
  expect_lt(stats::median(abs(tm$gain_timing - tm$truth)), 0.05)
})
