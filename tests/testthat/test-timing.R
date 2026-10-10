test_that("gain timing reproduces the worked example of Figure 6-13", {
  tm <- gain_timing_from_counts(states = c(1, 2, 3), counts = c(3, 2, 1))
  expect_equal(tm$gain_timing, c(1 / 3, 2 / 3))
  expect_equal(tm$copy_number_before_gain, c(1, 2))
  expect_equal(tm$copy_number_after_gain, c(2, 3))
})

test_that("a single gain is timed by early / all copy number adjusted mutations", {
  # 20 mutations before a 1 -> 2 gain (multiplicity 2), 60 after on 2 copies
  tm <- gain_timing_from_counts(c(1, 2), c(60, 20))
  expect_equal(tm$gain_timing, 20 / (20 + 60 / 2))
  # legacy adjustment divides by multiplicity
  tl <- gain_timing_from_counts(c(1, 2), c(60, 20), adjust = "multiplicity")
  expect_equal(tl$gain_timing, (20 / 2) / (20 / 2 + 60))
  expect_true(tm$loCI < tm$gain_timing && tm$hiCI > tm$gain_timing)
})

test_that("dropped boundaries pool mutations into one interval", {
  tm <- gain_timing_from_counts(c(1, 2, 3), c(3, 2, 1), keep = c(TRUE, FALSE))
  expect_equal(nrow(tm), 1)
  # states 2 and 1 are pooled and adjusted by the copies present after the gain (2)
  expect_equal(tm$gain_timing, 1 / (1 + (2 + 3) / 2))
  expect_equal(nrow(gain_timing_from_counts(3, 10)), 0)
})

make_muts <- function(states, counts, cluster = "1", cnA = max(states), seg = 1) {
  data.frame(chr = "chr1", pos = seq_len(sum(counts)) + seg * 1e6,
             seg_start = seg * 1e6, seg_end = seg * 1e6 + 1e5,
             segID = seg_id("chr1", seg * 1e6, seg * 1e6 + 1e5),
             phase = "A", cnA = cnA, cnB = 1,
             mut_cpn = rep(states, counts) + rep(c(-0.1, 0.1), length.out = sum(counts)),
             cluster = cluster)
}

test_that("mutations are classified as early or late", {
  m <- classify_mutation_timing(make_muts(c(1, 2), c(10, 5)))
  expect_equal(as.character(m$timing), rep(c("late", "early"), c(10, 5)))
  m$cnA <- 2.5
  expect_true(all(classify_mutation_timing(m)$timing == "not_possible"))
})

test_that("copy number states with fewer than min_muts mutations are removed", {
  m <- assign_cn_states(make_muts(c(1, 2, 3), c(10, 3, 10)))
  expect_equal(sort(unique(stats::na.omit(m$cn_state))), c(1, 3))
})

test_that("time_gains times clonal gains and excludes gains of descendant clusters", {
  tree <- list(edges = cbind(c("1"), c("2")), trunk = "1")
  # trunk: gain 1 -> 2 (states 4 and 2), cluster 2: doubling 2 -> 4 (states 2 and 1)
  m <- rbind(make_muts(c(2, 4), c(40, 20), "1", cnA = 4),
             make_muts(c(1, 2), c(40, 10), "2", cnA = 4))
  m$pos <- seq_len(nrow(m))
  clonal <- time_gains(m, cluster_col = "cluster", tree = tree)
  expect_equal(nrow(clonal), 1)
  expect_equal(clonal$cluster, "1")
  expect_equal(c(clonal$copy_number_before_gain, clonal$copy_number_after_gain), c(1, 2))
  # copies: state 4 -> 1 copy, state 2 -> 2 copies
  expect_equal(clonal$gain_timing, 20 / (20 + 40 / 2))
  all <- time_gains(m, cluster_col = "cluster", tree = tree, clonal_only = FALSE)
  sub <- all[all$cluster == "2", ]
  expect_equal(c(sub$copy_number_before_gain, sub$copy_number_after_gain), c(2, 4))
  expect_equal(sub$gain_timing, (10 / 2) / (10 / 2 + 40 / 4))
})

test_that("untimed segments are reported with a reason", {
  m <- make_muts(c(1), c(10), cnA = 1)
  out <- time_gains(m, include_untimed = TRUE)
  expect_true(all(out$status %in% c("no_gain")))
})

test_that("gain timing is compared between regions", {
  t1 <- data.frame(sample = c("T1_R1", "T1_R2"), chr = "chr1", seg_start = 1, seg_end = c(100, 95),
                   segID = c("chr1_1_100", "chr1_1_95"), haplotype = "A", gain_index = 1,
                   gain_timing = c(0.3, 0.35), se = c(0.05, 0.05))
  cmp <- compare_region_timing(t1)
  expect_equal(nrow(cmp), 1)
  expect_equal(cmp$difference, 0.05)
  expect_equal(cmp$p_value, 2 * pnorm(-0.05 / sqrt(2 * 0.05^2)))
})

test_that("CONIPHER style trees are converted", {
  tr <- as_phaseTime_tree(list(graph_pyclone = list(default_tree = matrix(c(1, 2, 2, 3), 2, byrow = TRUE),
                                                    trunk = 1)))
  expect_equal(tr$trunk, "1")
  expect_equal(unname(tr$edges[, 2]), c("2", "3"))
})
