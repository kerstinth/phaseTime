# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #    functions to analyse allele specific gain timing across regions    # # #  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# author: Kerstin Thol

################################################################################
#####################################################################  libraries

library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(cowplot)
library(BSgenome.Hsapiens.UCSC.hg38)
library(RColorBrewer)
library(GenomicRanges)

################################################################################
##################################################################### file paths

################################################################################
####################################################################   functions

# function to make mutations data frame continuous across several chromosomes
genome.continous.fun <- function(data, genome){
  
  data <- data[data$chr %in% paste0("chr", c(1:22)), ]
  chr_sizes <- seqlengths(genome)[paste0("chr", 1:22)]
  chr_sizes <- data.frame(chr_sizes)
  chr_sizes$chr <- c(1:22)
  chr_sizes <- chr_sizes[order(chr_sizes$chr, decreasing = F), ]
  chr_sizes <- chr_sizes$chr_sizes
  
  # adjust pos to continous value
  cumsum_chr_size <- cumsum(as.numeric(chr_sizes))
  
  data$adjust_pos <- data$pos
  data$chr_adjust <- sub("chr", "", data$chr)
  data$chr_adjust <- as.numeric(as.character(data$chr_adjust))
  
  for(i in c(2:22)){
    data$adjust_pos[data$chr_adjust == i] <- data$adjust_pos[data$chr_adjust == i] + cumsum_chr_size[i-1]
  }
  
  # also adjust seg start
  data$adjust_seg_start <- data$seg_start
  for(i in c(2:22)){
    data$adjust_seg_start[data$chr_adjust == i] <- data$adjust_seg_start[data$chr_adjust == i] + cumsum_chr_size[i-1]
  }
  
  # and seg end
  data$adjust_seg_end <- data$seg_end
  for(i in c(2:22)){
    data$adjust_seg_end[data$chr_adjust == i] <- data$adjust_seg_end[data$chr_adjust == i] + cumsum_chr_size[i-1]
  }
  
  return(data)
}

# function to get colors
get.colors <- function(max_cn = max(mutation_df$mut_cpn_round), opacity = 255){
  cols <- paste(brewer.pal(max_cn, name = "Dark2"), sep = "")
  nopaq <- round(opacity / ceiling(max_cn / 8))
  opacity_val <- c(sapply(1:ceiling(max_cn / 8), function(x) rep(opacity - nopaq * (x - 1), 8)))[1:max_cn]
  cols.opac <- paste(cols, as.hexmode(opacity_val), sep = "")
}


# function to compare timing estimates between two regions of the same tumour

plot.gaintiming.regions.fun <- function(timed_data = all_timed_data, regions, output_path = input_path, genome = BSgenome.Hsapiens.UCSC.hg38){
  
  print(paste0("plotting haplotype specific gain timing of regions: ", paste(regions, collapse = " & ")))
  
  # get only segment data so we don't have this duplicated for each mutation
  timed_data <- unique(timed_data[, c("sample", "segID", "chr", "seg_start", "seg_end",
                               "maj_cn", "min_cn", "tot_cn", "phase", "early", "late", "total_mut", "seg_timing", "seg_length",
                               "adjustment_factor", "late_adjusted", "seg_timing_adjusted", "adjusted_total_muts_on_segment",
                               "error", "loCI", "hiCI")])
  
  # subset for regions of interest if that wasnt already done before running this function
  timed_data_region <- timed_data[timed_data$sample %in% regions, ]
  
  timed_data_region$phase_info <- paste0(timed_data_region$early, "_", timed_data_region$late, "_", timed_data_region$total_mut, 
                                         "_", timed_data_region$seg_timing, "_", timed_data_region$seg_length, "_", timed_data_region$adjustment_factor,
                                         "_", timed_data_region$late_adjusted, "_", timed_data_region$seg_timing_adjusted, "_",
                                         timed_data_region$adjusted_total_muts_on_segment, "_", timed_data_region$error, "_", timed_data_region$loCI,
                                         "_", timed_data_region$hiCI)
  
  timing_mat <- reshape2::dcast(timed_data_region, sample + segID + chr + seg_start + seg_end + maj_cn + min_cn + tot_cn ~ phase, value.var = "phase_info")
  timing_mat <- timing_mat %>% separate(A, into = c("early_A", "late_A", "total_mut_A", "seg_timing_A", "seg_length_A",
                                                    "adjustment_factor_A", "late_adjusted_A", "adjusted_seg_timing_A", "adjusted_total_muts_on_segment_A",
                                                    "error_A", "loCI_A", "hiCI_A"), sep = "_")
  
  timing_mat <- timing_mat %>% separate(B, into = c("early_B", "late_B", "total_mut_B", "seg_timing_B", "seg_length_B",
                                                    "adjustment_factor_B", "late_Bdjusted_B", "adjusted_seg_timing_B", "adjusted_total_muts_on_segment_B",
                                                    "error_B", "loCI_B", "hiCI_B"), sep = "_")
  
  timing_mat$adjusted_seg_timing_A <- as.numeric(as.character(timing_mat$adjusted_seg_timing_A))
  timing_mat$adjusted_seg_timing_B <- as.numeric(as.character(timing_mat$adjusted_seg_timing_B))
  timing_mat$loCI_A <- as.numeric(as.character(timing_mat$loCI_A))
  timing_mat$loCI_B <- as.numeric(as.character(timing_mat$loCI_B))
  timing_mat$hiCI_A <- as.numeric(as.character(timing_mat$hiCI_A))
  timing_mat$hiCI_B <- as.numeric(as.character(timing_mat$hiCI_B))
  
  timing_mat[is.na(timing_mat$adjusted_seg_timing_A), "timing_possibleA"] <- "no"
  timing_mat[is.na(timing_mat$adjusted_seg_timing_B), "timing_possibleB"] <- "no"
  timing_mat[!is.na(timing_mat$adjusted_seg_timing_A), "timing_possibleA"] <- "yes"
  timing_mat[!is.na(timing_mat$adjusted_seg_timing_B), "timing_possibleB"] <- "yes"
  timing_mat[is.na(timing_mat$adjusted_seg_timing_A), "adjusted_seg_timing_A"] <- 0
  timing_mat[is.na(timing_mat$adjusted_seg_timing_B), "adjusted_seg_timing_B"] <- 0
  
  timing_mat$timing_possibleA <- factor(timing_mat$timing_possibleA, levels = c("yes", "no"))
  timing_mat$timing_possibleB <- factor(timing_mat$timing_possibleB, levels = c("yes", "no"))
  
  # first let's plot this for both regions
  # let's try to write this code in a way that it could take more than two regions into account even though we don't have that right now
  # first plot the timing per haplotype for all regions within one plot, leave our CIs for this first plot
  # we shall calculate if there is an actual difference between regions later
  
  colors <- get.colors(length(regions))
  
  # make a dummy haplotype collumn for plotting purposes
  timing_mat$haplotype <- "haplotype"
  
  # make the genome continuos
  timing_mat <- genome.continous.fun(timing_mat, genome)

  chr_sizes <- seqlengths(genome)[paste0("chr", 1:22)]
  chr_lines <- cumsum(as.numeric(as.character(chr_sizes)))
  chr_midpoint <- chr_lines - chr_sizes / 2
  xaxis_label <- data.frame(chr_midpoint)
  xaxis_label$label <- rownames(xaxis_label)
  
  p_genome1 <- ggplot() +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    geom_segment(data = timing_mat,
                 aes(x = adjust_seg_start, xend = adjust_seg_end, y = adjusted_seg_timing_A, yend = adjusted_seg_timing_A, color = sample, alpha = timing_possibleA),
                 linewidth = 2) +
    scale_color_manual(name = "tumour region",
                       values = c(colors)) +
    scale_alpha_manual(name = "", breaks = c("yes", "no"),
                       values = c("yes" = 1,
                                  "no" = 0.3),
                       labels = c("gain timing possible", "gain timing not possible")) +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    facet_wrap(~haplotype, labeller = labeller(haplotype = c(haplotype = "haplotype A")), strip.position = "right") +
    ylim(c(0, 1)) +
    ylab("gain timing") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm"))
  
  p_genome2 <- ggplot() +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    geom_segment(data = timing_mat,
                 aes(x = adjust_seg_start, xend = adjust_seg_end, y = adjusted_seg_timing_B, yend = adjusted_seg_timing_B, color = sample, alpha = timing_possibleB),
                 linewidth = 2) +
    scale_color_manual(name = "tumour region",
                       values = c(colors)) +
    scale_alpha_manual(name = "", breaks = c("yes", "no"),
                       values = c("yes" = 1,
                                  "no" = 0.3),
                       labels = c("gain timing possible", "gain timing not possible")) +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    facet_wrap(~haplotype, labeller = labeller(haplotype = c(haplotype = "haplotype B")), strip.position = "right") +
    ylim(c(0, 1)) +
    ylab("gain timing") +
    xlab("") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          plot.margin = unit(c(0, 0.1, 0.1, 0.1), "cm"))
  
  # combine the two plots
  p_genome_c <- plot_grid(p_genome1, p_genome2, ncol = 1, align = "v", axis = "rl", rel_heights = c(1, 1.1))
  
  if(!dir.exists(paste0(output_path, tumour, "/", "/plots/comparing_regions_phased_mutation_timing"))){
    dir.create(paste0(output_path, tumour, "/", "/plots/comparing_regions_phased_mutation_timing"), recursive = T)
  }
  
  pdf(paste0(output_path, tumour, "/plots/comparing_regions_phased_mutation_timing/", tumour, "_genome_phased_gain_timing_per_haplotype_regions.pdf"), width = 18, height = 7)
  print(p_genome_c)
  dev.off()
  
}

compare.gaintiming.regions.fun <- function(timed_data = all_timed_data, tumour, regions, output_path = input_path, genome = BSgenome.Hsapiens.UCSC.hg38){

  print(paste0("comparing haplotype specific gain timing of regions: ", paste(regions, collapse = " & ")))
  
  # get only segment data so we don't have this duplicated for each mutation
  timed_data <- unique(timed_data[, c("sample", "segID", "chr", "seg_start", "seg_end",
                                      "maj_cn", "min_cn", "tot_cn", "phase", "early", "late", "total_mut", "seg_timing", "seg_length",
                                      "adjustment_factor", "late_adjusted", "seg_timing_adjusted", "adjusted_total_muts_on_segment",
                                      "error", "loCI", "hiCI")])
  
  # subset for regions of interest if that wasnt already done before running this function
  timed_data_region <- timed_data[timed_data$sample %in% regions, ]
  
  # build minimum consistent regions
  segTable.gr   <- GenomicRanges::makeGRangesFromDataFrame(timed_data_region, seqnames.field = "chr", start.field = "seg_start", end.field = "seg_end")
  minConsReg.gr <- disjoin(segTable.gr)
  minConsReg    <- as.data.frame(minConsReg.gr)
  # minConsReg$segment <- seq(1, nrow(minConsReg), by = 1)
  
  # need to overlap the original data with the minimum consistent regions
  mcols(segTable.gr) <- data.frame(timed_data_region)
  overlap <- findOverlaps(segTable.gr, minConsReg.gr)
  
  timed_min_cons <- cbind(timed_data_region[queryHits(overlap), ], minConsReg[subjectHits(overlap), ])
  timed_min_cons$width <- NULL
  timed_min_cons$strand <- NULL
  
  # test for differences between regions
  timed_min_cons$inf_col <- paste0(timed_min_cons$early, "_", timed_min_cons$late_adjusted, "_", timed_min_cons$adjusted_total_muts_on_segment, "_", timed_min_cons$seg_timing_adjusted)
  timed_min_cons_reg <- unique(timed_min_cons[, c("sample", "segID", "phase", "inf_col", "seqnames", "start", "end")])
  timed_min_cons_reg <- dcast(timed_min_cons_reg, segID + phase + seqnames + start + end ~ sample, value.var = "inf_col")
  
  # separate the column with necessary info
  if(colnames(timed_min_cons_reg)[6] %in% regions[1]){
    region1 <- regions[1]
    region2 <- regions[2]
    colnames(timed_min_cons_reg) <- c("segID", "phase", "seqnames", "start", "end", "region1", "region2")
  } else {
    region1 <- regions[1]
    region2 <- regions[2]
    colnames(timed_min_cons_reg) <- c("segID", "phase", "seqnames", "start", "end", "region2", "region1")
  }
  
  timed_min_cons_reg <- timed_min_cons_reg %>% separate(region1, into = c("region1:early", "region1:late_adjusted", "region1:adjusted_total_muts_on_segment", "region1:seg_timing_adjusted"), sep = "_")
  timed_min_cons_reg <- timed_min_cons_reg %>% separate(region2, into = c("region2:early", "region2:late_adjusted", "region2:adjusted_total_muts_on_segment", "region2:seg_timing_adjusted"), sep = "_")
  
  timed_min_cons_reg$`region1:early` <- as.numeric(timed_min_cons_reg$`region1:early`)
  timed_min_cons_reg$`region1:late_adjusted` <- as.numeric(timed_min_cons_reg$`region1:late_adjusted`)
  timed_min_cons_reg$`region1:adjusted_total_muts_on_segment` <- as.numeric(timed_min_cons_reg$`region1:adjusted_total_muts_on_segment`)
  timed_min_cons_reg$`region1:seg_timing_adjusted` <- as.numeric(timed_min_cons_reg$`region1:seg_timing_adjusted`)
  timed_min_cons_reg$`region2:early` <- as.numeric(timed_min_cons_reg$`region2:early`)
  timed_min_cons_reg$`region2:late_adjusted` <- as.numeric(timed_min_cons_reg$`region2:late_adjusted`)
  timed_min_cons_reg$`region2:adjusted_total_muts_on_segment` <- as.numeric(timed_min_cons_reg$`region2:adjusted_total_muts_on_segment`)
  timed_min_cons_reg$`region2:seg_timing_adjusted` <- as.numeric(timed_min_cons_reg$`region2:seg_timing_adjusted`)
  
  # remove rows with NA because let's only test where timing was possible in both regions
  timed_min_cons_reg_test <- na.omit(timed_min_cons_reg)

  # do a prop test to test for differences between the two regions
  timed_min_cons_reg_test <- timed_min_cons_reg_test %>%
    rowwise() %>%
    mutate(prop_test_p = prop.test(x = c(`region1:early`, `region2:early`), n = c(`region1:adjusted_total_muts_on_segment`, `region2:adjusted_total_muts_on_segment`))$p.value)
  
  timed_min_cons_reg_test$p_adjust <- p.adjust(timed_min_cons_reg_test$prop_test_p, method = "fdr")
  
  # make a plot indicating where there is a significant difference
  # for that let's add this back to the big dataframe
  timed_min_cons_reg <- left_join(timed_min_cons_reg, timed_min_cons_reg_test)
  timed_min_cons_reg_p <- timed_min_cons_reg[, c("segID", "phase", "seqnames", "start", "end", "region1:seg_timing_adjusted",
                                                 "region2:seg_timing_adjusted", "prop_test_p", "p_adjust")]
  
  timed_min_cons_reg_p <- melt(timed_min_cons_reg_p, id.vars = c("segID", "phase", "seqnames", "start", "end", "prop_test_p", "p_adjust"))
  timed_min_cons_reg_p$variable <- sub(":seg_timing_adjusted", "", timed_min_cons_reg_p$variable)
  colnames(timed_min_cons_reg_p) <- c("segID", "phase", "seqnames", "start", "end", "prop_test_p", "p_adjust", "region", "seg_timing_adjusted")

  # now seperate this out for the two haplotypes
  timed_min_cons_reg_p$info_col <- paste0(timed_min_cons_reg_p$seg_timing_adjusted, "_", timed_min_cons_reg_p$prop_test_p, "_", timed_min_cons_reg_p$p_adjust)
  timed_min_cons_reg_p <- reshape2::dcast(timed_min_cons_reg_p, segID + seqnames + start + end + region ~ phase, value.var = "info_col" )
  
  timed_min_cons_reg_p <- timed_min_cons_reg_p %>% separate(A, into = c("A_seg_timing_adjusted", "A_prop_test_p", "A_p_adjust"), sep = "_")
  timed_min_cons_reg_p <- timed_min_cons_reg_p %>% separate(B, into = c("B_seg_timing_adjusted", "B_prop_test_p", "B_p_adjust"), sep = "_")
  
  timed_min_cons_reg_p[which(timed_min_cons_reg_p$A_p_adjust <= 0.05), "significanceA"] <- "significant"
  timed_min_cons_reg_p[which(timed_min_cons_reg_p$A_p_adjust > 0.05), "significanceA"] <- "not_significant"
  timed_min_cons_reg_p[is.na(timed_min_cons_reg_p$significanceA), "significanceA"] <- "not_significant"
  
  timed_min_cons_reg_p[which(timed_min_cons_reg_p$B_p_adjust <= 0.05), "significanceB"] <- "significant"
  timed_min_cons_reg_p[which(timed_min_cons_reg_p$B_p_adjust > 0.05), "significanceB"] <- "not_significant"
  timed_min_cons_reg_p[is.na(timed_min_cons_reg_p$significanceB), "significanceB"] <- "not_significant"
  
  timed_min_cons_reg_p$significanceA <- factor(timed_min_cons_reg_p$significanceA, levels = c("significant", "not_significant"))
  timed_min_cons_reg_p$significanceB <- factor(timed_min_cons_reg_p$significanceB, levels = c("significant", "not_significant"))
  
  colors <- get.colors(length(regions))
  
  # make a dummy haplotype collumn for plotting purposes
  timed_min_cons_reg_p$haplotype <- "haplotype"
  
  # make the genome continuos
  colnames(timed_min_cons_reg_p) <- c("segID", "chr", "seg_start", "seg_end", "region", "A_seg_timing_adjusted", "A_prop_test_p", "A_p_adjust",
                                      "B_seg_timing_adjusted", "B_prop_test_p", "B_p_adjust", "significanceA", "significanceB", "haplotype")
  
  timed_min_cons_reg_p <- genome.continous.fun(timed_min_cons_reg_p, genome)
  
  chr_sizes <- seqlengths(genome)[paste0("chr", 1:22)]
  chr_lines <- cumsum(as.numeric(as.character(chr_sizes)))
  chr_midpoint <- chr_lines - chr_sizes / 2
  xaxis_label <- data.frame(chr_midpoint)
  xaxis_label$label <- rownames(xaxis_label)
  
  timed_min_cons_reg_p$A_seg_timing_adjusted <- as.numeric(as.character(timed_min_cons_reg_p$A_seg_timing_adjusted))
  timed_min_cons_reg_p$B_seg_timing_adjusted <- as.numeric(as.character(timed_min_cons_reg_p$B_seg_timing_adjusted))
  timed_min_cons_reg_p[is.na(timed_min_cons_reg_p$A_seg_timing_adjusted), "A_seg_timing_adjusted"] <- 0
  timed_min_cons_reg_p[is.na(timed_min_cons_reg_p$B_seg_timing_adjusted), "B_seg_timing_adjusted"] <- 0
  
  timed_min_cons_reg_p$region <- factor(timed_min_cons_reg_p$region, levels = c("region1", "region2"))
  
  p_genome1 <- ggplot() +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    geom_segment(data = timed_min_cons_reg_p,
                 aes(x = adjust_seg_start, xend = adjust_seg_end, y = A_seg_timing_adjusted, yend = A_seg_timing_adjusted, color = region, alpha = significanceA),
                 linewidth = 2) +
    scale_color_manual(name = "tumour region",
                       breaks = c("region1", "region2"),
                       values = c(colors),
                       labels = c(region1, region2)) +
    scale_alpha_manual(name = "", breaks = c("significant", "not_significant"),
                       values = c("significant" = 1,
                                  "not_significant" = 0.3),
                       labels = c("significant", "not significant\n/gain timing not possible\nfor both regions")) +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    facet_wrap(~haplotype, labeller = labeller(haplotype = c(haplotype = "haplotype A")), strip.position = "right") +
    ylim(c(0, 1)) +
    ylab("gain timing") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm"))
  
  p_genome2 <- ggplot() +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    geom_segment(data = timed_min_cons_reg_p,
                 aes(x = adjust_seg_start, xend = adjust_seg_end, y = B_seg_timing_adjusted, yend = B_seg_timing_adjusted, color = region, alpha = significanceB),
                 linewidth = 2) +
    scale_color_manual(name = "tumour region",
                       values = c(colors)) +
    scale_alpha_manual(name = "", breaks = c("significant", "not_significant"),
                       values = c("significant" = 1,
                                  "not_significant" = 0.3),
                       labels = c("significant", "not significant\n/gain timing not possible\nfor both regions")) +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    facet_wrap(~haplotype, labeller = labeller(haplotype = c(haplotype = "haplotype B")), strip.position = "right") +
    ylim(c(0, 1)) +
    ylab("gain timing") +
    xlab("") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          plot.margin = unit(c(0, 0.1, 0.1, 0.1), "cm"))
  
  # combine the two plots
  p_genome_c <- plot_grid(p_genome1, p_genome2, ncol = 1, align = "v", axis = "rl", rel_heights = c(1, 1.1))
  
  if(!dir.exists(paste0(output_path, tumour, "/", "/plots/comparing_regions_phased_mutation_timing"))){
    dir.create(paste0(output_path, tumour, "/", "/plots/comparing_regions_phased_mutation_timing"), recursive = T)
  }
  
  pdf(paste0(output_path, tumour, "/plots/comparing_regions_phased_mutation_timing/", tumour, "_genome_phased_gain_timing_per_haplotype_regions_significance.pdf"), width = 18, height = 7)
  print(p_genome_c)
  dev.off()
  
  return(timed_min_cons_reg_p)
  
}

# function to get a sequence of timing across regions
timing.sequence.fun <- function(all_region_data = all_timed_data, output_path, genome){
  
  # keep only relevant data
  all_region_data <- unique(all_region_data[, c("sample", "chr", "seg_start", "seg_end", "segID", "seg_timing_adjusted", "phase")])
  
  # have to build minimum consistent regions of segments
  # build minimum consistent regions
  segTable.gr   <- GenomicRanges::makeGRangesFromDataFrame(all_region_data, seqnames.field = "chr", start.field = "seg_start", end.field = "seg_end")
  minConsReg.gr <- disjoin(segTable.gr)
  minConsReg    <- as.data.frame(minConsReg.gr)

  # need to overlap the original data with the minimum consistent regions
  mcols(segTable.gr) <- data.frame(all_region_data)
  overlap <- findOverlaps(segTable.gr, minConsReg.gr)
  
  timed_min_cons <- cbind(all_region_data[queryHits(overlap), ], minConsReg[subjectHits(overlap), ])
  timed_min_cons$width <- NULL
  timed_min_cons$strand <- NULL
  timed_min_cons$chr <- NULL
  timed_min_cons$seg_start <- NULL
  timed_min_cons$seg_end <- NULL
  timed_min_cons$segID <- NULL
  timed_min_cons <- unique(timed_min_cons)
  
  # if we have a timing estimate for both regions of one tumour for a segment take the mean of those twp
  timed_min_cons$min_cons_segID <- paste0(timed_min_cons$seqnames, "_", timed_min_cons$start, "_", timed_min_cons$end)
  timed_min_cons$tumour <- sub("_.*", "", timed_min_cons$sample)
  
  timed_min_cons <- timed_min_cons %>% group_by(tumour, phase, min_cons_segID) %>%
                      mutate(mean_timing = mean(seg_timing_adjusted, na.rm = T))
  
  timed_min_cons$sample <- NULL
  timed_min_cons$seg_timing_adjusted <- NULL
  timed_min_cons <- unique(timed_min_cons)
  
  # for each tumour, determine a sequence of events
  timed_min_cons_seq <- lapply(unique(timed_min_cons$tumour), function(x){
    
    df <- timed_min_cons[timed_min_cons$tumour == x, ]
    
    dfA <- df[df$phase == "A",]
    dfA <- dfA[order(dfA$mean_timing, decreasing = F), ]
    dfA <- na.omit(dfA)
    for(i in unique(dfA$mean_timing)){
      event_number <- which(unique(dfA$mean_timing) == i)
      dfA[dfA$mean_timing == i, "event_order"] <- event_number
    }
    
    dfB <- df[df$phase == "B",]
    dfB <- dfB[order(dfB$mean_timing, decreasing = F), ]
    dfB <- na.omit(dfB)
    for(i in unique(dfB$mean_timing)){
      event_number <- which(unique(dfB$mean_timing) == i)
      dfB[dfB$mean_timing == i, "event_order"] <- event_number
    }
    
    dfAB <- rbind(dfA, dfB)
    return(dfAB)
  }) 
  
  timed_min_cons_seq <- do.call(rbind, timed_min_cons_seq)
  
  # gotta make an order from all tumours
  timed_min_cons_seq$event_order <- as.numeric(as.character(timed_min_cons_seq$event_order))
  timed_min_cons_seq <- data.frame(timed_min_cons_seq)
  
  timed_min_cons_seq$tumour <- NULL
  timed_min_cons_seq$mean_timing <- NULL
  timed_min_cons_seq <- unique(timed_min_cons_seq)
  
  event_sequence <- timed_min_cons_seq %>% group_by(phase, min_cons_segID) %>%
                      mutate(median_event_order = median(event_order)) %>%
                      mutate(mean_event_order = mean(event_order)) %>%
                      mutate(iqr = IQR(event_order)) %>%
                      mutate(num_tumours = n())
  
  event_sequence$min_IQR <- event_sequence$median_event_order - (event_sequence$iqr/2)
  event_sequence$max_IQR <- event_sequence$median_event_order + (event_sequence$iqr/2)
  
  # where min IQR is negative set it to 0
  event_sequence[event_sequence$min_IQR < 0, "min_IQR"] <- 0
  
  event_sequence$event_order <- NULL
  event_sequence <- unique(event_sequence)
  
  # plot this
  colnames(event_sequence) <- c("phase", "chr", "seg_start", "seg_end", "min_cons_segID", "median_event_order",
                                "mean_event_order", "iqr", "num_tumours", "min_IQR", "max_IQR")
  
  event_sequence <- genome.continous.fun(event_sequence, genome)
  
  chr_sizes <- seqlengths(genome)[paste0("chr", 1:22)]
  chr_lines <- cumsum(as.numeric(as.character(chr_sizes)))
  chr_midpoint <- chr_lines - chr_sizes / 2
  xaxis_label <- data.frame(chr_midpoint)
  xaxis_label$label <- rownames(xaxis_label)
  
  p_genome <- ggplot() +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    geom_rect(data = event_sequence, aes(xmin = adjust_seg_start, xmax = adjust_seg_end, ymin = min_IQR, ymax = max_IQR, fill = phase), color = NA, alpha = 0.15) +
    scale_fill_manual(values = c("A" = "darkred",
                                  "B" = "#313695"), guide = "none") +
    geom_segment(data = event_sequence,
                 aes(x = adjust_seg_start, xend = adjust_seg_end, y = median_event_order, yend = median_event_order,
                     color = phase), linewidth = 2) +
    scale_color_manual(values = c("A" = "darkred",
                                  "B" = "#313695"), guide = "none") +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    facet_wrap(~phase, labeller = labeller(phase = c(A = "haplotype A", B = "haplotype B")), strip.position = "right", nrow = 2) +
    ylab("event history") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
  
  if(!dir.exists(paste0(output_path, "/plots/"))){
    dir.create(paste0(output_path, "/plots/"), recursive = T)
  }
  
  pdf(paste0(output_path, "/plots/", "genome_phased_gain_timing_per_haplotype_gain_sequence_across_tumours.pdf"), width = 18, height = 7)
  print(p_genome)
  dev.off()

}

# haplo timing heat map

timing.heatmap.fun <- function(all_region_data = all_timed_data, output_path, genome){
  
  # keep only relevant data
  all_region_data <- all_region_data %>% separate(segID, into = c("chr", "seg_start", "seg_end"), sep = "_", remove = F)
  all_region_data <- all_region_data[which(all_region_data$gain_index == 1), ]
 
  # have to build minimum consistent regions of segments
  # build minimum consistent regions
  segTable.gr   <- GenomicRanges::makeGRangesFromDataFrame(all_region_data, seqnames.field = "chr", start.field = "seg_start", end.field = "seg_end")
  minConsReg.gr <- disjoin(segTable.gr)
  minConsReg    <- as.data.frame(minConsReg.gr)
  
  # need to overlap the original data with the minimum consistent regions
  mcols(segTable.gr) <- data.frame(all_region_data)
  overlap <- findOverlaps(segTable.gr, minConsReg.gr)
  
  timed_min_cons <- cbind(all_region_data[queryHits(overlap), ], minConsReg[subjectHits(overlap), ])
  timed_min_cons$width <- NULL
  timed_min_cons$strand <- NULL
  timed_min_cons$chr <- NULL
  timed_min_cons$seg_start <- NULL
  timed_min_cons$seg_end <- NULL
  timed_min_cons$segID <- NULL
  timed_min_cons$copy_number_before_gain <- NULL
  timed_min_cons$copy_number_after_gain <- NULL
  timed_min_cons$number_early_muts_cn_adjusted <- NULL
  timed_min_cons$number_all_muts_cn_adjusted <- NULL
  timed_min_cons$gain_index <- NULL
  timed_min_cons$cluster <- NULL
  timed_min_cons <- unique(timed_min_cons)
  
  
  # plot this
  colnames(timed_min_cons) <- c("tumour", "region", "gain_timing", "loCI", "hiCI", "phase", "chr", "seg_start", "seg_end")
  
  timed_min_cons <- genome.continous.fun(timed_min_cons, genome)
  
  chr_sizes <- seqlengths(genome)[paste0("chr", 1:22)]
  chr_lines <- cumsum(as.numeric(as.character(chr_sizes)))
  chr_midpoint <- chr_lines - chr_sizes / 2
  xaxis_label <- data.frame(chr_midpoint)
  xaxis_label$label <- rownames(xaxis_label)
  
  sample_levels <- data.frame(region = unique(timed_min_cons$region))
  sample_levels$index <- c(1:length(sample_levels$region))
  
  timed_min_cons <- left_join(timed_min_cons, sample_levels)

  tumour_lines <- sample_levels
  tumour_lines$tumour <- sub("_.*", "", tumour_lines$region)
  tumour_lines <- tumour_lines %>% group_by(tumour) %>% mutate(n = n())
  tumour_lines <- data.frame(tumour_lines)
  tumour_lines$region <- NULL
  tumour_lines$index <- NULL
  tumour_lines <- unique(tumour_lines)
  tumour_lines$line_index <- cumsum(tumour_lines$n)
  tumour_lines$line_index <- tumour_lines$line_index + 1


  p_heat <- ggplot() +
    geom_rect(data = timed_min_cons, aes(xmin = adjust_seg_start, xmax = adjust_seg_end, ymin = index, ymax = index + 1, fill = gain_timing)) +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "black") +
    geom_hline(yintercept = tumour_lines$line_index, linetype = 1, color = "#ffffff") +
    scale_fill_gradientn(name = "gain timing", colours = c("#8e0152", "#c51b7d", "#de77ae", "#f1b6da", "#fde0ef", "#e6f5d0", "#b8e186", "#7fbc41", "#4d9221", "#276419"),
                         breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0 (early)", 0.25, 0.5, 0.75, "1 (late)"), na.value = "#bababa", limits = c(0, 1)) +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label) +
    scale_y_continuous(breaks = unique(timed_min_cons$index) + 0.5, labels = unique(timed_min_cons$region), expand = c(0,0)) +
    facet_wrap(~phase, labeller = labeller(phase = c(A = "haplotype A", B = "haplotype B")), strip.position = "right", nrow = 2) +
    ylab("") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 16),
          axis.text.y = element_text( size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
          panel.background = element_rect(fill = "#d9d9d9"))
  
  if(!dir.exists(paste0(output_path, "/plots/"))){
    dir.create(paste0(output_path, "/plots/"), recursive = T)
  }
  
  pdf(paste0(output_path, "/plots/", "genome_phased_gain_timing_per_haplotype_heatmap_across_tumours_clonaL_first_gain.pdf"), width = 18, height = 7)
  print(p_heat)
  dev.off()
  
}

#function to analyse parallel gains within the same tumour region
analyse.parallel.gains <- function(all_region_data = all_timed_data, output_path){

  para_gains <- lapply(unique(all_region_data$region), function(region){

    region_df <- all_region_data[all_region_data$region == region, ]
    region_df <- region_df[region_df$gain_index == 1, ]
    seg_count <- data.frame(table(region_df$segID))
    para_seg  <- as.character(seg_count[seg_count$Freq == 2, "Var1"])
    
    if(length(para_seg) > 0){
      para_df   <- region_df[region_df$segID %in% para_seg, ]

      # what is the difference in gain timing between the two haplotypes
      time_df <- para_df[, c("tumour", "region", "segID", "gain_timing", "phase")]
      time_df <- reshape2::dcast(time_df, tumour + region + segID ~ phase, value.var = "gain_timing")
      time_df$difference <- abs(time_df$A - time_df$B)

      loci_df <- para_df[, c("tumour", "region", "segID", "loCI", "phase")]
      loci_df <- reshape2::dcast(loci_df, tumour + region + segID ~ phase, value.var = "loCI")
      colnames(loci_df) <- c("tumour", "region", "segID", "loCI_A", "loCI_B")

      hici_df <- para_df[, c("tumour", "region", "segID", "hiCI", "phase")]
      hici_df <- reshape2::dcast(hici_df, tumour + region + segID ~ phase, value.var = "hiCI")
      colnames(hici_df) <- c("tumour", "region", "segID", "hiCI_A", "hiCI_B")
  
      df_out <- left_join(time_df, loci_df)
      df_out <- left_join(df_out, hici_df)
 
      #are the confidence intervals overlapping
      df_out <- df_out %>% rowwise() %>% mutate(CI_overlap = ifelse((loCI_B < loCI_A & hiCI_B > hiCI_A) |
                                                                    (loCI_A < loCI_B & hiCI_A > hiCI_B) |
                                                                    (loCI_A < loCI_B & hiCI_A > loCI_B) |
                                                                    (loCI_B < loCI_A & hiCI_B > loCI_A), TRUE, FALSE))
      df_out <- data.frame(df_out)

      return(df_out) 
    }
  })
  para_gains <- do.call(rbind, para_gains)
  para_gains$width_CI_A <- abs(para_gains$loCI_A - para_gains$hiCI_A)
  para_gains$width_CI_B <- abs(para_gains$loCI_B - para_gains$hiCI_B)
  para_gains <- para_gains %>% rowwise() %>% mutate(maxCI_with = max(width_CI_A, width_CI_B)) 
  para_gains <- data.frame(para_gains)

  #plot
  plot_df <- melt(para_gains, id.vars = c("tumour", "region", "segID", "A", "B", "loCI_A", "loCI_B", "hiCI_A", "hiCI_B", "CI_overlap", "width_CI_A", "width_CI_B"))
  plot_df$variable <- factor(plot_df$variable, levels = c("maxCI_with", "difference"))
  plot_df$CI_overlap <- factor(plot_df$CI_overlap, levels = c(TRUE, FALSE))

  p <- ggplot() +
        geom_point(data = plot_df, aes(segID, value, color = variable, shape = CI_overlap)) +
        scale_color_manual(values = c("red", "blue"), guide = "none") +
        scale_shape_manual(name = "CIs overlap", breaks = c(TRUE, FALSE), values = c(16, 17)) +
        facet_grid(variable~region, scales = "free",switch = "y", labeller = labeller(variable = c(maxCI_with = "maximum\ngain timing\nCI width", difference = "gain timing\ndifference"))) +
        ylab("") +
        xlab("segment") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              strip.background = element_blank(),
              strip.placement = "outside")

  if(!dir.exists(paste0(output_path, "/plots/"))){
    dir.create(paste0(output_path, "/plots/"), recursive = T)
  }
  
  pdf(paste0(output_path, "/plots/", "parallel_gains_comparison.pdf"), width = 16, height = 6)
  print(p)
  dev.off()

}

analyse.clonal.gains <- function(all_region_data = all_timed_data, output_path){

  #which ones are the tumours with several regions
  tum_df <- unique(all_region_data[, c("tumour", "region")])
  dups <- tum_df$tumour[duplicated(tum_df$tumour)]

  clonal_gains <- lapply(dups, function(tumour){
  
    tumour_df <- all_region_data[all_region_data$tumour == tumour, ]
    tumour_df <- tumour_df[which(tumour_df$gain_index == 1), ]
    seg_count <- data.frame(table(tumour_df$segID))
    para_segs <- as.character(seg_count[seg_count$Freq > 1, "Var1"])

    para_df <- tumour_df[tumour_df$segID %in% para_segs, ]

    #they need to be present in both regions not just both haplotypes
    clonal_df <- list()
    for(segID in unique(para_df$segID)){
      
      df <- para_df[para_df$segID == segID, ]

      if(length(unique(df$region)) == 2){
        clonal_df[[segID]] <- df
      }

    }
    clonal_df <- do.call(rbind, clonal_df)

    if(length(clonal_df) > 0){

      # we have to deal with parallel gains on top of clonal gains
      ct <- data.frame(table(clonal_df$segID))
      keep <- as.character(ct[ct$Freq == 2, "Var1"])
      
      clonal_df <- clonal_df[clonal_df$segID %in% keep, ]

      if(length(clonal_df) > 0){

        # what is the difference in gain timing between the two regions
        time_df <- clonal_df[, c("tumour", "region", "segID", "gain_timing")]
        time_df <- reshape2::dcast(time_df, tumour + segID ~ region, value.var = "gain_timing")
        region_key <- colnames(time_df[3:4])
        colnames(time_df) <- c("tumour", "segID", "timing_region1", "timing_region2")

        loci_df <- clonal_df[, c("tumour", "region", "segID", "loCI")]
        loci_df <- reshape2::dcast(loci_df, tumour + segID ~ region, value.var = "loCI")
        colnames(loci_df) <- c("tumour", "segID", "loCI_region1", "loCI_region2")

        hici_df <- clonal_df[, c("tumour", "region", "segID", "hiCI")]
        hici_df <- reshape2::dcast(hici_df, tumour + segID ~ region, value.var = "hiCI")
        colnames(hici_df) <- c("tumour", "segID", "hiCI_region1", "hiCI_region2")
    
        df_out <- left_join(time_df, loci_df)
        df_out <- left_join(df_out, hici_df)
        df_out$difference <- abs(df_out$timing_region1 - df_out$timing_region2)
  
        #are the confidence intervals overlapping
        df_out <- df_out %>% rowwise() %>% mutate(CI_overlap = ifelse((loCI_region2 < loCI_region1 & hiCI_region2 > hiCI_region1) |
                                                                      (loCI_region1 < loCI_region2 & hiCI_region1 > hiCI_region2) |
                                                                      (loCI_region1 < loCI_region2 & hiCI_region1 > loCI_region2) |
                                                                      (loCI_region2 < loCI_region1 & hiCI_region2 > loCI_region1), TRUE, FALSE))
        df_out <- data.frame(df_out)

        return(df_out) 
      }
    }

  })
  clonal_gains <- do.call(rbind, clonal_gains)

  clonal_gains$width_CI_region1 <- abs(clonal_gains$hiCI_region1 - clonal_gains$loCI_region1)
  clonal_gains$width_CI_region2 <- abs(clonal_gains$hiCI_region2 - clonal_gains$loCI_region2)
  clonal_gains <- clonal_gains %>% rowwise() %>% mutate(maxCI_with = max(width_CI_region1, width_CI_region2)) 
  clonal_gains <- data.frame(clonal_gains)

  #plot
  plot_df <- melt(clonal_gains, id.vars = c("tumour", "segID", "timing_region1", "timing_region2", "loCI_region1", "loCI_region2", "hiCI_region1", "hiCI_region2", "CI_overlap", "width_CI_region1", "width_CI_region2"))
  plot_df$variable <- factor(plot_df$variable, levels = c("maxCI_with", "difference"))
  plot_df$CI_overlap <- factor(plot_df$CI_overlap, levels = c(TRUE, FALSE))

  p <- ggplot() +
        geom_point(data = plot_df, aes(segID, value, color = variable, shape = CI_overlap)) +
        scale_color_manual(values = c("red", "blue"), guide = "none") +
        scale_shape_manual(name = "CIs overlap", breaks = c(TRUE, FALSE), values = c(16, 17)) +
        facet_grid(variable~tumour, scales = "free",switch = "y", labeller = labeller(variable = c(maxCI_with = "maximum\ngain timing\nCI width", difference = "gain timing\ndifference"))) +
        ylab("") +
        xlab("segment") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              strip.background = element_blank(),
              strip.placement = "outside")

  if(!dir.exists(paste0(output_path, "/plots/"))){
    dir.create(paste0(output_path, "/plots/"), recursive = T)
  }
  
  pdf(paste0(output_path, "/plots/", "clonal_gains_comparison.pdf"), width = 20, height = 6)
  print(p)
  dev.off()


}

################################################################################
####################################################################      MAIN

all_timed_data <- list()
for(region in all_regions){

  tumour = strsplit(region, "_")[[1]][1]

  filepath <- paste0(input_path, tumour, "/", region, "/data/", "clonal_gain_timing.txt")

  if(file.exists(filepath)){
    df <- read.table(filepath, head = T, sep = "\t")
    all_timed_data[[region]] <- df
  }
}

all_timed_data <- do.call(rbind, all_timed_data)
all_timed_data <- all_timed_data[!is.na(all_timed_data$gain_timing), ]
all_timed_data <- all_timed_data[all_timed_data$gain_index == 1, ]


diff <- list()

for(tumour in unique(all_timed_data$tumour)){
  
  df <- all_timed_data[all_timed_data$tumour == tumour, ]
  dfA <- df[df$haplotype == "A", ]
  dfB <- df[df$haplotype == "B", ]
  
  if(nrow(dfA) > 0){
    reg_compareA <- dcast(dfA, segID ~ region, value.var = "gain_timing")
  }
  
  if(nrow(dfB) > 0){
    reg_compareB <- dcast(dfB, segID ~ region, value.var = "gain_timing")
  }
  
  if(all(ncol(reg_compareA) == 3 & ncol(reg_compareB) == 3)){
    
    both <- rbind(reg_compareA, reg_compareB)
    
  } else {
    
    if(ncol(reg_compareA) == 3){
      reg_compareA <- dcast(dfA, segID ~ region, value.var = "gain_timing")
      both <- reg_compareA
    }
    
    if(ncol(reg_compareB) == 3){
      
      reg_compareB <- dcast(dfB, segID ~ region, value.var = "gain_timing")
      both <- reg_compareB
    }
    
  }
  
  colnames(both) <- c("segID", "region1", "region2")
  
  diff[[tumour]] <- both
  
}

diff <- do.call(rbind, diff)
diff$diff <- abs(diff$region1 - diff$region2)

p_hist <- ggplot(diff, aes(x = diff)) + 
            geom_histogram(color = "darkgrey", fill = "darkblue") +
            xlab("gain timing difference") +
            theme_bw() 

pdf(paste0(output_path, "plots/hist_timing_difference_regions.pdf"), width = 5, height = 5)
p_hist
dev.off()

# now check if CIs overlap

ci_comp <- list()

for(tumour in unique(all_timed_data$tumour)){
  
  df <- all_timed_data[all_timed_data$tumour == tumour, ]
  df$error <- abs(df$loCI - df$gain_timing)
  
  dfA <- df[df$haplotype == "A", ]
  dfB <- df[df$haplotype == "B", ]
  
  if(nrow(dfA) > 0){

    reg_compareA <- dcast(dfA, segID ~ region, value.var = "gain_timing")
    
    if(ncol(reg_compareA) == 3){
      colnames(reg_compareA) <- c("segID", "region1", "region2")
      reg_compareA$diff <- abs(reg_compareA$region1 - reg_compareA$region2)
      
      reg_errorA <- dcast(dfA, segID ~ region, value.var = "error")
      colnames(reg_errorA) <- c("segID", "region1", "region2")
      reg_errorA$se_diff <- sqrt(reg_errorA$region1^2 + reg_errorA$region2^2)
      reg_errorA$region1 <- NULL
      reg_errorA$region2 <- NULL
      
      reg_compareA <- left_join(reg_compareA, reg_errorA)
      reg_compareA$z <- reg_compareA$diff / reg_compareA$se_diff
      reg_compareA$p <- 2 * pnorm(-abs(reg_compareA$z))
    }
    
  }
  
  if(nrow(dfB) > 0){
    
    reg_compareB <- dcast(dfB, segID ~ region, value.var = "gain_timing")
    
    if(ncol(reg_compareB) == 3){
      colnames(reg_compareB) <- c("segID", "region1", "region2")
      reg_compareB$diff <- abs(reg_compareB$region1 - reg_compareB$region2)
      
      reg_errorB <- dcast(dfB, segID ~ region, value.var = "error")
      colnames(reg_errorB) <- c("segID", "region1", "region2")
      reg_errorB$se_diff <- sqrt(reg_errorB$region1^2 + reg_errorB$region2^2)
      reg_errorB$region1 <- NULL
      reg_errorB$region2 <- NULL
      
      reg_compareB <- left_join(reg_compareB, reg_errorB)
      reg_compareB$z <- reg_compareB$diff / reg_compareB$se_diff
      reg_compareB$p <- 2 * pnorm(-abs(reg_compareB$z))
    }
  }
  
  if(all(ncol(reg_compareA) == 3 & ncol(reg_compareB) == 3)){
    
    both <- rbind(reg_compareA, reg_compareB)
    
  } else {
    
    if(ncol(reg_compareA) == 7){
      both <- reg_compareA
    }
    
    if(ncol(reg_compareB) == 7){
      both <- reg_compareB
    }
  }
  
  ci_comp[[tumour]] <- both
  
}
ci_comp <- do.call(rbind, ci_comp)

p_hist_p <- ggplot(ci_comp, aes(x = p)) + 
  geom_histogram(color = "darkgrey", fill = "darkblue") +
  xlab("p-value") +
  geom_vline(xintercept = 0.05, color = "darkred") +
  theme_bw() 

pdf(paste0(output_path, "plots/hist_timing_pvalue_difference_regions.pdf"), width = 5, height = 5)
p_hist_p
dev.off()








