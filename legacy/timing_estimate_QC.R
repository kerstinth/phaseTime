# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #           check allele specific timing independence                   # # #  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# author: Kerstin Thol

################################################################################
#####################################################################  libraries
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(data.table)
library(vcfR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(GenomicRanges)

################################################################################
##################################################################### file paths


################################################################################
##################################################################### functions
source("haplo_timing_functions.R")

################################################################################
#####################################################################     MAIN

# read in data and run timing functions
qc_sheet <- read.table(battenberg_qc_path, head = T, sep = ",")
qc_sheet <- qc_sheet[qc_sheet$manual_call == "OK", ]

mutation_data_all <- fread(mutation_data_base_path, head = T, sep = "\t", data.table = FALSE)

all_timed_data <- list()
for(region in all_regions){
  
  print(region)
  
  tumour = strsplit(region, "_")[[1]][1]
  
  # get data
  mutation_data <- mutation_data_all[mutation_data_all$sample == region, ]
  
  #gotta load the qc'd copy number, this should eventuell be streamlined
  reg_qc <- qc_sheet[qc_sheet$region == region, ]
  region_bb_path <- list.files(path = paste0(battenberg_base_path, tumour, "/Battenberg/aggregated/"), pattern = paste0(reg_qc$solution, "_subclones.txt.gz"), full.names = T)
  region_bb_path <- region_bb_path[grep(region, region_bb_path)]
  cn <- read.table(region_bb_path, head = T, sep = "\t")
  
  #gotta merge the subclones
  cn$nAraw <- sapply(1:nrow(cn), function(x){
    ifelse(cn$frac1_A[x]==1,
           cn$nMaj1_A[x],
           cn$nMaj1_A[x]*cn$frac1_A[x]+cn$nMaj2_A[x]*cn$frac2_A[x])
  })
  
  cn$nBraw <- sapply(1:nrow(cn), function(x){
    ifelse(cn$frac1_A[x]==1,
           cn$nMin1_A[x],
           cn$nMin1_A[x]*cn$frac1_A[x]+cn$nMin2_A[x]*cn$frac2_A[x])
  })
  
  cn$nMajor   <- round(cn$nAraw)
  cn$nMinor   <- round(cn$nBraw)
  cn$cnTotal  <- cn$nMajor + cn$nMinor
  
  # flip this around already now
  #need to also read in the MSAI battenberg file
  msai_path <- paste0(battenberg_base_path, tumour, "/Battenberg/aggregated/multisample_MSAI.txt")
  msai <- read.table(msai_path, head = T, sep = "\t")
  
  # Step 1: Separate `frac_consensus` into rows
  df_separated <- msai %>%
    separate_rows(frac_consensus, sep = ";")
  
  # Step 2: Separate the key-value pairs
  df_separated <- df_separated %>%
    separate(frac_consensus, into = c("Key", "Value"), sep = "=") %>%
    mutate(Value = as.numeric(Value))
  
  df_separated$Key <- sub("--..........", "", df_separated$Key)
  
  # Step 3: Pivot the data to wide format
  msai_wide <- df_separated %>%
    pivot_wider(names_from = Key, values_from = Value)
  
  msai_wide <- data.frame(msai_wide)
  msai_long <- melt(msai_wide, id.vars = c("seqnames", "start", "end", "msai"))
  
  #where msai is TRUE and value is below 0.5 we have to swap
  msai_long[which(msai_long$value < 0.5), "flip"] <- TRUE
  msai_long[which(msai_long$value >= 0.5), "flip"] <- FALSE
  msai_long[is.na(msai_long$value), "flip"] <- FALSE
  
  #make one more collumn indicating if there is allelic imbalance present
  msai_long[is.na(msai_long$value), "allelic_imbalance"] <- FALSE
  msai_long[!is.na(msai_long$value), "allelic_imbalance"] <- TRUE
  
  #we also need to create a map from region specific CN segments to minimum consistent region segments
  mcr <- unique(msai_long[, c("seqnames", "start", "end")])
  colnames(mcr) <- c("chr", "seg_start", "seg_end")
  mcr$mcr <- TRUE
  mcr_gr <- GRanges(seqnames = mcr$chr, IRanges(start = mcr$seg_start, end = mcr$seg_end))
  mcols(mcr_gr) <- data.frame(mcr = TRUE, mcr_seg_start = mcr$seg_start, mcr_seg_end = mcr$seg_end)
  regionCN <- unique(cn[, c("chr", "startpos", "endpos")])
  regionCN <- regionCN[!is.na(regionCN$startpos), ]
  regionCN <- unique(regionCN)
  regionCN$regionCN <- TRUE
  regionCN_gr <- GRanges(seqnames = regionCN$chr, IRanges(start = regionCN$startpos, end = regionCN$endpos))
  mcols(regionCN_gr) <- data.frame(regionCN = TRUE, regionCN_seg_start = regionCN$startpos, regionCN_seg_end = regionCN$endpos)
  
  overlap         <-  findOverlaps(regionCN_gr, mcr_gr)
  intersect_gr    <- pintersect(mcr_gr[subjectHits(overlap)], regionCN_gr[queryHits(overlap)])
  mcols(intersect_gr) <- data.frame(mcols(intersect_gr), mcols(regionCN_gr[queryHits(overlap)]))
  intersect_df <- data.frame(intersect_gr)
  
  segment_map <- intersect_df[, c("seqnames", "regionCN_seg_start", "regionCN_seg_end", "mcr_seg_start", "mcr_seg_end")]
  
  cn <- left_join(cn, segment_map, by = c("chr" = "seqnames",
                                          "startpos" = "regionCN_seg_start",
                                          "endpos" = "regionCN_seg_end"))
  
  cn[!is.na(cn$mcr_seg_start), "startpos"] <-  cn[!is.na(cn$mcr_seg_start), "mcr_seg_start"]
  cn[!is.na(cn$mcr_seg_end), "endpos"] <-  cn[!is.na(cn$mcr_seg_end), "mcr_seg_end"]
  
  cn$mcr_seg_start <- NULL
  cn$mcr_seg_end <- NULL
  
  #we have to properly assign major and minor CN to haplotypes
  #separate the regions
  msai_long$region <- sub("_Major", "", msai_long$variable)
  msai_long$region <- sub("\\.", "-", msai_long$region)
  msai_long <- msai_long[msai_long$region == region, ]
  
  #map the start and end positions to the region specific start and end positions
  # for subsequent analysis let's use minimum consistent regions, might have to have a single regions and a multi region mode
  
  cn <- left_join(cn, msai_long[, c("seqnames", "start", "end", "msai", "flip", "allelic_imbalance")], by = c("chr" = "seqnames",
                                                                                                              "startpos" = "start",
                                                                                                              "endpos" = "end"))
  
  
  #sort out NAs in msai, flip and allelic_imbalance, these are porbably cases without AI in any region and is therefore not in MSAI table
  cn[is.na(cn$msai), "msai"] <- FALSE
  cn[is.na(cn$flip), "flip"] <- FALSE
  cn[is.na(cn$allelic_imbalance), "allelic_imbalance"] <- FALSE
  
  
  #now using the flip collumn we can assign copy number to A and B instead major and minor based on being the greater value
  cn[which(cn$flip == FALSE), "cnA"] <- cn[which(cn$flip == FALSE), "nAraw"]
  cn[which(cn$flip == FALSE), "cnB"] <- cn[which(cn$flip == FALSE), "nBraw"]
  cn[which(cn$flip == TRUE), "cnA"] <- cn[which(cn$flip == TRUE), "nBraw"]
  cn[which(cn$flip == TRUE), "cnB"] <- cn[which(cn$flip == TRUE), "nAraw"]
  
  #lets remove the old non-flipped copy number
  cn$nAraw <- NULL
  cn$nBraw <- NULL
  
  # and rename into major and minor, this makes no sense but prevents me having to change the code everywhere
  
  colnames(cn)[(ncol(cn) -1) : ncol(cn)] <- c("maj_cn", "min_cn")
  
  #have to change the copy number that's in the mutation_data file to this copy number
  muts_for_cn <- unique(mutation_data[, c("chr", "pos", "vaf")])
  muts_for_cn_gr <- GRanges(seqnames = muts_for_cn$chr, IRanges(start = muts_for_cn$pos, end = muts_for_cn$pos))
  mcols(muts_for_cn_gr)          <- data.frame(muts_for_cn[, c("vaf")])

  cn_for_muts <- cn[, c("chr", "startpos", "endpos", "maj_cn", "min_cn", "msai", "allelic_imbalance")]
  cn_for_muts_gr <- GRanges(seqnames = cn_for_muts$chr, IRanges(start = cn_for_muts$startpos, end = cn_for_muts$endpos))
  mcols(cn_for_muts_gr)          <- data.frame(cn_for_muts[, c("startpos", "endpos", "maj_cn", "min_cn", "msai", "allelic_imbalance")])
  
  overlap               <- findOverlaps(muts_for_cn_gr, cn_for_muts_gr)
  intersect_gr          <- pintersect(cn_for_muts_gr[subjectHits(overlap)], muts_for_cn_gr[queryHits(overlap)])
  mcols(intersect_gr)   <- data.frame(mcols(intersect_gr), mcols(muts_for_cn_gr[queryHits(overlap)]))
  intersect_df          <- data.frame(intersect_gr)
  intersect_df$hit      <- NULL
  intersect_df$width    <- NULL
  intersect_df$strand    <- NULL
  intersect_df$end      <- NULL
  colnames(intersect_df) <- c("chr", "pos", "seg_start", "seg_end", "maj_cn", "min_cn", "msai", "allelic_imbalance", "VAF")
  intersect_df$tot_cn   <- intersect_df$maj_cn + intersect_df$min_cn
  
  region_pp_path <- list.files(path = paste0(battenberg_base_path, tumour, "/Battenberg/aggregated/"), pattern = paste0(reg_qc$solution, "_purity_ploidy.txt.gz"), full.names = T)
  region_pp_path <- region_pp_path[grep(region, region_pp_path)]
  pp <- read.table(region_pp_path, head = T, sep = "\t")
  
  intersect_df$ploidy <- pp$ploidy
  intersect_df$purity <- pp$purity
  
  intersect_df$mut_cpn <- (intersect_df$VAF * 1/intersect_df$purity) * ((intersect_df$purity * intersect_df$tot_cn) + 2 * (1-intersect_df$purity))
  intersect_df$VAF <- NULL
  intersect_df$purity <- NULL
  intersect_df$ploidy <- NULL
  
  #add this back onto the mutation_data table
  mutation_data$maj_cn <- NULL
  mutation_data$min_cn <- NULL
  mutation_data$tot_cn <- NULL
  mutation_data$seg_start <- NULL
  mutation_data$seg_end <- NULL
  mutation_data$mut_cpn <- NULL
  mutation_data <- unique(mutation_data)
  
  mutation_data <- left_join(mutation_data, intersect_df)
  
  if(nrow(mutation_data_all[mutation_data_all$sample == region, ]) > 0){
    
    mutation_data <- mutation_data[, c("sample", "chr", "pos", "hap1_id", "hap1_bases", "hap2_bases", "hap_id_SNP_pos", "phase","maj_cn", "min_cn", "tot_cn", "seg_start", "seg_end", "mut_cpn", "allelic_imbalance", "msai")]
    
    snp_data <- fread(paste0(snp_data_base_path, region, "_phased_muts_long.tsv"), head = T, sep = "\t", data.table = FALSE)
    snp_data$mutation_snp_id <- paste0(snp_data$chr_mut, "_", snp_data$pos_mut, "_", snp_data$pos_hetsnp)
    
    haplo_paths <- list.files(path = paste0(battenberg_base_path, tumour, "/Battenberg/aggregated/"), pattern = "_multisample_haplotypes_", full.names = T)
    haplo_paths <- haplo_paths[grep("RDS", haplo_paths)]
    
    BAF_data <- lapply(haplo_paths, function(x){
      gr <- readRDS(x)
    })
    
    BAF_data <- do.call(c, BAF_data)
    colnames(mcols(BAF_data)) <- sub("--..........", "", colnames(mcols(BAF_data)))

    # run the functions
    print("getting data in shape")
    data <- wrangle.data.fun(region, tumour, snp_data, mutation_data, BAF_data)

    #read in the ploidy data
    ploidy <- pp$ploidy
    
    print("phasing mutations")
    phased_data <- phasing.phun(data, tumour, region, ploidy, input_path)
    
    print("plotting phased mutations")
    phun.plots(phased_data, tumour, region, data, cn, input_path, genome)
    
    print("reading in and formatting mutation cluster info")
    cluster_mut_path <- paste0(cluster_mut_base_path, tumour, "/", tumour, ".SCoutput.CLEAN.tsv")
    cluster_muts <- read.table(cluster_mut_path, head = T, sep = "\t")
    cluster_muts$CHR <- paste0("chr", cluster_muts$CHR)
    cluster_muts$SAMPLE <- sub("\\.", "-", cluster_muts$SAMPLE)
    
    phased_data <- left_join(phased_data, cluster_muts[, c("SAMPLE", "CHR", "POS", "CLUSTER", "CCF_PHYLO")], by = c("sample" = "SAMPLE", "chr" = "CHR", "pos" = "POS"))
    
    
    print("reading in tree object")
    tree_path <- paste0(tree_base_path, tumour, "/", tumour, ".tree.RDS")
    tree_object <- readRDS(tree_path)
    
    phased_timed_data <- segment.cluster.timing.fun(phased_data, tree_object, output_path)
    
    
    # read in SVs to be able to plot them 
    SVs <- get.gridds.SVs(tumour, region, battenberg_base_path)
    SVs <- SVs[SVs$chr1 %in% paste0("chr", 1:22), ]
    SVs <- SVs[SVs$chr2 %in% paste0("chr", 1:22), ]
    SVs[is.na(SVs$SVLEN), "SVLEN"] <- 10001
    SVs <- SVs[SVs$SVLEN >= 10000, ]
    
    print("plotting gain timing")

    plot.timing.fun(timing_df = phased_timed_data, mutation_df = phased_data, cn_df = cn, all_mut = data, output_path = input_path, genome = BSgenome.Hsapiens.UCSC.hg38, SV_data = SVs)

    
    #plot tree plots for each timable segment

    
    # for plotting of copy number we need to get all segments and sort the haplotypes
    
    
    # 
    # timing.tree.fun(tumour_id = tumour, data = phased_timed_data, fullTree = tree_object, output_path, region)
    # 
  
    #all_timed_data[[region]] <- phased_timed_data
    
    # print("estimating mutation timing")
    # phased_data <- mut.cpn.match.fun(phased_data)
    
    # print("plotting early and late mutations")
    # plot.early.late.fun(phased_data, genome = BSgenome.Hsapiens.UCSC.hg38, output_path = input_path)

    # print("plotting phased mutations by CN state")
    # copynumber.state.plot.fun(phased_data, opacity = 255, tumour, region, input_path, genome)
 
    # print("calculating allele specific timing estimates")
    # phased_timed_data <- calc.timing.fun(phased_data, input_path)
    

    # # cn_states_data <- compare.cn.states.fun(phased_data, region)
    
    # all_timed_data[[region]] <- phased_timed_data
    
  }
  
}

all_timed_data <- do.call(rbind, all_timed_data)

write.table(all_timed_data, paste0(output_path, "plots/timed_data_cluster_test.txt"), col.names = T, row.names = F, sep = "\t", quote = F)

all_timed_data <- list()
for(region in all_regions){

   tumour = strsplit(region, "_")[[1]][1]

   filepath <- paste0(output_path, tumour, "/", region, "/data/clonal_gain_timing.txt")

   if(file.exists(filepath)){
     df <- read.table(filepath, head = T, sep = "\t")
     all_timed_data[[region]] <- df
   }
}

all_timed_data <- do.call(rbind, all_timed_data)

# check if timing estimates are independent of segment lenght and numbers of mutations

# segment length
p_seglength <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "seg_length", "seg_timing", "sample")])), x = "seg_length", y = "seg_timing", color = "sample",
                         add = "reg.line",
                         conf.int = TRUE,
                         add.params = list(color = "blue",
                                           fill = "lightgray")) +
                    xlab("segment length") +
                    ylab("segment timing estimate") +
                    stat_cor(method = "pearson", label.x = 3, label.y = 1) +
                    theme(legend.position = "right")

all_timed_data$seg_length_log <- log10(all_timed_data$seg_length)
p_seglength_log <- ggscatter(na.omit(unique(all_timed_data[ , c("segID", "seg_length_log", "seg_timing", "sample")])), x = "seg_length_log", y = "seg_timing", color = "sample",
                         add = "reg.line",
                         conf.int = TRUE,
                         add.params = list(color = "blue",
                                           fill = "lightgray")) +
  xlab("segment length (log10)") +
  ylab("segment timing estimate") +
  stat_cor(method = "pearson", label.x = 5, label.y = 1) +
  theme(legend.position = "right")

all_timed_data$total_mut_log <- log10(all_timed_data$total_mut)
p_mut <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "total_mut", "seg_timing", "sample")])), x = "total_mut", y = "seg_timing", color = "sample",
                         add = "reg.line",
                         conf.int = TRUE,
                         add.params = list(color = "blue",
                                           fill = "lightgray")) +
  xlab("# muts on segment ") +
  ylab("segment timing estimate") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  theme(legend.position = "right")

all_timed_data$total_mut_log <- log10(all_timed_data$total_mut)
p_mut_log <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "total_mut_log", "seg_timing", "sample")])), x = "total_mut_log", y = "seg_timing", color = "sample",
                   add = "reg.line",
                   conf.int = TRUE,
                   add.params = list(color = "blue",
                                     fill = "lightgray")) +
  xlab("# muts on segment (log10)") +
  ylab("segment timing estimate") +
  stat_cor(method = "pearson", label.x = 2, label.y = 1) +
  theme(legend.position = "right")

pdf(paste0(output_path, "unadjusted_segtiming_QC_plots.pdf"), width = 7, height = 5)
p_seglength
p_seglength_log
p_mut
p_mut_log
dev.off()


# make these plots for every region
p_seglength_tumour <- ggscatter(na.omit(unique(all_timed_data[, c("sample", "segID", "seg_length", "seg_timing", "sample")])), x = "seg_length", y = "seg_timing", color = "sample",
                                add = "reg.line",
                                conf.int = TRUE,
                                add.params = list(color = "blue",
                                                  fill = "lightgray")) +
  xlab("segment length") +
  ylab("segment timing estimate") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  facet_wrap(~sample) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

p_seglength_log_tumour <- ggscatter(na.omit(unique(all_timed_data[, c("sample", "segID", "seg_length_log", "seg_timing", "sample")])), x = "seg_length_log", y = "seg_timing", color = "sample",
                                add = "reg.line",
                                conf.int = TRUE,
                                add.params = list(color = "blue",
                                                  fill = "lightgray")) +
  xlab("segment length (log10)") +
  ylab("segment timing estimate") +
  stat_cor(method = "pearson", label.x = 5, label.y = 1) +
  facet_wrap(~sample) +
  theme(legend.position = "none")

p_mut_tumour <- ggscatter(na.omit(unique(all_timed_data[, c("sample", "segID", "total_mut", "seg_timing", "sample")])), x = "total_mut", y = "seg_timing", color = "sample",
                   add = "reg.line",
                   conf.int = TRUE,
                   add.params = list(color = "blue",
                                     fill = "lightgray")) +
  xlab("# muts on segment ") +
  ylab("segment timing estimate") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  facet_wrap(~sample) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

p_mut_log_tumour <- ggscatter(na.omit(unique(all_timed_data[, c("sample", "segID", "total_mut_log", "seg_timing", "sample")])), x = "total_mut_log", y = "seg_timing", color = "sample",
                          add = "reg.line",
                          conf.int = TRUE,
                          add.params = list(color = "blue",
                                            fill = "lightgray")) +
  xlab("# muts on segment (log10)") +
  ylab("segment timing estimate") +
  stat_cor(method = "pearson", label.x = 1, label.y = 1) +
  facet_wrap(~sample) +
  theme(legend.position = "none")

pdf(paste0(output_path, "unadjusted_segtiming_QC_per_tumour_plots.pdf"), width = 7, height = 7)
p_seglength_tumour
p_seglength_log_tumour
p_mut_tumour
p_mut_log_tumour
dev.off()

# there is a relationship between length of segment and number of early and late muts
p_early <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "early", "seg_timing", "sample")])), x = "early", y = "seg_timing", color = "sample",
                     add = "reg.line",
                     conf.int = TRUE,
                     add.params = list(color = "blue",
                                       fill = "lightgray")) +
  xlab("# early muts per segment") +
  ylab("segment timing estimate") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  theme(legend.position = "right")

all_timed_data$early_log <- log10(all_timed_data$early)
p_early_log <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "early_log", "seg_timing", "sample")])), x = "early_log", y = "seg_timing", color = "sample",
                     add = "reg.line",
                     conf.int = TRUE,
                     add.params = list(color = "blue",
                                       fill = "lightgray")) +
  xlab("# early muts per segment (log10)") +
  ylab("segment timing estimate") +
  stat_cor(method = "pearson", label.x = 1, label.y = 1) +
  theme(legend.position = "right")


p_late <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "late", "seg_timing", "sample")])), x = "late", y = "seg_timing", color = "sample",
                     add = "reg.line",
                     conf.int = TRUE,
                     add.params = list(color = "blue",
                                       fill = "lightgray")) +
  xlab("# late muts per segment") +
  ylab("segment timing estimate") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  theme(legend.position = "right")

all_timed_data$late_log <- log10(all_timed_data$late)
p_late_log <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "late_log", "seg_timing", "sample")])), x = "late_log", y = "seg_timing", color = "sample",
                    add = "reg.line",
                    conf.int = TRUE,
                    add.params = list(color = "blue",
                                      fill = "lightgray")) +
  xlab("# late muts per segment (log10)") +
  ylab("segment timing estimate") +
  stat_cor(method = "pearson", label.x = 1, label.y = 1) +
  theme(legend.position = "right")

p_early_length <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "early", "seg_length", "sample")])), y = "early", x = "seg_length", color = "sample",
                     add = "reg.line",
                     conf.int = TRUE,
                     add.params = list(color = "blue",
                                       fill = "lightgray")) +
  ylab("# early muts per segment") +
  xlab("segment length") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  theme(legend.position = "right")

p_early_length_log <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "early_log", "seg_length_log", "sample")])), y = "early_log", x = "seg_length_log", color = "sample",
                            add = "reg.line",
                            conf.int = TRUE,
                            add.params = list(color = "blue",
                                              fill = "lightgray")) +
  ylab("# early muts per segment (log10)") +
  xlab("segment length (log10)") +
  stat_cor(method = "pearson", label.x = 5, label.y = 3) +
  theme(legend.position = "right")

p_late_length <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "late", "seg_length", "sample")])), y = "late", x = "seg_length", color = "sample",
                    add = "reg.line",
                    conf.int = TRUE,
                    add.params = list(color = "blue",
                                      fill = "lightgray")) +
  ylab("# late muts per segment") +
  xlab("segment length") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  theme(legend.position = "right")

p_late_length_log <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "late_log", "seg_length_log", "sample")])), y = "late_log", x = "seg_length_log", color = "sample",
                           add = "reg.line",
                           conf.int = TRUE,
                           add.params = list(color = "blue",
                                             fill = "lightgray")) +
  ylab("# late muts per segment (log10)") +
  xlab("segment length (log10)") +
  stat_cor(method = "pearson", label.x = 5, label.y = 3) +
  theme(legend.position = "right")

all_timed_data$early_length_ratio <- all_timed_data$early / all_timed_data$seg_length
all_timed_data$late_length_ratio <- all_timed_data$late / all_timed_data$seg_length

p_early_length_ratio <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "early_length_ratio", "seg_length", "sample")])), y = "early_length_ratio", x = "seg_length", color = "sample",
                            add = "reg.line",
                            conf.int = TRUE,
                            add.params = list(color = "blue",
                                              fill = "lightgray")) +
  ylab("# early muts per segment / segment length") +
  xlab("segment length") +
  stat_cor(method = "pearson", label.x = 3, label.y = 0) +
  theme(legend.position = "right")

p_late_length_ratio <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "late_length_ratio", "seg_length", "sample")])), y = "late_length_ratio", x = "seg_length", color = "sample",
                                  add = "reg.line",
                                  conf.int = TRUE,
                                  add.params = list(color = "blue",
                                                    fill = "lightgray")) +
  ylab("# late muts per segment / segment length") +
  xlab("segment length") +
  stat_cor(method = "pearson", label.x = 3, label.y = 0) +
  theme(legend.position = "right")

pdf(paste0(output_path, "early_late_relationship_seglenght.pdf"), width = 7, height = 5)
p_early
p_early_log
p_late
p_late_log
p_early_length
p_early_length_log
p_late_length
p_late_length_log
p_early_length_ratio
p_late_length_ratio
dev.off()


# check if there is a relationship with copy number and timing estimate
cn_timing_df <- unique(all_timed_data[, c("segID", "maj_cn", "min_cn", "phase", "seg_timing", "sample")])
p_cn_timing <- ggscatter(cn_timing_df[cn_timing_df$phase == "A", ], x = "maj_cn", y = "seg_timing", color = "sample",
                                 add = "reg.line",
                                 conf.int = TRUE,
                                 add.params = list(color = "blue",
                                                   fill = "lightgray")) +
  ylab("segment timing") +
  xlab("haplotype A copy number") +
  stat_cor(method = "pearson", label.y = 1, label.x = 40) +
  theme(legend.position = "right")

p_cn_timing_cut <- ggscatter(cn_timing_df[cn_timing_df$phase == "A" & cn_timing_df$maj_cn <= 10, ], x = "maj_cn", y = "seg_timing", color = "sample",
                         add = "reg.line",
                         conf.int = TRUE,
                         add.params = list(color = "blue",
                                           fill = "lightgray")) +
  ylab("segment timing") +
  xlab("haplotype A copy number") +
  stat_cor(method = "pearson", label.y = 1, label.x = 5) +
  theme(legend.position = "right")

p_cn_timingB <- ggscatter(cn_timing_df[cn_timing_df$phase == "B", ], x = "min_cn", y = "seg_timing", color = "sample",
                         add = "reg.line",
                         conf.int = TRUE,
                         add.params = list(color = "blue",
                                           fill = "lightgray")) +
  ylab("segment timing") +
  xlab("haplotype B copy number") +
  stat_cor(method = "pearson", label.y = 1, label.x = 3) +
  theme(legend.position = "right")

pdf(paste0(output_path, "timing_copynumber_relationship.pdf"), width = 7, height = 5)
p_cn_timing
p_cn_timing_cut
p_cn_timingB
dev.off()

# # # # # # # # 

# make the same plots but with the copy number adjusted timing estimates

# # # # # # # # 

# segment length
p_seglength <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "seg_length", "seg_timing_adjusted", "sample")])), x = "seg_length", y = "seg_timing_adjusted", color = "sample",
                         add = "reg.line",
                         conf.int = TRUE,
                         add.params = list(color = "blue",
                                           fill = "lightgray")) +
  xlab("segment length") +
  ylab("adjusted segment timing estimate") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  theme(legend.position = "right")

all_timed_data$seg_length_log <- log10(all_timed_data$seg_length)
p_seglength_log <- ggscatter(na.omit(unique(all_timed_data[ , c("segID", "seg_length_log", "seg_timing_adjusted", "sample")])), x = "seg_length_log", y = "seg_timing_adjusted", color = "sample",
                             add = "reg.line",
                             conf.int = TRUE,
                             add.params = list(color = "blue",
                                               fill = "lightgray")) +
  xlab("segment length (log10)") +
  ylab("adjusted segment timing estimate") +
  stat_cor(method = "pearson", label.x = 5, label.y = 1) +
  theme(legend.position = "right")

# remove the outlier
p_seglength_log_o <- ggscatter(na.omit(unique(all_timed_data[ , c("segID", "seg_length_log", "seg_timing_adjusted", "sample")])), x = "seg_length_log", y = "seg_timing_adjusted", color = "sample",
                             add = "reg.line",
                             conf.int = TRUE,
                             add.params = list(color = "blue",
                                               fill = "lightgray")) +
  xlab("segment length (log10)") +
  ylab("adjusted segment timing estimate") +
  stat_cor(method = "pearson", label.x = 5, label.y = 1) +
  theme(legend.position = "right")


all_timed_data$total_mut_log <- log10(all_timed_data$total_mut)
p_mut <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "total_mut", "seg_timing_adjusted", "sample")])), x = "total_mut", y = "seg_timing_adjusted", color = "sample",
                   add = "reg.line",
                   conf.int = TRUE,
                   add.params = list(color = "blue",
                                     fill = "lightgray")) +
  xlab("# muts on segment ") +
  ylab("adjusted segment timing estimate") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  theme(legend.position = "right")

all_timed_data$total_mut_log <- log10(all_timed_data$total_mut)
p_mut_log <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "total_mut_log", "seg_timing_adjusted", "sample")])), x = "total_mut_log", y = "seg_timing_adjusted", color = "sample",
                       add = "reg.line",
                       conf.int = TRUE,
                       add.params = list(color = "blue",
                                         fill = "lightgray")) +
  xlab("# muts on segment (log10)") +
  ylab("adjusted segment timing estimate") +
  stat_cor(method = "pearson", label.x = 2, label.y = 1) +
  theme(legend.position = "right")

pdf(paste0(output_path, "adjusted_segtiming_QC_plots.pdf"), width = 7, height = 5)
p_seglength
p_seglength_log
p_seglength_log_o
p_mut
p_mut_log
dev.off()

# make these plots for every region
p_seglength_tumour <- ggscatter(na.omit(unique(all_timed_data[, c("sample", "segID", "seg_length", "seg_timing_adjusted")])), x = "seg_length", y = "seg_timing_adjusted", color = "sample",
                                add = "reg.line",
                                conf.int = TRUE,
                                add.params = list(color = "blue",
                                                  fill = "lightgray")) +
  xlab("segment length") +
  ylab("adjusted segment timing estimate") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  facet_wrap(~sample) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

p_seglength_log_tumour <- ggscatter(na.omit(unique(all_timed_data[, c("sample", "segID", "seg_length_log", "seg_timing_adjusted")])), x = "seg_length_log", y = "seg_timing_adjusted", color = "sample",
                                    add = "reg.line",
                                    conf.int = TRUE,
                                    add.params = list(color = "blue",
                                                      fill = "lightgray")) +
  xlab("segment length (log10)") +
  ylab("adjusted segment timing estimate") +
  stat_cor(method = "pearson", label.x = 5, label.y = 1) +
  facet_wrap(~sample) +
  theme(legend.position = "none")

p_mut_tumour <- ggscatter(na.omit(unique(all_timed_data[, c("sample", "segID", "total_mut", "seg_timing_adjusted")])), x = "total_mut", y = "seg_timing_adjusted", color = "sample",
                          add = "reg.line",
                          conf.int = TRUE,
                          add.params = list(color = "blue",
                                            fill = "lightgray")) +
  xlab("# muts on segment ") +
  ylab("adjusted segment timing estimate") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  facet_wrap(~sample) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

p_mut_log_tumour <- ggscatter(na.omit(unique(all_timed_data[, c("sample", "segID", "total_mut_log", "seg_timing_adjusted")])), x = "total_mut_log", y = "seg_timing_adjusted", color = "sample",
                              add = "reg.line",
                              conf.int = TRUE,
                              add.params = list(color = "blue",
                                                fill = "lightgray")) +
  xlab("# muts on segment (log10)") +
  ylab("adjusted segment timing estimate") +
  stat_cor(method = "pearson", label.x = 1, label.y = 1) +
  facet_wrap(~sample) +
  theme(legend.position = "none")

pdf(paste0(output_path, "adjusted_segtiming_QC_per_tumour_plots.pdf"), width = 7, height = 7)
p_seglength_tumour
p_seglength_log_tumour
p_mut_tumour
p_mut_log_tumour
dev.off()

# there is a relationship between length of segment and number of early and late muts
p_early <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "early", "seg_timing_adjusted", "sample")])), x = "early", y = "seg_timing_adjusted", color = "sample",
                     add = "reg.line",
                     conf.int = TRUE,
                     add.params = list(color = "blue",
                                       fill = "lightgray")) +
  xlab("# early muts per segment") +
  ylab("adjusted segment timing estimate") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  theme(legend.position = "right")

all_timed_data$early_log <- log10(all_timed_data$early)
p_early_log <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "early_log", "seg_timing_adjusted", "sample")])), x = "early_log", y = "seg_timing_adjusted", color = "sample",
                         add = "reg.line",
                         conf.int = TRUE,
                         add.params = list(color = "blue",
                                           fill = "lightgray")) +
  xlab("# early muts per segment (log10)") +
  ylab("adjusted segment timing estimate") +
  stat_cor(method = "pearson", label.x = 1, label.y = 1) +
  theme(legend.position = "right")


p_late <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "late", "seg_timing_adjusted", "sample")])), x = "late", y = "seg_timing_adjusted", color = "sample",
                    add = "reg.line",
                    conf.int = TRUE,
                    add.params = list(color = "blue",
                                      fill = "lightgray")) +
  xlab("# late muts per segment") +
  ylab("adjusted segment timing estimate") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  theme(legend.position = "right")

all_timed_data$late_log <- log10(all_timed_data$late)
p_late_log <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "late_log", "seg_timing_adjusted", "sample")])), x = "late_log", y = "seg_timing_adjusted", color = "sample",
                        add = "reg.line",
                        conf.int = TRUE,
                        add.params = list(color = "blue",
                                          fill = "lightgray")) +
  xlab("# late muts per segment (log10)") +
  ylab("adjusted segment timing estimate") +
  stat_cor(method = "pearson", label.x = 1, label.y = 1) +
  theme(legend.position = "right")


p_late_adjusted <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "late_adjusted", "seg_timing_adjusted", "sample")])), x = "late_adjusted", y = "seg_timing_adjusted", color = "sample",
                    add = "reg.line",
                    conf.int = TRUE,
                    add.params = list(color = "blue",
                                      fill = "lightgray")) +
  xlab("# late adjusted muts per segment") +
  ylab("adjusted segment timing estimate") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  theme(legend.position = "right")

all_timed_data$late_adjusted_log <- log10(all_timed_data$late_adjusted)
p_late_adjusted_log <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "late_adjusted_log", "seg_timing_adjusted", "sample")])), x = "late_adjusted_log", y = "seg_timing_adjusted", color = "sample",
                        add = "reg.line",
                        conf.int = TRUE,
                        add.params = list(color = "blue",
                                          fill = "lightgray")) +
  xlab("# late adjusted  muts per segment (log10)") +
  ylab("adjusted segment timing estimate") +
  stat_cor(method = "pearson", label.x = 1, label.y = 1) +
  theme(legend.position = "right")


p_early_length <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "early", "seg_length", "sample")])), y = "early", x = "seg_length", color = "sample",
                            add = "reg.line",
                            conf.int = TRUE,
                            add.params = list(color = "blue",
                                              fill = "lightgray")) +
  ylab("# early muts per segment") +
  xlab("segment length") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  theme(legend.position = "right")

p_early_length_log <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "early_log", "seg_length_log", "sample")])), y = "early_log", x = "seg_length_log", color = "sample",
                                add = "reg.line",
                                conf.int = TRUE,
                                add.params = list(color = "blue",
                                                  fill = "lightgray")) +
  ylab("# early muts per segment (log10)") +
  xlab("segment length (log10)") +
  stat_cor(method = "pearson", label.x = 5, label.y = 3) +
  theme(legend.position = "right")

p_late_length <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "late", "seg_length", "sample")])), y = "late", x = "seg_length", color = "sample",
                           add = "reg.line",
                           conf.int = TRUE,
                           add.params = list(color = "blue",
                                             fill = "lightgray")) +
  ylab("# late muts per segment") +
  xlab("segment length") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  theme(legend.position = "right")

p_late_length_log <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "late_log", "seg_length_log", "sample")])), y = "late_log", x = "seg_length_log", color = "sample",
                               add = "reg.line",
                               conf.int = TRUE,
                               add.params = list(color = "blue",
                                                 fill = "lightgray")) +
  ylab("# late muts per segment (log10)") +
  xlab("segment length (log10)") +
  stat_cor(method = "pearson", label.x = 5, label.y = 3) +
  theme(legend.position = "right")

p_late_adjusted_length <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "late_adjusted", "seg_length", "sample")])), y = "late_adjusted", x = "seg_length", color = "sample",
                           add = "reg.line",
                           conf.int = TRUE,
                           add.params = list(color = "blue",
                                             fill = "lightgray")) +
  ylab("# late adjusted muts per segment") +
  xlab("segment length") +
  stat_cor(method = "pearson", label.x = 3, label.y = 1) +
  theme(legend.position = "right")

p_late_adjusted_length_log <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "late_adjusted_log", "seg_length_log", "sample")])), y = "late_adjusted_log", x = "seg_length_log", color = "sample",
                               add = "reg.line",
                               conf.int = TRUE,
                               add.params = list(color = "blue",
                                                 fill = "lightgray")) +
  ylab("# late adjusted muts per segment (log10)") +
  xlab("segment length (log10)") +
  stat_cor(method = "pearson", label.x = 5, label.y = 3) +
  theme(legend.position = "right")

all_timed_data$early_length_ratio <- all_timed_data$early / all_timed_data$seg_length
all_timed_data$late_length_ratio <- all_timed_data$late / all_timed_data$seg_length

p_early_length_ratio <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "early_length_ratio", "seg_length", "sample")])), y = "early_length_ratio", x = "seg_length", color = "sample",
                                  add = "reg.line",
                                  conf.int = TRUE,
                                  add.params = list(color = "blue",
                                                    fill = "lightgray")) +
  ylab("# early muts per segment / segment length") +
  xlab("segment length") +
  stat_cor(method = "pearson", label.x = 3, label.y = 0) +
  theme(legend.position = "right")

p_late_length_ratio <- ggscatter(na.omit(unique(all_timed_data[, c("segID", "late_length_ratio", "seg_length", "sample")])), y = "late_length_ratio", x = "seg_length", color = "sample",
                                 add = "reg.line",
                                 conf.int = TRUE,
                                 add.params = list(color = "blue",
                                                   fill = "lightgray")) +
  ylab("# late muts per segment / segment length") +
  xlab("segment length") +
  stat_cor(method = "pearson", label.x = 3, label.y = 0) +
  theme(legend.position = "right")

pdf(paste0(output_path, "adjusted_early_late_relationship_seglenght.pdf"), width = 7, height = 5)
p_early
p_early_log
p_late
p_late_log
p_late_adjusted
p_late_adjusted_log
p_early_length
p_early_length_log
p_late_length
p_late_length_log
p_late_adjusted_length
p_late_adjusted_length_log
p_early_length_ratio
p_late_length_ratio
dev.off()

# check if there is a relationship with copy number and timing estimate
cn_timing_df <- unique(all_timed_data[, c("segID", "maj_cn", "min_cn", "phase", "seg_timing_adjusted", "sample")])
p_cn_timing <- ggscatter(cn_timing_df[cn_timing_df$phase == "A", ], x = "maj_cn", y = "seg_timing_adjusted", color = "sample",
                         add = "reg.line",
                         conf.int = TRUE,
                         add.params = list(color = "blue",
                                           fill = "lightgray")) +
  ylab("adjusted segment timing") +
  xlab("haplotype A copy number") +
  stat_cor(method = "pearson", label.y = 1, label.x = 40) +
  theme(legend.position = "right")

p_cn_timing_cut <- ggscatter(cn_timing_df[cn_timing_df$phase == "A" & cn_timing_df$maj_cn <= 10, ], x = "maj_cn", y = "seg_timing_adjusted", color = "sample",
                             add = "reg.line",
                             conf.int = TRUE,
                             add.params = list(color = "blue",
                                               fill = "lightgray")) +
  ylab("adjusted segment timing") +
  xlab("haplotype A copy number") +
  stat_cor(method = "pearson", label.y = 1, label.x = 5) +
  theme(legend.position = "right")

p_cn_timingB <- ggscatter(cn_timing_df[cn_timing_df$phase == "B", ], x = "min_cn", y = "seg_timing_adjusted", color = "sample",
                          add = "reg.line",
                          conf.int = TRUE,
                          add.params = list(color = "blue",
                                            fill = "lightgray")) +
  ylab("adjusted segment timing") +
  xlab("haplotype B copy number") +
  stat_cor(method = "pearson", label.y = 1, label.x = 3) +
  theme(legend.position = "right")

pdf(paste0(output_path, "adjusted_timing_copynumber_relationship.pdf"), width = 7, height = 5)
p_cn_timing
p_cn_timing_cut
p_cn_timingB
dev.off()


#density plot
phased_data$segID <- paste0(phased_data$chr, "_", phased_data$seg_start, "_", phased_data$seg_end)
unique(phased_data[phased_data$chr == "chr12", "segID"])

dat <- phased_data[phased_data$segID == "chr12_72815038_133254685" & phased_data$phase == "A", ]

p <- ggplot(dat, aes(x = mut_cpn)) +
        geom_density() +
        theme_bw()

pdf(paste0(output_path, "QC/density_example.pdf"), width = 4, height = 4)
p
dev.off()
