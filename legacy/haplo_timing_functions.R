# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #           functions to phase mutations and copy number                # # #
# # #       using long read and short read whole genome sequencing          # # #
# # #                 and to perform timing analysis                        # # #  
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
library(patchwork)
library(igraph)
library(GenomicRanges)
library(IRanges)
library(ggnewscale)

################################################################################
##################################################################### file paths


# genome                     <- BSgenome.Hsapiens.UCSC.hg38

################################################################################
#################################################################### functions

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

# function to make ascat like copy number data frame continuous across several chromosomes
# genome.continous.cn.fun <- function(data, genome){
#   
#   chr_sizes <- seqlengths(genome)[paste0("chr", 1:22)]
#   chr_sizes <- data.frame(chr_sizes)
#   chr_sizes$chr <- c(1:22)
#   chr_sizes <- chr_sizes[order(chr_sizes$chr, decreasing = F), ]
#   chr_sizes <- chr_sizes$chr_sizes
#   
#   # adjust pos to continous value
#   cumsum_chr_size <- cumsum(as.numeric(chr_sizes))
#   
#   
#   data$chr_adjust <- data$chr
#   data$chr_adjust <- sub("chr", "", data$chr_adjust)
#   # also adjust seg start
#   data$adjust_seg_start <- data$seg_start
#   for(i in c(2:22)){
#     data$adjust_seg_start[data$chr_adjust == i] <- data$adjust_seg_start[data$chr_adjust == i] + cumsum_chr_size[i-1]
#   }
#   
#   # and seg end
#   data$adjust_seg_end <- data$seg_end
#   for(i in c(2:22)){
#     data$adjust_seg_end[data$chr_adjust == i] <- data$adjust_seg_end[data$chr_adjust == i] + cumsum_chr_size[i-1]
#   }
#   
#   return(data)
# }


# function to annotate parallel mutations
parallel.mut.fun <- function(snp_data){
  
  snp_data_par <- snp_data[, c("mutid", "SNV_MT_SNP_MT_N", "SNV_WT_SNP_MT_N", "SNV_MT_SNP_WT_N", "SNV_WT_SNP_WT_N")]
  snp_data_par <- snp_data_par %>% 
                    group_by(mutid) %>% 
                    mutate(sum_SNV_MT_SNP_MT_N = sum(SNV_MT_SNP_MT_N)) %>%
                    mutate(sum_SNV_WT_SNP_MT_N = sum(SNV_WT_SNP_MT_N)) %>%
                    mutate(sum_SNV_MT_SNP_WT_N = sum(SNV_MT_SNP_WT_N)) %>%
                    mutate(sum_SNV_WT_SNP_WT_N = sum(SNV_WT_SNP_WT_N))
  
  snp_data_par <- unique(snp_data_par[, c("mutid", "sum_SNV_MT_SNP_MT_N", "sum_SNV_WT_SNP_MT_N", "sum_SNV_MT_SNP_WT_N", "sum_SNV_WT_SNP_WT_N")])
    
  snp_data_par[snp_data_par$sum_SNV_MT_SNP_MT_N > 0 & snp_data_par$sum_SNV_MT_SNP_WT_N > 0, "parallel_mutation"] <- TRUE
  snp_data_par[is.na(snp_data_par$parallel_mutation), "parallel_mutation"] <- FALSE
  
  snp_data <- left_join(snp_data, snp_data_par)
  
  return(snp_data)
  
}

# function to get data ready for phasing
wrangle.data.fun <- function(region, tumour, snp_data, mutation_data, BAF_data){
  
  # let's get this in long format to get all the separate SNPs that are linked to a mutation
  mutation_data$mutation_id <- paste0(mutation_data$chr, "_", mutation_data$pos, ";", mutation_data$hap_id_SNP_pos) 
  
  mutation_data_hap1_id <- mutation_data[, c("hap_id_SNP_pos", "hap1_id", "mutation_id")]
  mutation_data_hap1_id <- mutation_data_hap1_id %>% 
  separate_rows(hap_id_SNP_pos, hap1_id, sep = ":")
  
  mutation_data_hap1_bases <- mutation_data[, c("hap_id_SNP_pos", "hap1_bases", "mutation_id")]
  mutation_data_hap1_bases <- mutation_data_hap1_bases %>% 
  separate_rows(hap_id_SNP_pos, hap1_bases, sep = ":")
  
  mutation_data_hap2_bases <- mutation_data[, c("hap_id_SNP_pos", "hap2_bases", "mutation_id")]
  mutation_data_hap2_bases <- mutation_data_hap2_bases %>% 
  separate_rows(hap_id_SNP_pos, hap2_bases, sep = ":")
  
  mutation_data_long <- cbind(mutation_data_hap1_id, mutation_data_hap1_bases)
  mutation_data_long$hap_id_SNP_pos <- NULL  
  mutation_data_long$mutation_id <- NULL 
  mutation_data_long <- cbind(mutation_data_long, mutation_data_hap2_bases)
  mutation_data_long$hap_id_SNP_pos <- NULL  
  mutation_data_long$mutation_id <- NULL 
  
  mutation_data_long <- mutation_data_long[, c("mutation_id", "hap_id_SNP_pos", "hap1_id", "hap1_bases", "hap2_bases")]
  mutation_data_long <- unique(mutation_data_long)
  mutation_data      <- unique(mutation_data)
  mutation_data_long <- left_join(mutation_data_long, mutation_data[, c("sample", "chr", "pos", "maj_cn", "min_cn", "tot_cn", "seg_start", "seg_end", "mut_cpn", "phase","mutation_id", "allelic_imbalance", "msai")])
  mutation_data_long$phase <- NULL
  
  mutation_data_long$chr <- paste0("chr", mutation_data_long$chr)
  mutation_data_long$mutation_snp_id <- paste0(mutation_data_long$chr, "_", mutation_data_long$pos, "_", mutation_data_long$hap_id_SNP_pos)
  
  snp_data <-  parallel.mut.fun(snp_data)
  
  snp_data_small <- snp_data[, c("mutation_snp_id", "ref_mut", "alt_mut", "counts_mut", "chr_hetsnp", "pos_hetsnp", "ref_hetsnp", "alt_hetsnp", "counts_hetsnp", 
                               "sum_SNV_MT_SNP_MT_N", "sum_SNV_MT_SNP_WT_N", "sum_SNV_WT_SNP_MT_N", "sum_SNV_WT_SNP_WT_N",  "parallel_mutation", "snp_vaf_on_mut")]
  
  snp_data_small <- unique(snp_data_small)
  
  snp_mutation_data <- left_join(mutation_data_long, snp_data_small)
  snp_mutation_data$snp_id <- paste0(snp_mutation_data$chr, "_", snp_mutation_data$pos_hetsnp)
  
  # add BAF data
  col_Major <- paste0(region, "_Major")
  col_BAF <- paste0(region, "_BAF")
  col_PS <-  paste0(region, "_PS")
  
  BAF_region <- BAF_data[, c("REF", "ALT", paste0(col_Major), paste0(col_BAF), paste0(col_PS), "multisample_haplo")]
  BAF_region$REF <- as.character(BAF_region$REF)
  BAF_region$ALT <- as.character(unlist(BAF_region$ALT))
  
  BAF_region <- data.frame(BAF_region)
  BAF_region$width <- NULL
  BAF_region$strand <- NULL
  
  # add this to the other dataframe
  BAF_region$seqnames <- paste0("chr", BAF_region$seqnames)
  BAF_region$snp_id <- paste0(BAF_region$seqnames, "_", BAF_region$start)
  colnames(BAF_region) <- c("bb_snp_chr", "bb_snp_start", "bb_snp_end", "bb_snp_REF", "bb_SNP_ALT", "bb_Major", "bb_BAF", "bb_PS", "bb_mutisample_haplo", "snp_id")

  BAF_region$bb_snp_chr <- NULL
  BAF_region$bb_snp_end <- NULL
  BAF_region$bb_snp_start <- NULL

  snp_mutation_data_chr <- left_join(snp_mutation_data, BAF_region)
  snp_mutation_df <- snp_mutation_data_chr[!is.na(snp_mutation_data_chr$bb_BAF), ]
  
  snp_mutation_df <- snp_mutation_df[!is.na(snp_mutation_df$snp_vaf_on_mut), ]

  # #we have to properly assign major and minor CN to haplotypes
  # #separate the regions
  # msai_long$region <- sub("_Major", "", msai_long$variable)
  # msai_long$region <- sub("\\.", "-", msai_long$region)
  # msai_long <- msai_long[msai_long$region == region, ]
  # msai_long$seqnames <- paste0("chr", msai_long$seqnames)
  # 
  # #map the start and end positions to the region specific start and end positions
  # # for subsequent analysis let's use minimum consistent regions, might have to have a single regions and a multi region mode
  # segment_map$seqnames <- paste0("chr", segment_map$seqnames)
  # 
  # snp_mutation_df$seg_start <- NULL
  # snp_mutation_df$seg_end <- NULL
  # 
  # #do an overlap for the muts with the minimum consisten regions
  # snp_mutation_gr <- GRanges(seqnames = snp_mutation_df$chr, IRanges(start = snp_mutation_df$pos, end = snp_mutation_df$pos))
  # segment_map_gr <- GRanges(seqnames = segment_map$seqnames, IRanges(start = segment_map$regionCN_seg_start, end = segment_map$regionCN_seg_end))
  # 
  # overlap         <-  findOverlaps(snp_mutation_gr, segment_map_gr)
  # snp_mutation_df <-  cbind(snp_mutation_df[queryHits(overlap), ], segment_map[subjectHits(overlap), ])
  # 
  # colnames(snp_mutation_df)[(ncol(snp_mutation_df) -1) : ncol(snp_mutation_df)] <- c("seg_start", "seg_end")
  # 
  # snp_mutation_df <- unique(snp_mutation_df)
  # 
  # snp_mutation_df <- left_join(snp_mutation_df, msai_long[, c("seqnames", "start", "end", "msai", "flip", "allelic_imbalance")], by = c("chr" = "seqnames",
  #                                                                                                                                       "seg_start" = "start",
  #                                                                                                                                       "seg_end" = "end"))
  # snp_mutation_df <- snp_mutation_df[!is.na(snp_mutation_df$maj_cn), ]
  # snp_mutation_df <- snp_mutation_df[!is.na(snp_mutation_df$min_cn), ]
  # 
  # #now using the flip collumn we can assign copy number to A and B instead major and minor based on being the greater value
  # snp_mutation_df[which(snp_mutation_df$flip == FALSE), "cnA"] <- snp_mutation_df[which(snp_mutation_df$flip == FALSE), "maj_cn"]
  # snp_mutation_df[which(snp_mutation_df$flip == FALSE), "cnB"] <- snp_mutation_df[which(snp_mutation_df$flip == FALSE), "min_cn"]
  # snp_mutation_df[which(snp_mutation_df$flip == TRUE), "cnA"] <- snp_mutation_df[which(snp_mutation_df$flip == TRUE), "min_cn"]
  # snp_mutation_df[which(snp_mutation_df$flip == TRUE), "cnB"] <- snp_mutation_df[which(snp_mutation_df$flip == TRUE), "maj_cn"]
  # 
  # #lets remove the old non-flipped copy number
  # snp_mutation_df$maj_cn <- NULL
  # snp_mutation_df$min_cn <- NULL
  # 
  # # and rename into major and minor, this makes no sense but prevents me having to change the code everywhere
  # 
  # colnames(snp_mutation_df)[(ncol(snp_mutation_df) -1) : ncol(snp_mutation_df)] <- c("maj_cn", "min_cn")

  return(snp_mutation_df)
}       

#function to properly assign multi-region haplotypes to copy number
assign.haplo.cn.fun <- function(cn_data = cn, msai_long){
  
  #we also need to create a map from region specific CN segments to minimum consistent region segments
  mcr <- unique(msai_long[, c("seqnames", "start", "end")])
  colnames(mcr) <- c("chr", "seg_start", "seg_end")
  mcr$mcr <- TRUE
  mcr_gr <- GRanges(seqnames = mcr$chr, IRanges(start = mcr$seg_start, end = mcr$seg_end))
  mcols(mcr_gr) <- data.frame(mcr = TRUE, mcr_seg_start = mcr$seg_start, mcr_seg_end = mcr$seg_end)
  
  regionCN <- unique(cn_data[, c("chr", "startpos", "endpos")])
  colnames(regionCN) <- c("chr", "seg_start", "seg_end")
  regionCN <- regionCN[!is.na(regionCN$seg_start), ]
  regionCN <- unique(regionCN)
  regionCN$regionCN <- TRUE
  regionCN_gr <- GRanges(seqnames = regionCN$chr, IRanges(start = regionCN$seg_start, end = regionCN$seg_end))
  mcols(regionCN_gr) <- data.frame(regionCN = TRUE, regionCN_seg_start = regionCN$seg_start, regionCN_seg_end = regionCN$seg_end)
  
  overlap         <-  findOverlaps(mcr_gr, regionCN_gr)
  intersect_gr    <- pintersect(regionCN_gr[subjectHits(overlap)], mcr_gr[queryHits(overlap)])
  mcols(intersect_gr) <- data.frame(mcols(intersect_gr), mcols(mcr_gr[queryHits(overlap)]))
  intersect_df <- data.frame(intersect_gr)
  
  segment_map <- intersect_df[, c("seqnames", "regionCN_seg_start", "regionCN_seg_end", "mcr_seg_start", "mcr_seg_end")]
  
  #we have to properly assign major and minor CN to haplotypes
  #separate the regions
  msai_long$region <- sub("_Major", "", msai_long$variable)
  msai_long$region <- sub("\\.", "-", msai_long$region)
  msai_long <- msai_long[msai_long$region == region, ]
  msai_long$seqnames <- paste0("chr", msai_long$seqnames)
  
  #map the start and end positions to the region specific start and end positions
  # for subsequent analysis let's use minimum consistent regions, might have to have a single regions and a multi region mode
  segment_map$seqnames <- paste0("chr", segment_map$seqnames)
  cn_data$chr <- paste0("chr", cn_data$chr)
  cn_data <- left_join(cn_data, segment_map, by = c("chr" = "seqnames",
                                                    "startpos" = "regionCN_seg_start",
                                                     "endpos" = "regionCN_seg_end"))
  
  cn_data$startpos <- NULL
  cn_data$endpos <- NULL
  
  colnames(cn_data)[(ncol(cn_data) -1) : ncol(cn_data)] <- c("seg_start", "seg_end")
  
  cn_data <- unique(cn_data)
  
  cn_data <- left_join(cn_data, msai_long[, c("seqnames", "start", "end", "msai", "flip", "allelic_imbalance")], by = c("chr" = "seqnames",
                                                                                                                                        "seg_start" = "start",
                                                                                                                                        "seg_end" = "end"))
  
  cn_data <- cn_data[!is.na(cn_data$nMajor), ]
  cn_data <- cn_data[!is.na(cn_data$nMinor), ]
  
  cn_data <- cn_data[, c("chr", "seg_start", "seg_end", "nAraw", "nBraw", "msai", "flip", "allelic_imbalance")]
  colnames(cn_data) <- c("chr", "seg_start", "seg_end", "maj_cn", "min_cn", "msai", "flip", "allelic_imbalance")
  
  #now using the flip collumn we can assign copy number to A and B instead major and minor based on being the greater value
  cn_data[which(cn_data$flip == FALSE), "cnA"] <- cn_data[which(cn_data$flip == FALSE), "maj_cn"]
  cn_data[which(cn_data$flip == FALSE), "cnB"] <- cn_data[which(cn_data$flip == FALSE), "min_cn"]
  cn_data[which(cn_data$flip == TRUE), "cnA"] <- cn_data[which(cn_data$flip == TRUE), "min_cn"]
  cn_data[which(cn_data$flip == TRUE), "cnB"] <- cn_data[which(cn_data$flip == TRUE), "maj_cn"]
  
  #lets remove the old non-flipped copy number
  cn_data$maj_cn <- NULL
  cn_data$min_cn <- NULL
  
  colnames(cn_data)[(ncol(cn_data) -1) : ncol(cn_data)] <- c("maj_cn", "min_cn")
  
  return(cn_data)
  
}



#function to determine whole genome doubling, written by Alex Frankell and adapted by KT
detect_gd <- function( seg ){
  
  # only keep the copy number info
  seg <- unique(seg[, c("chr", "maj_cn", "min_cn", "seg_start", "seg_end")])

  #lets unphase the copy number for these purposes
  seg[seg$maj_cn < seg$min_cn, "major_cn"] <- seg[seg$maj_cn < seg$min_cn, "min_cn"] 
  seg[seg$maj_cn < seg$min_cn, "minor_cn"] <- seg[seg$maj_cn < seg$min_cn, "maj_cn"] 

  seg[seg$maj_cn > seg$min_cn, "major_cn"] <- seg[seg$maj_cn > seg$min_cn, "maj_cn"] 
  seg[seg$maj_cn > seg$min_cn, "minor_cn"] <- seg[seg$maj_cn > seg$min_cn, "min_cn"] 

  seg[seg$maj_cn == seg$min_cn, "major_cn"] <- seg[seg$maj_cn == seg$min_cn, "maj_cn"] 
  seg[seg$maj_cn == seg$min_cn, "minor_cn"] <- seg[seg$maj_cn == seg$min_cn, "min_cn"]

  seg$maj_cn <- NULL
  seg$min_cn <- NULL
  
  seg <- seg[, c("chr", "major_cn", "minor_cn", "seg_start", "seg_end")]
  colnames(seg) <- c("chr", "maj_cn", "min_cn", "seg_start", "seg_end")

  # which segments are in autosomes
  autosomal_index <- seg$chr %in% paste0("chr", 1:22)
  
  # determine width of each segment
  segment_widths <- seg$seg_end - seg$seg_start
  
  # determine the % of genome (where CN is called) which each segment contains
  segment_genome_fraction <- segment_widths / sum( segment_widths )
  
  # determine fraction of genme at alleicbalenece (ie nMaj = nMin)
  genome_frac_balence <- sum( segment_genome_fraction[ seg$maj_cn == seg$min_cn ] ) 
  
  # determine mean differnce in alleles accross the genome
  mean_allelic_difference <- sum( (seg$maj_cn - seg$min_cn) * segment_genome_fraction ) 
  
  # determine fraction of autosomal genome with Major allele => 2, if > 50% >= 1 GD has occurred
  genome_frac_nMaj2 <- sum( segment_genome_fraction[ autosomal_index & seg$maj_cn >= 2 ] ) 
  first_gd <- genome_frac_nMaj2 > 0.5
  
  # determine fraction of autosomal genome with Major allele => 3, if > 50% >= 2 GDa have occurred
  genome_frac_nMaj3 <- sum( segment_genome_fraction[ autosomal_index & seg$maj_cn >= 3 ] ) 
  second_gd <- genome_frac_nMaj3 > 0.5
  
  # determine fraction of autosomal genome with Major allele => 4
  genome_frac_nMaj4 <- sum( segment_genome_fraction[ autosomal_index & seg$maj_cn >= 4 ] ) 
  
  # as.char required to allowing mixingof numeric and logical data types
  out <- data.frame( first_gd, genome_frac_nMaj2, second_gd, genome_frac_nMaj3, genome_frac_nMaj4, genome_frac_balence, mean_allelic_difference )
  
  return( out )
}
   
# function to phase mutations to haplotypes 
phasing.phun <- function(snp_mutation_df, tumour, region, ploidy, output_path){
  
    # snp_mutation_df[snp_mutation_df$bb_BAF < 0.5, "bb_phase"] <- "B"
    # snp_mutation_df[snp_mutation_df$bb_BAF >= 0.5, "bb_phase"] <- "A"
    
    snp_mutation_df[snp_mutation_df$bb_mutisample_haplo == 0, "bb_phase"] <- "B"
    snp_mutation_df[snp_mutation_df$bb_mutisample_haplo == 1, "bb_phase"] <- "A"
    
    snp_mutation_df[snp_mutation_df$bb_phase == "B", "snp_vaf_on_mut_bb_phase"] <- abs(1 - snp_mutation_df[snp_mutation_df$bb_phase == "B", "snp_vaf_on_mut"])
    snp_mutation_df[snp_mutation_df$bb_phase == "A", "snp_vaf_on_mut_bb_phase"] <- snp_mutation_df[snp_mutation_df$bb_phase == "A", "snp_vaf_on_mut"] 
    
    snp_mutation_df$mutation_id <- paste0(snp_mutation_df$chr, "_", snp_mutation_df$pos)                
    snp_mutation_df <- snp_mutation_df %>% group_by(mutation_id) %>% mutate(average_snp_vaf_on_mut_bb_phase = mean(snp_vaf_on_mut_bb_phase))                
    snp_mutation_df <- data.frame(snp_mutation_df)             
    
    mutation_df <- snp_mutation_df           
    mutation_df <- unique(mutation_df[, c("sample", "chr", "pos", "maj_cn", "min_cn", "tot_cn", "seg_start", "seg_end", "mut_cpn",
                                          "ref_mut", "alt_mut", "counts_mut", "average_snp_vaf_on_mut_bb_phase",
                                          "sum_SNV_MT_SNP_MT_N", "sum_SNV_MT_SNP_WT_N", "sum_SNV_WT_SNP_MT_N", "sum_SNV_WT_SNP_WT_N", "parallel_mutation", "allelic_imbalance", "msai")])
    
    mutation_df[mutation_df$average_snp_vaf_on_mut_bb_phase >= 0.5, "phase"] <- "A"
    mutation_df[mutation_df$average_snp_vaf_on_mut_bb_phase < 0.5, "phase"] <- "B"

    mutation_df <- mutation_df[!is.na(mutation_df$maj_cn), ]
    mutation_df <- mutation_df[!is.na(mutation_df$min_cn), ]
    mutation_df <- mutation_df[!is.na(mutation_df$mut_cpn), ]
    
    #determine if a gain is present, this needs to differ from diploid to wgd
    #so first, determine if the tumour is genome doubled
    # wgd <- detect_gd(mutation_df)$first_gd

    # if(wgd == FALSE){

      mutation_df[which(mutation_df$maj_cn > 1), "maj_cn_gain"] <- TRUE
      mutation_df[which(mutation_df$maj_cn <= 1), "maj_cn_gain"] <- FALSE
      mutation_df[which(mutation_df$min_cn > 1), "min_cn_gain"] <- TRUE
      mutation_df[which(mutation_df$min_cn <= 1), "min_cn_gain"] <- FALSE

    # }

    #  if(wgd){
    # 
    #   thresh <- ceiling(ploidy)
    #   if(thresh %% 2 == 0){
    #     mutation_df[which(mutation_df$maj_cn > thresh/2), "maj_cn_gain"] <- TRUE
    #     mutation_df[which(mutation_df$maj_cn <= thresh/2), "maj_cn_gain"] <- FALSE
    #     mutation_df[which(mutation_df$min_cn > thresh/2), "min_cn_gain"] <- TRUE
    #     mutation_df[which(mutation_df$min_cn <= thresh/2), "min_cn_gain"] <- FALSE
    #   }
    # 
    #   if(thresh %% 2 != 0){
    #     
    #     thresh <- thresh+1
    # 
    #     mutation_df[which(mutation_df$maj_cn > thresh/2), "maj_cn_gain"] <- TRUE
    #     mutation_df[which(mutation_df$maj_cn <= thresh/2), "maj_cn_gain"] <- FALSE
    #     mutation_df[which(mutation_df$min_cn > thresh/2), "min_cn_gain"] <- TRUE
    #     mutation_df[which(mutation_df$min_cn <= thresh/2), "min_cn_gain"] <- FALSE
    #   }
    # 
    # }
    
    # where there is no allelic imbalance lets remove the phase
    mutation_df[mutation_df$allelic_imbalance == FALSE, "phase"] <- NA

    # save this
    if(!dir.exists(paste0(output_path, tumour, "/", region, "/data"))){
      dir.create(paste0(output_path, tumour, "/", region, "/data"), recursive = T)
    }
    
    write.table(mutation_df, paste0(output_path, tumour, "/", region, "/data/", region, "_phased_mutations.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
    
    return(mutation_df)
  
}  

# function to make battenberg data into ascat like clonal format
# bb.into.ascat.form.fun <- function(seg.table){
#   
#   seg.table$nAraw <- sapply(1:nrow(seg.table), function(x){
#     ifelse(seg.table$frac1_A[x]==1, 
#            seg.table$nMaj1_A[x], 
#            seg.table$nMaj1_A[x]*seg.table$frac1_A[x]+seg.table$nMaj2_A[x]*seg.table$frac2_A[x])
#   })
#   seg.table$nBraw <- sapply(1:nrow(seg.table), function(x){
#     ifelse(seg.table$frac1_A[x]==1, 
#            seg.table$nMin1_A[x], 
#            seg.table$nMin1_A[x]*seg.table$frac1_A[x]+seg.table$nMin2_A[x]*seg.table$frac2_A[x])
#   })
#   
#   seg <- data.frame(chr = seg.table$chr,
#                     startpos = seg.table$startpos,
#                     endpos = seg.table$endpos,
#                     n.het = NA,
#                     cnTotal = NA,
#                     nMajor = round(seg.table$nAraw),
#                     nMinor = round(seg.table$nBraw),
#                     nAraw = seg.table$nAraw,
#                     nBraw = seg.table$nBraw)
#   
#   seg$cnTotal <- as.numeric(seg$nMajor) + as.numeric(seg$nMinor)
#   seg$chr <- paste0("chr", seg$chr)
#   
#   return(seg)
# 
# }
   
# function to plot mutations colored by haplotype across the genome and seperately for each chromosome
phun.plots <- function(mutation_df = phased_data, tumour, region, all_mut_data = data, cn, output_path, genome = BSgenome.Hsapiens.UCSC.hg38){
  
  #we've gotta add non phased mutations back on so we can plot those but in grey
  all_chr_pos <- unique(all_mut_data[, c("chr", "pos", "mut_cpn", "maj_cn", "min_cn")])
  
  all_chr_pos <- left_join(all_chr_pos, mutation_df)
  
  mutation_df <- all_chr_pos
  mutation_df[is.na(mutation_df$phase), "phase"] <- "no_allelic_imbalance"
  mutation_df$phase <- factor(mutation_df$phase, levels = c("A", "B", "no_allelic_imbalance"))
  
  # we need to make sure that the two copy number alleles are not plotted on top of each other
  # but we also don't want to move them too much so we can still distinguish subclonality

  mutation_df <- mutation_df[!is.na(mutation_df$maj_cn), ]
  mutation_df <- mutation_df[!is.na(mutation_df$min_cn), ]
  
  mutation_df[mutation_df$maj_cn == mutation_df$min_cn, "seg_y_pos_maj"] <- mutation_df[mutation_df$maj_cn == mutation_df$min_cn, "maj_cn"] + 0.025
  mutation_df[mutation_df$maj_cn != mutation_df$min_cn, "seg_y_pos_maj"] <- mutation_df[mutation_df$maj_cn != mutation_df$min_cn, "maj_cn"]
  mutation_df[mutation_df$maj_cn == mutation_df$min_cn, "seg_y_pos_min"] <- mutation_df[mutation_df$maj_cn == mutation_df$min_cn, "min_cn"] - 0.025
  mutation_df[mutation_df$maj_cn != mutation_df$min_cn, "seg_y_pos_min"] <- mutation_df[mutation_df$maj_cn != mutation_df$min_cn, "min_cn"]
  
  mutation_df$allelic_imbalance <- factor(mutation_df$allelic_imbalance, levels = c("TRUE", "FALSE"))
  
  # make sure our segments also have a legend
  mutation_df$seg_legendA <- "haplotype A"
  mutation_df$seg_legendB <- "haplotype B"
  
  # also plot it across the whole genome
  chr_df <- genome.continous.fun(mutation_df, genome)
  chr_sizes <- seqlengths(genome)[paste0("chr", 1:22)]
  chr_lines <- cumsum(as.numeric(as.character(chr_sizes)))
  chr_midpoint <- chr_lines - chr_sizes / 2
  
  xaxis_label <- data.frame(chr_midpoint)
  xaxis_label$label <- rownames(xaxis_label)

  # for the cross genome plot let's limit the highest cn because otherwise we can't see anything, let's set the limit to 10
  chr_df[chr_df$mut_cpn > 10, "mut_cpn_lim"] <- 10
  chr_df[chr_df$mut_cpn <= 10, "mut_cpn_lim"] <- chr_df[chr_df$mut_cpn <= 10, "mut_cpn"] 
  chr_df[chr_df$seg_y_pos_min > 10, "seg_y_pos_min_lim"] <- 10
  chr_df[chr_df$seg_y_pos_min <= 10, "seg_y_pos_min_lim"] <- chr_df[chr_df$seg_y_pos_min <= 10, "seg_y_pos_min"]
  chr_df[chr_df$seg_y_pos_maj > 10, "seg_y_pos_maj_lim"] <- 10
  chr_df[chr_df$seg_y_pos_maj <= 10, "seg_y_pos_maj_lim"] <- chr_df[chr_df$seg_y_pos_maj <= 10, "seg_y_pos_maj"]

  mut_cpn_ylim_genome <- max(c(chr_df$mut_cpn_lim, chr_df$seg_y_pos_maj_lim))
  mut_cpn_ylim_genome_min <- min(c(chr_df$mut_cpn_lim, chr_df$seg_y_pos_min_lim))
  
  p_genome <- ggplot(chr_df, aes(adjust_pos, mut_cpn_lim, fill = phase), color = "black") +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
    geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min_lim", "seg_y_pos_maj_lim", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_maj_lim, yend = seg_y_pos_maj_lim, color = 'haplotype A', alpha = allelic_imbalance), linewidth = 2) +
    geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min_lim", "seg_y_pos_maj_lim", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_min_lim, yend = seg_y_pos_min_lim, color = 'haplotype B', alpha = allelic_imbalance),linewidth = 2) +
    scale_fill_manual(name = "mutation haplotype", breaks = c("A", "B", "no_allelic_imbalance"), labels = c("allele A", "allele B", "no allelic imbalance"), values = c("#d73027", "#4575b4", "#bababa")) +
    scale_color_manual(name = "copy number haplotype", values=c("haplotype A" = "darkred", "haplotype B" = "#313695")) +
    scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10), limits = c(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome)) +
    ylab("mutation copy number") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  p_genome_red <- ggplot(chr_df[chr_df$phase == "A", ], aes(adjust_pos, mut_cpn_lim, fill = phase), color = "black") +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
    geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min_lim", "seg_y_pos_maj_lim", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_maj_lim, yend = seg_y_pos_maj_lim, color = 'haplotype A', alpha = allelic_imbalance), linewidth = 2) +
    scale_fill_manual(name = "mutation haplotype", breaks = c("A", "B", "no_allelic_imbalance"), labels = c("allele A", "allele B", "no allelic imbalance"), values = c("#d73027", "#4575b4", "#bababa")) +
    scale_color_manual(name = "copy number haplotype", values=c("haplotype A" = "darkred", "haplotype B" = "#313695")) +
    scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10), limits = c(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome)) +
    ylab("mutation copy number") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  p_genome_blue <- ggplot(chr_df[chr_df$phase == "B", ], aes(adjust_pos, mut_cpn_lim, fill = phase), color = "black") +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
    geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min_lim", "seg_y_pos_maj_lim", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_min_lim, yend = seg_y_pos_min_lim, color = 'haplotype B', alpha = allelic_imbalance),linewidth = 2) +
    scale_fill_manual(name = "mutation haplotype", breaks = c("A", "B", "no_allelic_imbalance"), labels = c("allele A", "allele B", "no allelic imbalance"), values = c("#d73027", "#4575b4", "#bababa")) +
    scale_color_manual(name = "copy number haplotype", values=c("haplotype A" = "darkred", "haplotype B" = "#313695")) +
    scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10), limits = c(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome)) +
    ylab("mutation copy number") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  p_genome_grey <- ggplot(chr_df, aes(adjust_pos, mut_cpn_lim, fill = phase), color = "black") +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
    geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min_lim", "seg_y_pos_maj_lim", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_maj_lim, yend = seg_y_pos_maj_lim, color = 'haplotype A', alpha = allelic_imbalance), linewidth = 2) +
    geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min_lim", "seg_y_pos_maj_lim", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_min_lim, yend = seg_y_pos_min_lim, color = 'haplotype B', alpha = allelic_imbalance),linewidth = 2) +
    scale_fill_manual(name = "mutation haplotype", breaks = c("A", "B", "no_allelic_imbalance"), labels = c("allele A", "allele B", "no allelic imbalance"), values = c("#878787", "#878787", "#878787")) +
    scale_color_manual(name = "copy number haplotype", values=c("haplotype A" = "darkred", "haplotype B" = "#313695")) +
    scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10), limits = c(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome)) +
    ylab("mutation copy number") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())

  p_genome_blank <- ggplot(chr_df, aes(adjust_pos, mut_cpn_lim, fill = phase), color = "black") +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min_lim", "seg_y_pos_maj_lim", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_maj_lim, yend = seg_y_pos_maj_lim, color = 'haplotype A', alpha = allelic_imbalance), linewidth = 2) +
    geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min_lim", "seg_y_pos_maj_lim", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_min_lim, yend = seg_y_pos_min_lim, color = 'haplotype B', alpha = allelic_imbalance),linewidth = 2) +
    scale_fill_manual(name = "mutation haplotype", breaks = c("A", "B", "no_allelic_imbalance"), labels = c("allele A", "allele B", "no allelic imbalance"), values = c("#d73027", "#4575b4", "#bababa")) +
    scale_color_manual(name = "copy number haplotype", values=c("haplotype A" = "darkred", "haplotype B" = "#313695")) +
    scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10), limits = c(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome)) +
    ylab("copy number") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())
  
  if(!dir.exists(paste0(output_path, tumour, "/", region, "/plots/phased"))){
    dir.create(paste0(output_path, tumour, "/", region, "/plots/phased"), recursive = T)
  }
  
  pdf(paste0(output_path, tumour, "/", region, "/plots/phased/", region, "_genome_phased_mutation_copynumber_red_blue.pdf"), width = 18, height = 3.5)
  print(p_genome)
  dev.off()
  
  pdf(paste0(output_path, tumour, "/", region, "/plots/phased/", region, "_genome_phased_mutation_copynumber_red.pdf"), width = 18, height = 3.5)
  print(p_genome_red)
  dev.off()
  
  pdf(paste0(output_path, tumour, "/", region, "/plots/phased/", region, "_genome_phased_mutation_copynumber_blue.pdf"), width = 18, height = 3.5)
  print(p_genome_blue)
  dev.off()
  
  pdf(paste0(output_path, tumour, "/", region, "/plots/phased/", region, "_genome_phased_mutation_copynumber_grey.pdf"), width = 18, height = 3.5)
  print(p_genome_grey)
  dev.off()
  
  pdf(paste0(output_path, tumour, "/", region, "/plots/phased/", region, "_genome_phased_mutation_copynumber_CN.pdf"), width = 18, height = 3.5)
  print(p_genome_blank)
  dev.off()
  
  
  # make one plot for each chromosome
  
  for(chr in unique(mutation_df$chr)){
    
    mutation_chr_df <- mutation_df[mutation_df$chr == chr, ]
    mutation_chr_df[mutation_chr_df$mut_cpn > 10, "mut_cpn_lim"] <- 10
    mutation_chr_df[mutation_chr_df$mut_cpn <= 10, "mut_cpn_lim"] <- mutation_chr_df[mutation_chr_df$mut_cpn <= 10, "mut_cpn"] 
    mutation_chr_df[mutation_chr_df$seg_y_pos_min > 10, "seg_y_pos_min_lim"] <- 10
    mutation_chr_df[mutation_chr_df$seg_y_pos_min <= 10, "seg_y_pos_min_lim"] <- mutation_chr_df[mutation_chr_df$seg_y_pos_min <= 10, "seg_y_pos_min"]
    mutation_chr_df[mutation_chr_df$seg_y_pos_maj > 10, "seg_y_pos_maj_lim"] <- 10
    mutation_chr_df[mutation_chr_df$seg_y_pos_maj <= 10, "seg_y_pos_maj_lim"] <- mutation_chr_df[mutation_chr_df$seg_y_pos_maj <= 10, "seg_y_pos_maj"]
    
    mut_cpn_ylim_chr <- max(mutation_chr_df$mut_cpn_lim, mutation_chr_df$seg_y_pos_min_lim, mutation_chr_df$seg_y_pos_maj_lim)

    p <- ggplot(mutation_chr_df, aes(pos, mut_cpn_lim, fill = phase), color = "black") +
      geom_segment(data = unique(mutation_chr_df[, c("seg_start", "seg_end", "seg_y_pos_maj_lim", "seg_y_pos_min_lim", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_maj_lim, yend = seg_y_pos_maj_lim, color = 'haplotype A', alpha = allelic_imbalance), linewidth = 2) +
      geom_segment(data = unique(mutation_chr_df[, c("seg_start", "seg_end", "seg_y_pos_maj_lim", "seg_y_pos_min_lim", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_min_lim, yend = seg_y_pos_min_lim, color = 'haplotype B', alpha = allelic_imbalance), linewidth = 2) +
      geom_point(alpha = 0.7, shape = 21, stroke = NA) +
      scale_fill_manual(name = "mutation phase", breaks = c("A", "B", "no_allelic_imbalance"), labels = c("haplotype A", "haplotype B", "no allelic imbalance"), values = c("#d73027", "#4575b4", "#bababa")) +
      scale_color_manual(name = "copy number haplotype", values=c("haplotype A" = "darkred", "haplotype B" = "#313695")) +
      scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
      ylab("mutation copy number") +
      ylim(0, mut_cpn_ylim_chr) +
      xlab(chr) +
      theme_bw() 
    
    p_red <- ggplot(mutation_chr_df[mutation_chr_df$phase == "A", ], aes(pos, mut_cpn_lim, fill = phase), color = "black") +
      geom_segment(data = unique(mutation_chr_df[, c("seg_start", "seg_end", "seg_y_pos_min_lim", "seg_y_pos_maj_lim", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_maj_lim, yend = seg_y_pos_maj_lim, color = 'haplotype A', alpha = allelic_imbalance), linewidth = 2) +
      geom_point(alpha = 0.7, shape = 21, stroke = NA) +
      scale_fill_manual(name = "mutation haplotype", breaks = c("A", "B", "no_allelic_imbalance"), labels = c("haplotype A", "haplotype B", "no allelic imbalance"), values = c("#d73027", "#4575b4", "#bababa")) +
      scale_color_manual(name = "copy number haplotype", values=c("haplotype A" = "darkred", "haplotype B" = "#313695")) +
      scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c(FALSE, TRUE), label = c("FALSE", "TRUE")) +
      ylab("mutation copy number") +
      ylim(0, mut_cpn_ylim_chr) +
      xlab(chr) +
      theme_bw() 
    
    p_blue <- ggplot(mutation_chr_df[mutation_chr_df$phase == "B", ], aes(pos, mut_cpn_lim, fill = phase), color = "black") +
      geom_segment(data = unique(mutation_chr_df[, c("seg_start", "seg_end", "seg_y_pos_min_lim", "seg_y_pos_maj_lim", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_min_lim, yend = seg_y_pos_min_lim, color = 'haplotype B', alpha = allelic_imbalance),linewidth = 2) +
      geom_point(alpha = 0.7, shape = 21, stroke = NA) +
      scale_fill_manual(name = "mutation haplotype", breaks = c("A", "B", "no_allelic_imbalance"), labels = c("haplotype A", "haplotype B", "no allelic imbalance"), values = c("#d73027", "#4575b4", "#bababa")) +
      scale_color_manual(name = "copy number haplotype", values=c("haplotype A" = "darkred", "haplotype B" = "#313695")) +
      scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
      ylab("mutation copy number") +
      ylim(0, mut_cpn_ylim_chr) +
      xlab(chr) +
      theme_bw() 
    
    p_grey <- ggplot(mutation_chr_df, aes(pos, mut_cpn_lim, fill = phase), color = "black") +
      geom_segment(data = unique(mutation_chr_df[, c("seg_start", "seg_end", "seg_y_pos_min_lim", "seg_y_pos_maj_lim", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_maj_lim, yend = seg_y_pos_maj_lim, color = 'haplotype A', alpha = allelic_imbalance), linewidth = 2) +
      geom_segment(data = unique(mutation_chr_df[, c("seg_start", "seg_end", "seg_y_pos_min_lim", "seg_y_pos_maj_lim", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_min_lim, yend = seg_y_pos_min_lim, color = 'haplotype B', alpha = allelic_imbalance),linewidth = 2) +
      geom_point(alpha = 0.7, shape = 21, stroke = NA) +
      scale_fill_manual(name = "mutation haplotype", breaks = c("A", "B", "no_allelic_imbalance"), labels = c("haplotype A", "haplotype B", "no allelic imbalance"), values = c("#878787", "#878787", "#878787")) +
      scale_color_manual(name = "copy number haplotype", values=c("haplotype A" = "darkred", "haplotype B" = "#313695")) +
      scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
      ylab("mutation copy number") +
      ylim(0, mut_cpn_ylim_chr) +
      xlab(chr) +
      theme_bw() 
    
    p_blank <- ggplot(mutation_chr_df, aes(pos, mut_cpn_lim, fill = phase), color = "black") +
      geom_segment(data = unique(mutation_chr_df[, c("seg_start", "seg_end", "seg_y_pos_min_lim", "seg_y_pos_maj_lim", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_maj_lim, yend = seg_y_pos_maj_lim, color = 'haplotype A', alpha = allelic_imbalance), linewidth = 2) +
      geom_segment(data = unique(mutation_chr_df[, c("seg_start", "seg_end", "seg_y_pos_min_lim", "seg_y_pos_maj_lim", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_min_lim, yend = seg_y_pos_min_lim, color = 'haplotype B', alpha = allelic_imbalance),linewidth = 2) +
      scale_color_manual(name = "copy number haplotype", values=c("haplotype A" = "darkred", "haplotype B" = "#313695")) +
      scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
      scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), limits = c(0, mut_cpn_ylim_chr)) +
      ylab("mutation copy number") +
      xlab(chr) +
      theme_bw() +
      theme(panel.grid.minor.y = element_blank())
    
    
    if(!dir.exists(paste0(output_path, tumour, "/", region, "/plots/phased/per_chromosome/"))){
      dir.create(paste0(output_path, tumour, "/", region, "/plots/phased/per_chromosome"), recursive = T)
    }
    
    pdf(paste0(output_path, tumour, "/", region, "/plots/phased/per_chromosome/", region, "_", chr, "_phased_mutation_copynumber_red_blue.pdf"), width = 9, height = 3.5)
    print(p)
    dev.off()
    
    pdf(paste0(output_path, tumour, "/", region, "/plots/phased/per_chromosome/", region, "_", chr, "_phased_mutation_copynumber_red.pdf"), width = 9, height = 3.5)
    print(p_red)
    dev.off()
    
    pdf(paste0(output_path, tumour, "/", region, "/plots/phased/per_chromosome/", region, "_", chr, "_phased_mutation_copynumber_blue.pdf"), width = 9, height = 3.5)
    print(p_blue)
    dev.off()
    
    pdf(paste0(output_path, tumour, "/", region, "/plots/phased/per_chromosome/", region, "_", chr, "_phased_mutation_copynumber_grey.pdf"), width = 9, height = 3.5)
    print(p_grey)
    dev.off()
    
    pdf(paste0(output_path, tumour, "/", region, "/plots/phased/per_chromosome/", region, "_", chr, "_phased_mutation_copynumber_CN.pdf"), width = 9, height = 3.5)
    print(p_blank)
    dev.off()

  }
  
}

#function that works out which cluster mutations are on that segment
segment.cluster.timing.fun <- function(data = phased_data, tree_object, output_path){
  
  #for each segment, let's see if we have mutations from several clusters on that segment
  data$segID <- paste0(data$chr, "_", data$seg_start, "_", data$seg_end)
  data <- data %>% group_by(segID) %>% mutate(num_clusters_on_seg = length(unique(CLUSTER)))
  data <- data.frame(data)
  
  #now we need to check if there both early and late mutations within one of the clusters
  # if a segment is on the integer value lets use integers to time mutations as early and late
  # if a segment is not on the integer it is a subclonal gain, this means two things:
    # we should only be using subclonal mutations here
    # and the integer thing won't work, so here use local optima to define early and late
  
  # seg = "chr1_231605523_234078809"
  # seg = chr3_93719007_123865333
  
  timed_data <- lapply(unique(data$segID), function(seg){
    
    print(seg)
    seg_df <- data[data$segID == seg, ]
    
    #need to do the next steps per haplotype
    
    dfA <- seg_df[which(seg_df$phase == "A"), ]
    dfB <- seg_df[which(seg_df$phase == "B"), ]
    
    #start with haplotype A
    #first, does this segment have a timeable copy number gain?

    if(all(nrow(dfA) > 0 && all(dfA$maj_cn_gain))){
      
      #and for now only use clonal copy number gains
      if(as.numeric(as.character(unique(dfA$maj_cn))) == round(as.numeric(as.character(unique(dfA$maj_cn))))){
        
        # which clusters have several copy number steps? need to define steps first
        # because this is a clonal gain we can use the simple definition
        dfA$local_optimum <- round(dfA$mut_cpn)
        dfA[dfA$local_optimum > dfA$maj_cn, "local_optimum"] <- unique(dfA$maj_cn)
        
        #in theory I should not have any muts at 0
        #but this might be the case due to some noise and then rounding down
        # let's remove those
        dfA <- dfA[dfA$local_optimum >= 1, ]
        
        if(nrow(dfA) > 0){
          #we need to adjust for some noise where there are less than 5 mutations on one state, these mutations get removed
          optima_count <- data.frame(table(dfA$local_optimum))
          colnames(optima_count) <- c("local_optimum", "muts_per_optimum")
          optima_count$local_optimum <- as.numeric(as.character(optima_count$local_optimum))
          dfA <- left_join(dfA, optima_count)
          dfA <- dfA[dfA$muts_per_optimum >= 5, ]
          dfA$muts_per_optimum <- NULL
          
          counts <- data.frame(table(dfA$local_optimum, dfA$CLUSTER))
          counts <- counts[counts$Freq > 0, ]
          
          interesting_cluster <- names(table(counts$Var2)[table(counts$Var2) > 1])
          
          # for this to be a real gain the mut cpn must be at 1, or an integer above 1  and at an integer above even that
          interesting_cluster_muts <- dfA[dfA$CLUSTER %in% interesting_cluster, ]
        } else { 
          interesting_cluster_muts <- data.frame(matrix(NA, 1, 1)) 
          interesting_cluster_muts <- na.omit(interesting_cluster_muts)
        }
        
        if(nrow(interesting_cluster_muts) > 0){
          #let's eliminate those clusters where this is not the case
          cpn_cl_states <- data.frame(table(interesting_cluster_muts$local_optimum, interesting_cluster_muts$CLUSTER))
          cpn_cl_states$Var1 <- as.numeric(as.character(cpn_cl_states$Var1))
          cpn_cl_states <- cpn_cl_states[cpn_cl_states$Freq > 0 & cpn_cl_states$Var1 >= 1, ]
          
          usable_cluster_muts <- interesting_cluster_muts[interesting_cluster_muts$CLUSTER %in% unique(cpn_cl_states$Var2), ]
          
          if(nrow(usable_cluster_muts) > 0){
            
            #let's try to create a flag that says whether a gain is only seen on one cluster or on several
            optimum_cluster <- data.frame(table(usable_cluster_muts$CLUSTER, usable_cluster_muts$local_optimum))
            colnames(optimum_cluster) <- c("CLUSTER", "local_optimum", "Freq")
            optimum_cluster <- optimum_cluster[optimum_cluster$Freq > 0, ]
            optimum_cluster <- optimum_cluster[optimum_cluster$local_optimum != 0, ]
            
            # we only want to keep the gain in the clonal cluster which is not also detected in a subclonal cluster
            trunk <- tree_object$graph_pyclone$trunk
            
            optimum_cluster[optimum_cluster$CLUSTER == trunk, "trunk"] <- "trunk"
            optimum_cluster[optimum_cluster$CLUSTER != trunk, "trunk"] <- "not_trunk"
            
            optimum_df <- data.frame(table(optimum_cluster$local_optimum))
            colnames(optimum_df) <- c("local_optimum", "Freq")
            optimum_df <- optimum_df[optimum_df$Freq != 0, ]
            optimum_to_use <- optimum_df[optimum_df$Freq == 1, ]
            optimum_to_use <- left_join(optimum_to_use, optimum_cluster[, c("CLUSTER", "local_optimum", "trunk")], by = c("local_optimum" = "local_optimum"))
            
            #now if this has the trunk in it we can use it
            optimum_to_use <- optimum_to_use[optimum_to_use$trunk == "trunk", ]
            
            #now for each cluster separately we can time the gains
            if(nrow(optimum_to_use) > 0){
              
              cluster_gain_time_df <- list()
              for(cluster in unique(optimum_to_use$CLUSTER)){
                

                df <- usable_cluster_muts[usable_cluster_muts$CLUSTER == cluster, ]
              
                
                #only keep the first gain
                cl_optimum <- optimum_to_use$local_optimum
                
                # late_optimum <- as.numeric(as.character(unique(optimum_cluster[optimum_cluster$CLUSTER == cluster,"local_optimum"])))
                # late_optimum <- sort(late_optimum, decreasing = F)
                # late_optimum <- late_optimum[length(late_optimum)-1]
                
                # optima_to_use <- c(as.numeric(as.character(cl_optimum)), late_optimum)
                optima_to_use <- c(as.numeric(as.character(cl_optimum)))
                
                # how many copy number steps are there
                number_steps <- length(optima_to_use)
                
                #for each copy number state count the numbers of mutations
                df$mutID <- paste0(df$chr, ":", df$pos)
                df <- df %>% group_by(local_optimum) %>% mutate(num_muts_on_cn_step = n())
                df <- data.frame(df)
                
                #only keep relevant info
                mutation_count_df <- unique(df[, c("sample", "chr", "seg_start", "seg_end", "phase", "maj_cn_gain", "min_cn_gain", "local_optimum", "num_muts_on_cn_step")])
                mutation_count_df <- mutation_count_df[order(mutation_count_df$local_optimum, decreasing = F), ]
                
                mutation_count_df$copy_number_before_gain <- mutation_count_df$local_optimum
                
                for(s in 1:nrow(mutation_count_df)){
                  mutation_count_df[s, "copy_number_after_gain"] <- ifelse(s < nrow(mutation_count_df), mutation_count_df[s+1, "copy_number_before_gain"], mutation_count_df[s, "copy_number_before_gain"])
                }
                
                mutation_count_df$adjustment_factor <- mutation_count_df$copy_number_before_gain
                
                #if adjustment factor isbelow one make it into one, although this should not happen anymore
                mutation_count_df[mutation_count_df$adjustment_factor < 1, "adjustment_factor"] <- 1
                
                mutation_count_df$number_muts_cn_adjusted <- mutation_count_df$num_muts_on_cn_step / mutation_count_df$adjustment_factor
                mutation_count_df <- mutation_count_df[mutation_count_df$local_optimum %in% optima_to_use, ]
                
                gain_time_df <- data.frame(tumour = NA,
                                           region = NA,
                                           segID = NA,
                                           gain_index = NA,
                                           copy_number_before_gain = NA,
                                           copy_number_after_gain = NA,
                                           number_early_muts_cn_adjusted = NA,
                                           number_all_muts_cn_adjusted = NA,
                                           gain_timing = NA,
                                           loCI = NA,
                                           hiCI = NA,
                                           haplotype = NA) 
                
                for(event in 1:(number_steps - 1)){
                  
                  gain_index                            <- event
                  copy_number_before_gain               <- mutation_count_df[event, "copy_number_before_gain"]
                  copy_number_after_gain                <- mutation_count_df[event, "copy_number_after_gain"]
                  
                  number_early_muts                     <- sum(mutation_count_df[(event + 1) : max(number_steps), "number_muts_cn_adjusted"])
                  # number_early_muts                     <- sum(mutation_count_df[1:event, "number_muts_cn_adjusted"])
                  number_all_muts                       <- sum(mutation_count_df[, "number_muts_cn_adjusted"])
                  
                  gain_timing                           <- number_early_muts / number_all_muts
                  error                                 <- qnorm(0.975) * sqrt(gain_timing * (1 - gain_timing) / number_all_muts)
                  loCI                                  <- gain_timing - error
                  hiCI                                  <- gain_timing + error
                  
                  gain_time_df[event, "tumour"]                                 <- tumour
                  gain_time_df[event, "region"]                                 <- region
                  gain_time_df[event, "segID"]                                  <- seg
                  gain_time_df[event, "cluster"]                                <- cluster
                  gain_time_df[event, "gain_index"]                             <- number_steps-gain_index
                  gain_time_df[event, "copy_number_before_gain"]                <- copy_number_before_gain
                  gain_time_df[event, "copy_number_after_gain"]                 <- copy_number_after_gain
                  gain_time_df[event, "number_early_muts_cn_adjusted"]          <- number_early_muts
                  gain_time_df[event, "number_all_muts_cn_adjusted"]            <- number_all_muts
                  gain_time_df[event, "gain_timing"]                            <- gain_timing
                  gain_time_df[event, "loCI"]                                   <- loCI
                  gain_time_df[event, "hiCI"]                                   <- hiCI
                  gain_time_df[event, "haplotype"]                              <- "A"
                  
                }  
                
                cluster_gain_time_df[[cluster]] <- gain_time_df
                
              }
              
              cluster_gain_time_df <- do.call(rbind, cluster_gain_time_df)
              
              haplotypeA_gaintiming <- cluster_gain_time_df
              
            } else {
              
              gain_time_df <- data.frame(tumour = tumour,
                                         region = region,
                                         segID = seg,
                                         cluster = NA,
                                         gain_index = NA,
                                         copy_number_before_gain = NA,
                                         copy_number_after_gain = NA,
                                         number_early_muts_cn_adjusted = NA,
                                         number_all_muts_cn_adjusted = NA,
                                         gain_timing = NA,
                                         loCI = NA,
                                         hiCI = NA,
                                         haplotype = "A") 
              
              
              haplotypeA_gaintiming <- gain_time_df
              
            }
            
          } else {
            
            gain_time_df <- data.frame(tumour = tumour,
                                       region = region,
                                       segID = seg,
                                       cluster = NA,
                                       gain_index = NA,
                                       copy_number_before_gain = NA,
                                       copy_number_after_gain = NA,
                                       number_early_muts_cn_adjusted = NA,
                                       number_all_muts_cn_adjusted = NA,
                                       gain_timing = NA,
                                       loCI = NA,
                                       hiCI = NA,
                                       haplotype = "A") 
            
            
            haplotypeA_gaintiming <- gain_time_df
            
          } 
        } else {
          
          gain_time_df <- data.frame(tumour = tumour,
                                     region = region,
                                     segID = seg,
                                     cluster = NA,
                                     gain_index = NA,
                                     copy_number_before_gain = NA,
                                     copy_number_after_gain = NA,
                                     number_early_muts_cn_adjusted = NA,
                                     number_all_muts_cn_adjusted = NA,
                                     gain_timing = NA,
                                     loCI = NA,
                                     hiCI = NA,
                                     haplotype = "A") 
          
          
          haplotypeA_gaintiming <- gain_time_df
        }
        
      } else {
        
        #need to deal with non clonal gains for now
        
        gain_time_df <- data.frame(tumour = tumour,
                                   region = region,
                                   segID = seg,
                                   cluster = NA,
                                   gain_index = NA,
                                   copy_number_before_gain = NA,
                                   copy_number_after_gain = NA,
                                   number_early_muts_cn_adjusted = NA,
                                   number_all_muts_cn_adjusted = NA,
                                   gain_timing = NA,
                                   loCI = NA,
                                   hiCI = NA,
                                   haplotype = "A") 
        
        
        haplotypeA_gaintiming <- gain_time_df
        
      }
      
    } else {
      
      #and gotta deal with segments that cant be time cause they aren't a gain
      
      gain_time_df <- data.frame(tumour = tumour,
                                 region = region,
                                 segID = seg,
                                 cluster = NA,
                                 gain_index = NA,
                                 copy_number_before_gain = NA,
                                 copy_number_after_gain = NA,
                                 number_early_muts_cn_adjusted = NA,
                                 number_all_muts_cn_adjusted = NA,
                                 gain_timing = NA,
                                 loCI = NA,
                                 hiCI = NA,
                                 haplotype = "A") 
      
      
      haplotypeA_gaintiming <- gain_time_df
      
      
    }
    
    
    #now do haplotype B
    
    if(all(nrow(dfB) > 0 && all(dfB$min_cn_gain))){
      
      #and for now only use clonal copy number gains
      if(as.numeric(as.character(unique(dfB$min_cn))) == round(as.numeric(as.character(unique(dfB$min_cn))))){
        
        # which clusters have several copy number steps? need to define steps first
        # because this is a clonal gain we can use the simple definition
        dfB$local_optimum <- round(dfB$mut_cpn)
        dfB[dfB$local_optimum > dfB$min_cn, "local_optimum"] <- unique(dfB$min_cn)
        
        #in theory I should not have any muts at 0
        #but this might be the case due to some noise and then rounding down
        # let's remove those
        dfB <- dfB[dfB$local_optimum >= 1, ]
        
        if(nrow(dfB) > 0){
        
          #we need to adjust for some noise where there are less than 5 mutations on one state, these mutations get removed
          optima_count <- data.frame(table(dfB$local_optimum))
          colnames(optima_count) <- c("local_optimum", "muts_per_optimum")
          optima_count$local_optimum <- as.numeric(as.character(optima_count$local_optimum))
          dfB <- left_join(dfB, optima_count)
          dfB <- dfB[dfB$muts_per_optimum >= 5, ]
          dfB$muts_per_optimum <- NULL
          
          counts <- data.frame(table(dfB$local_optimum, dfB$CLUSTER))
          counts <- counts[counts$Freq > 0, ]
          
          interesting_cluster <- names(table(counts$Var2)[table(counts$Var2) > 1])
          
          # for this to be a real gain the mut cpn must be at 1 and at an integer above 1
          interesting_cluster_muts <- dfB[dfB$CLUSTER %in% interesting_cluster, ]
          
        } else {
          interesting_cluster_muts <- data.frame(matrix(NA, 1, 1)) 
          interesting_cluster_muts <- na.omit(interesting_cluster_muts)
        }
        
        if(nrow(interesting_cluster_muts) > 0){
        
          #let's eliminate those clusters where this is not the case
          cpn_cl_states <- data.frame(table(interesting_cluster_muts$local_optimum, interesting_cluster_muts$CLUSTER))
          cpn_cl_states$Var1 <- as.numeric(as.character(cpn_cl_states$Var1))
          cpn_cl_states <- cpn_cl_states[cpn_cl_states$Freq > 0 & cpn_cl_states$Var1 >= 1, ]
          
          usable_cluster_muts <- interesting_cluster_muts[interesting_cluster_muts$CLUSTER %in% unique(cpn_cl_states$Var2), ]
          
          if(nrow(usable_cluster_muts) > 0){
          
            #let's try to create a flag that says whether a gain is only seen on one cluster or on several
            optimum_cluster <- data.frame(table(usable_cluster_muts$CLUSTER, usable_cluster_muts$local_optimum))
            colnames(optimum_cluster) <- c("CLUSTER", "local_optimum", "Freq")
            optimum_cluster <- optimum_cluster[optimum_cluster$Freq > 0, ]
            optimum_cluster <- optimum_cluster[optimum_cluster$local_optimum != 0, ]
            
            # we only want to keep the gain in the clonal cluster which is not also detected in a subclonal cluster
            trunk <- tree_object$graph_pyclone$trunk
            
            optimum_cluster[optimum_cluster$CLUSTER == trunk, "trunk"] <- "trunk"
            optimum_cluster[optimum_cluster$CLUSTER != trunk, "trunk"] <- "not_trunk"
            
            optimum_df <- data.frame(table(optimum_cluster$local_optimum))
            colnames(optimum_df) <- c("local_optimum", "Freq")
            optimum_df <- optimum_df[optimum_df$Freq != 0, ]
            optimum_to_use <- optimum_df[optimum_df$Freq == 1, ]
            optimum_to_use <- left_join(optimum_to_use, optimum_cluster[, c("CLUSTER", "local_optimum", "trunk")], by = c("local_optimum" = "local_optimum"))
            
            #now if this has the trunk in it we can use it
            optimum_to_use <- optimum_to_use[optimum_to_use$trunk == "trunk", ]
            
            
            #now for each cluster separately we can time the gains
            if(nrow(optimum_to_use) > 0){
            
              cluster_gain_time_df <- list()
              for(cluster in unique(optimum_to_use$CLUSTER)){
                
                df <- usable_cluster_muts[usable_cluster_muts$CLUSTER == cluster, ]
                
                #in theory I should not have any muts at 0
                #but this might be the case due to some noise and then rounding down
                # lets make all those mutations into 1
                df[df$local_optimum < 1, "local_optimum"] <- 1
                
                #only keep the first gain
                cl_optimum <- optimum_to_use$local_optimum
                
                # late_optimum <- as.numeric(as.character(unique(optimum_cluster[optimum_cluster$CLUSTER == cluster,"local_optimum"])))
                # late_optimum <- sort(late_optimum, decreasing = F)
                # late_optimum <- late_optimum[length(late_optimum)-1]
                
                optima_to_use <- c(as.numeric(as.character(cl_optimum)))
                # optima_to_use <- c(as.numeric(as.character(cl_optimum)), late_optimum)
                
                # how many copy number steps are there
                number_steps <- length(optima_to_use)
                
                #for each copy number state count the numbers of mutations
                df$mutID <- paste0(df$chr, ":", df$pos)
                df <- df %>% group_by(local_optimum) %>% mutate(num_muts_on_cn_step = n())
                df <- data.frame(df)
                
                #only keep relevant info
                mutation_count_df <- unique(df[, c("sample", "chr", "seg_start", "seg_end", "phase", "maj_cn_gain", "min_cn_gain", "local_optimum", "num_muts_on_cn_step")])
                mutation_count_df <- mutation_count_df[order(mutation_count_df$local_optimum, decreasing = F), ]
                
                mutation_count_df$copy_number_before_gain <- mutation_count_df$local_optimum
                
                for(s in 1:nrow(mutation_count_df)){
                  mutation_count_df[s, "copy_number_after_gain"] <- ifelse(s < nrow(mutation_count_df), mutation_count_df[s+1, "copy_number_before_gain"], mutation_count_df[s, "copy_number_before_gain"])
                }
                
                mutation_count_df$adjustment_factor <- mutation_count_df$copy_number_before_gain
                
                #if adjustment factor isbelow one make it into one, although this should not happen anymore
                mutation_count_df[mutation_count_df$adjustment_factor < 1, "adjustment_factor"] <- 1
                
                mutation_count_df$number_muts_cn_adjusted <- mutation_count_df$num_muts_on_cn_step / mutation_count_df$adjustment_factor
                mutation_count_df <- mutation_count_df[mutation_count_df$local_optimum %in% optima_to_use, ]
                
                gain_time_df <- data.frame(tumour = NA,
                                           region = NA,
                                           segID = NA,
                                           gain_index = NA,
                                           copy_number_before_gain = NA,
                                           copy_number_after_gain = NA,
                                           number_early_muts_cn_adjusted = NA,
                                           number_all_muts_cn_adjusted = NA,
                                           gain_timing = NA,
                                           loCI = NA,
                                           hiCI = NA,
                                           haplotype = NA) 
                
                for(event in 1:(number_steps - 1)){
                  
                  gain_index                            <- event
                  copy_number_before_gain               <- mutation_count_df[event, "copy_number_before_gain"]
                  copy_number_after_gain                <- mutation_count_df[event, "copy_number_after_gain"]
                  
                  number_early_muts                     <- sum(mutation_count_df[(event + 1) : max(number_steps), "number_muts_cn_adjusted"])
                  number_all_muts                       <- sum(mutation_count_df[, "number_muts_cn_adjusted"])
                  
                  gain_timing                           <- number_early_muts / number_all_muts
                  error                                 <- qnorm(0.975) * sqrt(gain_timing * (1 - gain_timing) / number_all_muts)
                  loCI                                  <- gain_timing - error
                  hiCI                                  <- gain_timing + error
                  
                  gain_time_df[event, "tumour"]                                 <- tumour
                  gain_time_df[event, "region"]                                 <- region
                  gain_time_df[event, "segID"]                                  <- seg
                  gain_time_df[event, "cluster"]                                <- cluster
                  gain_time_df[event, "gain_index"]                             <- number_steps-gain_index
                  gain_time_df[event, "copy_number_before_gain"]                <- copy_number_before_gain
                  gain_time_df[event, "copy_number_after_gain"]                 <- copy_number_after_gain
                  gain_time_df[event, "number_early_muts_cn_adjusted"]          <- number_early_muts
                  gain_time_df[event, "number_all_muts_cn_adjusted"]            <- number_all_muts
                  gain_time_df[event, "gain_timing"]                            <- gain_timing
                  gain_time_df[event, "loCI"]                                   <- loCI
                  gain_time_df[event, "hiCI"]                                   <- hiCI
                  gain_time_df[event, "haplotype"]                              <- "B"
                  
         
                }  
                
                cluster_gain_time_df[[cluster]] <- gain_time_df
                
              }
              
              cluster_gain_time_df <- do.call(rbind, cluster_gain_time_df)
              
              haplotypeB_gaintiming <- cluster_gain_time_df
            
              } else {
              
                gain_time_df <- data.frame(tumour = tumour,
                                           region = region,
                                           segID = seg,
                                           cluster = NA,
                                           gain_index = NA,
                                           copy_number_before_gain = NA,
                                           copy_number_after_gain = NA,
                                           number_early_muts_cn_adjusted = NA,
                                           number_all_muts_cn_adjusted = NA,
                                           gain_timing = NA,
                                           loCI = NA,
                                           hiCI = NA,
                                           haplotype = "B") 
                
                
                haplotypeB_gaintiming <- gain_time_df
                
              }
              
              } else {
                
                gain_time_df <- data.frame(tumour = tumour,
                                           region = region,
                                           segID = seg,
                                           cluster = NA,
                                           gain_index = NA,
                                           copy_number_before_gain = NA,
                                           copy_number_after_gain = NA,
                                           number_early_muts_cn_adjusted = NA,
                                           number_all_muts_cn_adjusted = NA,
                                           gain_timing = NA,
                                           loCI = NA,
                                           hiCI = NA,
                                           haplotype = "B") 
                
                
                haplotypeB_gaintiming <- gain_time_df
                
              }
              } else {
                gain_time_df <- data.frame(tumour = tumour,
                                           region = region,
                                           segID = seg,
                                           cluster = NA,
                                           gain_index = NA,
                                           copy_number_before_gain = NA,
                                           copy_number_after_gain = NA,
                                           number_early_muts_cn_adjusted = NA,
                                           number_all_muts_cn_adjusted = NA,
                                           gain_timing = NA,
                                           loCI = NA,
                                           hiCI = NA,
                                           haplotype = "B") 
                
                
                haplotypeB_gaintiming <- gain_time_df
              }
          
      } else {
        
        #need to deal with non clonal gains for now
        
        gain_time_df <- data.frame(tumour = tumour,
                                   region = region,
                                   segID = seg,
                                   cluster = NA,
                                   gain_index = NA,
                                   copy_number_before_gain = NA,
                                   copy_number_after_gain = NA,
                                   number_early_muts_cn_adjusted = NA,
                                   number_all_muts_cn_adjusted = NA,
                                   gain_timing = NA,
                                   loCI = NA,
                                   hiCI = NA,
                                   haplotype = "B") 
        
        
        haplotypeB_gaintiming <- gain_time_df
        
        
      }
      
    } else {
      
      #and gotta deal with segments that cant be time cause they aren't a gain
      
      gain_time_df <- data.frame(tumour = tumour,
                                 region = region,
                                 segID = seg,
                                 cluster = NA,
                                 gain_index = NA,
                                 copy_number_before_gain = NA,
                                 copy_number_after_gain = NA,
                                 number_early_muts_cn_adjusted = NA,
                                 number_all_muts_cn_adjusted = NA,
                                 gain_timing = NA,
                                 loCI = NA,
                                 hiCI = NA,
                                 haplotype = "B") 
      
      
      haplotypeB_gaintiming <- gain_time_df
      
      
    }
    
    
    both_haplos <- rbind(haplotypeA_gaintiming, haplotypeB_gaintiming)
    
  })
  
  timed_data <- do.call(rbind, timed_data)
  
  write.table(timed_data, paste0(output_path, tumour, "/", region, "/data/clonal_gain_timing.txt"), col.names = T, row.names = T, sep = "\t", quote = F)
  
  return(timed_data)

}


segment.cluster.timing.example.fun <- function(data = phased_data, cn_data = cn, tree_object){
  
  #for each segment, let's see if we have mutations from several clusters on that segment
  data$segID <- paste0(data$chr, "_", data$seg_start, "_", data$seg_end)
  data <- data %>% group_by(segID) %>% mutate(num_clusters_on_seg = length(unique(CLUSTER)))
  data <- data.frame(data)
  
  #now we need to check if there both early and late mutations within one of the clusters
  # if a segment is on the integer value lets use integers to time mutations as early and late
  # if a segment is not on the integer it is a subclonal gain, this means two things:
  # we should only be using subclonal mutations here
  # and the integer thing won't work, so here use local optima to define early and late
  
  seg = "chr14_47734531_83979663"
  
  timed_data <- lapply(unique(data$segID), function(seg){
    
    seg_df <- data[data$segID == seg, ]
    
    #need to do the next steps per haplotype
    
    dfA <- seg_df[which(seg_df$phase == "A"), ]
    dfB <- seg_df[which(seg_df$phase == "B"), ]
    
    #start with haplotype A
    #first, does this segment have a timeable copy number gain?
    
    if(all(nrow(dfA) > 0 && all(dfA$maj_cn_gain))){
      
      #and for now only use clonal copy number gains
      if(as.numeric(as.character(unique(dfA$min_cn))) == round(as.numeric(as.character(unique(dfA$min_cn))))){
        
        # which clusters have several copy number steps? need to define steps first
        # because this is a clonal gain we can use the simple definition
        dfA$local_optimum <- round(dfA$mut_cpn)
        dfA[dfA$local_optimum > dfA$maj_cn, "local_optimum"] <- unique(dfA$maj_cn)
        
        counts <- data.frame(table(dfA$local_optimum, dfA$CLUSTER))
        counts <- counts[counts$Freq > 0, ]
        
        interesting_cluster <- names(table(counts$Var2)[table(counts$Var2) > 1])
        
        # for this to be a real gain the mut cpn must be at 1 and at an integer above 1
        interesting_cluster_muts <- dfA[dfA$CLUSTER %in% interesting_cluster, ]
        
        #let's eliminate those clusters where this is not the case
        cpn_cl_states <- data.frame(table(interesting_cluster_muts$local_optimum, interesting_cluster_muts$CLUSTER))
        cpn_cl_states$Var1 <- as.numeric(as.character(cpn_cl_states$Var1))
        cpn_cl_states <- cpn_cl_states[cpn_cl_states$Freq > 0 & cpn_cl_states$Var1 > 1, ]
        
        usable_cluster_muts <- interesting_cluster_muts[interesting_cluster_muts$CLUSTER %in% unique(cpn_cl_states$Var2), ]
        
        #let's try to create a flag that says whether a gain is only seen on one cluster or on several
        optimum_cluster <- data.frame(table(usable_cluster_muts$CLUSTER, usable_cluster_muts$local_optimum))
        colnames(optimum_cluster) <- c("CLUSTER", "local_optimum", "Freq")
        optimum_cluster <- optimum_cluster[optimum_cluster$Freq > 0, ]
        optimum_cluster <- optimum_cluster[optimum_cluster$local_optimum != 0, ]
        optimum_cluster$local_optimum <- as.numeric(as.character(optimum_cluster$local_optimum))
        
        gain_index_df <- data.frame(local_optimum = sort(unique(as.numeric(as.character(optimum_cluster$local_optimum))), decreasing = T),
                                    gain_index = c(1:length(unique(optimum_cluster$local_optimum))))
        
        gain_index_df <- gain_index_df[gain_index_df$local_optimum != 1, ]
        
        optimum_cluster <- left_join(optimum_cluster, gain_index_df)
        
        cluster_gain_time_df <- list()
        for(cluster in unique(usable_cluster_muts$CLUSTER)){
          
          df <- usable_cluster_muts[usable_cluster_muts$CLUSTER == cluster, ]
          
          #in theory I should not have any muts at 0
          #but this might be the case due to some noise and then rounding down
          # lets make all those mutations into 1
          df[df$local_optimum < 1, "local_optimum"] <- 1
          
          # how many copy number steps are there which are truly just in this cluster
          optima_to_use <- optimum_cluster
          
          #which optima in this cluster is this highest !!! this will not work for all scenarios!!! Attention attention !!!
          optima_to_use <- optima_to_use[optima_to_use$CLUSTER == cluster, ]
          optima_to_use <- optima_to_use[optima_to_use$local_optimum == max(optima_to_use$local_optimum), ]
          optima_to_use <- optima_to_use$local_optimum
          
          #now i need to figure out which one is the step one below which will be the late muts
          # late_optimum <- unique(optimum_cluster$local_optimum)
          # late_optimum <- late_optimum[late_optimum[which(late_optimum %in% optima_to_use)]-1]
          
          # optima_to_use <- c(optima_to_use, late_optimum)
          optima_to_use <- c(optima_to_use)
          
          number_steps <- length(optima_to_use)
          
          #for each copy number state count the numbers of mutations
          df$mutID <- paste0(df$chr, ":", df$pos)
          df <- df %>% group_by(local_optimum) %>% mutate(num_muts_on_cn_step = n())
          df <- data.frame(df)
          
          #only keep relevant info
          mutation_count_df <- unique(df[, c("sample", "chr", "seg_start", "seg_end", "phase", "maj_cn_gain", "min_cn_gain", "local_optimum", "num_muts_on_cn_step")])
          mutation_count_df <- mutation_count_df[order(mutation_count_df$local_optimum, decreasing = F), ]
          
          mutation_count_df$copy_number_before_gain <- mutation_count_df$local_optimum
          
          for(s in 1:nrow(mutation_count_df)){
            mutation_count_df[s, "copy_number_after_gain"] <- ifelse(s < nrow(mutation_count_df), mutation_count_df[s+1, "copy_number_before_gain"], mutation_count_df[s, "copy_number_before_gain"])
          }
          
          mutation_count_df$adjustment_factor <- mutation_count_df$copy_number_before_gain
          
          #if adjustment factor isbelow one make it into one, although this should not happen anymore
          mutation_count_df[mutation_count_df$adjustment_factor < 1, "adjustment_factor"] <- 1
          
          mutation_count_df$number_muts_cn_adjusted <- mutation_count_df$num_muts_on_cn_step / mutation_count_df$adjustment_factor
          
          # how many copy number steps are there

          gain_time_df <- data.frame(tumour = NA,
                                     region = NA,
                                     segID = NA,
                                     gain_index = NA,
                                     copy_number_before_gain = NA,
                                     copy_number_after_gain = NA,
                                     number_early_muts_cn_adjusted = NA,
                                     number_all_muts_cn_adjusted = NA,
                                     gain_timing = NA,
                                     loCI = NA,
                                     hiCI = NA,
                                     haplotype = NA) 
          
          for(event in 1:(number_steps - 1)){
            
            gain_index                            <- event
            copy_number_before_gain               <- mutation_count_df[event, "copy_number_before_gain"]
            copy_number_after_gain                <- mutation_count_df[event, "copy_number_after_gain"]
            
            number_early_muts                     <- sum(mutation_count_df[event + 1, "number_muts_cn_adjusted"])
            number_all_muts                       <- sum(mutation_count_df[, "number_muts_cn_adjusted"])
            
            gain_timing                           <- number_early_muts / number_all_muts
            error                                 <- qnorm(0.975) * sqrt(gain_timing * (1 - gain_timing) / number_all_muts)
            loCI                                  <- gain_timing - error
            hiCI                                  <- gain_timing + error
            
            gain_time_df[event, "tumour"]                                 <- tumour
            gain_time_df[event, "region"]                                 <- region
            gain_time_df[event, "segID"]                                  <- seg
            gain_time_df[event, "cluster"]                                <- cluster
            gain_time_df[event, "gain_index"]                             <- gain_index
            gain_time_df[event, "copy_number_before_gain"]                <- copy_number_before_gain
            gain_time_df[event, "copy_number_after_gain"]                 <- copy_number_after_gain
            gain_time_df[event, "number_early_muts_cn_adjusted"]          <- number_early_muts
            gain_time_df[event, "number_all_muts_cn_adjusted"]            <- number_all_muts
            gain_time_df[event, "gain_timing"]                            <- gain_timing
            gain_time_df[event, "loCI"]                                   <- loCI
            gain_time_df[event, "hiCI"]                                   <- hiCI
            gain_time_df[event, "haplotype"]                              <- "A"
            
          }  
          
          cluster_gain_time_df[[cluster]] <- gain_time_df
          
        }
        
        cluster_gain_time_df <- do.call(rbind, cluster_gain_time_df)
        
        haplotypeA_gaintiming <- cluster_gain_time_df
        
      } else {
        
        #need to deal with non clonal gains for now
        
        gain_time_df <- data.frame(tumour = tumour,
                                   region = region,
                                   segID = seg,
                                   cluster = NA,
                                   gain_index = NA,
                                   copy_number_before_gain = NA,
                                   copy_number_after_gain = NA,
                                   number_early_muts_cn_adjusted = NA,
                                   number_all_muts_cn_adjusted = NA,
                                   gain_timing = NA,
                                   loCI = NA,
                                   hiCI = NA,
                                   haplotype = "A") 
        
        
        haplotypeA_gaintiming <- gain_time_df
        
      }
      
    } else {
      
      #and gotta deal with segments that cant be time cause they aren't a gain
      
      gain_time_df <- data.frame(tumour = tumour,
                                 region = region,
                                 segID = seg,
                                 cluster = NA,
                                 gain_index = NA,
                                 copy_number_before_gain = NA,
                                 copy_number_after_gain = NA,
                                 number_early_muts_cn_adjusted = NA,
                                 number_all_muts_cn_adjusted = NA,
                                 gain_timing = NA,
                                 loCI = NA,
                                 hiCI = NA,
                                 haplotype = "A") 
      
      
      haplotypeA_gaintiming <- gain_time_df
      
      
    }
    
    
    #now do haplotype B
    
    if(all(nrow(dfB) > 0 && all(dfB$min_cn_gain))){
      
      #and for now only use clonal copy number gains
      if(as.numeric(as.character(unique(dfB$min_cn))) == round(as.numeric(as.character(unique(dfB$min_cn))))){
        
        # which clusters have several copy number steps? need to define steps first
        # because this is a clonal gain we can use the simple definition
        dfB$local_optimum <- round(dfB$mut_cpn)
        dfB[dfB$local_optimum > dfB$min_cn, "local_optimum"] <- unique(dfB$maj_cn)
        
        counts <- data.frame(table(dfB$local_optimum, dfB$CLUSTER))
        counts <- counts[counts$Freq > 0, ]
        
        interesting_cluster <- names(table(counts$Var2)[table(counts$Var2) > 1])
        
        # for this to be a real gain the mut cpn must be at 1 and at an integer above 1
        interesting_cluster_muts <- dfB[dfB$CLUSTER %in% interesting_cluster, ]
        
        #let's eliminate those clusters where this is not the case
        cpn_cl_states <- data.frame(table(interesting_cluster_muts$local_optimum, interesting_cluster_muts$CLUSTER))
        cpn_cl_states$Var1 <- as.numeric(as.character(cpn_cl_states$Var1))
        cpn_cl_states <- cpn_cl_states[cpn_cl_states$Freq > 0 & cpn_cl_states$Var1 > 1, ]
        
        usable_cluster_muts <- interesting_cluster_muts[interesting_cluster_muts$CLUSTER %in% unique(cpn_cl_states$Var2), ]
        
        #let's try to create a flag that says whether a gain is only seen on one cluster or on several
        optimum_cluster <- data.frame(table(usable_cluster_muts$CLUSTER, usable_cluster_muts$local_optimum))
        colnames(optimum_cluster) <- c("CLUSTER", "local_optimum", "Freq")
        optimum_cluster <- optimum_cluster[optimum_cluster$Freq > 0, ]
        optimum_cluster <- optimum_cluster[optimum_cluster$local_optimum != 0, ]
        optimum_cluster$local_optimum <- as.numeric(as.character(optimum_cluster$local_optimum))
        
        gain_index_df <- data.frame(local_optimum = sort(unique(as.numeric(as.character(optimum_cluster$local_optimum))), decreasing = T),
                                    gain_index = c(1:length(unique(optimum_cluster$local_optimum))))
        
        gain_index_df <- gain_index_df[gain_index_df$local_optimum != 1, ]
        
        optimum_cluster <- left_join(optimum_cluster, gain_index_df)
        
        cluster_gain_time_df <- list()
        for(cluster in unique(usable_cluster_muts$CLUSTER)){
          
          df <- usable_cluster_muts[usable_cluster_muts$CLUSTER == cluster, ]
          
          #in theory I should not have any muts at 0
          #but this might be the case due to some noise and then rounding down
          # lets make all those mutations into 1
          df[df$local_optimum < 1, "local_optimum"] <- 1
          
          # how many copy number steps are there which are truly just in this cluster
          optima_to_use <- optimum_cluster
          
          #which optima in this cluster is this highest !!! this will not work for all scenarios!!! Attention attention !!!
          optima_to_use <- optima_to_use[optima_to_use$CLUSTER == cluster, ]
          optima_to_use <- optima_to_use[optima_to_use$local_optimum == max(optima_to_use$local_optimum), ]
          optima_to_use <- optima_to_use$local_optimum
          
          #now i need to figure out which one is the step one below which will be the late muts
          # late_optimum <- unique(optimum_cluster$local_optimum)
          # late_optimum <- late_optimum[late_optimum[which(late_optimum %in% optima_to_use)]-1]
          
          # optima_to_use <- c(optima_to_use, late_optimum)
          optima_to_use <- c(optima_to_use)
          
          # how many copy number steps are there
          number_steps <- length(unique(df$local_optimum))
          
          #for each copy number state count the numbers of mutations
          df$mutID <- paste0(df$chr, ":", df$pos)
          df <- df %>% group_by(local_optimum) %>% mutate(num_muts_on_cn_step = n())
          df <- data.frame(df)
          
          #only keep relevant info
          mutation_count_df <- unique(df[, c("sample", "chr", "seg_start", "seg_end", "phase", "maj_cn_gain", "min_cn_gain", "local_optimum", "num_muts_on_cn_step")])
          mutation_count_df <- mutation_count_df[order(mutation_count_df$local_optimum, decreasing = F), ]
          
          mutation_count_df$copy_number_before_gain <- mutation_count_df$local_optimum
          
          for(s in 1:nrow(mutation_count_df)){
            mutation_count_df[s, "copy_number_after_gain"] <- ifelse(s < nrow(mutation_count_df), mutation_count_df[s+1, "copy_number_before_gain"], mutation_count_df[s, "copy_number_before_gain"])
          }
          
          mutation_count_df$adjustment_factor <- mutation_count_df$copy_number_before_gain
          
          #if adjustment factor isbelow one make it into one, although this should not happen anymore
          mutation_count_df[mutation_count_df$adjustment_factor < 1, "adjustment_factor"] <- 1
          
          mutation_count_df$number_muts_cn_adjusted <- mutation_count_df$num_muts_on_cn_step / mutation_count_df$adjustment_factor
          
          gain_time_df <- data.frame(tumour = NA,
                                     region = NA,
                                     segID = NA,
                                     gain_index = NA,
                                     copy_number_before_gain = NA,
                                     copy_number_after_gain = NA,
                                     number_early_muts_cn_adjusted = NA,
                                     number_all_muts_cn_adjusted = NA,
                                     gain_timing = NA,
                                     loCI = NA,
                                     hiCI = NA,
                                     haplotype = NA) 
          
          for(event in 1:(number_steps - 1)){
            
            gain_index                            <- event
            copy_number_before_gain               <- mutation_count_df[event, "copy_number_before_gain"]
            copy_number_after_gain                <- mutation_count_df[event, "copy_number_after_gain"]
            
            number_early_muts                     <- sum(mutation_count_df[event + 1, "number_muts_cn_adjusted"])
            number_all_muts                       <- sum(mutation_count_df[, "number_muts_cn_adjusted"])
            
            gain_timing                           <- number_early_muts / number_all_muts
            error                                 <- qnorm(0.975) * sqrt(gain_timing * (1 - gain_timing) / number_all_muts)
            loCI                                  <- gain_timing - error
            hiCI                                  <- gain_timing + error
            
            gain_time_df[event, "tumour"]                                 <- tumour
            gain_time_df[event, "region"]                                 <- region
            gain_time_df[event, "segID"]                                  <- seg
            gain_time_df[event, "cluster"]                                <- cluster
            gain_time_df[event, "gain_index"]                             <- gain_index
            gain_time_df[event, "copy_number_before_gain"]                <- copy_number_before_gain
            gain_time_df[event, "copy_number_after_gain"]                 <- copy_number_after_gain
            gain_time_df[event, "number_early_muts_cn_adjusted"]          <- number_early_muts
            gain_time_df[event, "number_all_muts_cn_adjusted"]            <- number_all_muts
            gain_time_df[event, "gain_timing"]                            <- gain_timing
            gain_time_df[event, "loCI"]                                   <- loCI
            gain_time_df[event, "hiCI"]                                   <- hiCI
            gain_time_df[event, "haplotype"]                              <- "B"
            
            
          }  
          
          cluster_gain_time_df[[cluster]] <- gain_time_df
          
        }
        
        cluster_gain_time_df <- do.call(rbind, cluster_gain_time_df)
        
        haplotypeB_gaintiming <- cluster_gain_time_df
        
        # if there are subseqwuent gains on non-clonal clusters we've now timed them within the cluster
        # but it would be nice to also work out their timing in terms of the start of tumour formation
        
        # tree_structure <- tree_object$graph_pyclone$default_tree
        # clonal_cluster <- tree_object$graph_pyclone$trunk
        # 
        # #for the clonal cluster everything stays as is
        # haplotypeB_gaintiming[haplotypeB_gaintiming$cluster == clonal_cluster, "tree_time"]
        
        
      } else {
        
        #need to deal with non clonal gains for now
        
        gain_time_df <- data.frame(tumour = tumour,
                                   region = region,
                                   segID = seg,
                                   cluster = NA,
                                   gain_index = NA,
                                   copy_number_before_gain = NA,
                                   copy_number_after_gain = NA,
                                   number_early_muts_cn_adjusted = NA,
                                   number_all_muts_cn_adjusted = NA,
                                   gain_timing = NA,
                                   loCI = NA,
                                   hiCI = NA,
                                   haplotype = "B") 
        
        
        haplotypeB_gaintiming <- gain_time_df
        
        
      }
      
    } else {
      
      #and gotta deal with segments that cant be time cause they aren't a gain
      
      gain_time_df <- data.frame(tumour = tumour,
                                 region = region,
                                 segID = seg,
                                 cluster = NA,
                                 gain_index = NA,
                                 copy_number_before_gain = NA,
                                 copy_number_after_gain = NA,
                                 number_early_muts_cn_adjusted = NA,
                                 number_all_muts_cn_adjusted = NA,
                                 gain_timing = NA,
                                 loCI = NA,
                                 hiCI = NA,
                                 haplotype = "B") 
      
      
      haplotypeB_gaintiming <- gain_time_df
      
      
    }
    
    
    both_haplos <- rbind(haplotypeA_gaintiming, haplotypeB_gaintiming)
    
    timing.tree.fun(tumour_id = tumour, data = both_haplos, fullTree = tree_object, output_path, region)
    
  })
  
  timed_data <- do.call(rbind, timed_data)
  return(timed_data)
  
}










# functions to create segment timing on tree plots

getDataForTree <- function(tree, trunk, node.fill = "gray", node.shape = 21, node.size = 13, stroke.size = 1) {
  tree <- rbind(c("0", trunk), tree)
  
  g.tree    <- igraph::as.undirected(igraph::graph_from_data_frame(tree))
  node.indx <- igraph::V(g.tree)$name
  l.tree    <- igraph::layout_as_tree(g.tree, root = '0', flip.y = FALSE)
  
  tree.plot.nodes <- data.frame(Node = node.indx[-1], 
                                x = l.tree[-1, 1], 
                                y = l.tree[-1, 2],
                                size = node.size,
                                shape = node.shape,
                                fill = node.fill,
                                color = "black",
                                stroke = stroke.size,
                                node.label = "", 
                                stringsAsFactors = FALSE)
  
  
  tree.plot.edges <- data.frame(Edge = apply(tree, 1, function(x) paste0(x, collapse = "-")),
                                xstart = l.tree[match(tree[,1], node.indx), 1],
                                xend = l.tree[match(tree[,2], node.indx), 1],
                                ystart = l.tree[match(tree[,1], node.indx), 2],
                                yend = l.tree[match(tree[,2], node.indx), 2],
                                color = "black",
                                size = stroke.size, 
                                stringsAsFactors = FALSE)
  tree.plot.edges$edge.label <- ""
  tree.plot.edges$label.x <- (tree.plot.edges$xstart + tree.plot.edges$xend) / 2
  tree.plot.edges$label.y <- (tree.plot.edges$ystart + tree.plot.edges$yend) / 2
  tree.plot.edges$highlight.color <- "transparent"
  
  return(list(tree.plot.nodes = tree.plot.nodes, tree.plot.edges = tree.plot.edges))
}


plottingTxTree <- function(tree.plot.nodes, tree.plot.edges, tree.title = NULL, tree.subtitle = NULL, node.label.type = "none", add.edge.labels = FALSE, edge.highlights = NULL) {
  p <- ggplot()
  
  if (!is.null(edge.highlights)) {
    p <- p + geom_segment(data = tree.plot.edges, aes(x = xstart, y = ystart, xend = xend, yend = yend), color = tree.plot.edges$highlight.color, linewidth = tree.plot.edges$size * 5)
  }
  
  p <- p + geom_segment(data = tree.plot.edges, aes(x = xstart, y = ystart, xend = xend, yend = yend), color = tree.plot.edges$color, linewidth = tree.plot.edges$size) +
    geom_point(data = tree.plot.nodes, aes(x = x, y = y), fill = "white", color = "black", size = tree.plot.nodes$size, shape = tree.plot.nodes$shape, stroke = tree.plot.nodes$stroke) +
    geom_point(data = tree.plot.nodes, aes(x = x, y = y), fill = tree.plot.nodes$fill, color = tree.plot.nodes$color, size = tree.plot.nodes$size, shape = tree.plot.nodes$shape, stroke = tree.plot.nodes$stroke) +
    theme_void()
  
  if (!is.null(tree.title)) {
    p <- p + ggtitle(tree.title) + theme(plot.title = element_text(hjust = 0.5))
  }
  
  if (!is.null(tree.subtitle)) {
    p <- p + labs(subtitle = tree.subtitle) + theme(plot.subtitle = element_text(hjust = 0.5))
  }
  
  if (add.edge.labels) {
    p <- p + geom_text(data = tree.plot.edges, aes(x = label.x, y = label.y, label = stringr::str_wrap(edge.label, 20)), hjust = 0, nudge_x = 0.05, vjust = 1, size = 3)
    # p <- p + geom_text(data = tree.plot.edges, aes(x = label.x, y = label.y, label = stringr::str_wrap(label, 20)), hjust = 0.5, vjust = 0.5, size = 20)
  }
  
  if (node.label.type == "none") {
    return(p)
  } else if (node.label.type == "nodeNumber") {
    p <- p + geom_text(data = tree.plot.nodes, aes(x = x, y = y, label = Node), size = 7)
    return(p)
  } else if (node.label.type == "custom") {
    p <- p + geom_point(data = tree.plot.nodes %>% filter(node.label != ""), aes(x = x, y = y), fill = "white", color = "white", shape = 24, size = tree.plot.nodes %>% filter(node.label != "") %>% pull(size) * 0.5)
    return(p)
  } else {
    warning(paste0("Node label type ", node.label.type, " unknown. Returning basic plot."))
    return(p)
  }
}

plottingTimepointTreeNodeHighlight <- function(tumorID, tree, trunk, edgelength, colorTable, time_df, region, haplotype) {
  tmp <- getDataForTree(tree, trunk)
  tree.plot.nodes <- tmp$tree.plot.nodes
  tree.plot.edges <- tmp$tree.plot.edges
  rm(tmp)
  
  #sort out the labels
  tree.plot.edges$cluster <- sapply(tree.plot.edges$Edge, function(x){strsplit(x, "-")[[1]][2]})
  tree.plot.edges <- full_join(tree.plot.edges, time_df, multiple = "all")
  
  tree.plot.edges[!is.na(tree.plot.edges$gain_timing), "label.x"] <- (tree.plot.edges[!is.na(tree.plot.edges$gain_timing), "xend"] - tree.plot.edges[!is.na(tree.plot.edges$gain_timing), "xstart"]) * tree.plot.edges[!is.na(tree.plot.edges$gain_timing), "gain_timing"] + tree.plot.edges[!is.na(tree.plot.edges$gain_timing), "xstart"]
  tree.plot.edges[!is.na(tree.plot.edges$gain_timing), "label.y"] <- (tree.plot.edges[!is.na(tree.plot.edges$gain_timing), "yend"] - tree.plot.edges[!is.na(tree.plot.edges$gain_timing), "ystart"]) * tree.plot.edges[!is.na(tree.plot.edges$gain_timing), "gain_timing"] + tree.plot.edges[!is.na(tree.plot.edges$gain_timing), "ystart"]
  tree.plot.edges[!is.na(tree.plot.edges$gain_timing), "edge.label"] <- paste0("--", round(tree.plot.edges[!is.na(tree.plot.edges$gain_timing), "gain_timing"] * 100, 1), "%")
  tree.plot.edges$cluster <- NULL
  tree.plot.edges$gain_timing <- NULL
  
  plot.title.label <- paste0(region, " - ", strsplit(seg, "_")[[1]][1], ":", strsplit(seg, "_")[[1]][2], ":", strsplit(seg, "_")[[1]][3], " - haplotype ", haplotype)
  colorTable <- colorTable %>% dplyr::filter(tumour_id == tumorID) %>% dplyr::rename(color = colors)
  
  
  colorTable <- colorTable %>% mutate(line.color = "black") %>% mutate(node.size = 13)
  
  
  tree.plot.nodes <- tree.plot.nodes %>% 
    select(-fill, -color, -size, -stroke) %>% 
    left_join(colorTable %>% 
                select(Node = clones, fill = color, color = line.color, size = node.size, stroke) %>% 
                unique(), by = "Node") %>% 
    select(Node, x, y, size, shape, fill, color, stroke, node.label)
  
  tree.plot.edges <- tree.plot.edges %>% 
    mutate(daughter = sapply(strsplit(Edge, split = "-"), function(x) x[2])) %>%
    select(-color) %>%
    left_join(colorTable %>% 
                select(daughter = clones, color = line.color) %>%
                unique(), by = "daughter") %>%
    select(Edge, xstart, xend, ystart, yend, color, size, edge.label, label.x, label.y, highlight.color)
  
  
  p <- plottingTxTree(tree.plot.nodes, 
                      tree.plot.edges, 
                      tree.title = plot.title.label, 
                      tree.subtitle = NULL, 
                      node.label.type = "nodeNumber", 
                      add.edge.labels = TRUE, 
                      edge.highlights = TRUE)
  return(p)
}

timing.tree.fun <- function(tumour_id = tumour, data = phased_timed_data, fullTree = tree_object, output_path, region){
  
  plot_list <- list()
  for(seg in unique(data$segID)){
    
    df <- data[data$segID == seg, ]
    df$cluster <- as.character(df$cluster)
    
    seg_time_dfA <- df[df$haplotype == "A", c("cluster", "gain_timing")]
    seg_time_dfB <- df[df$haplotype == "B", c("cluster", "gain_timing")]
    
    cloneColorDF <- data.frame(tumour_id = tumour_id,
                               clones = as.character(sort(unique(as.numeric(fullTree$graph_pyclone$default_tree)))),
                               colors = "gray",
                               stroke = 1)
    
    cloneColorDF$colors <- brewer.pal(12, "Paired")[as.numeric(cloneColorDF$clones)]
    
    plotA <- plottingTimepointTreeNodeHighlight(tumorID = tumour_id, tree = fullTree$graph_pyclone$default_tree, trunk = fullTree$graph_pyclone$trunk, edgelength = fullTree$graph_pyclone$edgelength, colorTable = cloneColorDF, time_df = seg_time_dfA, region = region, haplotype = "A") + scale_y_reverse()
    plotA <- plotA + theme(plot.margin = unit(c(1,2,1,2), "cm"))
    plotB <- plottingTimepointTreeNodeHighlight(tumorID = tumour_id, tree = fullTree$graph_pyclone$default_tree, trunk = fullTree$graph_pyclone$trunk, edgelength = fullTree$graph_pyclone$edgelength, colorTable = cloneColorDF, time_df = seg_time_dfB, region = region, haplotype = "B") + scale_y_reverse()
    plotB <- plotB + theme(plot.margin = unit(c(1,2,1,2), "cm"))
    
    p <- plotA + plotB
    plot_list[[seg]] <- p
  }
  
  
  if(!dir.exists(paste0(output_path, tumour, "/", region, "/plots/timing_trees/"))){
    dir.create(paste0(output_path, tumour, "/", region, "/plots/timing_trees/"), recursive = T)
  }
  
  pdf(paste0(output_path, tumour, "/", region, "/plots/timing_trees/", "gain_timing_on_tree.pdf"), width = 12, height = 6)
  for(i in 1:length(seg)){
    print(plot_list[[i]])
  }
  dev.off()
  
}


# function that estimates if a mutation is early or late
# mut.cpn.match.fun <- function(mutation_df = phased_data){
#   
#   mutation_df$mut_cpn_round <- round(mutation_df$mut_cpn)
#   
#   # what is the maximum observed copy number for each segment, rounded
#   mutation_df$segID <- paste0(mutation_df$chr, "_", mutation_df$seg_start, "_", mutation_df$seg_end)
#   mutation_df <- mutation_df %>% 
#                     group_by(segID) %>%
#                     rowwise() %>%
#                     mutate(max_CN_A = round(max(maj_cn))) %>%
#                     mutate(max_CN_B = round(max(min_cn)))
#   
#   mutation_df <- data.frame(mutation_df)
#   
#   # which mutations are on the gain and which ones are not
#   # calculate the local maxima for each copy number segment
#   mutation_opt_df <- lapply(unique(mutation_df$segID), function(segID){
#     
#     df <- mutation_df[mutation_df$segID == segID, ]
# 
#     if(!all(is.na(df$phase))){
#       #need to do this seperately for the two haplotypes
#       dfA <- df[df$phase == "A", ]
#       dfB <- df[df$phase == "B", ]
# 
#       if(nrow(dfA) > 1){
#         dns_max <- ifelse(max(dfA$maj_cn, na.rm = T) < min(dfA$mut_cpn, na.rm = T), max(dfA$mut_cpn, na.rm = T), max(dfA$maj_cn, na.rm = T))
#         dns <- density(dfA$mut_cpn,from=min(dfA$mut_cpn)-0.2,to=dns_max+0.2)
#         hypercube.size = 4
#         localOptima = NULL
#         peak.indices = NULL
#         for(i in (1+hypercube.size):(length(dns$y)-hypercube.size)){
#           if(dns$y[i] == max(dns$y[(i-hypercube.size):(i+hypercube.size)])){
#             localOptima = c(localOptima,dns$x[i])
#             peak.indices = c(peak.indices,i)
#           }
#         }
#         
# 
#         # let's remove any local optima where there is no density, set this slightly higher
#         if (length(localOptima)>1){
#           localOptima <- localOptima[which(dns$y[dns$x%in%localOptima]>0.005)]
#           
#            #we need to merge any local optima below one
#           below <- localOptima[localOptima < 1]
#           below <- mean(below)
#           above <- localOptima[localOptima > 1]
#           localOptima <- c(below, above)
#           localOptima <- localOptima[!is.na(localOptima)]
# 
#           #now we need to assign the mutations to the optima based on which one is the closest
#           for(mut in 1:nrow(dfA)){
#             target.index <- which(abs(localOptima - dfA[mut, "mut_cpn"]) == min(abs(localOptima - dfA[mut, "mut_cpn"])))
# 
#             # if target index is longer than 1 then we have more local optima than mutations, this can't be
#             if(length(target.index) > 1){
#               localOptima <- mean(localOptima)
#               target.index <- which(abs(localOptima - dfA[mut, "mut_cpn"]) == min(abs(localOptima - dfB[mut, "mut_cpn"])))
#             }
# 
#             dfA[mut, "local_optimum"] <- localOptima[target.index]
#           }
#           
#           dfA$global_maximum <- max(localOptima)
#           dfA$global_minimum <- min(localOptima)
# 
#         } else {
#           dfA$local_optimum <- NA
#           dfA$global_maximum <- NA
#           dfA$global_minimum <- NA
#         }
# 
#       }
# 
#       if(nrow(dfB) > 1){
#         dns_max <- ifelse(max(dfB$min_cn, na.rm = T) < min(dfB$mut_cpn, na.rm = T), max(dfB$mut_cpn, na.rm = T), max(dfB$min_cn, na.rm = T))
# 
#         dns <- density(dfB$mut_cpn,from=min(dfB$mut_cpn)-0.2,to=dns_max+0.2)
#         hypercube.size = 4
#         localOptima = NULL
#         peak.indices = NULL
#         for(i in (1+hypercube.size):(length(dns$y)-hypercube.size)){
#           if(dns$y[i] == max(dns$y[(i-hypercube.size):(i+hypercube.size)])){
#             localOptima = c(localOptima,dns$x[i])
#             peak.indices = c(peak.indices,i)
#           }
#         }
#         
#         # let's remove any local optima where there is no density, set this slightly higher
#         if (length(localOptima)>1){
#           localOptima <- localOptima[which(dns$y[dns$x%in%localOptima]>0.005)]
# 
#            #we need to merge any local optima below one
#           below <- localOptima[localOptima < 1]
#           below <- mean(below)
#           above <- localOptima[localOptima > 1]
#           localOptima <- c(below, above)
#           localOptima <- localOptima[!is.na(localOptima)]
# 
#           #now we need to assign the mutations to the optima based on which one is the closest
#           for(mut in 1:nrow(dfB)){
#             target.index <- which(abs(localOptima - dfB[mut, "mut_cpn"]) == min(abs(localOptima - dfB[mut, "mut_cpn"])))
# 
#             # if target index is longer than 1 then we have more local optima than mutations, this can't be
#             if(length(target.index) > 1){
#               localOptima <- mean(localOptima)
#               target.index <- which(abs(localOptima - dfB[mut, "mut_cpn"]) == min(abs(localOptima - dfB[mut, "mut_cpn"])))
#             }
# 
#             dfB[mut, "local_optimum"] <- localOptima[target.index]
#           }
#           
#           dfB$global_maximum <- max(localOptima)
#           dfB$global_minimum <- min(localOptima)
#         } else {
#           dfB$local_optimum <- NA
#           dfB$global_maximum <- NA
#           dfB$global_minimum <- NA
#         }
# 
# 
#       }
# 
#       #if only one mut present, lets make this into NA cause we cant do anything about this one lil guy
#       if(nrow(dfA) == 1){
#         
#         dfA$local_optimum <- NA
#         dfA$global_maximum <- NA
#         dfA$global_minimum <- NA
# 
#       }
# 
#       # if only one local optima present let's also resport back to rounding approach
#       if(length(unique(dfA$local_optimum)) == 1 & all(!is.na(unique(dfA$local_optimum)))){
#         
#         dfA$local_optimum <- dfA$mut_cpn_round
#         dfA$global_maximum <- max(dfA$max_CN_A, dfA$mut_cpn_round)
#         dfA$global_minimum <- min(dfA$max_CN_A, dfA$mut_cpn_round)
# 
#       }
# 
# 
#       if(nrow(dfB) == 1){
#         
#         dfB$local_optimum <- NA
#         dfB$global_maximum <- NA
#         dfB$global_minimum <- NA
# 
#       }
# 
#       # if only one local optima present let's also resport back to rounding approach
#       if(length(unique(dfB$local_optimum)) == 1 & all(!is.na(unique(dfB$local_optimum)))){
#         
#         dfB$local_optimum <- dfB$mut_cpn_round
#         dfB$global_maximum <- max(dfB$maj_cn, dfB$min_cn, dfB$mut_cpn_round)
#         dfB$global_minimum <- min(dfB$maj_cn, dfB$min_cn, dfB$mut_cpn_round)
# 
#       }
# 
#       #make a timing suggestion
#       #need to introduce a special clause for gained segments which dont have any mutations
#       # plus a special clause when only one optima exists
#       if(nrow(dfB) > 0){
#         if(!any(dfB$mut_cpn >= max(dfB$max_CN_B)-1)){
#           dfB$timing <- "late"
#           dfB[is.na(dfB$local_optimum), "timing"] <- "not_possible"
#         } else {
#           dfB[which(dfB$local_optimum >= dfB$global_maximum), "timing"] <- "early"
#           dfB[which(dfB$local_optimum < dfB$global_maximum), "timing"] <- "late"
#           dfB[is.na(dfB$local_optimum), "timing"] <- "not_possible"
#         }
#       }
# 
#       if(nrow(dfA) > 0){
#         if(!any(dfA$mut_cpn >= max(dfA$max_CN_A)-1)){
#           dfA$timing <- "late"
#           dfA[is.na(dfA$local_optimum), "timing"] <- "not_possible"
#         } else {
#           dfA[which(dfA$local_optimum >= dfA$global_maximum), "timing"] <- "early"
#           dfA[which(dfA$local_optimum < dfA$global_maximum), "timing"] <- "late"
#           dfA[is.na(dfA$local_optimum), "timing"] <- "not_possible"
#         }
#       }
# 
#       df <- rbind(dfA, dfB)
#       return(df)
#     }
# 
#     if(all(is.na(df$phase))){
# 
#       df$local_optimum  <- NA
#       df$global_maximum <- NA
#       df$global_minimum <- NA
#       df$timing <- "not_possible"
# 
#       return(df)
#     
#     }
#     
#   })
#   
#   mutation_opt_df <- do.call(rbind, mutation_opt_df)
# 
#   return(mutation_opt_df)
#   
# } 

# function to get colors
get.colors <- function(max_cn = max(mutation_df$mut_cpn_round), opacity = 255){
  cols <- paste(brewer.pal(max_cn, name = "Paired"), sep = "")
  nopaq <- round(opacity / ceiling(max_cn / 12))
  opacity_val <- c(sapply(1:ceiling(max_cn / 12), function(x) rep(opacity - nopaq * (x - 1), 12)))[1:max_cn]
  cols.opac <- paste(cols, as.hexmode(opacity_val), sep = "")
}

# function to plot mutations coloured by rounded mutation copy number and seperately for each allele
copynumber.state.plot.fun <- function(mutation_df, opacity = 255, tumour, region, output_path, genome = BSgenome.Hsapiens.UCSC.hg38){
  
  mutation_df$haplotype <- paste0("haplotype", mutation_df$phase)
  mutation_df$local_optimum <- round(mutation_df$mut_cpn)
  mutation_df[mutation_df$phase == "A" & mutation_df$local_optimum > mutation_df$maj_cn, "local_optimum"] <- mutation_df[mutation_df$phase == "A" & mutation_df$local_optimum > mutation_df$maj_cn, "maj_cn"]
  mutation_df[mutation_df$phase == "B" & mutation_df$local_optimum > mutation_df$min_cn, "local_optimum"] <- mutation_df[mutation_df$phase == "B" & mutation_df$local_optimum > mutation_df$min_cn, "min_cn"]
  
  mutation_df$local_optimum <- as.character(mutation_df$local_optimum)
  
  
  #color every other step differently
  mutation_df_seg <- lapply(unique(mutation_df$segID), function(segID){

    df <- mutation_df[mutation_df$segID == segID, ]
    dfA <- df[which(df$phase == "A" & df$timing != "not_possible"), ]
    dfB <- df[which(df$phase == "B" & df$timing != "not_possible"), ]
    dfNA <- df[which(df$timing == "not_possible"), ]

    # how many steps are there in this segment

    if(nrow(dfA) > 0){
      dfA$number_steps <- length(unique(dfA$local_optimum))
      color_indexA <- data.frame(local_optimum = sort(unique(dfA$local_optimum)),
                               color_index = rep(c("one", "two"), length.out = length(unique(dfA$local_optimum))))
      dfA <- suppressMessages(left_join(dfA, color_indexA))
    }
    
    if(nrow(dfB) > 0){
      dfB$number_steps <- length(unique(dfB$local_optimum))
      color_indexB <- data.frame(local_optimum = sort(unique(dfB$local_optimum)),
                                color_index = rep(c("one", "two"), length.out = length(unique(dfB$local_optimum))))
      dfB <- suppressMessages(left_join(dfB, color_indexB))
    }

    if(nrow(dfNA) > 0){
      dfNA$color_index <- "not_possible"
      dfNA$number_steps <- NA
    }
    
    df_out <- rbind(dfA, dfB, dfNA)

  })
  mutation_df_seg <- do.call(rbind, mutation_df_seg)

  mutation_df_seg[mutation_df_seg$maj_cn_gain == FALSE & mutation_df_seg$min_cn_gain == FALSE, "timing"] <- "not_possible"
  mutation_df_seg[which(mutation_df_seg$maj_cn_gain == FALSE & mutation_df_seg$phase == "A"), "timing"] <- "not_possible"
  mutation_df_seg[which(mutation_df_seg$min_cn_gain == FALSE & mutation_df_seg$phase == "B"), "timing"] <- "not_possible"
  mutation_df_seg[is.na(mutation_df_seg$phase), "timing"] <- "not_possible"
  mutation_df_seg$timing <- factor(mutation_df_seg$timing, levels = c("early", "late", "not_possible"))

  mutation_df_seg[mutation_df_seg$timing == "not_possible", "color_index"] <- "not_possible"

  mutation_df_seg$color_index <- factor(mutation_df_seg$color_index, levels = c("one", "two", "not_possible"))
  mutation_df_seg$allelic_imbalance <- factor(mutation_df_seg$allelic_imbalance, levels = c("TRUE", "FALSE"))
  
  mutation_df_seg[mutation_df_seg$mut_cpn > 10, "mut_cpn_lim"] <- 10
  mutation_df_seg[mutation_df_seg$mut_cpn <= 10, "mut_cpn_lim"] <- mutation_df_seg[mutation_df_seg$mut_cpn <= 10, "mut_cpn"] 
  mutation_df_seg[mutation_df_seg$min_cn > 10, "seg_y_pos_min_lim"] <- 10
  mutation_df_seg[mutation_df_seg$min_cn <= 10, "seg_y_pos_min_lim"] <- mutation_df_seg[mutation_df_seg$min_cn <= 10, "min_cn"]
  mutation_df_seg[mutation_df_seg$maj_cn > 10, "seg_y_pos_maj_lim"] <- 10
  mutation_df_seg[mutation_df_seg$maj_cn <= 10, "seg_y_pos_maj_lim"] <- mutation_df_seg[mutation_df_seg$maj_cn <= 10, "maj_cn"]
  
  mut_cpn_ylim_genome <- max(c(mutation_df_seg$mut_cpn_lim, mutation_df_seg$seg_y_pos_maj_lim))
  mut_cpn_ylim_genome_min <- min(c(mutation_df_seg$mut_cpn_lim, mutation_df_seg$seg_y_pos_min_lim))

  # make genome continuous
  mutation_df_cont <- genome.continous.fun(mutation_df_seg, genome)
  chr_sizes <- seqlengths(genome)[paste0("chr", 1:22)]
  chr_lines <- cumsum(as.numeric(as.character(chr_sizes)))
  chr_midpoint <- chr_lines - chr_sizes / 2

  xaxis_label <- data.frame(chr_midpoint)
  xaxis_label$label <- rownames(xaxis_label)
  
  p_genome1 <- ggplot(mutation_df_cont[which(mutation_df_cont$phase == "A"), ], aes(adjust_pos, mut_cpn_lim, fill = color_index), color = "black") +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
    scale_fill_manual(name = "", breaks = c("one", "two", "not_possible"), labels = c("one", "two", "not possible"), values = c("#d8b365", "#5ab4ac", "#bababa"), guide = "none") +
    geom_segment(data = unique(mutation_df_cont[, c("adjust_seg_start", "adjust_seg_end", "allelic_imbalance", "phase", "seg_y_pos_maj_lim", "seg_y_pos_min_lim", "color_index")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_maj_lim, yend = seg_y_pos_maj_lim, color = 'allele A', alpha = allelic_imbalance), linewidth = 2) +
    scale_color_manual(name = "copy number phase", values=c("allele A" = "darkred")) +
    scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    facet_wrap(~haplotype, ncol = 1, strip.position = "right") +
    ylab("mutation copy number") +
    xlab("") +
    ylim(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.background = element_rect(fill = "white"))
  
  p_genome2 <- ggplot(mutation_df_cont[which(mutation_df_cont$phase == "B"), ], aes(adjust_pos, mut_cpn_lim, fill = color_index), color = "black") +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
    scale_fill_manual(name = "", breaks = c("one", "two", "not_possible"), labels = c("one", "two", "not possible"), values = c("#d8b365", "#5ab4ac", "#bababa"), guide = "none") +
    geom_segment(data = unique(mutation_df_cont[, c("adjust_seg_start", "adjust_seg_end", "allelic_imbalance", "phase", "seg_y_pos_maj_lim", "seg_y_pos_min_lim", "color_index")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_min_lim, yend = seg_y_pos_min_lim, color = 'allele B', alpha = allelic_imbalance), linewidth = 2) +
    scale_color_manual(name = "copy number phase", values=c("allele B" = "#313695")) +
    scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    facet_wrap(~haplotype, ncol = 1, strip.position = "right") +
    ylab("mutation copy number") +
    xlab("") +
    ylim(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.background = element_rect(fill = "white"))
  
  
  # get the legend
  theme_margin <- theme(legend.box.margin = margin(0, 0, 0, 0))
  
  legend_genome  <- cowplot::get_legend(p_genome1 + theme_margin)
  
  # remove the legend
  p_genome1 <- p_genome1 + theme(legend.position = "none", plot.margin =  unit(c(0.5, 0.5, 0, 0.5), "cm"))
  p_genome2 <- p_genome2 + theme(legend.position = "none", plot.margin =  unit(c(0, 0.5, 0.5, 0.5), "cm"))
  
  # combine the two plots
  p_genome_c <- plot_grid(p_genome1, p_genome2, ncol = 1, align = "v", axis = "rl")
  
  # add the legend back on
  p_genome_cl <- plot_grid(p_genome_c, legend_genome, nrow = 1, rel_widths = c(10, 1))
  
  
  if(!dir.exists(paste0(output_path, tumour, "/", region, "/plots/phased_mutcpn"))){
    dir.create(paste0(output_path, tumour, "/", region, "/plots/phased_mutcpn"), recursive = T)
  }
  
  pdf(paste0(output_path, tumour, "/", region, "/plots/phased_mutcpn/", region, "_genome_phased_mutation_copynumber_per_haplotype_mutcnp.pdf"), width = 18, height = 7)
  print(p_genome_cl)
  dev.off()
  
  # now plot this for every chromosome
  for(chr in unique(mutation_df$chr)){
  
    mutation_chr_df <- mutation_df_cont[mutation_df_cont$chr == chr, ]
      
    mut_cpn_ylim_chr <- max(c(mutation_chr_df$mut_cpn_lim, mutation_chr_df$seg_y_pos_maj_lim))
    mut_cpn_ylim_chr_min <- min(c(mutation_chr_df$mut_cpn_lim, mutation_chr_df$seg_y_pos_min_lim))

    
    if(nrow(mutation_chr_df[which(mutation_chr_df$phase == "A"), ]) > 0){
      p_chr_genome1 <- ggplot(mutation_chr_df[which(mutation_chr_df$phase == "A"), ], aes(pos, mut_cpn_lim, fill = color_index), color = "black") +
        geom_segment(data = unique(mutation_chr_df[, c("seg_start", "seg_end", "allelic_imbalance", "phase", "seg_y_pos_maj_lim", "seg_y_pos_min_lim", "color_index")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_maj_lim, yend = seg_y_pos_maj_lim, color = 'allele A', alpha = allelic_imbalance), linewidth = 2) +
        scale_color_manual(name = "copy number phase", values=c("allele A" = "darkred")) +
        scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
        geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
        scale_fill_manual(name = "", breaks = c("one", "two", "not_possible"), labels = c("one", "two", "not possible"), values = c("#d8b365", "#5ab4ac", "#bababa"), guide = "none") +
        facet_wrap(~haplotype, ncol = 1, strip.position = "right") +
        ylab("mutation copy number") +
        xlab("") +
        ylim(mut_cpn_ylim_chr_min, mut_cpn_ylim_chr) +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              strip.background = element_rect(fill = "white"))
    }

    if(nrow(mutation_chr_df[which(mutation_chr_df$phase == "B"), ]) > 0){
    p_chr_genome2 <- ggplot(mutation_chr_df[which(mutation_chr_df$phase == "B"), ], aes(pos, mut_cpn_lim, fill = color_index), color = "black") +
      geom_segment(data = unique(mutation_chr_df[, c("seg_start", "seg_end", "allelic_imbalance", "phase", "seg_y_pos_maj_lim", "seg_y_pos_min_lim", "color_index")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_min_lim, yend = seg_y_pos_min_lim, color = 'allele B', alpha = allelic_imbalance), linewidth = 2) +
      scale_color_manual(name = "copy number phase", values=c("allele B" = "#313695")) +
      scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
      geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
      scale_fill_manual(name = "", breaks = c("one", "two", "not_possible"), labels = c("one", "two", "not possible"), values = c("#d8b365", "#5ab4ac", "#bababa"), guide = "none") +
      facet_wrap(~haplotype, ncol = 1, strip.position = "right") +
      ylab("mutation copy number") +
      xlab(chr) +
      ylim(mut_cpn_ylim_chr_min, mut_cpn_ylim_chr) +
      theme_bw() +
      theme(legend.position = "none",
            strip.background = element_rect(fill = "white"))
    }
    
    
    if(all(c(nrow(mutation_chr_df[which(mutation_chr_df$phase == "A"), ]) > 0, nrow(mutation_chr_df[which(mutation_chr_df$phase == "B"), ]) > 0))){
      # get the legend
      theme_margin <- theme(legend.box.margin = margin(0, 0, 0, 0))
      
      legend_chr  <- cowplot::get_legend(p_chr_genome1 + theme_margin)
      
      # remove the legend
      p_chr_genome1 <- p_chr_genome1 + theme(legend.position = "none", plot.margin =  unit(c(0.5, 0.5, 0, 0.5), "cm"))
      p_chr_genome2 <- p_chr_genome2 + theme(legend.position = "none", plot.margin =  unit(c(0, 0.5, 0.5, 0.5), "cm"))
      
      # combine the two plots
      p_chr_c <- plot_grid(p_chr_genome1, p_chr_genome2, ncol = 1, align = "v", axis = "rl")
      
      # add the legend back on
      p_chr_cl <- plot_grid(p_chr_c, legend_chr, nrow = 1, rel_widths = c(10, 1.5))
      
      
      if(!dir.exists(paste0(output_path, tumour, "/", region, "/plots/phased_mutcpn/per_chromosome/"))){
        dir.create(paste0(output_path, tumour, "/", region, "/plots/phased_mutcpn/per_chromosome"), recursive = T)
      }
      
      pdf(paste0(output_path, tumour, "/", region, "/plots/phased_mutcpn/per_chromosome/", region, "_", chr, "_genome_phased_mutation_copynumber_per_haplotype_mutcnp.pdf"), width = 10, height = 7)
      print(p_chr_cl)
      dev.off()
      
    } else {
      
      if(nrow(mutation_chr_df[which(mutation_chr_df$phase == "A"), ]) > 0){
        p_chr_genome1 <- ggplot(mutation_chr_df[which(mutation_chr_df$phase == "A"), ], aes(pos, mut_cpn_lim, fill = color_index), color = "black") +
          geom_segment(data = unique(mutation_chr_df[, c("seg_start", "seg_end", "allelic_imbalance", "phase", "seg_y_pos_maj_lim", "seg_y_pos_min_lim", "color_index")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_maj_lim, yend = seg_y_pos_maj_lim, color = 'allele A', alpha = allelic_imbalance), linewidth = 2) +
          scale_color_manual(name = "copy number phase", values=c("allele A" = "darkred")) +
          scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
          geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
          scale_fill_manual(name = "", breaks = c("one", "two", "not_possible"), labels = c("one", "two", "not possible"), values = c("#d8b365", "#5ab4ac", "#bababa"), guide = "none") +
          facet_wrap(~haplotype, ncol = 1, strip.position = "right") +
          ylab("mutation copy number") +
          xlab(chr) +
          ylim(mut_cpn_ylim_chr_min, mut_cpn_ylim_chr) +
          theme_bw() +
          theme(strip.background = element_rect(fill = "white"))
        
        pdf(paste0(output_path, tumour, "/", region, "/plots/phased_mutcpn/per_chromosome/", region, "_", chr, "_genome_phased_mutation_copynumber_per_haplotype_mutcnp.pdf"), width = 10, height = 7)
        print(p_chr_genome1)
        dev.off()
      }  
      
      if(nrow(mutation_chr_df[which(mutation_chr_df$phase == "B"), ]) > 0){
        p_chr_genome2 <- ggplot(mutation_chr_df[which(mutation_chr_df$phase == "B"), ], aes(pos, mut_cpn_lim, fill = color_index), color = "black") +
          geom_segment(data = unique(mutation_chr_df[, c("seg_start", "seg_end", "allelic_imbalance", "phase", "seg_y_pos_maj_lim", "seg_y_pos_min_lim", "color_index")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_min_lim, yend = seg_y_pos_min_lim, color = 'allele B', alpha = allelic_imbalance), linewidth = 2) +
          scale_color_manual(name = "copy number phase", values=c("allele B" = "#313695")) +
          scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
          geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
          scale_fill_manual(name = "", breaks = c("one", "two", "not_possible"), labels = c("one", "two", "not possible"), values = c("#d8b365", "#5ab4ac", "#bababa"), guide = "none") +
          facet_wrap(~haplotype, ncol = 1, strip.position = "right") +
          ylab("mutation copy number") +
          xlab(chr) +
          ylim(mut_cpn_ylim_chr_min, mut_cpn_ylim_chr) +
          theme_bw() +
          theme(strip.background = element_rect(fill = "white"))
        
        pdf(paste0(output_path, tumour, "/", region, "/plots/phased_mutcpn/per_chromosome/", region, "_", chr, "_genome_phased_mutation_copynumber_per_haplotype_mutcnp.pdf"), width = 10, height = 7)
        print(p_chr_genome2)
        dev.off()
      }
      
    }
  }
}

# function to compare mutations on different copy number states of different alleles
compare.cn.states.fun <- function(phased_data, region){
  
  # for each haplotype seperately gte the number of mutations for each copy number state
  # also keep the actual copy number so we can know which mutations are showing us an "ancestral" copy number state
  
  mutcn_df <- unique(phased_data[, c("chr", "pos", "maj_cn", "min_cn", "seg_start", "seg_end", "mut_cpn", "mut_cpn_round", "phase")])
  mutcn_df$segID <- paste0(mutcn_df$chr, "_", mutcn_df$seg_start, "_", mutcn_df$seg_end)
  
  mutcn_count <- mutcn_df %>% 
    group_by(segID, mut_cpn_round, phase) %>%
    summarize(n = n())
  
  mutcn_count <- left_join(mutcn_count, unique(mutcn_df[, c("segID", "maj_cn", "min_cn")]))
  
  # do a prop test between each haplotype for each CN state
  mutcn_count <- dcast(mutcn_count, segID + mut_cpn_round + maj_cn + min_cn ~ phase, value.var = "n")
  mutcn_count[is.na(mutcn_count$A), "A"] <- 0
  mutcn_count[is.na(mutcn_count$B), "B"] <- 0
  mutcn_count$total_mutations <- mutcn_count$A + mutcn_count$B
  
  mutcn_count <- mutcn_count %>%
    rowwise() %>%
    mutate(prop_test_p = prop.test(x = c(A, B), n = c(total_mutations, total_mutations))$p.value)
  
  # calculate ratio between A and B
  mutcn_count$haplo_ratio <- mutcn_count$A / mutcn_count$B
  mutcn_count$region <- region
  
  return(mutcn_count)
  
}

# calc.timing.fun <- function(mutation_df = phased_data, output_path){
#   
#   #for each segment and each haplotype separately
#   mutation_dfA <- mutation_df[which(mutation_df$phase == "A"), ]
#   mutation_dfB <- mutation_df[which(mutation_df$phase == "B"), ]
#   
#   haplotypeA_gaintiming <- list()
#   for(segID in unique(mutation_dfA$segID)){
#     
#     df <- mutation_dfA[mutation_dfA$segID == segID, ]
#     
#     #is timing possible?
#     gain_present <- ifelse(any(df$maj_cn_gain), TRUE, FALSE)
#     two_optima <-  ifelse(length(unique(df$local_optimum)) > 1, TRUE, FALSE)
# 
#     timing_possible <- ifelse(all(gain_present, two_optima), TRUE, FALSE)
# 
# 
#     if(timing_possible){
#       
#       # how many copy number steps are there
#       number_steps <- length(unique(df$local_optimum))
#       
#       #for each copy number state count the numbers of mutations
#       df$mutID <- paste0(df$chr, ":", df$pos)
#       df <- df %>% group_by(local_optimum) %>% mutate(num_muts_on_cn_step = n())
#       df <- data.frame(df)
# 
#       #only keep relevant info
#       mutation_count_df <- unique(df[, c("sample", "chr", "seg_start", "seg_end", "phase", "maj_cn_gain", "min_cn_gain", "local_optimum", "num_muts_on_cn_step")])
#       mutation_count_df <- mutation_count_df[order(mutation_count_df$local_optimum, decreasing = F), ]
# 
#       mutation_count_df$copy_number_before_gain <- mutation_count_df$local_optimum
# 
#       for(s in 1:nrow(mutation_count_df)){
#         mutation_count_df[s, "copy_number_after_gain"] <- ifelse(s < nrow(mutation_count_df), mutation_count_df[s+1, "copy_number_before_gain"], mutation_count_df[s, "copy_number_before_gain"])
#       }
#       
#       mutation_count_df$adjustment_factor <- mutation_count_df$copy_number_before_gain
# 
#       #if adjustment factor isbelow one make it into one
#       mutation_count_df[mutation_count_df$adjustment_factor < 1, "adjustment_factor"] <- 1
# 
#       mutation_count_df$number_muts_cn_adjusted <- mutation_count_df$num_muts_on_cn_step / mutation_count_df$adjustment_factor
# 
#       # how many copy number steps are there
#       number_steps <- length(unique(mutation_count_df$local_optimum))
#       steps <- sort(unique(mutation_count_df$local_optimum))
#       
#       gain_time_df <- data.frame(tumour = NA,
#                                 region = NA,
#                                 segID = NA,
#                                 gain_index = NA,
#                                 copy_number_before_gain = NA,
#                                 copy_number_after_gain = NA,
#                                 number_early_muts_cn_adjusted = NA,
#                                 number_all_muts_cn_adjusted = NA,
#                                 gain_timing = NA,
#                                 loCI = NA,
#                                 hiCI = NA) 
# 
#       for(event in 1:(number_steps - 1)){
#         
#         gain_index                            <- event
#         copy_number_before_gain               <- mutation_count_df[event, "copy_number_before_gain"]
#         copy_number_after_gain                <- mutation_count_df[event, "copy_number_after_gain"]
#         
#         number_early_muts                     <- sum(mutation_count_df[1:event, "number_muts_cn_adjusted"])
#         number_all_muts                       <- sum(mutation_count_df[, "number_muts_cn_adjusted"])
#         
#         gain_timing                           <- number_early_muts / number_all_muts
#         error                                 <- qnorm(0.975) * sqrt(gain_timing * (1 - gain_timing) / number_all_muts)
#         loCI                                  <- gain_timing - error
#         hiCI                                  <- gain_timing + error
#         
#         gain_time_df[event, "tumour"]                                 <- tumour
#         gain_time_df[event, "region"]                                 <- region
#         gain_time_df[event, "segID"]                                  <- segID
#         gain_time_df[event, "gain_index"]                             <- gain_index
#         gain_time_df[event, "copy_number_before_gain"]                <- copy_number_before_gain
#         gain_time_df[event, "copy_number_after_gain"]                 <- copy_number_after_gain
#         gain_time_df[event, "number_early_muts_cn_adjusted"]          <- number_early_muts
#         gain_time_df[event, "number_all_muts_cn_adjusted"]            <- number_all_muts
#         gain_time_df[event, "gain_timing"]                            <- gain_timing
#         gain_time_df[event, "loCI"]                                   <- loCI
#         gain_time_df[event, "hiCI"]                                   <- hiCI
# 
#       }  
#         
#         haplotypeA_gaintiming[[segID]] <- gain_time_df
#         
#     } else {
#       
#       gain_time_df <- data.frame(tumour = tumour,
#                                 region = region,
#                                 segID = segID,
#                                 gain_index = NA,
#                                 copy_number_before_gain = NA,
#                                 copy_number_after_gain = NA,
#                                 number_early_muts_cn_adjusted = NA,
#                                 number_all_muts_cn_adjusted = NA,
#                                 gain_timing = NA,
#                                 loCI = NA,
#                                 hiCI = NA) 
# 
#       
#       haplotypeA_gaintiming[[segID]] <- gain_time_df
#     
#     }
# 
#   }
#   haplotypeA_gaintiming <- do.call(rbind, haplotypeA_gaintiming)
#   haplotypeA_gaintiming$phase <- "A"
#   
#   haplotypeB_gaintiming <- list()
#   for(segID in unique(mutation_dfB$segID)){
#     
#     df <- mutation_dfB[mutation_dfB$segID == segID, ]
#     
#     #is timing possible?
#     gain_present <- ifelse(any(df$min_cn_gain), TRUE, FALSE)
#     two_optima <-  ifelse(length(unique(df$local_optimum)) > 1, TRUE, FALSE)
# 
#     timing_possible <- ifelse(all(gain_present, two_optima), TRUE, FALSE)
# 
#     if(timing_possible){
#       
#        # how many copy number steps are there
#       number_steps <- length(unique(df$local_optimum))
#       
#       #for each copy number state count the numbers of mutations
#       df$mutID <- paste0(df$chr, ":", df$pos)
#       df <- df %>% group_by(local_optimum) %>% mutate(num_muts_on_cn_step = n())
#       df <- data.frame(df)
# 
#       #only keep relevant info
#       mutation_count_df <- unique(df[, c("sample", "chr", "seg_start", "seg_end", "phase", "maj_cn_gain", "min_cn_gain", "local_optimum", "num_muts_on_cn_step")])
#       mutation_count_df <- mutation_count_df[order(mutation_count_df$local_optimum, decreasing = F), ]
# 
#       mutation_count_df$copy_number_before_gain <- mutation_count_df$local_optimum
# 
#       for(s in 1:nrow(info)){
#         mutation_count_df[s, "copy_number_after_gain"] <- ifelse(s < nrow(mutation_count_df), mutation_count_df[s+1, "copy_number_before_gain"], mutation_count_df[s, "copy_number_before_gain"])
#       }
#       
#       mutation_count_df$adjustment_factor <- mutation_count_df$copy_number_before_gain
# 
#       #if adjustment factor isbelow one make it into one
#       mutation_count_df[mutation_count_df$adjustment_factor < 1, "adjustment_factor"] <- 1
# 
#       mutation_count_df$number_muts_cn_adjusted <- mutation_count_df$num_muts_on_cn_step / mutation_count_df$adjustment_factor
# 
#       # how many copy number steps are there
#       number_steps <- length(unique(mutation_count_df$local_optimum))
#       steps <- sort(unique(mutation_count_df$local_optimum))
#       
#       gain_time_df <- data.frame(tumour = NA,
#                                 region = NA,
#                                 segID = NA,
#                                 gain_index = NA,
#                                 copy_number_before_gain = NA,
#                                 copy_number_after_gain = NA,
#                                 number_early_muts_cn_adjusted = NA,
#                                 number_all_muts_cn_adjusted = NA,
#                                 gain_timing = NA,
#                                 loCI = NA,
#                                 hiCI = NA) 
# 
#       for(event in 1:(number_steps - 1)){
#         
#         gain_index                            <- event
#         copy_number_before_gain               <- mutation_count_df[event, "copy_number_before_gain"]
#         copy_number_after_gain                <- mutation_count_df[event, "copy_number_after_gain"]
#         
#         number_early_muts                     <- sum(mutation_count_df[1:event, "number_muts_cn_adjusted"])
#         number_all_muts                       <- sum(mutation_count_df[, "number_muts_cn_adjusted"])
#         
#         gain_timing                           <- number_early_muts / number_all_muts
#         error                                 <- qnorm(0.975) * sqrt(gain_timing * (1 - gain_timing) / number_all_muts)
#         loCI                                  <- gain_timing - error
#         hiCI                                  <- gain_timing + error
#         
#         gain_time_df[event, "tumour"]                                 <- tumour
#         gain_time_df[event, "region"]                                 <- region
#         gain_time_df[event, "segID"]                                  <- segID
#         gain_time_df[event, "gain_index"]                             <- gain_index
#         gain_time_df[event, "copy_number_before_gain"]                <- copy_number_before_gain
#         gain_time_df[event, "copy_number_after_gain"]                 <- copy_number_after_gain
#         gain_time_df[event, "number_early_muts_cn_adjusted"]          <- number_early_muts
#         gain_time_df[event, "number_all_muts_cn_adjusted"]            <- number_all_muts
#         gain_time_df[event, "gain_timing"]                            <- gain_timing
#         gain_time_df[event, "loCI"]                                   <- loCI
#         gain_time_df[event, "hiCI"]                                   <- hiCI
# 
#       }  
# 
#       haplotypeB_gaintiming[[segID]] <- gain_time_df
# 
#     } else {
#       
#       gain_time_df <- data.frame(tumour = tumour,
#                                 region = region,
#                                 segID = segID,
#                                 gain_index = NA,
#                                 copy_number_before_gain = NA,
#                                 copy_number_after_gain = NA,
#                                 number_early_muts_cn_adjusted = NA,
#                                 number_all_muts_cn_adjusted = NA,
#                                 gain_timing = NA,
#                                 loCI = NA,
#                                 hiCI = NA) 
# 
#       haplotypeB_gaintiming[[segID]] <- gain_time_df
#       
#     }
#     
#   }
#   haplotypeB_gaintiming <- do.call(rbind, haplotypeB_gaintiming)
#   haplotypeB_gaintiming$phase <- "B"
#   
#   # deal with the segments we can't time
#   mutation_dfNA <- mutation_df[is.na(mutation_df$phase), ]
#   
#   
#   haplotypeNA_gaintiming <- data.frame(tumour = NA,
#                                       region = NA,
#                                       segID = NA,
#                                       gain_index = NA,
#                                       copy_number_before_gain = NA,
#                                       copy_number_after_gain = NA,
#                                       number_early_muts_cn_adjusted = NA,
#                                       number_all_muts_cn_adjusted = NA,
#                                       gain_timing = NA,
#                                       loCI = NA,
#                                       hiCI = NA,
#                                       phase = NA) 
# 
#   for(segID in 1:length(unique(mutation_dfNA$segID))){
#     
#     haplotypeNA_gaintiming[segID, "tumour"]                                 <- tumour
#     haplotypeNA_gaintiming[segID, "region"]                                 <- region
#     haplotypeNA_gaintiming[segID, "segID"]                                  <- unique(mutation_dfNA$segID)[segID]
#     haplotypeNA_gaintiming[segID, "gain_index"]                             <- NA
#     haplotypeNA_gaintiming[segID, "copy_number_before_gain"]                <- NA
#     haplotypeNA_gaintiming[segID, "copy_number_after_gain"]                 <- NA
#     haplotypeNA_gaintiming[segID, "number_early_muts_cn_adjusted"]          <- NA
#     haplotypeNA_gaintiming[segID, "number_all_muts_cn_adjusted"]            <- NA
#     haplotypeNA_gaintiming[segID, "gain_timing"]                            <- NA
#     haplotypeNA_gaintiming[segID, "loCI"]                                   <- NA
#     haplotypeNA_gaintiming[segID, "hiCI"]                                   <- NA
#     haplotypeNA_gaintiming[segID, "phase"]                                  <- "not_possible"
#     
#   }
# 
#   #combine the gain timing data frames
#   gaintiming <- rbind(haplotypeA_gaintiming, haplotypeB_gaintiming, haplotypeNA_gaintiming)
#   
#   if(!dir.exists(paste0(output_path, tumour, "/", region, "/data/"))){
#     dir.create(paste0(output_path, tumour, "/", region, "/data/"), recursive = T)
#   }
#   
#   write.table(gaintiming, paste0(output_path, tumour, "/", region, "/data/", region, "_copynumber_gain_timing.txt"),
#               col.names = T,
#               row.names = F,
#               quote = F,
#               sep = "\t")
#   
#   return(gaintiming)
# 
# }

# function to plot each copy number segment on the x axis and timing estimate on the y axis
plot.timing.fun <- function(timing_df = timed_data, mutation_df = phased_data, cn_df = cn, all_mut = data, output_path, genome = BSgenome.Hsapiens.UCSC.hg38, SV_data = SVs){
  

  # separate the segIDs
  timing_df <- timing_df %>% separate(segID, into = c("chr", "seg_start", "seg_end"), sep = "_", remove = FALSE)
  timing_df$seg_start <- as.numeric(as.character(timing_df$seg_start))
  timing_df$seg_end <- as.numeric(as.character(timing_df$seg_end))

  # first, make the genome continous
  timing_df <- genome.continous.fun(timing_df, genome)

  chr_sizes <- seqlengths(genome)[paste0("chr", 1:22)]
  chr_lines <- cumsum(as.numeric(as.character(chr_sizes)))
  chr_midpoint <- chr_lines - chr_sizes / 2
  
  xaxis_label <- data.frame(chr_midpoint)
  xaxis_label$label <- rownames(xaxis_label)
  
  timing_df <- unique(timing_df[, c("segID", "chr", "cluster", "gain_index", "adjust_seg_start", "adjust_seg_end", "haplotype", "gain_timing", "loCI", "hiCI", "seg_start", "seg_end")])
  colnames(timing_df) <- c("segID", "chr", "cluster", "gain_index", "adjust_seg_start", "adjust_seg_end", "phase", "gain_timing", "loCI", "hiCI", "seg_start", "seg_end")
  
  # where confidence intervals are greater than 1 or smaller than 0 make them 1 or 0 because anything higher or lower makes no sense
  timing_df[which(timing_df$hiCI > 1), "hiCI"] <- 1
  timing_df[which(timing_df$loCI < 0), "loCI"] <- 0
  
  # make dummy collumn just for plotting purposes cause I can't figure out how to give the plot a title on the rigth without faceting
  timing_df$haplotype <- "haplotype"
  
  #where timing is not possible make timing into 0 so that we can still plot it
  timing_df[is.na(timing_df$gain_timing), "phase"] <- "not_possible"
  timing_df[is.na(timing_df$gain_timing), "gain_timing"] <- 0
  timing_df[is.na(timing_df$loCI), "loCI"] <- 0
  timing_df[is.na(timing_df$hiCI), "hiCI"] <- 0
    
  # for an initial plot only keep the first clonal gain
  timing_df[timing_df$phase == "not_possible", "gain_index"] <- 0

  
  p_genome1 <- ggplot() +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    geom_rect(data = timing_df[timing_df$phase %in% c("A", "not_possible") & timing_df$gain_index %in% c(0, 1), ], aes(xmin = adjust_seg_start, xmax = adjust_seg_end, ymin = loCI, ymax = hiCI), fill = "darkred", color = NA, alpha = 0.15) +
    geom_segment(data = timing_df[timing_df$phase %in% c("A", "not_possible") & timing_df$gain_index %in% c(0, 1), ],
                 aes(x = adjust_seg_start, xend = adjust_seg_end, y = gain_timing, yend = gain_timing, color = phase, alpha = phase),
                 linewidth = 2) +
    scale_color_manual(name = "", breaks = c("A", "not_possible"),
                       values = c("A" = "darkred",
                                  "not_possible" = "white"),
                       labels = c("gain timing possible", "gain timing not possible"),
                       guide = "none") +
    scale_alpha_manual(breaks = c("A", "not_possible"), values = c(1, 0), guide = "none") +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.00\n(early)", "0.25", "0.50", "0.75", "1.00\n(late)"), limits = c(0, 1)) + 
    facet_wrap(~haplotype, labeller = labeller(haplotype = c(haplotype = "haplotype A")), strip.position = "right") +
    ylab("copy number state timing") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_rect(fill = "white"),
          plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"))

  p_genome2 <- ggplot() +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    geom_rect(data = timing_df[timing_df$phase %in% c("B", "not_possible") & timing_df$gain_index %in% c(0, 1), ], aes(xmin = adjust_seg_start, xmax = adjust_seg_end, ymin = loCI, ymax = hiCI), fill = "#313695", color = NA, alpha = 0.15) +
    geom_segment(data = timing_df[timing_df$phase %in% c("B", "not_possible") & timing_df$gain_index %in% c(0, 1), ],
                 aes(x = adjust_seg_start, xend = adjust_seg_end, y = gain_timing, yend = gain_timing, color = phase, alpha = phase),
                 linewidth = 2) +
    scale_color_manual(name = "", breaks = c("B", "not_possible"),
                       values = c("B" = "#313695",
                                  "not_possible" = "white"),
                       labels = c("gain timing possible", "gain timing not possible"),
                       guide = "none") +
    scale_alpha_manual(breaks = c("A", "not_possible"), values = c(1, 0), guide = "none") +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.00\n(early)", "0.25", "0.50", "0.75", "1.00\n(late)"), limits = c(0, 1)) + 
    facet_wrap(~haplotype, labeller = labeller(haplotype = c(haplotype = "haplotype B")), strip.position = "right") +
    ylab("copy number state timing") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_rect(fill = "white"),
          plot.margin = unit(c(0, 0.1, 0.1, 0.1), "cm"))

  # make a plot with copy number
  # also plot the ancestral copy number on here
  
  #make genome continuous
  cn_df$chr <- paste0("chr", cn_df$chr)
  colnames(cn_df)[2:3] <- c("seg_start", "seg_end")
  cn_df <- genome.continous.fun(cn_df, genome)

  #make max cn into 10
  cn_df[which(cn_df$maj_cn > 10), "maj_cn_lim"] <- 10
  cn_df[which(cn_df$maj_cn <= 10), "maj_cn_lim"] <- cn_df[which(cn_df$maj_cn <= 10), "maj_cn"]
  cn_df[which(cn_df$min_cn > 10), "min_cn_lim"] <- 10
  cn_df[which(cn_df$min_cn <= 10), "min_cn_lim"] <- cn_df[which(cn_df$min_cn <= 10), "min_cn"]

  # for the cross genome copy number profile plot we want to indicate where SVs are located
  # need to reformat the SV table for thos
  SVs_long <- SV_data
  SVs_1 <- SVs_long[, c("chr1", "pos1")]
  colnames(SVs_1) <- c("chr", "pos")
  SVs_2 <- SVs_long[, c("chr2", "pos2")]
  colnames(SVs_2) <- c("chr", "pos")
  SVs_long <- rbind(SVs_1, SVs_2)
  
  # make continous
  SVs_long <- genome.continous.fun(SVs_long, genome)
    

  p_cn <- ggplot() +
            geom_vline(xintercept = chr_lines, linetype = 2, color = "black") +
            geom_vline(xintercept = SVs_long$adjust_pos, linewidth = 0.05, color = "#bababa") +
            geom_segment(data = cn_df,
                        aes(x = adjust_seg_start, xend = adjust_seg_end, y = maj_cn_lim, yend = maj_cn_lim, color = 'allele A'), linewidth = 2, alpha = 0.7) +
            geom_segment(data = cn_df,
                        aes(x = adjust_seg_start, xend = adjust_seg_end, y = min_cn_lim, yend = min_cn_lim, color = 'allele B'), linewidth = 2, alpha = 0.7) +
            # geom_segment(data = cn_df[cn_df$timing == "early" & cn_df$phase == "A", ],
            #             aes(x = adjust_seg_start, xend = adjust_seg_end, y = local_optimum_lim, yend = local_optimum_lim, color = 'allele A'), linewidth = 1, alpha = 0.2) +
            # geom_segment(data = cn_df[cn_df$timing == "early" & cn_df$phase == "B", ],
            #             aes(x = adjust_seg_start, xend = adjust_seg_end, y = local_optimum_lim, yend = local_optimum_lim, color = 'allele B'), linewidth = 1, alpha = 0.2) +
            scale_color_manual(name = "haplotype", values = c("allele A" = "darkred",
                                                            "allele B" = "#313695"),
                                                    labels = c("haplotype A", "haplotype B")) +
            scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
            scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
            ylab("copy number") +
            xlab("") +
            theme_bw() +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  strip.background = element_rect(fill = "white"),
                  plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm"),
                  panel.grid.minor.y = element_blank())   

    # also plot the mutations
    mut_df <- unique(mutation_df[, c("chr", "seg_start", "seg_end", "pos", "phase", "mut_cpn", "maj_cn_gain", "min_cn_gain", "maj_cn", "min_cn", "allelic_imbalance", "CLUSTER")])

    #only keep clonal cluster
    mut_df <- mut_df[mut_df$CLUSTER %in% unique(timing_df$cluster), ]
    mut_df <- mut_df[!is.na(mut_df$CLUSTER), ]
    
    #make genome continuous
    mut_df <- genome.continous.fun(mut_df, genome)
    
    mut_df$local_optimum <- round(mut_df$mut_cpn)
    mut_df <- mut_df[mut_df$local_optimum >= 1, ]
    
    mut_df[which(mut_df$phase == "A" & mut_df$local_optimum > mut_df$maj_cn), "local_optimum"] <- mut_df[which(mut_df$phase == "A" & mut_df$local_optimum > mut_df$maj_cn), "maj_cn"]
    mut_df[which(mut_df$phase == "B" & mut_df$local_optimum > mut_df$min_cn), "local_optimum"] <- mut_df[which(mut_df$phase == "B" & mut_df$local_optimum > mut_df$min_cn), "min_cn"]
    
    #we need to adjust for some noise where there are less than 5 mutations on one state, these mutations get removed
    mut_df$segID <- paste0(mut_df$chr, ":", mut_df$seg_start, ":", mut_df$seg_end, ":", mut_df$phase)
    mut_df <- mut_df %>% group_by(segID, local_optimum) %>% mutate(num_muts_per_optimum_seg = n())
    mut_df <- data.frame(mut_df)
    mut_df <- mut_df[mut_df$num_muts_per_optimum_seg >= 5, ]
    mut_df$num_muts_per_optimum_seg <- NULL
    
    #for  each segment color the top most steps
    mut_df$segID <- paste0(mut_df$chr, ":", mut_df$seg_start, ":", mut_df$seg_end)
    
    mut_df_seg <- list()
    for(seg in unique(mut_df$segID)) {
    
      df <- mut_df[mut_df$segID == seg, ]
      
      dfD <- colnames(df)
      
      dfA <- df[df$phase == "A", ]
      dfB <- df[df$phase == "B", ]
      
      A_optima <- unique(dfA$local_optimum)
      A_optima <- sort(A_optima, decreasing = T)
      A_top <- A_optima[1]
      A_second <- A_optima[2]

      if(all(nrow(dfA) > 0 & length(A_optima) > 1)){
        dfA[dfA$local_optimum == A_top, "color"] <- "top"
        dfA[dfA$local_optimum == A_second, "color"] <- "second"
        dfA[is.na(dfA$color), "color"] <- "irrelevant"
      } else {
        dfA <- data.frame(matrix(NA, ncol = 18, nrow = 1))
        colnames(dfA) <- dfD
        dfA$color <- "irrelevant"
      }

      
      B_optima <- unique(dfB$local_optimum)
      B_optima <- sort(B_optima, decreasing = T)
      B_top <- B_optima[1]
      B_second <- B_optima[2]
      
      if(all(nrow(dfB) > 0 & length(B_optima) > 1)){
        dfB[which(dfB$local_optimum == B_top), "color"] <- "top"
        dfB[dfB$local_optimum == B_second, "color"] <- "second"
        dfB[is.na(dfB$color), "color"] <- "irrelevant"
      } else {
        dfB <- data.frame(matrix(NA, ncol = 18, nrow = 1))
        colnames(dfB) <- dfD
        dfB$color <- "irrelevant"
      }
      
      df_out <- rbind(dfA, dfB)
      df_out <- df_out[!is.na(df_out$segID), ]
      
      mut_df_seg[[seg]] <- df_out 
    }
    mut_df_seg <- do.call(rbind, mut_df_seg)
    
    #color every other step differently
  #   mut_df_seg <- lapply(unique(mut_df$segID), function(segID){
  # 
  #     df <- mut_df[mut_df$segID == segID, ]
  #     dfA <- df[which(df$phase == "A" & df$timing != "not_possible"), ]
  #     dfB <- df[which(df$phase == "B" & df$timing != "not_possible"), ]
  #     dfNA <- df[df$timing == "not_possible", ]
  # 
  #     # how many steps are there in this segment
  # 
  #     if(nrow(dfA) > 0){
  #       dfA$number_steps <- length(unique(dfA$local_optimum))
  #       color_indexA <- data.frame(local_optimum = sort(unique(dfA$local_optimum)),
  #                               color_index = rep(c("one", "two"), length.out = length(unique(dfA$local_optimum))))
  #       dfA <- suppressMessages(left_join(dfA, color_indexA))
  #     }
  #     
  #     if(nrow(dfB) > 0){
  #       dfB$number_steps <- length(unique(dfB$local_optimum))
  #       color_indexB <- data.frame(local_optimum = sort(unique(dfB$local_optimum)),
  #                                 color_index = rep(c("one", "two"), length.out = length(unique(dfB$local_optimum))))
  #       dfB <- suppressMessages(left_join(dfB, color_indexB))
  #     }
  # 
  #     if(nrow(dfNA) > 0){
  #       dfNA$color_index <- "not_possible"
  #       dfNA$number_steps <- NA
  #     }
  #     
  #     df_out <- rbind(dfA, dfB, dfNA)
  # 
  # })
  # mut_df_seg <- do.call(rbind, mut_df_seg)

  mut_df_seg[mut_df_seg$maj_cn_gain == FALSE & mut_df_seg$min_cn_gain == FALSE, "color"] <- "irrelevant"
  mut_df_seg[which(mut_df_seg$maj_cn_gain == FALSE & mut_df_seg$phase == "A"), "color"] <- "irrelevant"
  mut_df_seg[which(mut_df_seg$min_cn_gain == FALSE & mut_df_seg$phase == "B"), "color"] <- "irrelevant"
  mut_df_seg[is.na(mut_df_seg$phase), "color"] <- "irrelevant"
  
  mut_df_seg[mut_df_seg$phase == "A" & mut_df_seg$maj_cn != round(mut_df_seg$maj_cn), "color"] <- "irrelevant"
  mut_df_seg[mut_df_seg$phase == "B" & mut_df_seg$min_cn != round(mut_df_seg$min_cn), "color"] <- "irrelevant"
  
  
  # mut_df_seg[mut_df_seg$timing == "not_possible", "color_index"] <- "not_possible"
  # mut_df_seg$color_index <- factor(mut_df_seg$color_index, levels = c("one", "two", "not_possible"))

  #make max cn into 10
  mut_df_seg[which(mut_df_seg$mut_cpn > 10), "mut_cpn_lim"] <- 10
  mut_df_seg[which(mut_df_seg$mut_cpn <= 10), "mut_cpn_lim"] <- mut_df_seg[which(mut_df_seg$mut_cpn <= 10), "mut_cpn"]
  mut_df_seg[which(mut_df_seg$maj_cn > 10), "maj_cn_lim"] <- 10
  mut_df_seg[which(mut_df_seg$maj_cn <= 10), "maj_cn_lim"] <- mut_df_seg[which(mut_df_seg$maj_cn <= 10), "maj_cn"]
  mut_df_seg[which(mut_df_seg$min_cn > 10), "min_cn_lim"] <- 10
  mut_df_seg[which(mut_df_seg$min_cn <= 10), "min_cn_lim"] <- mut_df_seg[which(mut_df_seg$min_cn <= 10), "min_cn"]

  mut_cpn_ylim_genome_min <- 0
  mut_cpn_ylim_genome_max <- max(mut_df_seg$mut_cpn_lim, mut_df_seg$maj_cn_lim, mut_df_seg$min_cn_lim)

  mut_df_seg$haplotype <- paste0("haplotype", mut_df_seg$phase)
  
  mut_df_seg$color <- factor(mut_df_seg$color, levels = c("top", "second", "irrelevant"))

  if(nrow(mut_df_seg[which(mut_df_seg$phase == "A"), ]) > 0){
    p_mutA <- ggplot() +
      geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
      geom_point(data = mut_df_seg[which(mut_df_seg$phase == "A"), ], aes(adjust_pos, mut_cpn_lim, fill = color), color = "black", alpha = 0.7, shape = 21, stroke = NA, size = 1) +
      scale_fill_manual(name = "", breaks = c("top", "second", "irrelevant"), labels = c("one", "two", "not possible"), values = c("#d8b365", "#5ab4ac", "#bababa"), guide = "none") +
      # geom_segment(data = unique(mut_df_seg[, c("adjust_seg_start", "adjust_seg_end", "allelic_imbalance", "phase", "maj_cn_lim", "min_cn_lim", "color")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = maj_cn_lim, yend = maj_cn_lim, color = 'allele A', alpha = allelic_imbalance), linewidth = 1) +
      geom_segment(data = cn_df,
                   aes(x = adjust_seg_start, xend = adjust_seg_end, y = maj_cn_lim, yend = maj_cn_lim, color = 'allele A'), linewidth = 1, alpha = 0.7) +
      scale_color_manual(name = "copy number phase", values = c("allele A" = "darkred"), guide = "none") +
      scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
      scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
      facet_wrap(~haplotype, labeller = labeller(haplotype = c("haplotypeA" = "haplotype A")), strip.position = "right") +
      ylab("mutation copy number") +
      xlab("") +
      ylim(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome_max) +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            strip.background = element_rect(fill = "white"))
  } else {
    cn_dfA <- cn_df
    cn_dfA$haplotype <- "haplotypeA"
    
    p_mutA <- ggplot() +
      geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
      # geom_segment(data = unique(mut_df_seg[, c("adjust_seg_start", "adjust_seg_end", "allelic_imbalance", "phase", "maj_cn_lim", "min_cn_lim", "color")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = maj_cn_lim, yend = maj_cn_lim, color = 'allele A', alpha = allelic_imbalance), linewidth = 1) +
      geom_segment(data = cn_dfA,
                   aes(x = adjust_seg_start, xend = adjust_seg_end, y = maj_cn_lim, yend = maj_cn_lim, color = 'allele A'), linewidth = 1, alpha = 0.7) +
      scale_color_manual(name = "copy number phase", values = c("allele A" = "darkred"), guide = "none") +
      scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
      scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
      facet_wrap(~haplotype, labeller = labeller(haplotype = c("haplotypeA" = "haplotype A")), strip.position = "right") +
      ylab("mutation copy number") +
      xlab("") +
      ylim(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome_max) +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            strip.background = element_rect(fill = "white"))
    
  }
  

  if(nrow(mut_df_seg[which(mut_df_seg$phase == "B"), ]) > 0){
    p_mutB <- ggplot() +
      geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
      geom_point(data = mut_df_seg[which(mut_df_seg$phase == "B"), ], aes(adjust_pos, mut_cpn_lim, fill = color), color = "black", alpha = 0.7, shape = 21, stroke = NA, size = 1) +
      scale_fill_manual(name = "", breaks = c("top", "second", "irrelevant"), labels = c("one", "two", "not possible"), values = c("#d8b365", "#5ab4ac", "#bababa"), guide = "none") +
      #geom_segment(data = unique(mut_df_seg[, c("adjust_seg_start", "adjust_seg_end", "allelic_imbalance", "phase", "maj_cn_lim", "min_cn_lim", "color")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = min_cn_lim, yend = min_cn_lim, color = 'allele B', alpha = allelic_imbalance), linewidth = 1) +
      geom_segment(data = cn_df,
                   aes(x = adjust_seg_start, xend = adjust_seg_end, y = min_cn_lim, yend = min_cn_lim, color = 'allele B'), linewidth = 1, alpha = 0.7) +
      scale_color_manual(name = "copy number phase", values = c("allele B" = "#313695"), guide = "none") +
      scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
      
      scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
      facet_wrap(~haplotype, labeller = labeller(haplotype = c("haplotypeB" = "haplotype B")), strip.position = "right") +
      ylab("mutation copy number") +
      xlab("") +
      ylim(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome_max) +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            strip.background = element_rect(fill = "white"))
  }   else {
    
    cn_dfB <- cn_df
    cn_dfB$haplotype <- "haplotypeB"
    p_mutB <- ggplot() +
      geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
      #geom_segment(data = unique(mut_df_seg[, c("adjust_seg_start", "adjust_seg_end", "allelic_imbalance", "phase", "maj_cn_lim", "min_cn_lim", "color")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = min_cn_lim, yend = min_cn_lim, color = 'allele B', alpha = allelic_imbalance), linewidth = 1) +
      geom_segment(data = cn_dfB,
                   aes(x = adjust_seg_start, xend = adjust_seg_end, y = min_cn_lim, yend = min_cn_lim, color = 'allele B'), linewidth = 1, alpha = 0.7) +
      scale_color_manual(name = "copy number phase", values = c("allele B" = "#313695"), guide = "none") +
      scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
      
      scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
      facet_wrap(~haplotype, labeller = labeller(haplotype = c("haplotypeB" = "haplotype B")), strip.position = "right") +
      ylab("mutation copy number") +
      xlab("") +
      ylim(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome_max) +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            strip.background = element_rect(fill = "white"))
  }



    # combine the plots
    p_genome_c <- p_cn / p_mutA / p_genome1 / p_mutB / p_genome2 + plot_annotation(title = region, theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
    
    if(!dir.exists(paste0(output_path, tumour, "/", region, "/plots/phased_timing/"))){
      dir.create(paste0(output_path, tumour, "/", region, "/plots/phased_timing/"), recursive = T)
    }
    
    pdf(paste0(output_path, tumour, "/", region, "/plots/phased_timing/", region, "_genome_phased_copynumber_timing.pdf"), width = 18, height = 10)
    print(p_genome_c)
    dev.off()
    
    
    # make a plot for every chromosome
    for(chr in unique(timing_df$chr)){
      
      print(chr)
      # first, timing plots
      timing_df_chr <- timing_df[timing_df$chr == chr, ]

      p_chr1 <- ggplot() +
        geom_rect(data = timing_df_chr[timing_df_chr$phase %in% c("A", "not_possible") & timing_df_chr$gain_index %in% c(0, 1), ], aes(xmin = seg_start, xmax = seg_end, ymin = loCI, ymax = hiCI), fill = "darkred", color = NA, alpha = 0.15) +
        geom_segment(data = timing_df_chr[timing_df_chr$phase %in% c("A", "not_possible") & timing_df_chr$gain_index %in% c(0, 1), ],
                    aes(x = seg_start, xend = seg_end, y = gain_timing, yend = gain_timing, color = phase, alpha = phase),
                    linewidth = 2) +
        scale_color_manual(name = "", breaks = c("A", "not_possible"),
                          values = c("A" = "darkred",
                                      "not_possible" = "white"),
                          labels = c("gain timing possible", "gain timing not possible"),
                          guide = "none") +
        scale_alpha_manual(breaks = c("A", "not_possible"), values = c(1, 0), guide = "none") +
        scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.00\n(early)", "0.25", "0.50", "0.75", "1.00\n(late)"), limits = c(0, 1)) + 
        scale_x_continuous(expand = c(0,0)) +
        facet_wrap(~haplotype, labeller = labeller(haplotype = c(haplotype = "haplotype A")), strip.position = "right") +
        ylab("copy number state timing") +
        xlab("") +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              strip.background = element_rect(fill = "white"),
              plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"))

      p_chr2 <- ggplot() +
        geom_rect(data = timing_df_chr[timing_df_chr$phase %in% c("B", "not_possible") & timing_df_chr$gain_index %in% c(0, 1), ], aes(xmin = seg_start, xmax = seg_end, ymin = loCI, ymax = hiCI), fill = "#313695", color = NA, alpha = 0.15) +
        geom_segment(data = timing_df_chr[timing_df_chr$phase %in% c("B", "not_possible") & timing_df_chr$gain_index %in% c(0, 1), ],
                    aes(x = seg_start, xend = seg_end, y = gain_timing, yend = gain_timing, color = phase, alpha = phase),
                    linewidth = 2) +
        scale_color_manual(name = "", breaks = c("B", "not_possible"),
                          values = c("B" = "#313695",
                                      "not_possible" = "white"),
                          labels = c("gain timing possible", "gain timing not possible"),
                          guide = "none") +
        scale_alpha_manual(breaks = c("A", "not_possible"), values = c(1, 0), guide = "none") +
        scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.00\n(early)", "0.25", "0.50", "0.75", "1.00\n(late)"), limits = c(0, 1)) + 
        scale_x_continuous(expand = c(0,0)) +
        facet_wrap(~haplotype, labeller = labeller(haplotype = c(haplotype = "haplotype B")), strip.position = "right") +
        ylab("copy number state timing") +
        xlab("") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              strip.background = element_rect(fill = "white"),
              plot.margin = unit(c(0, 0.1, 0.1, 0.1), "cm"))

      # and plot the copy number
      cn_df_chr <- cn_df[cn_df$chr == chr, ]
      
      # here we want SVs linekd to each other 
      
      SVs_chr            <- SVs[which(SVs$chr1 == chr & SVs$chr2 == chr), ]
      SVs_chr_long       <- SVs[which(SVs$chr1 == chr | SVs$chr2 == chr), ]
      SVs_chr_long       <- SVs_chr_long[SVs_chr_long$chr1 != SVs_chr_long$chr2, ]
      SVs_chr_long1      <- SVs_chr_long[, c("chr1", "pos1")]
      colnames(SVs_chr_long1) <- c("chr", "pos")
      SVs_chr_long2      <- SVs_chr_long[, c("chr2", "pos2")]
      colnames(SVs_chr_long2) <- c("chr", "pos")
      SVs_chr_long       <- rbind(SVs_chr_long1, SVs_chr_long2)

      # where should the SVS curves go?
      maxcn <- round(max(cn_df_chr$maj_cn, cn_df_chr$min_cn) + 0.55)

        
      p_chr_cn <- ggplot() +
                geom_segment(data = cn_df_chr,
                            aes(x = seg_start, xend = seg_end, y = maj_cn_lim, yend = maj_cn_lim, color = 'allele A'), linewidth = 2, alpha = 0.7) +
                geom_segment(data = cn_df_chr,
                            aes(x = seg_start, xend = seg_end, y = min_cn_lim, yend = min_cn_lim, color = 'allele B'), linewidth = 2, alpha = 0.7) +
                # geom_segment(data = cn_df[cn_df$timing == "early" & cn_df$phase == "A", ],
                #             aes(x = adjust_seg_start, xend = adjust_seg_end, y = local_optimum_lim, yend = local_optimum_lim, color = 'allele A'), linewidth = 1, alpha = 0.2) +
                # geom_segment(data = cn_df[cn_df$timing == "early" & cn_df$phase == "B", ],
                #             aes(x = adjust_seg_start, xend = adjust_seg_end, y = local_optimum_lim, yend = local_optimum_lim, color = 'allele B'), linewidth = 1, alpha = 0.2) +
                scale_color_manual(name = "haplotype", values = c("allele A" = "darkred",
                                                                "allele B" = "#313695"),
                                                        labels = c("haplotype A", "haplotype B")) +
                scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
                scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
                ylab("copy number") +
                xlab("") +
                theme_bw() +
                theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      panel.grid.major.x = element_blank(),
                      panel.grid.minor.x = element_blank(),
                      strip.background = element_rect(fill = "white"),
                      plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm"),
                      panel.grid.minor.y = element_blank())   
          
      # and plot the mutations
      mut_df_seg_chr <- mut_df_seg[mut_df_seg$chr == chr, ]
      
      if(nrow(mut_df_seg_chr) > 0){
        mut_cpn_ylim_chr_min <- 0
        mut_cpn_ylim_chr_max <- max(mut_df_seg_chr$mut_cpn_lim, mut_df_seg_chr$maj_cn_lim, mut_df_seg_chr$min_cn_lim)
        
        if(nrow(mut_df_seg_chr[which(mut_df_seg_chr$phase == "A"), ]) > 0){
          
          cn_df_A <- cn_df_chr
          cn_df_A$phase <- "A"
          
          p_chr_mutA <- ggplot() +
            # geom_vline(xintercept = SVs_chr_long$pos, linewidth = 0.05, alpha = 0.7, colour = "#DDCC77") +
            scale_colour_manual(name = "Structural\nVariant\nType", values = c("#88CCEE", "#CC6677", "#117733","#44AA99", "#DDCC77"), labels = c("deletion", "duplication", "h2h inversion", "t2t inversion", "single breakend")) +
            geom_curve(data = SVs_chr[SVs_chr$SVTYPE == "DEL", ], aes(x = pos1, y = maxcn, xend = pos2, yend = maxcn), colour = "#88CCEE", curvature = -0.2, linewidth = 1) +
            geom_segment(data = SVs_chr[SVs_chr$SVTYPE == "DEL", ], aes(x = pos1, xend = pos1, y = 0, yend = maxcn), linewidth = 0.5, colour = "#88CCEE") +
            geom_segment(data = SVs_chr[SVs_chr$SVTYPE == "DEL", ], aes(x = pos2, xend = pos2, y = 0, yend = maxcn), linewidth = 0.5, colour = "#88CCEE") +
            geom_curve(data = SVs_chr[SVs_chr$SVTYPE == "h2hINV", ], aes(x = pos1, y = maxcn, xend = pos2, yend = maxcn), colour = "#117733", curvature = -0.2, linewidth = 1) +
            geom_segment(data = SVs_chr[SVs_chr$SVTYPE == "h2hINV", ], aes(x = pos1, xend = pos1, y = 0, yend = maxcn), linewidth = 0.5, colour = "#117733") +
            geom_segment(data = SVs_chr[SVs_chr$SVTYPE == "h2hINV", ], aes(x = pos2, xend = pos2, y = 0, yend = maxcn), linewidth = 0.5, colour = "#117733") +
            geom_curve(data = SVs_chr[SVs_chr$SVTYPE == "t2tINV", ], aes(x = pos1, y = maxcn, xend = pos2, yend = maxcn), colour = "#44AA99", curvature = -0.2, linewidth = 1) +
            geom_segment(data = SVs_chr[SVs_chr$SVTYPE == "t2tINV", ], aes(x = pos1, xend = pos1, y = 0, yend = maxcn), linewidth = 0.5, colour = "#44AA99") +
            geom_segment(data = SVs_chr[SVs_chr$SVTYPE == "t2tINV", ], aes(x = pos2, xend = pos2, y = 0, yend = maxcn), linewidth = 0.5, colour = "#44AA99") +
            geom_curve(data = SVs_chr[SVs_chr$SVTYPE == "DUP", ], aes(x = pos1, y = maxcn, xend = pos2, yend = maxcn), colour = "#CC6677", curvature = -0.2, linewidth = 1) +
            geom_segment(data = SVs_chr[SVs_chr$SVTYPE == "DUP", ], aes(x = pos1, xend = pos1, y = 0, yend = maxcn), linewidth = 0.5, colour = "#CC6677") +
            geom_segment(data = SVs_chr[SVs_chr$SVTYPE == "DUP", ], aes(x = pos2, xend = pos2, y = 0, yend = maxcn), linewidth = 0.5, colour = "#CC6677") +
            new_scale_color() +
            geom_point(data = mut_df_seg_chr[which(mut_df_seg_chr$phase == "A"), ], aes(pos, mut_cpn_lim, fill = color), color = "black", alpha = 0.7, shape = 21, stroke = NA, size = 1) +
            scale_fill_manual(name = "", breaks = c("top", "second", "irrelevant"), labels = c("one", "two", "not possible"), values = c("#d8b365", "#5ab4ac", "#bababa"), guide = "none") +
            geom_segment(data = cn_df_A,
                         aes(x = seg_start, xend = seg_end, y = maj_cn_lim, yend = maj_cn_lim, color = 'allele A'), linewidth = 1, alpha = 0.7) +
            scale_color_manual(name = "copy number phase", values = c("allele A" = "darkred"), guide = "none") +
            scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
            scale_x_continuous(expand = c(0,0)) +
            facet_wrap(~phase, labeller = labeller(phase = c("A" = "haplotype A")), strip.position = "right") +
            ylab("mutation copy number") +
            xlab("") +
            ylim(mut_cpn_ylim_chr_min, max(mut_cpn_ylim_chr_max, maxcn)) +
            theme_bw() +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  strip.background = element_rect(fill = "white"))
        } else {
          
          cn_df_A <- cn_df_chr
          cn_df_A$phase <- "A"
          
          p_chr_mutA <- ggplot() +
            geom_segment(data = cn_df_A,
                         aes(x = adjust_seg_start, xend = adjust_seg_end, y = maj_cn_lim, yend = maj_cn_lim, color = 'allele A'), linewidth = 1, alpha = 0.7) +
            scale_color_manual(name = "copy number phase", values = c("allele A" = "darkred"), guide = "none") +
            scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
            scale_x_continuous(expand = c(0,0)) +
            facet_wrap(~phase, labeller = labeller(phase = c("A" = "haplotype A")), strip.position = "right") +
            ylab("mutation copy number") +
            xlab("") +
            ylim(mut_cpn_ylim_chr_min, mut_cpn_ylim_chr_max) +
            theme_bw() +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  strip.background = element_rect(fill = "white"))
        }
        
        if(nrow(mut_df_seg_chr[which(mut_df_seg_chr$phase == "B"), ]) > 0){
          
          cn_df_B <- cn_df_chr
          cn_df_B$phase <- "B"
          
          p_chr_mutB <- ggplot() +
            geom_point(data = mut_df_seg_chr[which(mut_df_seg_chr$phase == "B"), ], aes(adjust_pos, mut_cpn_lim, fill = color), color = "black", alpha = 0.7, shape = 21, stroke = NA, size = 1) +
            scale_fill_manual(name = "", breaks = c("top", "second", "irrelevant"), labels = c("one", "two", "not possible"), values = c("#d8b365", "#5ab4ac", "#bababa"), guide = "none") +
            geom_segment(data = cn_df_B,
                         aes(x = adjust_seg_start, xend = adjust_seg_end, y = min_cn_lim, yend = min_cn_lim, color = 'allele B'), linewidth = 1, alpha = 0.7) +
            scale_color_manual(name = "copy number phase", values = c("allele B" = "#313695"), guide = "none") +
            scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
            scale_x_continuous(expand = c(0,0)) +
            facet_wrap(~phase, labeller = labeller(phase = c("B" = "haplotype B")), strip.position = "right") +
            ylab("mutation copy number") +
            xlab("") +
            ylim(mut_cpn_ylim_chr_min, mut_cpn_ylim_chr_max) +
            theme_bw() +
            theme(legend.position = "none",
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  strip.background = element_rect(fill = "white"))
        } else {
          
          cn_df_B <- cn_df_chr
          cn_df_B$phase <- "B"
          
          p_chr_mutB <- ggplot() +
            geom_segment(data = cn_df_B,
                         aes(x = adjust_seg_start, xend = adjust_seg_end, y = min_cn_lim, yend = min_cn_lim, color = 'allele B'), linewidth = 1, alpha = 0.7) +
            scale_color_manual(name = "copy number phase", values = c("allele B" = "#313695"), guide = "none") +
            scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
            scale_x_continuous(expand = c(0,0)) +
            facet_wrap(~phase, labeller = labeller(phase = c("B" = "haplotype B")), strip.position = "right") +
            ylab("mutation copy number") +
            xlab("") +
            ylim(mut_cpn_ylim_chr_min, mut_cpn_ylim_chr_max) +
            theme_bw() +
            theme(legend.position = "none",
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  strip.background = element_rect(fill = "white"))
          
          
        }
        
        p_chr_c <- p_chr_cn / p_chr_mutA / p_chr1 / p_chr_mutB / p_chr2 + plot_annotation(title = paste0(region, " : ", chr), theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
        
        if(!dir.exists(paste0(output_path, tumour, "/", region, "/plots/phased_timing/per_chromosome/"))){
          dir.create(paste0(output_path, tumour, "/", region, "/plots/phased_timing/per_chromosome/"), recursive = T)
        }
        
        pdf(paste0(output_path, tumour, "/", region, "/plots/phased_timing/per_chromosome/", region, "_", chr, "_phased_copynumber_timing.pdf"), width = 10, height = 10)
        print(p_chr_c)
        dev.off()
        
        
      }

    }
    
    # we also want to plot the various timing of the different steps
    # lets get some colors, gainindex 0 is always white, although that doesnt really matter cause we make it transparent anyway
    
    # colorsA <- get.colors(max_cn = max(as.numeric(as.character(timing_df[timing_df$phase %in% c("A", "not_possible"), "gain_index"]))))
    # colorsA <- c("white", colorsA)
    # 
    # timing_df$gain_index <- factor(timing_df$gain_index, levels = sort(unique(timing_df$gain_index)))
    # 
    # p_genome1_steps <- ggplot() +
    #         geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    #         geom_rect(data = timing_df[timing_df$phase %in% c("A", "not_possible"), ], aes(xmin = adjust_seg_start, xmax = adjust_seg_end, ymin = loCI, ymax = hiCI, fill = gain_index), color = NA, alpha = 0.15) +
    #         scale_fill_manual(name = "gain order",
    #                           values = colorsA,
    #                           labels = c("", as.character(unique(sort(timing_df[timing_df$phase %in% c("A", "not_possible"), "gain_index"]))[-1])), guide = "none") +
    #         geom_segment(data = timing_df[timing_df$phase %in% c("A", "not_possible"), ],
    #                     aes(x = adjust_seg_start, xend = adjust_seg_end, y = gain_timing, yend = gain_timing, color = gain_index, alpha = phase),
    #                     linewidth = 2) +
    #         scale_color_manual(name = "gain order",
    #                           values = colorsA,
    #                           labels = c("", as.character(unique(sort(timing_df[timing_df$phase %in% c("A", "not_possible"), "gain_index"]))[-1]))) +
    #         scale_alpha_manual(breaks = c("A", "not_possible"), values = c(1, 0), guide = "none") +
    #         scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    #         scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.00\n(early)", "0.25", "0.50", "0.75", "1.00\n(late)"), limits = c(0, 1)) + 
    #         facet_wrap(~haplotype, labeller = labeller(haplotype = c(haplotype = "haplotype A")), strip.position = "right") +
    #         ylab("gain timing") +
    #         xlab("") +
    #         theme_bw() +
    #         theme(axis.text.x = element_blank(),
    #               axis.ticks.x = element_blank(),
    #               panel.grid.major.x = element_blank(),
    #               panel.grid.minor.x = element_blank(),
    #               panel.grid.minor.y = element_blank(),
    #               strip.background = element_rect(fill = "white"),
    #               plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"))
    # 
    # colorsB <- get.colors(max_cn = max(as.numeric(as.character(timing_df[timing_df$phase %in% c("B", "not_possible"), "gain_index"]))))
    # colorsB <- c("white", colorsB)
    # 
    # p_genome2_steps <- ggplot() +
    #         geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    #         geom_rect(data = timing_df[timing_df$phase %in% c("B", "not_possible"), ], aes(xmin = adjust_seg_start, xmax = adjust_seg_end, ymin = loCI, ymax = hiCI, fill = gain_index), color = NA, alpha = 0.15) +
    #         scale_fill_manual(name = "gain order",
    #                           values = colorsB,
    #                           labels = c("", as.character(unique(sort(timing_df[timing_df$phase %in% c("B", "not_possible"), "gain_index"]))[-1])), guide = "none") +
    #         geom_segment(data = timing_df[timing_df$phase %in% c("B", "not_possible"), ],
    #                     aes(x = adjust_seg_start, xend = adjust_seg_end, y = gain_timing, yend = gain_timing, color = gain_index, alpha = phase),
    #                     linewidth = 2) +
    #         scale_color_manual(name = "gain order",
    #                           values = colorsB,
    #                           labels = c("", as.character(unique(sort(timing_df[timing_df$phase %in% c("B", "not_possible"), "gain_index"]))[-1]))) +
    #         scale_alpha_manual(breaks = c("B", "not_possible"), values = c(1, 0), guide = "none") +
    #         scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    #         scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.00\n(early)", "0.25", "0.50", "0.75", "1.00\n(late)"), limits = c(0, 1)) + 
    #         facet_wrap(~haplotype, labeller = labeller(haplotype = c(haplotype = "haplotype B")), strip.position = "right") +
    #         ylab("gain timing") +
    #         xlab("") +
    #         theme_bw() +
    #         theme(panel.grid.major.x = element_blank(),
    #               panel.grid.minor.x = element_blank(),
    #               panel.grid.minor.y = element_blank(),
    #               strip.background = element_rect(fill = "white"),
    #               plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"))
    # 
    # cn_df[cn_df$maj_cn_gain == FALSE & cn_df$min_cn_gain == FALSE, "timing"] <- "not_possible"
    # cn_df[which(cn_df$maj_cn_gain == FALSE & cn_df$phase == "A"), "timing"] <- "not_possible"
    # cn_df[which(cn_df$min_cn_gain == FALSE & cn_df$phase == "B"), "timing"] <- "not_possible"
    # cn_df[is.na(cn_df$phase), "timing"] <- "not_possible"
    # 
    # 
    # # add copy number
    # p_cn <- ggplot() +
    #       geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    #       geom_segment(data = cn_df,
    #                   aes(x = adjust_seg_start, xend = adjust_seg_end, y = maj_cn_lim, yend = maj_cn_lim, color = 'allele A'), linewidth = 2, alpha = 0.7) +
    #       geom_segment(data = cn_df,
    #                   aes(x = adjust_seg_start, xend = adjust_seg_end, y = min_cn_lim, yend = min_cn_lim, color = 'allele B'), linewidth = 2, alpha = 0.7) +
    #       geom_segment(data = cn_df[cn_df$timing == "late" & cn_df$phase == "A", ],
    #                   aes(x = adjust_seg_start, xend = adjust_seg_end, y = local_optimum_lim, yend = local_optimum_lim, color = 'allele A'), linewidth = 1, alpha = 0.2) +
    #       geom_segment(data = cn_df[cn_df$timing == "late" & cn_df$phase == "B", ],
    #                   aes(x = adjust_seg_start, xend = adjust_seg_end, y = local_optimum_lim, yend = local_optimum_lim, color = 'allele B'), linewidth = 1, alpha = 0.2) +
    #       scale_color_manual(name = "haplotype", values = c("allele A" = "darkred",
    #                                                       "allele B" = "#313695"),
    #                                               labels = c("haplotype A", "haplotype B")) +
    #       scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    #       scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
    #       ylab("copy number") +
    #       xlab("") +
    #       theme_bw() +
    #       theme(axis.text.x = element_blank(),
    #             axis.ticks.x = element_blank(),
    #             panel.grid.major.x = element_blank(),
    #             panel.grid.minor.x = element_blank(),
    #             strip.background = element_rect(fill = "white"),
    #             plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm"),
    #             panel.grid.minor.y = element_blank())   
    # 
    # # combine the plots
    # p_genome_steps_c <- p_cn / p_mutA / p_genome1_steps / p_mutB / p_genome2_steps + plot_annotation(title = region, theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
    # 
    # if(!dir.exists(paste0(output_path, tumour, "/", region, "/plots/phased_timing_steps/"))){
    #   dir.create(paste0(output_path, tumour, "/", region, "/plots/phased_timing_steps/"), recursive = T)
    # }
    # 
    # pdf(paste0(output_path, tumour, "/", region, "/plots/phased_timing_steps/", region, "_genome_phased_copynumber_timing_steps.pdf"), width = 18, height = 10)
    # print(p_genome_steps_c)
    # dev.off()
    # 
    # 
    # # and now also do this plot per chromosome
    # for(chr in unique(timing_df$chr)){
    #   
    #   # first, timing plots
    #   timing_df_chr <- timing_df[timing_df$chr == chr, ]
    # 
    #   p_chrA_steps <- ggplot() +
    #                     geom_rect(data = timing_df_chr[timing_df_chr$phase %in% c("A", "not_possible"), ], aes(xmin = seg_start, xmax = seg_end, ymin = loCI, ymax = hiCI, fill = gain_index), color = NA, alpha = 0.15) +
    #                     scale_fill_manual(name = "gain order",
    #                                       values = colorsA,
    #                                       labels = c("", as.character(unique(sort(timing_df_chr[timing_df_chr$phase %in% c("A", "not_possible"), "gain_index"]))[-1])), guide = "none") +
    #                     geom_segment(data = timing_df_chr[timing_df_chr$phase %in% c("A", "not_possible"), ],
    #                                 aes(x = seg_start, xend = seg_end, y = gain_timing, yend = gain_timing, color = gain_index, alpha = phase),
    #                                 linewidth = 2) +
    #                     scale_color_manual(name = "gain order",
    #                                       values = colorsA,
    #                                       labels = c("", as.character(unique(sort(timing_df_chr[timing_df_chr$phase %in% c("A", "not_possible"), "gain_index"]))[-1]))) +
    #                     scale_alpha_manual(breaks = c("A", "not_possible"), values = c(1, 0), guide = "none") +
    #                     scale_x_continuous(expand = c(0,0)) +
    #                     scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.00\n(early)", "0.25", "0.50", "0.75", "1.00\n(late)"), limits = c(0, 1)) + 
    #                     facet_wrap(~haplotype, labeller = labeller(haplotype = c(haplotype = "haplotype A")), strip.position = "right") +
    #                     ylab("gain timing") +
    #                     xlab("") +
    #                     theme_bw() +
    #                     theme(axis.text.x = element_blank(),
    #                           axis.ticks.x = element_blank(),
    #                           panel.grid.major.x = element_blank(),
    #                           panel.grid.minor.x = element_blank(),
    #                           panel.grid.minor.y = element_blank(),
    #                           strip.background = element_rect(fill = "white"),
    #                           plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"))
    # 
    #   p_chrB_steps <- ggplot() +
    #                     geom_rect(data = timing_df_chr[timing_df_chr$phase %in% c("B", "not_possible"), ], aes(xmin = seg_start, xmax = seg_end, ymin = loCI, ymax = hiCI, fill = gain_index), color = NA, alpha = 0.15) +
    #                     scale_fill_manual(name = "gain order",
    #                                       values = colorsB,
    #                                       labels = c("", as.character(unique(sort(timing_df_chr[timing_df_chr$phase %in% c("B", "not_possible"), "gain_index"]))[-1])), guide = "none") +
    #                     geom_segment(data = timing_df_chr[timing_df_chr$phase %in% c("B", "not_possible"), ],
    #                                 aes(x = seg_start, xend = seg_end, y = gain_timing, yend = gain_timing, color = gain_index, alpha = phase),
    #                                 linewidth = 2) +
    #                     scale_color_manual(name = "gain order",
    #                                       values = colorsB,
    #                                       labels = c("", as.character(unique(sort(timing_df_chr[timing_df_chr$phase %in% c("B", "not_possible"), "gain_index"]))[-1]))) +
    #                     scale_alpha_manual(breaks = c("B", "not_possible"), values = c(1, 0), guide = "none") +
    #                     scale_x_continuous(expand = c(0,0)) +
    #                     scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.00\n(early)", "0.25", "0.50", "0.75", "1.00\n(late)"), limits = c(0, 1)) + 
    #                     facet_wrap(~haplotype, labeller = labeller(haplotype = c(haplotype = "haplotype B")), strip.position = "right") +
    #                     ylab("gain timing") +
    #                     xlab("") +
    #                     theme_bw() +
    #                     theme(panel.grid.major.x = element_blank(),
    #                           panel.grid.minor.x = element_blank(),
    #                           panel.grid.minor.y = element_blank(),
    #                           strip.background = element_rect(fill = "white"),
    #                           plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"))
    # 
    #    # and plot the copy number
    #   cn_df_chr <- cn_df[cn_df$chr == chr, ]
    #     
    #   p_chr_cn <- ggplot() +
    #             geom_segment(data = cn_df_chr,
    #                         aes(x = seg_start, xend = seg_end, y = maj_cn_lim, yend = maj_cn_lim, color = 'allele A'), linewidth = 2, alpha = 0.7) +
    #             geom_segment(data = cn_df_chr,
    #                         aes(x = seg_start, xend = seg_end, y = min_cn_lim, yend = min_cn_lim, color = 'allele B'), linewidth = 2, alpha = 0.7) +
    #             geom_segment(data = cn_df_chr[cn_df_chr$timing == "late" & cn_df_chr$phase == "A", ],
    #                         aes(x = seg_start, xend = seg_end, y = local_optimum_lim, yend = local_optimum_lim, color = 'allele A'), linewidth = 1, alpha = 0.2) +
    #             geom_segment(data = cn_df_chr[cn_df_chr$timing == "late" & cn_df_chr$phase == "B", ],
    #                         aes(x = seg_start, xend = seg_end, y = local_optimum_lim, yend = local_optimum_lim, color = 'allele B'), linewidth = 1, alpha = 0.2) +
    #             scale_color_manual(name = "haplotype", values = c("allele A" = "darkred",
    #                                                             "allele B" = "#313695"),
    #                                                     labels = c("haplotype A", "haplotype B")) +
    #             scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    #             scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) +
    #             ylab("copy number") +
    #             xlab("") +
    #             theme_bw() +
    #             theme(axis.text.x = element_blank(),
    #                   axis.ticks.x = element_blank(),
    #                   panel.grid.major.x = element_blank(),
    #                   panel.grid.minor.x = element_blank(),
    #                   strip.background = element_rect(fill = "white"),
    #                   plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm"),
    #                   panel.grid.minor.y = element_blank())   
    # 
    #   # and plot the mutations
    #   mut_df_seg_chr <- mut_df_seg[mut_df_seg$chr == chr, ]
    # 
    #   mut_cpn_ylim_chr_min <- 0
    #   mut_cpn_ylim_chr_max <- max(mut_df_seg_chr$mut_cpn_lim, mut_df_seg_chr$maj_cn_lim, mut_df_seg_chr$min_cn_lim)
    # 
    #   if(nrow(mut_df_seg_chr[which(mut_df_seg_chr$phase == "A"), ]) > 0){
    #     p_chr_mutA <- ggplot(mut_df_seg_chr[which(mut_df_seg_chr$phase == "A"), ], aes(adjust_pos, mut_cpn_lim, fill = color_index), color = "black") +
    #                 geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
    #                 scale_fill_manual(name = "", breaks = c("one", "two", "not_possible"), labels = c("one", "two", "not possible"), values = c("#d8b365", "#5ab4ac", "#bababa"), guide = "none") +
    #                 geom_segment(data = unique(mut_df_seg_chr[, c("adjust_seg_start", "adjust_seg_end", "allelic_imbalance", "phase", "maj_cn_lim", "min_cn_lim", "color_index")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = maj_cn_lim, yend = maj_cn_lim, color = 'allele A', alpha = allelic_imbalance), linewidth = 1) +
    #                 scale_color_manual(name = "copy number phase", values = c("allele A" = "darkred"), guide = "none") +
    #                 scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
    #                 scale_x_continuous(expand = c(0,0)) +
    #                 facet_wrap(~haplotype, labeller = labeller(haplotype = c("haplotypeA" = "haplotype A")), strip.position = "right") +
    #                 ylab("mutation copy number") +
    #                 xlab("") +
    #                 ylim(mut_cpn_ylim_chr_min, mut_cpn_ylim_chr_max) +
    #                 theme_bw() +
    #                 theme(axis.text.x = element_blank(),
    #                       axis.ticks.x = element_blank(),
    #                       panel.grid.major.x = element_blank(),
    #                       panel.grid.minor.x = element_blank(),
    #                       strip.background = element_rect(fill = "white"))
    #   } else {
    #     p_chr_mutA <- ggplot(unique(mut_df_seg_chr[, c("adjust_seg_start", "adjust_seg_end", "allelic_imbalance", "phase", "maj_cn_lim", "min_cn_lim", "color_index")])) +
    #         geom_segment(data = unique(mut_df_seg_chr[, c("adjust_seg_start", "adjust_seg_end", "allelic_imbalance", "phase", "maj_cn_lim", "min_cn_lim", "color_index")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = maj_cn_lim, yend = maj_cn_lim, color = 'allele A', alpha = allelic_imbalance), linewidth = 1) +
    #         scale_color_manual(name = "copy number phase", values = c("allele A" = "darkred"), guide = "none") +
    #         scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
    #         scale_x_continuous(expand = c(0,0)) +
    #         facet_wrap(~phase, labeller = labeller(phase = c("A" = "haplotype A")), strip.position = "right") +
    #         ylab("mutation copy number") +
    #         xlab("") +
    #         ylim(mut_cpn_ylim_chr_min, mut_cpn_ylim_chr_max) +
    #         theme_bw() +
    #         theme(axis.text.x = element_blank(),
    #               axis.ticks.x = element_blank(),
    #               panel.grid.major.x = element_blank(),
    #               panel.grid.minor.x = element_blank(),
    #               strip.background = element_rect(fill = "white"))
    #   }
    #   
    #   if(nrow(mut_df_seg_chr[which(mut_df_seg_chr$phase == "B"), ]) > 0){
    #     p_chr_mutB <- ggplot(mut_df_seg_chr[which(mut_df_seg_chr$phase == "B"), ], aes(adjust_pos, mut_cpn_lim, fill = color_index), color = "black") +
    #                 geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
    #                 scale_fill_manual(name = "", breaks = c("one", "two", "not_possible"), labels = c("one", "two", "not possible"), values = c("#d8b365", "#5ab4ac", "#bababa"), guide = "none") +
    #                 geom_segment(data = unique(mut_df_seg_chr[, c("adjust_seg_start", "adjust_seg_end", "allelic_imbalance", "phase", "maj_cn_lim", "min_cn_lim", "color_index")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = min_cn_lim, yend = min_cn_lim, color = 'allele B', alpha = allelic_imbalance), linewidth = 1) +
    #                 scale_color_manual(name = "copy number phase", values = c("allele B" = "#313695"), guide = "none") +
    #                 scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
    #                 scale_x_continuous(expand = c(0,0)) +
    #                 facet_wrap(~haplotype, labeller = labeller(haplotype = c("haplotypeB" = "haplotype B")), strip.position = "right") +
    #                 ylab("mutation copy number") +
    #                 xlab("") +
    #                 ylim(mut_cpn_ylim_chr_min, mut_cpn_ylim_chr_max) +
    #                 theme_bw() +
    #                 theme(legend.position = "none",
    #                       axis.text.x = element_blank(),
    #                       axis.ticks.x = element_blank(),
    #                       panel.grid.major.x = element_blank(),
    #                       panel.grid.minor.x = element_blank(),
    #                       strip.background = element_rect(fill = "white"))                
    #   } else {
    #     p_chr_mutB <- ggplot(unique(mut_df_seg_chr[, c("adjust_seg_start", "adjust_seg_end", "allelic_imbalance", "phase", "maj_cn_lim", "min_cn_lim", "color_index")])) +
    #         geom_segment(data = unique(mut_df_seg_chr[, c("adjust_seg_start", "adjust_seg_end", "allelic_imbalance", "phase", "maj_cn_lim", "min_cn_lim", "color_index")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = min_cn_lim, yend = min_cn_lim, color = 'allele B', alpha = allelic_imbalance), linewidth = 1) +
    #         scale_color_manual(name = "copy number phase", values = c("allele B" = "#313695"), guide = "none") +
    #         scale_alpha_manual(name = "allelic imbalance", values = c(0.7, 0.3), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
    #         scale_x_continuous(expand = c(0,0)) +
    #         facet_wrap(~phase, labeller = labeller(phase = c("B" = "haplotype B")), strip.position = "right") +
    #         ylab("mutation copy number") +
    #         xlab("") +
    #         ylim(mut_cpn_ylim_chr_min, mut_cpn_ylim_chr_max) +
    #         theme_bw() +
    #         theme(legend.position = "none",
    #               axis.text.x = element_blank(),
    #               axis.ticks.x = element_blank(),
    #               panel.grid.major.x = element_blank(),
    #               panel.grid.minor.x = element_blank(),
    #               strip.background = element_rect(fill = "white"))    
    #   }
    # 
    #   p_chr_steps_c <- p_chr_cn / p_chr_mutA / p_chrA_steps / p_chr_mutB / p_chrB_steps + plot_annotation(title = paste0(region, " : ", chr), theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
    #   
    #   if(!dir.exists(paste0(output_path, tumour, "/", region, "/plots/phased_timing_steps/per_chromosome/"))){
    #     dir.create(paste0(output_path, tumour, "/", region, "/plots/phased_timing_steps/per_chromosome/"), recursive = T)
    #   }
    #   
    #   pdf(paste0(output_path, tumour, "/", region, "/plots/phased_timing_steps/per_chromosome/", region, "_", chr, "_phased_copynumber_timing_steps.pdf"), width = 10, height = 10)
    #   print(p_chr_steps_c)
    #   dev.off()
    # 
    # } 
    #   
}

plot.early.late.fun <- function(mutation_df = phased_data, genome = BSgenome.Hsapiens.UCSC.hg38, output_path){
  
  mutation_df <- data.frame(mutation_df)

  # first, make the genome continous
  mutation_df <- genome.continous.fun(mutation_df, genome)
  
  chr_sizes <- seqlengths(genome)[paste0("chr", 1:22)]
  chr_lines <- cumsum(as.numeric(as.character(chr_sizes)))
  chr_midpoint <- chr_lines - chr_sizes / 2
  
  xaxis_label <- data.frame(chr_midpoint)
  xaxis_label$label <- rownames(xaxis_label)
  
  timing_df <- unique(mutation_df[, c("segID", "chr", "pos", "adjust_pos", "seg_start", "seg_end", "adjust_seg_start", "adjust_seg_end", "maj_cn", "min_cn", "mut_cpn", "phase", "timing", "allelic_imbalance", "maj_cn_gain", "min_cn_gain")])
  timing_df$haplotype <- paste0("haplotype ", timing_df$phase)
  
  # grey out the mutations that we can't time, ie those where there is no gain or where there is no allelic imbalance
  timing_df[timing_df$maj_cn_gain == FALSE & timing_df$min_cn_gain == FALSE, "timing"] <- "not_possible"
  timing_df[which(timing_df$maj_cn_gain == FALSE & timing_df$phase == "A"), "timing"] <- "not_possible"
  timing_df[which(timing_df$min_cn_gain == FALSE & timing_df$phase == "B"), "timing"] <- "not_possible"
  timing_df[is.na(timing_df$phase), "timing"] <- "not_possible"
  timing_df$timing <- factor(timing_df$timing, levels = c("early", "late", "not_possible"))
  
  timing_df[timing_df$mut_cpn > 10, "mut_cpn_lim"] <- 10
  timing_df[timing_df$mut_cpn <= 10, "mut_cpn_lim"] <- timing_df[timing_df$mut_cpn <= 10, "mut_cpn"] 
  timing_df[timing_df$min_cn > 10, "seg_y_pos_min_lim"] <- 10
  timing_df[timing_df$min_cn <= 10, "seg_y_pos_min_lim"] <- timing_df[timing_df$min_cn <= 10, "min_cn"]
  timing_df[timing_df$maj_cn > 10, "seg_y_pos_maj_lim"] <- 10
  timing_df[timing_df$maj_cn <= 10, "seg_y_pos_maj_lim"] <- timing_df[timing_df$maj_cn <= 10, "maj_cn"]
  
  mut_cpn_ylim_genome <- max(c(timing_df$mut_cpn_lim, timing_df$seg_y_pos_maj_lim))
  mut_cpn_ylim_genome_min <- min(c(timing_df$mut_cpn_lim, timing_df$seg_y_pos_min_lim))

  p_genome1 <- ggplot(timing_df[which(timing_df$phase == "A"), ], aes(adjust_pos, mut_cpn_lim, fill = timing), color = "black") +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
    scale_fill_manual(name = "mutation timing", breaks = c("early", "late", "not_possible"), values = c("#e94e1b", "#662483", "#bababa"), labels = c("early", "late", "not possible")) +
    geom_segment(data = unique(timing_df[, c("adjust_seg_start", "adjust_seg_end", "phase", "seg_y_pos_maj_lim", "seg_y_pos_min_lim", "timing")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_maj_lim, yend = seg_y_pos_maj_lim, color = 'haplotype A'), alpha = 0.7, linewidth = 2) +
    scale_color_manual(name = "", values = c("haplotype A" = "darkred"), label = "CN haplotype A") +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    facet_wrap(~haplotype, ncol = 1, strip.position = "right") +
    #ylim(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome) +
    ylab("mutation copy number") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm"))
  
  p_genome2 <- ggplot(timing_df[which(timing_df$phase == "B"), ], aes(adjust_pos, mut_cpn_lim, fill = timing), color = "black") +
    geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
    geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
    scale_fill_manual(name = "mutation timing", breaks = c("early", "late", "not_possible"), values = c("#e94e1b", "#662483", "#bababa"), labels = c("early", "late", "not possible"), guide = "none") +
    geom_segment(data = unique(timing_df[, c("adjust_seg_start", "adjust_seg_end", "phase", "seg_y_pos_maj_lim", "seg_y_pos_min_lim", "timing")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_min_lim, yend = seg_y_pos_min_lim, color = 'haplotype B'), alpha = 0.7, linewidth = 2) +
    scale_color_manual(name = "", values=c("haplotype B" = "#313695"), label = "CN haplotype B") +
    scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
    facet_wrap(~haplotype, ncol = 1, strip.position = "right") +
    #ylim(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome) +
    ylab("mutation copy number") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.background = element_rect(fill = "white"),
          plot.margin = unit(c(0, 0.1, 0.1, 0.1), "cm"))
  
  # combine the two plots
  p_genome_c <- plot_grid(p_genome1, p_genome2, ncol = 1, align = "v", axis = "rl", rel_heights = c(1, 1.1))
  
  if(!dir.exists(paste0(output_path, tumour, "/", region, "/plots/phased_mutation_timing"))){
    dir.create(paste0(output_path, tumour, "/", region, "/plots/phased_mutation_timing"), recursive = T)
  }
  
  pdf(paste0(output_path, tumour, "/", region, "/plots/phased_mutation_timing/", region, "_genome_phased_mutation_timing_per_haplotype_mutcnp.pdf"), width = 18, height = 7)
  print(p_genome_c)
  dev.off()
  
  
  # now make one plot for each chromosome
  for(chr in unique(timing_df$chr)){
    
    chr_timing_df <- timing_df[timing_df$chr == chr, ]
    
    mut_cpn_ylim_chr <- max(c(chr_timing_df$mut_cpn_lim, chr_timing_df$seg_y_pos_maj_lim))
    mut_cpn_ylim_chr_min <- min(c(chr_timing_df$mut_cpn_lim, chr_timing_df$seg_y_pos_min_lim))
    
    if(nrow(chr_timing_df[which(chr_timing_df$phase == "A"), ]) > 0){
      p_chr1 <- ggplot(chr_timing_df[which(chr_timing_df$phase == "A"), ], aes(pos, mut_cpn_lim, fill = timing), color = "black") +
        geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
        scale_fill_manual(name = "mutation timing", breaks = c("early", "late", "not_possible"), values = c("#e94e1b", "#662483", "#bababa"), labels = c("early", "late", "not possible")) +
        geom_segment(data = unique(chr_timing_df[, c("seg_start", "seg_end", "phase", "seg_y_pos_maj_lim", "seg_y_pos_min_lim", "timing")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_maj_lim, yend = seg_y_pos_maj_lim, color = 'haplotype A'), alpha = 0.7, linewidth = 2) +
        scale_color_manual(name = "", values = c("haplotype A" = "darkred"), label = "CN haplotype A") +
        facet_wrap(~haplotype, ncol = 1, strip.position = "right") +
        ylim(mut_cpn_ylim_chr_min, mut_cpn_ylim_chr) +
        ylab("mutation copy number") +
        xlab("") +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              strip.background = element_rect(fill = "white"),
              plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm"))
    }
    
    if(nrow(chr_timing_df[which(chr_timing_df$phase == "B"), ]) > 0){
      p_chr2 <- ggplot(chr_timing_df[which(chr_timing_df$phase == "B"), ], aes(pos, mut_cpn_lim, fill = timing), color = "black") +
        geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
        scale_fill_manual(name = "mutation timing", breaks = c("early", "late", "not_possible"), values = c("#e94e1b", "#662483", "#bababa"), labels = c("early", "late", "not possible"), guide = "none") +
        geom_segment(data = unique(chr_timing_df[, c("seg_start", "seg_end", "phase", "seg_y_pos_maj_lim", "seg_y_pos_min_lim", "timing")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_min_lim, yend = seg_y_pos_min_lim, color = 'haplotype B'), alpha = 0.7, linewidth = 2) +
        scale_color_manual(name = "", values=c("haplotype B" = "#313695"), label = "CN haplotype B") +
        facet_wrap(~haplotype, ncol = 1, strip.position = "right") +
        ylim(mut_cpn_ylim_chr_min, mut_cpn_ylim_chr) +
        ylab("mutation copy number") +
        xlab(chr) +
        theme_bw() +
        theme(strip.background = element_rect(fill = "white"),
              plot.margin = unit(c(0, 0.1, 0.1, 0.1), "cm"))
    }
    
    # combine the two plots
    if(all(c(nrow(chr_timing_df[which(chr_timing_df$phase == "A"), ]) > 0, nrow(chr_timing_df[which(chr_timing_df$phase == "B"), ]) > 0))){
      
      p_chr_c <- plot_grid(p_chr1, p_chr2, ncol = 1, align = "v", axis = "rl", rel_heights = c(1, 1.1))
      
      if(!dir.exists(paste0(output_path, tumour, "/", region, "/plots/phased_mutation_timing/per_chromosome"))){
        dir.create(paste0(output_path, tumour, "/", region, "/plots/phased_mutation_timing/per_chromosome"), recursive = T)
      }
      
      pdf(paste0(output_path, tumour, "/", region, "/plots/phased_mutation_timing/per_chromosome/", region, "_", chr, "_phased_mutation_timing_per_haplotype_mutcnp.pdf"), width = 10, height = 7)
      print(p_chr_c)
      dev.off()
    
    } else {
      
      if(nrow(chr_timing_df[which(chr_timing_df$phase == "A"), ]) > 0){
       
        p_chr1 <- ggplot(chr_timing_df[which(chr_timing_df$phase == "A"), ], aes(pos, mut_cpn_lim, fill = timing), color = "black") +
          geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
          scale_fill_manual(name = "mutation timing", breaks = c("early", "late", "not_possible"), values = c("#e94e1b", "#662483", "#bababa"), labels = c("early", "late", "not possible")) +
          geom_segment(data = unique(chr_timing_df[, c("seg_start", "seg_end", "phase", "seg_y_pos_maj_lim", "seg_y_pos_min_lim", "timing")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_maj_lim, yend = seg_y_pos_maj_lim, color = 'haplotype A'), alpha = 0.7, linewidth = 2) +
          scale_color_manual(name = "", values = c("haplotype A" = "darkred"), label = "CN haplotype A") +
          facet_wrap(~haplotype, ncol = 1, strip.position = "right") +
          ylim(mut_cpn_ylim_chr_min, mut_cpn_ylim_chr) +
          ylab("mutation copy number") +
          xlab(chr) +
          theme_bw() +
          theme(strip.background = element_rect(fill = "white"),
                plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
        
        if(!dir.exists(paste0(output_path, tumour, "/", region, "/plots/phased_mutation_timing/per_chromosome"))){
          dir.create(paste0(output_path, tumour, "/", region, "/plots/phased_mutation_timing/per_chromosome"), recursive = T)
        }
        
        pdf(paste0(output_path, tumour, "/", region, "/plots/phased_mutation_timing/per_chromosome/", region, "_", chr, "_phased_mutation_timing_per_haplotype_mutcnp.pdf"), width = 10, height = 7)
        print(p_chr1)
        dev.off()
      }  
      
      if(nrow(chr_timing_df[which(chr_timing_df$phase == "B"), ]) > 0){
        
        p_chr2 <- ggplot(chr_timing_df[which(chr_timing_df$phase == "B"), ], aes(pos, mut_cpn_lim, fill = timing), color = "black") +
          geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
          scale_fill_manual(name = "mutation timing", breaks = c("early", "late", "not_possible"), values = c("#e94e1b", "#662483", "#bababa"), labels = c("early", "late", "not possible")) +
          geom_segment(data = unique(chr_timing_df[, c("seg_start", "seg_end", "phase", "seg_y_pos_maj_lim", "seg_y_pos_min_lim", "timing")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_min_lim, yend = seg_y_pos_min_lim, color = 'haplotype B'), alpha = 0.7, linewidth = 2) +
          scale_color_manual(name = "", values=c("haplotype B" = "#313695"), label = "CN haplotype B") +
          facet_wrap(~haplotype, ncol = 1, strip.position = "right") +
          ylim(mut_cpn_ylim_chr_min, mut_cpn_ylim_chr) +
          ylab("mutation copy number") +
          xlab(chr) +
          theme_bw() +
          theme(strip.background = element_rect(fill = "white"),
                plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
        
        if(!dir.exists(paste0(output_path, tumour, "/", region, "/plots/phased_mutation_timing/per_chromosome"))){
          dir.create(paste0(output_path, tumour, "/", region, "/plots/phased_mutation_timing/per_chromosome"), recursive = T)
        }
        
        pdf(paste0(output_path, tumour, "/", region, "/plots/phased_mutation_timing/per_chromosome/", region, "_", chr, "_phased_mutation_timing_per_haplotype_mutcnp.pdf"), width = 10, height = 7)
        print(p_chr2)
        dev.off()
      }  
      
    }
    
  }
  
}

################################################################################
#####################################################################     MAIN

# test 
# mutation_data_all <- read.table(mutation_data_base_path, head = T, sep = "\t")
#       
# all_timed_data <- list()
# for(region in all_regions){
#   
#   print(region)
#   
#   # get data
#   mutation_data <- mutation_data_all[mutation_data_all$sample == region, ]
#   
#   if(nrow(mutation_data_all[mutation_data_all$sample == region, ]) > 0){
#     
#     mutation_data <- mutation_data[, c("sample", "chr", "pos", "hap1_id", "hap1_bases", "hap2_bases", "hap_id_SNP_pos", "phase","maj_cn", "min_cn", "tot_cn", "seg_start", "seg_end", "mut_cpn")]
#     
#     snp_data = read.table(paste0(snp_data_base_path, region, "_phased_muts_long.tsv"), head = T, sep = "\t")
#     snp_data$mutation_snp_id <- paste0(snp_data$chr_mut, "_", snp_data$pos_mut, "_", snp_data$pos_hetsnp)
#     
#     tumour = strsplit(region, "_")[[1]][1]
#     BAF_data <- lapply(paste0("chr", c(1:22, "X")), function(chr){
#       gr <- readRDS(paste0(battenberg_base_path, tumour, "/Battenberg/aggregated/", tumour, "_BS_GL_multisample_haplotypes_", chr, "_loci.RDS"))
#     })
#     
#     BAF_data <- do.call(c, BAF_data)
#     
#     # run the functions
#     print("getting data in shape")
#     data <- wrangle.data.fun(region, tumour, snp_data, mutation_data, BAF_data)
#     
#     print("phasing mutations")
#     phased_data <- phasing.phun(data, tumour, region, output_path)
#     
#     # print("plotting phased mutations")
#     # phun.plots(phased_data, tumour, region, output_path, genome)
#     # 
#     print("estimating mutation timing")
#     phased_data <- mut.cpn.match.fun(phased_data)
#     
#     print("calculating allele specific copy number timing estimates")
#     phased_timed_data <- calc.timing.fun(phased_data)
#     
#     # print("plotting phased mutations by CN state")
#     # copynumber.state.plot.fun(phased_data, opacity = 255, tumour, region, output_path, genome)
#     # 
#     # cn_states_data <- compare.cn.states.fun(phased_data, region)
#     
#     all_timed_data[[region]] <- phased_timed_data
#     
#   }
# 
# }
# 
# all_timed_data <- do.call(rbind, all_timed_data)


# TO DO: make this function work
# write a function that runs all functions sequentially
# phunction <- function(tumourID, regionID, mutation_df, snp_df, BAF_df, genome_object = BSgenome.Hsapiens.UCSC.hg38, output_dir){
#   
#   print(paste0("phasing mutations for tumour: ", tumourID, " - region ", regionID))
#   # tumour <- tumourID
#   # region <- regionID
#   # mutation_data <- mutation_df
#   # snp_data <- snp_df
#   # BAF_data <- BAF_df
#   # genome <- genome_object
#   # output_path <- output_dir
#   # 
#   print("getting data in shape")
#   data <- wrangle.data.fun(region = regionID, tumour = tumourID, snp_data = snp_df, mutation_data = mutation_df, BAF_data = BAF_df)
#   
#   print("phasing mutations")
#   phased_data <- phasing.phun(data, tumour, region, output_path)
#   print("plotting phased mutations")
#   phun.plots(phased_data, tumour, region, output_path, genome)
#   
#   print("plotting phased mutations by CN state")
#   phased_data <- mut.cpn.match.fun(phased_data)
#   copynumber.state.plot.fun(phased_data, opacity = 255, tumour, region, output_path, genome)
#   
# }
# 

