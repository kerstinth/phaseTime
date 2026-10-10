# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #                   allele specific mutation timing                     # # #  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################################################################################
#####################################################################  libraries
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(data.table)
library(vcfR)
library(BSgenome.Hsapiens.UCSC.hg38)

options(scipen=999)
################################################################################
##################################################################### file paths


genome                     <- BSgenome.Hsapiens.UCSC.hg38

################################################################################
#################################################################### functions

genome.continous.fun <- function(data, genome){
  
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

chr_sizes <- seqlengths(genome)[paste0("chr", 1:22)]
chr_lines <- cumsum(as.numeric(as.character(chr_sizes)))
chr_midpoint <- chr_lines - chr_sizes / 2

################################################################################
#####################################################################     MAIN

mutation_data_all <- read.table(mutation_data_base_path, head = T, sep = "\t")

for(region in all_regions){
  
  tumour <- sub("_.*", "", region)
  print(tumour)
  print(region)
  
  snp_data_path <- paste0(snp_data_base_path, region, "_phased_muts_long.tsv")
  
  if(file.exists(snp_data_path)){
    
    if(file.exists(paste0(battenberg_base_path, tumour, "/Battenberg/aggregated/", tumour, "_BS_GL_multisample_haplotypes_chr1_loci.RDS"))){
      
      if(region %in% unique(mutation_data_all$sample)){
        
        snp_data <- read.table(snp_data_path, head = T, sep = "\t")
        snp_data$mutation_snp_id <- paste0(snp_data$chr_mut, "_", snp_data$pos_mut, "_", snp_data$pos_hetsnp)
        
        mutation_data <- mutation_data_all[mutation_data_all$sample == region, ]
        mutation_data <- mutation_data[, c("sample", "chr", "pos", "hap1_id", "hap1_bases", "hap2_bases", "hap_id_SNP_pos", "phase","maj_cn", "min_cn", "tot_cn", "seg_start", "seg_end", "mut_cpn")]
        
        # now let's get this in long format to get all the separate SNPs that are linked to a mutation
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
        mutation_data_long <- left_join(mutation_data_long, mutation_data[, c("sample", "chr", "pos", "maj_cn", "min_cn", "tot_cn", "seg_start", "seg_end", "mut_cpn", "phase","mutation_id")])
        mutation_data_long$phase <- NULL
        
        mutation_data_long$chr <- paste0("chr", mutation_data_long$chr)
        mutation_data_long$mutation_snp_id <- paste0(mutation_data_long$chr, "_", mutation_data_long$pos, "_", mutation_data_long$hap_id_SNP_pos)
        
        snp_data_small <- snp_data[, c("mutation_snp_id", "ref_mut", "alt_mut", "counts_mut", "chr_hetsnp", "pos_hetsnp", "ref_hetsnp", "alt_hetsnp", "counts_hetsnp", 
                                       "snp_vaf_on_mut")]
        
        snp_data_small <- unique(snp_data_small)
        
        snp_mutation_data <- left_join(mutation_data_long, snp_data_small)
        snp_mutation_data$snp_id <- paste0(snp_mutation_data$chr, "_", snp_mutation_data$pos_hetsnp)
        
        chr_df <- list()
        for(chr in c(paste0("chr", 1:22))){
          print(chr)
          
          # load battenberg BAF data
          BAF_data <- readRDS(paste0(battenberg_base_path, tumour, "/Battenberg/aggregated/", tumour, "_BS_GL_multisample_haplotypes_", chr, "_loci.RDS"))
          BAF_region <- BAF_data[, c("REF", "ALT", paste0(region, "_Major"), paste0(region, "_BAF"), paste0(region, "_PS"), "multisample_haplo")]
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
          
          if(nrow(snp_mutation_df) > 0){
            snp_mutation_df <- snp_mutation_df[!is.na(snp_mutation_df$snp_vaf_on_mut), ]
            
            snp_mutation_df[snp_mutation_df$bb_BAF < 0.5, "bb_phase"] <- "B"
            snp_mutation_df[snp_mutation_df$bb_BAF >= 0.5, "bb_phase"] <- "A"
            
            snp_mutation_df[snp_mutation_df$bb_phase == "B", "snp_vaf_on_mut_bb_phase"] <- abs(1 - snp_mutation_df[snp_mutation_df$bb_phase == "B", "snp_vaf_on_mut"])
            snp_mutation_df[snp_mutation_df$bb_phase == "A", "snp_vaf_on_mut_bb_phase"] <- snp_mutation_df[snp_mutation_df$bb_phase == "A", "snp_vaf_on_mut"] 
            
            snp_mutation_df$mutation_id <- paste0(snp_mutation_df$chr, "_", snp_mutation_df$pos)                
            snp_mutation_df <- snp_mutation_df %>% group_by(mutation_id) %>% mutate(average_snp_vaf_on_mut_bb_phase = mean(snp_vaf_on_mut_bb_phase))                
            snp_mutation_df <- data.frame(snp_mutation_df)             
            
            mutation_df <- snp_mutation_df           
            mutation_df <- unique(mutation_df[, c("sample", "chr", "pos", "maj_cn", "min_cn", "tot_cn", "seg_start", "seg_end", "mut_cpn",
                                                  "ref_mut", "alt_mut", "counts_mut", "average_snp_vaf_on_mut_bb_phase")])
            
            mutation_df[mutation_df$average_snp_vaf_on_mut_bb_phase >= 0.5, "phase"] <- "A"
            mutation_df[mutation_df$average_snp_vaf_on_mut_bb_phase < 0.5, "phase"] <- "B"
            mutation_df$phase <- factor(mutation_df$phase, levels = c("A", "B"))
            
            mutation_df <- mutation_df[!is.na(mutation_df$maj_cn), ]
            mutation_df <- mutation_df[!is.na(mutation_df$min_cn), ]
            mutation_df <- mutation_df[!is.na(mutation_df$mut_cpn), ]
            
            # we need to make sure that the two copy number alleles are not plotted on top of each other
            # but we also don't want to move them too much so we can still distinguish subclonality
            
            mutation_df[mutation_df$maj_cn == mutation_df$min_cn, "seg_y_pos_maj"] <- mutation_df[mutation_df$maj_cn == mutation_df$min_cn, "maj_cn"] + 0.025
            mutation_df[mutation_df$maj_cn != mutation_df$min_cn, "seg_y_pos_maj"] <- mutation_df[mutation_df$maj_cn != mutation_df$min_cn, "maj_cn"]
            mutation_df[mutation_df$maj_cn == mutation_df$min_cn, "seg_y_pos_min"] <- mutation_df[mutation_df$maj_cn == mutation_df$min_cn, "min_cn"] - 0.025
            mutation_df[mutation_df$maj_cn != mutation_df$min_cn, "seg_y_pos_min"] <- mutation_df[mutation_df$maj_cn != mutation_df$min_cn, "min_cn"]
            
            # indicate where there is allelic imbalance
            mutation_df[mutation_df$min_cn != mutation_df$maj_cn, "allelic_imbalance"] <- "TRUE"
            mutation_df[mutation_df$min_cn == mutation_df$maj_cn, "allelic_imbalance"] <- "FALSE"
            
            mutation_df$allelic_imbalance <- factor(mutation_df$allelic_imbalance, levels = c("TRUE", "FALSE"))
            
            # make sure our segments also have a legend
            mutation_df$seg_legendA <- "allele A"
            mutation_df$seg_legendB <- "allele B"
            
            colnames(mutation_df)[c((ncol(mutation_df) - 1) : ncol(mutation_df))] <- c("allele A", "allele B")
            
            chr_df[[chr]] <- mutation_df
            
            # now, for every segment we can time the gain, or the mutations 
            # mutation_df$segID <- paste0(mutation_df$chr, "_", mutation_df$seg_start, "_", mutation_df$seg_end)
            # 
            # on top of the phase, for every mutation we need to determine whether the mutation sits on the major or minor copy number allele
            # is the multiplicity closer to the major or to the minor allele?
            # mutation_df$multiplicity_maj_distance <- abs(mutation_df$maj_cn - mutation_df$mut_cpn)
            # mutation_df$multiplicity_min_distance <- abs(mutation_df$min_cn - mutation_df$mut_cpn)
            # mutation_df[mutation_df$multiplicity_maj_distance < mutation_df$multiplicity_min_distance, "copy_number_allele"] <- "major"
            # mutation_df[mutation_df$multiplicity_min_distance < mutation_df$multiplicity_maj_distance, "copy_number_allele"] <- "minor"
            # mutation_df[is.na(mutation_df$copy_number_allele), "copy_number_allele"] <- "unknown"
            # 
            # mutation_df$segID_allele <- paste0(mutation_df$segID, "_", mutation_df$copy_number_allele)
            
            # for every segment and each phase, calculate the ratio of mutations on the major and the minor allele
            # counts <- mutation_df %>% group_by(segID) %>% count(phase, copy_number_allele) %>% complete(phase, copy_number_allele, fill = list(n = 0))
            # counts <- dcast(counts, segID + phase ~ copy_number_allele, value.var = "n")
            # counts$ratio <- counts$minor / counts$major
            
            # counts <- mutation_df %>% group_by(segID) %>% count(phase) %>% complete(phase, fill = list(n = 0))
            # counts <- reshape2::dcast(counts, segID ~ phase, value.var = "n")
            # counts$ratio <- counts$A / counts$B
            # 
            # chr_timing_df[[chr]] <- counts
            
            p <- ggplot(mutation_df, aes(pos, mut_cpn, fill = phase), color = "black") +
              geom_segment(data = unique(mutation_df[, c("seg_start", "seg_end", "seg_y_pos_maj", "seg_y_pos_min", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_maj, yend = seg_y_pos_maj, color = 'allele A', alpha = allelic_imbalance), linewidth = 2) +
              geom_segment(data = unique(mutation_df[, c("seg_start", "seg_end", "seg_y_pos_maj", "seg_y_pos_min", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_min, yend = seg_y_pos_min, color = 'allele B', alpha = allelic_imbalance),linewidth = 2) +
              geom_point(alpha = 0.7, shape = 21, stroke = NA) +
              scale_fill_manual(name = "mutation phase", breaks = c("A", "B"), labels = c("allele A", "allele B"), values = c("#d73027", "#4575b4")) +
              scale_color_manual(name = "copy number phase", values=c("allele A" = "darkred", "allele B" = "#313695")) +
              scale_alpha_manual(name = "allelic imbalance", values = c(1, 0.5), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
              ylab("mutation copy number") +
              # ylim(0, mut_cpn_ylim) +
              xlab(chr) +
              theme_bw() 
            
            p_red <- ggplot(mutation_df[mutation_df$phase == "A", ], aes(pos, mut_cpn, fill = phase), color = "black") +
              geom_segment(data = unique(mutation_df[, c("seg_start", "seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_maj, yend = seg_y_pos_maj, color = 'allele A', alpha = allelic_imbalance), linewidth = 2) +
              geom_segment(data = unique(mutation_df[, c("seg_start", "seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_min, yend = seg_y_pos_min, color = 'allele B', alpha = allelic_imbalance),linewidth = 2) +
              geom_point(alpha = 0.7, shape = 21, stroke = NA) +
              scale_fill_manual(name = "mutation phase", breaks = c("A", "B"), labels = c("allele A", "allele B"), values = c("#d73027", "#4575b4")) +
              scale_color_manual(name = "copy number phase", values=c("allele A" = "darkred", "allele B" = "#313695")) +
              scale_alpha_manual(name = "allelic imbalance", values = c(1, 0.5), breaks = c(FALSE, TRUE), label = c("FALSE", "TRUE")) +
              ylab("mutation copy number") +
              # ylim(0, mut_cpn_ylim) +
              xlab(chr) +
              theme_bw() 
            
            p_blue <- ggplot(mutation_df[mutation_df$phase == "B", ], aes(pos, mut_cpn, fill = phase), color = "black") +
              geom_segment(data = unique(mutation_df[, c("seg_start", "seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_maj, yend = seg_y_pos_maj, color = 'allele A', alpha = allelic_imbalance), linewidth = 2) +
              geom_segment(data = unique(mutation_df[, c("seg_start", "seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_min, yend = seg_y_pos_min, color = 'allele B', alpha = allelic_imbalance),linewidth = 2) +
              geom_point(alpha = 0.7, shape = 21, stroke = NA) +
              scale_fill_manual(name = "mutation phase", breaks = c("A", "B"), labels = c("allele A", "allele B"), values = c("#d73027", "#4575b4")) +
              scale_color_manual(name = "copy number phase", values=c("allele A" = "darkred", "allele B" = "#313695")) +
              scale_alpha_manual(name = "allelic imbalance", values = c(1, 0.5), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
              ylab("mutation copy number") +
              # ylim(0, mut_cpn_ylim) +
              xlab(chr) +
              theme_bw() 
            
            p_grey <- ggplot(mutation_df, aes(pos, mut_cpn, fill = phase), color = "black") +
              geom_segment(data = unique(mutation_df[, c("seg_start", "seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_maj, yend = seg_y_pos_maj, color = 'allele A', alpha = allelic_imbalance), linewidth = 2) +
              geom_segment(data = unique(mutation_df[, c("seg_start", "seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_min, yend = seg_y_pos_min, color = 'allele B', alpha = allelic_imbalance),linewidth = 2) +
              geom_point(alpha = 0.7, shape = 21, stroke = NA) +
              scale_fill_manual(name = "mutation phase", breaks = c("A", "B"), labels = c("allele A", "allele B"), values = c("#878787", "#878787")) +
              scale_color_manual(name = "copy number phase", values=c("allele A" = "darkred", "allele B" = "#313695")) +
              scale_alpha_manual(name = "allelic imbalance", values = c(1, 0.5), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
              ylab("mutation copy number") +
              # ylim(0, mut_cpn_ylim) +
              xlab(chr) +
              theme_bw() 
            
            p_blank <- ggplot(mutation_df, aes(pos, mut_cpn, fill = phase), color = "black") +
              geom_segment(data = unique(mutation_df[, c("seg_start", "seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_maj, yend = seg_y_pos_maj, color = 'allele A', alpha = allelic_imbalance), linewidth = 2) +
              geom_segment(data = unique(mutation_df[, c("seg_start", "seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = seg_start, xend = seg_end, y = seg_y_pos_min, yend = seg_y_pos_min, color = 'allele B', alpha = allelic_imbalance),linewidth = 2) +
              scale_color_manual(name = "copy number phase", values=c("allele A" = "darkred", "allele B" = "#313695")) +
              scale_alpha_manual(name = "allelic imbalance", values = c(1, 0.5), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
              ylab("mutation copy number") +
              # ylim(0, mut_cpn_ylim) +
              xlab(chr) +
              theme_bw() 
            
            
            if(!dir.exists(paste0(output_path, tumour, "/", region))){
              dir.create(paste0(output_path, tumour, "/", region), recursive = T)
            }
            
            pdf(paste0(output_path, tumour, "/", region, "/", region, "_", chr, "_phased_mutation_copynumber_red_blue.pdf"), width = 9, height = 3.5)
            print(p)
            dev.off()
            
            pdf(paste0(output_path, tumour, "/", region, "/", region, "_", chr, "_phased_mutation_copynumber_red.pdf"), width = 9, height = 3.5)
            print(p_red)
            dev.off()
            
            pdf(paste0(output_path, tumour, "/", region, "/", region, "_", chr, "_phased_mutation_copynumber_blue.pdf"), width = 9, height = 3.5)
            print(p_blue)
            dev.off()
            
            pdf(paste0(output_path, tumour, "/", region, "/", region, "_", chr, "_phased_mutation_copynumber_grey.pdf"), width = 9, height = 3.5)
            print(p_grey)
            dev.off()
            
            pdf(paste0(output_path, tumour, "/", region, "/", region, "_", chr, "_phased_mutation_copynumber_CN.pdf"), width = 9, height = 3.5)
            print(p_blank)
            dev.off()
            
          }
        }
        
        chr_df <- do.call(rbind, chr_df)
        write.table(chr_df, paste0(output_path, tumour, "/", region, "/", region, "_", "segment_mutation_phase_info.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
        
        # add some metrics to this data frame for every segment etc
        metrics_chr_df <- chr_df
        metrics_chr_df$segID <- paste0(metrics_chr_df$chr, "_", metrics_chr_df$seg_start, "_", metrics_chr_df$seg_end)
        
        # get overall number of mutations on each segment
        # get number of A allele mutations on each segment
        
        metrics_chr_df <- metrics_chr_df %>% 
                          group_by(segID) %>% 
                          mutate(number_mutations_on_seg = n(),
                                 number_A_mutations_on_seg = sum(phase == "A"),
                                 number_B_mutations_on_seg = sum(phase == "B"))
                          
        #need to determine if mutation is gained or not
        # is the multiplicity closer to the major or to the minor allele?
        metrics_chr_df$multiplicity_maj_distance <- abs(metrics_chr_df$maj_cn - metrics_chr_df$mut_cpn)
        metrics_chr_df$multiplicity_min_distance <- abs(metrics_chr_df$min_cn - metrics_chr_df$mut_cpn)
        metrics_chr_df[metrics_chr_df$multiplicity_maj_distance < metrics_chr_df$multiplicity_min_distance, "copy_number_allele"] <- "major"
        metrics_chr_df[metrics_chr_df$multiplicity_min_distance < metrics_chr_df$multiplicity_maj_distance, "copy_number_allele"] <- "minor"
        metrics_chr_df[is.na(metrics_chr_df$copy_number_allele), "copy_number_allele"] <- "unknown"
        
        # now for A and B mutations count how many are on major and on minor allele
        metrics_chr_df <- metrics_chr_df %>% 
                            group_by(segID) %>% 
                            mutate(number_A_mutations_on_major = sum(phase == "A" & copy_number_allele == "major"),
                                   number_B_mutations_on_major = sum(phase == "B" & copy_number_allele == "major"),
                                   number_A_mutations_on_minor = sum(phase == "A" & copy_number_allele == "minor"),
                                   number_B_mutations_on_minor = sum(phase == "B" & copy_number_allele == "minor"))
                            
        
        
        metrics_seg_df <- metrics_chr_df[, c("sample", "segID", "chr", "seg_start", "seg_end", "maj_cn", "min_cn", "allelic_imbalance", "number_mutations_on_seg", "number_A_mutations_on_seg", 
                                             "number_B_mutations_on_seg", "number_A_mutations_on_major", "number_B_mutations_on_major", "number_A_mutations_on_minor", "number_B_mutations_on_minor")]
        
        metrics_seg_df <- unique(metrics_seg_df)
        
        

        # also plot it across the whole genome
        chr_df <- genome.continous.fun(chr_df, genome)
        chr_df$allelic_imbalance <- factor(chr_df$allelic_imbalance, levels = c(TRUE, FALSE))
        
        mut_cpn_ylim_genome <- max(chr_df$mut_cpn)
        mut_cpn_ylim_genome_min <- min(chr_df$mut_cpn)
        
        xaxis_label <- data.frame(chr_midpoint)
        xaxis_label$label <- rownames(xaxis_label)

        p_genome <- ggplot(chr_df, aes(adjust_pos, mut_cpn, fill = phase), color = "black") +
          geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
          geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
          geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_maj, yend = seg_y_pos_maj, color = 'allele A', alpha = allelic_imbalance), linewidth = 2) +
          geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_min, yend = seg_y_pos_min, color = 'allele B', alpha = allelic_imbalance),linewidth = 2) +
          scale_fill_manual(name = "mutation phase", breaks = c("A", "B"), labels = c("allele A", "allele B"), values = c("#d73027", "#4575b4")) +
          scale_color_manual(name = "copy number phase", values=c("allele A" = "darkred", "allele B" = "#313695")) +
          scale_alpha_manual(name = "allelic imbalance", values = c(1, 0.5), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
          scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
          ylab("mutation copy number") +
          xlab("") +
          ylim(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank())
        
        p_genome_red <- ggplot(chr_df[chr_df$phase == "A", ], aes(adjust_pos, mut_cpn, fill = phase), color = "black") +
          geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
          geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
          geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_maj, yend = seg_y_pos_maj, color = 'allele A', alpha = allelic_imbalance), linewidth = 2) +
          geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_min, yend = seg_y_pos_min, color = 'allele B', alpha = allelic_imbalance),linewidth = 2) +
          scale_fill_manual(name = "mutation phase", breaks = c("A", "B"), labels = c("allele A", "allele B"), values = c("#d73027", "#4575b4")) +
          scale_color_manual(name = "copy number phase", values=c("allele A" = "darkred", "allele B" = "#313695")) +
          scale_alpha_manual(name = "allelic imbalance", values = c(1, 0.5), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
          scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
          ylab("mutation copy number") +
          xlab("") +
          ylim(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank())
        
        p_genome_blue <- ggplot(chr_df[chr_df$phase == "B", ], aes(adjust_pos, mut_cpn, fill = phase), color = "black") +
          geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
          geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
          geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_maj, yend = seg_y_pos_maj, color = 'allele A', alpha = allelic_imbalance), linewidth = 2) +
          geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_min, yend = seg_y_pos_min, color = 'allele B', alpha = allelic_imbalance),linewidth = 2) +
          scale_fill_manual(name = "mutation phase", breaks = c("A", "B"), labels = c("allele A", "allele B"), values = c("#d73027", "#4575b4")) +
          scale_color_manual(name = "copy number phase", values=c("allele A" = "darkred", "allele B" = "#313695")) +
          scale_alpha_manual(name = "allelic imbalance", values = c(1, 0.5), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
          scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
          ylab("mutation copy number") +
          xlab("") +
          ylim(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank())
        
        p_genome_grey <- ggplot(chr_df, aes(adjust_pos, mut_cpn, fill = phase), color = "black") +
          geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
          geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
          geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_maj, yend = seg_y_pos_maj, color = 'allele A', alpha = allelic_imbalance), linewidth = 2) +
          geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_min, yend = seg_y_pos_min, color = 'allele B', alpha = allelic_imbalance),linewidth = 2) +
          scale_fill_manual(name = "mutation phase", breaks = c("A", "B"), labels = c("allele A", "allele B"), values = c("#878787", "#878787")) +
          scale_color_manual(name = "copy number phase", values=c("allele A" = "darkred", "allele B" = "#313695")) +
          scale_alpha_manual(name = "allelic imbalance", values = c(1, 0.5), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
          scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
          ylab("mutation copy number") +
          xlab("") +
          ylim(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank())
        
        p_genome_blank <- ggplot(chr_df, aes(adjust_pos, mut_cpn, fill = phase), color = "black") +
          geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
          geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_maj, yend = seg_y_pos_maj, color = 'allele A', alpha = allelic_imbalance), linewidth = 2) +
          geom_segment(data = unique(chr_df[, c("adjust_seg_start", "adjust_seg_end", "seg_y_pos_min", "seg_y_pos_maj", "allelic_imbalance", "phase")]), aes(x = adjust_seg_start, xend = adjust_seg_end, y = seg_y_pos_min, yend = seg_y_pos_min, color = 'allele B', alpha = allelic_imbalance),linewidth = 2) +
          scale_fill_manual(name = "mutation phase", breaks = c("A", "B"), labels = c("allele A", "allele B"), values = c("#d73027", "#4575b4")) +
          scale_color_manual(name = "copy number phase", values=c("allele A" = "darkred", "allele B" = "#313695")) +
          scale_alpha_manual(name = "allelic imbalance", values = c(1, 0.5), breaks = c("TRUE", "FALSE"), label = c("TRUE", "FALSE")) +
          scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
          ylab("copy number") +
          xlab("") +
          ylim(mut_cpn_ylim_genome_min, mut_cpn_ylim_genome) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank())
        
        if(!dir.exists(paste0(output_path, tumour, "/", region))){
          dir.create(paste0(output_path, tumour, "/", region), recursive = T)
        }
        
        pdf(paste0(output_path, tumour, "/", region, "/", region, "_genome_phased_mutation_copynumber_red_blue.pdf"), width = 18, height = 3.5)
        print(p_genome)
        dev.off()
        
        pdf(paste0(output_path, tumour, "/", region, "/", region, "_genome_phased_mutation_copynumber_red.pdf"), width = 18, height = 3.5)
        print(p_genome_red)
        dev.off()
        
        pdf(paste0(output_path, tumour, "/", region, "/", region, "_genome_phased_mutation_copynumber_blue.pdf"), width = 18, height = 3.5)
        print(p_genome_blue)
        dev.off()
        
        pdf(paste0(output_path, tumour, "/", region, "/", region, "_genome_phased_mutation_copynumber_grey.pdf"), width = 18, height = 3.5)
        print(p_genome_grey)
        dev.off()
        
        pdf(paste0(output_path, tumour, "/", region, "/", region, "_genome_phased_mutation_copynumber_CN.pdf"), width = 18, height = 3.5)
        print(p_genome_blank)
        dev.off()
        
      }
    }
  }
}


