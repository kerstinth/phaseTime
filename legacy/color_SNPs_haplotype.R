# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # #                       color SNPs by haplotype                       # # #   
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###############################################################################
##################################################################    libraries
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)

source("haplo_timing_functions.R")

###############################################################################
##################################################################    file path

genome                     <- BSgenome.Hsapiens.UCSC.hg38

###############################################################################
##################################################################    functions


###############################################################################
##################################################################         MAIN

for(region in all_regions){
 
  print(region)
  
  tumour = strsplit(region, "_")[[1]][1]

  haplo_paths <- list.files(path = paste0(battenberg_base_path, tumour, "/Battenberg/aggregated/"), pattern = "_multisample_haplotypes_", full.names = T)
  haplo_paths <- haplo_paths[grep("RDS", haplo_paths)]
  
  BAF_data <- lapply(haplo_paths, function(x){
    gr <- readRDS(x)
  })

  BAF_data <- do.call(c, BAF_data)
  colnames(mcols(BAF_data)) <- sub("--..........", "", colnames(mcols(BAF_data)))
  
  col_Major <- paste0(region, "_Major")
  col_BAF   <- paste0(region, "_BAF")
  col_PS    <-  paste0(region, "_PS")
  
  metadata_columns <- mcols(BAF_data)
  column_names <- colnames(metadata_columns)

  if(any(grepl(region, column_names))){

    BAF_region <- BAF_data[, c("REF", "ALT", paste0(col_Major), paste0(col_BAF), paste0(col_PS), "multisample_haplo")]
    BAF_region$REF <- as.character(BAF_region$REF)
    BAF_region$ALT <- as.character(unlist(BAF_region$ALT))
    
    BAF_region <- data.frame(BAF_region)
    BAF_region$width <- NULL
    BAF_region$strand <- NULL
    BAF_region$end <- NULL
    
    colnames(BAF_region) <- c("chr", "pos", "REF", "ALT", "Major", "BAF", "PS", "multisample_haplo")
    
    #cerate a collumn making a phase with BAF phase only
    BAF_region[BAF_region$BAF < 0.5 & BAF_region$Major %in% c("0|1", "1|0"), "bb_phase"] <- "B"
    BAF_region[BAF_region$BAF >= 0.5 & BAF_region$Major %in% c("0|1", "1|0"), "bb_phase"] <- "A"
    BAF_region[is.na(BAF_region$bb_phase), "bb_phase"] <- "not_phasable"

    # make a plot with the BAF phase and with the multisample haplo phase
    
    #make genome continous
    BAF_region$chr <- paste0("chr", BAF_region$chr)
    BAF_region <- genome.continous.fun(BAF_region, BSgenome.Hsapiens.UCSC.hg38)
    
    chr_sizes <- seqlengths(genome)[paste0("chr", 1:22)]
    chr_lines <- cumsum(as.numeric(as.character(chr_sizes)))
    chr_midpoint <- chr_lines - chr_sizes / 2
    
    xaxis_label <- data.frame(chr_midpoint)
    xaxis_label$label <- rownames(xaxis_label)
    
    BAF_region$bb_phase <- factor(BAF_region$bb_phase, levels = c("A", "B", "not_phasable"))
    
    BAF_region[BAF_region$multisample_haplo == 0, "multisample_haplo_phase"] <- "A"
    BAF_region[BAF_region$multisample_haplo == 1, "multisample_haplo_phase"] <- "B"
    
    BAF_region$multisample_haplo_phase <- factor(BAF_region$multisample_haplo_phase, levels = c("A", "B"))
    
    p_genome_bb <- ggplot(BAF_region, aes(adjust_pos, BAF, fill = bb_phase), color = "black") +
                        geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
                        geom_point(shape = 21, stroke = NA, size = 1) +
                        scale_fill_manual(name = "SNP single region haplotype", breaks = c("A", "B", "not_phasable"), labels = c("allele A", "allele B", "not phasable"), values = c("#d73027", "#4575b4", "#bababa")) +
                        scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
                        ylab("B-allele frequency (BAF)") +
                        xlab("") +
                        ylim(c(0, 1)) +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
                              panel.grid.major.x = element_blank(),
                              panel.grid.minor.x = element_blank())
    
    p_genome_multi <- ggplot(BAF_region, aes(adjust_pos, BAF, fill = multisample_haplo_phase), color = "black") +
                        geom_vline(xintercept = chr_lines, linetype = 2, color = "#bababa") +
                        geom_point(shape = 21, stroke = NA, size = 1) +
                        scale_fill_manual(name = "SNP multi region haplotype", breaks = c("A", "B"), labels = c("allele A", "allele B"), values = c("#d73027", "#4575b4")) +
                        scale_x_continuous(breaks = xaxis_label$chr_midpoint, labels = xaxis_label$label, expand = c(0,0)) +
                        ylab("B-allele frequency (BAF)") +
                        xlab("") +
                        ylim(c(0, 1)) +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
                              panel.grid.major.x = element_blank(),
                              panel.grid.minor.x = element_blank())
    
    #make these into one plot
    p <- p_genome_bb / p_genome_multi + plot_annotation(title = region)
    
    
    if(!dir.exists(paste0(output_path, tumour, "/", region, "/plots/phased_SNPs"))){
      dir.create(paste0(output_path, tumour, "/", region, "/plots/phased_SNPs"), recursive = T)
    }

    write.table(BAF_region, paste0(out_temp, region, "_SNPs_phased_single_multi_haplotype.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
    
    # png(paste0(output_path, tumour, "/", region, "/plots/phased_SNPs/", region, "_SNPs_phased_single_multi_haplotype.png"), width = 24, height = 12, units = "cm", res = 300)
    # print(p)
    # dev.off()
    
    
    # make one plot for each chromosome
    
    # for(chr in unique(BAF_region$chr)){
      
    #   print(chr)
      
    #   BAF_region_chr <- BAF_region[BAF_region$chr == chr, ]
      
    #   BAF_region_chr$bb_phase <- factor(BAF_region_chr$bb_phase, levels = c("A", "B", "not_phasable"))
    #   BAF_region_chr$multisample_haplo_phase <- factor(BAF_region_chr$multisample_haplo_phase, levels = c("A", "B"))
      
    #   p_chr_genome_bb <- ggplot(BAF_region_chr, aes(pos, BAF, fill = bb_phase), color = "black") +
    #     geom_point(shape = 21, stroke = NA, size = 1) +
    #     scale_fill_manual(name = "SNP single region haplotype", breaks = c("A", "B", "not_phasable"), labels = c("allele A", "allele B", "not phasable"), values = c("#d73027", "#4575b4", "#bababa")) +
    #     ylab("B-allele frequency (BAF)") +
    #     xlab(chr) +
    #     ylim(c(0, 1)) +
    #     theme_bw() +
    #     theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
      
    #   p_chr_genome_multi <- ggplot(BAF_region_chr, aes(pos, BAF, fill = multisample_haplo_phase), color = "black") +
    #     geom_point(shape = 21, stroke = NA, size = 1) +
    #     scale_fill_manual(name = "SNP multi region haplotype", breaks = c("A", "B"), labels = c("allele A", "allele B"), values = c("#d73027", "#4575b4")) +
    #     ylab("B-allele frequency (BAF)") +
    #     xlab(chr) +
    #     ylim(c(0, 1)) +
    #     theme_bw() +
    #     theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
      
    #   #make these into one plot
    #   p_chr <- p_chr_genome_bb / p_chr_genome_multi + plot_annotation(title = paste0(region, " ", chr))
      
    #   if(!dir.exists(paste0(output_path, tumour, "/", region, "/plots/phased_SNPs/per_chromosome/"))){
    #     dir.create(paste0(output_path, tumour, "/", region, "/plots/phased_SNPs/per_chromosome"), recursive = T)
    #   }
      
    #   png(paste0(output_path, tumour, "/", region, "/plots/phased_SNPs/per_chromosome/", region, "_", chr, "_SNPs_phased_single_multi_haplotype.png"), width = 24, height = 12, units = "cm", res = 300)
    #   print(p_chr)
    #   dev.off()

    # }
  }
}




