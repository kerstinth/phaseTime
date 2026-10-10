################################################################################
####                  simulate phasing SNPs to mutations                    ####  
################################################################################

################################################################################
###################################################################### libraries
library(dplyr)
library(ggplot2)

################################################################################
######################################################################file paths


################################################################################
######################################################################      MAIN

#let's assume a genome has 4.5million SNP (average form 1000 genomes) https://www.nature.com/articles/nature15393
#but acraully we only want to phase heterozygous SNPs which is about 2.25 million: https://web.stanford.edu/group/pritchardlab/HGbook/Release-2023-09/HGBook-2023-09-chapters/HGBook-2023-09-23-ch1.3.pdf
# and let's assume the gemome has 3,054,815,472 nucleotides https://www.science.org/doi/10.1126/science.abj6987

#lets iterate over different length of sequencing reads
read_lenghts <- c(100, 150, 200, 300, 500, 1000, 2000, 5000, 7500, 10000)

#at over different numbers of mutations per genome
mutations <- c(50000)

all_out <- list()
for(read_length in read_lenghts){
  
  print(read_length)
  
  mut_out <- list()
  for(muts in mutations){
    print(muts)
    
    #randomly sample the SNP positions
    SNP_pos <- sample(3054815472, 2250000, replace = FALSE)
    
    #and now randomly sample the mut positions
    mut_pos <- sample(3054815472, muts, replace = FALSE)

    #gotta make sure the muts and the SNP pos are not overlapping
    if(any(mut_pos %in% SNP_pos)){
      
      overlaps <- mut_pos[mut_pos %in% SNP_pos]
      overlap_length <- length(overlaps)

      mut_pos <- mut_pos[!mut_pos %in% SNP_pos]

      new_muts <- sample(3054815472, overlap_length, replace = FALSE)

      mut_pos <- c(mut_pos, new_muts)

       if(any(mut_pos %in% SNP_pos)){

        overlaps <- mut_pos[mut_pos %in% SNP_pos]
        overlap_length <- length(overlaps)

        mut_pos <- mut_pos[!mut_pos %in% SNP_pos]

        new_muts <- sample(3054815472, overlap_length, replace = FALSE)

        mut_pos <- c(mut_pos, new_muts)

          if(any(mut_pos %in% SNP_pos)){

            print("still overlaps in 3rd round")

          }

       }

    }
    
    #now bin the genome into read length and check if any muts and SNPs are ever in the same bin
    bins <- data.frame(bin_start = seq(0, 3054815472, by = read_length))
    bins$bin_end <- bins$bin_start + (read_length-1)
    bins$index <- 1:nrow(bins)

    #gotta bin the muts
    mut_bins <- data.frame(mutations = mut_pos)
    mut_bins <- mut_bins %>% mutate(bin_index = findInterval(mutations, bins$bin_start))
    mut_bins <- left_join(mut_bins, bins, by = c("bin_index" = "index"))
    mut_bins$binID <- paste0(mut_bins$bin_start, ":", mut_bins$bin_end)

    #and bin the SNPs
    snp_bins <- data.frame(snps = SNP_pos)
    snp_bins <- snp_bins %>% mutate(bin_index = findInterval(snps, bins$bin_start))
    snp_bins <- left_join(snp_bins, bins, by = c("bin_index" = "index"))
    snp_bins$binID <- paste0(snp_bins$bin_start, ":", snp_bins$bin_end)

    # now, how many muts are in the same bin as a SNP
    muts_in_SNP_bin <- length(mut_bins$binID[mut_bins$binID %in% snp_bins$binID])
    
    df <- data.frame(read_length = read_length,
                    number_mutations = muts,
                    number_phased_muts = muts_in_SNP_bin)


    mut_out[[as.character(muts)]] <- df
  }
  mut_out <- do.call(rbind, mut_out)
  
  all_out[[as.character(read_length)]] <- mut_out
}

all_out <- do.call(rbind, all_out)

#calculate propoprtion
all_out$prop <- (as.numeric(as.character(all_out$number_phased_muts)) / as.numeric(as.character(all_out$number_mutations)))*100

# make a plot
all_out$number_mutations <- as.character(all_out$number_mutations)
all_out$read_length <- as.character(all_out$read_length)

all_out$number_mutations <- factor(all_out$number_mutations, levels = unique(all_out$number_mutations))
all_out$read_length <- factor(all_out$read_length, levels = unique(all_out$read_length))


p <- ggplot(all_out, aes(read_length, prop)) +
      geom_bar(stat = "identity", fill = "darkblue") +
      # scale_fill_gradient(name = "# phased mutations", low="white", high="#253494", trans = "log10") +
      ylab("% mutations phased to heterozygous SNPs") +
      xlab("read length") +
      theme_bw() 



pdf(paste0(output_path, "simulated_SNP_mutation_phasing.pdf"), width = 4, height = 4)
p
dev.off()
