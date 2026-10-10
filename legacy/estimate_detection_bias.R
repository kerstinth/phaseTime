# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #               long-read mutatino timing detection bias              # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###############################################################################
#################################################################     libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

###############################################################################
#################################################################    file paths

###############################################################################
#################################################################          MAIN

# mutation_data_all <- fread(mutation_data_base_path, head = T, sep = "\t", data.table = FALSE)
# 
# muts <- mutation_data_all[, c("sample", "chr", "pos", "ref", "alt", "varcount_tumour", "refcount_tumour", "vaf", "maj_cn", "min_cn", "tot_cn", "mut_cpn"), ]
# muts$cov <- muts$varcount_tumour + muts$refcount_tumour

# min_varcount <- min(muts$varcount_tumour, na.rm = T)
# min_vaf <- min(muts$vaf, na.rm = T)
# max_cov <- max(muts$cov, na.rm = T)
# min_cov <- min(muts$cov, na.rm = T)

min_vaf <- 0.005
max_cov <- 5000
min_cov <- 1

# how many reads do I need to at least get a vaf of target_vaf up until a coverage of max_cov
#do this for several vafs
# target_vaf = 0.005

# df_reads <- data.frame(coverage = seq(min_cov, max_cov, 50))
# df_vafs <- data.frame(vaf = seq(0.005, 1, 0.005))
# 
# df_reads <- cross_join(df_reads, df_vafs)
# df_reads$mutated_reads <- df_reads$coverage * df_reads$vaf


# but we need to take mutation multiplicity into account
# and we need to do this for various puities

purities <- c(seq(0, 1, by = 0.05))
purities <- purities[purities != 0.0]

# and we want to do this up until the maximum copy number that we detect
# max_maj_cn <- max(muts$maj_cn, na.rm = T)
max_maj_cn <- 65
normal_cn <- 1

#let's iterate over several VAFs

multiplicities <- list()
for(p in purities){
  
  print(p)
  
  vaf_list <- list()
  for(v in seq(0.005, 1, 0.005)){
     print(v)
    
    cn_list <- list()
    for(cn in 1:max_maj_cn){
      
      print(cn)
      multi <- (v * 1/p) * ((p * cn) + (normal_cn * (1 -p)))
      
      df_out <- data.frame(tumour_copynumber = cn,
                           purity = p,
                           normal_cn = normal_cn,
                           multiplicity = multi,
                           vaf = v)
      cn_list[[cn]] <- df_out
      
    }
    
    cn_list <- do.call(rbind, cn_list)
    vaf_list[[as.character(v)]] <- cn_list
    
  }
  vaf_list <- do.call(rbind, vaf_list)
  
  multiplicities[[as.character(p)]] <- vaf_list

}
    
multiplicities <- do.call(rbind, multiplicities)  

# multiplicities_cov <- left_join(multiplicities, df_reads)


##############  let's try to plot this

# first of all, show that in regions with really high copy number we have very high coverage
# cn_cov <- muts[, c("maj_cn", "cov")]
# cn_cov <- na.omit(cn_cov)
# cn_cov <- unique(cn_cov)
# 
# p1 <- ggscatter(cn_cov, x = "maj_cn", y = "cov",
#                 add = "reg.line",
#                 conf.int = TRUE,
#                 add.params = list(color = "blue",
#                                   fill = "lightgray")) +
#           stat_cor(method = "pearson", label.x = 0, label.y = 5000) +
#           xlab("major copynumber") +
#           ylab("coverage")
# 
# pdf(paste0(output_path, "cnmaj_cov_cor.pdf"), width = 4, height = 4)
# p1
# dev.off()


# for each purity we want to calculate the cut off copy number value where we cannot pick up a mutation of one copy with a 5% chance anymore

# at a purity of x which is the copy number value where I can still get a multiplicity of 1 at a vaf of 0.05


# df <- multiplicities[multiplicities$purity == 0.4, ]
df <- multiplicities
df <- df[df$multiplicity >= 1, ]
df <- df %>% group_by(purity, tumour_copynumber) %>% slice_min(multiplicity)

p <- ggplot(df, aes(tumour_copynumber, vaf)) +
        geom_point() +
        geom_hline(yintercept = 0.05, colour = "blue") +
        facet_wrap(~purity) +
        theme_bw()


cut_offs <- df[df$vaf >= 0.05, ]
cut_offs <- cut_offs %>% group_by(purity) %>% slice_min(vaf)
cut_offs <- cut_offs %>% group_by(purity) %>% slice_max(tumour_copynumber)


