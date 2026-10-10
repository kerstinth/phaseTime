# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # #                simulations for copy number timing                   # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###############################################################################
#################################################################### libraries

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggnewscale)
library(patchwork)
library(ggborderline)

###############################################################################
#################################################################### file path

###############################################################################
####################################################################  functions

# function to get colors
get.colors <- function(max_cn = max(mutation_df$mut_cpn_round), opacity = 255){
  cols <- paste(brewer.pal(max_cn, name = "Paired"), sep = "")
  nopaq <- round(opacity / ceiling(max_cn / 12))
  opacity_val <- c(sapply(1:ceiling(max_cn / 12), function(x) rep(opacity - nopaq * (x - 1), 12)))[1:max_cn]
  cols.opac <- paste(cols, as.hexmode(opacity_val), sep = "")
}

#function to simulate mutations per cell division of a copy number segment with copy number gains at certain time points
division_simulator <- function(timepoints = c(1000, 2000, 5000, 10000), segment_size = 5000000, cn_steps = c(1, 2, 3, 5, 5), mutation_rate_bp = 0.0000000133){
  
  df_out <- data.frame(interval_index = NA,
                       mutation_rate_bp = NA,
                       segment_size = NA,
                       copy_number_before_gain = NA,
                       copy_number_after_gain = NA,
                       number_cell_divisions_since_last_gain = NA,
                       number_muts_accumulated = NA)
  
  for(i in 1:length(timepoints)){
    
    number_divisions_since_last_gain  <- ifelse(i == 1, timepoints[i], timepoints[i] - timepoints[i-1]) 
    
    # our segment size grows with copy number
    segment_size_cn_adjusted <- ifelse(i == 1, segment_size, segment_size * cn_steps[i])
    
    #whats the mutation rate now
    total_mut_rate <- segment_size_cn_adjusted*mutation_rate_bp
    
    #now let's sample some mutations
    mutations_accumulated <- rpois(n = number_divisions_since_last_gain, lambda = total_mut_rate)
    interval_muts  <- sum(mutations_accumulated) #add these up
    
    # store results
    copy_number_before_gain <- cn_steps[i]
    copy_number_after_gain  <- cn_steps[i + 1]
    
    number_muts_cn_adjusted <- interval_muts / copy_number_before_gain
    
    df_out[i, "interval_index"] <- i
    df_out[i, "mutation_rate_bp"] <- mutation_rate_bp
    df_out[i, "copy_number_before_gain"] <- copy_number_before_gain
    df_out[i, "copy_number_after_gain"] <- copy_number_after_gain
    df_out[i, "segment_size"] <- segment_size
    df_out[i, "present_dna_length_before_gain"] <- segment_size_cn_adjusted
    df_out[i, "number_cell_divisions_since_last_gain"] <- number_divisions_since_last_gain
    df_out[i, "number_muts_accumulated"] <- interval_muts
    df_out[i, "adjustment_factor"] <- copy_number_before_gain
    df_out[i, "number_muts_cn_adjusted"] <- number_muts_cn_adjusted
    
  }
  return(df_out)
}

# function to calculate and plot copy number gain timing
calc_gain_timing <- function(mutation_count_df = df, events = c(1:(length(c(1000, 2000, 5000, 10000))- 1)), output_path = output_path, plot_name = "sim_gaintiming_1.pdf", plot = FALSE, write_results = FALSE){
  
  df_out <- data.frame(gain_index = NA,
                       mutation_rate_bp = NA,
                       segment_size = NA,
                       copy_number_before_gain = NA,
                       copy_number_after_gain = NA,
                       number_cell_divisions_since_last_gain = NA,
                       number_early_muts_cn_adjusted = NA,
                       number_all_muts_cn_adjusted = NA,
                       gain_timing = NA,
                       loCI = NA,
                       hiCI = NA,
                       theoretical_gain_timing = NA)
  
  for(event in events){
    
    gain_index                            <- event
    mutation_rate_bp                      <- unique(mutation_count_df$mutation_rate_bp)
    segment_size                          <- unique(mutation_count_df$segment_size)
    copy_number_before_gain               <- mutation_count_df[event, "copy_number_before_gain"]
    copy_number_after_gain                <- mutation_count_df[event, "copy_number_after_gain"]
    number_cell_divisions_since_last_gain <- mutation_count_df[event, "number_cell_divisions_since_last_gain"]
    
    number_early_muts                     <- sum(mutation_count_df[1:event, "number_muts_cn_adjusted"])
    number_all_muts                       <- sum(mutation_count_df[, "number_muts_cn_adjusted"])
    
    gain_timing                           <- number_early_muts / number_all_muts
    error                                 <- qnorm(0.975) * sqrt(gain_timing * (1 - gain_timing) / number_all_muts)
    loCI                                  <- gain_timing - error
    hiCI                                  <- gain_timing + error
    theoretical_gain_timing               <- ifelse(event == 1, number_cell_divisions_since_last_gain / sum(mutation_count_df$number_cell_divisions_since_last_gain), (number_cell_divisions_since_last_gain / sum(mutation_count_df$number_cell_divisions_since_last_gain)) + df_out[event-1, "theoretical_gain_timing"])
    
    df_out[event, "gain_index"]                             <- gain_index
    df_out[event, "mutation_rate_bp"]                       <- mutation_rate_bp
    df_out[event, "segment_size"]                           <- segment_size
    df_out[event, "copy_number_before_gain"]                <- copy_number_before_gain
    df_out[event, "copy_number_after_gain"]                 <- copy_number_after_gain
    df_out[event, "number_cell_divisions_since_last_gain"]  <- number_cell_divisions_since_last_gain
    df_out[event, "number_early_muts_cn_adjusted"]          <- number_early_muts
    df_out[event, "number_all_muts_cn_adjusted"]            <- number_all_muts
    df_out[event, "gain_timing"]                            <- gain_timing
    df_out[event, "loCI"]                                   <- loCI
    df_out[event, "hiCI"]                                   <- hiCI
    df_out[event, "theoretical_gain_timing"]                <- theoretical_gain_timing
    
  }
  

  #also create a plot to show the result
  df_out[df_out$loCI < 0, "loCI"] <- 0
  df_out[df_out$loCI < 0, "hiCI"] <- 0
  df_out[df_out$loCI > 1, "loCI"] <- 1
  df_out[df_out$loCI > 1, "hiCI"] <- 1
  
  colors <- get.colors(max_cn = max(df_out$gain_index))
  
  df_out$seg_start  <- 1
  df_out$seg_end    <- df_out$segment_size
  df_out$gain_index <- factor(df_out$gain_index, levels = sort(df_out$gain_index))
  
  
  #make a dummy dataframe for the theoretical timing
  theo_time <- df_out[, c("gain_index", "theoretical_gain_timing")]
  theo_time <- theo_time[rep(seq_len(nrow(theo_time)), each = 2), ]
  theo_time$pos <- rep(c(1, segment_size), nrow(theo_time)/2)
  
  p <- ggplot() +
    geom_rect(data = df_out, aes(xmin = seg_start, xmax = seg_end, ymin = loCI, ymax = hiCI, fill = gain_index), color = NA, alpha = 0.15) +
    geom_segment(data = df_out,
                 aes(x = seg_start, xend = seg_end, y = gain_timing, yend = gain_timing, color = gain_index), linewidth = 3) +
    scale_fill_manual(name = "gain order", values = colors, guide = "none") +
    scale_color_manual(name = "calculated gain timing", values = colors) +
    new_scale_color() +
    geom_borderline(data = theo_time,
                    aes(x = pos, y = theoretical_gain_timing, group = gain_index, color = gain_index), bordercolour = "#525252", linewidth = 1) +
    scale_color_manual(name = "ground truth gain timing", values = colors) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.00\n(early)", "0.25", "0.50", "0.75", "1.00\n(late)"), limits = c(0, 1)) + 
    ylab("gain timing") +
    xlab("copy number segment") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())
  
  # also make a plot what the mutations on a copy number profile would look like for visualization purposes
  df_muts <- data.frame(mut_cpn = mutation_count_df$copy_number_before_gain,
                        number_muts_cn_adjusted = round(mutation_count_df$number_muts_cn_adjusted))
  
  #let's sample our mutations
  mutations <- list()
  for(i in 1:nrow(df_muts)){
    
    num <- df_muts[i, "number_muts_cn_adjusted"]
    pos <- sample(segment_size, num, replace = FALSE)
    
    mut_cpn <- df_muts[i, "mut_cpn"]

    multi  <- rep(mut_cpn, num)
    multi  <- jitter(multi)
    
    if(length(multi) == 0){
      multi <- NA
      pos <- NA
    }
    
    
    mut_out_df <- data.frame(mut_cpn = df_muts[i, "mut_cpn"],
                             number_muts_cn_adjusted = num,
                             pos = pos,
                             multiplicity = multi)
    mutations[[i]] <- mut_out_df
    
  }
  
  mutations <- do.call(rbind, mutations)
  
  #create plot
  color_index <- data.frame(mut_cpn = sort(unique(mutations$mut_cpn)),
                            color_index = rep(c("one", "two"), length.out = length(unique(mutations$mut_cpn))))
  
  mutations <- suppressMessages(left_join(mutations, color_index))
  mutations$color_index <- factor(mutations$color_index, levels = c("one", "two"))
  
  p2 <- ggplot(mutations, aes(pos, multiplicity, fill = color_index), color = "black") +
            geom_point(alpha = 0.7, shape = 21, stroke = NA, size = 1) +
            scale_fill_manual(name = "", breaks = c("one", "two"), labels = c("one", "two"), values = c("#d8b365", "#5ab4ac"), guide = "none") +
            ylab("mutation copy number") +
            xlab("copy number segment") +
            scale_y_continuous(breaks = seq(0, round(max(mutations$multiplicity, na.rm = T)), by = 2), limits = c(0, max(mutations$multiplicity, na.rm = T))) +
            theme_bw() +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank())
  
  
  #join the two plots
  p_all <- p2 + p + plot_annotation(title = paste0("mutation rate: ", unique(mutation_count_df$mutation_rate_bp)))
  
  if(plot){
    pdf(paste0(output_path, plot_name), width = 10, height = 4)
    print(p_all)
    dev.off()
  }

  if(write_results){
    write.table(df_out, paste0(output_path, sub(".pdf", ".txt", plot_name)), row.names = F, col.names = T, sep = "\t", quote = F)
  }

  return(df_out)
  
}

###############################################################################
####################################################################       MAIN

#### run

mutation_rate_bp    <- 0.0000000133
cn_steps            <- c(1, 2, 3, 5, 5)
timepoints          <- c(1000, 2000, 5000, 10000)
segment_size        <- 5000000

df <- division_simulator(timepoints = timepoints, segment_size = segment_size, cn_steps = cn_steps, mutation_rate_bp = mutation_rate_bp)
timing <- calc_gain_timing(mutation_count_df = df, events = c(1:(length(timepoints)- 1)), output_path = output_path, plot_name = "sim_gaintiming_1.pdf", plot = TRUE, write_results = F)

mutation_rate_bp    <- 0.0000000133
cn_steps            <- c(1, 3, 10, 10)
timepoints          <- c(7000, 8000, 10000)
segment_size        <- 5000000

df <- division_simulator(timepoints = timepoints, segment_size = segment_size, cn_steps = cn_steps, mutation_rate_bp = mutation_rate_bp)
timing <- calc_gain_timing(mutation_count_df = df, events = c(1:(length(timepoints)- 1)), output_path = output_path, plot_name = "sim_gaintiming_2.pdf", plot = TRUE, write_results = F)

mutation_rate_bp    <- 0.0000000133
cn_steps            <- c(1, 2, 3, 5, 8, 12)
timepoints          <- c(1000, 2000, 2300, 5000, 10000)
segment_size        <- 5000000

df <- division_simulator(timepoints = timepoints, segment_size = segment_size, cn_steps = cn_steps, mutation_rate_bp = mutation_rate_bp)
timing <- calc_gain_timing(mutation_count_df = df, events = c(1:(length(timepoints)- 1)), output_path = output_path, plot_name = "sim_gaintiming_3.pdf", plot = TRUE, write_results = F)

#iterate over all sorts of scenarios to say soomething about how accurate this is

all_res <- list()
for(i in 1:1000){
  
  mutation_rate_bp    <- 0.0000000133
  segment_size        <- 5000000
  
  #randomly determine how many CN gains there are
  num_gains <- round(runif(1, min = 1, max = 5))
  
  cn_steps            <- sample(2:10, num_gains)
  cn_steps            <- sort(cn_steps)
  cn_steps            <- c(1, cn_steps, max(cn_steps))
  
  timepoints          <- sample(1:10000, (length(cn_steps)-2))
  timepoints          <- sort(timepoints)
  timepoints          <- c(timepoints, 10000)
  
  df <- division_simulator(timepoints = timepoints, segment_size = segment_size, cn_steps = cn_steps, mutation_rate_bp = mutation_rate_bp)
  timing <- calc_gain_timing(mutation_count_df = df, events = c(1:(length(timepoints)- 1)), output_path = output_path, plot_name = "sim_gaintiming_2.pdf", plot = FALSE, write_results = FALSE)
  timing$iteration <- i
  
  all_res[[i]] <- timing
  
}

all_res <- do.call(rbind, all_res)

all_res$diff <- abs(all_res$gain_timing - all_res$theoretical_gain_timing)

p_hist <- ggplot(all_res, aes(x = diff)) + 
            geom_histogram(color = "darkgrey", fill = "darkblue") +
            xlab("absolute difference") +
            theme_bw()

pdf(paste0(output_path, "hist_abs_diff.pdf"), width = 5, height = 5)
p_hist
dev.off()


#run this with decreasing mutation rates
for(mut_rate in c(0.0000000133, 0.00000000133, 0.000000000133)){
  
  mutation_rate_bp    <- mut_rate
  cn_steps            <- c(1, 2, 2)
  timepoints          <- c(5000, 10000)
  segment_size        <- 5000000
  
  df <- division_simulator(timepoints = timepoints, segment_size = segment_size, cn_steps = cn_steps, mutation_rate_bp = mutation_rate_bp)
  timing <- calc_gain_timing(mutation_count_df = df, events = c(1:(length(timepoints)- 1)), output_path = output_path, plot_name = paste0("sim_gaintiming_1_", "mutation_rate_", mut_rate, ".pdf"), plot = TRUE, write_results = TRUE)
  
}

#let's read the stats back in an figure out how many mutations i need to have to be able to work out the gain timing
all_mutrate <- lapply(c(0.0000000133, 0.00000000133, 0.000000000133), function(mut_rate){
  
  read.table(paste0(output_path, "sim_gaintiming_1_", "mutation_rate_", mut_rate, ".txt"), head = T, sep = "\t")
  
})

all_mutrate <- do.call(rbind, all_mutrate)

# indicator if theoretical gain timing is within confidence intervals

all_mutrate[all_mutrate$loCI <= all_mutrate$theoretical_gain_timing & all_mutrate$hiCI >= all_mutrate$theoretical_gain_timing, "fits"] <- TRUE
all_mutrate[is.na(all_mutrate$fits), "fits"] <- FALSE


all_res_mutrate <- list()
for(i in 1:1000){
  
  segment_size        <- 5000000
  
  mut_rate_df <- list()
  for(mutation_rate_bp in c(0.0000000133, 0.00000000133, 0.000000000133)){
    
    #randomly determine how many CN gains there are
    num_gains <- round(runif(1, min = 1, max = 5))
    
    cn_steps            <- sample(2:10, num_gains)
    cn_steps            <- sort(cn_steps)
    cn_steps            <- c(1, cn_steps, max(cn_steps))
    
    timepoints          <- sample(1:10000, (length(cn_steps)-2))
    timepoints          <- sort(timepoints)
    timepoints          <- c(timepoints, 10000)
    
    df <- division_simulator(timepoints = timepoints, segment_size = segment_size, cn_steps = cn_steps, mutation_rate_bp = mutation_rate_bp)
    timing <- calc_gain_timing(mutation_count_df = df, events = c(1:(length(timepoints)- 1)), output_path = output_path, plot_name = "sim_gaintiming_2.pdf", plot = FALSE, write_results = FALSE)
    timing$iteration <- i
    
    mut_rate_df[[as.character(mutation_rate_bp)]] <- timing
  }
  mut_rate_df <- do.call(rbind, mut_rate_df)
  
  all_res_mutrate[[i]] <- mut_rate_df
  
}

all_res_mutrate <- do.call(rbind, all_res_mutrate)
all_res_mutrate$CI_size <- all_res_mutrate$hiCI - all_res_mutrate$loCI
all_res_mutrate$mutation_rate_bp <- factor(all_res_mutrate$mutation_rate_bp, levels = c(0.0000000133, 0.00000000133, 0.000000000133))

p_hist_ci <- ggplot(all_res_mutrate, aes(x = CI_size)) + 
  geom_histogram(color = "darkgrey", fill = "darkblue") +
  xlab("CI width") +
  facet_wrap(~ mutation_rate_bp, nrow = 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))

pdf(paste0(output_path, "hist_CI_width.pdf"), width = 10, height = 5)
p_hist_ci
dev.off()





# calculate some gain timing
# calc_gain_timing <- function(mutation_count_df = df, events = 1:length(c(1000, 2000, 5000, 10000))- 1){
#   
#   df_out <- data.frame(gain_index = NA,
#                        mutation_rate_bp = NA,
#                        segment_size = NA,
#                        copy_number_before_gain = NA,
#                        copy_number_after_gain = NA,
#                        number_cell_divisions_since_last_gain = NA,
#                        number_early_muts = NA,
#                        number_late_muts = NA,
#                        segment_size_adjustment_factor = NA,
#                        gain_timing = NA,
#                        loCI = NA,
#                        hiCI = NA,
#                        theoretical_gain_timing = NA)
#   
#   for(event in events){
#     
#     gain_index                            <- event
#     mutation_rate_bp                      <- unique(mutation_count_df$mutation_rate_bp)
#     segment_size                          <- unique(mutation_count_df$segment_size)
#     copy_number_before_gain               <- mutation_count_df[event, "copy_number_before_gain"]
#     copy_number_after_gain                <- mutation_count_df[event, "copy_number_after_gain"]
#     number_cell_divisions_since_last_gain <- mutation_count_df[event, "number_cell_divisions_since_last_gain"]
#     number_early_muts                     <- mutation_count_df[1, "number_muts_accumulated"]
#     number_late_muts                      <- mutation_count_df[event + 1, "number_muts_accumulated"]
#     segment_size_adjustment_factor        <- mutation_count_df[event, "copy_number_after_gain"]
#     number_late_muts_cn_adjusted          <- number_late_muts / segment_size_adjustment_factor
#     gain_timing                           <- number_early_muts / (number_early_muts + number_late_muts_cn_adjusted)
#     error                                 <- qnorm(0.975)*sqrt(gain_timing*(1-gain_timing)/(number_late_muts_cn_adjusted + number_early_muts))
#     loCI                                  <- gain_timing - error
#     hiCI                                  <- gain_timing + error
#     theoretical_gain_timing               <- ifelse(event == 1, number_cell_divisions_since_last_gain / sum(mutation_count_df$number_cell_divisions_since_last_gain), (number_cell_divisions_since_last_gain / sum(mutation_count_df$number_cell_divisions_since_last_gain)) + df_out[event-1, "theoretical_gain_timing"])
#     
#     df_out[event, "gain_index"]                             <- gain_index
#     df_out[event, "mutation_rate_bp"]                       <- mutation_rate_bp
#     df_out[event, "segment_size"]                           <- segment_size
#     df_out[event, "copy_number_before_gain"]                <- copy_number_before_gain
#     df_out[event, "copy_number_after_gain"]                 <- copy_number_after_gain
#     df_out[event, "number_cell_divisions_since_last_gain"]  <- number_cell_divisions_since_last_gain
#     df_out[event, "number_early_muts"]                      <- number_early_muts
#     df_out[event, "number_late_muts"]                       <- number_late_muts
#     df_out[event, "segment_size_adjustment_factor"]         <- segment_size_adjustment_factor
#     df_out[event, "number_late_muts_cn_adjusted"]           <- number_late_muts_cn_adjusted
#     df_out[event, "gain_timing"]                            <- gain_timing
#     df_out[event, "loCI"]                                   <- loCI
#     df_out[event, "hiCI"]                                   <- hiCI
#     df_out[event, "theoretical_gain_timing"]                <- theoretical_gain_timing
#       
#   }
#   
#   return(df_out)
#   
# }
# 
# calc_gain_timing(mutation_count_df = df, events = 1:length(c(1000, 2000, 5000, 10000))- 1)


# division_simulator_repeats <- function(timepoints = c(1000, 2000, 5000, 10000), segment_size = 5000000, cn_steps = c(1, 2, 3, 5, 5), mutation_rate_bp = 0.0000000133, repeats = 10000){
#  
#   all_interval_muts <- list()
#   for(r in 1:repeats){
#     all_interval_muts[[r]] <- division_simulator(timepoints, segment_size, cn_steps, mutation_rate_bp)
#   } 
#   
#   all_interval_muts$iteration <- r
#   all_interval_muts <- do.call(rbind, all_interval_muts)
#   
#   return(all_interval_muts)
# 
# }

# #function to calculate ratios across all iterations
# all_iterations <- division_simulator_repeats(timepoints = c(1000, 2000, 5000, 10000), segment_size = 5000000, cn_steps = c(1, 2, 3, 5, 5), mutation_rate_bp = 0.0000000133, repeats = 10000)


# calculate gain timing maybe fixed
# calc_gain_timing <- function(mutation_count_df = df, events = c(1:(length(c(1000, 2000, 5000, 10000))- 1))){
#   
#   df_out <- data.frame(gain_index = NA,
#                        mutation_rate_bp = NA,
#                        segment_size = NA,
#                        copy_number_before_gain = NA,
#                        copy_number_after_gain = NA,
#                        number_cell_divisions_since_last_gain = NA,
#                        number_early_muts = NA,
#                        number_late_muts = NA,
#                        segment_size_adjustment_factor = NA,
#                        gain_timing = NA,
#                        loCI = NA,
#                        hiCI = NA,
#                        theoretical_gain_timing = NA)
#   
#   for(event in events){
#     
#     gain_index                            <- event
#     mutation_rate_bp                      <- unique(mutation_count_df$mutation_rate_bp)
#     segment_size                          <- unique(mutation_count_df$segment_size)
#     copy_number_before_gain               <- mutation_count_df[event, "copy_number_before_gain"]
#     copy_number_after_gain                <- mutation_count_df[event, "copy_number_after_gain"]
#     number_cell_divisions_since_last_gain <- mutation_count_df[event, "number_cell_divisions_since_last_gain"]
#     number_early_muts                     <- mutation_count_df[1, "number_muts_accumulated"]
#     
#     number_late_muts                      <- mutation_count_df[c(event + 1 : nrow(mutation_count_df)), "number_muts_accumulated"]
#     number_late_muts                      <- number_late_muts[!is.na(number_late_muts)]
#     number_late_muts_sum                  <- sum(number_late_muts)
#     
#     segment_size_adjustment_factor        <- mutation_count_df[event : (nrow(mutation_count_df) - 1), "copy_number_after_gain"]
#     number_late_muts_cn_adjusted          <- number_late_muts / segment_size_adjustment_factor
#     number_late_muts_cn_adjusted_cumu     <- sum(number_late_muts_cn_adjusted)
#     gain_timing                           <- number_early_muts / (number_early_muts + number_late_muts_cn_adjusted_cumu)
#     error                                 <- qnorm(0.975)*sqrt(gain_timing*(1-gain_timing)/(number_late_muts_cn_adjusted_cumu + number_early_muts))
#     loCI                                  <- gain_timing - error
#     hiCI                                  <- gain_timing + error
#     theoretical_gain_timing               <- ifelse(event == 1, number_cell_divisions_since_last_gain / sum(mutation_count_df$number_cell_divisions_since_last_gain), (number_cell_divisions_since_last_gain / sum(mutation_count_df$number_cell_divisions_since_last_gain)) + df_out[event-1, "theoretical_gain_timing"])
#     
#     df_out[event, "gain_index"]                             <- gain_index
#     df_out[event, "mutation_rate_bp"]                       <- mutation_rate_bp
#     df_out[event, "segment_size"]                           <- segment_size
#     df_out[event, "copy_number_before_gain"]                <- copy_number_before_gain
#     df_out[event, "copy_number_after_gain"]                 <- copy_number_after_gain
#     df_out[event, "number_cell_divisions_since_last_gain"]  <- number_cell_divisions_since_last_gain
#     df_out[event, "number_early_muts"]                      <- number_early_muts
#     df_out[event, "number_late_muts"]                       <- number_late_muts_sum
#     df_out[event, "number_late_muts_cn_adjusted"]           <- number_late_muts_cn_adjusted_cumu
#     df_out[event, "gain_timing"]                            <- gain_timing
#     df_out[event, "loCI"]                                   <- loCI
#     df_out[event, "hiCI"]                                   <- hiCI
#     df_out[event, "theoretical_gain_timing"]                <- theoretical_gain_timing
#     
#   }
#   
#   return(df_out)
#   
# }
# 
# calc_gain_timing(mutation_count_df = df, events = c(1:(length(c(1000, 2000, 5000, 10000))- 1)))









# for(iter in 1:10000){
#   print(iter)
#   
#   df_out <- data.frame(interval_index = NA,
#                        mutation_rate_bp = NA,
#                        segment_size = NA,
#                        copy_number_before_gain = NA,
#                        copy_number_after_gain = NA,
#                        number_cell_divisions_since_last_gain = NA,
#                        number_muts_accumulated = NA)
#   
#   for(i in 1:length(timepoints)){
#     
#     print(i)
#     
#     number_divisions_since_last_gain  <- ifelse(i == 1, timepoints[i], timepoints[i] - timepoints[i-1]) 
#     
#     # our segment size grows with copy number
#     segment_size_cn_adjusted <- ifelse(i == 1, segment_size, segment_size * cn_steps[i])
#     
#     #whats the mutation rate now
#     total_mut_rate <- segment_size_cn_adjusted*mutation_rate_bp
#     
#     #now let's sample some mutations
#     mutations_accumulated <- rpois(n = number_divisions_since_last_gain, lambda = total_mut_rate)
#     interval_muts  <- sum(mutations_accumulated) #add these up
#     
#     # store results
#     copy_number_before_gain <- ifelse(i == 1, 1, cn_steps[i])
#     copy_number_after_gain  <- ifelse(i == 1, 2, cn_steps[i + 1])
#     
#     df_out[i, "interval_index"] <- i
#     df_out[i, "mutation_rate_bp"] <- mutation_rate_bp
#     df_out[i, "copy_number_before_gain"] <- copy_number_before_gain
#     df_out[i, "copy_number_after_gain"] <- copy_number_after_gain
#     df_out[i, "segment_size"] <- segment_size
#     df_out[i, "present_dna_length_before_gain"] <- segment_size_cn_adjusted
#     df_out[i, "number_cell_divisions_since_last_gain"] <- number_divisions_since_last_gain
#     df_out[i, "number_muts_accumulated"] <- interval_muts
#     
#   }
#   
#   df_out$iteration <- iter
#   all_out[[iter]] <- df_out
# }
# 
# all_out <- do.call(rbind, all_out)






