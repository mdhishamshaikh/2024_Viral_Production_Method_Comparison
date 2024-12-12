# The aim of the script is to run calculate collision rates for  sampling station sin the coastal Dutch North Sea (CNS), open North Sea (ONS), and Curacao Caribbean Sea (CCS).
# 0. Setting up ####
source("./scripts/0_source.R")


# 1.0 Importing dataframes ####
counts_df <- read.csv("./results/vp_assays/viral_counts/filtered_counts_without_outliers.csv") %>%
  select(Location, Station_Number, Sample_Type, Timepoint, Replicate, c_Bacteria, c_Viruses, VBR)

# 2.0 Data wrangling and calculating collision rates
# Calculating  mean c_Bacteria and c_Viruses across replicates, then computing BV
bv_counts_df <- counts_df  %>%
  group_by(Location, Station_Number, Sample_Type, Timepoint) %>%
  summarise(
    mean_c_Bacteria = mean(c_Bacteria, na.rm = TRUE),
    mean_c_Viruses = mean(c_Viruses, na.rm = TRUE)
  ) %>%
  mutate(BV = mean_c_Bacteria * mean_c_Viruses) %>%
  ungroup()

# Initializing an empty list to store results
cr_list <- list()

# Calculating rates based on the mean BV values
for (loc in unique(bv_counts_df$Location)) {
  for (expt in unique(bv_counts_df$Station_Number)) {
    for (type in unique(bv_counts_df$Sample_Type)) {
      
    
      df_filtered <- bv_counts_df %>%
        filter(Location == loc, Station_Number == expt, Sample_Type == type)
      
      # Checking if Timepoint 0 exists to avoid empty divisions
      if (nrow(df_filtered %>% filter(Timepoint == 0)) > 0) {
        for (time in unique(df_filtered$Timepoint)) {
          # Only proceed if both the current timepoint and timepoint 0 data exist
          if (nrow(df_filtered %>% filter(Timepoint == time)) > 0) {
            cr_value <- (df_filtered %>% filter(Timepoint == time))$BV / 
              (df_filtered %>% filter(Timepoint == 0))$BV
            
            # Store the result as a list
            cr_list[[length(cr_list) + 1]] <- list(
              Location = loc,
              Station_Number = expt,
              Sample_Type = type,
              Timepoint = time,
              rate = cr_value
            )
          }
        }
      }
    }
  }
}

# Converting the list of lists to a data frame
cr_df_result <- do.call(rbind, lapply(cr_list, as.data.frame))

# Converting rate to numeric and handle any potential NA values
cr_df_result$rate <- as.numeric(cr_df_result$rate)


# Summarizing 
cr_df_summary <- cr_df_result %>%
  group_by(Location, Station_Number, Sample_Type, Timepoint) %>%
  summarise(rate_mean = mean(rate, na.rm = TRUE),
            rate_se = sd(rate, na.rm = TRUE) / sqrt(n())) %>%
  ungroup() %>%
  mutate(Location_Station = paste(Location, Station_Number, sep = "_"))

# Loading bacterial production data for endpoints
bep_df <- read_csv("./results/viralprod/vp_results_BP.csv") %>%
  select(#Location_Station, 
         Location, Station_Number, Time_Range) %>%
  unique() %>%
  mutate(endpoint = as.numeric(sub(".*\\T0_T?\\s*(\\d+)", "\\1", Time_Range)))


# Correcting timepoints
cr_opt <- cr_df_summary %>%
  mutate(Timepoint = case_when(
    Station_Number == 1 & Timepoint == 9 ~ 17,
    Station_Number == 2 & Timepoint == 9 ~ 17.5,
    Station_Number == 1 & Timepoint == 12 ~ 20,
    Station_Number == 2 & Timepoint == 12 ~ 20.5,
    TRUE ~ Timepoint
  )) %>%
  mutate(Timepoint = as.numeric(Timepoint))

#Merging collision rate and bacterial production data
cr_bep_df <- merge(cr_opt, bep_df) 

# 3.0 Visualization #####

# Setting colors and shapes
cols <- c("VP" = '#d32f27ff', "VPC" = '#4888a2ff', "Diff" = '#26547c')
shapes <- c("VP" = 15, "VPC" = 17, "Diff" = 19)

# Defining a plotting function for reusability
plot_cr <- function(data, cols, shapes) {
  ggplot(data) +
    geom_hline(yintercept = 1, color = 'gray') +
    geom_vline(data = data, aes(xintercept = endpoint), colour = '#7E6148FF', linewidth = 1.5, alpha = 0.5) +
    geom_point(aes(x = as.numeric(Timepoint), y = rate_mean, fill = Sample_Type, shape = Sample_Type, colour = Sample_Type), show.legend = TRUE) +
    theme_bw() +
    facet_wrap(. ~ Station_Number, nrow = 2) +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shapes) +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold', size = 10),
          #panel.border = element_rect(linewidth = 1),
          panel.background = element_rect(fill = NA),
          legend.title = element_text(face = 'bold', size = 10),
          legend.text = element_text(size = 9),
          axis.title = element_text(face = 'bold', size = 10),
          axis.text = element_text(size = 10)) +
    labs(x = NULL,
           #"Timepoint (in hours)", 
         y = NULL,
         #"Mean Relative Collision Rate"
         ) +
    guides(color = guide_legend(title = "Treatment"), shape = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment"))
}

# Plotting for all stations
cr_plot <- plot_cr(cr_bep_df, cols, shapes)
cr_plot
ggsave(cr_plot + theme(legend.position = 'none'), filename = "./figures/collision_rates_ngte.png", dpi = 800, width = 200, height = 80, units = "mm")


# Plotting for CNS stations
cr_cns_plot <- plot_cr(cr_bep_df %>%
                        dplyr::filter(Location == "CNS"), cols, shapes)
cr_cns_plot
ggsave(cr_cns_plot + theme(legend.position = 'none'), filename = "./figures/collision_rates_ngte_CNS.png", dpi = 800, width = 130, height = 80, units = "mm")


