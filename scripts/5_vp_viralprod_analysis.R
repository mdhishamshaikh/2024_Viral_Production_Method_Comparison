# The aim of the script is to remove outliers and run viralprod R package on DNS (coastal Dutch North Sea) data

# 0. Setting up ####
source("./scripts/0_source.R")

# 1.0 Importing per mL counts ####

nj_counts <- read.csv("./results/nj2020_fcm/nj2020_fcm_corrected_counts.csv")
pe_counts <- read.csv("./data/vp_assays/counts_per_ml/pe_cruises_vp_assay_fcm_counts_per_ml.csv")
cr_counts <- readxl::read_xlsx("./data/vp_assays/counts_per_ml/caribbean_vp_assay_fcm_counts_per_ml.xlsx")

# Combining all counts df
counts <- rbind(nj_counts, pe_counts, cr_counts) %>% 
  mutate(
    Original_Location = Location,
    Original_Station_Number = Station_Number,
    Location = case_when(
    Location == "NJ2020" ~ "DNS",
    Location == "PE477" ~ "ONS",
    Location == "PE486" ~ "ONS",
    Location == "CR" ~ "CCS",
    TRUE ~ Location  # Keep the original value if it doesn't match any condition
  ))  %>%
  group_by(Location) %>%
  mutate(
    Station_Number = dense_rank(Station_Number)
  ) %>%
  ungroup()


filtered_counts <- counts %>%
  filter(Sample_Type %in% c("VP", "VPC"),
         Staining_Protocol == "Viruses")  %>%
  mutate(Location_Station = paste(Location, Station_Number, sep = "_"))


# 2.0 Visualizing and filtering viral counts ####

plot_viral_abundance(file_path = "./results/vp_assays/viral_counts/vdc_plot_per_location_station.pdf",
                     data = filtered_counts,
                     location_column = "Location_Station")

# identifying potential outliers 
# Making an exclusion csv in Excel with colnames - Location, Station_Number, Sample_Type, Timepoint, Replicate
exclusion_df <- read.csv("./results/vp_assays/viral_counts/excluded_samples.csv")

# Filtering based on exclusion key
filtered_counts <- filtered_counts %>%
  mutate(Exclusion_Key = paste(Location, Station_Number, Timepoint, Replicate, Sample_Type, sep = "_"))
exclusion_df <- exclusion_df %>%
  mutate(Exclusion_Key = paste(Location, Station_Number, Timepoint, Replicate, Sample_Type, sep = "_"))

filtered_counts_without_outliers <- filtered_counts %>%
  mutate(c_Viruses = ifelse(Exclusion_Key %in% exclusion_df$Exclusion_Key, NA, c_Viruses)) %>%
  select(-Exclusion_Key)

# Plotting filtered counts
plot_viral_abundance(file_path = "./results/vp_assays/viral_counts/vdc_plot_per_location_station_without_outliers.pdf",
                     data = filtered_counts_without_outliers,
                     location_column = "Location_Station")



# Visualizing difference curves with VP and VPC
average_counts_diff_with_errors <- calculate_mean_viral_count_diff_with_errors(filtered_counts_without_outliers)

write.csv(average_counts_diff_with_errors, "./results/vp_assays/viral_counts/average_counts_diff_with_errors.csv", row.names = F)

plot_mean_viral_abundance(file_path = "./results/vp_assays/viral_counts/vdc_plot_with_diff_per_location_station_without_outliers.pdf",
                     data = average_counts_diff_with_errors,
                     location_column = "Location_Station")


# Visualizing difference curves only
diff_curves_plot <- ggplot(data = average_counts_diff_with_errors %>% dplyr::filter(Sample_Type == "Diff"),
       aes(x = Timepoint, y = avg_c_Viruses / 1e+6)) +
  geom_line(size = 0.7) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = (avg_c_Viruses - sem_c_Viruses) / 1e+6,
                    ymax = (avg_c_Viruses + sem_c_Viruses) / 1e+6), width = 0.2) +
  scale_color_manual(values = c("VP" = "#ef476f", "VPC" = "#ffd166", "Diff" = "#26547c")) +
  facet_wrap(Location_Station ~ ., scale = "free") +
  theme_bw() +
  labs(
    #title = paste("Average Viral Abundance for Location Station:", loc_station),
    x = "Timepoint",
    y = "Average Viral Count (c_Viruses in millions)",
    color = "Sample Type"
  )
diff_curves_plot
ggsave(filename = "difference_curves.png", plot = diff_curves_plot, path = "./results/vp_assays/viral_counts", dpi = 1000)



# Saving the dataframe with outliers removed 

write.csv(filtered_counts_without_outliers, "./results/vp_assays/viral_counts/filtered_counts_without_outliers.csv", row.names = F)



# 3.0 Importing Abundance data ####

abundance <- read.csv("./data/metadata/bacterial_viral_abundances_abiotic.csv")


# 4.0 Running viralprod #####

vp_counts <- read.csv("./results/vp_assays/viral_counts/filtered_counts_without_outliers.csv")

#Only keeping c_Bacteria and c_viruses
vp_counts <-  vp_counts %>%
  dplyr::select(-c("c_HNA", "c_LNA", "c_V1", "c_V2", "c_V3"))

# Checking populations
vp_check_populations(vp_counts)
class(abundance)
abundance <- vp_class_ori_abu(abundance)
class(abundance)

class(vp_counts) #failed
vp_counts <- vp_class_count_data(vp_counts) #passed
class(vp_counts)

#Running viralprod

vp_end_to_end(data = vp_counts ,
              original_abundances = abundance,
              methods = c(2,9,10),
              write_output = T,
              output_dir = './results/viralprod2')

