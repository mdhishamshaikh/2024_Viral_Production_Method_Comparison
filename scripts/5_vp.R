# The aim of the script is to remove outliers and run viralprod R package on 
# DNS (coastal Dutch North Sea), ONS (Open North Sea), and CCS ( Curacao Caribbean Sea)


# 0. Setting up ####
source("./scripts/0_source.R")

# 1.0 Importing per mL counts ####

nj_counts <- read.csv("./results/nj2020_fcm/nj2020_fcm_corrected_counts.csv")
pe_counts <- read.csv("./data/vp_assays/counts_per_ml/pe_cruises_vp_assay_fcm_counts_per_mL.csv")
cr_counts <- readxl::read_xlsx("./data/vp_assays/counts_per_ml/caribbean_vp_assay_fcm_counts_per_mL.xlsx")

# Combining all counts df
counts <- rbind(nj_counts, pe_counts, cr_counts) %>% 
  mutate(Location2 = case_when(
  Location == "NJ2020" ~ "DNS",
  Location == "PE477" ~ "ONS",
  Location == "PE486" ~ "ONS",
  Location == "CR" ~ "CCS",
  TRUE ~ Location  # Keep the original value if it doesn't match any condition
))  %>%
  group_by(Location2) %>%
  mutate(
    Station_Number2 = dense_rank(Station_Number)
  ) %>%
  ungroup()

filtered_counts <- counts %>%
  filter(Sample_Type %in% c("VP", "VPC"),
         Staining_Protocol == "Viruses")  %>%
  mutate(Location_Station = paste(Location2, Station_Number2, sep = "_"))

# 2.0 Visualizing viral counts ####

vdc_plot <- ggplot(filtered_counts, aes(x = Timepoint, y = c_Viruses/1e+6)) +
  geom_line(aes(group = as.factor(Replicate))) +  # Add lines to connect points for each sample
  geom_point(aes(color = as.factor(Replicate)), size = 2) +  # Color points by Sample_Type
  facet_grid(Sample_Type ~ Location_Station, scales = "free") +  # Facet by Location_Station and Sample_Type
  scale_color_npg() +
  theme_bw(base_size = 15) +
  labs(#title = "Viral Abundance (c_Viruses) over Timepoints",
    x = "Timepoint",
    y = "Viral Count (c_Viruses)",
    color = "Sample Type")
vdc_plot
ggsave(vdc_plot, filename = "./figures/vdc_plot.svg", width = 20, height = 5)
