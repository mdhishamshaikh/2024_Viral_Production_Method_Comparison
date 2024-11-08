# Loading essential libraries.

source("./scripts/0_source.R")
packages_to_load <- c("tidyverse",
                      "ggsci")
lapply(packages_to_load, load_or_install)


# 1.0 Importing the literature survey sheet ####

lit_df <- readxl::read_xlsx("./data/vp_literature_survey/2024_vp_assay_literature_survey.xlsx")

# 2.0 Data wrangling ####

# Filtering for 'Keep'
lit_df <- lit_df %>%
  dplyr::filter(Keep == 1)

# Creating a tally dataframe
tally_df <- lit_df %>%
  pivot_longer(cols = c(Lytic, Lysogenic), names_to = "Study_Type", values_to = "Assay_Count") %>%
  pivot_longer(cols = c(LR, VIPCAL), names_to = "Method_Used", values_to = "Method_Count") %>%
  dplyr::filter(Assay_Count > 0, Method_Count > 0) %>%
  #select(-Method_Value) %>%      
  group_by(Year_Published, Study_Type, Method_Used) %>%
  summarise(Study_Count = sum(Method_Count, na.rm = TRUE), .groups = 'drop')

tally_df2 <- lit_df %>%
  pivot_longer(cols = c(Lytic, Lysogenic), names_to = "Study_Type", values_to = "Assay_Count") %>%
  pivot_longer(cols = c(LR, VIPCAL), names_to = "Method_Used", values_to = "Method_Count") %>%
  filter(Assay_Count > 0, Method_Count > 0) %>%
  #select(-Method_Value) %>%      
  group_by(Year_Published, Study_Type, Method_Used) %>%
  summarise(Study_Count = sum(Method_Count, na.rm = TRUE), .groups = 'drop')


# Finding the range of years
year_range <- seq(min(tally_df$Year_Published), max(tally_df$Year_Published))

# Expand the df to include all years for each Study_Type and Method_Used combination
expanded_df <- expand.grid(
  Year_Published = year_range,
  Study_Type = unique(tally_df$Study_Type),
  Method_Used = unique(tally_df$Method_Used)
)

# Merging with the original data to fill missing years with NA
tally_df_expanded <- expanded_df %>%
  dplyr::left_join(tally_df, by = c("Year_Published", "Study_Type", "Method_Used")) %>%
  # Replacing NA Study_Count with 0 (no studies in those years)
  dplyr::mutate(Study_Count = replace_na(Study_Count, 0)) %>%
  # Grouping again to calculate cumulative count, ensuring all years are accounted for
  dplyr::group_by(Study_Type, Method_Used) %>%
  dplyr::arrange(Year_Published) %>%
  dplyr::mutate(Cumulative_Count = cumsum(Study_Count)) %>%
  dplyr::ungroup()

print(head(tally_df_expanded))

# Combining 'Method_Used' and 'Study_Type' into a single variable in the expanded DataFrame
tally_df_expanded <- tally_df_expanded %>%
  mutate(Combined = paste(Method_Used, Study_Type, sep = "_")) 
  



# 3.0 Visualization ####

# Ensure 'Combined' is a factor for consistent plotting
tally_df_expanded$Combined <- factor(tally_df_expanded$Combined)
colors <- c("LR" = "#DC2828", "VIPCAL" = "#003049")
shapes <- c("Lytic" = 15, "Lysogenic" = 17)



# Plot the cumulative graph using the expanded data
vp_lit_plot<- ggplot(tally_df_expanded, aes(x = Year_Published, y = Cumulative_Count)) +
  geom_line(aes(color = Method_Used, group = interaction(Method_Used, Study_Type)), size = 0.5) +
  geom_point(aes(shape = Study_Type, color = Method_Used), size = 2.0) +
  scale_shape_manual(values = shapes, name = "Viral production") +
  scale_color_manual(values = colors, name = "Analytical approach") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


vp_lit_plot

# Saving plots
ggsave(plot = vp_lit_plot, filename = "./figures/vp_lit_plot_legend.svg", dpi = 800, width = 160, height = 70, units = "mm")
ggsave(plot = vp_lit_plot + theme(legend.position = "none"), filename = "./figures/vp_lit_plot.png", dpi = 800, width = 160, height = 70, units = "mm")

