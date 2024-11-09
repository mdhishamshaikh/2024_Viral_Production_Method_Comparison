library(tidyverse)
library(ggsci)
library(patchwork)

# Here I will create individual plots for each scenario

# Plot function ####

generate_plots <- function(data, save_folder) {
  
  # Create the save folder if it does not exist
  if (!dir.exists(save_folder)) {
    dir.create(save_folder, recursive = TRUE)
  }
  
  # Summarize data
  summary_df <- data %>%
    group_by(Timepoint, Treatment) %>%
    summarise(
      Mean = mean(Count),
      SEM = sd(Count) / sqrt(n())
    ) %>%
    mutate(Treatment = factor(Treatment, levels = c("VP", "VPC", "Diff"))) %>%
    ungroup()
  
  # Prepare data for plotting
  lytic_df <- summary_df %>% dplyr::filter(Treatment == 'VP')
  lytic_lysogenic_df <- summary_df %>% dplyr::filter(Treatment == 'VPC')
  combined_df <- inner_join(lytic_df, lytic_lysogenic_df, by = "Timepoint", suffix = c("_Lytic", "_Lytic_Lysogenic"))
  combined_df <- combined_df %>%
    mutate(
      Diff = Mean_Lytic_Lysogenic - Mean_Lytic,
      Diff_SEM = sqrt(SEM_Lytic_Lysogenic^2 + SEM_Lytic^2)
    )
  
  summary_df <- bind_rows(
    summary_df,
    combined_df %>%
      select(Timepoint, Treatment = Diff, Mean = Diff, SEM = Diff_SEM) %>%
      mutate(Treatment = "Diff")
  ) %>%
    mutate(Treatment = factor(Treatment, levels = c("VP", "VPC", "Diff")))
  
  # Define plot function
  plot_and_save <- function(plot, filename) {
    ggsave(plot = plot, filename = file.path(save_folder, filename), dpi = 800, width = 70, height = 70, units = "mm")
  }
  
  # Generate individual plots
  lm_vp_plot <- ggplot(data %>% dplyr::filter(Treatment == 'VP'), aes(x = Timepoint, y = Count)) +
    geom_point() +
    geom_smooth(aes(group = Replicate), method = "lm", se = FALSE, color = "#ef476f") +
    xlab(NULL) +
    ylab(NULL) +
    ylim(c(4, 8)) +
    theme_classic(base_size = 18)
  
  lm_vpc_plot <- ggplot(data %>% dplyr::filter(Treatment == 'VPC'), aes(x = Timepoint, y = Count)) +
    geom_point() +
    geom_smooth(aes(group = Replicate), method = "lm", se = FALSE, color = "#ffd166") +
    xlab(NULL) +
    ylab(NULL) +
    ylim(c(4, 8)) +
    theme_classic(base_size = 18)
  
  lm_diff_plot <- ggplot(summary_df %>% dplyr::filter(Treatment == 'Diff'), aes(x = Timepoint, y = Mean)) +
    geom_point(alpha = 0.3) +
    geom_line(data = summary_df %>% dplyr::filter(Treatment == 'Diff'), aes(x = Timepoint, y = Mean, group = 1), color = "#26547c", inherit.aes = F) +
    xlab(NULL) +
    ylab(NULL) +
    ylim(c(-1.5, 1.5)) +
    theme_classic(base_size = 18)
  
  vipcal_vp_plot <- ggplot(data %>% dplyr::filter(Treatment == 'VP'), aes(x = Timepoint, y = Count)) +
    geom_point(alpha = 0.3) +
    geom_line(data = summary_df %>% dplyr::filter(Treatment == 'VP'), aes(x = Timepoint, y = Mean, group = 1), color = "#ef476f", inherit.aes = F) +
    geom_point(data = summary_df %>% dplyr::filter(Treatment == 'VP'), aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 15) +
    xlab(NULL) +
    ylab(NULL) +
    ylim(c(4, 8)) +
    theme_classic(base_size = 18)
  
  vipcal_vpc_plot <- ggplot(data %>% dplyr::filter(Treatment == 'VPC'), aes(x = Timepoint, y = Count)) +
    geom_point(alpha = 0.3) +
    geom_line(data = summary_df %>% dplyr::filter(Treatment == 'VPC'), aes(x = Timepoint, y = Mean, group = 1), color = "#ffd166", inherit.aes = F) +
    geom_point(data = summary_df %>% dplyr::filter(Treatment == 'VPC'), aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 19) +
    xlab(NULL) +
    ylab(NULL) +
    ylim(c(4, 8)) +
    theme_classic(base_size = 18)
  
  vipcal_diff_plot <- ggplot(data %>% dplyr::filter(Treatment == 'Diff'), aes(x = Timepoint, y = Count)) +
    geom_point(alpha = 0.3) +
    geom_line(data = summary_df %>% dplyr::filter(Treatment == 'Diff'), aes(x = Timepoint, y = Mean, group = 1), color = "#26547c", inherit.aes = F) +
    geom_point(data = summary_df %>% dplyr::filter(Treatment == 'Diff'), aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 17) +
    xlab(NULL) +
    ylab(NULL) +
    ylim(c(-1.5, 1.5)) +
    theme_classic(base_size = 18)
  
  vipcal_se_vp_plot <- ggplot(data %>% dplyr::filter(Treatment == 'VP'), aes(x = Timepoint, y = Count)) +
    geom_point(alpha = 0.3) +
    geom_line(data = summary_df %>% dplyr::filter(Treatment == 'VP'), aes(x = Timepoint, y = Mean, group = 1), color = "#ef476f", inherit.aes = F) +
    geom_point(data = summary_df %>% dplyr::filter(Treatment == 'VP'), aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 15) +
    geom_pointrange(data = summary_df %>% dplyr::filter(Treatment == 'VP'), aes(x = Timepoint, y = Mean, ymin = Mean - SEM, ymax = Mean + SEM), color = "black", inherit.aes = F) +
    xlab(NULL) +
    ylab(NULL) +
    ylim(c(4, 8)) +
    theme_classic(base_size = 18)
  
  vipcal_se_vpc_plot <- ggplot(data %>% dplyr::filter(Treatment == 'VPC'), aes(x = Timepoint, y = Count)) +
    geom_point(alpha = 0.3) +
    geom_line(data = summary_df %>% dplyr::filter(Treatment == 'VPC'), aes(x = Timepoint, y = Mean, group = 1), color = "#ffd166", inherit.aes = F) +
    geom_point(data = summary_df %>% dplyr::filter(Treatment == 'VPC'), aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 19) +
    geom_pointrange(data = summary_df %>% dplyr::filter(Treatment == 'VPC'), aes(x = Timepoint, y = Mean, ymin = Mean - SEM, ymax = Mean + SEM), color = "black", inherit.aes = F) +
    xlab(NULL) +
    ylab(NULL) +
    ylim(c(4, 8)) +
    theme_classic(base_size = 18)
  
  vipcal_se_diff_plot <- ggplot(data %>% dplyr::filter(Treatment == 'Diff'), aes(x = Timepoint, y = Count)) +
    geom_point(alpha = 0.3) +
    geom_line(data = summary_df %>% dplyr::filter(Treatment == 'Diff'), aes(x = Timepoint, y = Mean, group = 1), color = "#26547c", inherit.aes = F) +
    geom_point(data = summary_df %>% dplyr::filter(Treatment == 'Diff'), aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 17) +
    geom_pointrange(data = summary_df %>% dplyr::filter(Treatment == 'Diff'), aes(x = Timepoint, y = Mean, ymin = Mean - SEM, ymax = Mean + SEM), color = "black", shape = 17, inherit.aes = F) +
    xlab(NULL) +
    ylab(NULL) +
    ylim(c(-1.5, 1.5)) +
    theme_classic(base_size = 18) 
  
  # Create a cumulative plot layout using patchwork
  cumulative_plot <- (lm_vp_plot | lm_vpc_plot | lm_diff_plot) /
    (vipcal_vp_plot | vipcal_vpc_plot | vipcal_diff_plot) /
    (vipcal_se_vp_plot | vipcal_se_vpc_plot | vipcal_se_diff_plot)
  
  
  # Saving individual plots
  plot_and_save(lm_vp_plot, "lm_vp_plot.png")
  plot_and_save(lm_vpc_plot, "lm_vpc_plot.png")
  plot_and_save(lm_diff_plot, "lm_diff_plot.png")
  plot_and_save(vipcal_vp_plot, "vipcal_vp_plot.png")
  plot_and_save(vipcal_vpc_plot, "vipcal_vpc_plot.png")
  plot_and_save(vipcal_diff_plot, "vipcal_diff_plot.png")
  plot_and_save(vipcal_se_vp_plot, "vipcal_se_vp_plot.png")
  plot_and_save(vipcal_se_vpc_plot, "vipcal_se_vpc_plot.png")
  plot_and_save(vipcal_se_diff_plot, "vipcal_se_diff_plot.png")
  
  # Save the cumulative plot
  ggsave(
    plot = cumulative_plot, 
    filename = file.path(save_folder, "cumulative_plot.png"), 
    dpi = 800, 
    width = 3 * 70,  # Triple the width
    height = 3 * 70, # Triple the height
    units = "mm"
  )
}


# Creating a skeleton dataframe ####

mock_df <-  tibble(
  Treatment = rep(c("VP", "VPC"), each = 18),    # 2 treatments, each repeated 18 times (3 replicates * 6 timepoints)
  Replicate = rep(rep(1:3, each = 6), times = 2), # 3 replicates, each repeated 6 times, for each treatment
  Timepoint = rep(c("T1", "T2", "T3", "T4", "T5", "T6"), times = 6)                 # 6 timepoints, repeated for each replicate and treatment
)


#### Scenario 1: Perfect replicates, steady increase in VP and VPC, VPC is higher than VP ####

s1_count <- c(5.0, 5.5, 6.1, 6.4, 6.8, 7.3,
              5.0, 5.5, 6.1, 6.4, 6.8, 7.3,
              5.0, 5.5, 6.1, 6.4, 6.8, 7.3,
              5.0, 5.6, 6.3, 6.8, 7.1, 7.8,
              5.0, 5.6, 6.3, 6.8, 7.1, 7.8,
              5.0, 5.6, 6.3, 6.8, 7.1, 7.8)

s1_mock_df <- mock_df %>% mutate(Count = s1_count)

generate_plots(s1_mock_df, "./figures/scenarios/s1")


#### Scenario 2: Perfect replicates, steady increase in VP and VPC, VP is higher than VPC ####

s2_count <- c(5.0, 5.6, 6.3, 6.8, 7.1, 7.8,
              5.0, 5.6, 6.3, 6.8, 7.1, 7.8,
              5.0, 5.6, 6.3, 6.8, 7.1, 7.8,
              5.0, 5.5, 6.1, 6.4, 6.8, 7.3,
              5.0, 5.5, 6.1, 6.4, 6.8, 7.3,
              5.0, 5.5, 6.1, 6.4, 6.8, 7.3)

s2_mock_df <- mock_df %>% mutate(Count = s2_count)

generate_plots(s2_mock_df, "./figures/scenarios/s2")


