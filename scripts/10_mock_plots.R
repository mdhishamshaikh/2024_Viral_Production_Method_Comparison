# The aim of this script is to generate explanatory figures for LM, VIPCAL, and VIPCAL-SE

library(tidyverse)
library(ggsci)
#  Let's make a function to create a set of mockplots based on a given dataset. I will update the dataset later #


mock_df <-  tibble(
  Treatment = rep(c("VP", "VPC"), each = 18),    # 2 treatments, each repeated 18 times (3 replicates * 6 timepoints)
  Replicate = rep(rep(1:3, each = 6), times = 2), # 3 replicates, each repeated 6 times, for each treatment
  Timepoint = rep(c("T1", "T2", "T3", "T4", "T5", "T6"), times = 6)                 # 6 timepoints, repeated for each replicate and treatment
)

# Now adding viral counts

mock_df$Count <- c(5.0, 4.2, 5.9, 5.4, 6.5, 6.1,
                   4.9, 4.1, 6.3, 5.3, 6.7, 6.8,
                   5.3, 4.3, 6.4, 5.8, 6.1, 7.9,
                   4.7, 4.2, 5.8, 5.9, 6.4, 7.3,
                   4.6, 4.5, 6.2, 5.1, 6.3, 7.8,
                   5.0, 4.4, 6.3, 6.2, 5.7, 7.9)



{
summary_df <-  mock_df %>%
  group_by(Timepoint, Treatment) %>%
  summarise(
    Mean = mean(Count),
    SEM = sd(Count) / sqrt(n())
  ) %>%
    mutate(Treatment = factor(Treatment, levels = c("VP", "VPC", "Diff")))%>%
  ungroup()


lytic_df <- summary_df %>% filter(Treatment == 'VP')
lytic_lysogenic_df <- summary_df %>% filter(Treatment == 'VPC')
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



summary_df
}

# Base plot

ggplot(mock_df, aes(x = Timepoint, y = Count, color = as.factor(Replicate))) +
  geom_point() +
  geom_smooth(aes(group = Replicate), method = "lm", se = FALSE) + 
  geom_point(data = summary_df, aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 17) +
  geom_errorbar(data = summary_df, aes(x = Timepoint, ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2, color = "black", inherit.aes = F) +
  facet_grid(cols = vars(Treatment), scales = "fixed") +
 # ylim(c(3, 8)) +
  theme_classic(base_size = 18)





# Viral counts over time

vp_plot<- ggplot(mock_df %>% dplyr::filter(Treatment == 'VP'), aes(x = Timepoint, y = Count)) +
  geom_point() +
  #geom_smooth(aes(group = Replicate), method = "lm", se = FALSE, color = "#ef476f") + 
  #geom_point(data = summary_df, aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 17) +
  #geom_errorbar(data = summary_df, aes(x = Timepoint, ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2, color = "black", inherit.aes = F) +
  #facet_grid(cols = vars(Treatment), scales = "fixed") +
  ylim(c(4, 8)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_classic(base_size = 18)

vp_plot

ggsave(plot = vp_plot, filename = "./figures/vp_plot.png", dpi = 800, width = 70, height = 70, units = "mm")

vpc_plot<- ggplot(mock_df %>% dplyr::filter(Treatment == 'VPC'), aes(x = Timepoint, y = Count)) +
  geom_point() +
  #geom_smooth(aes(group = Replicate), method = "lm", se = FALSE, color = "#ffd166") + 
  #geom_point(data = summary_df, aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 17) +
  #geom_errorbar(data = summary_df, aes(x = Timepoint, ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2, color = "black", inherit.aes = F) +
  #facet_grid(cols = vars(Treatment), scales = "fixed") +
  ylim(c(4, 8)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_classic(base_size = 18)

vpc_plot

ggsave(plot = vpc_plot, filename = "./figures/vpc_plot.png", dpi = 800, width = 70, height = 70, units = "mm")






# LM plots ####
# LM Lytic 


lm_vp_plot<- ggplot(mock_df %>% dplyr::filter(Treatment == 'VP') %>%dplyr::mutate(Count = if_else(Timepoint == 'T1', NA_real_, Count)), aes(x = Timepoint, y = Count)) +
  geom_point() +
  geom_smooth(aes(group = Replicate), method = "lm", se = FALSE, color = "#ef476f") + 
  #geom_point(data = summary_df, aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 17) +
  #geom_errorbar(data = summary_df, aes(x = Timepoint, ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2, color = "black", inherit.aes = F) +
  #facet_grid(cols = vars(Treatment), scales = "fixed") +
  ylim(c(4, 8)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_classic(base_size = 18)

lm_vp_plot

ggsave(plot = lm_vp_plot, filename = "./figures/lm_vp_plot.png", dpi = 800, width = 70, height = 70, units = "mm")

# LM Lytic + Lysogenic
lm_vpc_plot<- ggplot(mock_df %>% dplyr::filter(Treatment == 'VPC') %>%dplyr::mutate(Count = if_else(Timepoint == 'T1', NA_real_, Count)), aes(x = Timepoint, y = Count)) +
  geom_point() +
  geom_smooth(aes(group = Replicate), method = "lm", se = FALSE, color = "#ffd166") + 
  #geom_point(data = summary_df, aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 17) +
  #geom_errorbar(data = summary_df, aes(x = Timepoint, ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2, color = "black", inherit.aes = F) +
  #facet_grid(cols = vars(Treatment), scales = "fixed") +
  ylim(c(4, 8)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_classic(base_size = 18)

lm_vpc_plot

ggsave(plot = lm_vpc_plot, filename = "./figures/lm_vpc_plot.png", dpi = 800, width = 70, height = 70, units = "mm")



# vipcal plots ####
# VIPCAL Lytic 


vipcal_vp_plot<- ggplot(mock_df %>% dplyr::filter(Treatment == 'VP'), aes(x = Timepoint, y = Count)) +
  geom_point(alpha = 0.3) +
  #geom_smooth(aes(group = Replicate), method = "lm", se = FALSE, color = "#ef476f") + 
  geom_line(data = summary_df %>% dplyr::filter(Treatment == 'VP'), aes(x = Timepoint, y = Mean, group = 1), color = "#ef476f", inherit.aes = F, linewidth = 1) +
  geom_point(data = summary_df %>% dplyr::filter(Treatment == 'VP'), aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 15) +
  #geom_errorbar(data = summary_df, aes(x = Timepoint, ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2, color = "black", inherit.aes = F) +
  #facet_grid(cols = vars(Treatment), scales = "fixed") +
  ylim(c(4, 8)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_classic(base_size = 18)

vipcal_vp_plot

ggsave(plot = vipcal_vp_plot, filename = "./figures/vipcal_vp_plot.png", dpi = 800, width = 70, height = 70, units = "mm")


vipcal_vpc_plot<- ggplot(mock_df %>% dplyr::filter(Treatment == 'VPC'), aes(x = Timepoint, y = Count)) +
  geom_point(alpha = 0.3) +
  #geom_smooth(aes(group = Replicate), method = "lm", se = FALSE, color = "#ef476f") + 
  geom_line(data = summary_df %>% dplyr::filter(Treatment == 'VPC'), aes(x = Timepoint, y = Mean, group = 1), color = "#ffd166", inherit.aes = F, linewidth = 1) +
  geom_point(data = summary_df %>% dplyr::filter(Treatment == 'VPC'), aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 19) +
  #geom_errorbar(data = summary_df, aes(x = Timepoint, ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2, color = "black", inherit.aes = F) +
  #facet_grid(cols = vars(Treatment), scales = "fixed") +
  ylim(c(4, 8)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_classic(base_size = 18)

vipcal_vpc_plot

ggsave(plot = vipcal_vpc_plot, filename = "./figures/vipcal_vpc_plot.png", dpi = 800, width = 70, height = 70, units = "mm")


vipcal_diff_plot<- ggplot(mock_df %>% dplyr::filter(Treatment == 'Diff'), aes(x = Timepoint, y = Count)) +
  geom_point(alpha = 0.3) +
  #geom_smooth(aes(group = Replicate), method = "lm", se = FALSE, color = "#ef476f") + 
  geom_line(data = summary_df %>% dplyr::filter(Treatment == 'Diff'), aes(x = Timepoint, y = Mean, group = 1), color = "#26547c", inherit.aes = F, linewidth = 1) +
  geom_point(data = summary_df %>% dplyr::filter(Treatment == 'Diff'), aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 17) +
  #geom_errorbar(data = summary_df, aes(x = Timepoint, ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2, color = "black", inherit.aes = F) +
  #facet_grid(cols = vars(Treatment), scales = "fixed") +
  ylim(c(-1.5, 1.5)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_classic(base_size = 18)

vipcal_diff_plot

ggsave(plot = vipcal_diff_plot, filename = "./figures/vipcal_diff_plot.png", dpi = 800, width = 70, height = 70, units = "mm")



# vipcal-se plots ####
# VIPCAL-SE Lytic 


vipcal_se_vp_plot<- ggplot(mock_df %>% dplyr::filter(Treatment == 'VP'), aes(x = Timepoint, y = Count)) +
  geom_point(alpha = 0.3) +
  #geom_smooth(aes(group = Replicate), method = "lm", se = FALSE, color = "#ef476f") + 
  geom_line(data = summary_df %>% dplyr::filter(Treatment == 'VP'), aes(x = Timepoint, y = Mean, group = 1), color = "#ef476f", inherit.aes = F, linewidth = 1) +
  geom_point(data = summary_df %>% dplyr::filter(Treatment == 'VP'), aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 15) +
  geom_pointrange(data = summary_df %>% dplyr::filter(Treatment == 'VP'), aes(x = Timepoint, y = Mean, ymin = Mean - SEM, ymax = Mean + SEM),  color = "black", inherit.aes = F) +
  #facet_grid(cols = vars(Treatment), scales = "fixed") +
  ylim(c(4, 8)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_classic(base_size = 18)

vipcal_se_vp_plot

ggsave(plot = vipcal_se_vp_plot, filename = "./figures/vipcal_se_vp_plot.png", dpi = 800, width = 70, height = 70, units = "mm")


vipcal_se_vpc_plot<- ggplot(mock_df %>% dplyr::filter(Treatment == 'VPC'), aes(x = Timepoint, y = Count)) +
  geom_point(alpha = 0.3) +
  #geom_smooth(aes(group = Replicate), method = "lm", se = FALSE, color = "#ef476f") + 
  geom_line(data = summary_df %>% dplyr::filter(Treatment == 'VPC'), aes(x = Timepoint, y = Mean, group = 1), color = "#ffd166", inherit.aes = F, linewidth = 1) +
  geom_point(data = summary_df %>% dplyr::filter(Treatment == 'VPC'), aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 19) +
  geom_pointrange(data = summary_df %>% dplyr::filter(Treatment == 'VPC'), aes(x = Timepoint, y = Mean, ymin = Mean - SEM, ymax = Mean + SEM),  color = "black", inherit.aes = F) +
  #facet_grid(cols = vars(Treatment), scales = "fixed") +
  ylim(c(4, 8)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_classic(base_size = 18)

vipcal_se_vpc_plot

ggsave(plot = vipcal_se_vpc_plot, filename = "./figures/vipcal_se_vpc_plot.png", dpi = 800, width = 70, height = 70, units = "mm")


vipcal_se_diff_plot<- ggplot(mock_df %>% dplyr::filter(Treatment == 'Diff'), aes(x = Timepoint, y = Count)) +
  geom_point(alpha = 0.3) +
  #geom_smooth(aes(group = Replicate), method = "lm", se = FALSE, color = "#ef476f") + 
  geom_line(data = summary_df %>% dplyr::filter(Treatment == 'Diff'), aes(x = Timepoint, y = Mean, group = 1), color = "#26547c", inherit.aes = F, linewidth = 1) +
  geom_point(data = summary_df %>% dplyr::filter(Treatment == 'Diff'), aes(x = Timepoint, y = Mean), color = "black", size = 3, shape = 17) +
  geom_pointrange(data = summary_df %>% dplyr::filter(Treatment == 'Diff'), aes(x = Timepoint, y = Mean, ymin = Mean - SEM, ymax = Mean + SEM),  color = "black", shape = 17, inherit.aes = F) +
  #facet_grid(cols = vars(Treatment), scales = "fixed") +
  ylim(c(-1.5, 1.5)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_classic(base_size = 18) 

vipcal_se_diff_plot

ggsave(plot = vipcal_se_diff_plot, filename = "./figures/vipcal_se_diff_plot.png", dpi = 800, width = 70, height = 70, units = "mm")




