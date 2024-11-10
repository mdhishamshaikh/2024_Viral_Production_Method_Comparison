source("./scripts/0_source.R")


# Importing data
# LM supervised
# VIPCAL T24
#VIPCAL-SE GTE

# 
lr_24<- readxl::read_xlsx("./data/vp_assays/lm_vp_assays/dns_ons_ccs_vp_assays_lm.xlsx") %>%
  dplyr::select(-c(bac_efficiency, dplyr::starts_with("corrected"))) %>%
  rename_with(~ str_replace_all(., c("Lytic" = "VP", "Lysogenic" = "Diff"))) %>%
  pivot_longer(
    cols = c(VP, Diff),
    names_to = "Sample_Type",
    values_to = "VP"
  )%>%
  mutate(Time_Range = case_when(
    Sample_Type == "VP" ~ VP_Time_Range,
    Sample_Type == "Diff" ~ Diff_Time_Range
  )) %>%
  select(-VP_Time_Range, -Diff_Time_Range, -station_code)%>%
  mutate(VP_Method = 'LR')

vpcl_24<- read.csv("./results/viralprod/vp_results_24.csv") %>%
  mutate(Timepoint = 'T24') %>%
  dplyr::filter(VP_Method %in% c("VPCL_AR_DIFF"),
                Population == "c_Viruses",
                Sample_Type != "VPC") %>%
  select(c("Location", "Station_Number", "VP_Method", "Sample_Type", "VP", "VP_SE", "Time_Range")) %>%
  mutate(VP_Method = 'VIPCAL')

vpcl_se_24<- read.csv("./results/viralprod/vp_results_24.csv") %>%
  mutate(Timepoint = 'T24') %>%
  dplyr::filter(VP_Method %in% c("VPCL_AR_DIFF_SE"),
                Population == "c_Viruses",
                Sample_Type != "VPC") %>%
  select(c("Location", "Station_Number", "VP_Method", "Sample_Type", "VP", "VP_SE", "Time_Range")) %>%
  mutate(VP_Method = 'VIPCAL-SE')

vpcl_se_gte<- read.csv("./results/viralprod/vp_results_BP.csv") %>%
  mutate(Timepoint = 'GTE') %>%
  dplyr::filter(VP_Method == "VPCL_AR_DIFF_SE",
                Population == "c_Viruses",
                Sample_Type != "VPC") %>%
  select(c("Location", "Station_Number", "VP_Method", "Sample_Type", "VP", "VP_SE", "Time_Range"))%>%
  mutate(VP_Method = 'VIPCAL-SE-GTE')





all_vp<- dplyr:: bind_rows(lr_24, vpcl_24, vpcl_se_24, vpcl_se_gte)  %>%
  mutate(Sample_Type = recode(Sample_Type,
                              "VP" = "Lytic",
                              "Diff" = "Lysogenic"),
         Sample_Type = factor(Sample_Type, levels = c("Lytic", "Lysogenic")),
         Time_Range = factor(Time_Range, levels = c("T0_T3", "T0_T6",
                                                    "T0_T9", "T0_T12",
                                                    "T0_T17", "T0_T20",
                                                    "T0_T17.5", "T0_T20.5",
                                                    "T0_T24")))

# Adding station codes
station_legends <- read.csv("./results/vp_assays/viral_counts/station_legends.csv") %>%
  dplyr::select(-c(Location, Station_rank))  %>%
  arrange(Station_Number) %>%
  mutate(Station_code = factor(Station_code, levels = unique(Station_code)))

all_vp <- all_vp %>%
  dplyr::left_join(station_legends, by = "Station_Number") %>%
  mutate(
    VP = if_else(VP < 0, 0, VP),
    SE = if_else(VP == 0, 0, VP_SE)
  )

# Replacing negative values with zero



lm_vs_vpcl_vs_vpcl_se<- ggplot(all_vp  %>% dplyr::filter(Sample_Type != "VPC") ,  
       aes(x = as.factor(Station_code), y = VP/1e+6, fill = VP_Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = "black") +
 geom_errorbar(aes(ymin = (VP - VP_SE)/1e+6, ymax = (VP + VP_SE)/1e+6),
               position = position_dodge(width = 0.7), width = 0.25) +
  geom_hline(yintercept = 0) +
  labs(#title = "NJ2020 Viral Production",
    y = expression("Viral production rate" ~ (x ~ 10^6 ~ "VLPs" ~ mL^-1 ~ h^-1)),
    x = "Sampling event",
    x = NULL,
    y = NULL,
    fill = "Analytical approach") +
  scale_fill_manual(values = c("LR" = "#DC2828", "VIPCAL" = "#003049", "VIPCAL-SE" ="#2A9D8F", "VIPCAL-SE-GTE" ="#F6CD61")) +
  #scale_y_continuous(labels = scales::label_math(expr = 10^.x)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
  facet_wrap(~ Sample_Type, ncol = 1, scales = "free") +  # Use 'facet_wrap' with free y-axis scaling
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 90))



lm_vs_vpcl_vs_vpcl_se

ggsave(plot = lm_vs_vpcl_vs_vpcl_se + theme(legend.position = 'none'), filename = "./figures/lm_vs_vpcl_vs_vpcl_se.png", dpi = 800, width = 220, height = 100, units = "mm")
ggsave(plot = lm_vs_vpcl_vs_vpcl_se, filename = "./figures/lm_vs_vpcl_vs_vpcl_se_legend.svg", dpi = 800, width = 200, height = 80, units = "mm")

























#### Trying to define the geom_points Y values by making a a new column ####


all_vp<- dplyr:: bind_rows(lr_24, vpcl_24, vpcl_se_24, vpcl_se_gte)    %>%
  mutate(Sample_Type = recode(Sample_Type,
                              "VP" = "Lytic",
                              "Diff" = "Lysogenic"),
         Sample_Type = factor(Sample_Type, levels = c("Lytic", "Lysogenic")),
         Time_Range = factor(Time_Range, levels = c("T0_T3", "T0_T6",
                                                    "T0_T9", "T0_T12",
                                                    "T0_T17", "T0_T20",
                                                    "T0_T24")))%>%
  dplyr::filter(!is.na(Sample_Type)) %>%
  mutate(VP_SE = replace_na(VP_SE, 0)) %>% #substituting NA with O
  mutate(geom_point_position = ifelse(VP > 0, ((VP + VP_SE)/1e+6) + 0.2, (((VP/1e+6) - (2 *(VP_SE/1e+6))) - 0.2)))


lm_vs_vpcl_vs_vpcl_se_time <- ggplot(all_vp %>% dplyr::filter(Sample_Type != "VPC"),  
                                     aes(x = as.factor(Station_Number), y = VP / 1e+6, fill = VP_Method)) +
  # Bar plot
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
  
  # Error bars
  geom_errorbar(aes(ymin = (VP - VP_SE) / 1e+6, ymax = (VP + VP_SE) / 1e+6),
                position = position_dodge(width = 0.7), width = 0.25) +
  
  # Labels
  labs(
    x = NULL,
    y = NULL,
    fill = "Assay Duration"
  ) +
  
  # Manual fill colors for bars
  scale_fill_manual(values = c("LM_S" = "#DC2828", "VIPCAL" = "#003049", "VIPCAL_SE" = "#F6CD61")) +
  
  # New scale for fill
  new_scale_fill() +
  
  # Points for LM_S, nudged to the left
  geom_point(data = all_vp %>% dplyr::filter(Sample_Type != "VPC", VP_Method == "LR"),
             aes(x = as.factor(Station_Number), y = geom_point_position, fill = Time_Range),
             position = position_nudge(x = -0.24), shape = 21, size = 3, stroke = 0.5) +
  
  # Points for VIPCAL, no adjustment
  geom_point(data = all_vp %>% dplyr::filter(Sample_Type != "VPC", VP_Method == "VIPCAL"),
             aes(x = as.factor(Station_Number), y = geom_point_position, fill = Time_Range),
             position = position_nudge(x = 0), shape = 21, size = 3, stroke = 0.5) +
  
  # Points for VIPCAL_SE, nudged to the right
  geom_point(data = all_vp %>% dplyr::filter(Sample_Type != "VPC", VP_Method == "VIPCAL-SE"),
             aes(x = as.factor(Station_Number), y = geom_point_position, fill = Time_Range),
             position = position_nudge(x = 0.24), shape = 21, size = 3, stroke = 0.5) +
  
  # Scale for Time Range fill colors
  scale_fill_manual("Time Range", values = c("T0_T12" = "#E41A1C", "T0_T17" = "#377EB8", 
                                             "T0_T20" = "#4DAF4A", "T0_T24" = "#984EA3", 
                                             "T0_T3" = "#FF7F00", "T0_T6" = "#FFFF33", 
                                             "T0_T9" = "#A65628")) +
  
  
  # Additional plot settings
  scale_y_continuous(breaks = seq(-3, 6, 1), expand = expansion(add = c(0.5, 0.5))) +
  facet_wrap(~ Sample_Type, ncol = 2, scales = "fixed") +  
  theme_classic(base_size = 18) +
  theme(strip.background = element_blank())

lm_vs_vpcl_vs_vpcl_se_time

ggsave(plot = lm_vs_vpcl_vs_vpcl_se_time + theme(legend.position = 'none'), filename = "./figures/lm_vs_vpcl_vs_vpcl_se_time.png", dpi = 800, width = 240, height = 120, units = "mm")
ggsave(plot = lm_vs_vpcl_vs_vpcl_se_time, filename = "./figures/lm_vs_vpcl_vs_vpcl_se_time_legend.svg", dpi = 800, width = 240, height = 120, units = "mm")



