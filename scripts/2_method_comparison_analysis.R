#Aim: Perform statistics and genertae plots for comparison.

#### 0.0 Setting up
set.seed(2023)
library(tidyverse)
library(ggsci)
library(dunn.test)


simu_vp_all<- read.csv("./results/simulation_viral_production/vp_results_ALL.csv")

unique(simu_vp_all$VP_Method)
str(simu_vp_all)

simu_vp_all <- simu_vp_all %>%
  mutate(Station_Number = as.numeric(Station_Number))


#### 1.0 LM methods ####

#### 1.1 LM methods - dataframe ####

#Filtering only the linear regression model approaches used in literature
simu_vp_lm <- simu_vp_all %>%
  filter(VP_Method %in% c("LM_SR_AVG", "LM_AR", "LM_ALLPOINTS"),
         Time_Range == "T1_T6",
         Sample_Type != "VPC")

#making the labels simpler
replacements <- c("LM_SR_AVG" = "LM_1", "LM_AR" = "LM_2", "LM_ALLPOINTS" = "LM_3",
                  "VP" = "Lytic", "Diff" = "Lysogenic")
simu_vp_lm <- simu_vp_lm %>%
  mutate(
    VP_Method = case_when(
      VP_Method %in% names(replacements) ~ replacements[VP_Method],
      TRUE ~ VP_Method
    ),
    Sample_Type = case_when(
      Sample_Type %in% names(replacements) ~ replacements[Sample_Type],
      TRUE ~ Sample_Type
    )
  ) %>%
  mutate(
    Sample_Type = factor(Sample_Type, levels = c("Lytic", "Lysogenic"))
  )

rm(replacements)

#### 1.2 LM methods - Statistical tests ####

kruskal.test(VP ~ VP_Method, data = simu_vp_lm) # no difference
kruskal.test(VP_SE ~ VP_Method, data = simu_vp_lm) # significantly different.
dunn_test_results <- dunn.test(simu_vp_lm$VP_SE, simu_vp_lm$VP_Method, method="bonferroni")
print(dunn_test_results)


#### 1.3 LM methods - basic plots ####

ggplot(simu_vp_lm, aes(x = VP_Method, y = VP, fill = Sample_Type)) +
  geom_split_violin(nudge = 0.02) +
  scale_fill_manual(values = c( "Lytic" = "#d43028", "Lysogenic" = "#4d778b")) +
  theme_classic()
ggplot(simu_vp_lm, aes(x = VP_Method, y = VP_SE, fill = Sample_Type)) +
  geom_split_violin(nudge = 0.02)+
  scale_fill_manual(values = c( "Lytic" = "#d43028", "Lysogenic" = "#4d778b")) +
  theme_classic()


#### 2.0 LM vs VIPCAL ####

#### 2.1 LM vs VIPCAL - data frame ####

#Filtering the linear regression model and VIPCAL approach
simu_vp_lm_vs_vpcl <- simu_vp_all %>%
  filter(VP_Method %in% c("LM_SR_AVG", "VPCL_AR_DIFF"),
         Time_Range == "T1_T6",
         Sample_Type != "VPC")

#making the labels simpler
replacements <- c("LM_SR_AVG" = "LM", "VPCL_AR_DIFF" = "VIPCAL",
                  "VP" = "Lytic", "Diff" = "Lysogenic")
simu_vp_lm_vs_vpcl <- simu_vp_lm_vs_vpcl %>%
  mutate(
    VP_Method = case_when(
      VP_Method %in% names(replacements) ~ replacements[VP_Method],
      TRUE ~ VP_Method
    ),
    Sample_Type = case_when(
      Sample_Type %in% names(replacements) ~ replacements[Sample_Type],
      TRUE ~ Sample_Type
    )
  ) %>%
  mutate(
    Sample_Type = factor(Sample_Type, levels = c("Lytic", "Lysogenic"))
  )

rm(replacements) 

# Percenatge comparisons between the
comp_lm_vs_vpcl <- simu_vp_lm_vs_vpcl %>%
  select(!c(abs_VP, VP_SE, VP_R_Squared)) %>%
  pivot_wider(names_from = VP_Method, values_from = VP) %>%
  mutate(Comparison = case_when(
    LM == VIPCAL ~ "Same",
    LM > VIPCAL ~ "LM_Higher",
    LM < VIPCAL ~ "VIPCAL_Higher",
    TRUE ~ "Data Missing"  # Catching cases where data might be missing or NA
  ),
  ZeroStatus = case_when(
    LM == 0 & VIPCAL == 0 ~ "Both_Zero",
    LM == 0 & VIPCAL != 0 ~ "LM_Zero",
    LM != 0 & VIPCAL == 0 ~ "VIPCAL_Zero",
    TRUE ~ "None_Zero"
  ))

#### 2.2 LM vs VIPCAL - statistical analysis ####

kruskal.test(VP ~ VP_Method, data = simu_vp_lm_vs_vpcl) # significant difference
kruskal.test(VP ~ VP_Method, data = simu_vp_lm_vs_vpcl %>%
               filter(Sample_Type == 'Lytic')) # significant difference
kruskal.test(VP ~ VP_Method, data = simu_vp_lm_vs_vpcl %>%
               filter(Sample_Type == 'Lysogenic')) # significant difference

# Comparing the distribution 
comp_lm_vs_vpcl %>% 
  filter(Sample_Type == 'Lytic') %>%  #also tried with lysogenic
  count(Comparison) %>%
  mutate(Percentage = n / sum(n) * 100)

# VIPCAL-SE provides zero estimations more readily
comp_lm_vs_vpcl %>% 
  filter(Sample_Type == 'Lytic') %>%
  count(ZeroStatus) %>%
  mutate(Percentage = n / sum(n) * 100)


#### 3.0 VIPCAL vs VIPCAL_SE ####

#### 3.1 VIPCAL vs VIPCAL_SE - data frame ####

#Filtering VIPCAL and VIPCAL_SE approaches
simu_vp_vpcl_vs_vpclse <- simu_vp_all %>%
  filter(VP_Method %in% c("VPCL_AR_DIFF", "VPCL_AR_DIFF_SE"),
         Time_Range == "T1_T6",
         Sample_Type != "VPC")

#making the labels simpler
replacements <- c("VPCL_AR_DIFF" = "VIPCAL", "VPCL_AR_DIFF_SE" = "VIPCAL_SE",
                  "VP" = "Lytic", "Diff" = "Lysogenic")
simu_vp_vpcl_vs_vpclse <- simu_vp_vpcl_vs_vpclse %>%
  mutate(
    VP_Method = case_when(
      VP_Method %in% names(replacements) ~ replacements[VP_Method],
      TRUE ~ VP_Method
    ),
    Sample_Type = case_when(
      Sample_Type %in% names(replacements) ~ replacements[Sample_Type],
      TRUE ~ Sample_Type
    )
  ) %>%
  mutate(
    Sample_Type = factor(Sample_Type, levels = c("Lytic", "Lysogenic"))
  )

rm(replacements)

comp_vpcl_vs_vpclse <- simu_vp_vpcl_vs_vpclse %>%
  select(!c(abs_VP, VP_SE, VP_R_Squared)) %>%
  pivot_wider(names_from = VP_Method, values_from = VP) %>%
    mutate(Comparison = case_when(
      VIPCAL == VIPCAL_SE ~ "Same",
      VIPCAL > VIPCAL_SE ~ "VIPCAL_Higher",
      VIPCAL < VIPCAL_SE ~ "VIPCAL_SE_Higher",
    TRUE ~ "Data Missing"  # Catching cases where data might be missing or NA
  ),
  ZeroStatus = case_when(
    VIPCAL == 0 & VIPCAL_SE == 0 ~ "Both_Zero",
    VIPCAL == 0 & VIPCAL_SE != 0 ~ "VIPCAL_Zero",
    VIPCAL != 0 & VIPCAL_SE == 0 ~ "VIPCAL_SE_Zero",
    TRUE ~ "None_Zero"
  ))

#### 3.2 VIPCAL vs VIPCAL_SE - statistical analysis ####

kruskal.test(VP ~ VP_Method, data = simu_vp_vpcl_vs_vpclse) # significant difference
kruskal.test(VP ~ VP_Method, data = simu_vp_vpcl_vs_vpclse %>%
               filter(Sample_Type == 'Lytic')) # significant difference
kruskal.test(VP ~ VP_Method, data = simu_vp_vpcl_vs_vpclse %>%
               filter(Sample_Type == 'Lysogenic')) # significant difference


# Comparing the distribution 
comp_vpcl_vs_vpclse %>% 
  filter(Sample_Type == 'Lytic') %>%
  count(Comparison) %>%
  mutate(Percentage = n / sum(n) * 100)

# VIPCAL-SE provides zero estimations more readily - lytic
comp_vpcl_vs_vpclse %>% 
  filter(Sample_Type == 'Lytic') %>%
  count(ZeroStatus) %>%
  mutate(Percentage = n / sum(n) * 100)


#Now do the same for lysogenic
