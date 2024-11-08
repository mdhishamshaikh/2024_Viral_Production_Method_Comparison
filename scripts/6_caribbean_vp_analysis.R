#### 0.0 Setting up ####
library(tidyverse)
library(viralprod)

#### 1.0 Calculating viral production ####


caribbean<- readxl::read_xlsx("./data/caribbean_vp_assays/caribbean_vp_assay.xlsx")
caribbean <- caribbean %>% select(-c(c_HNA, c_LNA, c_V1, c_V2, c_V3))


#Checks before viralprod
try(vp_class_count_data(caribbean)) #failed

# I have to change the station numbers here to numeric. So I will com up with a key 

caribbean <- caribbean %>%
  mutate(Station_Number = recode(Station_Number, !!!c("13-1-15" = "1", 
                                                      "30-1-16" = "2", 
                                                      "40-1-9" = "3", 
                                                      "52-1-13" = "4", 
                                                      "62-1-20" = "5")))

str(caribbean)

#Adjusting columns
caribbean <- caribbean %>%
  mutate(across(all_of(c("Timepoint", "Replicate")), ~ as.numeric(as.character(.))),
         across(all_of(c("Location", "Sample_Type")), as.character),
         across(all_of(c("Depth", "Station_Number")), ~as.integer(as.character(.))))
caribbean$Depth <- 1
str(caribbean)

vp_class_count_data(caribbean) #passed

#Assigning the class
caribbean<- vp_class_count_data(caribbean)

# caribbean original bacterial and viral abundance

caribbean_abundance<- read.csv("./data/metadata/caribbean_original_abundances.csv")

#Running viralprod calculate function to extract viral production rate


try(viralprod::vp_end_to_end(caribbean ,
                             original_abundances = caribbean_abundance,
                             output_dir = "./results/caribbean_viral_production/",
                             SR_calc = T,
                             BP_endpoint = T))


caribbean_vp_all<- read.csv("./results/caribbean_viral_production/vp_results_ALL.csv")

unique(caribbean_vp_all$VP_Method)
unique(caribbean_vp_all$Station_Number)
str(caribbean_vp_all)
summarise(caribbean_vp_all)

caribbean_vp_all <- caribbean_vp_all %>%
  mutate(Station_Number = as.numeric(Station_Number))
str(caribbean_vp_all)

