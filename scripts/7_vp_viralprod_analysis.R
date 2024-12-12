# The aim of the script is to run viralprod R package on DNS (coastal Dutch North Sea) data

# 0. Setting up ####
source("./scripts/0_source.R")


# 1.0 Importing Abundance data ####

abundance <- read.csv("./data/metadata/bacterial_viral_abundances_abiotic.csv")


# 2.0 Running viralprod #####

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
              output_dir = './results/viralprod')

