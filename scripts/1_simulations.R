#Aim: To generate simulation dataset to comparatively assess viral production analyses methods.

####0.0 Setting up the environment####

# Installing packages
library(tidyverse)
library(devtools)
#install_github("mdhishamshaikh/ViralProduction_R") #viralprod package
library(viralprod)


####1.0 Creating a simulation dataframe####

# #Set number of dataframes you'd like to create
simu_length<- 1000
{
  set.seed(2023) #Setting seed for reproducibility
  simu_df<- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(simu_df)<- c("Timepoint", "c_Viruses", "Replicate", "Station_Number", "Sample_Type")
  simu_df<- simu_df %>%  mutate(across(!Sample_Type, as.numeric)) %>%
    mutate(across(Sample_Type, as.character))
  
  
  for (df_no in 1:simu_length){
    
    A<- data.frame(Count = runif(n=12, min = 2, max = 3)) 
    A$Timepoint<- rep(1:6,2)
    A$SD_percent<- runif(n = 12, min = 0, max = 100) #Up to 100% error 
    A$SD<- (A$Count*A$SD_percent)/100
    
    B<- apply(A, 1, function(x){x[1]+x[4]*abs(scale(rnorm(3)))} ) #Ensuring positive values
    B<- as.data.frame(B)
    colnames(B)<- rep(1:6,2)
    B<- B%>% pivot_longer(cols = everything(), names_to = 'Timepoint',
                          values_to = 'c_Viruses') %>%
      mutate(c_Viruses = as.numeric(c_Viruses))
    B<- B%>%
      arrange(Timepoint) %>%
      mutate(Replicate = rep(1:6, 6))%>%
      mutate(Sample_Type = rep(c(rep('VP',3), rep('VPC',3)),6))
    
    B$Station_Number = df_no
    
    
    B<- B%>% mutate(across(!Sample_Type, as.numeric))
    
    try(simu_df<- simu_df %>%full_join(B))
    rm(A,B, df_no)
  }
}

#Adding additional column to fit viralprod tool
cols_to_add<- c("Location", "Depth")

for(col in cols_to_add) {
  simu_df[[col]] <- 1
}

#Writing csv
write.csv(simu_df, file = "data/simulation_df.csv", row.names = F)



####2.0 Running viralprod on simulated dataset"

#Importing simulation dataset
simu_df<- read.csv("./data/simulation_df.csv")

#Checks before viralprod
try(vp_class_count_data(simu_df)) #failed

str(simu_df)

#Adjusting columns
simu_df <- simu_df %>%
  mutate(across(all_of(c("Timepoint", "Replicate")), ~ as.numeric(as.character(.))),
         across(all_of(c("Location", "Sample_Type")), as.character),
         across(all_of(c("Depth", "Station_Number")), ~as.integer(as.character(.))))

str(simu_df)

vp_class_count_data(simu_df) #passed

#Assigning the class
simu_df<- vp_class_count_data(simu_df)


#Running viralprod calculate function to extract viral production rate

try(viralprod::vp_calculate(simu_df,
        output_dir = "./results/simulation_viral_production/",
        SR_calc = F,
        bp_endpoint = F))


simu_vp_all<- read.csv("./results/simulation_viral_production/vp_results_ALL.csv")

unique(simu_vp_all$VP_Method)
str(simu_vp_all)

simu_vp_all <- simu_vp_all %>%
  mutate(Station_Number = as.numeric(Station_Number))



