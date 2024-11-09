
# The aim fo the script is to extract viral counts from FCM data of NJ2020 viral production assays

#### 0.0 Setting up ####

work_dir<- paste0(getwd(), "/") # work_dir is an important variable for further calculations
project_title<- "nj2020_fcm" # project title
source("./scripts/4
       _nj2020_fcm_count_source.R")
set_up_vp_count()

#### 1.0 Importing metadata ####
metadata_processing("./data/metadata/nj2020_fcm_metadata.xlsx")


#### 2.0 FCS file shave been moved to `data/raw_data` ####

#if not, define origin and the following function will read `Sample_Name` variable
# from themetadata file to move the fcs files to `data/raw_data`

#origin<- "./data/raw_data"
#import_fcs(origin)


#### 3.0 Setting up gates ####
#This is already a part of the source script. Although, if anything needs to be chnaged, it can be done here
{
  polycut<- matrix(c(1.1, 2.0, 2.5, 3.2, 3.7, 3.7, 2.0, 0.6,
                     1.9, 1.5, 1.5, 1.8, 2.8, 3.7, 3.2, 2.75), nrow = 8, ncol=2)
  colnames(polycut)<- c("SSC-H", "FL1-H")
  bgate<- polygonGate(.gate = polycut, filterId = "Bacteria" )
  vgate<- rectangleGate(filterId= "Viruses", "SSC-H" = c(-0.1,1.35), "FL1-H" = c(-0.1, 1.7)) 
  #Different bacterial gates for viral vs bacterial stained samples. We need different gates cause sometimes the bacterial
  #and viral samples have a slight shift.
  HNA_Bacteria_b<- rectangleGate(filterId="HNA_Bacteria",
                                 "SSC-H"=c(0.6, 3.5), "FL1-H"=c(2.15, 3.5))
  LNA_Bacteria_b<- rectangleGate(filterId="LNA_Bacteria",
                                 "SSC-H"=c(1.0, 3.5), "FL1-H"=c(1.5, 2.15))
  
  HNA_Bacteria_v<- rectangleGate(filterId="HNA_Bacteria",
                                 "SSC-H"=c(0.6, 3.5), "FL1-H"=c(2.3, 3.5))
  LNA_Bacteria_v<- rectangleGate(filterId="LNA_Bacteria",
                                 "SSC-H"=c(1.0, 3.5), "FL1-H"=c(1.5, 2.3))
  
  #in case it doesn't, we could define the gate as following
  HNA_Bacteria_bv<- rectangleGate(filterId="HNA_Bacteria",
                                  "SSC-H"=c(0.6, 3.7), "FL1-H"=c(2.3, 3.6))
  LNA_Bacteria_bv<- rectangleGate(filterId="LNA_Bacteria",
                                  "SSC-H"=c(1.0, 3.7), "FL1-H"=c(1.5, 2.3))
  
  #same viral gates as we don't utilise the viral info from bacterial samples
  v1<- rectangleGate(filterId="V1", 
                     "SSC-H"= c(-0.1, 0.90), "FL1-H"= c(-0.1, 0.8)) 
  v2<- rectangleGate(filterId="V2", 
                     "SSC-H"= c(-0.1, 0.90), "FL1-H"= c(0.8, 1.25)) 
  v3<- rectangleGate(filterId="V3", 
                     "SSC-H"= c(-0.1, 1.3), "FL1-H"= c(1.25, 1.7))
  
  detectors<- c("FSC-H", "SSC-H", "FL1-H", "FL2-H", "FL3-H")
  
  translist_bv<- transformList(detectors, logTransform())
}

#### 4.0 Performing a small test ####
get_bv_stats(test = T)
get_bv_plots(test = T)

#examine the test files and remove

#### 5.0 Extracting FCM count ####
get_bv_plots()
get_bv_stats()

#### 6.0 Tris-EDTA (TE) control examination ####

#Combining counts and metadata
counts <- as.data.frame(read_csv("./results/nj2020_fcm/counts/nj2020_fcm_counts.csv")) %>%
  mutate(file_name = str_replace_all(file_name, "\\.1$", ""))

combine_metadata_counts()

#Separate the dataframe containing TE
TE<- counts_metadata[counts_metadata$Sample_Type == 'TE',]
plotly::ggplotly(ggplot(data = TE[TE$Staining_Protocol == 'Viruses',], aes(x = Sample_Name ,y = V1V2V3))+
                   
                   geom_point(shape = 1))

{
#Use this to identify off TEs. Perhaps, we could get rid of all the TEs that are above 2000. Ideally TEs should eb below 500, but the FCM was acting up.
TE<- TE%>% dplyr::filter(V1V2V3 <2000)
plotly::ggplotly(ggplot(data = TE[TE$Staining_Protocol == 'Viruses',], aes(x = Sample_Name ,y = V1V2V3))+
                   geom_point(shape = 1))
#Looking at this i could also get rid of values above 1000 as there aren't any consecutive TEs above 1000
TE<- TE %>% dplyr::filter( V1V2V3 <1000)
plotly::ggplotly(ggplot(data = TE[TE$Staining_Protocol == 'Viruses',], aes(x = Sample_Name ,y = V1V2V3))+
                   geom_point(shape = 1))

#Outlier analysis 1

z_scores <- scale(TE$V1V2V3)  #compute the z-score
outliers <- (TE$V1V2V3)[abs(z_scores) > 3]
outliers

#Looking at this i could also get rid of values above 900 as there aren't any consecutive TEs above 1000
TE<- TE %>% dplyr::filter( V1V2V3 <900)
plotly::ggplotly(ggplot(data = TE[TE$Staining_Protocol == 'Viruses',], aes(x = Sample_Name ,y = V1V2V3))+
                   geom_point(shape = 1))

#Outlier analysis 2
z_scores <- scale(TE$V1V2V3)  #compute the z-score
outliers <- (TE$V1V2V3)[abs(z_scores) > 3]
outliers

#Looking at this i could also get rid of values above 900 as there aren't any consecutive TEs above 1000
TE<- TE %>% dplyr::filter( V1V2V3 >31)
plotly::ggplotly(ggplot(data = TE[TE$Staining_Protocol == 'Viruses',], aes(x = Sample_Name ,y = V1V2V3))+
                   geom_point(shape = 1))
}

#Now we need to remove all the TEs we got rid of

calc_TE()


#Let's examine the output of this script.
plotly::ggplotly(ggplot(data = counts_metadata, aes(x = Sample_Name ,y = TE_Vi))+
                   geom_point(shape = 1))
#The TE values look good.Except for the 800 something values
#We can now adjust the TE values

adjust_TE()

#### 7.0 Visualization ####

nj2020<- read.csv("./results/nj2020_fcm/nj2020_fcm_corrected_counts.csv")

viral_count_overview_plots(nj2020)
bacterial_count_overview_plots(nj2020)


#### 8.0 Calculating viral production ####


#Importing nj2020 dataset
nj2020<- read.csv("./results/nj2020_fcm/nj2020_fcm_corrected_counts.csv")

#Checks before viralprod
try(vp_class_count_data(nj2020)) #failed

str(nj2020)

#Adjusting columns
nj2020 <- nj2020 %>%
  mutate(across(all_of(c("Timepoint", "Replicate")), ~ as.numeric(as.character(.))),
         across(all_of(c("Location", "Sample_Type")), as.character),
         across(all_of(c("Depth", "Station_Number")), ~as.integer(as.character(.))))

str(nj2020)

vp_class_count_data(nj2020) #passed

#Assigning the class
nj2020<- vp_class_count_data(nj2020)

# NJ2020 original bacterial and viral abundance

nj2020_abundance<- read.csv("./data/metadata/nj2020_original_abundances.csv")

#Running viralprod calculate function to extract viral production rate


try(viralprod::vp_end_to_end(nj2020,
                             original_abundances = nj2020_abundance,
                            output_dir = "./results/nj2020_viral_production/",
                            SR_calc = T,
                            BP_endpoint = T))
# 
# # Some checks
# 
# nj2020_vp_all<- read.csv("./results/nj2020_viral_production/vp_results_ALL.csv")
# 
# unique(nj2020_vp_all$VP_Method)
# unique(nj2020_vp_all$Station_Number)
# str(nj2020_vp_all)
# 
# nj2020_vp_all <- nj2020_vp_all %>%
#   mutate(Station_Number = as.numeric(Station_Number))
# str(nj2020_vp_all)
# 
# 
