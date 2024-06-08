
# The aim fo the script is to extract viral counts from FCM data of NJ2020 viral production assays

#### 0.0 Setting up ####

work_dir<- paste0(getwd(), "/") # work_dir is an important variable for further calculations
project_title<- "nj2020_fcm" # project title
source("./scripts/3_nj2020_fcm_count_source.R")
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
