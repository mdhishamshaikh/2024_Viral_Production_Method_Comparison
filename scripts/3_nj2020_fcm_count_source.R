# 1.0 Set up functions ####

{
  #List of packages to install
  packages_to_load<- c("tidyverse", 
                       "flowWorkspace",
                       "flowCore",
                       "scales",
                       "readxl",
                       "ggcyto",
                       "ggsci",
                       "svglite")
  
  {
    if (!requireNamespace("BiocManager"))
      install.packages("BiocManager")
  }
  
  #Checking if packages are already present. If absent, then installing packages from BiocManager 
  for (pack in packages_to_load){
    if(!requireNamespace(pack))
      BiocManager::install(pack, force = T)
  }
  
  #Loading libraries  
  for (pack in packages_to_load){
    library(pack, character.only = T)
  }
  
  
  
  rm(pack)
  print(paste0("Package loaded: ", packages_to_load))
  rm(packages_to_load)
}

set_up_vp_count<- function(){
  { #create directories
    dir_to_create<- c("./data", "./results", "./data/raw_data", 
                      "./data/metadata")
    
    for (dir in 1: length(dir_to_create)){
      
      if("FALSE" %in% file.exists(dir_to_create[dir])){
        dir.create(dir_to_create[dir])
        
        print(paste0("The ", dir_to_create[dir], " directory was created"))
      }
      
    }
    
    rm(dir_to_create, dir)
    
  }
  
  
}

metadata_processing<- function(file_path, extension = ".xlsx", sheet =NULL, format = "%y%m%d", ...){
  if (extension == ".xlsx"){
    metadata<- read_excel(file_path, sheet)
  } else if (extension == ".csv") {
    metadata<- read.csv(file_path)
  } else {
    print("File format not supported. Please convert the metadata file to .xlsx or .csv")
  }
  
  metadata<- metadata[, c('Sample_Name', 'Staining_Protocol', 'Date_Measurement', 
                          'Location', 'Expt_Date', 'Expt_No', 'Depth', 
                          'Sample_Type', 'Timepoint', 'Replicate', 'Acquisition_Duration',
                          'Dilution', 'Flowrate')] %>%
    mutate(Expt_Date = as.Date(as.character(Expt_Date), format)) %>% #converting Expt_Date and Date_Measurement into Date format
    mutate(Date_Measurement= as.Date(as.character(Date_Measurement), format))
  .GlobalEnv$metadata<- metadata
  
  write.csv(metadata, paste0("data/metadata/",project_title,"_metadata.csv"), row.names=F)
  print("Metadata processed and stored under `data/metadata`")
  return(metadata)
  
  
}

import_fcs<- function(fcs_dir, ...){
  
  print("It'll take some time before all the files are transferred")
  file.copy(from = paste0(fcs_dir, metadata$Sample_Name),
            to = "./data/raw_data")
  
  if( "FALSE" %in% (metadata$Sample_Name %in% list.files("./data/raw_data"))) { #files are missing in the raw_data folder
    print("Missing files. Check `missing_fcs_file`")
    print(setdiff(metadata$Sample_Name,list.files("./data/raw_data"))) #what to do with the missing files
    missing_fcs_file<- setdiff(metadata$Sample_Name,list.files("./data/raw_data"))
  } else {
    print("All files were transferred.")
  }
  fcs_data<- paste0(work_dir,"data/raw_data/", 
                    list.files(paste0(work_dir, "data/raw_data/")))
}



# 2.0 Gating Functions ####
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

bacterial_gate<- function(gate = "same", df = metadata, fcs_file, ...){
  
  if (gate == "different"){
    if(df[df$Sample_Name == fcs_file,]$Staining_Protocol == "Bacteria") { #MAKE THIS INTO A FUNCTION
      print("Bacteria")
      .GlobalEnv$HNA_Bacteria<- HNA_Bacteria_b
      .GlobalEnv$LNA_Bacteria<- LNA_Bacteria_b
    } else if (df[df$Sample_Name == fcs_file,]$Staining_Protocol == "Viruses"){
      print("Viruses")
      .GlobalEnv$HNA_Bacteria<- HNA_Bacteria_v
      .GlobalEnv$LNA_Bacteria<- LNA_Bacteria_v
      
    } else {
      print("Staining protocol niether Bacteria or Viruses")
    }
  } else if (gate == "same"){
    .GlobalEnv$HNA_Bacteria<- HNA_Bacteria_bv
    .GlobalEnv$LNA_Bacteria<- LNA_Bacteria_bv
  }
  
  
}

read_transform_fs_bv <- function(x){ #function to read and transform fcs files
  flowCore::read.flowSet(c(x, x)) %>%
    flowCore::transform(translist_bv)
  
}

gatingset_bv_stats<- function(flowset){ #flowset here is already transformed, cleaned, and compensated
  gsbv_fs<- flowWorkspace::GatingSet(flowset)
  gs_pop_add(gsbv_fs, bgate, parent="root")
  gs_pop_add(gsbv_fs, HNA_Bacteria, parent= "Bacteria")
  gs_pop_add(gsbv_fs, LNA_Bacteria, parent= "Bacteria")
  gs_pop_add(gsbv_fs, vgate, parent= "root")
  gs_pop_add(gsbv_fs, v1, parent="Viruses")
  gs_pop_add(gsbv_fs, v2, parent="Viruses")
  gs_pop_add(gsbv_fs, v3, parent="Viruses")
  recompute(gsbv_fs)
  stats<- gs_pop_get_stats(gsbv_fs[[2]], c("root", "Bacteria", "HNA_Bacteria",
                                           "LNA_Bacteria", "Viruses",
                                           "V1", "V2", "V3"), "count")
  write.table(stats, file = csv_file_name, sep = ",", append = T, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  
}

get_bv_stats <- function(df = metadata, gate = "same", write_csv = TRUE, test = FALSE, ...) {
  # Base directory as per test
  base_dir <- if (test) {
    df <- df[1:10,]  # Subsetting the dataframe for test
    "./results/test/test_stats/"
  } else {
    "./results/"
  }
  
  # Checking for project_title and setting directory paths
  if ("project_title" %in% ls(.GlobalEnv)) {
    if (test) {
      results_dir <- base_dir  # Using test directory if test condition
      dir.create(results_dir, recursive = TRUE)
    } else {
      # Creating directory with project title under results if not a test
      results_dir <- paste0(base_dir, project_title, "/stats/")
      dir.create(results_dir, recursive = TRUE)
    }
  } else {
    results_dir <- paste0(base_dir, "stats/")  # Default stats directory
    dir.create(results_dir, recursive = TRUE)
  }
  
  # Setting the file name for the CSV file
  .GlobalEnv$csv_file_name <- paste0(results_dir, ifelse(test, paste0(project_title, "_test"),  project_title), "_counts.csv")
  
  # Creating a CSV file to store all the stats
  if (write_csv) {
    x <- data.frame(file_name = character(), pop = character(), count = numeric())
    write.table(x, file = csv_file_name, sep = ",", col.names = TRUE, row.names = FALSE)
    print(paste0(csv_file_name, " created"))
  }
  
  # Processing samples
  t <- 0  # Counter for tracking progress
  for (i in 1:length(df$Sample_Name)) {
    fcs_file <- df$Sample_Name[i]
    fcs_data <- paste0(work_dir, "data/raw_data/", fcs_file)
    
    bacterial_gate(fcs_file = fcs_file)
    print(df$Sample_Name[i])
    
    try({
      read_transform_fs_bv(fcs_data) %>% 
        gatingset_bv_stats()
    })
    t <- t + 1
    print(paste0(t, "/", length(df$Sample_Name)))
    
    print(HNA_Bacteria)
    print(LNA_Bacteria)
  }
}

gatingset_bv_plots <- function(flowset, position, fcs_file, bins = 600, ...){ 
  library(flowWorkspace)
  library(ggcyto)
  
  gsbv_fs <- flowWorkspace::GatingSet(flowset)
  # Adding population gates 
  gs_pop_add(gsbv_fs, bgate, parent="root")
  gs_pop_add(gsbv_fs, HNA_Bacteria, parent="Bacteria")
  gs_pop_add(gsbv_fs, LNA_Bacteria, parent="Bacteria")
  gs_pop_add(gsbv_fs, vgate, parent="root")
  gs_pop_add(gsbv_fs, v1, parent="Viruses")
  gs_pop_add(gsbv_fs, v2, parent="Viruses")
  gs_pop_add(gsbv_fs, v3, parent="Viruses")
  recompute(gsbv_fs)
  
  # Configuring the plot
  p <- ggcyto::ggcyto(gsbv_fs[[2]], aes(x = `SSC-H`, y = `FL1-H`), subset = "root") +
    geom_hex(bins = bins) +  
    theme_bw() +
    labs(title= paste0(position, "   ", fcs_file), x = "Side scatter (a.u.)", y = "Green Fluorescence (a.u.)") +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "white")) +
    ggcyto_par_set(limits = list(x = c(-0.75,3.8), y = c(0, 3.6)))
  
  p <- p + geom_gate("Viruses", colour = "black", size = 1) + geom_stats("Viruses", type = c("gate_name", "count"), adjust = c(1.25, 0), colour = "black") +
    geom_gate("Bacteria", colour = "black", size = 1) + geom_stats("Bacteria", type = c("gate_name", "count"), adjust = c(-0, 1), colour = "black") +
    geom_gate(c("HNA_Bacteria", "LNA_Bacteria"), colour = "red", size = 0.8) + geom_stats(c("HNA_Bacteria", "LNA_Bacteria"), type = c("gate_name", "count", "percent"), adjust = c(-0.25, 0.75), colour = "red") +
    geom_gate(c("V1", "V2", "V3"), colour = "blue", size = 0.8) + geom_stats(c("V1", "V2", "V3"), type = c("gate_name", "count", "percent"), adjust = c(-0.25, 0.75), colour = "blue")
  
  print(p)
}

get_bv_plots <- function(df = metadata, gate = "same", write_pdf = TRUE, test = FALSE, ...) {
  # Define the base directory based on the test condition
  base_dir <- if (test) {
    df <- df[1:10,]  # Subsetting the dataframe for test
    "./results/test/test_plots/"
  } else {
    "./results/"
  }
  
  # Set up the directory paths based on project_title presence
  if ("project_title" %in% ls(.GlobalEnv)) {
    if (test) {
      results_dir <- base_dir  # Using test directory if test condition
      dir.create(results_dir, recursive = TRUE)
    } else {
      # Creating a directory with the project title under results if not a test
      results_dir <- paste0(base_dir, project_title, "/plots/")
      dir.create(results_dir, recursive = TRUE)
    }
  } else {
    results_dir <- paste0(base_dir, "plots/")  # Default plots directory
    dir.create(results_dir, recursive = TRUE)
  }
  
  # Define the file name for the PDF based on conditions
  pdf_file_name <- paste0(results_dir, ifelse(test, paste0(project_title, "_test"), project_title), "_plots.pdf")
  
  # Create a PDF file to store all the plots
  if (write_pdf) {
    pdf(pdf_file_name, onefile = TRUE)
    print(paste0(pdf_file_name, " created"))
  }
  
  t <- 0 # Initialize counter to 0
  p <- list() # Initialize an empty list for storing plot objects
  
  for(position in 1:length(df$Sample_Name)) {
    fcs_file <- df$Sample_Name[position]
    fcs_data <- paste0(work_dir, "data/raw_data/", fcs_file)
    
    # Assume bacterial_gate is a function that prepares or filters data
    bacterial_gate(fcs_file = fcs_file)
    print(df$Sample_Name[position])
    
    # Generate and store plots using a try block to catch errors
    p[[position]] <- try({
      read_transform_fs_bv(fcs_data) %>% 
        gatingset_bv_plots(position, fcs_file)  # Call plotting function with position
    })
    
    t <- t + 1
    print(paste0(t, "/", length(df$Sample_Name)))
    
    # Assume these are diagnostic prints to output information
    print(HNA_Bacteria)
    print(LNA_Bacteria)
  }
  
  graphics.off()  # Close the graphics device once plotting is complete
}


# 3.0 Wrangling ####
combine_metadata_counts <- function(metadata_df = metadata, counts_df = counts){
  
  counts<- pivot_wider(counts_df, names_from = "pop", values_from = "count")
  #dim(counts)
  colnames(counts)[colnames(counts) == 'root'] <- "Total"
  colnames(counts)[colnames(counts) == 'file_name'] <- "Sample_Name"
  #head(counts)
  
  
  if (length(setdiff(counts$Sample_Name, metadata_df$Sample_Name)) > 0) {
    print("There are samples that were counted but doesn't have associated metadata")
  } else {
    print("All samples processed have associated metadata")
  }
  
  {
    counts_metadata<- base::merge(counts, metadata_df, by = "Sample_Name") 
    #counts_metadata<- dplyr::mutate(counts_metadata,Expt_Date = as.Date(as.character(Expt_Date), format= "%y%m%d"))
    #counts_metadata<- dplyr::mutate(counts_metadata,Date_Measurement = as.Date(as.character(Date_Measurement), format= "%y%m%d"))
    #flowrate<- dplyr::mutate(flowrate, Date_Measurement = as.Date(as.character(Date_Measurement), format= "%y%m%d"))
    #we could have done flowrate per date of measurement, but because one of our sample is from 2020, the flowrate wasn't noted. So let's go with average.
    # dim(counts_metadata)
    # head(counts_metadata)
    #counts_metadata[,'Flowrate'] <- mean(flowrate$Flow_rate)
  }
  
  #Adding Boolean columns for bacterial and viral total counts
  {
    counts_metadata$HNALNA <- counts_metadata$HNA_Bacteria + counts_metadata$LNA_Bacteria
    counts_metadata$V1V2V3 <- counts_metadata$V1 + counts_metadata$V2 + counts_metadata$V3
    
    #We can perform a simple linear regression to see if we could use Boolean 
    # addition of HNA/LNA and V1/V2/V3 as replacements for total bacterial and viral 
    # counts. Alternatively, these could have been added as Boolean gates during processing
    # scatter.smooth(counts_metadata$Bacteria, counts_metadata$HNALNA, main = "Bacterial Counts")
    #  summary(lm(counts_metadata$Bacteria ~ counts_metadata$HNALNA))  #Gives an R-square value of 1.0
    #Perhaps the next few lines aren't as important. 
    # plot(fitted(lm(counts_metadata$Bacteria ~ counts_metadata$HNALNA)), resid(lm(counts_metadata$Bacteria ~ counts_metadata$HNALNA)))
    # abline(0,0)
    # qqnorm(resid(lm(counts_metadata$Bacteria ~ counts_metadata$HNALNA)))
    # qqline(resid(lm(counts_metadata$Bacteria ~ counts_metadata$HNALNA)))
    
    # scatter.smooth(counts_metadata$Viruses, counts_metadata$V1V2V3, main = "Viral Counts")
    #  summary(lm(counts_metadata$Viruses ~ counts_metadata$V1V2V3))  #Gives an R-square value of 0.999
    
  }
  .GlobalEnv$counts_metadata<- counts_metadata
  return(counts_metadata)
  
}

calc_TE<- function(df = counts_metadata, TE_df = TE ){ 
  
  counts_metadata<- df
  TE<- TE_df
  counts_metadata<- rbind(counts_metadata%>% dplyr::filter(Sample_Type != 'TE'), TE)  %>%
    arrange(Sample_Name)
  
  counts_metadata[,'TE_Ba']<- NA
  counts_metadata[,'TE_Vi']<- NA
  
  #We can now run TE calculation script
  for (name in counts_metadata$Sample_Name){ #Here I added TE value for both viruses and bacteria
    if (counts_metadata[counts_metadata$Sample_Name == name,]$Staining_Protocol == 'Viruses'
        #selecting rows with the specified file name that also had viral staining protocol
    ) {
      if (counts_metadata[counts_metadata$Sample_Name == name,]$Sample_Type  == 'TE') {
        print("TE") #we don't want this. so we move on and look for a count file
      } else if (counts_metadata[counts_metadata$Sample_Name == name,]$Sample_Type  != 'TE') {
        if (counts_metadata[which(counts_metadata$Sample_Name == name) + c(-1), ]$Sample_Type  == 'TE'
            #looking to see if there is a TE above this count file. if yes, we record it as 'a'
        ){
          a<- counts_metadata[which(counts_metadata$Sample_Name == name) + c(-1), ]$HNALNA #HNALNA and not V1V2V3 because we might use viral samples for bacterial measurements in VP assay
          if (counts_metadata[which(counts_metadata$Sample_Name == name) + c(-2), ]$Sample_Type  == 'TE'
              #here we see if there is a TE above our first TE. we record this as 'b'
          ){
            b<- counts_metadata[which(counts_metadata$Sample_Name == name) + c(-2), ]$HNALNA
          }
        }
        print(paste(name, a, b, mean(c(a,b)))) #too see what the output is
        counts_metadata[counts_metadata$Sample_Name == name,]$TE_Ba <- mean(c(a,b)) ##We are recording TE for bacteria in viral samples
        #adding the output to TE_value column
      }
      
      if (counts_metadata[counts_metadata$Sample_Name == name,]$Sample_Type  == 'TE') {
        print("TE") #we don't want this. so we move on and look for a count file
      } else if (counts_metadata[counts_metadata$Sample_Name == name,]$Sample_Type  != 'TE') {
        if (counts_metadata[which(counts_metadata$Sample_Name == name) + c(-1), ]$Sample_Type  == 'TE'
            #looking to see if there is a TE above this count file. if yes, we record it as 'a'
        ){
          c<- counts_metadata[which(counts_metadata$Sample_Name == name) + c(-1), ]$V1V2V3
          if (counts_metadata[which(counts_metadata$Sample_Name == name) + c(-2), ]$Sample_Type  == 'TE'
              #here we see if there is a TE above our first TE. we record this as 'b'
          ){
            d<- counts_metadata[which(counts_metadata$Sample_Name == name) + c(-2), ]$V1V2V3
          }
        }
        print(paste(name, a, b, mean(c(a,b)))) #too see what the output is
        counts_metadata[counts_metadata$Sample_Name == name,]$TE_Vi <- mean(c(c,d))
      }
    } else { #THIS PART WILL NOT RUN AS ALL MY SAMPLES ARE VIRAL 
      if (counts_metadata[counts_metadata$Sample_Name == name,]$Staining_Protocol == 'Bacteria') {
        if (counts_metadata[counts_metadata$Sample_Name == name,]$Sample_Type  == 'TE') {
          print("yes")
        } else if (counts_metadata[counts_metadata$Sample_Name == name,]$Sample_Type != 'TE') {
          if (counts_metadata[which(counts_metadata$Sample_Name == name) + c(-1), ]$Sample_Type  == 'TE'){
            a<- counts_metadata[which(counts_metadata$Sample_Name == name) + c(-1), ]$HNALNA
            if (counts_metadata[which(counts_metadata$Sample_Name == name) + c(-2), ]$Sample_Type  == 'TE'){
              b<- counts_metadata[which(counts_metadata$Sample_Name == name) + c(-2), ]$HNALNA
            }
          }
          print(paste(name, a, b, mean(c(a,b))))
          counts_metadata[counts_metadata$Sample_Name == name,]$TE_Ba <- mean(c(a,b)) ##We are recording TE for bacteria in bacterial samples. There is no need to do it for viruses
        }
      }
    }
  }
  .GlobalEnv$counts_metadata<- counts_metadata
  return(counts_metadata)
}

adjust_TE<- function(counts_metadata_df = counts_metadata, write_csv = T){
  
  counts_metadata<- counts_metadata_df
  for (cols in c( "c_Bacteria", "c_HNA", "c_LNA", "c_Viruses", "c_V1", "c_V2", "c_V3")){
    counts_metadata[, cols] <- NA
  }
  
  {
    counts_metadata$c_Viruses<- with(counts_metadata, ((V1V2V3-((V1V2V3/V1V2V3)*TE_Vi))*Dilution*60*1000)/(Flowrate*Acquisition_Duration))
    counts_metadata$c_V1<- with(counts_metadata, ((V1-((V1/V1V2V3)*TE_Vi))*Dilution*60*1000)/(Flowrate*Acquisition_Duration))
    counts_metadata$c_V2<- with(counts_metadata, ((V2-((V2/V1V2V3)*TE_Vi))*Dilution*60*1000)/(Flowrate*Acquisition_Duration))
    counts_metadata$c_V3<- with(counts_metadata, ((V3-((V3/V1V2V3)*TE_Vi))*Dilution*60*1000)/(Flowrate*Acquisition_Duration))
    counts_metadata$c_Bacteria<- with(counts_metadata, ((HNALNA-((HNALNA/HNALNA)*TE_Ba))*Dilution*60*1000)/(Flowrate*Acquisition_Duration))
    counts_metadata$c_HNA<- with(counts_metadata, ((HNA_Bacteria-((HNA_Bacteria/HNALNA)*TE_Ba))*Dilution*60*1000)/(Flowrate*Acquisition_Duration)) 
    counts_metadata$c_LNA<- with(counts_metadata, ((LNA_Bacteria-((LNA_Bacteria/HNALNA)*TE_Ba))*Dilution*60*1000)/(Flowrate*Acquisition_Duration))
    counts_metadata$VBR<- with(counts_metadata, c_Viruses/c_Bacteria) #calculated with bacteria from viral samples. Did an LR and the R2 was 0.95
    counts_metadata$HNAperLNA<- with(counts_metadata, c_HNA/c_LNA)
    #scatter.smooth(S22[S22$Staining_Protocol == 'Bacteria',]$c_Bacteria ~ S22[S22$Staining_Protocol == 'Viruses',]$c_Bacteria )
    #summary(lm(S22[S22$Staining_Protocol == 'Bacteria',]$c_Bacteria ~ S22[S22$Staining_Protocol == 'Viruses',]$c_Bacteria ))
  }
  
  {
    
    otpt_df<- counts_metadata[counts_metadata$Sample_Type != 'TE',]
    otpt_df<- otpt_df[, c('Sample_Name', 'Staining_Protocol', 'Expt_Date', 'Date_Measurement',
                          'Location', 'Expt_No', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate', 'c_Bacteria', 'c_HNA', 'c_LNA', 
                          'c_Viruses', 'c_V1', 'c_V2', 'c_V3', 'VBR', 'HNAperLNA')]
    otpt_df<- otpt_df[
      with(otpt_df,
           order(
             otpt_df[, 'Staining_Protocol'],
             otpt_df[, 'Location'],
             otpt_df[, 'Expt_No'],
             otpt_df[, 'Depth'],
             otpt_df[, 'Timepoint'],
             otpt_df[, 'Replicate'],
             otpt_df[, 'Sample_Type']
             
           )),
    ]
  }
  if (write_csv == T){
    write.csv(otpt_df, paste0("./results/", project_title, ".csv"), row.names = F)
  }
  print("Counts were adjusted with TE")
  return(otpt_df)
  
}

bacterial_count_overview_plots<- function(data){
  plot1<- ggplot(data = NJ2020, aes(x = Sample_Type, y = c_Bacteria, col = as.factor(Expt_No), group = interaction(Sample_Type,Expt_No)))+
    geom_violin(aes(group = Sample_Type, fill = Sample_Type), alpha = 0.15, color = NA)+
    labs(fill = "Sample Type")+
    scale_fill_lancet()+
    geom_boxplot( fill = 'white', outlier.alpha =0)+
    guides(color = guide_legend(title = "Expt No"))+
    geom_jitter(alpha = 1.0, aes(shape = as.factor(Timepoint)), position = position_jitterdodge(jitter.width = 2, dodge.width = 0.9))+
    labs(shape = "Timepoint")+
    scale_color_lancet()+
    theme_classic()+
    labs(title = "Overview of Bacterial Counts per Sample Type")+
    xlab("Sample Type")+
    ylab("Total Virus Count")+
    scale_shape_manual(values = c(19, 0, 1, 2, 3, 4, 7, 8, 10, 15))
  
  print(plot1)
  ggsave("./results/TotalBacteria_perSampleType.svg", width = 20, height = 20, units = "cm")
  
  
  plot2<- ggplot(data = NJ2020, aes(x = as.factor(Expt_No), y = c_Bacteria, col = as.factor(Sample_Type), group = interaction(Sample_Type,Expt_No)))+
    geom_violin(aes(group = as.factor(Expt_No), fill = as.factor(Expt_No)), alpha = 0.25, color = NA)+
    labs(fill = "Expt No")+
    scale_fill_lancet()+
    geom_boxplot( fill = 'white', outlier.alpha =0)+
    guides(color = guide_legend(title = "Expt No"))+
    geom_jitter(alpha = 1.0, aes(shape = as.factor(Timepoint)), position = position_jitterdodge(jitter.width = 2, dodge.width = 0.9))+
    labs(shape = "Timepoint")+
    scale_color_lancet()+
    theme_classic()+
    labs(title = "Overview of Bacterial Counts per Expt No")+
    xlab("Sample Type")+
    ylab("Total Virus Count")+
    scale_shape_manual(values = c(19, 0, 1, 2, 3, 4, 7, 8, 10, 15))
  
  print(plot2)
  ggsave("./results/TotalBacteria_perExptNo.svg", width = 30, height = 20, units = "cm")
  
}

viral_count_overview_plots<- function(data){
  plot1<- ggplot(data = NJ2020, aes(x = Sample_Type, y = c_Viruses, col = as.factor(Expt_No), group = interaction(Sample_Type,Expt_No)))+
    geom_violin(aes(group = Sample_Type, fill = Sample_Type), alpha = 0.15, color = NA)+
    labs(fill = "Sample Type")+
    scale_fill_lancet()+
    geom_boxplot( fill = 'white', outlier.alpha =0)+
    guides(color = guide_legend(title = "Expt No"))+
    geom_jitter(alpha = 1.0, aes(shape = as.factor(Timepoint)), position = position_jitterdodge(jitter.width = 2, dodge.width = 0.9))+
    labs(shape = "Timepoint")+
    scale_color_lancet()+
    theme_classic()+
    labs(title = "Overview of Viral Counts per Sample Type")+
    xlab("Sample Type")+
    ylab("Total Virus Count")+
    scale_shape_manual(values = c(19, 0, 1, 2, 3, 4, 7, 8, 10, 15))
  
  print(plot1)
  ggsave("./results/TotalViruses_perSampleType.svg", width = 20, height = 20, units = "cm")
  
  plot2<- ggplot(data = NJ2020, aes(x = as.factor(Expt_No), y = c_Viruses, col = as.factor(Sample_Type), group = interaction(Sample_Type,Expt_No)))+
    geom_violin(aes(group = as.factor(Expt_No), fill = as.factor(Expt_No)), alpha = 0.25, color = NA)+
    labs(fill = "Expt No")+
    scale_fill_lancet()+
    geom_boxplot( fill = 'white', outlier.alpha =0)+
    guides(color = guide_legend(title = "Expt No"))+
    geom_jitter(alpha = 1.0, aes(shape = as.factor(Timepoint)), position = position_jitterdodge(jitter.width = 2, dodge.width = 0.9))+
    labs(shape = "Timepoint")+
    scale_color_lancet()+
    theme_classic()+
    labs(title = "Overview of Viral Counts per Expt No")+
    xlab("Sample Type")+
    ylab("Total Virus Count")+
    scale_shape_manual(values = c(19, 0, 1, 2, 3, 4, 7, 8, 10, 15))
  
  
  print(plot2)
  ggsave("./results/TotalViruses_perExptNo.svg", width = 30, height = 20, units = "cm")
  
}
