# 0. Loading essential packages  ####

set.seed(2024)
load_or_install <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

packages_to_load <- c("tidyverse",
                      "ggsci",
                      "data.table",
                      "readxl",
                      "devtools",
                      "dunn.test",
                      "dplyr",
                      "pracma",
                      "cowplot",
                      "ggpubr",
                      "colorspace",
                      "plotrix",
                      "ggnewscale",
                      "patchwork")

lapply(packages_to_load, load_or_install)


# Installing viralprod R package
{
  if (!requireNamespace("viralprod", quietly = TRUE)) {
  devtools::install_github("mdhishamshaikh/ViralProduction_R")
3}
library(viralprod)
}
#lme4 needs to be installed from source
{
  if (!requireNamespace("lme4", quietly = TRUE)) {
    install.packages("lme4", type = "source") 
  }
  library(lme4)
}

#### Plot viral counts ####
plot_viral_abundance <- function(file_path, data, location_column = "Location_Station") {
  # Open a PDF device
  pdf(file = file_path, width = 10, height = 5)
  
  # Loop through each unique location station
  for (loc_station in unique(data[[location_column]])) {
    # Filter data for the current location station
    data_subset <- data %>% filter(!!sym(location_column) == loc_station)
    
    # Create the plot
    vdc_plot <- ggplot(data_subset, aes(x = Timepoint, y = c_Viruses / 1e+6)) +
      geom_line(aes(group = as.factor(Replicate))) +  # Add lines to connect points for each replicate
      geom_point(aes(color = as.factor(Replicate)), size = 2) +  # Color points by Replicate
      facet_grid(. ~ Sample_Type, scales = "fixed") +  # Facet by Sample_Type only
      scale_color_npg() +
      theme_bw(base_size = 15) +
      labs(
        title = paste("Viral Abundance for Location Station:", loc_station),
        x = "Timepoint",
        y = "Viral Count (c_Viruses in millions)",
        color = "Replicate"
      )
    
    # Print the plot to the PDF
    print(vdc_plot)
  }
  
  # Close the PDF device
  dev.off()
}

# Averaging dataframe to get difference curve ####
calculate_mean_viral_count_diff_with_errors <- function(data) {
  # Calculating average and standard error for each Sample_Type
  average_counts <- data %>%
    group_by(Location_Station, Timepoint, Sample_Type) %>%
    summarize(
      avg_c_Viruses = mean(c_Viruses, na.rm = TRUE),
      sem_c_Viruses = sd(c_Viruses, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Calculating Diff (VPC - VP) and its standard error per Timepoint
  diff_counts <- average_counts %>%
    pivot_wider(names_from = Sample_Type, values_from = c(avg_c_Viruses, sem_c_Viruses)) %>%
    mutate(avg_c_Viruses = avg_c_Viruses_VPC - avg_c_Viruses_VP,
      sem_c_Viruses = sqrt(sem_c_Viruses_VPC^2 + sem_c_Viruses_VP^2)
    ) %>%
    select(Location_Station, Timepoint, avg_c_Viruses, sem_c_Viruses) %>%
    mutate(
      Sample_Type = "Diff")
  
  # Combining averaged data with Diff data
  final_data <- average_counts %>%
    bind_rows(diff_counts)
  
  return(final_data)
}


plot_mean_viral_abundance <- function(file_path, data, location_column = "Location_Station") {
 
  pdf(file = file_path, width = 10, height = 5)
  
  # Looping through each unique location station
  for (loc_station in unique(data[[location_column]])) {
    # Filter data for the current location station
    data_subset <- data %>% filter(!!sym(location_column) == loc_station)
    
    # Creating the plot
    vdc_plot <- ggplot(data_subset, aes(x = Timepoint, y = avg_c_Viruses / 1e+6, color = Sample_Type, group = Sample_Type)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = (avg_c_Viruses - sem_c_Viruses) / 1e+6,
                        ymax = (avg_c_Viruses + sem_c_Viruses) / 1e+6), width = 0.2) +
      scale_color_manual(values = c("VP" = "#ef476f", "VPC" = "#ffd166", "Diff" = "#26547c")) +
      theme_bw(base_size = 15) +
      labs(
        title = paste("Average Viral Abundance for Location Station:", loc_station),
        x = "Timepoint",
        y = "Average Viral Count (c_Viruses in millions)",
        color = "Sample Type"
      )
    
    print(vdc_plot)
  }
  

  dev.off()
}


#### gEOM SPLIT VIOLIN
{ #https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
  GeomSplitViolin <- ggplot2::ggproto(
    "GeomSplitViolin",
    ggplot2::GeomViolin,
    draw_group = function(self,
                          data,
                          ...,
                          # add the nudge here
                          nudge = 0,
                          draw_quantiles = NULL) {
      data <- transform(data,
                        xminv = x - violinwidth * (x - xmin),
                        xmaxv = x + violinwidth * (xmax - x))
      grp <- data[1, "group"]
      newdata <- plyr::arrange(transform(data,
                                         x = if (grp %% 2 == 1) xminv else xmaxv),
                               if (grp %% 2 == 1) y else -y)
      newdata <- rbind(newdata[1, ],
                       newdata,
                       newdata[nrow(newdata), ],
                       newdata[1, ])
      newdata[c(1, nrow(newdata)-1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
      
      # now nudge them apart
      newdata$x <- ifelse(newdata$group %% 2 == 1,
                          newdata$x - nudge,
                          newdata$x + nudge)
      
      if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
        
        stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
        
        quantiles <- ggplot2:::create_quantile_segment_frame(data,
                                                             draw_quantiles)
        aesthetics <- data[rep(1, nrow(quantiles)),
                           setdiff(names(data), c("x", "y")),
                           drop = FALSE]
        aesthetics$alpha <- rep(1, nrow(quantiles))
        both <- cbind(quantiles, aesthetics)
        quantile_grob <- ggplot2::GeomPath$draw_panel(both, ...)
        ggplot2:::ggname("geom_split_violin",
                         grid::grobTree(ggplot2::GeomPolygon$draw_panel(newdata, ...),
                                        quantile_grob))
      }
      else {
        ggplot2:::ggname("geom_split_violin",
                         ggplot2::GeomPolygon$draw_panel(newdata, ...))
      }
    }
  )
  geom_split_violin <- function(mapping = NULL,
                                data = NULL,
                                stat = "ydensity",
                                position = "identity",
                                # nudge param here
                                nudge = 0,
                                ...,
                                draw_quantiles = NULL,
                                trim = TRUE,
                                scale = "area",
                                na.rm = FALSE,
                                show.legend = NA,
                                inherit.aes = TRUE) {
    
    ggplot2::layer(data = data,
                   mapping = mapping,
                   stat = stat,
                   geom = GeomSplitViolin,
                   position = position,
                   show.legend = show.legend,
                   inherit.aes = inherit.aes,
                   params = list(trim = trim,
                                 scale = scale,
                                 # don't forget the nudge
                                 nudge = nudge,
                                 draw_quantiles = draw_quantiles,
                                 na.rm = na.rm,
                                 ...))
  }
  
}
