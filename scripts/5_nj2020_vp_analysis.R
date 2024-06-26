#### 0.0 Setting up ####

library(tidyverse)


nj2020_vp_all<- read.csv("./results/nj2020_viral_production/vp_results_ALL.csv")

unique(nj2020_vp_all$VP_Method)
unique(nj2020_vp_all$Station_Number)
str(nj2020_vp_all)

nj2020_vp_all <- nj2020_vp_all %>%
  mutate(Station_Number = as.numeric(Station_Number),
         Time_Range = as.character(Time_Range))
str(nj2020_vp_all)

nj2020_vp_filtered <- nj2020_vp_all %>%
  dplyr::filter(Population == "c_Viruses",
                Sample_Type != "VPC")

#### 1.0 Comparing lytic and lysogenic viral production between VIPCAL and VIPCAL-SE at T0-T24 ####

# 1,1 Wrangling data frame
nj2020_vp_24<- nj2020_vp_filtered %>%
  dplyr::filter(Time_Range == 'T0_T24',
                VP_Method %in% c("VPCL_AR_DIFF", "VPCL_AR_DIFF_SE"))

unique(nj2020_vp_24$Time_Range)
unique(nj2020_vp_24$VP_Method)

# 1.2 Kruskal tests
 
kruskal.test(VP ~ VP_Method, data = nj2020_vp_24 %>%
               dplyr::filter(Sample_Type == 'VP')) # not significant
kruskal.test(VP ~ VP_Method, data = nj2020_vp_24 %>%
               dplyr::filter(Sample_Type == 'Diff')) #barely not significant

# 1.3 Linear regression
replacements <- c("VPCL_AR_DIFF" = "VIPCAL", "VPCL_AR_DIFF_SE" = "VIPCAL_SE",
                  "VP" = "Lytic", "Diff" = "Lysogenic")
nj2020_vp_plot_df<- nj2020_vp_24 %>%
  select(!c("abs_VP", "VP_R_Squared")) %>%
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
    Sample_Type = factor(Sample_Type, levels = c("Lytic", "Lysogenic")),
    Station_Number = factor(Station_Number)
  )


ggplot(nj2020_vp_plot_df %>%
         dplyr::filter(Sample_Type == 'Lytic'), aes(x = Station_Number, y = VP, fill = VP_Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
  geom_errorbar(aes(ymin = VP - VP_SE, ymax = VP + VP_SE), 
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(title = "NJ2020 Lytic Viral Production",
       x = "Station Number",
       y = "Viral Production (VP)",
       fill = "VP Methods") + 
  scale_fill_manual(values = c( "VIPCAL" = "#d43028", "VIPCAL_SE" = "#4d778b")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.0))) + 
  theme_classic()

ggplot(nj2020_vp_plot_df %>%
         dplyr::filter(Sample_Type == 'Lysogenic'), aes(x = Station_Number, y = VP, fill = VP_Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
  geom_errorbar(aes(ymin = VP - VP_SE, ymax = VP + VP_SE), 
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(title = "NJ2020 Lysogenic Viral Production",
       x = "Station Number",
       y = "Viral Production (VP)",
       fill = "VP Methods") +
  scale_fill_manual(values = c( "VIPCAL" = "#d43028", "VIPCAL_SE" = "#4d778b")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
  theme_classic()

#Facet wrap
ggplot(nj2020_vp_plot_df ,  # Filter for both conditions
       aes(x = Station_Number, y = VP, fill = VP_Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
  geom_errorbar(aes(ymin = VP - VP_SE, ymax = VP + VP_SE),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(title = "NJ2020 Viral Production",
       x = "Station Number",
       y = "Viral Production (VP)",
       fill = "VP Methods") +
  scale_fill_manual(values = c("VIPCAL" = "#d43028", "VIPCAL_SE" = "#4d778b")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
  facet_wrap(~ Sample_Type, ncol = 2, scales = "fixed") +  # Use 'facet_wrap' with free y-axis scaling
  theme_classic() +
  theme(strip.background = element_blank() )


#### 2.0 Comparing bacterial endpoint vs T24 ####

nj2020_vp_24<- read.csv("./results/nj2020_viral_production/vp_results_24.csv") %>%
  mutate(Timepoint = 'T24') 
nj2020_vp_bp<- read.csv("./results/nj2020_viral_production/vp_results_BP.csv") %>%
  mutate(Timepoint = 'BP')

nj2020_24_bp <- bind_rows(nj2020_vp_24, nj2020_vp_bp) %>%
  dplyr::filter(VP_Method == "VPCL_AR_DIFF_SE",
                Population == 'c_Viruses',
                Sample_Type != 'VPC')


replacements <- c("VPCL_AR_DIFF_SE" = "VIPCAL_SE",
                  "VP" = "Lytic", "Diff" = "Lysogenic")
nj2020_24_bp <- bind_rows(nj2020_vp_24, nj2020_vp_bp) %>%
  dplyr::filter(VP_Method == "VPCL_AR_DIFF_SE",
                Population == 'c_Viruses',
                Sample_Type != 'VPC') %>%
  select(!c("abs_VP", "VP_R_Squared")) %>%
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
    Sample_Type = factor(Sample_Type, levels = c("Lytic", "Lysogenic")),
    Station_Number = factor(Station_Number)
  )

ggplot(nj2020_24_bp ,  
       aes(x = Station_Number, y = VP, fill = Timepoint)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7, color = 'black') +
  geom_errorbar(aes(ymin = VP - VP_SE, ymax = VP + VP_SE),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(title = "NJ2020 Viral Production",
       x = "Station Number",
       y = "Viral Production (VP)",
       fill = "Assay Duration") +
  scale_fill_manual(values = c("T24" = "#d43028", "BP" = "#4d778b")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
  facet_wrap(~ Sample_Type, ncol = 2, scales = "fixed") +  # Use 'facet_wrap' with free y-axis scaling
  theme_classic() +
  theme(strip.background = element_blank() )



#### 3.0 Contact rates ####

cr_df<- read.csv("./results/nj2020_fcm/nj2020_fcm_corrected_counts.csv")
library(tidyverse)
library(data.table)
library(cowplot)
library(ggpubr)
library(colorspace)







cr_df<- cr_df %>% select(c("Location", "Station_Number", "Sample_Type", "Timepoint",
                           "Replicate", "c_Bacteria", "c_Viruses", "VBR"))

abundance<- read.csv("./data/metadata/nj2020_original_abundances.csv")

abundance1<- abundance %>%
  select(Location, Station_Number, Total_Bacteria, Total_Viruses, VBR) %>%
  rename(c_Bacteria = Total_Bacteria,
         c_Viruses = Total_Viruses )

loc<- 'NJ2020'
expt_no<- 1:7
type<- c('VP', 'VPC', '0.22')
time<- -3
rep<- 1:3
abundance2<- data.frame(expand.grid(loc, expt_no, type, time, rep))
colnames(abundance2)<- c('Location', 'Station_Number', 'Sample_Type',
                         'Timepoint', 'Replicate')




abundance3<- full_join(abundance2, abundance1)


cr_df<- full_join(abundance3, cr_df)


cr_df$BV<-  cr_df$c_Bacteria * cr_df$c_Viruses

cr_list<- list()

for (loc in unique(cr_df$Location)){
  for (expt in unique(cr_df$Station_Number)){
    for(type in unique(cr_df$Sample_Type)){
      for (rep in unique(cr_df$Replicate)){
        
        df10<- cr_df %>%
          dplyr::filter(Location == loc,
                 Station_Number == expt,
                 Sample_Type == type,
                 Replicate == rep)
        
        for(time in unique(cr_df$Timepoint)){
          
          cr<-  (df10 %>% dplyr::filter(Timepoint == time))$BV/(df10 %>% dplyr::filter(Timepoint == 0))$BV
          print(cr)
          cr_value<- c(loc, expt, type, rep, time, cr)
          
          cr_list[[length(cr_list)+1]]<- cr_value
        }
        
      }
    }
  }
}


cr_opt<- data.table(data.table::transpose(as.data.table(cr_list)))
colnames(cr_opt)<- c("Location", "Station_Number", "Sample_Type","Replicate", "Timepoint",
                     "rate")
cr_opt$rate <- as.numeric(cr_opt$rate)
#Warnings are because of o.22 samples that have negative amount of bacteria. This could be due to TE corection.
#We either get call the bacterial count as 0 and calculate for 0.22 or exclude it from further analyses.
#For now, we'll exclude it.
cr_opt <- cr_opt %>% dplyr::filter(Sample_Type != '0.22')


#Let's visualise this data #####

#First we summarise the replicates
cr_opt<- cr_opt %>% dplyr::filter(Timepoint != -3)

cr_opt <- cr_opt %>% 
  group_by(Location, Station_Number, Sample_Type, Timepoint) %>%
  summarise(rate_mean = mean(rate), rate_se = plotrix::std.error(rate))

cr_opt2<- cr_opt %>% ungroup() %>%
  select(-rate_se) %>%
  pivot_wider(names_from = Sample_Type,
              values_from = rate_mean) %>%
  mutate(cr_diff = VP - VPC)


cr_opt[cr_opt$Station_Number == 1 & cr_opt$Timepoint == 9,]$Timepoint<- '17'
cr_opt[cr_opt$Station_Number == 2 & cr_opt$Timepoint == 9,]$Timepoint<- '17.5'
cr_opt[cr_opt$Station_Number == 1 & cr_opt$Timepoint == 12,]$Timepoint<- '20'
cr_opt[cr_opt$Station_Number == 2 & cr_opt$Timepoint == 12,]$Timepoint<- '20.5'
cr_opt$Timepoint<- as.numeric(cr_opt$Timepoint)

cr_opt2[cr_opt2$Station_Number == 1 & cr_opt2$Timepoint == 9,]$Timepoint<- '17'
cr_opt2[cr_opt2$Station_Number == 2 & cr_opt2$Timepoint == 9,]$Timepoint<- '17.5'
cr_opt2[cr_opt2$Station_Number == 1 & cr_opt2$Timepoint == 12,]$Timepoint<- '20'
cr_opt2[cr_opt2$Station_Number == 2 & cr_opt2$Timepoint == 12,]$Timepoint<- '20.5'
cr_opt2$Timepoint<- as.numeric(cr_opt2$Timepoint)
cols<- c("VP" = '#d32f27ff',
         "VPC" = '#4888a2ff',
         "Diff" = '#edb81dff')
shapes<- c("VP" = 21,
           "VPC" = 24,
           "Diff" = 19)
shapes<- c("VP" = 15,
           "VPC" = 17,
           "Diff" = 19)


# Now let's overlap bacterial production values here
bep_df <- read_csv("./results/nj2020_viral_production/vp_results_BP.csv") 
bep_df <- bep_df %>%
  select(Location, Station_Number,  Time_Range) %>%
  unique() %>%
  mutate(endpoint = as.numeric(sub(".*\\T0_T?\\s*(\\d+)", "\\1", Time_Range)) )#function to extract the endpoint

bep_df <- merge(cr_opt2, bep_df)



cr_plot1<- ggplot()+
  geom_vline(data = bep_df,
             aes(xintercept = endpoint),
             colour = '#7E6148FF',
             linewidth =1.5,
             alpha = 0.5,
             show.legend = F)+
  geom_point(data = cr_opt,
             aes(x = as.numeric(Timepoint),
                 y = rate_mean,
                 fill = Sample_Type,
                 shape = Sample_Type,
                 colour = Sample_Type),
             show.legend =  T)+
  geom_line(data = cr_opt2,
            aes(x = as.numeric(Timepoint),
                y = cr_diff),
            size = 1,
            alpha = 0.5,
            show.legend = F) +
  theme_bw()+
  facet_wrap(. ~ Station_Number, ncol = 4
             # scales = 'free',
  )+
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols#, breaks = 'VP', labels = 'Bacterial Endpoint'
  )+
  scale_shape_manual(values = shapes)+
  theme(strip.background = element_lank,
        strip.text = element_text(face = 'bold',
                                  color = 'white',
                                  size = 10),
        panel.border = element_rect(linewidth = 2),
        panel.background = element_rect(fill = NA),
        #legend.position = c(0.75, 0.1),
        legend.title = element_text(face = 'bold',
                                    size = 10),
        legend.text = element_text(size = 9),
        axis.title = element_text(face = 'bold',
                                  size = 10),
        axis.text = element_text(size = 10))+
  labs(x = 'Timepoint (in hours)',
       y = 'Mean Relative Collision Rate')  +
  guides(color = guide_legend(title="Treatment"),
         shape = guide_legend(title="Treatment"),
         fill= guide_legend(title="Treatment"))

cr_plot1 

legend1<- ggpubr::as_ggplot( ggpubr::get_legend(cr_plot1))

cr_plot<- cr_plot1 +theme(legend.position = 'none')

cr_plot2<- ggplot()+
  geom_vline(data = bep_df,
             aes(xintercept = endpoint),
             colour = '#7E6148FF',
             linewidth =1.5,
             alpha = 0.5,
             show.legend = T)+
  geom_point(data = cr_opt,
             aes(x = as.numeric(Timepoint),
                 y = rate_mean,
                 #fill = Sample_Type,
                 #shape = Sample_Type,
                 colour = Sample_Type
             ),
             show.legend =  F)+
  # geom_line(data = cr_opt2,
  #           aes(x = as.numeric(Timepoint),
  #               y = cr_diff),
  #           size = 1,
  #           alpha = 0.5,
  #           show.legend = F) +
  theme_bw()+
  facet_wrap(. ~ Station_Number, ncol = 2
             # scales = 'free',
  )+
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols, breaks = 'VP', labels = 'Bacterial\nEndpoint'
  )+
  scale_shape_manual(values = shapes)+
  theme(strip.background = element_rect(fill = 'gray',
                                        color = NA),
        strip.text = element_text(face = 'bold'),
        panel.border = element_rect(linewidth = 2),
        panel.background = element_rect(fill = NA),
        #legend.position = c(0.75, 0.1),
        legend.title = element_text(face = 'bold',
                                    size = 10),
        legend.text = element_text(size = 9),
        axis.title = element_text(face = 'bold',
                                  size = 15),
        axis.text = element_text(size = 15))+
  labs(x = 'Timepoint (in hours)',
       y = 'Mean Relative Collision Rate')  +
  guides(color = guide_legend(title = ''),
         #shape = guide_legend(title="Treatment"),
         #fill= guide_legend(title="Treatment")
  )

cr_plot2

legend2<- ggpubr::as_ggplot( ggpubr::get_legend(cr_plot2))


cr_plot3<- ggplot()+
  # geom_vline(data = bep_df,
  #            aes(xintercept = endpoint),
  #            colour = '#7E6148FF',
  #            linewidth =1.5,
  #            alpha = 0.5,
  #            show.legend = T)+
  geom_point(data = cr_opt,
             aes(x = as.numeric(Timepoint),
                 y = rate_mean,
                 #fill = Sample_Type,
                 #shape = Sample_Type,
                 colour = Sample_Type
             ),
             show.legend =  F)+
  geom_line(data = cr_opt2,
            aes(x = as.numeric(Timepoint),
                y = cr_diff),
            size = 1,
            alpha = 0.5,
            show.legend = T) +
  theme_bw()+
  facet_wrap(. ~ Station_Number, ncol = 2
             # scales = 'free',
  )+
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols, breaks = 'VP', labels = 'Differences in\ncollision rates\nbetween VP\nand VPC\ntreatments'
  )+
  scale_shape_manual(values = shapes)+
  theme(strip.background = element_rect(fill = lighten("#323D5E",
                                                       amount = 0.1) ,
                                        color = NA),
        strip.text = element_text(face = 'bold',
                                  color = 'white'),
        panel.border = element_rect(linewidth = 2),
        panel.background = element_rect(fill = NA),
        #legend.position = c(0.75, 0.1),
        legend.title = element_text(face = 'bold',
                                    size = 10),
        legend.text = element_text(size = 9),
        axis.title = element_text(face = 'bold',
                                  size = 15),
        axis.text = element_text(size = 15))+
  labs(x = 'Timepoint (in hours)',
       y = 'Mean Relative Collision Rate')  +
  guides(color = guide_legend(title = ''),
         #shape = guide_legend(title="Treatment"),
         #fill= guide_legend(title="Treatment")
  )

cr_plot3

legend3<- ggpubr::as_ggplot( ggpubr::get_legend(cr_plot3))

legends<- cowplot::plot_grid(plotlist = list(legend1, legend2, legend3),
                             ncol= 1,
                             rel_widths = c(1, 1,  1))
legends2<- plot_grid( NULL, legends,
                      rel_widths = c(2,3),
                      nrow = 2)



cr_plot_final<- cr_plot + NULL+ legends +
  plot_layout(width = c(8, 2,2))
cr_plot_final



cr_plot_final<- plot_grid(cr_plot, legends,
                          nrow = 1,
                          rel_widths = c( 2, 1))


cr_plot_final


ggsave(file = 'collision_rate_ww.png', width = 10,
       height = 18, unit = 'cm')









