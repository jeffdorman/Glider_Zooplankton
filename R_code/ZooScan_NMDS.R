###NMDS for Krill
## Load Packages
library(tidyverse) 
library(lubridate)
library(vegan)
library(ggrepel)
library(ggpubr)
library(nlme)
library(MuMIn)
library(sjPlot)
library(multcomp)
library(agricolae)
library(car)
library(Rmisc)
library(scales)
library(pairwiseAdonis)


###Load Data####################################################################
#import updated data
data <- read_csv('/Users/gammonkoval/Library/CloudStorage/GoogleDrive-gkoval@faralloninstitute.org/Shared drives/Commons/PROJECTS/IOOS_GliderZooplankton/data/calcofi/ZooScan/ZooScan_combined.csv')

#split cruise column
data <- separate(data, Cruise, into = c("Year", "Season", "Ship_Code"), sep = c(4, 6))
data <- transform(data, Season = as.numeric(Season))

#add column for season
seasons = function(x){
  if(x %in% 1:3) return("Winter")
  if(x %in% 4:6) return("Spring")
  if(x %in% 7:9) return("Summer")
  if(x %in% c(10:12)) return("Fall")
}

data$Season = sapply(data$Season, seasons)

#add column for day or night
data$Hour <- hour(data$DateTime)
data$Day_Night <- ifelse(data$Hour<7 | 19<data$Hour, "Night","Day")

#convert data from wide to long
data_data <- data %>% 
  dplyr::select(-contains("Abundance")) %>% 
  pivot_longer(
    cols = `appendicularia Estimated C Biomass (mgC m-2)`:`salps Estimated C Biomass (mgC m-2)`,
    names_to = 'Species',
    values_to = 'biomass'
  )

#create dictionary of groupings
spp_group_dict <- data.frame('Species' = c('appendicularia Estimated C Biomass (mgC m-2)', "bryozoan_larvae Estimated C Biomass (mgC m-2)", 
                                           "chaetognatha Estimated C Biomass (mgC m-2)", "cnidaria_ctenophores Estimated C Biomass (mgC m-2)", 
                                           "copepoda_calanoida Estimated C Biomass (mgC m-2)", "copepoda_eucalanids Estimated C Biomass (mgC m-2)", 
                                           "copepoda_harpacticoida Estimated C Biomass (mgC m-2)", "copepoda_oithona_like Estimated C Biomass (mgC m-2)", 
                                           "copepoda_others Estimated C Biomass (mgC m-2)", "copepoda_poecilostomatoids Estimated C Biomass (mgC m-2)", 
                                           "crustacea_others Estimated C Biomass (mgC m-2)", "doliolids Estimated C Biomass (mgC m-2)", 
                                           "eggs Estimated C Biomass (mgC m-2)", "euphausiids Estimated C Biomass (mgC m-2)", 
                                           "multiples Estimated C Biomass (mgC m-2)", "nauplii Estimated C Biomass (mgC m-2)", 
                                           "ostracods Estimated C Biomass (mgC m-2)", "others Estimated C Biomass (mgC m-2)",
                                           "polychaete Estimated C Biomass (mgC m-2)", "pteropoda Estimated C Biomass (mgC m-2)", 
                                           "pyrosomes Estimated C Biomass (mgC m-2)", "rhizaria Estimated C Biomass (mgC m-2)", 
                                           "salps Estimated C Biomass (mgC m-2)"),
                             'Spp_Group' = c('others', 'others', 'others', 
                                         'gelatinous', 'copepods', 'copepods', 
                                         'copepods', 'copepods', 'copepods', 
                                         'copepods', 'others', 'gelatinous', 
                                         'gelatinous', 'euphausiids', 'others', 
                                         'others', 'others', 'others', 'others', 
                                         'others', 'gelatinous', 'others', 'gelatinous'))

#join back to data to get the group
data_data <- left_join(data_data, spp_group_dict)

#summarize the data by groups and convert back to wide
data_data <- data_data %>% 
  group_by(Year, Season, Ship_Code, DateTime, station, season,
           line, latitude, longitude, Day_Night, Spp_Group) %>% 
  dplyr::summarise(biomass = mean(biomass)) %>% 
  ungroup() %>% 
  spread(key = Spp_Group, value = biomass) 

#extract groups
data_group <- data_data %>% dplyr::select(c(Year, Season, Ship_Code, DateTime, station, 
                                      season, line, latitude, longitude, Day_Night))
data_data <- data_data %>% dplyr::select(c(copepods, euphausiids, gelatinous, others))

#run nMDS
krill_nmds <- metaMDS(as.matrix(data_data), distance = "bray", 
                      autotransform = TRUE, model = 'global', weakties = TRUE,
                      k = 2, maxit = 999, try = 20, trymax = 10000)

krill_nmds #stress = 0.08600059
stressplot(krill_nmds)

#extract stores from nMDS
data.scores <- as.data.frame(scores(krill_nmds, display = c('sites')))
data.scores$site <- rownames(data.scores)
data.scores$Year <- data_group$Year
data.scores$Season <- data_group$Season
data.scores$Ship_Code <- data_group$Ship_Code
data.scores$DateTime <- data_group$DateTime
data.scores$station <- data_group$station
data.scores$line <- data_group$line
data.scores$latitude <- data_group$latitude
data.scores$longitude <- data_group$longitude
data.scores$Day_Night <- data_group$Day_Night

#convert data to factors
data.scores$Day_Night <- factor(data.scores$Day_Night, levels = c("Day", "Night"))
data.scores$Season <- factor(data.scores$Season, levels = c("Spring", "Summer",
                                                            "Fall", "Winter"))
data.scores$line <- factor(data.scores$line, levels = c(80, 90))

#add species vectors 
vec.sp <- envfit(krill_nmds, as.matrix(sqrt(data_data)), perm = 1000)
vec.sp.df <- as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$species <- rownames(vec.sp.df)
# vec.sp.df <- vec.sp.df |> separate_wider_delim(species, 
#                                   delim = "Estimated", 
#                                   names = c("species_name", "units"))

#nMDS by season with species vectors
ggplot(data = data.scores) + 
  geom_hline(yintercept = 0, color = 'grey', linetype = 'dashed') + 
  geom_vline(xintercept = 0, color = 'grey', linetype = 'dashed') + 
  geom_point(aes(x = NMDS1, y = NMDS2), size = 1.5, alpha = 0.2) + 
  annotate(geom = "text", x = 0.75, y = 1.5, label = "Stress = 0.086", size = 6) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Season), level = 0.95, size = 1,
               linetype = 'dashed') +
  geom_segment(data = vec.sp.df, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.5, "cm")), color = "grey50",
               inherit.aes = FALSE, linewidth = 0.65) +
  # scale_color_manual(name = "Season", values = c("dodgerblue","orchid","forestgreen","darkorange")) +
  geom_text_repel(data = vec.sp.df, aes(x = NMDS1, y = NMDS2, label = species),
                  size = 6, segment.color = "white") +
  theme_bw() + 
  # ylim(-0.8, 0.8) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  labs(title = 'nMDS of ZooScan') + 
  # labs(title = "nMDS of SHAPES Krill by Survey Direction") + 
  theme(plot.title = element_text(size = 24),
        plot.subtitle = element_text(size = 22),
        plot.caption = element_text(size = 12),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 18))


#run PERMANOVA by Day/Night
adon.results <- adonis2(wisconsin(sqrt(as.matrix(data_data))) ~ data_group$Day_Night, 
                        method = "euclidian", perm = 999)
print(adon.results) #p = 0.001 

#run PERMANOVA by Line
adon.results <- adonis2(wisconsin(sqrt(as.matrix(data_data)))  ~ as.factor(data_group$line),
                        method = "euclidian", perm = 999)
print(adon.results) #p = 0.053 

#run PERMANOVA by Season
adon.results <- adonis2(wisconsin(sqrt(as.matrix(data_data))) ~ data_group$Season, 
                        method = "euclidian", perm = 999)
print(adon.results) #p = 0.001 

#pairwise test for Season
pairwise.adonis2(wisconsin(sqrt(as.matrix(data_data)))  ~ Season, dat = data_group, 
                 sim.method='bray', p.adjust.m='tukey')

#Summer_vs_Fall p = 0.001 ***
#Summer_vs_Winter p = 0.001 ***
#Summer_vs_Spring p = 0.001 ***
#Fall_vs_Winter p = 0.001 ***
#Fall_vs_Spring p = 0.001 ***
#Winter_vs_Spring p = 0.001 ***

#run PERMANOVA by Season
adon.results <- adonis2(wisconsin(sqrt(as.matrix(data_data))) ~ data_group$Season *
                        as.factor(data_group$line), 
                        by = 'terms',
                        method = "euclidian", perm = 999)
print(adon.results) #p = 0.001 
#pairwise test for Season
pairwise.adonis2(wisconsin(sqrt(as.matrix(data_data)))  ~ Season * Line, 
                 dat = data_group, 
                 sim.method='bray', 
                 p.adjust.m='tukey')




###Stacked Bar Plot#############################################################
data_bar <- cbind(data_group, data_data)
data_bar <- data_bar %>% 
  pivot_longer(
    # cols = `appendicularia Estimated C Biomass (mgC m-2)`:`salps Estimated C Biomass (mgC m-2)`,
    cols = copepods:others,
    names_to = 'Species',
    values_to = 'biomass'
  )

data_bar <- separate_wider_delim(data_bar, cols = Species, delim = "Estimated", names = c("Species", "units"))
data_bar$Year<- substr(data_bar$Year, 3, 4)

data_bar %>% 
  filter(season != 'nan') %>% 
ggplot(., aes(fill=Species, y=biomass, x=season)) + 
  geom_bar(position="fill", stat="identity") + 
  theme_bw() + 
  # ylim(-0.8, 0.8) +
  ylab("Percent of Biomass (%)") +
  # facet_wrap(.~season) + 
  # labs(title = 'Annual Distribution of Biomass from ZooScan') + 
  # labs(title = "nMDS of SHAPES Krill by Survey Direction") + 
  theme(plot.title = element_text(size = 24),
        plot.subtitle = element_text(size = 22),
        plot.caption = element_text(size = 12),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 18))




