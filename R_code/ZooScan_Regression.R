###NMDS for Krill
## Load Packages
library(tidyverse) 
library(lubridate)
library(ggrepel)
library(ggpubr)
library(nlme)
library(MuMIn)
library(sjPlot)
library(multcomp)
library(agricolae)
library(car)
# library(Rmisc)
library(scales)
library(mgcv)

#import updated data
line80 <- read_csv('/Users/gammonkoval/Documents/GitHub/gammon_sol/IRA_Krill/data/ZooScan_Colocated_Line80_refine.csv')
line80 <- line80 %>% drop_na(acoustic_backscatter)

line90 <- read_csv('/Users/gammonkoval/Documents/GitHub/gammon_sol/IRA_Krill/data/ZooScan_Colocated_Line90_refine.csv')
line90 <- line90 %>% drop_na(acoustic_backscatter)

data <- rbind(line80, line90)
data <- data %>%
  filter(depth <= 210)
#rename biomass columns
colnames(data)[colnames(data) == 'appendicularia Estimated C Biomass (mgC m-2)'] <- 'appendicularia'
colnames(data)[colnames(data) == 'bryozoan_larvae Estimated C Biomass (mgC m-2)'] <- 'bryozoan_larvae'
colnames(data)[colnames(data) == 'chaetognatha Estimated C Biomass (mgC m-2)'] <- 'chaetognatha'
colnames(data)[colnames(data) == 'cnidaria_ctenophores Estimated C Biomass (mgC m-2)'] <- 'cnidaria_ctenophores'
colnames(data)[colnames(data) == 'copepoda_calanoida Estimated C Biomass (mgC m-2)'] <- 'copepoda_calanoida'
colnames(data)[colnames(data) == 'copepoda_eucalanids Estimated C Biomass (mgC m-2)'] <- 'copepoda_eucalanids'
colnames(data)[colnames(data) == 'copepoda_harpacticoida Estimated C Biomass (mgC m-2)'] <- 'copepoda_harpacticoida'
colnames(data)[colnames(data) == 'copepoda_oithona_like Estimated C Biomass (mgC m-2)'] <- 'copepoda_oithona_like'
colnames(data)[colnames(data) == 'copepoda_others Estimated C Biomass (mgC m-2)'] <- 'copepoda_others'
colnames(data)[colnames(data) == 'copepoda_poecilostomatoids Estimated C Biomass (mgC m-2)'] <- 'copepoda_poecilostomatoids'
colnames(data)[colnames(data) == 'crustacea_others Estimated C Biomass (mgC m-2)'] <- 'crustacea_others'
colnames(data)[colnames(data) == 'doliolids Estimated C Biomass (mgC m-2)'] <- 'doliolids'
colnames(data)[colnames(data) == 'euphausiids Estimated C Biomass (mgC m-2)'] <- 'euphausiids'
colnames(data)[colnames(data) == 'eggs Estimated C Biomass (mgC m-2)'] <- 'eggs'
colnames(data)[colnames(data) == 'multiples Estimated C Biomass (mgC m-2)'] <- 'multiples'
colnames(data)[colnames(data) == 'nauplii Estimated C Biomass (mgC m-2)'] <- 'nauplii'
colnames(data)[colnames(data) == 'ostracods Estimated C Biomass (mgC m-2)'] <- 'ostracods'
colnames(data)[colnames(data) == 'others Estimated C Biomass (mgC m-2)'] <- 'others'
colnames(data)[colnames(data) == 'polychaete Estimated C Biomass (mgC m-2)'] <- 'polychaete'
colnames(data)[colnames(data) == 'pteropoda Estimated C Biomass (mgC m-2)'] <- 'pteropoda'
colnames(data)[colnames(data) == 'pyrosomes Estimated C Biomass (mgC m-2)'] <- 'pyrosomes'
colnames(data)[colnames(data) == 'rhizaria Estimated C Biomass (mgC m-2)'] <- 'rhizaria'
colnames(data)[colnames(data) == 'salps Estimated C Biomass (mgC m-2)'] <- 'salps'

#remove abundance columns
data <- data[,!grepl("Abundance", colnames(data))]

#integrate on depth
data_summary <- data %>%
  pivot_longer(cols = 'appendicularia':'salps',
               names_to = 'spp',
               values_to = 'biomass') %>%
  dplyr::group_by(Cruise, DateTime, line, station, spp) %>%
  summarize(biomass = mean(log10(biomass + 1)),
            acoustic_backscatter = mean(acoustic_backscatter)) %>%
  spread(key = spp, value = biomass) %>%
  ungroup()

###Full Species GAM
gam.control(maxit = 999)
m1 <- gam(acoustic_backscatter ~ s(appendicularia) + s(bryozoan_larvae) + s(chaetognatha) +
            s(cnidaria_ctenophores) + s(copepoda_calanoida) + s(copepoda_eucalanids) +
            s(copepoda_harpacticoida) + s(copepoda_oithona_like) + s(copepoda_others) +
            s(copepoda_poecilostomatoids) + s(crustacea_others) + s(doliolids) +
            s(euphausiids) + s(eggs) + s(multiples) +
            s(nauplii) + s(ostracods) + s(others) +
            s(polychaete) + s(pteropoda) + s(pyrosomes) +
            s(rhizaria) + s(salps), data = data_summary, select = TRUE)

m1 <- gam(acoustic_backscatter ~ s(copepoda_calanoida) + s(copepoda_eucalanids) +
            s(copepoda_harpacticoida) + s(copepoda_oithona_like) + s(copepoda_others) +
            s(copepoda_poecilostomatoids), data = data_summary, select = TRUE)

summary(m1)
plot(m1)


###Grouped GAM
data <- rbind(line80, line90)
data <- data %>%
  filter(depth <= 210)

spp_group_dict <- data.frame('spp' = c('appendicularia Estimated C Biomass (mgC m-2)', "bryozoan_larvae Estimated C Biomass (mgC m-2)",
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
data_data <- data %>%
  pivot_longer(cols = 'appendicularia Estimated C Biomass (mgC m-2)':"salps Estimated C Biomass (mgC m-2)",
               names_to = 'spp',
               values_to = 'biomass')
data_data <- left_join(data_data, spp_group_dict)

#summarize the data by groups and convert back to wide
data_data <- data_data %>%
  drop_na(Spp_Group) %>%
  group_by(Cruise, DateTime, line, station, Spp_Group) %>%
  dplyr::summarise(biomass = mean(biomass),
                   acoustic_backscatter = mean(acoustic_backscatter)) %>%
  ungroup() %>%
  spread(key = Spp_Group, value = biomass)

data_data[is.na(data_data)] <- 0

m1 <- gam(acoustic_backscatter ~ s(euphausiids) + s(others) + s(copepods) +
            s(gelatinous), data = data_data, select = FALSE)
summary(m1)
par(mfrow = c(2, 2))
plot(m1)



