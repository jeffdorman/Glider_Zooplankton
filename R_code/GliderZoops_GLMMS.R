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
library(lme4)
library(merTools)


#open data
zooscan <- read_csv('/Users/gammonkoval/Documents/GitHub/gammon_sol/IRA_Krill/data/ZooScan_Model_Data.csv')
biovol <- read_csv('/Users/gammonkoval/Documents/GitHub/gammon_sol/IRA_Krill/data/BioVolume_Model_Data.csv', na = 'nan')
glider <- read_csv('/Users/gammonkoval/Documents/GitHub/gammon_sol/IRA_Krill/data/Glider_Model_Data.csv', na = 'NaN')

#clean data
zooscan <- zooscan %>% 
  drop_na(season) %>% 
  dplyr::select(c("line", "season", "Inshore_Offshore", "total_biomass", "Day_Night",
                  "total_biomass_log", "year")) %>% 
  group_by(line, season, Inshore_Offshore, Day_Night, year) %>% 
  summarize(total_biomass = mean(total_biomass, na.rm = TRUE),
            total_biomass_log = mean(total_biomass_log, na.rm = TRUE)) %>% 
  ungroup()

biovol <- biovol %>% 
  drop_na(season) %>% 
  rename(line = St_Line) %>% 
  dplyr::filter(line %in% c(80, 90)) %>% 
  dplyr::select(c("line", "season", "Inshore/Offshore", "Norm_Sml_PVolC3", "year")) %>% 
  group_by(line, season, `Inshore/Offshore`, year) %>% 
  summarize(Norm_Sml_PVolC3 = mean(Norm_Sml_PVolC3, na.rm = TRUE)) %>% 
  ungroup()

glider <- glider %>% 
  rename(line = Line) %>% 
  drop_na(season) %>% 
  dplyr::select(c("line", "season", "Inshore_Offshore", "Day_Night", 
                  "big_sv", "little_sv", "year")) %>% 
  mutate(line = as.numeric(str_extract_all(line, "\\d+"))) %>%
  group_by(line, season, Inshore_Offshore, year) %>% 
  summarize(big_sv = mean(big_sv, na.rm = TRUE),
            little_sv = mean(little_sv, na.rm = TRUE)) %>% 
  ungroup() 


#join together
# data <- full_join(zooscan, biovol)
data <- full_join(zooscan, glider)
#total_biomass = sum of biomass from ZooScan 
#total_biomass_log = log10(biomass + 1) transformed biomass from ZooScan
#Norm_Sml_PVolC3 = normalized small volume from biovolume (small volume / water filtered)
#big_sv = log scaled acoustic backscatter from glider
#little_sv = linear scaled acoustic backscatter from glider 


data$line <- as.factor(data$line)
data$season <- factor(data$season, levels = c('Winter', 'Spring', 'Summer', 'Fall'))
data$`Inshore/Offshore` <- as.factor(data$`Inshore/Offshore`)

ggplot(data, aes(x = total_biomass_log)) + 
  geom_histogram(aes(y = ..density..), color = 'black', fill = 'white') + 
  geom_vline(aes(xintercept=mean(total_biomass_log, na.rm = TRUE)), 
             color = 'blue', linetype = 'dashed', linewidth = 1) + 
  theme_bw() + 
  labs(title = 'Log of Total Biomass from ZooScan') +
  stat_function(fun = dnorm, args = list(mean = mean(data$total_biomass_log),
                                         sd = sd(data$total_biomass_log)))

ggplot(data, aes(x = big_sv)) + 
  geom_histogram(aes(y = ..density..), color = 'black', fill = 'white') + 
  geom_vline(aes(xintercept=mean(big_sv, na.rm = TRUE)), 
             color = 'red', linetype = 'dashed', size = 1) + 
  theme_bw() + 
  labs(title = 'Acoustic Backscatter') + 
  stat_function(fun = dnorm, args = list(mean = mean(data$big_sv),
                                         sd = sd(data$big_sv)))



ggplot(data, aes(x = year, y = big_sv)) + 
  geom_point() + 
  geom_smooth(method = 'lm')


data %>% 
  drop_na(big_sv, total_biomass_log) %>% 
  group_by(year) %>% 
  summarize(n())
  group_by(line, season, Inshore_Offshore, Day_Night) %>% 
  summarize(n())

shapiro.test(data$big_sv) # 0.005285
shapiro.test(data$total_biomass_log) #0.001373

f.sbs <- lme(big_sv ~ total_biomass_log + season + `Inshore/Offshore`, 
             random = ~ total_biomass_log | line/year, 
             method = 'REML',
             data = data)




##geom_histogram()###Glider vs ZooScan
# data_90 <- dplyr::filter(data, line == 90)
m2 <- lmer(big_sv ~ total_biomass_log + (1|line) + (1|season) + (1|`Inshore/Offshore`),
           data = data, REML = TRUE)

summary(m2)
r.squaredGLMM(m2)
Anova(m2, type = 'III')

predictInterval(m2)
REsim(m2)
plotREsim((REsim(m2)))

