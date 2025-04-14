###Script for testing how many years are required for relationship
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
library(lme4)
library(ggplot2)
library(gridExtra)
library(suncalc)
library(lubridate)
library(lattice)

#open data
zooscan <- read_csv('/Users/gammonkoval/Documents/GitHub/gammon_sol/IRA_Krill/data/ZooScan_Model_Data.csv')
# biovol <- read_csv('/Users/gammonkoval/Documents/GitHub/gammon_sol/IRA_Krill/data/BioVolume_Model_Data.csv', na = 'nan')
glider <- read_csv('/Users/gammonkoval/Documents/GitHub/gammon_sol/IRA_Krill/data/Glider_Model_Data_newInshoreOffshore_standardized.csv', na = 'NaN')

#group taxa as needed
zooscan$Copepod <- zooscan[, 15] + zooscan[, 17] + zooscan[, 19] + zooscan[, 21] + zooscan[, 23] + zooscan[, 25]
zooscan$Krill <- zooscan[, 37]
zooscan$Gelatinous <- zooscan[, 47] + zooscan[, 51]

#clean data
zooscan <- zooscan %>% 
  drop_na(season) %>% 
  dplyr::select(c("line", "season", "Inshore_Offshore", "total_biomass", 
                  "total_biomass_log", "year","Day_Night","Copepod","Krill","Gelatinous")) %>% 
  group_by(line, season, Inshore_Offshore, year,Day_Night) %>% 
  summarize(total_biomass = mean(total_biomass, na.rm = TRUE),
            total_Cop = mean(Copepod, na.rm = TRUE),
            total_Krill = mean(Krill, na.rm = TRUE),
            total_Gel = mean(Gelatinous, na.rm = TRUE),
            total_biomass_log = mean(total_biomass_log, na.rm = TRUE)) %>% 
  ungroup()

# biovol <- biovol %>% 
#   drop_na(season) %>% 
#   rename(line = St_Line) %>% 
#   dplyr::filter(line %in% c(80, 90)) %>% 
#   dplyr::select(c("line", "season", "Inshore.Offshore", "Norm_Sml_PVolC3", "year")) %>% 
#   group_by(line, season, `Inshore.Offshore`, year) %>% 
#   summarize(Norm_Sml_PVolC3 = mean(Norm_Sml_PVolC3, na.rm = TRUE)) %>% 
#   ungroup()

#glider data needs to be formatted and assigned Day or Night given sunrise/sunset times
#Step 1 convert gmt to local time
glider <- glider %>% 
  rename(line = Line) %>% 
  drop_na(season) %>% 
  mutate(
    time_gmt = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "GMT"),
    time_local = with_tz(time_gmt, tzone = "America/Los_Angeles"),
    local_date = if_else(
      date(time_local) != date(time_gmt), 
      date(time_local), 
      date(time_gmt)
    )
  )

#Step 2 assign Day or Night given Sunset/Sunrise at Los Angeles
glider <- glider %>%
  rowwise() %>%
  mutate(
    sun_times = list(getSunlightTimes(
      date = as.Date(local_date), 
      lat = lat, 
      lon = lon, 
      keep = c("sunrise", "sunset"), 
      tz = "America/Los_Angeles"
    )),
    sunrise = sun_times$sunrise,
    sunset = sun_times$sunset,
    Day_Night = if_else(time_local >= sunrise & time_local < sunset, "Day", "Night")
  ) %>%
  ungroup()

#Step3 group and produce summary df.
glider <- glider %>%
  dplyr::select(c("line", "season", "Inshore_Offshore", "big_sv", "little_sv", "year","Day_Night")) %>% 
  mutate(line = as.numeric(str_extract_all(line, "\\d+"))) %>%
  group_by(line, season, Inshore_Offshore, year, Day_Night) %>% 
  summarize(big_sv = mean(big_sv, na.rm = TRUE),
            little_sv = mean(little_sv, na.rm = TRUE)) %>% 
  ungroup() 


###OLD GAMMON NOTES
##### I have questions about this - does sum biomass include multiples?
#total_biomass = sum of biomass from ZooScan 
#total_biomass_log = log10(biomass + 1) transformed biomass from ZooScan
#Norm_Sml_PVolC3 = normalized small volume from biovolume (small volume / water filtered)
#big_sv = log scaled acoustic backscatter from glider
#little_sv = linear scaled acoustic backscatter from glider 


#join data
# data <- full_join(zooscan, biovol)
data <- full_join(zooscan, glider)

#clean joined data
data$line <- as.factor(data$line)
data$season <- factor(data$season, levels = c('Winter', 'Spring', 'Summer', 'Fall'))
data$Inshore_Offshore <- as.factor(data$Inshore_Offshore)
data$`Day_Night` <- as.factor(data$`Day_Night`)
data$`year` <- as.factor(data$`year`)
data <- data %>% filter(season!="NA") %>% filter(big_sv!="NaN") %>% filter(total_biomass_log!="NA") %>% droplevels()



f.sbs.fin <- lme(total_biomass_log ~ big_sv + season + Inshore_Offshore,
                 random = list(line = ~big_sv, year = ~big_sv),
                 control = list(maxIter = 10000, niterEM = 10000),
                 method = "REML",
                 data = data)

r.squaredGLMM(f.sbs.fin)

test <- data %>% filter(year == 2021)
test$year = 2025
test$predict <- predict(f.sbs.fin, newdata = test)

ggplot(test, aes(x = total_biomass_log, y = predict)) + 
  geom_point()

###Bootstrap random years (without replacement)#################################
#create dataframe to store rsquared results
rsquared_m_mean = NULL
rsquared_c_mean = NULL
rsquared_m_std = NULL
rsquared_c_std = NULL
na_count = NULL
rand_rsquared_m = NULL
rand_rsquared_c = NULL

# n = 17
# r = 2
# 
# factorial(n) / (factorial(r) *  factorial(n - r))

for (sub in c(2:15)) {
  
  for (i in c(1:100)) {
    
    tryCatch(
      expr = {
        years_sel <- sample(x = c(2008:2024), size = sub, replace = FALSE)
        
        data_temp <- dplyr::filter(data, year %in% years_sel)
        f.sbs.fin <- lme(total_biomass_log ~ big_sv + season + Inshore_Offshore,
                         random = list(line = ~big_sv, year = ~big_sv),
                         control = list(maxIter = 10000, niterEM = 10000),
                         method = "REML",
                         data = data_temp)
        
        rand_rsquared_m[i] <- r.squaredGLMM(f.sbs.fin)[1]
        rand_rsquared_c[i] <- r.squaredGLMM(f.sbs.fin)[2]
      },
      error = function(e){
        rand_rsquared_m[i] <- NA
        rand_rsquared_c[i] <- NA
      }
    )
  }
  
  rsquared_m_mean[sub - 1] = mean(rand_rsquared_m, na.rm = TRUE)
  rsquared_c_mean[sub - 1] = mean(rand_rsquared_c, na.rm = TRUE)
  rsquared_m_std[sub - 1] = sd(rand_rsquared_m, na.rm = TRUE)
  rsquared_c_std[sub - 1] = sd(rand_rsquared_c, na.rm = TRUE)
  na_count[sub - 1] = sum(is.na(rand_rsquared_c))
}

rsquared <- do.call(rbind, Map(data.frame, 
                               sample=c(2:15),
                               rsquared_m_mean=rsquared_m_mean,
                               rsquared_m_std=rsquared_m_std,
                               rsquared_c_mean=rsquared_c_mean,
                               rsquared_c_std=rsquared_c_std,
                               na_count=na_count))

colors <- c("Conditional" = "blue", "Marginal" = "red")

ggplot(rsquared, aes(x = sample)) + 
  geom_line(aes(y = rsquared_m_mean, color = "Marginal")) +
  geom_point(aes(y = rsquared_m_mean, color = "Marginal")) +
  geom_ribbon(aes(y = rsquared_m_mean, 
                  ymin = rsquared_m_mean - rsquared_m_std, 
                  ymax = rsquared_m_mean + rsquared_m_std), 
              alpha = .2, fill = 'red') +
  geom_line(aes(y = rsquared_c_mean, color = "Conditional")) +
  geom_point(aes(y = rsquared_c_mean, color = "Conditional")) +
  geom_ribbon(aes(y = rsquared_c_mean, 
                  ymin = rsquared_c_mean - rsquared_c_std, 
                  ymax = rsquared_c_mean + rsquared_c_std), 
              alpha = .2, fill = 'blue') +
  labs(title = 'Random Sampling of Years for ZooScan vs Glider Model', 
       x = 'Number of Years Selected',
       y = expression (~R^2~"from 100 Model Runs"),
       color = 'Legend') + 
  scale_color_manual(values = colors) + 
  theme_bw() 
  
###Bootstrap sequential years###################################################
#create dataframe to store rsquared results
rsquared_m_mean = NULL
rsquared_c_mean = NULL
rsquared_m_std = NULL
rsquared_c_std = NULL
na_count = NULL
rand_rsquared_m = NULL
rand_rsquared_c = NULL

# n = 17
# r = 2
# 
# factorial(n) / (factorial(r) *  factorial(n - r))

for (sub in c(2:15)) {
  
  for (year in c(2008:2024)) {
    years_sel <- c(year:(year+sub))
    if (max(years_sel) <= 2024) {
      print(years_sel)
    }
  }
}

for (sub in c(2:15)) {
  
  for (year in c(2008:2024)) {
    tryCatch(
      expr = {
        years_sel <- c(year:(year+sub))
        
        if (max(years_sel) <= 2024) {
        
          data_temp <- dplyr::filter(data, year %in% years_sel)
          f.sbs.fin <- lme(total_biomass_log ~ big_sv + season + Inshore_Offshore,
                           random = list(line = ~big_sv, year = ~big_sv),
                           control = list(maxIter = 10000, niterEM = 10000),
                           method = "REML",
                           data = data_temp)
          
          rand_rsquared_m[i] <- r.squaredGLMM(f.sbs.fin)[1]
          rand_rsquared_c[i] <- r.squaredGLMM(f.sbs.fin)[2]
        }
      },
      error = function(e){
        rand_rsquared_m[i] <- NA
        rand_rsquared_c[i] <- NA
      }
    )
  }
  
  rsquared_m_mean[sub - 1] = mean(rand_rsquared_m, na.rm = TRUE)
  rsquared_c_mean[sub - 1] = mean(rand_rsquared_c, na.rm = TRUE)
  rsquared_m_std[sub - 1] = sd(rand_rsquared_m, na.rm = TRUE)
  rsquared_c_std[sub - 1] = sd(rand_rsquared_c, na.rm = TRUE)
  na_count[sub - 1] = sum(is.na(rand_rsquared_c))
}

rsquared <- do.call(rbind, Map(data.frame, 
                               sample=c(2:15),
                               rsquared_m_mean=rsquared_m_mean,
                               rsquared_m_std=rsquared_m_std,
                               rsquared_c_mean=rsquared_c_mean,
                               rsquared_c_std=rsquared_c_std,
                               na_count=na_count))

colors <- c("Conditional" = "blue", "Marginal" = "red")

ggplot(rsquared, aes(x = sample)) + 
  geom_line(aes(y = rsquared_m_mean, color = "Marginal")) +
  geom_point(aes(y = rsquared_m_mean, color = "Marginal")) +
  geom_ribbon(aes(y = rsquared_m_mean, 
                  ymin = rsquared_m_mean - rsquared_m_std, 
                  ymax = rsquared_m_mean + rsquared_m_std), 
              alpha = .2, fill = 'red') +
  geom_line(aes(y = rsquared_c_mean, color = "Conditional")) +
  geom_point(aes(y = rsquared_c_mean, color = "Conditional")) +
  geom_ribbon(aes(y = rsquared_c_mean, 
                  ymin = rsquared_c_mean - rsquared_c_std, 
                  ymax = rsquared_c_mean + rsquared_c_std), 
              alpha = .2, fill = 'blue') +
  labs(title = 'Sequential Sampling of Years for ZooScan vs Glider Model', 
       x = 'Number of Years Selected',
       y = expression (~R^2~"from Maximum Model Runs"),
       color = 'Legend') + 
  scale_color_manual(values = colors) + 
  theme_bw() 



