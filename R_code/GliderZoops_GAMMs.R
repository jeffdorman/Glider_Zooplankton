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
library(gridExtra)
library(suncalc)
library(lubridate)
library(lattice)
library(gratia)


#open data
zooscan <- read_csv('/Users/gammonkoval/Documents/GitHub/gammon_sol/IRA_Krill/data/ZooScan_Model_Data.csv')
glider <- read_csv('/Users/gammonkoval/Documents/GitHub/gammon_sol/IRA_Krill/data/Glider_Model_Data_newInshoreOffshore_standardized.csv', na = 'NaN')
#group taxa as needed
zooscan$Copepod <- zooscan[, 15] + zooscan[, 17] + zooscan[, 19] + zooscan[, 21] + zooscan[, 23] + zooscan[, 25]
zooscan$Krill <- zooscan[, 37]
zooscan$Gelatinous <- zooscan[, 47] + zooscan[, 51]

ggplot(zooscan, aes(x = sqrt(bathymetry), y = total_biomass_log)) + 
  geom_point() + 
  geom_smooth()

zooscan$Inshore_Offshore <- ifelse((zooscan$line == 80) & (zooscan$longitude < -121), 
                                   'Offshore', 'Inshore')
zooscan$Inshore_Offshore <- ifelse((zooscan$line == 90) & (zooscan$longitude < -119.1),
                                   'Offshore', zooscan$Inshore_Offshore)
zooscan$Inshore_Offshore <- ifelse((zooscan$line == 90) & (zooscan$longitude >= -119.1),
                                   'Inshore', zooscan$Inshore_Offshore)

#clean data
zooscan <- zooscan %>% 
  drop_na(season) %>% 
  rename(lat = latitude) %>% 
  dplyr::select(c("line", "season", "Inshore_Offshore", "total_biomass", 
                  "total_biomass_log", "year", "month",
                  "Day_Night","Copepod","Krill","Gelatinous")) %>% 
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

# #Step 2 assign Day or Night given Sunset/Sunrise at Los Angeles
# glider <- glider %>%
#   rowwise() %>%
#   mutate(
#     sun_times = list(getSunlightTimes(
#       date = as.Date(local_date), 
#       lat = lat, 
#       lon = lon, 
#       keep = c("sunrise", "sunset"), 
#       tz = "America/Los_Angeles"
#     )),
#     sunrise = sun_times$sunrise,
#     sunset = sun_times$sunset,
#     Day_Night = if_else(time_local >= sunrise & time_local < sunset, "Day", "Night")
#   ) %>%
#   ungroup()

#Step3 group and produce summary df.
glider <- glider %>%
  dplyr::select(c("line", "season", "Inshore_Offshore", "big_sv", "little_sv",
                  "year", "month", "Day_Night")) %>% 
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

ggplot(data, aes(x = year)) + 
  geom_point(aes(y = big_sv), color = 'red') + 
  geom_smooth(aes(y = big_sv), color = 'red') + 
  geom_point(aes(x = year, y = total_biomass_log*20), color = 'blue') +
  geom_smooth(aes(x = year, y = total_biomass_log* 20), color = 'blue') +
  scale_y_continuous(name = "Backscatter (dB)",
                     sec.axis = sec_axis( trans=~./20, name="Total Biomass")) 


data %>% 
  drop_na(total_biomass_log) %>% 
  mutate(line = as.factor(line)) %>% 
  ggplot(., aes(x = big_sv, y = total_biomass_log, group = line)) + 
    geom_point() + 
    geom_smooth(aes(color = line), method = 'lm', se = FALSE, linetype = 'dashed') + 
    facet_wrap(.~year)
  

#clean joined data
# data$line <- as.factor(data$line)
data$season <- factor(data$season, levels = c('Winter', 'Spring', 'Summer', 'Fall'))
data$Inshore_Offshore <- as.factor(data$Inshore_Offshore)
data$`Day_Night` <- as.factor(data$`Day_Night`)
data$`year` <- as.factor(data$`year`)
data <- data %>% filter(season!="NA") %>% filter(big_sv!="NaN") %>% filter(total_biomass_log!="NA") %>% droplevels()


f.sbs <- gam(total_biomass_log ~ s(big_sv, by = season) + 
               s(big_sv, by = Inshore_Offshore) + s(big_sv),
             method = 'REML',
             data = data)

draw(f.sbs, residuals = TRUE, type = "response")

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

