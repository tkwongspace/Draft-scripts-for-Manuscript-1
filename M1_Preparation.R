# M1 Preparation
#
# This script loads data tables from directory
# and gets everything ready for further analysis
# (e.g., loads packages, fills up missing data, scales data)

# ---- Directory ----
baseDir = ""

# ---- Packages ----
library(dplyr)
library(ggplot2)
library(lubridate)
library(ggpubr)
library(ggalt)
library(zoo)

# ---- Colors ----
base.Colors = list('White' = '#F2F2F2',
                   'LightGrey' = '#9BAFB6',
                   'Blue' = '#2E6171',
                   'SeaBlue' = '#27187E',
                   'Cyan' = '#00798C',
                   'Orange' = '#F79824',
                   'Red' = '#772014',
                   'Purple' = '#42253B',
                   'Violet' = '#A05C7B',
                   'Eggplant' = '#673C4F',
                   'Pink' = '#E26D5A',
                   'Aqua' = '#4CE0B3',
                   'LightGreen' = '#B6C649',
                   'Green' = '#488B49',
                   'DarkLava' = '#524632',
                   'Black' = '#07020D')

# ---- Tools ----
source(paste(baseDir, "M1_Toolbox.R", sep = '/'))

# ---- Data: Litter ----
## Read table
tp0 = read.csv(paste(baseDir, 'DataFiles/litterByTrap.csv', sep = '')) %>%
  mutate(species = paste(Age, "a ", Species, sep = ""),
         litter = paste(Type, "DM", sep = ""),
         date = ym(paste(Year, Month, sep = "-")) %m+%
           months(1) - 1) %>%
  select(date, species, litter, trapID = TrapID, DM, 
         Year, Month) %>%
  subset(species != "28a Sc")
## Fill missing data
tp = lapply(unique(tp0$species), function(sp){
  lapply(unique(tp0$litter), function(lt){
    df = tp0 %>% subset(species == sp & litter == lt)
    res = lapply(unique(df$trapID), function(ID){
      gapCheck = data.frame(date = seq(as.Date('1999-01-01'),
                                       as.Date('2019-10-01'),
                                       by = 'month')) %>%
        mutate(date = date %m+% months(1) - 1) %>%
        left_join(df %>% subset(trapID == ID),
                  by = 'date')
      gapFilled = data.frame(
        date = gapCheck$date,
        species = sp,
        litter = lt,
        trapID = ID,
        DM = na.approx(gapCheck$DM)
      )
      return(gapFilled)
    }) %>% do.call(rbind.data.frame, args = .) 
    return(res)
  }) %>% do.call(rbind.data.frame, args = .)
}) %>% do.call(rbind.data.frame, args = .) %>%
  ltCleaning() %>%
  select(date, species, litter, trapID, DM)
rm(tp0)

# ---- Data: Environment ----
## Temperature and Rainfall (LFS)
weather = read.csv(file = paste(baseDir, "/DataFiles/LFS_weather.csv", sep = "")) %>%
  mutate(date = ymd(paste(year, month, '01', sep = '-')) %m+% months(1) -1) %>%
  select(Date = date,
         AVGT = temp, 
         PRCP = prcp)
wt.rmExt = data.frame(Date = weather$Date, lapply(weather[,-1], rmExtm, lb = .05, ub = .95))
rm(weather)

## Sea Levels (TBT)
tide = read.csv(file = paste(baseDir, "/DataFiles/TBT_tides.csv", sep = "")) %>%
  mutate(date = ymd(paste(year, month, '01', sep = '-')) %m+% months(1) -1) %>%
  select(Date = date,
         SL = msl,
         MinLev = minlevel,
         MaxLev = maxlevel)
tide = data.frame(Date = tide$Date, na.approx(tide[,-1]))
tide.rmExt = data.frame(Date = tide$Date, lapply(tide[,-1], rmExtm, lb = .05, ub = .95))
rm(tide)

## Salinity and Nutrients
wq = read.csv(file = paste(baseDir, "/DataFiles/HKEPD_quality.csv", sep = "")) %>%
  mutate(SampleDate = as_date(Date),
         Date = ymd(paste(Year, Month, '01', sep = '-')) %m+% months(1) -1)
dateRange = data.frame(date = seq(as.Date('1998-01-01'), 
                                  as.Date('2020-10-01'), by = 'month')) %>%
  mutate(Date = date %m+% months(1) - 1) %>% select(Date)
gapFilled = dateRange %>%
  left_join(wq %>% subset(Station == 'DM2'), by = 'Date') %>%
  select(-Station, -Date, -Year, -Month, -SampleDate) %>%
  na.approx()
wqGapCheck = data.frame(Date = dateRange, gapFilled)
wq.rmExt = data.frame(Date = wqGapCheck$Date, 
                      lapply(wqGapCheck %>% select(-Date),
                             rmExtm, lb = .05, ub = .95))
rm(wq, dateRange, gapFilled, wqGapCheck)

## Join all environment
wt = left_join(wt.rmExt, tide.rmExt, by = 'Date') %>% 
  left_join(wq.rmExt, by = 'Date')
rm(wt.rmExt, tide.rmExt, wq.rmExt)

# ---- Data: Arrangement ----
## Extract environmental data in target temporal period
climBase = wt %>% subset(Date >= ym('1999-01') & Date <= ym('2019-11')) %>%
  select(Date, AVGT, PRCP, SL, SAL, TN, TP)

## Scaled version of environment
climScale = climBase %>% 
  mutate_at(c('AVGT', 'PRCP', 'SL', 'SAL', 'TN', 'TP'),
            function(x){as.numeric(scale(x))}) %>% 
  as.data.frame()

## Join litter biomass and scaled environment
ltc.full = tp %>% left_join(climScale, by = c('date' = 'Date')) %>%
  mutate(Year = year(date), Month = month(date))

## Lookup table from traps to species 
ltc.match = ltc.full %>% select(species, trapID) %>% unique()