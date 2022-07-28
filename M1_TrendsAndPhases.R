# M1 Trends and Phases
#
# This script proceeds estimations of temporal trends at monthly scale
# and a comparison across mangrove assemblages on annual production
# and visualization of pheno-phases

# ---- Temporal Trends ----

library(runner)
library(nlme)
library(GGally)

# Litter biomass differentials
edm.dmTbl = lapply(unique(tp$trapID), function(tr){
  byLT = lapply(c('FDM', 'LDM'), function(lt){
    dfBase = data.frame(
      time = subset(tp, trapID == tr & litter == lt)$date,
      DM = c(NA, diff(subset(tp, trapID == tr & litter == lt)$DM))
    )
    dfBase = dfBase[-1,]
    return(dfBase)
  }); names(byLT) = c('FDM', 'LDM')
  return(byLT)
}); names(edm.dmTbl) = unique(tp$trapID)

## NOTE:
## For stationary tests on differentials refers to 'M1_StationarityTest.R"

# Moving windows
# Slice time series with a 36-month-long moving window

## Litter
mw.sp = lapply(ltc.match$trapID, function(ID){
  res = lapply(c("FDM", "LDM"), function(lt){
    lb = edm.dmTbl[[as.character(ID)]][[lt]][['DM']]
    res = runner(lb, f = c, k = 36)[-c(1:35)]
    return(res)
  }); names(res) = c("FDM", "LDM"); return(res)
}); names(mw.sp) = ltc.match$trapID

## Calculate running mean & PV along with time
litRunning = applyCal(slices = mw.sp, 
                      fun = c('mean', 'PV'), 
                      type = 'litter',
                      version = 'dif') %>% 
  mutate(date = ym('1999-02') %m+% months(date + 18)) %>%
  left_join(ltc.match, by = 'trapID')

## Temporal trends for mean
litMeanTrend = litRunning %>% subset(param == 'mean') %>%
  mutate(time = scale(date)) %>%
  nest_by(species, litter) %>%
  mutate(mod = list(lme(value~time, data = data, 
                        random = ~1|trapID, method = 'ML',
                        correlation = corAR1()))) %>%
  summarise(slope = summary(mod)$coefficients$fixed[[2]],
            p = round(anova(mod)$`p-value`[2], digits = 4)) %>%
  mutate(sig = ifelse(p < 0.05, 'solid', 'dashed'),
         sigMark = ifelse(p>0.05, 'NS',
                          ifelse(p>0.01, '*', 
                                 ifelse(p>0.001, '**', '***'))))
## Temporal trends for PV
litPVtrend = litRunning %>% subset(param == 'PV') %>%
  mutate(time = scale(date)) %>%
  nest_by(species, litter) %>%
  mutate(mod = list(lme(value~time, data = data, 
                        random = ~1|trapID, method = 'ML',
                        correlation = corAR1()))) %>%
  summarise(slope = summary(mod)$coefficients$fixed[[2]],
            p = round(anova(mod)$`p-value`[2], digits = 4)) %>%
  mutate(sig = ifelse(p < 0.05, 'solid', 'dashed'),
         sigMark = ifelse(p>0.05, 'NS',
                          ifelse(p>0.01, '*', 
                                 ifelse(p>0.001, '**', '***'))))

# Environment
mw.wt = lapply(names(climBase)[-1], function(v){
  res = runner(climBase[[v]], f = c, k = 36)[-c(1:35)]
}); names(mw.wt) = names(climBase)[-1]

## Calculate running mean & PV along with time
clmRunning = applyCal(slices = mw.wt, fun = c('mean', 'PV'), type = 'env')  %>% 
  mutate(date = ym('1999-01') %m+% months(date + 18))
## Temporal trends for mean
clmMeanTrend = clmRunning %>% subset(param == 'mean') %>%
  mutate(time = scale(date)) %>%
  nest_by(var) %>%
  mutate(mod = list(gls(value ~ time, data = data,
                        correlation = corAR1()))) %>%
  summarise(slope = summary(mod)$tTable[2],
            p = round(summary(mod)$tTable[8], digits = 4)) %>%
  mutate(sig = ifelse(p < 0.05, 'solid', 'dashed'),
         sigMark = ifelse(p>0.05, 'NS',
                          ifelse(p>0.01, '*', 
                                 ifelse(p>0.001, '**', '***'))))
## Temporal trends for PV
clmPVTrend = clmRunning %>% subset(param == 'PV') %>%
  mutate(time = scale(date)) %>%
  nest_by(var) %>%
  mutate(mod = list(gls(value ~ time, data = data,
                        correlation = corAR1()))) %>%
  summarise(slope = summary(mod)$tTable[2],
            p = round(summary(mod)$tTable[8], digits = 4)) %>%
  mutate(sig = ifelse(p < 0.05, 'solid', 'dashed'),
         sigMark = ifelse(p>0.05, 'NS',
                          ifelse(p>0.01, '*', 
                                 ifelse(p>0.001, '**', '***'))))

## 20-year PV for each corresponding calendar month
clmMonthly = climBase %>% 
  tidyr::gather(key = 'Var', value = 'Value', AVGT:PO4) %>%
  mutate(Month = month(Date)) %>%
  group_by(Month, Var) %>%
  summarise(PV = getPV(Value)) %>%
  select(Var, Month, PV)

## Visualization (Figure S3)
{
  plotTbl.envPV = clmMonthly %>% 
    subset(Var %in% c('AVGT', 'PRCP', 'SL', 'SAL', 'TN', 'TP')) %>%
    mutate(Var = factor(Var, 
                        levels = c('AVGT', 'PRCP', 'SL', 
                                   'SAL', 'TN', 'TP'),
                        labels = c('Temperature', 'Rainfall', 
                                   'Mean Sea Level', 'Salinity',
                                   'Nitrogen', 'Phosphorus'))) %>%
    as.data.frame()
  ## P1-6: AVGT, PRCP, SL, SAL, N, P
  varList.envPV = c('Temperature', 'Rainfall','Mean Sea Level', 'Salinity',
                    'Nitrogen', 'Phosphorus')
  plotList.envPV = lapply(varList.envPV, function(v){
    p = ggplot(plotTbl.envPV %>% subset(Var == v),
               aes(x = Month, y = PV)) +
      geom_col(fill = NA, color = base.Colors[['Black']], width = .5) +
      scale_x_continuous(breaks = 1:12) +
      scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
      facet_wrap(~Var) +
      theme_classic() +
      theme(text = element_text(size = 16),
            strip.background = element_blank(),
            strip.text = element_text(size = 16),
            axis.title = element_blank())
    return(p)
  })
  ggarrange(plotList.envPV[[1]], plotList.envPV[[2]], plotList.envPV[[3]],
            plotList.envPV[[4]], plotList.envPV[[5]], plotList.envPV[[6]],
            labels = letters[1:6]) %>%
    annotate_figure(left = text_grob('Proportional Variation',
                                     rot = 90, size = 16),
                    bottom = text_grob('Month', size = 16))
}

# ---- Annual Production ----

## Summary on annual production
## Sum by years
litAnnual = tp %>% mutate(Year = year(date)) %>%
  group_by(species, trapID, litter, Year) %>%
  summarise(annSum = sum(DM) * 10) %>%  # x10 to convert unit from g/m2 -> kg/ha
  as.data.frame()
## 20-year mean with SE
litAnnSummary = litAnnual %>%
  group_by(species, litter) %>%
  summarise(biomass = mean(annSum), SE = sd(annSum)/sqrt(n())) %>%
  as.data.frame()
## 20-year PV
litAnnPV = tp %>% mutate(Year = year(date)) %>%
  group_by(species, trapID, litter, Year) %>%
  summarise(PV = getPV(DM)) %>%
  as.data.frame()

## T-test (or Welch's ANOVA)
compare_means(annSum~species, data = litAnnual,
              group.by = 'litter', method = 't.test')
compare_means(PV~species, data = litAnnPV, 
              group.by = 'litter', method = 't.test')

## Visualization (Figure S1)
{
  ## P1/2 - Annual production of vegetation / reproduction
  plotTbl.litAnn = litAnnual %>%
    mutate(species = factor(species, 
                            levels = c('27a Sa', '25a Sc', '27a Ko',
                                       '82a Ko', '82a Am'),
                            labels = c('28a Sa', '26a Sc', '28a Ko',
                                       '83a Ko', '83a Am')),
           litter = factor(litter, 
                           levels = c('LDM', 'FDM'),
                           labels = c('Vegetative Growth',
                                      'Reproduction')))
  ## P1
  p1 = ggplot(plotTbl.litAnn %>% subset(litter == 'Vegetative Growth'), 
              aes(x = species, y = annSum)) +
    geom_jitter(shape = 20, width = .2, 
                color = base.Colors[['Green']], alpha = .5) +
    geom_boxplot(width = .5, fill = NA, color = base.Colors[['Black']]) +
    geom_text(data = data.frame(
      species = factor(x = c('28a Sa', '26a Sc', '28a Ko',
                             '83a Ko', '83a Am'),
                       levels = c('28a Sa', '26a Sc', '28a Ko',
                                  '83a Ko', '83a Am')),
      sigMark = c('a', 'b', 'c', 'd', 'c')
    ), aes(label = sigMark), y = 19000, size = 5) +
    scale_y_continuous(limits = c(0, 20000), breaks = seq(0, 15000, 5000),
                       labels = c('0', '5', '10', '15')) +
    labs(x = '', y = '') +
    facet_wrap(~litter) +
    theme_classic() +
    theme(text = element_text(size = 16),
          strip.background = element_blank(),
          strip.text = element_text(size = 16),
          axis.title.x = element_blank())
  p2 = ggplot(plotTbl.litAnn %>% subset(litter == 'Reproduction'), 
              aes(x = species, y = annSum)) +
    geom_jitter(shape = 20, width = .2, 
                color = base.Colors[['Green']], alpha = .5) +
    geom_boxplot(width = .5, fill = NA, color = base.Colors[['Black']]) +
    geom_text(data = data.frame(
      species = factor(x = c('28a Sa', '26a Sc', '28a Ko',
                             '83a Ko', '83a Am'),
                       levels = c('28a Sa', '26a Sc', '28a Ko',
                                  '83a Ko', '83a Am')),
      sigMark = c('bc', 'b', 'c', 'a', 'd')
    ), aes(label = sigMark), y = 19000, size = 5) +
    scale_y_continuous(limits = c(0, 20000), breaks = seq(0, 15000, 5000),
                       labels = c('0', '5', '10', '15')) +
    labs(x = '', y = '') +
    facet_wrap(~litter) +
    theme_classic() +
    theme(text = element_text(size = 16),
          strip.background = element_blank(),
          strip.text = element_text(size = 16),
          axis.title.x = element_blank())
  ## P3/4 - PV of vegetative biomass / reproduction
  plotTbl.PV = litAnnPV %>%
    mutate(species = factor(species, 
                            levels = c('27a Sa', '25a Sc', '27a Ko',
                                       '82a Ko', '82a Am'),
                            labels = c('28a Sa', '26a Sc', '28a Ko',
                                       '83a Ko', '83a Am')),
           litter = factor(litter, 
                           levels = c('LDM', 'FDM'),
                           labels = c('Vegetative Growth',
                                      'Reproduction')))
  p3 = ggplot(plotTbl.PV %>% subset(litter == 'Vegetative Growth'), 
              aes(x = species, y = PV)) +
    geom_jitter(shape = 20, width = .2, 
                color = base.Colors[['Green']], alpha = .5) +
    geom_boxplot(width = .5, fill = NA, color = base.Colors[['Black']]) +
    geom_text(data = data.frame(
      species = factor(x = c('28a Sa', '26a Sc', '28a Ko',
                             '83a Ko', '83a Am'),
                       levels = c('28a Sa', '26a Sc', '28a Ko',
                                  '83a Ko', '83a Am')),
      sigMark = c('a', 'c', 'c', 'd', 'b')
    ), aes(label = sigMark), y = 1, size = 5) +
    scale_y_continuous(limits = c(0, 1.1),
                       breaks = c(0, 0.5, 1)) +
    labs(x = '', y = '') +
    theme_classic() +
    theme(text = element_text(size = 16),
          axis.title.x = element_blank())
  p4 = ggplot(plotTbl.PV %>% subset(litter == 'Reproduction'), 
              aes(x = species, y = PV)) +
    geom_jitter(shape = 20, width = .2, 
                color = base.Colors[['Green']], alpha = .5) +
    geom_boxplot(width = .5, fill = NA, color = base.Colors[['Black']]) +
    geom_text(data = data.frame(
      species = factor(x = c('28a Sa', '26a Sc', '28a Ko',
                             '83a Ko', '83a Am'),
                       levels = c('28a Sa', '26a Sc', '28a Ko',
                                  '83a Ko', '83a Am')),
      sigMark = c('a', 'b', 'b', 'b', 'c')
    ), aes(label = sigMark), y = 1, size = 5) +
    scale_y_continuous(limits = c(0, 1.1),
                       breaks = c(0, 0.5, 1)) +
    labs(x = '', y = '') +
    theme_classic() +
    theme(text = element_text(size = 16),
          axis.title.x = element_blank())
  ## Arrange all subplots
  ggarrange(p1, p2, p3, p4,
            labels = letters[1:4],
            nrow = 2, ncol = 2,
            heights = c(0.53, 0.47)) %>%
    annotate_figure(left = text_grob(paste("Proportional Variation",
                                           "Annual Biomass (Mg/ha)",
                                           sep = "           "),
                                     rot = 90,
                                     size = 16),
                    bottom = text_grob("Assemblages", size = 16))
}

# ---- Seasonal Phases ----

# Phase table
phaseSubTbl = tp %>% 
  group_by(trapID, litter) %>%
  mutate(DM = range01(DM),
         Month = month(date)) %>%
  group_by(Month, species, trapID, litter) %>%
  summarise(DM = mean(DM)) %>%
  group_by(Month, species, litter) %>%
  summarise(Mean = mean(DM), SE = sd(DM)/sqrt(n())) %>%
  mutate(upper = Mean + SE, lower = Mean - SE) %>%
  as.data.frame()
phaseTbl = full_join(phaseSubTbl %>% subset(litter == 'FDM') %>%
                       select(species, Month, Mean, upper, lower),
                     phaseSubTbl %>% subset(litter == 'LDM') %>%
                       select(species, Month, Mean, upper, lower),
                     by = c('species', 'Month'),
                     suffix = c('F', 'L')) %>%
  select(species, Month, MeanF, upperF, lowerF, MeanL, upperL, lowerL) %>%
  mutate(Season = ifelse(Month %in% 3:5, 'Spring',
                         ifelse(Month %in% 6:8, 'Summer',
                                ifelse(Month %in% 9:11, 'Autumn', 'Winter'))),
         Season = factor(Season, 
                         levels = c('Spring', 'Summer', 
                                    'Autumn', 'Winter'))) %>%
  arrange(Month)

# Visualization (Figure 2)
{
  ## Create list of subplots
  plotList = lapply(1:6, function(SPID){
    if (SPID == 1) {  # 28a Sa
      pointTbl = phaseTbl %>% subset(species == '27a Sa') %>% 
        mutate(species = as.character(species),
               species = 'Exotic: 28a Sa')
      labelTbl = pointTbl %>% subset(Month %in% c(1, 12)) %>%
        mutate(labelX = c(0.0396 + 0.025, 0), 
               labelY = c(0.00926 + 0.015, 0)) %>%
        select(species, Month, labelX, labelY)
      segmentTbl = bind_rows(
        pointTbl %>% subset(Month %in% 3:6) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 3:6),
                 Season = 'Spring'),
        pointTbl %>% subset(Month %in% 6:9) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 6:9),
                 Season = 'Summer'),
        pointTbl %>% subset(Month %in% 9:12) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 9:12),
                 Season = 'Autumn'),
        pointTbl %>% subset(Month %in% c(12, 1, 2, 3)) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = c(12, 1, 2, 3)),
                 Season = 'Winter') %>% arrange(Month)
      ) %>% mutate(Season = factor(Season, 
                                   levels = c('Spring', 'Summer',
                                              'Autumn', 'Winter'),
                                   labels = c('Spring (Mar - May)',
                                              'Summer (Jun - Aug)',
                                              'Autumn (Sep - Nov)',
                                              'Winter (Dec - Feb)')))
    } else if (SPID == 2) {  # 26a Sc
      pointTbl = phaseTbl %>% subset(species == '25a Sc') %>% 
        mutate(species = as.character(species),
               species = 'Exotic: 26a Sc')
      labelTbl = pointTbl %>% subset(Month %in% c(1, 12)) %>%
        mutate(labelX = c(0.299 + 0.015, 0.206 - 0.015), 
               labelY = c(0.139 + 0.015, 0.163 - 0.015)) %>%
        select(species, Month, labelX, labelY)
      segmentTbl = bind_rows(
        pointTbl %>% subset(Month %in% 3:6) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 3:6),
                 Season = 'Spring'),
        pointTbl %>% subset(Month %in% 6:9) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 6:9),
                 Season = 'Summer'),
        pointTbl %>% subset(Month %in% 9:12) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 9:12),
                 Season = 'Autumn'),
        pointTbl %>% subset(Month %in% c(12, 1, 2, 3)) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = c(12, 1, 2, 3)),
                 Season = 'Winter') %>% arrange(Month)
      ) %>% mutate(Season = factor(Season, 
                                   levels = c('Spring', 'Summer',
                                              'Autumn', 'Winter'),
                                   labels = c('Spring (Mar - May)',
                                              'Summer (Jun - Aug)',
                                              'Autumn (Sep - Nov)',
                                              'Winter (Dec - Feb)')))
    } else if (SPID == 3) {  # 28a Ko
      pointTbl = phaseTbl %>% subset(species == '27a Ko') %>% 
        mutate(species = as.character(species),
               species = 'Indigenous: 28a Ko')
      labelTbl = pointTbl %>% subset(Month %in% c(1, 12)) %>%
        mutate(labelX = c(0.117 + 0.015, 0.127 - 0.015), 
               labelY = c(0.046 + 0.01, 0.0234 - 0.015)) %>%
        select(species, Month, labelX, labelY)
      segmentTbl = bind_rows(
        pointTbl %>% subset(Month %in% 3:6) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 3:6),
                 Season = 'Spring'),
        pointTbl %>% subset(Month %in% 6:9) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 6:9),
                 Season = 'Summer'),
        pointTbl %>% subset(Month %in% 9:12) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 9:12),
                 Season = 'Autumn'),
        pointTbl %>% subset(Month %in% c(12, 1, 2, 3)) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = c(12, 1, 2, 3)),
                 Season = 'Winter') %>% arrange(Month)
      ) %>% mutate(Season = factor(Season, 
                                   levels = c('Spring', 'Summer',
                                              'Autumn', 'Winter'),
                                   labels = c('Spring (Mar - May)',
                                              'Summer (Jun - Aug)',
                                              'Autumn (Sep - Nov)',
                                              'Winter (Dec - Feb)')))
    } else if (SPID == 4) {  # 83a Ko
      pointTbl = phaseTbl %>% subset(species == '82a Ko') %>% 
        mutate(species = as.character(species),
               species = 'Indigenous: 83a Ko')
      labelTbl = pointTbl %>% subset(Month %in% c(1, 12)) %>%
        mutate(labelX = c(0.165 + 0.015, 0.154 - 0.015), 
               labelY = c(0.0324 + 0.01, 0.0186 - 0.02)) %>%
        select(species, Month, labelX, labelY)
      segmentTbl = bind_rows(
        pointTbl %>% subset(Month %in% 3:6) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 3:6),
                 Season = 'Spring'),
        pointTbl %>% subset(Month %in% 6:9) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 6:9),
                 Season = 'Summer'),
        pointTbl %>% subset(Month %in% 9:12) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 9:12),
                 Season = 'Autumn'),
        pointTbl %>% subset(Month %in% c(12, 1, 2, 3)) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = c(12, 1, 2, 3)),
                 Season = 'Winter') %>% arrange(Month)
      ) %>% mutate(Season = factor(Season, 
                                   levels = c('Spring', 'Summer',
                                              'Autumn', 'Winter'),
                                   labels = c('Spring (Mar - May)',
                                              'Summer (Jun - Aug)',
                                              'Autumn (Sep - Nov)',
                                              'Winter (Dec - Feb)')))
    } else {  # 83a Am
      pointTbl = phaseTbl %>% subset(species == '82a Am') %>% 
        mutate(species = as.character(species),
               species = 'Indigenous: 83a Am')
      labelTbl = pointTbl %>% subset(Month %in% c(1, 12)) %>%
        mutate(labelX = c(0.0926 + 0.015, 0.1364 + 0.02), 
               labelY = c(0.00755 + 0.005, 0.00337 + 0.01)) %>%
        select(species, Month, labelX, labelY)
      segmentTbl = bind_rows(
        pointTbl %>% subset(Month %in% 3:6) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 3:6),
                 Season = 'Spring'),
        pointTbl %>% subset(Month %in% 6:9) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 6:9),
                 Season = 'Summer'),
        pointTbl %>% subset(Month %in% 9:12) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = 9:12),
                 Season = 'Autumn'),
        pointTbl %>% subset(Month %in% c(12, 1, 2, 3)) %>%
          select(species, Month, MeanF, MeanL) %>%
          mutate(Month = factor(Month, levels = c(12, 1, 2, 3)),
                 Season = 'Winter') %>% arrange(Month)
      ) %>% mutate(Season = factor(Season, 
                                   levels = c('Spring', 'Summer',
                                              'Autumn', 'Winter'),
                                   labels = c('Spring (Mar - May)',
                                              'Summer (Jun - Aug)',
                                              'Autumn (Sep - Nov)',
                                              'Winter (Dec - Feb)')))
    }
    ## PLOTTING
    ## SPID == 6 is specially for legend
    if (SPID %in% 1:5) {
      p = ggplot(pointTbl, aes(x = MeanL, y = MeanF)) + 
        geom_path(data = segmentTbl, aes(color = Season), 
                  size = 1, alpha = .7) + 
        ## Line ranges of F and L
        geom_linerange(aes(xmin = lowerL, xmax = upperL)) +
        geom_linerange(aes(ymin = lowerF, ymax = upperF)) +
        ## Labels
        geom_text(data = labelTbl, aes(label = Month, x = labelX, y = labelY), 
                  check_overlap = T) +
        ## Assign special point marks and colors
        geom_point(data = pointTbl %>% subset(Month == 1),
                   color = base.Colors[['Purple']], size = 3, shape = 12) +
        geom_point(data = pointTbl %>% subset(Month == 12),
                   color = base.Colors[['Purple']], size = 3, shape = 5) +
        scale_color_manual(values = c(base.Colors[['Pink']], 
                                      base.Colors[['Green']], 
                                      base.Colors[['Orange']], 
                                      base.Colors[['LightGrey']])) +
        scale_y_continuous(breaks = seq(0, 0.5, 0.1)) +
        scale_x_continuous(breaks = seq(0, 0.5, 0.1)) +
        facet_wrap(~species) +
        theme_classic() +
        theme(text = element_text(size = 16),
              strip.text = element_text(size = 16),
              strip.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              panel.background = element_rect(color = 'black'),
              legend.position = 'none')
    } else { ## Get legend
      p = ggplot(pointTbl, aes(x = MeanL, y = MeanF)) + 
        geom_path(data = segmentTbl, aes(color = Season), 
                  size = 1, alpha = .7) +
        scale_color_manual(values = c(base.Colors[['Pink']], 
                                      base.Colors[['Green']], 
                                      base.Colors[['Orange']], 
                                      base.Colors[['LightGrey']])) +
        labs(color = 'Seasons') +
        theme_classic() +
        theme(text = element_text(size = 16))
      p = get_legend(p)
    }
    return(p)
  })
  # Get legend from subplots
  ggarrange(plotList[[1]], plotList[[2]], plotList[[3]],
            plotList[[4]], plotList[[5]], plotList[[6]], 
            labels = c(letters[1:5], '')) %>%
    annotate_figure(left = text_grob('Reproduction', size = 16, rot = 90),
                    bottom = text_grob('Vegetative Growth', size = 16))
}
