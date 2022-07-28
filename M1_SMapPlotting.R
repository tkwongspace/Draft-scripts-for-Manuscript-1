# M1 S-Map Plotting
#
# This script plots sensitivity of mangrove assemblages to environmental changes
# 
# NOTE: This script is after 'M1_SMap.R'

# # Directory of result from S-Map
# resultDir = ""
# load(file = resultDir)
# uvsmap.scneExp = uvsmap$scne

# ---- Summary by Factors ----

plotTbl.single = uvsmap.scneExp %>% 
  subset(var %in% c('AVGT', 'PRCP', 'SL', 'SAL', 'TN', 'TP') & !is.na(sens)) %>%
  group_by(species, litter, var) %>% 
  summarise(avg = mean(sens), sd = sd(sens), count = n(),
            margin = qt(0.975, df = count-1)*sd/sqrt(count)) %>% 
  as.data.frame() %>% 
  mutate(upper = avg + margin, lower = avg - margin,
         var = factor(var,
                      levels = c('AVGT', 'PRCP', 'SL', 'SAL', 'TN', 'TP'),
                      labels = c('Temperature', 'Rainfall', 'Sea Level',
                                 'Salinity', 'Nitrogen', 'Phosphorus')),
         species = factor(species, 
                          levels = c('25a Sc', '27a Sa',
                                     '27a Ko', '82a Ko', '82a Am'),
                          labels = c('26a Sc', '28a Sa',
                                     '28a Ko', '83a Ko', '83a Am')),
         litter = factor(litter,
                         levels = c('LDM', 'FDM'),
                         labels = c('Vegetative Growth', 'Reproduction')))

## Reproduction
ggplot(plotTbl.single %>% subset(litter == 'Reproduction'), 
       aes(x = species, y = avg, color = species)) +
  geom_hline(yintercept = 0, 
             color = base.Colors[['LightGrey']],
             linetype = 'dashed') +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  scale_y_continuous(limits = c(-50, 40)) +
  scale_x_discrete(limits = rev) +
  labs(x = 'Mangrove Assemblages', 
       y = 'Reproductive Sensitivity') +
  facet_wrap(~var, nrow = 2, ncol = 3, scales = 'free_x') +
  coord_flip() +
  theme_classic() +
  theme(text = element_text(size = 16),
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        legend.position = 'none')

## Vegetative growth (Figure S5)
ggplot(plotTbl.single %>% subset(litter == 'Vegetative Growth'), 
       aes(x = species, y = avg, color = species)) +
  geom_hline(yintercept = 0, 
             color = base.Colors[['LightGrey']],
             linetype = 'dashed') +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  scale_y_continuous(limits = c(-63, 25)) +
  scale_x_discrete(limits = rev) +
  labs(x = 'Mangrove Assemblages', 
       y = 'Vegetative Sensitivity') +
  facet_wrap(~var, nrow = 2, ncol = 2, scales = 'free_x') +
  coord_flip() +
  theme_classic() +
  theme(text = element_text(size = 16),
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        legend.position = 'none')

# ---- Main Plots Preparations ----
#
# For all plots about seasonal sensitivity across assemblage
# Each assemblage is assigned to a fixed color in base.Colors (from M1_Preparation.R)
# 26a Sc - Blue
# 28a Sa - Green
# 28a Ko - LightGreen 
# 83a Ko - Violet
# 83a Am - Pink

# Legend for all plots
pLegend = ggplot(sensSum %>% subset(litter == 'LDM' & var == 'MaxLev'),
                 aes(x = Month, y = avg, color = species)) +
  geom_line(aes(y = avg), alpha = .5, linetype = 'dotted') + 
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  geom_point(aes(y = avg, shape = species), fill = 'white', size = 2) +
  scale_color_manual(values = c(base.Colors[['Blue']],
                                base.Colors[['Green']],
                                base.Colors[['LightGreen']],
                                base.Colors[['Violet']],
                                base.Colors[['Pink']])) +
  scale_shape_manual(values = c(24, 24, 21, 21, 21)) +
  labs(color = 'Assemblage',
       shape = 'Assemblage') +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.position = 'bottom') 
colLegend = get_legend(pLegend)

# ---- >> Reproduction ----

# Figure 4
{
  # P1 - Temperature | P2 - Rainfall | P3 - Salinity
  lt = 'FDM'
  
  p1 = ggplot(sensSum %>% subset(litter == lt & var == 'AVGT') %>%
                mutate(var = 'Temperature'),
              aes(x = Month, color = species, shape = species)) + 
    geom_hline(yintercept = 0, color = base.Colors[['LightGrey']],
               linetype = 'dashed') +
    geom_line(aes(y = avg), position = position_dodge(width = 0.5), 
              alpha = .5, linetype = 'dotted') + 
    geom_linerange(aes(ymin = lower, ymax = upper),
                   position = position_dodge(width = 0.5)) + 
    geom_point(aes(y = avg), size = 2,
               fill = 'white',
               position = position_dodge(width = 0.5)) +
    scale_x_continuous(breaks = 1:12) +
    scale_color_manual(values = c(base.Colors[['Blue']],
                                  base.Colors[['LightGreen']],
                                  base.Colors[['Violet']],
                                  base.Colors[['Pink']])) +
    scale_shape_manual(values = c(24, 21, 21, 21)) +
    labs(color = 'Assemblage',
         shape = 'Assemblage') +
    facet_grid(~var) +
    theme_classic() +
    theme(text = element_text(size = 16),
          axis.title = element_blank(),
          strip.text = element_text(size = 16),
          strip.background = element_blank(),
          legend.position = 'top')
  p2 = ggplot(sensSum %>% subset(litter == lt & var == 'PRCP') %>%
                mutate(var = 'Rainfall'),
              aes(x = Month, color = species, shape = species)) + 
    geom_hline(yintercept = 0, color = base.Colors[['LightGrey']],
               linetype = 'dashed') +
    geom_line(aes(y = avg), position = position_dodge(width = 0.5), 
              alpha = .5, linetype = 'dotted') + 
    geom_linerange(aes(ymin = lower, ymax = upper),
                   position = position_dodge(width = 0.5)) + 
    geom_point(aes(y = avg), fill = 'white', size = 2,
               position = position_dodge(width = 0.5)) +
    scale_x_continuous(breaks = 1:12) +
    scale_color_manual(values = c(base.Colors[['LightGreen']],
                                  base.Colors[['Violet']],
                                  base.Colors[['Pink']])) +
    scale_shape_manual(values = c(21, 21, 21)) +
    labs(color = 'Assemblage',
         shape = 'Assemblage') +
    facet_grid(~var) +
    theme_classic() +
    theme(text = element_text(size = 16),
          axis.title = element_blank(),
          strip.text = element_text(size = 16),
          strip.background = element_blank(),
          legend.position = 'top')
  p3 = ggplot(sensSum %>% subset(litter == lt & var == 'SAL') %>%
                mutate(var = 'Salinity'),
              aes(x = Month, color = species, shape = species)) + 
    geom_hline(yintercept = 0, color = base.Colors[['LightGrey']],
               linetype = 'dashed') +
    geom_line(aes(y = avg), position = position_dodge(width = 0.5), 
              alpha = .5, linetype = 'dotted') + 
    geom_linerange(aes(ymin = lower, ymax = upper),
                   position = position_dodge(width = 0.5)) + 
    geom_point(aes(y = avg), fill = 'white', size = 2,
               position = position_dodge(width = 0.5)) +
    scale_x_continuous(breaks = 1:12) +
    scale_color_manual(values = c(base.Colors[['Blue']],
                                  base.Colors[['LightGreen']],
                                  base.Colors[['Violet']],
                                  base.Colors[['Pink']])) +
    scale_shape_manual(values = c(24, 21, 21, 21)) +
    labs(color = 'Assemblage',
         shape = 'Assemblage') +
    facet_grid(~var) +
    theme_classic() +
    theme(text = element_text(size = 16),
          axis.title = element_blank(),
          strip.text = element_text(size = 16),
          strip.background = element_blank(),
          legend.position = 'top')
  
  ggarrange(p1, p2, p3,
            nrow = 3, ncol = 1,
            labels = letters[1:3],
            common.legend = T) %>%
    annotate_figure(bottom = text_grob('Month', size = 16),
                    left = text_grob(expression(atop('Reproduction Sensitivity',
                                                     '(Standarized Change in Biomass Partitioning)')),
                                     size = 16, rot = 90))
}

# ---- >> Vegetative Growth ----

# Figure 3
{
  # P1 - Temperature | P2 - Rainfall | P3 - Salinity
  lt = 'LDM'
  
  p1 = ggplot(sensSum %>% subset(litter == lt & var == 'AVGT') %>%
                mutate(var = 'Temperature'),
              aes(x = Month, color = species, shape = species)) + 
    geom_hline(yintercept = 0, color = base.Colors[['LightGrey']],
               linetype = 'dashed') +
    geom_line(aes(y = avg), position = position_dodge(width = 0.5), 
              alpha = .5, linetype = 'dotted') + 
    geom_linerange(aes(ymin = lower, ymax = upper),
                   position = position_dodge(width = 0.5)) + 
    geom_point(aes(y = avg), fill = 'white', size = 2,
               position = position_dodge(width = 0.5)) +
    scale_x_continuous(breaks = 1:12) +
    scale_color_manual(values = c(base.Colors[['Green']],
                                  base.Colors[['LightGreen']],
                                  base.Colors[['Violet']],
                                  base.Colors[['Pink']])) +
    scale_shape_manual(values = c(24, 21, 21, 21)) +
    labs(color = 'Assemblages',
         shape = 'Assemblages') +
    facet_grid(~var) +
    theme_classic() +
    theme(text = element_text(size = 16),
          axis.title = element_blank(),
          strip.text = element_text(size = 16),
          strip.background = element_blank(),
          legend.position = 'bottom')
  p2 = ggplot(sensSum %>% subset(litter == lt & var == 'PRCP') %>%
                mutate(var = 'Rainfall'),
              aes(x = Month, color = species, shape = species)) + 
    geom_hline(yintercept = 0, color = base.Colors[['LightGrey']],
               linetype = 'dashed') +
    geom_line(aes(y = avg), position = position_dodge(width = 0.5), 
              alpha = .5, linetype = 'dotted') + 
    geom_linerange(aes(ymin = lower, ymax = upper),
                   position = position_dodge(width = 0.5)) + 
    geom_point(aes(y = avg), fill = 'white', size = 2,
               position = position_dodge(width = 0.5)) +
    scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(breaks = c(-150, -100, -50, 0, 50)) +
    scale_color_manual(values = c(base.Colors[['Green']],
                                  base.Colors[['Violet']],
                                  base.Colors[['Pink']])) +
    scale_shape_manual(values = c(24, 21, 21)) +
    labs(color = 'Assemblages',
         shape = 'Assemblages') +
    facet_wrap(~var) +
    theme_classic() +
    theme(text = element_text(size = 16),
          axis.title = element_blank(),
          strip.text = element_text(size = 16),
          strip.background = element_blank(),
          legend.position = 'bottom')
  
  p3 = ggplot(sensSum %>% subset(litter == lt & var == 'SAL') %>%
                mutate(var = 'Salinity'),
              aes(x = Month, color = species, shape = species)) + 
    geom_hline(yintercept = 0, color = base.Colors[['LightGrey']],
               linetype = 'dashed') +
    geom_line(aes(y = avg), position = position_dodge(width = 0.5), 
              alpha = .5, linetype = 'dotted') + 
    geom_linerange(aes(ymin = lower, ymax = upper),
                   position = position_dodge(width = 0.5)) + 
    geom_point(aes(y = avg), fill = 'white', size = 2,
               position = position_dodge(width = 0.5)) +
    scale_x_continuous(breaks = 1:12) +
    scale_color_manual(values = c(base.Colors[['Green']],
                                  base.Colors[['Violet']],
                                  base.Colors[['Pink']])) +
    scale_shape_manual(values = c(24, 21, 21)) +
    labs(color = 'Assemblages',
         shape = 'Assemblages') +
    facet_grid(~var) +
    theme_classic() +
    theme(text = element_text(size = 16),
          axis.title = element_blank(),
          strip.text = element_text(size = 16),
          strip.background = element_blank(),
          legend.position = 'bottom')
  
  ggarrange(p1, p2, p3, 
            nrow = 3, ncol = 1,
            labels = c(letters[1:3]),
            common.legend = T) %>%
    annotate_figure(bottom = text_grob('Month', size = 16),
                    left = text_grob(expression(atop('Sensitivity of Vegetative Growth',
                                                     '(Standarized Change in Biomass Partitioning)')),
                                     size = 16, rot = 90))
}
