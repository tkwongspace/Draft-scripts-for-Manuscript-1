# M1 Wet Season Visualization
#
# This script visualizes seasonal dynamics of environment in FMNR

library(GGally)
library(tidyr)

# ---- Rearrange data ----
wetTbl = climBase %>% 
  mutate_at(names(climBase)[-1], range01()) %>%
  gather(key = 'Var', value = 'Value', AVGT:TP) %>%
  mutate(Month = month(Date)) %>%
  group_by(Month, Var) %>%
  summarise(Env = mean(Value), SE = sd(Value) / sqrt(n())) %>%
  mutate(upper = Env + SE, lower = Env - SE,
         Var = factor(Var, 
                      levels = c('AVGT', 'PRCP', 'SL', 'SAL', 'TN', 'TP'),
                      labels = c('Temperature', 'Rainfall', 'Mean Sea Level',
                                 'Salinity', 'Nitrogen', 'Phosphorus'))) %>%
  as.data.frame()

# ---- Visualization ----
ggbarplot(wetTbl %>% subset(Var == 'Rainfall'),
          x = 'Month', y = 'Env', fill = 'Var', color = base.Colors[['Blue']]) +
  geom_xspline(data = wetTbl %>% subset(Var != 'Rainfall'),
               aes(color = Var), size = 1) +
  scale_fill_manual(values = base.Colors[['White']]) +
  labs(x = 'Month', 
       y = 'Environment') +
  theme(text = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 5)),
        legend.position = 'top',
        legend.title = element_blank())