# M1 Visualization of CCM Results
#
# This script summarizes data to see if CCM performance is significant 
# at the largest sampling size

# CCM directory
ccmLocation = ""

# Data to plot
sp = ""
lt = ""

# Arrange table of prediction skills from CCM
skillTbl = lapply(levels(ltc.match$species), function(sp){
  byLT = lapply(c('FDM', 'LDM'), function(lt){
    # Target traps
    trapIDs = ltc.match[ltc.match$species == sp,]$trapID
    # CCM Data loading
    ccmList = lapply(trapIDs, function(tID){
      fileName = paste('CCM_', tID, '_', lt, '.Rdata', sep = '')
      load(paste(ccmLocation, fileName, sep = ''))
      return(ccmPerform)
    })
    
    # Extract and Join CCM performance table
    ccmTbl = lapply(1:length(ccmList), function(i){
      ccmList[[i]][['tbl']]
    }) %>% do.call(rbind.data.frame, args = .)
    
    # Summary variables at maximum library size & their significance
    sigTbl = ccmTbl %>% subset(Size == unique(ccmTbl$Size)[length(unique(ccmTbl$Size))]) %>%
      nest_by(Var) %>%
      mutate(sig = list(unique(data$sigEBI))) %>%
      summarise(state = ifelse(length(sig) == 1 & sig == 'NS', 'NS',
                               ifelse(length(sig) != 1 & 'NS' %in% sig, 'Partly', 'Sig'))) %>%
      unique() %>%
      as.data.frame() %>%
      mutate(species = sp, litter = lt)
    
    # Join significance with rho and its SE
    ## Marking of significance in plots 
    ##  - [Solid round (#16) in BLACK] All traps are significant
    ##  - [Solid triangle (#17) in BLACK] Part of traps is/are significant
    ##  - [Cross (#4) in GREY] None of traps are significant
    rhoTbl = ccmTbl %>% subset(Size == unique(ccmTbl$Size)[length(unique(ccmTbl$Size))]) %>%
      group_by(Var) %>% 
      summarise(Skill = mean(Rho),
                SE = sd(Rho)/sqrt(n())) %>% as.data.frame() %>%
      left_join(sigTbl, by = 'Var') %>%
      mutate(species = factor(species, 
                              levels = c('82a Am', '82a Ko', '27a Ko',
                                         '27a Sa', '25a Sc')),
             Var = factor(Var, levels = c('AVGT', 'PRCP', 'SAL',
                                          'SL', 'MinLev', 'MaxLev',
                                          'NO3', 'TP'),
                          labels = c('Temperature', 'Rainfall', 'Salinity',
                                     'Sea Level', 'Min Sea Level', 'Max Sea Level',
                                     'Nitrate', 'Phosphorus')),
             upper = Skill + SE, lower = Skill - SE,
             pointShape = ifelse(state == 'Sig', 16,
                                 ifelse(state == 'NS', 4, 17)),
             pointCol = ifelse(state == 'NS', 
                               base.Colors[['LightGrey']],
                               base.Colors[['Black']]))
    return(rhoTbl)
  }) %>% do.call(rbind.data.frame, args = .)
}) %>% do.call(rbind.data.frame, args = .)

# Visualization

ggplot(skillTbl, aes(x = Var)) +
  geom_hline(yintercept = 0.5, color = base.Colors[['LightGrey']],
             linetype = 'dashed', alpha = .5) +
  geom_linerange(aes(ymin = lower, ymax = upper, color = pointCol)) +
  geom_point(aes(y = Skill, shape = pointShape, color = pointCol),
             size = 3) +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_color_identity() +
  scale_shape_identity() +
  labs(x = '',
       y = 'Prediction Skill') +
  coord_flip() +
  facet_grid(cols = vars(species),
             rows = vars(litter)) +
  theme_classic() +
  theme(text = element_text(size = 16),
        panel.grid.major.y = element_line(linetype = 'dotted'),
        axis.title.y = element_blank())
