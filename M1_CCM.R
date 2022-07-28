# M1 CCM
#
# This script proceeds state-space reconstruction and CCM

library(rEDM)
library(hash)

# ---- Preparations ----

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

# ---- Optimization of Embedding Size ----

# Identify optimal E
edm.uw.findE = lapply(unique(ltc.match$trapID), function(tr){
  lapply(c('FDM', 'LDM'), function(lt){
    findE = EmbedDimension(dataFrame = edm.dmTbl[[as.character(tr)]][[lt]],
                           lib = '1 179',  # 1999-01 to 2013-12
                           pred = '180 249',  # 2014-01 to 2019-10
                           maxE = 13,
                           columns = 'DM',
                           target = 'DM',
                           embedded = F,
                           showPlot = F)
    # Add species and litter marks
    findE = mutate(findE, 
                   trap = tr,
                   species = ltc.match[ltc.match$trapID == tr,]$species, 
                   litter = lt)
    return(findE)
  }) %>%
    do.call(rbind.data.frame, args = .)
}) %>%
  do.call(rbind.data.frame, args = .)

# Find rows with max rho in each species (potential E in each traps)
edm.E.candidates = lapply(levels(edm.uw.findE$species), function(sp){
  byLt = lapply(c('FDM', 'LDM'), function(lt){
    Esummary = edm.uw.findE %>%
      group_by(trap, litter) %>%
      mutate(maxRho = max(rho),
             optimal = ifelse(rho == maxRho, T, F)) %>%
      subset(optimal)
    output = unique(subset(Esummary,
                           species == sp & litter == lt)$E)
    return(output)
  }); names(byLt) = c('FDM', 'LDM')
  return(byLt)
}); names(edm.E.candidates) = levels(edm.uw.findE$species)

# Compare and take E when maximum prediction skill is reached
## Simplex models are built on traps and averaged to get rho for the species
edm.E.selection = lapply(levels(edm.uw.findE$species), function(sp){
  byLt = lapply(c('FDM', 'LDM'), function(lt){
    # Find E candidates
    EList = edm.E.candidates[[sp]][[lt]]
    byE = lapply(EList, function(testE){
      # Find traps
      trapList = ltc.match[(ltc.match$species == sp), 2]
      byTrap = lapply(trapList, function(tr){
        # Input to Simplex
        inputdf = edm.dmTbl[[as.character(tr)]][[lt]]
        models = Simplex(dataFrame = inputdf,
                         lib = '1 179',
                         pred = '180 249',
                         E = testE,
                         columns = 'DM',
                         target  = 'DM',
                         embedded = F, 
                         showPlot = F)
        rho = ComputeError(models$Observations, 
                           models$Predictions)$rho
        return(rho)
      }) #byTrap
      
      # Create data frame
      resTbl = data.frame(species = sp,
                          litter = lt,
                          E = testE,
                          rho = mean(unlist(byTrap)),
                          se = sd(unlist(byTrap))/sqrt(length(byTrap)))
      return(resTbl)
    }) %>% do.call(rbind.data.frame, args = .) #byE
  }) %>% do.call(rbind.data.frame, args = .)
}) %>% do.call(rbind.data.frame, args = .) %>%
  group_by(species, litter) %>%
  mutate(maxRho = max(rho),
         optimal = ifelse(rho == maxRho, T, F)) %>%
  ungroup() %>%
  select(species, litter, E, rho, se, optimal)

# Create list for optimal E (easier to further use)
edm.E = lapply(levels(edm.uw.findE$species), function(sp){
  byLT = lapply(c('FDM', 'LDM'), function(lt){
    optimalRow = edm.E.selection %>% 
      subset(optimal & species == sp & litter == lt)
    output = list('E' = optimalRow$E,
                  'rho' = optimalRow$rho)
    return(output)
  }); names(byLT) = c('FDM', 'LDM')
  return(byLT)
}); names(edm.E) = levels(edm.uw.findE$species)

# save(edm.E, file = paste(baseDir, 'DataFiles/edmE_diff.Rdata', sep = ''))
rm(edm.uw.findE, edm.E.selection, edm.E.candidates)

# ---- CCM on Biomass Differentials ----

ccmList = lapply(levels(ltc.match$species), function(sp){
  trapList = ltc.match[ltc.match$species == sp,]$trapID
  byLT = lapply(c('FDM', 'LDM'), function(lt){
    byTrap = lapply(trapList, function(tr){
      message(paste('> Process starts for ', sp, ' ', lt, 
                    ' (trap #', tr, ') at ', Sys.time(), sep = ''))
      ccmPerform = edm.ccm(tr = tr,
                           lt = lt,
                           lib = '36 236 3',
                           nullMod = 'ebi',
                           nullModRep = 50)
      # Save to disk
      fileName = paste('CCM', tr, lt, sep = '_')
      save(ccmPerform, file = paste(baseDir, 'DataFiles/edmE_diff/',
                                    fileName, '.Rdata', sep = ''))
      message(paste('> File output of ', sp, ' ', lt, 
                    ' (trap #', tr, ') finished at ', Sys.time(), '.', sep = ''))
      return(ccmPerform)
    }); names(byTrap) = trapList
    return(byTrap)
  }); names(byLT) = c('FDM', 'LDM')
  return(byLT)
}); names(ccmList) = levels(ltc.match$species)

# Visualization of CCM results please refer to 'M1_CCMPlotting.R'
