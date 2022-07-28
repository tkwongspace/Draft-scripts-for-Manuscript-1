# M1 S-Map on litter biomass ~ each abiotic factor
#
# This script proceeds univariate S-Map regression on litter biomass
# with each environmental factor that is proved significance in CCM
#
# NOTE:
# This script is after 'M1_CCM.R'

# ---- Preparations ----

library(rEDM)
library(ggExtra)

# Factors that are significant (from CCM, Table 2)
edm.Var.base = list('82a Am FDM' = c('AVGT', 'PRCP', 'SL', 'SAL'),
                    '82a Ko FDM' = c('AVGT', 'PRCP', 'SL', 'SAL', 'TN'),
                    '27a Ko FDM' = c('AVGT', 'PRCP', 'SAL', 'TN'),
                    '27a Sa FDM' = c('TN', 'TP'),
                    '25a Sc FDM' = c('AVGT', 'SAL', 'TN'),
                    '82a Am LDM' = c('AVGT', 'PRCP', 'SAL'),
                    '82a Ko LDM' = c('AVGT', 'PRCP', 'SL', 'SAL'),
                    '27a Ko LDM' = c('AVGT'),
                    '27a Sa LDM' = c('AVGT', 'PRCP', 'SL', 'SAL'),
                    '25a Sc LDM' = c('SL'))

# ---- Optimal E ----

# Basic table for S-Map
uvsmap.baseTbl = lapply(names(edm.dmTbl), function(tr){
  lapply(c('FDM', 'LDM'), function(lt){
    edm.dmTbl[[as.character(tr)]][[lt]] %>% 
      mutate(trap = as.numeric(tr), litter = lt,
             species = ltc.match[ltc.match$trapID == as.numeric(tr), ]$species)
  }) %>% do.call(rbind.data.frame, args = .)
}) %>% do.call(rbind.data.frame, args = .) %>%
  left_join(climScale, by = c('time' = 'Date'))

# Look for optimal E by Simplex
uvsmap.eTbl = pbapply::pblapply(levels(uvsmap.baseTbl$species), function(sp){
  lapply(c('FDM', 'LDM'), function(lt){
    # Trap list
    trapList = ltc.match[ltc.match$species == sp,]$trapID
    # Identify variables
    varList = edm.Var.base[[paste(sp, lt)]]
    # Process by trap & variables
    byTrap = lapply(trapList, function(tr){
      byV = lapply(varList, function(targetV){
        # Basic table for E finding (remove useless columns)
        inputBase = uvsmap.baseTbl %>% subset(trap == tr & litter == lt) %>%
          select(time, DM, Env = all_of(targetV))
        eToTest = 1:13
        # Perform Simplex on each E candidate
        findE = lapply(eToTest, function(optE){
          # Embed tables
          if (optE >1) {
            embBase = inputBase %>% select(time, Env)
            inputEmb = Embed(dataFrame = embBase,
                             E = optE,
                             columns = 'Env')
            colInput = c('Env', paste('Env', 1:(optE-1), sep = ''))
            names(inputEmb) = colInput
            inputTbl = inputBase %>% select(time, DM) %>%
              cbind(inputEmb) %>% slice(-(1:(optE-1)))
          } else {
            inputTbl = inputBase
            colInput = 'Env'
          }
          # Simplex
          trainLib = '1 179'; testLib = paste(180, dim(inputTbl)[1], collapse = ' ')
          mod = Simplex(dataFrame = inputTbl,
                        lib = trainLib,
                        pred = testLib,
                        E = optE,
                        columns = paste(colInput, collapse = ' '),
                        target = 'DM',
                        embedded = T,
                        showPlot = F)
          rho = ComputeError(mod$Observations, mod$Predictions)$rho
          resTbl = data.frame(species = sp,
                              litter = lt,
                              trapID = tr,
                              Var = targetV, 
                              E = optE, 
                              rho = rho)
          return(resTbl)
        }) %>% do.call(rbind.data.frame, args = .)
        return(findE)
      }) %>% do.call(rbind.data.frame, args = .) 
    }) %>% do.call(rbind.data.frame, args = .) 
  }) %>% do.call(rbind.data.frame, args = .) 
}) %>% do.call(rbind.data.frame, args = .)

# Best E for each trap
uvsmap.optEbyTrap = uvsmap.eTbl %>%
  group_by(litter, trapID, Var) %>%
  mutate(maxRho = max(rho), optimal = ifelse(rho == maxRho, T, F)) %>%
  subset(optimal) %>%
  as.data.frame() %>%
  group_by(species, litter, Var) %>%
  summarise(candE = paste(unique(E), collapse = ', '),
            rho = mean(rho))

# Identify optimal E 
# that creates the highest mean rho on all traps of same assemblage
uvsmap.optE = lapply(unique(uvsmap.optEbyTrap$species), function(sp){
  lapply(c('FDM', 'LDM'), function(lt){
    baseInfo = uvsmap.optEbyTrap %>% subset(species == sp & litter == lt)
    trapList = ltc.match[ltc.match$species == sp,]$trapID
    varToLook = baseInfo$Var
    # Search each variable
    byV = lapply(varToLook, function(targetV){
      paramRow = baseInfo %>% subset(Var == targetV)
      eToLook = as.numeric(strsplit(paramRow$candE, 
                                    split = ', ', fixed = T)[[1]])
      if (length(eToLook) >1){
        byE = lapply(eToLook, function(testE){
          # Perform Simplex under same E on each trap
          meanRho = lapply(trapList, function(tr) {
            inputBase = uvsmap.baseTbl %>% 
              subset(trap == tr & litter == lt) %>%
              select(time, DM, Env = all_of(targetV))
            if (testE > 1){
              embBase = inputBase %>% select(time, Env)
              inputEmb = Embed(dataFrame = embBase,
                               E = testE,
                               columns = 'Env')
              colInput = c('Env', paste('Env', 1:(testE-1), sep = ''))
              names(inputEmb) = colInput
              inputTbl = inputBase %>% select(time, DM) %>%
                cbind(inputEmb) %>% slice(-(1:(testE-1)))
            } else {
              inputTbl = inputBase; colInput = 'Env'
            }
            # Simplex
            trainLib = '1 179'; testLib = paste(180, 
                                                dim(inputTbl)[1], 
                                                collapse = ' ')
            mod = Simplex(dataFrame = inputTbl,
                          lib = trainLib,
                          pred = testLib,
                          E = testE,
                          columns = paste(colInput, collapse = ' '),
                          target = 'DM',
                          embedded = T,
                          showPlot = F)
            rho = ComputeError(mod$Observations, mod$Predictions)$rho
            return(rho)
          }) %>% unlist() %>% mean()
          resTbl = data.frame(species = sp, litter = lt, Var = targetV,
                              optE = testE, rho = meanRho)
          return(resTbl)
        }) %>% do.call(rbind.data.frame, args = .)
      } else {
        byE = data.frame(species = sp, litter = lt, Var = targetV,
                         optE = paramRow$candE, rho = paramRow$rho)
      }
      return(byE)
    }) %>% do.call(rbind.data.frame, args = .) #byV
  }) %>% do.call(rbind.data.frame, args = .)
}) %>% do.call(rbind.data.frame, args = .) %>%
  group_by(species, litter, Var) %>%
  mutate(maxRho = max(rho), optimal = ifelse(rho == maxRho, T, F)) %>%
  ungroup() %>% subset(optimal) %>% 
  select(species, litter, Var, optE, rho) %>% as.data.frame()

# ---- Optimal Theta ----

# Best theta for each trap
uvsmap.paramList = pbapply::pblapply(levels(ltc.match$species), function(sp){
  trapList = ltc.match[ltc.match$species == sp,]$trapID
  byLT = lapply(c('FDM', 'LDM'), function(lt){
    varList = edm.Var.base[[paste(sp, lt)]]
    byV = lapply(varList, function(targetV){
      paramRow = uvsmap.optE %>% subset(species == sp & 
                                          litter == lt & 
                                          Var == targetV)
      optE = as.numeric(paramRow$optE)
      # Process by trap
      byTrap = lapply(trapList, function(tr){
        baseTbl = uvsmap.baseTbl %>% subset(trap == tr & litter == lt) %>%
          select(time, DM, Env = all_of(targetV))
        if (optE > 1){
          embBase = baseTbl %>% select(time, Env)
          inputEmb = Embed(dataFrame = embBase,
                           E = optE,
                           columns = 'Env')
          colInput = c('Env', paste('Env', 1:(optE-1), sep = ''))
          names(inputEmb) = colInput
          inputTbl = baseTbl %>% select(time, DM) %>%
            cbind(inputEmb) %>% slice(-(1:(optE-1)))
        } else {
          inputTbl = baseTbl; colInput = 'Env'
        }
        # Identify optimal theta
        thetaTbl = PredictNonlinear(dataFrame = inputTbl,
                                    lib = '1 179',
                                    pred = paste(180, dim(inputTbl)[1], sep = ' '),
                                    E = optE,
                                    columns = paste(colInput, collapse = ' '),
                                    target = 'DM',
                                    embedded = T,
                                    showPlot = F)
        optTheta = thetaTbl[which.max(thetaTbl$rho),]$Theta
        optRho = thetaTbl[which.max(thetaTbl$rho),]$rho
        # Arrange to output
        outputList = list('data' = inputTbl,
                          'col' = paste(colInput, collapse = ' '),
                          'E' = optE,
                          'theta' = optTheta,
                          'rho' = optRho)
        return(outputList)
      }); names(byTrap) = trapList; return(byTrap)
    }); names(byV) = varList; return(byV)
  }); names(byLT) = c('FDM', 'LDM'); return(byLT)
}); names(uvsmap.paramList) = levels(ltc.match$species)

# Extract parameters from list to table
uvsmap.paramTbl = lapply(levels(ltc.match$species), function(sp){
  lapply(c('FDM', 'LDM'), function(lt){
    # Retrieve optimal parameters for all traps of the species
    listToCheck = uvsmap.paramList[[sp]][[lt]]
    varToCheck = names(listToCheck)
    tableConverted = lapply(varToCheck, function(V){
      lapply(1:length(listToCheck[[V]]), function(tID){
        baseInfo = listToCheck[[V]][[tID]]
        outTbl = data.frame(species = sp,
                            litter = lt,
                            var = V,
                            trap = as.numeric(names(listToCheck[[V]])[tID]),
                            colInput = baseInfo[['col']],
                            E = baseInfo[['E']],
                            theta = baseInfo[['theta']],
                            rho = baseInfo[['rho']])
        return(outTbl)
      }) %>% do.call(rbind.data.frame, args = .)
    }) %>% do.call(rbind.data.frame, args = .)
    return(tableConverted)
  }) %>% do.call(rbind.data.frame, args = .)
}) %>% do.call(rbind.data.frame, args = .) %>%
  group_by(species, litter, var, colInput, E) %>%
  summarise(avg = mean(rho),
            theta = paste(unique(theta), collapse = ', ')) %>%
  ungroup() %>% as.data.frame()

# Identify optimal Theta 
# that creates the highest mean rho on all traps of same assemblage
uvsmap.perform = lapply(levels(ltc.match$species), function(sp){
  byLT = lapply(c('FDM', 'LDM'), function(lt){
    trapList = ltc.match[ltc.match$species == sp,]$trapID
    varList = edm.Var.base[[paste(sp, lt)]]
    byVar = lapply(varList, function(V){
      paramRow = uvsmap.paramTbl %>% subset(species == sp & litter == lt & var == V)
      optE = paramRow$E; colInput = paramRow$colInput
      thetaToLook = as.numeric(strsplit(paramRow$theta, split = ', ', fixed = T)[[1]])
      #### Perform by theta
      byTheta = lapply(thetaToLook, function(optT){
        #### Perform S-Map by trap
        byTrap = lapply(trapList, function(tr){
          inputTbl = uvsmap.paramList[[sp]][[lt]][[V]][[as.character(tr)]][['data']]
          mod = SMap(dataFrame = inputTbl,
                     lib = '1 179',
                     pred = paste(180, dim(inputTbl)[1], collapse = ' '),
                     E = optE,
                     theta = optT,
                     columns = colInput,
                     target = 'DM',
                     embedded = T)
          smapTbl = data.frame(trap = tr, mod$predictions)
          return(smapTbl)
        }); names(byTrap) = trapList
        #### Summary rho under given theta
        thetaSummary = lapply(byTrap, function(eachTrap){
          outTbl = data.frame(theta = optT,
                              rho = ComputeError(eachTrap$Observations,
                                                 eachTrap$Predictions)$rho)
          return(outTbl)
        }) %>% do.call(rbind.data.frame, args = .) %>%
          group_by(theta) %>% summarise(rho = mean(rho))
        #### Save to output
        outList = list('base' = byTrap, 
                       'rho' = thetaSummary)
      }); names(byTheta) = thetaToLook
      ### Evaluate optimal theta
      evalTheta = lapply(1:length(byTheta), function(ID){
        byTheta[[ID]][['rho']]
      }) %>% do.call(rbind.data.frame, args = .) %>%
        mutate(maxRho = max(rho), optimal = ifelse(rho == maxRho, T, F)) %>%
        subset(optimal)
      ### Optimal theta & Smap prediction table
      outList = list('theta' = evalTheta$theta,
                     'rho' = evalTheta$rho,
                     'base' = byTheta[[as.character(evalTheta$theta)]][['base']])
      return(outList)
    }); names(byVar) = varList; return(byVar)
  }); names(byLT) = c('FDM', 'LDM'); return(byLT)
}); names(uvsmap.perform) = levels(ltc.match$species)

# ---- Scenario Exploration ----
# Algorithm modified from Nova et al 2020 (DOI: 10.1111/ele.13652)

uvsmap.scneExp = pbapply::pblapply(levels(ltc.match$species), function(sp){
  byLT = lapply(c('FDM', 'LDM'), function(lt){
    trapList = ltc.match[ltc.match$species == sp,]$trapID
    # Parameters of the system
    baseParam = uvsmap.paramTbl %>% subset(species == sp & litter == lt)
    varToLook = baseParam$var
    # Perform by variable
    byVar = lapply(varToLook, function(V){
      optE = subset(baseParam, var == V)$E
      optTheta = uvsmap.perform[[sp]][[lt]][[V]][['theta']]
      colInput = subset(baseParam, var == V)$colInput
      # By trap
      byTrap = lapply(trapList, function(tr){
        baseTbl = uvsmap.baseTbl %>% subset(trap == tr & litter == lt) %>%
          select(time, DM, Env = all_of(V))
        # Modify environmental factor
        valueToChange = baseTbl$Env
        dX = 0.1 * sd(valueToChange)
        xPlus = xMinus = valueToChange; maxLength = length(valueToChange)
        xPlus[(180+(optE-1)):maxLength] = xPlus[(180+(optE-1)):maxLength] + dX
        xMinus[(180+(optE-1)):maxLength] = xMinus[(180+(optE-1)):maxLength] - dX
        plusTbl = data.frame(time = baseTbl$time, Env = xPlus)
        minusTbl = data.frame(time = baseTbl$time, Env = xMinus)
        # Embed table
        plusEmb = Embed(dataFrame = plusTbl, E = optE, columns = 'Env')
        minusEmb = Embed(dataFrame = minusTbl, E = optE, columns = 'Env')
        names(plusEmb) = names(minusEmb) = strsplit(colInput, split = ' ', fixed = T)[[1]]
        plusTbl = bind_cols(baseTbl %>% select(time, DM), plusEmb) %>%
          slice(-(1:(optE-1)))
        minusTbl = bind_cols(baseTbl %>% select(time, DM), minusEmb) %>%
          slice(-(1:(optE-1)))
        rm(plusEmb, minusEmb)
        # S-Map
        dxPlus = SMap(dataFrame = plusTbl,
                      lib = '1 179',
                      pred = paste(180, dim(plusTbl)[1], collapse = ' '),
                      E = optE,
                      theta = optTheta,
                      columns = colInput,
                      target = 'DM',
                      embedded = T)
        dxMinus = SMap(dataFrame = minusTbl,
                       lib = '1 179',
                       pred = paste(180, dim(minusTbl)[1], collapse = ' '),
                       E = optE,
                       theta = optTheta,
                       columns = colInput,
                       target = 'DM',
                       embedded = T)
        scneOut = bind_rows(data.frame(mode = 'plus',
                                       dxPlus$predictions),
                            data.frame(mode = 'minus',
                                       dxMinus$predictions)) %>%
          mutate(dX = dX) %>%
          select(time, dX, Y1 = Predictions, mode)
        # Join SMap predictions with unchanged environment
        smapBase = uvsmap.perform[[sp]][[lt]][[V]][['base']][[as.character(tr)]]
        scneRes = smapBase %>% select(time, Y0 = Predictions) %>% 
          left_join(scneOut, by = 'time') %>%
          subset(!is.na(Y0)) %>%
          mutate(species = sp, litter = lt, trap = tr, var = V,
                 dY = Y1 - Y0,
                 sens = ifelse(mode == 'plus', dY/dX, -dY/dX)) %>%
          select(species, litter, trap, time, var, sens)
        return(scneRes)
      }) %>% do.call(rbind.data.frame, args = .)
    }) %>% do.call(rbind.data.frame, args = .) 
  }) %>% do.call(rbind.data.frame, args = .) 
}) %>% do.call(rbind.data.frame, args = .)

# ---- Visualization ----

# For visualization of S-Map results
# please refer to 'M1_SMapPlotting.R'