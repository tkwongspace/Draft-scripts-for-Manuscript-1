# M1 Stationary Tests

# This script proceed stationary test for litter biomass differentials
# Two approaches are involved:
# 1. The Kwiatkowski et al. Unit Root Test (ur.kpss)
# 2. Augmented Dickey-Fuller Test (adf.test)
# 
# This script follows the idea of 
# https://stats.stackexchange.com/questions/239360/contradictory-results-of-adf-and-kpss-unit-root-tests/239367
# where @Ferdi suggested A UNIT ROOT TEST SHOULD BE APPLIED DEPENDING ON KPSS TEST
# 
# There are four cases of results that may be output:
# i. ADF-test non-stationary + KPSS non-staionary = Unit root exists
# ii. ADF-test stationary + KPSS stationary = Stationary time series
# iii. ADF-test non-stationary + KPSS stationary = not enough observations
# iv. ADF-test stationary + KPSS non-stationary = potential heteroskedasticity
#
# For case [iv] a more profound approach shall be applied, i.e., Variance Ratio Test
# that measures if the data is 'between stationary and a unit root'.

library(tseries)
library(urca)

# Function to evaluate test results
testResSummary = function(sigTau, sigMu, sigADF){
  if (!(sigTau & sigMu & sigADF)){
    # Case i, all bad
    res = 'Non-stationary'
  } else if (sigTau & sigMu & sigADF){
    # Case ii, all good
    res = 'Stationary'
  } else if (!sigTau & sigMu & sigADF) {
    # Case iv
    res = 'Trends exist'
  } else if (sigTau & !sigMu & sigADF){
    # Case iv
    res = 'Shifts exist'
  } else if (!sigTau & !sigMu & sigADF){
    # Case iv
    res = 'Heteroskedasticity'
  } else if (sigTau & sigMu & !sigADF) {
    res = 'Obs. Not enough'
  }
  return(res)
}

# Litter biomass differentials
tpDif = lapply(unique(tp$trapID), function(tr){
  lapply(c('FDM', 'LDM'), function(lt){
    oldTbl = tp %>% subset(trapID == tr & litter == lt)
    newDM = c(NA, diff(oldTbl$DM))
    newDF = data.frame(date = oldTbl$date,
                       species = oldTbl$species,
                       litter = lt,
                       trapID = tr,
                       DM = newDM)
    return(newDF)
  }) %>% do.call(rbind.data.frame, args = .)
}) %>% do.call(rbind.data.frame, args = .) %>% na.omit()

# Apply tests
testDifParamTbl = tpDif %>% nest_by(trapID, litter) %>%
  mutate(kptau = list(summary(ur.kpss(data$DM, type = 'tau'))),
         kpmu = list(summary(ur.kpss(data$DM, type = 'mu'))),
         adf = list(adf.test(data$DM))) %>%
  summarise(estTau = kptau@teststat,
            criTau = kptau@cval[2],
            estMu = kpmu@teststat,
            criMu = kpmu@cval[2],
            adfp = adf$p.value) %>%
  mutate(sigTau = ifelse(estTau < criTau, T, F),
         sigMu = ifelse(estMu < criMu, T, F),
         sigADF = ifelse(adfp < 0.05, T, F),
         result = testResSummary(sigTau, sigMu, sigADF)) %>%
  select(-starts_with('sig'))
