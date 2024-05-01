## MR analysis for two-sample datasets 
##### Load Library ########
set.seed(100)
library(MRInstruments) 
library(readr)
library(MendelianRandomization)
library(ggplot2)
library(plotly)
library(TwoSampleMR)
library(MRPRESSO)
source("mr_simex.R")
source("IVO.R")

### EAST ASIA Cholesterol CAD #######
et <- extract_instruments("bbj-a-54") # 44  15
ot <- extract_outcome_data(snps=et$SNP, outcomes="bbj-a-159") # 43  23
e_LD <- clump_data(et, pop = "EAS") # 44 15
o_LD <- clump_data(ot, pop = "EAS") # 43 24
ht <- harmonise_data(e_LD, o_LD[1:13]) # 43
mr_default <- mr(ht) # Default Parameter
sim_mr <- mr_simex(harmo_dat = ht)


# Get instruments --> Total Cholesterol 
exposure_dat <- extract_instruments("bbj-a-54", p1 = 5e-03, r2 = 0.01) # 2474 16

# Get effects of instruments on outcome --> CAD
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="bbj-a-159") # 2343 23

# IVO apply
exp_t <- IVO(exposure_dat) # 349  16
out_t <- IVO(outcome_dat) # 457  16

expo_LD <- clump_data(exp_t, pop = "EAS") # 200 16
outc_LD <- clump_data(out_t , pop = "EAS") # 252 16

# Harmonise the exposure and outcome data
dat_t <- harmonise_data(expo_LD , outc_LD[1:13] ) # 41 33

mr_res_t <- mr(dat_t)
sim_mr <- mr_simex(harmo_dat = dat_t)
# Apply all MR methods

# Directionality
dir_mr <- directionality_test(dat_t)
#  Heterogeneity Statistics
het_mr <- mr_heterogeneity(dat_t)
# horizontal pleiotropy in MR analysis
ple_mr <- mr_pleiotropy_test(dat_t)

##### MR presso- for horizontal pleiotropy 
mr_p <- mr_presso(
  BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
  SdOutcome = "se.outcome", SdExposure = "se.exposure", 
  OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
  data = dat_t, NbDistribution = 3000, SignifThreshold = 0.05)


# Scatter Plot
mr_scatter_plot(mr_res_t, dat_t)

#Leave-one-out analysis (mr_ivw)
mr_loo_tc <- mr_leaveoneout(dat_t)
mr_leaveoneout_plot(mr_loo_tc)

# single SNP analyses (Forest plot)
res1_single <- mr_singlesnp(dat_t, all_method=c("mr_ivw", "mr_egger_regression"))
mr_forest_plot(res1_single)

########################## Liver ########################
# liver iron and liver cancer
et <- extract_instruments("ebi-a-GCST90016674") # 10  15
ot <- extract_outcome_data(snps=et$SNP, outcomes="ieu-b-4953") # 8  23
e_LD <- clump_data(et, pop = "EUR") # 7 15
o_LD <- clump_data(ot, pop = "EUR") # 5 24
ht <- harmonise_data(e_LD, o_LD[1:13]) # 5
mr_default <- mr(ht)
sim_mr <- mr_simex(harmo_dat = ht)
all_LIC_LCC <- mr_wrapper(dat_t)
sum(all_LIC_LCC$`ebi-a-GCST90016674.ieu-b-4953`$snps_retained$steiger==T) # 14
# Proposed Pipeline

exposure_dat <- extract_instruments("ebi-a-GCST90016674", p1 = 5e-03, r2 = 0.01) # 4233 15

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-b-4953") # 1624 23

# IVO apply
exp_t <- IVO(exposure_dat) # 1070   16
out_t <- IVO(outcome_dat) #  355  16

expo_LD <- clump_data(exp_t, pop = "EUR") # 473  16
outc_LD <- clump_data(out_t , pop = "EUR") # 221 25

expo_LD$F_test <- expo_LD$t_value**2
outc_LD$F_test <- outc_LD$t_value**2

# Harmonise the exposure and outcome data
dat_t <- harmonise_data(expo_LD , outc_LD[1:13] )  # 14  33

#  Heterogeneity Statistics
het_mr <- mr_heterogeneity(dat_t)

# horizontal pleiotropy in MR analysis
ple_mr <- mr_pleiotropy_test(dat_t)

####### MR presso- for horizontal pleiotropy 
mr_p <- mr_presso(
  BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
  SdOutcome = "se.outcome", SdExposure = "se.exposure", 
  OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
  data = dat_t, NbDistribution = 3000, SignifThreshold = 0.05)

# MR test
mr_res_t <- mr(dat_t)
sim_mr <- mr_simex(harmo_dat = dat_t)

# Directionality
out <- directionality_test(dat_t)

# Scatter Plot
mr_scatter_plot(mr_res_t, dat_t)

#Leave-one-out analysis (mr_ivw)
mr_loo_tc <- mr_leaveoneout(dat_t)
mr_leaveoneout_plot(mr_loo_tc)

# single SNP analyses (Forest plot)
res1_single <- mr_singlesnp(dat_t, all_method=c("mr_ivw", "mr_egger_regression"))
mr_forest_plot(res1_single)

########################## Obesity & osteoarthritis ########################
# Childhood obesity and Knee and hip osteoarthritis
# Default pipeline
et <- extract_instruments("ieu-a-1096") # 5  15
ot <- extract_outcome_data(snps=et$SNP, outcomes="ieu-a-1170") # 5  23
e_LD <- clump_data(et, pop = "EUR") # 7 15
o_LD <- clump_data(ot, pop = "EUR") # 5 24
ht <- harmonise_data(e_LD, o_LD[1:13]) # 5
mr_default <- mr(ht)
sim_mr <- mr_simex(harmo_dat = ht)

# Proposed pipeline

exposure_dat <- extract_instruments("ieu-a-1096", p1 = 5e-03, r2 = 0.01) # 1707 15

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-1170") # 1267 23

# IVO apply
exp_t <- IVO(exposure_dat) # 376  16
out_t <- IVO(outcome_dat)  # 291  16


expo_LD <- clump_data(exp_t, pop = "EUR") 
outc_LD <- clump_data(out_t , pop = "EUR") 
# Harmonise the exposure and outcome data
dat_t <- harmonise_data(expo_LD , outc_LD[1:13] ) # 27  33

mr_res_t <- mr(dat_t)
sim_mr <- mr_simex(harmo_dat = dat_t)
gc()
#  Heterogeneity Statistics
het_mr <- mr_heterogeneity(dat_t)

# horizontal pleiotropy in MR analysis
ple_mr <- mr_pleiotropy_test(dat_t)

##### MR presso- for horizontal pleiotropy 
mr_p <- mr_presso(
  BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
  SdOutcome = "se.outcome", SdExposure = "se.exposure", 
  OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
  data = dat_t, NbDistribution = 3000, SignifThreshold = 0.05)

# Directionality
dir_mr <- directionality_test(dat_t)

# Scatter Plot
mr_scatter_plot(mr_res_t, dat_t)

#Leave-one-out analysis (mr_ivw)
mr_loo_tc <- mr_leaveoneout(dat_t)
mr_leaveoneout_plot(mr_loo_tc)

# single SNP analyses (Forest plot)
res1_single <- mr_singlesnp(dat_t, all_method=c("mr_ivw", "mr_egger_regression"))
mr_forest_plot(res1_single)


############ TC-CAD (European) #############
# Get instruments --> Total Cholesterol (European)
exposure_dat <- extract_instruments("met-c-933", p1 = 5e-03, r2 = 0.01) # 4186 16

# Get effects of instruments on outcome --> CAD (European)
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST000998") # 1344 24

# IVO apply
exp_t <- IVO(exposure_dat) # 787  16
out_t <- IVO(outcome_dat) #  258  24

expo_LD <- clump_data(exp_t, pop = "EUR") # 393 16
outc_LD <- clump_data(out_t , pop = "EUR") # 179 25
# Harmonise the exposure and outcome data
dat_t <- harmonise_data(expo_LD , outc_LD[1:13] )  # 12

# Perform MR
mr_res_t <- mr(dat_t)
# MR Egger SIMEX
sim_mr <- mr_simex(harmo_dat = dat_t)

# Directional test to check that exposure causes outcome is valid
direction_t <- directionality_test(dat_t)

#  Heterogeneity Statistics
het_mr <- mr_heterogeneity(dat_t)

# horizontal pleiotropy in MR analysis
ple_mr <- mr_pleiotropy_test(dat_t)

##### MR presso- for horizontal pleiotropy 
mr_p <- mr_presso(
  BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
  SdOutcome = "se.outcome", SdExposure = "se.exposure", 
  OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
  data = dat_t, NbDistribution = 3000, SignifThreshold = 0.05)

# Directionality
dir_mr <- directionality_test(dat_t)

# Scatter Plot
mr_scatter_plot(mr_res_t, dat_t)

#Leave-one-out analysis (mr_ivw)
mr_loo_tc <- mr_leaveoneout(dat_t)
mr_leaveoneout_plot(mr_loo_tc)

# single SNP analyses (Forest plot)
res1_single <- mr_singlesnp(dat_t, all_method=c("mr_ivw", "mr_egger_regression"))
mr_forest_plot(res1_single)

## With Default parameters
et <- extract_instruments("met-c-933") # 23  15
ot <- extract_outcome_data(snps=et$SNP, outcomes="ebi-a-GCST000998") # 11  23
#cad <- extract_outcome_data(snps=et$SNP, "ebi-a-GCST005195") # 20   23
e_LD <- clump_data(et, pop = "EUR") # 23 15
o_LD <- clump_data(ot, pop = "EUR") # 11 24
ht <- harmonise_data(e_LD, o_LD[1:13]) # 11
mr_default <- mr(ht)
sim_mr <- mr_simex(harmo_dat = ht)

##################################
gc()
