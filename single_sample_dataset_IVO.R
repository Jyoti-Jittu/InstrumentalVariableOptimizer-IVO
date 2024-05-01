## Single-sample MR analysis

#### Load library ####
set.seed(100)
library(TwoSampleMR)
library(MRPRESSO)
library(snpStats)
library(ggplot2)
source("LinearR.R")
source("logisticR.R")
source("mr_simex.R")
source("IVO.R")

for(pkg in c("snpStats", "doParallel", "SNPRelate","GenABEL" )){
  if(!require(pkg, character.only = T)) {
    stop("At least one pckg is required for this script. Please install it on your system.")
  }
}

############# Read Raw  Genotype Data ########
path <- paste("processed_BED_ld_.3",
              c(".bed", ".bim", ".fam"), sep = "")
SNP_Detail <- read.plink(path[1], path[2], path[3])

dim(SNP_Detail$genotypes) # 1277 41,802
dim(SNP_Detail$map)      # 41802     6
dim(SNP_Detail$fam)     # 1277    6
gc()


############ Read Phenotype Data ################
Phenotype = read.csv("Total_Chl_ld_.3",  row.names = 1) 
dim(Phenotype)  # 1277    2
colnames(Phenotype)<-c("id", "phenotype")
gc()
phenodata <- data.frame("id" = Phenotype$id,"phenotype" = Phenotype$phenotype, stringsAsFactors = F)

############ Linear Regression (Total Chol) ############
target <- "chl_ld_.3"

start <- Sys.time()
GWAA(genodata = SNP_Detail$genotypes, phenodata = phenodata,
     filename = paste(target, ".txt",sep = ""))
Sys.time() - start


########### Logistic Regression (CAD) ###########
target <- "CAD_ld_.3"

trait <- read.csv("GWAStutorial_clinical.csv")
trait <- trait[trait$FamID %in% Phenotype$id,] # 1277 7

phenodata <- data.frame("id" = trait$FamID,
                        "phenotype" = trait$CAD, stringsAsFactors = F)

start <- Sys.time()
GWAA_logi(genodata = SNP_Detail$genotypes, phenodata = phenodata,
          filename = paste(target, ".txt",sep = ""))
Sys.time() - start
########## Summary-Stat-Data-Preparation #########
exposure_data <- read.table("chl_ld_.3.txt", header = TRUE)
outcome_data <- read.table("CAD_ld_.3.txt", header = TRUE)

sumStat <- col.summary(SNP_Detail$genotypes)
#write.csv(sumStat, "gwasColSum")
sumStat <-  read.csv("gwasColSum")

exp_dat <- cbind.data.frame(exposure_data,SNP_Detail$map$allele.1, SNP_Detail$map$allele.2, 
                            sumStat$MAF, SNP_Detail$map$chromosome, SNP_Detail$map$position)
out_dat <- cbind.data.frame(outcome_data,SNP_Detail$map$allele.1, SNP_Detail$map$allele.2, 
                            sumStat$MAF, SNP_Detail$map$chromosome, SNP_Detail$map$position)

exp_dat$samplesize <- 1277 # Effect Size
out_dat$samplesize <- 1277 # Effect Size

exp_dat$id <- "chl"
out_dat$id <- "CAD"

colnames(exp_dat)<- c("SNP","beta.exposure","se.exposure","t.value","p.value",
                      "effect_allele.exposure","other_allele.exposure",
                      "eaf.exposure","chromosome","position","samplesize",
                      "id")
colnames(out_dat)<- c("SNP","beta.outcome","se.outcome","t.value","p.value",
                      "effect_allele.outcome","other_allele.outcome",
                      "eaf.outcome","position","chromosome","samplesize",
                      "id")
# 41802 SNPs

write.csv(exp_dat, "chl_SummaryStatistics")
write.csv(out_dat, "CAD_SummaryStatistics")




###### Read exposure data ########
expo <- format_data(
  exp_dat,
  type = "exposure",
  header = TRUE,
  phenotype_col = "chl",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "p.value",
  id_col = "id",
  min_pval = 1e-200,
  chr_col = "chromosome",
  samplesize_col = "samplesize",
  pos_col = "chromosome",
  log_pval = FALSE
)
gc()

########### Read outcome data ##############
outc <- format_data(
  out_dat,
  type = "outcome",
  header = TRUE,
  phenotype_col = "CAD",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  eaf_col = "eaf.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "p.value",
  id_col = "id",
  min_pval = 1e-200,
  chr_col = "chromosome",
  samplesize_col = "samplesize",
  pos_col = "chromosome",
  log_pval = FALSE
)

gc()
write.csv(expo, "expo_before_LD")
write.csv(outc, "outcome_before_LD")

######## LD clumping on Both Data ##########
expo_LD <- clump_data(expo, clump_kb = 10000, clump_r2 = 0.01, pop = "AMR") # 4080   13

outc_LD <- clump_data(outc, clump_kb = 10000, clump_r2 = 0.01, pop = "AMR") # 4081   14


####### Default parameter for MR comparison ####
default_har <- harmonise_data(expo_LD, outc_LD[1:13]) # 679 30
dafault_mr <- mr(default_har)

######## Apply InstrumentalVariableOptimizer ####
exp_t <- IVO(expo_LD) # 988  14
out_t <- IVO(outc_LD) #  959  15


# Harmonised data
dat_t <- harmonise_data(exp_t , out_t[1:13]) # 26 31


## MR Egger SIMEX ######
sim_mr <- mr_simex(harmo_dat = dat_t)

##### Mendelian randomisation ######
mr_res_t <- mr(dat_t)
# Apply all MR methods
all_ss <- mr_wrapper(dat_t)
all_ss$chl.CAD$snps_retained$SNP[ss_moe$chl.CAD$snps_retained$steiger==T] # 74
########### directionality ######
directionality_test(dat_t) # TRUE

# Scatter Plot
mr_scatter_plot(mr_res_t, dat_t)

#Leave-one-out analysis (mr_ivw)
mr_loo_tc <- mr_leaveoneout(dat_t)
mr_leaveoneout_plot(mr_loo_tc)

# single SNP analyses (Forest plot)
res1_single <- mr_singlesnp(dat_t, all_method=c("mr_ivw", "mr_egger_regression"))
mr_forest_plot(res1_single)

gc()
