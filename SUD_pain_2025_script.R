
# Title: Genetic interrelationships between chronic pain and substance use disorders and functional enrichment
# Author: Lydia Rader
# Purpose: Munging summary statistics, cross-trait LDSC, Genomic SEM, and stratified genomic SEM for diss


###### Get set up##############
setwd("/Users/lydiarader/Documents/pain_research/GWAS_summary_statistics/dissertation/data")

require(data.table)
require(R.utils)
require(readr)
require(dplyr)
require(GenomicSEM)
library(ggplot2)
library(reshape2)

###########PAU########

# Action steps:
# 1. Make sure that our columns all look good
# 2. Rename Sample size to Neff

#read in the alcohol use 2023 file
pau<-fread("/Users/lydiarader/Documents/pain_research/GWAS_summary_statistics/dissertation/data/mvp_PAU_OUD_2023/PAU_2023/PAU_EUR_Aug2023.txt.gz",data.table=FALSE)
head(pau) #check that our columns are as expected; 
pau$EA <- NULL
# We have the EAF and Andrew says to use this one
#output datafile
write.table(pau, file = "PAU_withrsID.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE) # This is the munging file for PAU
rm(pau) 



###########CUD ########
#Original file name: dbGAP_EUR_CanMVP_META
#Levey 2023
#File to munge:
#Sumstats: CUD.sumstats

cud<-fread("~/Documents/pain_research/GWAS_summary_statistics/dissertation/data/CUD_Levey_2023/dbGAP_EUR_CanMVP_META.gz",data.table=FALSE)
head(cud)
# Preparation of CUD
#Column names = SNP  MarkerName CHR     BP A1  A2    BETA     SE P.value


cud$MarkerName<-NULL # drop this because it is duplicate of SNP
cud$Neff <- 161053
#cud$SNP <- sub("12:", "", cud$SNP) # remove the 12: since we have this info in CHR file

# Check the result
head(cud)

#output datafile
write.table(cud, file = "CUD_withrsID.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)
rm(cud)

###########OUD ########
#Original file name: OUD_EUR_Sumstats_Deak2022.gz
#File info: Deak 2022
#ReadME notes:
#For sample size of each of the respective GWAS included in this manuscript, we recommend using 
#the effective sample size: AFR OUD effective N = 20,032; EUR OUD effective N = 56,994; OUD-MTAG effective N = 128,748.
#For analyses requiring alternative sample size calculations, the N of all included cohorts can be found in Table 1.
# For allele frequency we recommend using ancestry-specific 1000 Genomes Project reference panels if required for your analysis.

#read in the OUD file 
oud<-fread("~/Documents/pain_research/GWAS_summary_statistics/dissertation/data/Deak_2022_OUD_sumstats/OUD_EUR_Sumstats_Deak2022.gz",data.table=FALSE)
head(oud)
#Column names =  snp chr      bp a1 a2      z p.value
# Action steps:
# 1. Add in effective sample size: N = 56,994
# 2. Add in MAF

# Keep naming consistent
head(oud)
colnames(oud)[1] <-c("SNP")
head(oud)

#output datafile
write.table(oud, file = "OUD_withrsID.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)
rm(oud)

###########CPD ########
cpd<-fread("~/Documents/pain_research/GWAS_summary_statistics/dissertation/data/gscan_DPW_CPD/EUR_stratified/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt.gz",data.table=FALSE)
head(cpd)
colnames(cpd)[3] <-c("SNP")
head(cpd)

write.table(cpd, file = "CPD_withrsID.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)

###########DPW########
dpw<-fread("~/Documents/pain_research/GWAS_summary_statistics/dissertation/data/gscan_DPW_CPD/EUR_stratified/GSCAN_DrnkWk_2022_GWAS_SUMMARY_STATS_EUR.txt.gz",data.table=FALSE)
head(dpw)
colnames(dpw)[3] <-c("SNP")
head(dpw)

write.table(dpw, file = "DPW_withrsID.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)

########AUDIT-C######
# Summary statistics of the Genome-Wide Association Study (GWAS) of AUDIT total score,
# AUDIT-C and AUDIT-P limited to UKBiobank (UKB) research participants of European
# ancestry. We are only needing the AUDIT-C for consumption so drop the total score and the problem score. AUDIT-P is already included in PAU.
audc <-fread("~/Documents/pain_research/GWAS_summary_statistics/dissertation/data/AUDIT_PGC/AUDIT_UKB_2018_AJP.txt.gz",data.table=FALSE)
head(audc)
colnames(audc)[2] <-c("SNP")
head(audc)

audc<-audc %>% select(chr, rsid, a_0, a_1, info, beta_C, se_C, p_C, N) %>%
  rename(SNP = rsid,
         A1 = a_1,
         A2 = a_0,
         beta = beta_C,
         SE = se_C,
         P = p_C
         )

head(audc)

write.table(audc, file = "AUDC_withrsID.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)

#### Multisite Pain ####


mcpain <-fread("~/Documents/pain_research/GWAS_summary_statistics/dissertation/data/chronic_pain-bgen.stats",data.table=FALSE)
head(mcpain)

mcpain<-mcpain %>% select(SNP, CHR, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM_INF) %>%
  rename(A1 = ALLELE1,
         A2 = ALLELE0,
         P = P_BOLT_LMM_INF
  )


####NEED TO ADD P VALUE

write.table(mcpain, file = "MCPAIN_withrsID.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)


setwd("/Users/lydiarader/Documents/pain_research/GWAS_summary_statistics/dissertation/data")
ld<-"~/Documents/pain_research/GWAS_summary_statistics/dissertation/data/eur_w_ld_chr/"
wld<-"~/Documents/pain_research/GWAS_summary_statistics/dissertation/data/eur_w_ld_chr/"
files<-c("MCPAIN_withrsID.txt")
trait.names<-c("MCPAIN") # We are not munging TUD because AH provided munged sumstats already
N=(387649) 
info.filter=0.9
maf.filter=0.01
hm3<-"~/Documents/pain_research/GWAS_summary_statistics/dissertation/data/eur_w_ld_chr/w_hm3.snplist"
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)


########## Do the munge ############
#the folder of LD scores
setwd("/Users/lydiarader/Documents/pain_research/GWAS_summary_statistics/dissertation/data")

ld<-"~/Documents/pain_research/GWAS_summary_statistics/dissertation/data/eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld<-"~/Documents/pain_research/GWAS_summary_statistics/dissertation/data/eur_w_ld_chr/"

# Traits
files<-c("PAU_withrsID.txt", "CUD_withrsID.txt", "OUD_withrsID.txt","CPD_withrsID.txt", "DPW_withrsID.txt", "AUDC_withrsID.txt")

#name the traits
trait.names<-c("PAU", "CUD", "OUD", "CPD", "DPW", "AUDc") # We are not munging TUD because AH provided munged sumstats already

#sample size for traits in file already; providing the effective sample size for CUD and OUD
N=c(NA, 161053, 56994, NA, NA, NA) 

#define the imputation quality filter
info.filter=0.9

#define the MAF filter
maf.filter=0.01


#define the reference file being used to allign alleles across summary stats
#here we are using hapmap3
hm3<-"~/Documents/pain_research/GWAS_summary_statistics/dissertation/data/eur_w_ld_chr/w_hm3.snplist"

# MUNGE 
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)
#manually move sumstats into sumstats folder



################LDSC #######################
setwd("/Users/lydiarader/Documents/pain_research/GWAS_summary_statistics/dissertation/data/sumstats")
#vector of munged summary stats
# 4 SUDs
# 3 SU consumption traits
# 24 CP conditions
# 1 MS pain (this is in output2)
traits<-c("OUD.sumstats.gz","PAU.sumstats.gz","TUD.sumstats.gz","CUD.sumstats.gz",
          "AUDc.sumstats.gz", "DPW.sumstats.gz", "CPD.sumstats.gz",
          "hdch.sumstats.gz", "mgrn.sumstats.gz", "neck.sumstats.gz",
          "back.sumstats.gz", "chDs.sumstats.gz","chPh.sumstats.gz",
          "IBS.sumstats.gz", "gast.sumstats.gz","oesp.sumstats.gz",
          "stmP.sumstats.gz", "crpl.sumstats.gz","cyst.sumstats.gz",
          "hipP.sumstats.gz", "kneP.sumstats.gz", "legP.sumstats.gz",
          "gout.sumstats.gz", "enLL.sumstats.gz", "hipA.sumstats.gz", 
          "kneA.sumstats.gz", "enth.sumstats.gz", "otRA.sumstats.gz",
          "arth.sumstats.gz", "pnjt.sumstats.gz", "genP.sumstats.gz")

# Output 2 includes: "MCPAIN.sumstats.gz"
#enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
sample.prev<-c(0.5, NA, 0.5, 0.5,
               NA, NA, NA,
               0.104, 0.102, 0.181,
               0.473, 0.167, 0.046,
               0.134, 0.188, 0.062,
               0.051, 0.027, 0.075,
               0.099, 0.190, 0.277,
               0.073, 0.035, 0.084, 
               0.145, 0.141, 0.042,
               0.339, 0.028, 0.014)

# Output 2 inlcudes: multisite pain 0.30
#vector of population prevalence
population.prev<-c(0.028, NA, 0.125, 0.048,
                   NA, NA, NA,
                   0.104, 0.102, 0.181,
                   0.473, 0.167, 0.046,
                   0.134, 0.188, 0.062,
                   0.051, 0.027, 0.075,
                   0.099, 0.190, 0.277,
                   0.073, 0.035, 0.084, 
                   0.145, 0.141, 0.042,
                   0.339, 0.028, 0.014)

# Output 2 inlcudes: multisite pain 0.30                 
#the folder of LD scores
ld<-"/Users/lydiarader/Documents/pain_research/GWAS_summary_statistics/dissertation/data/eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld<-"/Users/lydiarader/Documents/pain_research/GWAS_summary_statistics/dissertation/data/eur_w_ld_chr/"

#name the traits
trait.names<-c("OUD","PAU","TUD", "CanUD", 
               "AUDITc", "DPW", "CPD",
               "hdch","mgrn", "nksh", 
               "back","chDs", "chPh",
               "IBS", "gast","oesp",
               "stmP", "crpl","cyst",
               "hipP","kneP", "legP",
               "gout","enLL", "hipA", 
               "kneA","enth", "otRA",
               "arth",  "pnjt", "genP")

#Output 2 includes: "multiP"
#run LDSC
#test121324<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names, n.blocks = 600)



output<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names)

#output <- test121324
#Save output
save(output,file="output021325.RData")

####Load the LDSC output for GenomicSEM analyses####

setwd("~/Documents/pain_research/GWAS_summary_statistics/dissertation/data")
load("output021325.RData")
covstruc<-output #define covstruc (the covariance structure output from LDSC)

####Correlation plot for LDSC matrix####
#the S matrix
output$S # snp h2 on diagonal; covariances on off diagonal
cov2cor(output$S)

#the V matrix
output$V # the order of the V matrix is the order of the S matrix; squared standard error on diagonal but the OFF DIAGONAL is the sampling dependencies

#the LDCS intercept matrix
output$I #Diagonal are the univariate intercept and the OFF diagonal are bivariate intercepts
# In other words the OFF diagonal is = phenotypic correlation * the proportion of overlapping samples
# Important for multivariate GWAS
plot<-output
rownames(plot$S) <-colnames(plot$S)

#look at genetic correlation matrix
require(corrplot)
fig1<- corrplot(corr = cov2cor(plot$S),
         method = "color",
         addrect = 6,
         rect.col = 'red', 
         rect.lwd = 4,
         addCoef.col = "darkgrey",
         add = F,
         bg = "white",
         diag = T, # Autocorrelations
         outline = T,
         mar = c(0,0,0,5),
         cl.pos = "b",
         cl.ratio = 0.125,
         cl.align.text = "l",
         cl.offset = 0.2,
         tl.srt=45,
         tl.pos = 'full', #position of text label
         tl.offset=0.2,
         tl.col = "black",
         pch.col = "white",
         addgrid.col = "black",
         xpd = T,
         tl.cex=1,
         number.cex = 15/ncol(plot$S),
         is.corr=TRUE)

####Pain Model Measurement Model####
pain <-
"
# Specify Bi-factor Pain measurement model 
pain1 =~ hdch + mgrn + nksh + back + chDs + chPh + IBS + gast +
oesp + stmP + crpl + cyst + hipP + kneP + legP + gout + enLL + 
hipA + kneA + enth + otRA + arth + pnjt + genP

pain2 =~ crpl + hipP + kneP + legP + enLL + hipA + kneA + enth + otRA + arth + pnjt

# Specify residual covariances cross-pain conditions
chDs ~~ chPh
hipA ~~ hipP
kneA ~~ kneP
mgrn ~~ hdch


# Specify condition-level residual variances
arth ~~ r1*arth
r1 > 0.001


# Fix Pain1 and Pain2 factor correlations to 0 as this is a bi-factor model
pain1 ~~ 0*pain2
"
pain_results <-usermodel(covstruc=covstruc,estimation="DWLS",model=pain, std.lv = T)
pain_results


####Pain Model Measurement Model with Multisite Pain####
pain_test <-
  "
# Specify Bi-factor Pain measurement model 
pain1 =~ hdch + mgrn + nksh + back + chDs + chPh + IBS + gast +
oesp + stmP + crpl + cyst + hipP + kneP + legP + gout + enLL + 
hipA + kneA + enth + otRA + arth + pnjt + genP

pain2 =~ crpl + hipP + kneP + legP + enLL + hipA + kneA + enth + otRA + arth + pnjt

# Specify residual covariances cross-pain conditions
chDs ~~ chPh
hipA ~~ hipP
kneA ~~ kneP
mgrn ~~ hdch


# Specify condition-level residual variances
arth ~~ r1*arth
r1 > 0.001


# Fix Pain1 and Pain2 factor correlations to 0 as this is a bi-factor model
pain1 ~~ 0*pain2

pain1 ~~ multiP
pain2 ~~ multiP
"
pain_test_results <-usermodel(covstruc=covstruc,estimation="DWLS",model=pain_test, std.lv = T)
pain_test_results


####Addiction Measurement Model####
SUD <-
"
SUD=~ OUD + PAU + TUD + CanUD
"

SUD_results <-usermodel(covstruc=covstruc,estimation="DWLS",model=SUD,std.lv=T)
SUD_results


####Individual correlations between substance use disorder traits and 2 pain factors ####
pain_corr_su <-
"
# Specify Bi-factor Pain measurement model 
pain1 =~ hdch + mgrn + nksh + back + chDs + chPh + IBS + gast +
oesp + stmP + crpl + cyst + hipP + kneP + legP + gout + enLL + 
hipA + kneA + enth + otRA + arth + pnjt + genP

pain2 =~ crpl + hipP + kneP + legP + enLL + hipA + kneA + enth + otRA + arth + pnjt

# Specify residual covariances cross-pain conditions
chDs ~~ chPh
hipA ~~ hipP
kneA ~~ kneP
mgrn ~~ hdch

# Specify condition-level residual variances
arth ~~ r1*arth
r1 > 0.001

# Fix Pain1 and Pain2 factor correlations to 0 as this is a bi-factor model
pain1 ~~ 0*pain2

# Descriptive correlations between pain factors and substance use traits
pain1 ~~ OUD + PAU + TUD + CanUD + AUDITc + DPW + CPD
pain2 ~~ OUD + PAU + TUD + CanUD + AUDITc + DPW + CPD
"
pain_corr_su_results <-usermodel(covstruc=covstruc,estimation="DWLS",model=pain_corr_su, std.lv = T)
pain_corr_su_results

#### Model 1: Does OUD share genetic variance with General Pain & MSK pain factor over and above the Addiction factor ####
model1_oud <- "
SUD=~ OUD + PAU + TUD + CanUD

# OUD paths added to pain factors
pain1 =~ hdch + mgrn + nksh + back + chDs + chPh + IBS + gast +
  oesp + stmP + crpl + cyst + hipP + kneP + legP + gout + enLL + 
  hipA + kneA + enth + otRA + arth + pnjt + genP + OUD

pain2 =~ crpl + hipP + kneP + legP + enLL + hipA + kneA + enth + otRA + arth + pnjt + OUD 

chDs ~~ chPh
hipA ~~ hipP
kneA ~~ kneP
mgrn ~~ hdch

arth ~~ r1*arth
r1 > 0.001

pain1 ~~ 0*pain2

# Controlling for common Addiction pathways
SUD ~~ pain1
SUD ~~ pain2
"

model1oud_results <-usermodel(covstruc=output,estimation="DWLS",model=model1_oud,std.lv=T)
model1oud_results

#### Model 1: Does PAU share genetic variance with the General Pain & MSK factor over and above the Addiction factor ####
model1_pau <- "
SUD=~ OUD + PAU + TUD + CanUD

pain1 =~ hdch + mgrn + nksh + back + chDs + chPh + IBS + gast +
  oesp + stmP + crpl + cyst + hipP + kneP + legP + gout + enLL + 
  hipA + kneA + enth + otRA + arth + pnjt + genP + PAU

pain2 =~ crpl + hipP + kneP + legP + enLL + hipA + kneA + enth + otRA + arth + pnjt + PAU

chDs ~~ chPh
hipA ~~ hipP
kneA ~~ kneP
mgrn ~~ hdch

arth ~~ r1*arth
r1 > 0.001

pain1 ~~ 0*pain2

SUD ~~ pain1
SUD ~~ pain2
"
model1pau_results <-usermodel(covstruc=output,estimation="DWLS",model=model1_pau,std.lv=T)
model1pau_results

#### Model 1: Does TUD share genetic variance with the General Pain & MSK factor over and above the Addiction factor ####
model1_tud <- "
SUD=~ OUD + PAU + TUD + CanUD

pain1 =~ hdch + mgrn + nksh + back + chDs + chPh + IBS + gast +
  oesp + stmP + crpl + cyst + hipP + kneP + legP + gout + enLL + 
  hipA + kneA + enth + otRA + arth + pnjt + genP + TUD

pain2 =~ crpl + hipP + kneP + legP + enLL + hipA + kneA + enth + otRA + arth + pnjt + TUD

chDs ~~ chPh
hipA ~~ hipP
kneA ~~ kneP
mgrn ~~ hdch

arth ~~ r1*arth
r1 > 0.001

pain1 ~~ 0*pain2

SUD ~~ pain1
SUD ~~ pain2
"

model1tud_results <-usermodel(covstruc=output,estimation="DWLS",model=model1_tud,std.lv=T)
model1tud_results

#### Model 1: Does CanUD share genetic variance with the General Pain & MSK factor over and above the Addiction factor ####
model1_cud <- "
SUD=~ OUD + PAU + TUD + CanUD

pain1 =~ hdch + mgrn + nksh + back + chDs + chPh + IBS + gast +
  oesp + stmP + crpl + cyst + hipP + kneP + legP + gout + enLL + 
  hipA + kneA + enth + otRA + arth + pnjt + genP + CanUD

pain2 =~ crpl + hipP + kneP + legP + enLL + hipA + kneA + enth + otRA + arth + pnjt + CanUD

chDs ~~ chPh
hipA ~~ hipP
kneA ~~ kneP
mgrn ~~ hdch

arth ~~ r1*arth
r1 > 0.001

pain1 ~~ 0*pain2

SUD ~~ pain1
SUD ~~ pain2
"
model1cud_results <-usermodel(covstruc=output,estimation="DWLS",model=model1_cud,std.lv=T)
model1cud_results

####FDR-correction of p-values####
# Sample p-values
p_values <- c(7.976409e-87, #Add-GP
              7.761685e-02, #Add-MSK
              3.481367e-01, #OUD-GP
              1.905051e-01, #OUD-MSK
              7.040298e-84, #Add-GP
              2.431723e-02, #Add-MSK
              1.136622e-02, #PAU-GP
              1.343752e-02, #PAU-MSK
              3.412841e-65, #Add-GP
              8.213891e-01, #Add-MSK
              9.426570e-02, #TUD-GP
              1.141013e-03, #TUD-MSK
              1.256227e-80, #Add-GP
              1.315742e-01, #Add-MSK
              9.122435e-01, #CanUD-GP
              7.385403e-01) #CanUD-MSK
# Apply FDR correction
adjusted_p_values <- p.adjust(p_values, method = "fdr")
data.frame(
  Original_p_values = p_values,
  Adjusted_p_values = adjusted_p_values)

#### Model 2: Final model of common path and unique substances ####
model2 <- "
SUD=~ OUD + PAU + TUD + CanUD

pain1 =~ hdch + mgrn + nksh + back + chDs + chPh + IBS + gast +
  oesp + stmP + crpl + cyst + hipP + kneP + legP + gout + enLL + 
  hipA + kneA + enth + otRA + arth + pnjt + genP + PAU 

pain2 =~ crpl + hipP + kneP + legP + enLL + hipA + kneA + enth + otRA + arth + pnjt + PAU + TUD

chDs ~~ chPh
hipA ~~ hipP
kneA ~~ kneP
mgrn ~~ hdch

arth ~~ r1*arth
r1 > 0.001

pain1 ~~ 0*pain2

SUD ~~ pain1
SUD ~~ pain2
"

model2_results <-usermodel(covstruc=output,estimation="DWLS",model=model2,std.lv=T)
model2_results

####FDR-correction of p-values####
# Sample p-values
p_values <- c(2.435436e-85, #Add-GP
              9.877768e-01, #Add-MSK
              1.214363e-02, #PAU-GP
              4.317331e-01, #PAU-MSK
              1.790257e-03) #TUD-MSK
# Apply FDR correction
adjusted_p_values <- p.adjust(p_values, method = "fdr")
data.frame(
  Original_p_values = p_values,
  Adjusted_p_values = adjusted_p_values)


#### Model 3: Addiction factor expanded with Choleskys Alcohol ####
model3_alc <- "
# Cholesky
A1=~ PAU + AUDITc + DPW
A2 =~ AUDITc + DPW
A3 =~ DPW

PAU ~~ 0*PAU
AUDITc ~~ 0*AUDITc
DPW ~~ 0*DPW

A1 ~~ 0*A2
A1 ~~ 0*A3
A2 ~~ 0*A3

pain1 =~ hdch + mgrn + nksh + back + chDs + chPh + IBS + gast +
  oesp + stmP + crpl + cyst + hipP + kneP + legP + gout + enLL + 
  hipA + kneA + enth + otRA + arth + pnjt + genP 

pain2 =~ crpl + hipP + kneP + legP + enLL + hipA + kneA + enth + otRA + arth + pnjt + TUD

chDs ~~ chPh
hipA ~~ hipP
kneA ~~ kneP
mgrn ~~ hdch

arth ~~ r1*arth
r1 > 0.001

SUD=~ OUD + PAU + AUDITc + DPW + TUD + CanUD 

pain1 ~~ 0*pain2

A1 ~~ pain1 + pain2
A2 ~~ pain1 + pain2
A3 ~~ pain1 + pain2

SUD ~~ pain1
SUD ~~ pain2

SUD ~~ 0*A1
SUD ~~ 0*A2
SUD ~~ 0*A3
"
model3alc_results <-usermodel(covstruc=output,estimation="DWLS",model=model3_alc,std.lv=T)
model3alc_results

#save(model3alc_results, file="~/Documents/pain_research/GWAS_summary_statistics/dissertation/data/output/gsem_cholesky_alcohol.Rdata")
#### Model 3: Addiction factor expanded with Choleskys Tobacco####
model3_tob <- "
T1=~ TUD + CPD
T2=~ CPD
TUD ~~ 0*TUD
CPD ~~ 0*CPD

T1 ~~ 0*T2

pain1 =~ hdch + mgrn + nksh + back + chDs + chPh + IBS + gast +
  oesp + stmP + crpl + cyst + hipP + kneP + legP + gout + enLL + 
  hipA + kneA + enth + otRA + arth + pnjt + genP + PAU
  
T1 ~~ pain1 + pain2
T2 ~~ pain1 + pain2

pain2 =~ crpl + hipP + kneP + legP + enLL + hipA + kneA + enth + otRA + arth + pnjt + PAU

chDs ~~ chPh
hipA ~~ hipP
kneA ~~ kneP
mgrn ~~ hdch

arth ~~ r1*arth
r1 > 0.001

SUD=~ OUD + PAU + TUD + CPD + CanUD 

pain1 ~~ 0*pain2

SUD ~~ pain1
SUD ~~ pain2

SUD ~~ 0*T1
SUD ~~ 0*T2
"
model3tob_results <-usermodel(covstruc=output,estimation="DWLS",model=model3_tob,std.lv=T)
model3tob_results

#save(model3tob_results, file="~/Documents/pain_research/GWAS_summary_statistics/dissertation/data/output/gsem_cholesky_tobacco.Rdata")

####Model 3: Cholesky decomposition FDR-correction of p-values####
p_values <- c(1.977047e-02, #A1P1
              4.295728e-01, #A1P2
              1.112624e-07, #A2P1
              2.458824e-02, #A2P2
              4.083000e-01, #A3P1
              2.295114e-01, #A3P2
              4.626273e-01, #T1P1
              5.569711e-03, #T1P2
              9.922356e-17, #T2P1
              1.133185e-01) #T2P2
adjusted_p_values <- p.adjust(p_values, method = "fdr")
data.frame(
  Original_p_values = p_values,
  Adjusted_p_values = adjusted_p_values)

#### Model 3: Addiction factor expanded with Choleskys Tobacco AND Alcohol ####
model3_total <- "
# Alc Cholesky
A1=~ PAU + AUDITc + DPW
A2 =~ AUDITc + DPW
A3 =~ DPW

PAU ~~ 0*PAU
AUDITc ~~ 0*AUDITc
DPW ~~ 0*DPW

A1 ~~ 0*A2
A1 ~~ 0*A3
A2 ~~ 0*A3

# Tob cholesky
T1=~ TUD + CPD
T2=~ CPD
TUD ~~ 0*TUD
CPD ~~ 0*CPD

T1 ~~ 0*T2

A1 ~~ 0*T1
A1 ~~ 0*T2
A2 ~~ 0*T1
A2 ~~ 0*T2
A3 ~~ 0*T1
A3 ~~ 0*T2

pain1 =~ hdch + mgrn + nksh + back + chDs + chPh + IBS + gast +
  oesp + stmP + crpl + cyst + hipP + kneP + legP + gout + enLL + 
  hipA + kneA + enth + otRA + arth + pnjt + genP 

pain2 =~ crpl + hipP + kneP + legP + enLL + hipA + kneA + enth + otRA + arth + pnjt

chDs ~~ chPh
hipA ~~ hipP
kneA ~~ kneP
mgrn ~~ hdch

arth ~~ r1*arth
r1 > 0.001

T1 ~~ pain1 + pain2
T2 ~~ pain1 + pain2
A1 ~~ pain1 + pain2
A2 ~~ pain1 + pain2
A3 ~~ pain1 + pain2

SUD=~ OUD + PAU + AUDITc + DPW + TUD + CPD + CanUD 

pain1 ~~ 0*pain2

SUD ~~ pain1
SUD ~~ pain2

SUD ~~ 0*T1
SUD ~~ 0*T2
SUD ~~ 0*A1
SUD ~~ 0*A2
SUD ~~ 0*A3
"
model3all_results <-usermodel(covstruc=output,estimation="DWLS",model=model3_total,std.lv=T)
model3all_results

####Load Output for Stratified Genomic SEM####

# Conducted LDSC on a cluster and saved the RData as GSEM_sldsc.RData
require(GenomicSEM)
strat <- load("/Users/lydiarader/Documents/pain_research/GWAS_summary_statistics/dissertation/GSEM_sldsc.RData")

####Stratified GSEM: Main model of common and substance-specific associations between pain and addiction####
strat.pain.sud <- "
SUD=~ OUD + PAU + TUD + CUD

pain1 =~ hdch + mgrn + nksh + back + chDs + chPh + IBS + gast + oesp + stmP + crpl + cyst + hipP + kneP + legP + gout + enLL + hipA + kneA + enth + otRA + arth + pnjt + genP #moving PAU down below as stratified Genomic SEM likes this syntax better

pain2 =~ crpl + hipP + kneP + legP + enLL + hipA + kneA + enth + otRA + arth + pnjt + PAU #moving TUD down below as stratified Genomic SEM likes this syntax better

chDs ~~ chPh
hipA ~~ hipP
kneA ~~ kneP
mgrn ~~ hdch

arth ~~ r1*arth
r1 > 0.001

pain1 ~~ 0*pain2

SUD ~~ pain1
SUD ~~ pain2
PAU ~ pain1
TUD ~ pain2
"

#specify and fix parameters
params<-c("SUD ~~ pain1", "PAU ~ pain1", "TUD ~ pain2")
fixparam<-c("pain1~~pain1", "pain2~~pain2", "pain1~~pain2", "SUD~~pain2","SUD~~SUD","OUD~~OUD","PAU~~PAU", "TUD~~TUD", "CUD~~CUD","arth~~arth", "back ~~ back", "chDs ~~ chDs", "chPh ~~ chPh", "crpl ~~ cprl", "cyst ~~ cyst", "enLL ~~ enLL", "enth ~~ enth", "gast ~~ gast", "genP ~~ genP", "gout ~~ gout", "hdch ~~ hdch", "mgrn ~~ mgrn", "hipA ~~ hipA", "hipP ~~ hipP", "IBS ~~ IBS", "kneA ~~ kneA", "kneP ~~ kneP", "legP ~~ legP", "nksh ~~ nksh", "oesp ~~ oesp", "otRA ~~ otRA", "pnjt ~~ pnjt", "stmP ~~ stmP","crpl ~~ crpl", "hipP ~~ hipA", "kneP ~~ kneA","chDs ~~ chPh","hdch ~~ mgrn")
#run enrichment analysis
Enrich_results_pain_sud<-enrich(s_covstruc = GSEM_sldsc, model=strat.pain.sud,
                              params=params, fixparam=fixparam, std.lv=TRUE)
Enrich_results_pain_sud

####Visualize the General Pain - Addiction correlation####
#get top annotation names for SUD ~~ P1 relationship 
sud1p1_10_enrich<-Enrich_results_pain_sud[[1]] 
sud1p1_10_enrich$fdr_p <- p.adjust(sud1p1_10_enrich$Enrichment_p_value, method = "fdr") # Estimate FDR-corrected p-values
sud1p1_10_save<-sud1p1_10_enrich # save whole df for supp results
sud1p1_10_enrich <- sud1p1_10_enrich %>%
  filter(Warning==0) %>% # 75
  filter(fdr_p <= 5e-02) %>% # 58
  arrange(., fdr_p) %>% # Arrange by fdr_p significance
  head(10) # For the figures we will use the top ten enrichments
#names <- sud1p1_10_enrich %>% select(Annotation)

sud1p1_10_enrich$Annotation <- case_when(
  sud1p1_10_enrich$Annotation == "NewAnnotations/GABA1." ~ "GABAergic Neurons [Subset 1]",
  sud1p1_10_enrich$Annotation == "NewAnnotations/END." ~ "Endothelial Cells",
  sud1p1_10_enrich$Annotation == "NewAnnotations/exPFC1." ~ "Excitatory Prefrontal Cortex Neurons [Subset 1]",
  sud1p1_10_enrich$Annotation == "Conserved_Primate_phastCons46wayL2" ~ "Conserved Primate",
  sud1p1_10_enrich$Annotation == "NewAnnotations/PI_genes." ~ "PI Genes",
  sud1p1_10_enrich$Annotation == "NewAnnotations/OPC." ~ "Oligodendrocyte Precursor Cells",
  sud1p1_10_enrich$Annotation == "NewAnnotations/ASC1." ~ "Astrocytic Transporter 1",
  sud1p1_10_enrich$Annotation == "NewAnnotations/exCA1." ~ "Excitatory CA1 Hippocampal Neurons",
  sud1p1_10_enrich$Annotation == "NewAnnotations/exPFC2." ~ "PI x Excitatory Prefrontal Cortex Neurons [Subset 2]",
  sud1p1_10_enrich$Annotation == "NewAnnotations/PIxNSC." ~ "PI x Neuronal Stem Cells",
  TRUE ~ sud1p1_10_enrich$Annotation,)

ggplot(sud1p1_10_enrich, aes(x = Enrichment, y = reorder(Annotation, Enrichment))) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +  # Bar aesthetics
  labs(
    title = "Top Ten Significant Enrichments of the Addiction and General Chronic Pain Association",
    x = "Enrichment (Beta)",
    y = "Annotation"
  ) +
  theme_minimal() +  # Clean and modern theme
  theme(
    axis.text.y = element_text(size = 14, hjust = 1),  # Ensures readable annotation names
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  # Center and style the title
  ) +
  geom_errorbarh(aes(xmin = Enrichment - 1.96 * Enrichment_SE, 
                     xmax = Enrichment + 1.96 * Enrichment_SE, 
                     y = reorder(Annotation, Enrichment)),
                 linewidth = 0.5, height = 0.5, color = "black") +  # Error bar settings
  geom_vline(xintercept = 1, lty = 2) +  # Vertical reference line
  guides(alpha = FALSE) + 
  labs(color = NULL, fill = NULL) +  
  ylab("") +  
  xlab("Enrichment")  +
  geom_text(
    data = subset(sud1p1_10_enrich, fdr_p < 0.05),  # Filter only rows where 'fdr_p' < 0.05
    aes(
      x = Enrichment + 1.96 * Enrichment_SE,  # Position at the end of the error bar
      y = reorder(Annotation, Enrichment),
      label = "*"
    ),
    color = "red",  # Red color for the asterisk
    size = 5,  # Adjust the size of the asterisk
    hjust = -1  # Slightly offset the horizontal alignment
  )

#### Visualize the PAU - General Pain specific association ####
# [2] PAU ~ pain1 relationship
PAUp1_10_enrich<-Enrich_results_pain_sud[[2]] 
PAUp1_10_enrich$fdr_p <- p.adjust(PAUp1_10_enrich$Enrichment_p_value, method = "fdr") # Estimate FDR-corrected p-values
PAUp1_10_save <- PAUp1_10_enrich
PAUp1_10_enrich <- PAUp1_10_enrich %>%
  filter(Warning==0) %>%
  filter(fdr_p <= 5e-02) %>% # 52 significant hits, if we wanted to filter by FDR
  arrange(., fdr_p) %>% # Arrange by fdr_p significance
  head(10)
#names <- PAUp1_10_enrich %>% select(Annotation)
PAUp1_10_enrich$Annotation <- case_when(
  PAUp1_10_enrich$Annotation == "NewAnnotations/END." ~ "Endothelial Cells",
  PAUp1_10_enrich$Annotation == "NewAnnotations/exPFC1." ~ "Excitatory Prefrontal Cortex Neurons [Subset 1]",
  PAUp1_10_enrich$Annotation == "NewAnnotations/GABA1." ~ "GABAergic Neurons [Subset 1]",
  PAUp1_10_enrich$Annotation == "NewAnnotations/exPFC2." ~ "PI x Excitatory Prefrontal Cortex Neurons [Subset 2]",
  PAUp1_10_enrich$Annotation == "NewAnnotations/PI_genes." ~ "PI Genes",
  PAUp1_10_enrich$Annotation == "NewAnnotations/ASC1." ~ "Astrocytic Transporter 1",
  PAUp1_10_enrich$Annotation == "NewAnnotations/ASC2." ~ "Astrocytic Transporter 2",
  PAUp1_10_enrich$Annotation == "NewAnnotations/OPC." ~ "Oligodendrocyte Precursor Cells",
  PAUp1_10_enrich$Annotation == "NewAnnotations/exCA3." ~ "Excitatory CA3 Hippocampal Neurons",
  PAUp1_10_enrich$Annotation == "NewAnnotations/exCA1." ~ "Excitatory CA1 Hippocampal Neurons",
  TRUE ~ PAUp1_10_enrich$Annotation,)

ggplot(PAUp1_10_enrich, aes(x = Enrichment, y = reorder(Annotation, Enrichment))) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +  # Bar aesthetics
  labs(
    title = "Top Ten Significant Enrichments of the PAU-specific and General Chronic Pain Association",
    x = "Enrichment (Beta)",
    y = "Annotation"
  ) +
  theme_minimal() +  # Clean and modern theme
  theme(
    axis.text.y = element_text(size = 14, hjust = 1),  # Ensures readable annotation names
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  # Center and style the title
  ) +
  geom_errorbarh(aes(xmin = Enrichment - 1.96 * Enrichment_SE, 
                     xmax = Enrichment + 1.96 * Enrichment_SE, 
                     y = reorder(Annotation, Enrichment)),
                 linewidth = 0.5, height = 0.5, color = "black") +  # Error bar settings
  geom_vline(xintercept = 1, lty = 2) +  # Vertical reference line
  guides(alpha = FALSE) + 
  labs(color = NULL, fill = NULL) +  
  ylab("") +  
  xlab("Enrichment")  +
  geom_text(
    data = subset(PAUp1_10_enrich, fdr_p < 0.05),  # Filter only rows where 'fdr_p' < 0.05
    aes(
      x = Enrichment + 1.96 * Enrichment_SE,  # Position at the end of the error bar
      y = reorder(Annotation, Enrichment),
      label = "*"
    ),
    color = "red",  # Red color for the asterisk
    size = 5,  # Adjust the size of the asterisk
    hjust = -1  # Slightly offset the horizontal alignment
  )


#### Visualize the TUD - MSK specific association ####
# [3] TUD ~ pain2 relationship
TUDp2_10_enrich<-Enrich_results_pain_sud[[3]] 
TUDp2_10_enrich$fdr_p <- p.adjust(TUDp2_10_enrich$Enrichment_p_value, method = "fdr") # Estimate FDR-corrected p-values
TUDp2_10_save <- TUDp2_10_enrich
TUDp2_10_enrich <- TUDp2_10_enrich %>%
  filter(Warning==0) %>%
  filter(fdr_p <= 5e-02) %>% # 17 significant hits, if we wanted to filter by FDR
  arrange(., fdr_p) %>% # Arrange by fdr_p significance
  head(10) 
#names <- TUDp2_10_enrich %>% select(Annotation)
TUDp2_10_enrich$Annotation <- case_when(
  TUDp2_10_enrich$Annotation == "NewAnnotations/PIxEND." ~ "PI x Endothelial Cells",
  TUDp2_10_enrich$Annotation == "NewAnnotations/PIxNSC." ~ "PI x Neuronal Stem Cells",
  TUDp2_10_enrich$Annotation == "NewAnnotations/PIxASC2." ~ "PI x Astrocytic Transporter 2",
  TUDp2_10_enrich$Annotation == "NewAnnotations/PIxODC1." ~ "PI x Oligodendrocytes",
  TUDp2_10_enrich$Annotation == "UTR_5_UCSCL2" ~ "UTR 5 USC",
  TUDp2_10_enrich$Annotation == "NewAnnotations/PIxMG." ~ "PI x Microglia",
  TUDp2_10_enrich$Annotation == "NewAnnotations/PIxexCA1." ~ "PI x Excitatory CA1 Hippocampal Neurons",
  TUDp2_10_enrich$Annotation == "NewAnnotations/PIxOPC." ~ "PI x Oligodendrocyte Precursor Cells",
  TUDp2_10_enrich$Annotation == "Human_Promoter_Villar_ExACL2" ~ "Human Promoter ExAC",
  TUDp2_10_enrich$Annotation == "synonymousL2" ~ "Synonymous",
  TRUE ~ TUDp2_10_enrich$Annotation,)

ggplot(TUDp2_10_enrich, aes(x = Enrichment, y = reorder(Annotation, Enrichment))) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +  # Bar aesthetics
  labs(
    title = "Top Ten Significant Enrichments of the TUD-specific and Musculoskeletal Pain Association",
    x = "Enrichment (Beta)",
    y = "Annotation"
  ) +
  theme_minimal() +  # Clean and modern theme
  theme(
    axis.text.y = element_text(size = 14, hjust = 1),  # Ensures readable annotation names
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  # Center and style the title
  ) +
  geom_errorbarh(aes(xmin = Enrichment - 1.96 * Enrichment_SE, 
                     xmax = Enrichment + 1.96 * Enrichment_SE, 
                     y = reorder(Annotation, Enrichment)),
                 linewidth = 0.5, height = 0.5, color = "black") +  # Error bar settings
  geom_vline(xintercept = 1, lty = 2) +  # Vertical reference line
  guides(alpha = FALSE) + 
  labs(color = NULL, fill = NULL) +  
  ylab("") +  
  xlab("Enrichment")  +
  geom_text(
    data = subset(TUDp2_10_enrich, fdr_p < 0.05),  # Filter only rows where 'fdr_p' < 0.05
    aes(
      x = Enrichment + 1.96 * Enrichment_SE,  # Position at the end of the error bar
      y = reorder(Annotation, Enrichment),
      label = "*"
    ),
    color = "red",  # Red color for the asterisk
    size = 5,  # Adjust the size of the asterisk
    hjust = -1  # Slightly offset the horizontal alignment
  )

#### Stratified GSEM: Cholesky Decomposition of alcohol and consumption ####
strat.cd.alc <-
"
# Addiction factor w additional alcohol consumption indicators
SUD=~ OUD + PAU + AUDc + DPW + TUD + CUD 

# Alcohol Cholesky decomposition
A1=~ PAU + AUDc + DPW # Problematic Alcohol Use factor
A2 =~ AUDc + DPW      # Alcohol Consumption factor
A3 =~ DPW             # Alcohol Frequency Residual factor

PAU ~~ 0*PAU          # Set residuals to 0 so fully captured by the factors
AUDc ~~ 0*AUDc        # Note the indicators have variance left NOT captured by Addiction factor
DPW ~~ 0*DPW

A1 ~~ 0*A2            # No covariances among cholesky factors
A1 ~~ 0*A3
A2 ~~ 0*A3

SUD ~~ 0*A1           # No covariances among Addiction factor and cholesky factors
SUD ~~ 0*A2
SUD ~~ 0*A3

# Pain factor model
pain1 =~ hdch + mgrn + nksh + back + chDs + chPh + IBS + gast +
  oesp + stmP + crpl + cyst + hipP + kneP + legP + gout + enLL + 
  hipA + kneA + enth + otRA + arth + pnjt + genP 

# Add TUD to pain2 based on our final model
pain2 =~ crpl + hipP + kneP + legP + enLL + hipA + kneA + enth + otRA + arth + pnjt + TUD

chDs ~~ chPh
hipA ~~ hipP
kneA ~~ kneP
mgrn ~~ hdch

arth ~~ r1*arth
r1 > 0.001

pain1 ~~ 0*pain2

# Covariances controlling for common path through Addiction
SUD ~~ pain1
SUD ~~ pain2

# Covariances of interest: How much do the pain factors associate w problematic alcohol use and consumption AFTER controlling for common path
A1 ~~ pain1 + pain2     
A2 ~~ pain1 + pain2
A3 ~~ pain1 + pain2
"

#specify and fix parameters
params<-c("A1 ~~ pain1", "A2 ~~ pain1", "A2 ~~ pain2") # These 3 are the FDR-corrected significant parameters in Genomic SEM

fixparam<-c("pain1~~pain1", "pain2~~pain2", "pain1 ~~ pain2", "A1~~A1", "A2~~A2", "A3~~A3","A1~~A2","A1~~A3", "A2~~A3", "A3~~pain1","A3~~pain2","PAU~~PAU", "AUDc~~AUDc", "DPW~~DPW", "OUD~~OUD", "TUD~~TUD", "CUD~~CUD", "SUD~~SUD", "arth~~arth", "back ~~ back", "chDs ~~ chDs", "chPh ~~ chPh", "crpl ~~ cprl", "cyst ~~ cyst","enLL ~~ enLL", "enth ~~ enth", "gast ~~ gast", "genP ~~ genP", "gout ~~ gout", "hdch ~~ hdch", "mgrn ~~ mgrn", "hipA ~~ hipA", "hipP ~~ hipP", "IBS ~~ IBS", "kneA ~~ kneA", "kneP ~~ kneP", "legP ~~ legP", "nksh ~~ nksh", "oesp ~~ oesp", "otRA ~~ otRA", "pnjt ~~ pnjt", "stmP ~~ stmP","crpl ~~ crpl", "hipP ~~ hipA", "kneP ~~ kneA","chDs ~~ chPh","hdch ~~ mgrn", "SUD ~~ pain1", "SUD ~~ pain2", "SUD ~~ A1", "SUD ~~ A2", "SUD ~~ A3", "A1 ~~ pain2")

#run enrichment analysis
Enrich_results_cd_alc<-enrich(s_covstruc = GSEM_sldsc, model=strat.cd.alc,
                              params=params, fixparam=fixparam, std.lv=TRUE)
Enrich_results_cd_alc

####Visualize the Cholesky PAU - General Pain correlation####
#get top annotation names for
# [1] A1 ~~ pain 1
#filter by fdr.adjusted significant p=values
a1p1_10_enrich<-Enrich_results_cd_alc[[1]] 
a1p1_10_enrich$fdr_p <- p.adjust(a1p1_10_enrich$Enrichment_p_value, method = "fdr")
a1p1_10_save <- a1p1_10_enrich # saving full results
a1p1_10_enrich <- a1p1_10_enrich %>% # 80
  filter(Warning==0) %>% # 56 after filtering for models that failed to estimate
  filter(fdr_p <= 5e-02) %>% # 55 after filtering for fdr p
  arrange(., fdr_p) %>% 
  head(10)

# Rename the significant annotations to something more readable for the A1~~P1 association
a1p1_10_enrich$Annotation <- case_when(
  a1p1_10_enrich$Annotation =="Enhancer_AnderssonL2" ~ "Enhancer Andersson",
  a1p1_10_enrich$Annotation =="MAFbin1L2" ~ "MAF Bin 1: 5% - 7.1%",
  a1p1_10_enrich$Annotation =="synonymousL2" ~ "Synonymous",
  a1p1_10_enrich$Annotation =="non_synonymousL2" ~ "Non Synonymous",
  a1p1_10_enrich$Annotation =="Ancient_Sequence_Age_Human_EnhancerL2" ~ "Ancient Sequence Age Human Enhancer",
  a1p1_10_enrich$Annotation =="Ancient_Sequence_Age_Human_PromoterL2" ~ "Ancient Sequence Age Human Promoter",
  a1p1_10_enrich$Annotation =="Human_Promoter_Villar_ExACL2" ~ "Human Promoter ExAC",
  a1p1_10_enrich$Annotation =="NewAnnotations/exPFC2." ~ "Excitatory Prefrontal Cortex Neurons [Subset 2]",
  a1p1_10_enrich$Annotation =="NewAnnotations/NSC." ~ "Neuronal Stem Cells",
  a1p1_10_enrich$Annotation =="NewAnnotations/PIxexCA1." ~ "PI x Excitatory CA1 Hippocampal Neurons",
  TRUE ~ a1p1_10_enrich$Annotation,)

# Plot the top ten enrichments for the A1~~P1 association
ggplot(a1p1_10_enrich, aes(x = Enrichment, y = reorder(Annotation, Enrichment))) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +  # Bar aesthetics
  labs(
    title = "Top Ten Significant Enrichments of the Cholesky Problematic Alcohol Use factor and General Pain Association",
    x = "Enrichment (Beta)",
    y = "Annotation"
  ) +
  theme_minimal() +  # Clean and modern theme
  theme(
    axis.text.y = element_text(size = 14, hjust = 1),  # Ensures readable annotation names
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  # Center and style the title
  ) +
  geom_errorbarh(aes(xmin = Enrichment - 1.96 * Enrichment_SE, 
                     xmax = Enrichment + 1.96 * Enrichment_SE, 
                     y = reorder(Annotation, Enrichment)),
                 linewidth = 0.5, height = 0.5, color = "black") +  # Error bar settings
  geom_vline(xintercept = 1, lty = 2) +  # Vertical reference line
  guides(alpha = FALSE) + 
  labs(color = NULL, fill = NULL) +  
  ylab("") +  
  xlab("Enrichment")  +
  geom_text(
    data = subset(a1p1_10_enrich, fdr_p < 0.05),  # Filter only rows where 'fdr_p' < 0.05
    aes(
      x = Enrichment + 1.96 * Enrichment_SE,  # Position at the end of the error bar
      y = reorder(Annotation, Enrichment),
      label = "*"
    ),
    color = "red",  # Red color for the asterisk
    size = 5,  # Adjust the size of the asterisk
    hjust = -1  # Slightly offset the horizontal alignment
  )

# [2] A2 ~~ pain 1; alcohol consumption ~~ general pain
#get top 10 annotation names for A2 ~~ Pain 1 relationship
a2p1_10_enrich<-Enrich_results_cd_alc[[2]] 
a2p1_10_enrich$fdr_p <- p.adjust(a2p1_10_enrich$Enrichment_p_value, method = "fdr")
a2p1_10_save <- a2p1_10_enrich
a2p1_10_enrich <- a2p1_10_enrich %>%
  filter(Warning == 0) %>% # 56
  #filter(fdr_p <= 5e-02) %>% # 0
  arrange(., fdr_p) %>% # Let's get top 10 even though not significant for supp figure
  head(10)

# rename the top 10 annotations
a2p1_10_enrich$Annotation <- case_when(
  a2p1_10_enrich$Annotation == "baseL2" ~ "Base L2",
  a2p1_10_enrich$Annotation == "Conserved_LindbladTohL2" ~ "Conserved LinbladToh",
  a2p1_10_enrich$Annotation == "CTCF_HoffmanL2" ~ "CTCF Hoffman",
  a2p1_10_enrich$Annotation == "DGF_ENCODEL2" ~ "DGF Encode",
  a2p1_10_enrich$Annotation == "DHS_peaks_TrynkaL2" ~ "DHS Peaks",
  a2p1_10_enrich$Annotation == "DHS_TrynkaL2" ~ "DHS Trynka",
  a2p1_10_enrich$Annotation == "Enhancer_AnderssonL2" ~ "Enhancer Andersson",  
  a2p1_10_enrich$Annotation == "Enhancer_HoffmanL2" ~ "Enhancer Hoffman",
  a2p1_10_enrich$Annotation == "FetalDHS_TrynkaL2" ~ "Fetal DHS",
  a2p1_10_enrich$Annotation == "H3K27ac_HniszL2" ~ "H3K27ac",
  TRUE ~ a2p1_10_enrich$Annotation,)

# Plot the alcohol consumption - general pain association
ggplot(a2p1_10_enrich, aes(x = Enrichment, y = reorder(Annotation, Enrichment))) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +  # Bar aesthetics
  labs(
    title = "The Top Ten Enrichments of the Cholesky Alcohol Consumption Factor and General Pain Association",
    x = "Enrichment (Beta)",
    y = "Annotation"
  ) +
  theme_minimal() +  
  theme(
    axis.text.y = element_text(size = 14, hjust = 1),  
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold") 
  ) +
  geom_errorbarh(aes(xmin = Enrichment - 1.96 * Enrichment_SE, 
                     xmax = Enrichment + 1.96 * Enrichment_SE, 
                     y = reorder(Annotation, Enrichment)),
                 linewidth = 0.5, height = 0.5, color = "black") +  
  geom_vline(xintercept = 1, lty = 2) + 
  guides(alpha = FALSE) +  
  labs(color = NULL, fill = NULL) +  
  ylab("") + 
  xlab("Enrichment") +
  geom_text(
    data = subset(a2p1_10_enrich, fdr_p < 0.05), 
    aes(
      x = Enrichment + 1.96 * Enrichment_SE,  
      y = reorder(Annotation, Enrichment),
      label = "*"
    ),
    color = "red",  # Red color for the asterisk
    size = 5,  
    hjust = -1  
  )

# [3] A2 ~~ pain 2; alcohol consumption ~~ MSK pain
#get top 10 annotation names for A2 ~~ Pain 2 relationship
a2p2_10_enrich<-Enrich_results_cd_alc[[3]] 
a2p2_10_enrich$fdr_p <- p.adjust(a2p2_10_enrich$Enrichment_p_value, method = "fdr")
a2p2_10_save <- a2p2_10_enrich
a2p2_10_enrich <- a2p2_10_enrich %>%
  filter(Warning == 0) %>% # 56
  #filter(fdr_p <= 5e-02) %>% # 0
  arrange(., fdr_p) %>% # Let's get top 10 even though not significant for supp figure
  head(10)

# rename the annotation names
a2p2_10_enrich$Annotation <- case_when(
  a2p2_10_enrich$Annotation == "baseL2" ~ "Base L2",
  a2p2_10_enrich$Annotation == "Conserved_LindbladTohL2" ~ "Conserved LinbladToh",
  a2p2_10_enrich$Annotation == "CTCF_HoffmanL2" ~ "CTCF Hoffman",
  a2p2_10_enrich$Annotation == "DGF_ENCODEL2" ~ "DGF Encode",
  a2p2_10_enrich$Annotation == "DHS_peaks_TrynkaL2" ~ "DHS Peaks",
  a2p2_10_enrich$Annotation == "DHS_TrynkaL2" ~ "DHS Trynka",
  a2p2_10_enrich$Annotation == "Enhancer_AnderssonL2" ~ "Enhancer Andersson", 
  a2p2_10_enrich$Annotation == "Enhancer_HoffmanL2" ~ "Enhancer Hoffman",
  a2p2_10_enrich$Annotation == "FetalDHS_TrynkaL2" ~ "Fetal DHS",
  a2p2_10_enrich$Annotation == "H3K27ac_HniszL2" ~ "H3K27ac",
  TRUE ~ a2p2_10_enrich$Annotation,)

# plot the alcohol consumption musculoskeletal pain enrichments
ggplot(a2p2_10_enrich, aes(x = Enrichment, y = reorder(Annotation, Enrichment))) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +  # Bar aesthetics
  labs(
    title = "The Top Ten Enrichments of the Chlesky Alcohol Consumption Factor and Musculoskeletal Pain Association",
    x = "Enrichment (Beta)",
    y = "Annotation"
  ) +
  theme_minimal() +  
  theme(
    axis.text.y = element_text(size = 14, hjust = 1),  # Ensures readable annotation names
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  # Center and style the title
  ) +
  geom_errorbarh(aes(xmin = Enrichment - 1.96 * Enrichment_SE, 
                     xmax = Enrichment + 1.96 * Enrichment_SE, 
                     y = reorder(Annotation, Enrichment)),
                 linewidth = 0.5, height = 0.5, color = "black") +  # Error bar settings
  geom_vline(xintercept = 1, lty = 2) +  # Vertical reference line
  guides(alpha = FALSE) +  # Removes the alpha guide
  labs(color = NULL, fill = NULL) +  # Clean up legends
  ylab("") +  # Removes y-axis label
  xlab("Enrichment") +
  geom_text(
    data = subset(a2p2_10_enrich, fdr_p < 0.05),  # Filter only rows where 'fdr_p' < 0.05
    aes(
      x = Enrichment + 1.96 * Enrichment_SE,  # Position at the end of the error bar
      y = reorder(Annotation, Enrichment),
      label = "*"
    ),
    color = "red",  # Red color for the asterisk
    size = 5,  # Adjust the size of the asterisk
    hjust = -1  # Slightly offset the horizontal alignment
  )

#### Stratified GSEM: Cholesky Decomposition of smoking and consumption ####
strat.cd.tob<-
"
# Addiction factor w smoking frequency
SUD=~ OUD + PAU + TUD + CPD + CUD 

# Tobacco Cholesky decomposition
T1=~ TUD + CPD         # Tobacco use disorder factor
T2=~ CPD               # Tobacco consumption factor

TUD ~~ 0*TUD           # Set residual variance to 0
CPD ~~ 0*CPD           # Note: these indicators are all the variance AFTER controlling for addiction factor variance

T1 ~~ 0*T2

# No covariances with Addiction and cholesky factors
SUD ~~ 0*T1
SUD ~~ 0*T2

# Pain Factor
# Include PAU on pain1 since it was significant
pain1 =~ hdch + mgrn + nksh + back + chDs + chPh + IBS + gast +
  oesp + stmP + crpl + cyst + hipP + kneP + legP + gout + enLL + 
  hipA + kneA + enth + otRA + arth + pnjt + genP + PAU
  
pain2 =~ crpl + hipP + kneP + legP + enLL + hipA + kneA + enth + otRA + arth + pnjt + PAU

chDs ~~ chPh
hipA ~~ hipP
kneA ~~ kneP
mgrn ~~ hdch

arth ~~ r1*arth
r1 > 0.001

pain1 ~~ 0*pain2

# Covariance controlling for common path
SUD ~~ pain1
SUD ~~ pain2

# Covariances of interest
T1 ~~ pain1 + pain2 
T2 ~~ pain1 + pain2
"

#specify and fix parameters
params<-c("T1 ~~ pain2", "T2 ~~ pain1" ) # These were the 2 sig params in Genomic SEM
fixparam<-c("T1 ~~ pain1", "pain1~~pain1", "pain2~~pain2", "T1~~T1", "T2~~T2", "PAU~~PAU", "OUD~~OUD", "CUD~~CUD", "SUD~~SUD","arth ~~ arth", "back ~~ back", "chDs ~~ chDs", "chPh ~~ chPh", "chDs ~~ chPh", "CPD ~~ CPD", "crpl ~~ cprl", "cyst ~~ cyst", "enLL ~~ enLL", "enth ~~ enth", "gast ~~ gast", "genP ~~ genP", "gout ~~ gout", "hdch ~~ hdch", "hdch ~~ mgrn", "mgrn ~~ mgrn", "hipA ~~ hipA", "hipP ~~ hipP", "IBS ~~ IBS", "kneA ~~ kneA", "kneP ~~ kneP", "legP ~~ legP", "nksh ~~ nksh", "oesp ~~ oesp", "otRA ~~ otRA", "pain1 ~~ pain2", "pnjt ~~ pnjt", "stmP ~~ stmP", "TUD~~TUD", "CPD~~CPD", "T2~~pain2", "crpl ~~ crpl", "hipP ~~ hipA", "kneP ~~ kneA", "T1 ~~ T2", "SUD ~~ T1", "SUD ~~ T2", "SUD ~~ pain1", "SUD ~~ pain2")

#run enrichment analysis
Enrich_results_cd_tob<-enrich(s_covstruc = GSEM_sldsc, model=strat.cd.tob,
                              params=params, fixparam=fixparam, std.lv=TRUE)
Enrich_results_cd_tob


####Visualization of Stratified GSEM Cholesky Decomposition of smoking and consumption####
#get top annotation names for
# [1] T1 ~~ pain2
#filter by fdr.adjusted significant p=values
t1p2_10_enrich<-Enrich_results_cd_tob[[1]] 
t1p2_10_enrich$fdr_p <- p.adjust(t1p2_10_enrich$Enrichment_p_value, method = "fdr")
t1p2_10_save <- t1p2_10_enrich
t1p2_10_enrich <- t1p2_10_enrich %>% 
  filter(Warning==0) %>% # 69
  #filter(fdr_p <= 5e-02) %>%  #0
  arrange(., fdr_p) %>% 
  head(10)

# Rename the top ten annotations to something more readable for the T1~~P2 association
t1p2_10_enrich$Annotation <- case_when(
  t1p2_10_enrich$Annotation=="baseL2" ~ "Base L2",
  t1p2_10_enrich$Annotation == "Coding_UCSCL2" ~ "Coding USC",
  t1p2_10_enrich$Annotation == "Conserved_LindbladTohL2" ~ "Conserved LinbladToh",
  t1p2_10_enrich$Annotation == "CTCF_HoffmanL2" ~ "CTCF Hoffman",
  t1p2_10_enrich$Annotation == "DGF_ENCODEL2" ~ "DGF Encode",
  t1p2_10_enrich$Annotation == "DHS_peaks_TrynkaL2" ~ "DHS Peaks",
  t1p2_10_enrich$Annotation == "DHS_TrynkaL2" ~ "DHS Trynka",
  t1p2_10_enrich$Annotation == "Enhancer_AnderssonL2" ~ "Enhancer Andersson",
  t1p2_10_enrich$Annotation == "Enhancer_HoffmanL2" ~ "Enhancer Hoffman",
  t1p2_10_enrich$Annotation == "FetalDHS_TrynkaL2" ~ "Fetal DHS",
  TRUE ~ t1p2_10_enrich$Annotation,)

# plot the top enrichments of tobacco use disorder and MSK pain
ggplot(t1p2_10_enrich, aes(x = Enrichment, y = reorder(Annotation, Enrichment))) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +  # Bar aesthetics
  labs(
    title = "Top Ten Enrichments of the Cholesky Tobacco Use Disorder Factor and Musculoskeletal Pain Association",
    x = "Enrichment (Beta)",
    y = "Annotation"
  ) +
  theme_minimal() +  # Clean and modern theme
  theme(
    axis.text.y = element_text(size = 14, hjust = 1),  # Ensures readable annotation names
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  # Center and style the title
  ) +
  geom_errorbarh(aes(xmin = Enrichment - 1.96 * Enrichment_SE, xmax = Enrichment + 1.96 * Enrichment_SE, y = Annotation),
                 position = position_dodge(width = 0.8),  # Dodge by the factor variable
                 linewidth = 0.5, height = 0.5, color = "black")  +
  geom_vline(xintercept = 1, lty = 2) +
  guides(alpha = FALSE) +
  labs(color = NULL, fill = NULL)+
  ylab("")+
  xlab("Enrichment") +
  geom_text(
    data = subset(t1p2_10_enrich, fdr_p < 0.05),  # Filter only rows where 'fdr_p' < 0.05
    aes(
      x = Enrichment + 1.96 * Enrichment_SE,  # Position at the end of the error bar
      y = reorder(Annotation, Enrichment),
      label = "*"
    ),
    color = "red",  # Red color for the asterisk
    size = 5,  # Adjust the size of the asterisk
    hjust = -1  # Slightly offset the horizontal alignment
  )


# [2] T2 ~~ pain1
#get top annotation names for T2 ~~ P1 relationship
t2p1_10_enrich<-Enrich_results_cd_tob[[2]] 
t2p1_10_enrich$fdr_p <- p.adjust(t2p1_10_enrich$Enrichment_p_value, method = "fdr")
t2p1_10_save <- t2p1_10_enrich
t2p1_10_enrich <- t2p1_10_enrich %>% 
  filter(Warning == 0) %>% # 69
  #filter(fdr_p <= 5e-02) %>% # 0
  arrange(., fdr_p) %>% 
  head(10)


# Rename the significant annotations to something more readable for the T2~~P1 association
t2p1_10_enrich$Annotation <- case_when(
  t2p1_10_enrich$Annotation=="baseL2" ~ "Base L2",
  t2p1_10_enrich$Annotation == "Coding_UCSCL2" ~ "Coding USC",
  t2p1_10_enrich$Annotation == "Conserved_LindbladTohL2" ~ "Conserved LinbladToh",
  t2p1_10_enrich$Annotation == "CTCF_HoffmanL2" ~ "CTCF Hoffman",
  t2p1_10_enrich$Annotation == "DGF_ENCODEL2" ~ "DGF Encode",
  t2p1_10_enrich$Annotation == "DHS_peaks_TrynkaL2" ~ "DHS Peaks",
  t2p1_10_enrich$Annotation == "DHS_TrynkaL2" ~ "DHS Trynka",
  t2p1_10_enrich$Annotation == "Enhancer_AnderssonL2" ~ "Enhancer Andersson",
  t2p1_10_enrich$Annotation == "Enhancer_HoffmanL2" ~ "Enhancer Hoffman",
  t2p1_10_enrich$Annotation == "FetalDHS_TrynkaL2" ~ "Fetal DHS",
  TRUE ~ t2p1_10_enrich$Annotation,)

# plot the top enrichments for tob consumption and general pain
ggplot(t2p1_10_enrich, aes(x = Enrichment, y = reorder(Annotation, Enrichment))) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +  # Bar aesthetics
  labs(
    title = "Top Ten Enrichments of the Cholesky Tobacco Consumption Factor and General Pain Association",
    x = "Enrichment (Beta)",
    y = "Annotation"
  ) +
  theme_minimal() +  # Clean and modern theme
  theme(
    axis.text.y = element_text(size = 14, hjust = 1),  # Ensures readable annotation names
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  # Center and style the title
  ) +
  geom_errorbarh(aes(xmin = Enrichment - 1.96 * Enrichment_SE, xmax = Enrichment + 1.96 * Enrichment_SE, y = Annotation),
                 position = position_dodge(width = 0.8),  # Dodge by the factor variable
                 linewidth = 0.5, height = 0.5, color = "black")  +
  geom_vline(xintercept = 1, lty = 2) +
  guides(alpha = FALSE) +
  labs(color = NULL, fill = NULL)+
  ylab("")+
  xlab("Enrichment") +
  geom_text(
    data = subset(t2p1_10_enrich, fdr_p < 0.05),  # Filter only rows where 'fdr_p' < 0.05
    aes(
      x = Enrichment + 1.96 * Enrichment_SE,  # Position at the end of the error bar
      y = reorder(Annotation, Enrichment),
      label = "*"
    ),
    color = "red",  # Red color for the asterisk
    size = 5,  # Adjust the size of the asterisk
    hjust = -1  # Slightly offset the horizontal alignment
  )


#### Write spreadsheet for supplementary results ####
library(openxlsx)
setwd("~/Documents/pain_research/GWAS_summary_statistics/SUD_pain_genomicsem_2025/biopsych_drafts")
# For this, make sure all the stratified dataframes are NOT selected for only significant. Include ALL annotation results, then save.
# Create a new workbook
wb <- createWorkbook()
# Add worksheets
addWorksheet(wb, "AddictionPain") # Addiction and General pain
writeData(wb, "AddictionPain", sud1p1_10_save, colNames = TRUE)

addWorksheet(wb, "PAUPain") # PAU and General Pain
writeData(wb, "PAUPain", PAUp1_10_save, colNames = TRUE)

addWorksheet(wb, "TUDMSK") # TUD and MSK
writeData(wb, "TUDMSK", TUDp2_10_save, colNames = TRUE)

addWorksheet(wb, "CholeskyPAUPain") # PAU and General Pain
writeData(wb, "CholeskyPAUPain", a1p1_10_save, colNames = TRUE)

addWorksheet(wb, "CholeskyAlcConsPain") # Alc Consumption and General Pain
writeData(wb, "CholeskyAlcConsPain", a2p1_10_enrich, colNames = TRUE)

addWorksheet(wb, "CholeskyAlcConsMSK") # Alc Consumption and MSK
writeData(wb, "CholeskyAlcConsMSK", a2p2_10_enrich, colNames = TRUE)

addWorksheet(wb, "CholeskyTUDMSK") # TUD and MSK
writeData(wb, "CholeskyTUDMSK", t1p2_10_save, colNames = TRUE)

addWorksheet(wb, "CholeskyTobConsPain") # Tob Consumption and General Pain
writeData(wb, "CholeskyTobConsPain", t2p1_10_enrich, colNames = TRUE)

# Save the workbook 
saveWorkbook(wb, "supplementary_results.xlsx", overwrite = TRUE)


#### ARCHIVE SECTION ####

