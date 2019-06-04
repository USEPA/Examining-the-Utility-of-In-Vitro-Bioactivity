#-----------------------------------------------------------------------------------#
# APCRA case study: Make integrated dataset
# 
# Katie Paul Friedman, paul-friedman.katie@epa.gov
# Original 2 Apr 2018
# Updated 6 Apr 2018
#-----------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#
# loading libraries and sources
#-----------------------------------------------------------------------------------#

rm(list = ls())

library(ggplot2)
library(ggthemes)
library(ggpubr)
library(httk) #version 1.8 used
library(openxlsx)
library(plyr)
library(reshape2)
library(stringr)
library(tcpl) #version 2.0.1 used

#getwd()
#setwd("/APCRA_final_retrospective/R") # set directory and all others will be relative

load(file="../expocast/expocast_448.RData") # load ExpoCast data
load(file="../other_data/hc_exposure_total.RData") # load HC exposure data
load(file='../toxcast/source_invitro_data.RData') # load in vitro data
load(file='../toxval/pod_summary_master.RData') # load summary of pod master
load(file='source_invivo_data.RData')
ttc <- as.data.table(read.xlsx('../other_data/TTC/Copy of APCRA_448_v2000_TTC_2018-09-28.xlsx', sheet=1)) #TTC information

#-----------------------------------------------------------------------------------#
# Calculate the human-based administered equivalent doses (AEDs) using httk
#-----------------------------------------------------------------------------------#

load('../httk/new-httk-2018-03-26.RData') # new RData from John Wambaugh 3-26-2018

chem.physical_and_invitro.data <- add_chemtable(new.data2,
                                                current.table=chem.physical_and_invitro.data,
                                                data.list=list(Compound="Compound",
                                                               CAS="CAS",
                                                               DSSTox.GSID="DSSTox_Substance_Id",
                                                               Clint="Human.Clint",
                                                               Clint.pValue="Human.Clint.pValue",
                                                               Funbound.plasma="Human.Funbound.plasma",
                                                               LogP="logP",MW="MW",
                                                               pKa_Accept="pKa_Accept",
                                                               pKa_Donor="pKa_Donor"),
                                                reference="Unpublished Ceetox",
                                                species="Human",
                                                overwrite=T)



#need to first confirm which casrn's have httk information otherwise a loop will fail
invitro.dist.subset<-as.data.frame(subset(invitro.dist, chemcas %in% get_cheminfo(species='Human')))

length(unique(invitro.dist.subset$chemcas)) # 448
colnames(invitro.dist.subset)
aed.df=data.frame()
for (casn in invitro.dist.subset[,'chemcas'])
{
  AED_5<-(calc_mc_oral_equiv(conc=subset(invitro.dist.subset,chemcas==casn)[['AC50p5']], chem.cas=casn, 
                             which.quantile=c(0.05), species='Human', method='dr', 
                             well.stirred.correction=T, restrictive.clearance=T, output.units='mg'))
  AED_50<-(calc_mc_oral_equiv(conc=subset(invitro.dist.subset,chemcas==casn)[['AC50p5']], chem.cas=casn, 
                              which.quantile=c(0.5), species='Human', method='dr', 
                              well.stirred.correction=T, restrictive.clearance=T, output.units='mg'))
  AED_95<-(calc_mc_oral_equiv(conc=subset(invitro.dist.subset,chemcas==casn)[['AC50p5']], chem.cas=casn, 
                              which.quantile=c(0.95), species='Human', method='dr', 
                              well.stirred.correction=T, restrictive.clearance=T, output.units='mg'))
  AED_hipptox_50 <-(calc_mc_oral_equiv(conc=subset(invitro.dist.subset,chemcas==casn)[['min.ec10.um']], chem.cas=casn, 
                                       which.quantile=c(0.5), species='Human', method='dr', 
                                       well.stirred.correction=T, restrictive.clearance=T, output.units='mg'))
  AED_hipptox_95 <-(calc_mc_oral_equiv(conc=subset(invitro.dist.subset,chemcas==casn)[['min.ec10.um']], chem.cas=casn, 
                                    which.quantile=c(0.95), species='Human', method='dr', 
                                    well.stirred.correction=T, restrictive.clearance=T, output.units='mg'))
  max_AED <- (calc_mc_oral_equiv(conc=100, chem.cas=casn, 
                                 which.quantile=c(0.95), species='Human', method='dr', 
                                 well.stirred.correction=T, restrictive.clearance=T, output.units='mg'))
  aed.df<-rbind(aed.df,cbind(casn, AED_5, AED_50, AED_95, AED_hipptox_50, AED_hipptox_95,max_AED))
}

save(aed.df, file='AED_from_httk_18Dec18.RData')

#-----------------------------------------------------------------------------------#
# Build wide data integration
#-----------------------------------------------------------------------------------#
load('AED_from_httk_18Dec18.RData')
# get administered dose equivalents from ToxCast-httk
colnames(aed.df)
pod.ratio.master <- as.data.table(aed.df)

setnames(pod.ratio.master,
            c('casn'),
            c('CASRN'))

# get POD summary information for all sources

colnames(pod.summary)
pod.ratio.master <- merge(pod.ratio.master,
                          pod.summary,
                          by='CASRN')

# get exposure information from ExpoCast

colnames(expocast.448) #expocast seem2 data

pod.ratio.master$expocast.total.med <- expocast.448$Total.median[match(pod.ratio.master$CASRN,
                                                                       expocast.448$CASRN)]
pod.ratio.master$expocast.total.95 <- expocast.448$'Total.95%-ile'[match(pod.ratio.master$CASRN,
                                                                       expocast.448$CASRN)]
pod.ratio.master$Name <- toxcast.master$chnm[match(pod.ratio.master$CASRN,
                                                   toxcast.master$casn)]

# put data into log10 mg/kg-bw/day units

str(pod.ratio.master)
pod.ratio.master <- as.data.frame(pod.ratio.master)

pod.ratio.master[,'AED_5'] <- as.numeric(as.character(pod.ratio.master[,'AED_5']))
pod.ratio.master[,'AED_50'] <- as.numeric(as.character(pod.ratio.master[,'AED_50']))
pod.ratio.master[,'AED_95'] <- as.numeric(as.character(pod.ratio.master[,'AED_95']))
pod.ratio.master[,'AED_hipptox_50'] <- as.numeric(as.character(pod.ratio.master[,'AED_hipptox_50']))
pod.ratio.master[,'AED_hipptox_95'] <- as.numeric(as.character(pod.ratio.master[,'AED_hipptox_95']))
pod.ratio.master[,'max_AED'] <- as.numeric(as.character(pod.ratio.master[,'max_AED']))

pod.ratio.master <- as.data.table(pod.ratio.master)
colnames(pod.ratio.master)
pod.ratio.master[,log.AED.5 := log10(AED_5)]
pod.ratio.master[,log.AED.50 := log10(AED_50)]
pod.ratio.master[,log.AED.95 := log10(AED_95)]
pod.ratio.master[,log.AED.hipptox.50 := log10(AED_hipptox_50)]
pod.ratio.master[,log.AED.hipptox.95 := log10(AED_hipptox_95)]
pod.ratio.master[,log.max.AED := log10(max_AED)]
pod.ratio.master[,log.expo.50 := log10(expocast.total.med)]
pod.ratio.master[,log.expo.95 := log10(expocast.total.95)]
pod.ratio.master[,log.p5.POD := log10(p5.POD)]
pod.ratio.master[,log.p95.POD := log10(p95.POD)]
pod.ratio.master[,log.min.pod := log10(min.POD)]
pod.ratio.master[,log.max.pod := log10(max.POD)]
pod.ratio.master[,log.med.pod := log10(med.POD)]

# define POD-NAM
pod.ratio.master[,pod.nam.95 := min(log.AED.95,log.AED.hipptox.95, na.rm=TRUE), by=list(CASRN)]
pod.ratio.master[,pod.nam.50 := min(log.AED.50,log.AED.hipptox.50, na.rm=TRUE), by=list(CASRN)]

# define POD ratio
pod.ratio.master[,pod.ratio.95 := log.p5.POD - pod.nam.95, by=list(CASRN)]
length(unique(pod.ratio.master[pod.ratio.95<0]$CASRN)) #48 CASRN out of 448 have pod.ratio.95 <0

pod.ratio.master[,pod.ratio.50 := log.p5.POD - pod.nam.50, by=list(CASRN)]
length(unique(pod.ratio.master[pod.ratio.50<0]$CASRN)) #90 CASRN out of 448 have pod.ratio.50 <0

(48/448)*100 #10.7%
(90/448)*100 # 20%

# define bioactivity-exposure ratio (BER)

pod.ratio.master[, ber.95 := pod.nam.95 - log.expo.95]
pod.ratio.master[, ber.50 := pod.nam.95 - log.expo.50]
pod.ratio.master[, ber.95.pod.50 := pod.nam.50 - log.expo.95]
pod.ratio.master[, ber.50.pod.50 := pod.nam.50 - log.expo.50]
count(pod.ratio.master[ber.95<0]$CASRN) #11 CASRN out of 448 have ber <0
count(pod.ratio.master[ber.95.pod.50 < 0]$CASRN) #3 CASRN out of 448 have ber <0
#x freq
#1 219714-96-2    1 #penoxsulam, naphthalene, and coumarin
#2     91-20-3    1
#3     91-64-5    1
pod.ratio.master[ber.95.pod.50<0]

apcra.448 <- as.data.table(read.xlsx('APCRA_Retrospective_List_448.xlsx', sheet=1))

pod.ratio.master$DTXSID <- apcra.448$DTXSID[match(pod.ratio.master$CASRN,
                                                  apcra.448$CASRN)]

colnames(pod.ratio.master)
setcolorder(pod.ratio.master,c(38,1,16,30:37,17:29,2:15))
                               
#-----------------------------------------------------------------------------------#
# Add TTC information
#-----------------------------------------------------------------------------------#

colnames(ttc)
head(ttc)
setnames(ttc,
         c('TTC_ASSIGNED.(µg/.kg.bw/day)'),
         c('Toxtree.TTC.ug.kg.bw.d'))
str(ttc)

ttc[,ttc.ug.kg.d := as.numeric(as.character(Toxtree.TTC.ug.kg.bw.d))]
ttc[is.na(ttc.ug.kg.d)] #12 exclusion and 3 no structure
ttc[,ttc.mkd := `ttc.ug.kg.d`/1000]
ttc[,log.ttc.mkd := log10(ttc.mkd)]

table(ttc$ttc.mkd)
# 2.5e-06   3e-04  0.0015   0.009    0.03 
# 141      36     212       5      39 

pod.ratio.master[ttc, log.ttc.mkd := log.ttc.mkd, on='DTXSID']
pod.ratio.master[ttc, ttc.label := Toxtree.TTC.ug.kg.bw.d, on='DTXSID']
colnames(pod.ratio.master)
nrow(pod.ratio.master[log.ttc.mkd < pod.nam.95])#389 out of 448
389/448 #86.8% of the time, ttc value is less than the pod.nam.95

nrow(pod.ratio.master[log.ttc.mkd < pod.nam.50]) #413 out of 448
(413/448)*100 #92% of the time, the ttc value is less than the pod.nam.50

pod.ratio.master[,ttc.ratio.95 := ifelse(!is.na(log.ttc.mkd), pod.nam.95 - log.ttc.mkd, NA)]
pod.ratio.master[,ttc.ratio.50 := ifelse(!is.na(log.ttc.mkd), pod.nam.50 - log.ttc.mkd, NA)]

table(pod.ratio.master$ttc.label)
str(pod.ratio.master)
pod.ratio.master[ttc.label=='2.5000000000000001E-3', ttc.label :='0.0025']

# see if using TTC instead of POD-NAM would eliminate the BER < 0 priority chemicals

pod.ratio.master[ber.95 <0 & log.ttc.mkd < pod.nam.95] # 2 of the BER <0 chemicals have a TTC that's even lower (coumarin, tribuytl phosphate)
pod.ratio.master[ber.95 <0 & log.ttc.mkd > pod.nam.95] # 8 substances
ttc.solution <- pod.ratio.master[ber.95 <0]
ttc.solution[,ttc.ber.95 := log.ttc.mkd - log.expo.95] #napthalene is the only chemical for which the ttc.ber.95 <0
ttc.solution

#-----------------------------------------------------------------------------------#
# allometrically scale PODs to human equivalent doses as an approach to evaluate interspecies differences
#-----------------------------------------------------------------------------------#

#standardize species first

colnames(pod.master)
count(pod.master$species)

# rat
pod.master[species %like% c('rat'), allo.spec := 'rat']
pod.master[species %like% c('Rat'), allo.spec := 'rat']
count(pod.master[allo.spec=='rat']$species)

# mouse
pod.master[species %like% c('mice'), allo.spec := 'mouse']
pod.master[species %like% c('mouse'), allo.spec := 'mouse']
count(pod.master[allo.spec=='mouse']$species)

#dog
pod.master[species %like% c('dog'), allo.spec := 'dog']
pod.master[species %like% c('Dog'), allo.spec := 'dog']

#guinea pig
pod.master[species %like% c('guinea'), allo.spec := 'guinea pig']

#hamster
pod.master[species %like% c('amster'), allo.spec := 'hamster']

#rabbit
pod.master[species %like% c('abbit'), allo.spec := 'rabbit']

# for these species, allometrically scale per Nair AB, Jacob S. A simple practice guide for dose conversion between
# animals and human. Journal of basic and clinical pharmacy 2016, 7: 27-31.

pod.master[ allo.spec == "mouse", HED := toxval_numeric_mkd_only * 0.081]
pod.master[ allo.spec == "rat", HED := toxval_numeric_mkd_only * 0.162]
pod.master[ allo.spec == "guinea pig", HED := toxval_numeric_mkd_only * 0.216]
pod.master[ allo.spec == "rabbit", HED := toxval_numeric_mkd_only * 0.324]
pod.master[ allo.spec == "dog", HED := toxval_numeric_mkd_only * 0.541]
pod.master[ allo.spec == "hamster", HED := toxval_numeric_mkd_only * 0.135]
pod.master.hed <- pod.master[!is.na(HED)]

#-----------------------------------------------------------------------------------#
# create HED summary
#-----------------------------------------------------------------------------------#

pod.summary.hed <- ddply(pod.master.hed, .(CASRN), summarize, 
                         min.POD=min(HED),
                         mean.POD=mean(HED),
                         max.POD=max(HED),
                         med.POD=quantile(HED,0.5,type=2,na.rm=TRUE),
                         p5.POD=quantile(HED,0.05,type=2,na.rm=TRUE),
                         p95.POD=quantile(HED,0.95,type=2,na.rm=TRUE))

#-----------------------------------------------------------------------------------#
# amend data integration
#-----------------------------------------------------------------------------------#

pod.summary.hed <-as.data.table(pod.summary.hed)
pod.summary.hed[,log.hed.5 := log10(p5.POD)]
pod.summary.hed[,log.hed.95 := log10(p95.POD)]

pod.ratio.master$log.hed.5 <- pod.summary.hed$log.hed.5[match(pod.ratio.master$CASRN,
                                                              pod.summary.hed$CASRN)]

pod.ratio.master[,hed.pod.ratio.95 := log.hed.5 - pod.nam.95, by=list(CASRN)]
length(unique(pod.ratio.master[hed.pod.ratio.95<0]$CASRN)) #82
82/447 #18%

pod.ratio.master[,hed.pod.ratio.50 := log.hed.5 - pod.nam.50, by=list(CASRN)]
length(unique(pod.ratio.master[hed.pod.ratio.50<0]$CASRN)) #158
158/447 #35%

#-----------------------------------------------------------------------------------#
# add binary pod ratio field
#-----------------------------------------------------------------------------------#

pod.ratio.master[pod.ratio.95 < 0, pod.ratio.95.binary := 0]
pod.ratio.master[pod.ratio.95 > 0, pod.ratio.95.binary := 1]
table(pod.ratio.master$pod.ratio.95.binary)
#0   1 
#48 400 

pod.ratio.master[pod.ratio.50 < 0, pod.ratio.50.binary := 0]
pod.ratio.master[pod.ratio.50 > 0, pod.ratio.50.binary := 1]
table(pod.ratio.master$pod.ratio.95.binary)
#0   1 

colnames(pod.ratio.master)

#-----------------------------------------------------------------------------------#
# save final files
#-----------------------------------------------------------------------------------#

save(pod.ra.uniq,
     pod.summary.hed,
     pod.master,
     file='source_invivo_data.RData')


colnames(pod.ratio.master)
setcolorder(pod.ratio.master,
            c(1:40,46:47,41:45))
write.xlsx(pod.ratio.master, file='../toxval/pod_ratio_master_final.xlsx')
save(pod.ratio.master, file='../toxval/pod_ratio_master_final.RData')


