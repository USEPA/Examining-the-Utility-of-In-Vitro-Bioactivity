#-----------------------------------------------------------------------------------#
# This script brings in other in vivo and exposure data sources from EPA partners: 
# Health Canada, EFSA, ECHA, etc.to make the POD master

# Katie Paul Friedman, paul-friedman.katie@epa.gov

# Updated 6 Dec 2018
#-----------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#
# Loading libraries
#-----------------------------------------------------------------------------------#
rm(list = ls())

library(cowplot)
library(data.table)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(httk) # version 1.8 used
library(plyr)
library(RColorBrewer)
library(reshape2)
library(scales)
library(stringr)
library(tcpl) # version 2.0.1 used
library(openxlsx)

tcplConf(user='', pass='', db='invitrodb_v3', drvr='MySQL', host='')
setwd()# to folder APCRA_final_retrospective/R

#-----------------------------------------------------------------------------------#
# Combining the POD data from all sources into POD.master
#-----------------------------------------------------------------------------------#
load('../in_vivo/source/toxval_master.RData')
load('../in_vivo/source/echa_pod_master.RData')
load('../in_vivo/source/efsa_pod_master.RData')
load('../in_vivo/source/hc_POD_master.RData')
load('../in_vitro/source_invitro_data.RData')

pod.master<-rbind(toxval.master,
                  hc.pod.master,
                  efsa.pod.master,
                  echa.pod.master)

table(pod.master$risk_assessment_class)

pod.master[grepl('acute', risk_assessment_class), risk_assessment_class := 'acute']
pod.master[risk_assessment_class %in% c('growth:chronic',
                                        'mortality:chronic'), risk_assessment_class := 'chronic']
pod.master[grepl('develop', risk_assessment_class), risk_assessment_class := 'developmental / reproductive']
pod.master[grepl('repro', risk_assessment_class), risk_assessment_class := 'developmental / reproductive']

pod.master[risk_assessment_class %in% c('study with volunteers',
                                        'human'), risk_assessment_class := 'human']
pod.master[risk_assessment_class %in% c('short-term', 'other', 'range-finding'), risk_assessment_class := 'repeat dose']

pod.master <- pod.master[!(risk_assessment_class=='ecotox')] # decided to eliminate all terrestrial animal studies for ecotox in consultation with Rusty
pod.master <- pod.master[!(source=='ECOTOX')]
pod.master[risk_assessment_class %in% c('toxicokinetic study'), risk_assessment_class := 'developmental / reproductive']
pod.master <- pod.master[!(CASRN=='1162-65-8')] # aflatoxin B1 not in toxcast

#-----------------------------------------------------------------------------------#
# Filtering the POD master for httk availability
#-----------------------------------------------------------------------------------#

length(unique(pod.master$CASRN)) #510

load('../httk/new-httk-2018-03-26.RData')

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

colnames(chem.physical_and_invitro.data)

ade.get.pod <-as.data.frame(subset(pod.master, CASRN %in% get_cheminfo(species='Human')))
length(unique(ade.get.pod$CASRN)) #448 CASRN

# check to see addition of CASRN by source
table(pod.master$source)

length(setdiff(pod.master$CASRN, toxval.master$CASRN)) #68 CASRN are in pod.master and not toxval.master (added by outside partners)
length(setdiff(pod.master[source=='Health.Canada']$CASRN, toxval.master$CASRN)) #2 CASRN from Health.Canada
length(setdiff(pod.master[source=='ECHA_Rasenberg']$CASRN, toxval.master$CASRN)) #51 CASRN from ECHA
length(setdiff(pod.master[source=='EFSA.httk.common']$CASRN, toxval.master$CASRN)) #0

length(intersect(invitro.dist$chemcas, pod.master$CASRN)) #448 CASRN unique, confirmed as also in httk +Toxcast
pod.master.casn <- intersect(invitro.dist$chemcas, pod.master$CASRN)

pod.master <- pod.master[CASRN %in% pod.master.casn]
length(unique(pod.master$CASRN)) #confirmed 448

save(pod.master, file='../in_vivo/pod_master.RData')

#---------------------------------------------------------------------#
# Summarize POD information by substance
#---------------------------------------------------------------------#
load('../toxval/pod_master.RData')

pod.summary <- ddply(pod.master, .(CASRN), summarize, 
                     min.POD=min(toxval_numeric_mkd_only),
                     mean.POD=mean(toxval_numeric_mkd_only),
                     med.POD=quantile(toxval_numeric_mkd_only,0.5,type=2,na.rm=TRUE),
                     max.POD=max(toxval_numeric_mkd_only),
                     p5.POD=quantile(toxval_numeric_mkd_only,0.05,type=2,na.rm=TRUE),
                     p95.POD=quantile(toxval_numeric_mkd_only,0.95,type=2, na.rm=TRUE))

save(pod.summary, file='../in_vivo/pod_summary_master.RData')

#---------------------------------------------------------------------#
# Prepare in vivo data for hypergeometric enrichment tests
# to determine if certain study types are enriched when pod.ratio <0
#---------------------------------------------------------------------#
load('../in_vivo/pod_master.RData')
colnames(pod.master)

count(pod.master$risk_assessment_class)
#x  freq
#1                        acute    66
#2                      chronic 12464
#3 developmental / reproductive  2548
#4               immunotoxicity     1
#5             immunotoxicology    89
#6                neurotoxicity    80
#7                  repeat dose   487
#8                   subchronic  6892

pod.master[toxval_numeric_mkd_only==0, toxval_numeric_mkd_only:=NA]
pod.master[, min.POD := min(toxval_numeric_mkd_only, na.rm=TRUE), by=list(CASRN)]
pod.master[,log10.toxval.numeric.mkd.only := log10(toxval_numeric_mkd_only)]
pod.master[, min.log10.POD := min(log10.toxval.numeric.mkd.only, na.rm=TRUE), by=list(CASRN)]
pod.master[toxval_numeric_mkd_only == min.POD, min.pod.ra.class := risk_assessment_class, by=list(CASRN,min.POD)]

pod.ra <- pod.master[!is.na(min.pod.ra.class)]
length(unique(pod.ra$CASRN)) #448

colnames(pod.ra)
pod.ra.uniq <- unique(pod.ra[,c("CASRN",
                                "Name",
                                "toxval_numeric",
                                "risk_assessment_class",
                                "toxval_numeric_mkd_only",
                                "min.POD",
                                "log10.toxval.numeric.mkd.only",
                                "min.log10.POD",
                                "min.pod.ra.class")]) #gives a table of 530 obs, suggesting more than one row for each chem

# have to decide how to handle possible disagreements
# mostly the repetition is due to multiple values from the same study and ra class
count(pod.ra.uniq$min.pod.ra.class)
#x freq
#1                        acute   11
#2                      chronic  328
#3 developmental / reproductive   46
#4             immunotoxicology    2
#5                neurotoxicity    1
#6                  repeat dose   36
#7                   subchronic  106

pod.ra[min.pod.ra.class=='human' & Name=='Triclosan'] #Allmyr et al 2009 subchronic duration class
pod.ra.uniq[min.pod.ra.class %in% c('repeat dose','subchronic'), min.pod.ra.class := 'repeat']
pod.ra.uniq[min.pod.ra.class=='human'& Name=='Triclosan', min.pod.ra.class := 'repeat']

casrns.dup <- pod.ra.uniq[duplicated(CASRN)]$CASRN
length(unique(casrns.dup)) #69
pod.ra.uniq.dup <- pod.ra.uniq[CASRN %in% casrns.dup]
setDT(pod.ra.uniq.dup)
RA.types <- pod.ra.uniq.dup[,.(min.pod.ra.class = paste(min.pod.ra.class, collapse=",")), by=CASRN]
pod.ra.uniq.dup$RA.types <- RA.types$min.pod.ra.class[match(pod.ra.uniq.dup$CASRN,
                                                            RA.types$CASRN)]
count(pod.ra.uniq.dup$RA.types)

# order of replacement is important to give preference in this order: 
# 1-neurotoxicity, 2-developmental / reproductive, 3-acute, 4-repeat, 5-chronic
# decided to drop immunotoxicity because there was only 2 studies

pod.ra.uniq.dup[grepl('chronic', RA.types), min.pod.ra.class := 'chronic', by = CASRN]
pod.ra.uniq.dup[grepl('repeat', RA.types), min.pod.ra.class := 'repeat', by = CASRN]
pod.ra.uniq.dup[grepl('acute', RA.types), min.pod.ra.class := 'acute', by = CASRN]
pod.ra.uniq.dup[grepl('develop', RA.types), min.pod.ra.class := 'developmental / reproductive', by = CASRN]
pod.ra.uniq.dup[grepl('neurotoxicity', RA.types), min.pod.ra.class := 'neurotoxicity', by = CASRN]

pod.ra.uniq.dup[,c('RA.types'):= NULL]
colnames(pod.ra.uniq.dup)

pod.ra.uniq.dup <- unique(pod.ra.uniq.dup, by=c('CASRN', 'min.pod.ra.class'))
pod.ra.uniq <- pod.ra.uniq[!CASRN %in% casrns.dup]
pod.ra.uniq <- rbind(pod.ra.uniq,pod.ra.uniq.dup)

length(unique(pod.ra.uniq$CASRN))

count(pod.ra.uniq$min.pod.ra.class)

#x freq
#1                        acute   11
#2                      chronic  272
#3 developmental / reproductive   44
#4                neurotoxicity    1
#5                       repeat  120

#11+272+44+1+120 = 448
