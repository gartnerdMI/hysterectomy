############################################
# Project: Hyst rates                      #
# Date(s): 9/24/18-9/13/19                 #
# Programmer: D.R. Gartner                 #
#                                          #
# Input: inpt_2014.sas7bdat                #
#        outpt_2014.sas7bdat               #
#        inpt_2013.sas7bdat                #
#        outpt_2013.sas7bdat               #
#        inpt_2012.sas7bdat                #
#        outpt_2012.sas7bdat               #
#        inpt_2011.sas7bdat                #
#        outpt_2011.sas7bdat               #
# Years: 2011-2014                         #  
#                                          #  
# Associated Files:                        #  
# 1) EpidemiologyPaperCode.R (R Script)    #
# 2) nc.1990_2016.singleages.txt (census   #
#    estimates from SEER downloaded here:  #
# https://seer.cancer.gov/popdata/download.html)  #
# 3) NCrace.pooled.prev.csv                #
#     (BRFSS hyst estimates using combo of #
#     2010, 2012, 2014, est), by race/eth  #
# 4) NC.prev.overall.csv (hyst prev est)   #
# 5) county_helper.csv (county names)      #
# 6) NCrace.pooled.prev_95CI.csv (pooled   #
#     BRFSS hyst est w/ 95% CI)            #
#                                          #  
############################################
# Description: calculate race/ethnicity specific hysterectomy rates for NC residents
# w/ non-cancerous gyn conditions

# Analysis Steps:
# Pt 1: Load packages and ready data
# Pt 2: Apply exclusion criteria
# Pt 3: Recode procedure variables & apply further exclusions
# Pt 4: Load & prep census data
# Pt 5: Calculation of hysterectomy incidence rates by race/ethnicity
# Pt 6: Bootstrap rate 95% CI
# Pt 7: Bootstrap rate difference 95% CI
# Pt 8: CCI calculation
# Pt 9: Sensitivity Analyses
# Prologue: BRFSS Prevalence calculations


# NOTES: Data include all hysterectomy and oophorectomy occurrences and limiting to
# cases of hysterectomy is needed (i.e., exclude cases of ooph if there was no
# concommitant hyst)

#....................................................................
# Pt 1: Load packages, read hyst count files, make smaller datasets, merge ####
#....................................................................

library(tidyverse) #for data management
library(Hmisc) #for quick descriptives
library(sas7bdat) #for reading SAS data files
library(survey) #for 
library(srvyr)
library(icd) #for Charlson Comorbidity calculations
library(boot) #to bootstrap estimates
library(knitr)

setwd("savefile") #for saving smaller datasets, output, etc.
dir() #view contents

inpat14 <- read.sas7bdat("savefile/inpt_2014.sas7bdat")
outpat14 <- read.sas7bdat("savefile/outpt_2014.sas7bdat")
inpat13 <- read.sas7bdat("savefile/inpt_2013.sas7bdat")
outpat13 <- read.sas7bdat("savefile/outpt_2013.sas7bdat")
inpat12 <- read.sas7bdat("savefile/inpt_2012.sas7bdat")
outpat12 <- read.sas7bdat("savefile/outpt_2012.sas7bdat")
inpat11 <- read.sas7bdat("savefile/inpt_2011.sas7bdat")
outpat11 <- read.sas7bdat("savefile/outpt_2011.sas7bdat")

#view inpatient variable names by year
names(inpat14)
names(inpat13)
names(inpat12)
names(inpat11)

# to keep (inpat): "shepsid", "fyear", "ptstate/patst", "ptcnty", "ptzip", "agey", "sex", "race", "admitdx", "ethnicity", 
#"diag1-diag25", "proccd1-proccd20" requires different subset for each data set (thanks a lot sheps)
inpat14 <- inpat14[c(2:9, 12, 13, 16:40, 66:85, 101)] #obs=6660, vars=56
inpat13 <- inpat13[c(1:8, 11, 12, 15:39, 65:84, 91)] #obs=6871, vars=56
inpat12 <- inpat12[c(1:8, 11, 12, 15:39, 65:84, 91)] #obs=7730, vars=56
inpat11 <- inpat11[c(1:8, 11, 12, 15:59, 66)] #obs=9482, vars=56

#view outpatient variable names by year
names(outpat14)
names(outpat13)
names(outpat12)
names(outpat11)

# to keep (outpat): "fyear", "ptstate", "ptcnty", "ptzip", "agey", "sex", "race", "type", "admitdx", "ethnicity", 
#"diag1-diag25", "cpxcd1-cpxcd20", "shepsid"
outpat14 <- outpat14[c(2:9, 11, 12, 16:60, 14)] #obs=15122, vars=56
outpat13 <- outpat13[c(1:8, 10, 11, 15:59, 13)] #obs=14658, vars=56
outpat12 <- outpat12[c(1:8, 10, 11, 15:59, 13)] #obs=14952, vars=56
outpat11 <- outpat11[c(1:8, 11, 12, 15:39, 60:79, 86)] #obs=14856, vars=56

# add inpatient and outpatient indicator to each
inpat14$setting <- 1
outpat14$setting <- 2
inpat13$setting <- 1
outpat13$setting <- 2
inpat12$setting <- 1
outpat12$setting <- 2
inpat11$setting <- 1
outpat11$setting <- 2

#transform shepsid into character
inpat14$shepsid <- as.character(inpat14$shepsid)
inpat13$shepsid <- as.character(inpat13$shepsid)
inpat12$shepsid <- as.character(inpat12$shepsid)
inpat11$shepsid <- as.character(inpat11$shepsid)
outpat14$shepsid <- as.character(outpat14$shepsid)
outpat13$shepsid <- as.character(outpat13$shepsid)
outpat12$shepsid <- as.character(outpat12$shepsid)
outpat11$shepsid <- as.character(outpat11$shepsid)

#rename variables that mismatch (patst into ptstate)
inpat12 <- rename(inpat12, ptstate=patst)
inpat11 <- rename(inpat11, ptstate=patst)
outpat12 <- rename(outpat12, ptstate=patst)
outpat11 <- rename(outpat11, ptstate=patst)

#merge inpt files together
inpat1 <- dplyr::union_all(inpat14,inpat13) #obs=13531, vars=57
inpat2 <- dplyr::union_all(inpat1, inpat12) #obs=21261, vars=57
fullinpat <- dplyr::union_all(inpat2, inpat11) #obs=30743, vars=57

#merge outpt files together
outpat1 <- dplyr::union_all(outpat14,outpat13) #obs=29780, vars=57
outpat2 <- dplyr::union_all(outpat1,outpat12) #obs=44732, vars=57
fulloutpat <- dplyr::union_all(outpat2, outpat11) #obs=59588, vars=57

#merge in/out together
full <- dplyr::union_all(fullinpat, fulloutpat) #obs=90331, vars=77

#create factor indicator for setting
full$setting_f <- factor(full$setting, levels=1:2,
                         labels = c("inpatient", "outpatient")) #vars 78

#.......................................................
# NOTES: all variable types were coerced into character
# type when the in/out datasets were unioned.

#.......................................................
# Pt 2: Apply exclusion criteria to hyst procedure data ####
#.......................................................

#Restrict to females, NC residents, ages 18-44, non-trauma
# deleted 10 males
# deleted 4033 non-NC residents
# deleted 43114 ages >44
# deleted 348 ages < 18
# deleted 8 truamas

full <- full[ which(full$sex == "F"),] #obs=90321
fullNC <- full[ which(full$ptstate == "NC"),] #obs=86288
fullpremen <- fullNC[ which(fullNC$agey < 45),] #obs=43174 
fullpremen2 <- fullpremen[ which(fullpremen$ agey > 17),] #obs=42826
fullnotrauma <- fullpremen2[ which(fullpremen2$type != "5"),] #obs=42818

#Restrict to benign conditions
#excluded 1539 w/ cancer dx

#dataframe w/ cancers included = fullnotrauma (analysis script elsewhere)
#dataframe w/ cancers excluded = fullnocancer

#create function to see if cancer codes exist in diag1-25 or admitdx
cancer <- function(x) {
  x %in% c('1746', '1741', '1743', '1745', '174', '1748', '1749', '1742', '1744',
           '1809', '1843', '182', '180', '1801', '1832', '182', '1821', '1841', '1842',
           '182','183','1808','1953','184','1844','1986','19882','1986','181',
           '2361', '2362','17351','17359','1735','1588','1589','158','1588',
           '17352','1820','1830','179','1976','1541','1800','1533','V5049',
           '1889','1540','2395') 
}

#apply function
vars_of_interest <- c(paste("diag", 1:25, sep=""), "admitdx")
cancerlist <- as.data.frame(lapply(fullnotrauma[,vars_of_interest], FUN = cancer))
fullnotrauma$cancer <- apply(cancerlist, MARGIN = 1, function(x) {
  ifelse(all(x==0),0, 1)
}
)
table(fullnotrauma$cancer, useNA = "always")

fullnocancer <- fullnotrauma[ which(fullnotrauma$cancer == 0),] #obs=41279, var=79

# restrict to hyst (w/ or w/o ooph) = delete ooph only (regardless of uni or bilateral
# or unclear if uni or bi removal) - Inpatient first
# excluded 2145 ooph only cases
# double checked codes w/ gyn provider on 4/19/19 (see emails for details)

ooph_icd9 <- function(x) {
  x %in% c('655','6551','6552','6553','6554','656','6561','6562','6563','6564'
           ,'653','6531','6539','654','6541','6549')
}

anyhyst_icd9 <- function(x) {
  x %in% c('68','683','6831','6839','684','6841','6849','685','6851','6859',
           '686','6861','6869','687','6871','6879','689')
}

#apply functions
vars_of_interest2 <- c(paste("proccd", 1:20, sep=""))

oophlist <- as.data.frame(lapply(fullnocancer[,vars_of_interest2], FUN = ooph_icd9))
hystlist <- as.data.frame(lapply(fullnocancer[,vars_of_interest2], FUN = anyhyst_icd9))

#create indicator for all oophs
fullnocancer$ooph_icd_ind <- apply(oophlist, MARGIN = 1, function(x) {
  ifelse(all(x==0),0, 1)
}
)
table(fullnocancer$ooph_icd_ind, useNA = "always")

#create indicator for all hysts
fullnocancer$hyst_icd_ind <- apply(hystlist, MARGIN = 1, function(x) {
  ifelse(all(x==0),0, 1)
}
)
table(fullnocancer$hyst_icd_ind, useNA = "always")

# restrict to hyst (w/ or w/o ooph) = delete ooph only (regardless of uni or bilateral
# or unclear if uni or bi removal) - Outpatient
# excluded 7874 ooph only cases
# double checked codes w/ gyn provider on 4/19/19 (see emails for details)

ooph_cpt <- function(x){
  x %in% c('58291','58542','58544','58548','58552','58554','58571','58573',
           '58951','58953','58954','58956','58700','58720','58925','58940',
           '58943','58950','58951','58952','58661')
}

anyhyst_cpt <- function(x){
  x %in% c('58291','58542','58544','58548','58552','58554','58571','58573',
           '56308','58150','58152',
           '58180','58200','58210','58260','58285','58290','58541','58543',
           '58550','58553','58570','58572','58262','58270','58294','58267',
           '58263','58275','58280','58292','58293'
           #,'58951','58953','58954','58956'cancer specific codes
           )
}

#apply functions
vars_of_interest3 <- c(paste("cpxcd", 1:20, sep=""))

cpt_oophlist <- as.data.frame(lapply(fullnocancer[,vars_of_interest3], FUN = ooph_cpt))
cpt_hystlist <- as.data.frame(lapply(fullnocancer[,vars_of_interest3], FUN = anyhyst_cpt))

#create indicator for all oophs
fullnocancer$ooph_cpt_ind <- apply(cpt_oophlist, MARGIN = 1, function(x) {
  ifelse(all(x==0),0, 1)
}
)
table(fullnocancer$ooph_cpt_ind, useNA = "always")

#create indicator for all hysts
fullnocancer$hyst_cpt_ind <- apply(cpt_hystlist, MARGIN = 1, function(x) {
  ifelse(all(x==0),0, 1)
}
)
table(fullnocancer$hyst_cpt_ind, useNA = "always")

#exclude ooph only cases
fullnocancer$oophonly <- ifelse(((fullnocancer$ooph_icd_ind==1 & fullnocancer$hyst_icd_ind==0) |
                                   (fullnocancer$ooph_cpt_ind==1 & fullnocancer$hyst_cpt_ind==0)), 1, 0)
table(fullnocancer$oophonly, useNA = "always") #10019 w/ ooph only

oophcheck <- ifelse((fullnocancer$ooph_icd_ind==1 & fullnocancer$hyst_icd_ind==0), 1,0)
table(oophcheck) #2145 had ooph only (inpatient)
oophcheck2 <- ifelse((fullnocancer$ooph_cpt_ind==1 & fullnocancer$hyst_cpt_ind==0), 1,0)
table(oophcheck2) #7874 had ooph only (oupatient)

#create indicator for those 2011 procedures that aren't hyst or ooph
fullnocancer$nohystnoooph_ind <- ifelse((fullnocancer$fyear == "2011" &
                                           fullnocancer$setting_f == "outpatient" &
                                           fullnocancer$hyst_cpt_ind == 0 &
                                           fullnocancer$ooph_cpt_ind == 0), 1,0)
table(fullnocancer$nohystnoooph_ind)

fullfinal <- fullnocancer[ which(fullnocancer$oophonly == 0 & 
                                   fullnocancer$nohystnoooph_ind == 0),] #obs=31235, var=85
#.......................................................
# NOTES: 3.4% w/ a cancer dx...does that seem low? Checked with
# old output from 2011-2013 and it looks right on par.

# One surgery (shepsid = 1013700601515304657) where there are no
# other info: 32 years, NH Black, total abdominal hyst, 2011, outpatient
# included in numerator estimates

# 25 with no hyst and no ooph codes; have been excluded.

# Exclusion totals:
# males = 10
# non-NC residents = 4033
# <18 or >44 = 43462
# traumas = 8
# cancers = 1539
# ooph only (no hyst) = 10019
# no ooph & no hyst = 25
# Final N = 31235

#.......................................................
# Pt 3: Recode hyst procedure variables & apply further exclusions ####
#.......................................................

#recode ethnicity variable (to "hisp" and make factor "hisp_f")
describe(fullfinal$ethnicity)

fullfinal$hisp <- ifelse(fullfinal$ethnicity == 2, 1, 0) #recode missing as non-Hisp
table(fullfinal$hisp, useNA = "always")

fullfinal$hisp_f <- factor(fullfinal$hisp, levels = 0:1,
                           labels = c("Non-Hispanic", "Hispanic"))
table(fullfinal$hisp_f, useNA = "always")

#recode race for
describe(fullfinal$race)

fullfinal$raceeth[fullfinal$hisp==1] <- 3 #Hispanic
fullfinal$raceeth[fullfinal$race==1 & fullfinal$hisp==0] <- 4 #NH AmInd
fullfinal$raceeth[fullfinal$race==2 & fullfinal$hisp==0] <- 5 #NH Asian
fullfinal$raceeth[fullfinal$race==3 & fullfinal$hisp==0] <- 2 #NH Black
fullfinal$raceeth[fullfinal$race==4 & fullfinal$hisp==0] <- 5 #NH Asian
fullfinal$raceeth[fullfinal$race==5 & fullfinal$hisp==0] <- 1 #NH White
fullfinal$raceeth[fullfinal$race==6 & fullfinal$hisp==0] <- 6 # NH Other
fullfinal$raceeth[fullfinal$race==9 & fullfinal$hisp==0] <- NA #Missing

table(fullfinal$raceeth, useNA = "always")

fullfinal$raceeth_f <- factor(fullfinal$raceeth, levels = 1:6,
                              labels = c("NH White", "NH Black", "Hispanic",
                                         "NH American Indian", "NH Asian/PI", "NH Other"))
table(fullfinal$raceeth_f, useNA = "always")

#subset to non-NH Other (to match census denominators)
fullfinal <- fullfinal[ which(fullfinal$raceeth_f != "NH Other" | fullfinal$raceeth != NA),]
#recode raceeth_f variable to to have 5 levels
fullfinal$raceeth_f <- factor(fullfinal$raceeth, levels = 1:5,
                              labels = c("NH White", "NH Black", "Hispanic",
                                         "NH American Indian", "NH Asian/PI"))
table(fullfinal$raceeth_f, useNA = "always")

#Excluded also: 
# 355 w/ missing Race
# 451 that are NH Other
# Analytic N = 30429
#.......................................................
#NOTE: be careful about 2013 - codebook says 9=other race; all other years say 6=other race
# Subset to all non-NH Other due to denominators not existing for NH Other
# Analytic N = 30429

#.......................................................
# Pt 4: Load & Prep Census Data ####
#.......................................................

census <- read.table(".savefile/nc.1990_2016.singleages.txt", stringsAsFactors = FALSE)

census2 <- census %>%
  separate(V1, into = c("Year", "State", "St fips", "Cnty fips", "Registry", 
                        "Race", "Origin", "Sex", "Age", "Population"), 
           sep = c(4,6,8,11,13,14,15,16,18))

#Restrict to 2011-2014, Females, ages 18-44
census2$Agenum <- as.numeric(census2$Age) #recode to be numeric
census2$Popsize <- as.numeric(census2$Population) #recode to be numeric

census3 <- census2[ which((census2$Year == "2014" | census2$Year == "2013" | 
                             census2$Year == "2012" | census2$Year == "2011") & 
                            census2$Sex =="2" & census2$Agenum >17 & 
                            census2$Agenum < 45),]

#Code race/eth variable and create factor variable (raceeth and raceeth_f)
describe(census3$Race)
describe(census3$Origin)
census3$raceeth[census3$Race == 1 & census3$Origin == 0] <- 1
census3$raceeth[census3$Race == 2 & census3$Origin == 0] <- 2
census3$raceeth[census3$Race == 3 & census3$Origin == 0] <- 4
census3$raceeth[census3$Race == 4 & census3$Origin == 0] <- 5
census3$raceeth[census3$Origin == 1] <- 3
table(census3$raceeth)

census3$raceeth_f <- factor(census3$raceeth, levels = 1:5,
                            labels = c("NH White", "NH Black", "Hispanic",
                                       "NH American Indian", "NH Asian/PI"))
table(census3$raceeth_f, useNA = "always") #check

#Code 5 year age groups: 18-19, 20-24, 25-29, 30-34, 35-39, 40-44 (agegrp)
#and factor version (agegrp_f)
summary(census3$Agenum)
census3$agegrp[census3$Agenum < 20] <- 1 #18-19, group 1
census3$agegrp[census3$Agenum > 19 & census3$Agenum < 25] <- 2 #20-24, grp 2
census3$agegrp[census3$Agenum > 24 & census3$Agenum < 30] <- 3 #25-29, grp 3
census3$agegrp[census3$Agenum > 29 & census3$Agenum < 35] <- 4 #30-34, grp 4
census3$agegrp[census3$Agenum > 34 & census3$Agenum < 40] <- 5 #35-39, grp 5
census3$agegrp[census3$Agenum > 39 & census3$Agenum < 45] <- 6 #40-44, grp 6
describe(census3$agegrp)

census3$agegrp_f <- factor(census3$agegrp, levels = 1:6,
                           labels = c("18-19 yrs", "20-24 yrs", "25-29 yrs", "30-34 yrs",
                                      "35-39 yrs", "40-44 yrs"))
table(census3$agegrp_f, useNA = "always") #check

#merge county names data
names <- read.csv(".savefile/county_helper.csv",
                  header = TRUE, stringsAsFactors = FALSE)
names <- names[c(1:100),] #delete last observation (out of state)
names <- names %>%
  separate(FIPS, into=c("State", "cnty"), sep=2) %>%
  select(-State)

census3$cnty <- census3$`Cnty fips`
census3 <- full_join(census3, names, by=c("cnty"))

#.......................................................
# NOTES: Deleted the NH other for the hyst data


#.......................................................
# Pt 5: Calc race/eth specific rates - avg over 4 years ####
#.......................................................

#CRUDE RATES
#create numerator for each race/eth by year, then average over the 4 years
# to create the average number of hysts by race eth during the 4 year period
nums <- fullfinal %>%
  group_by(raceeth_f, fyear) %>%
  dplyr::summarise(total=n())

combo.nums <- nums %>%
  group_by(raceeth_f) %>%
  dplyr::summarise(all=sum(total))

combo.nums$avgtotal <- combo.nums$all/4
combo.nums <- rename(combo.nums, total4yr=all)

#create denominator for each race/eth by year, then average over the 4 years
# to create the average population by race eth during the 4 year period
dens <- census3 %>%
  group_by(raceeth_f, Year) %>%
  dplyr::summarise(Popsize=sum(Popsize))

combo.dens <- dens %>%
  group_by(raceeth_f) %>%
  dplyr::summarise(all=sum(Popsize))

combo.dens$avgPopsize <- combo.dens$all/4
combo.dens <- rename(combo.dens, Popsize4yr=all)

#merge by raceeth and calc rates
crude <- left_join(combo.dens, combo.nums, by=c("raceeth_f"="raceeth_f"))

#calculate rates
crude$rate <- crude$avgtotal/crude$avgPopsize
crude$avgrateby10K <- crude$rate*10000
crude$avgrateby10K

#AGE ADJUSTED RATES
#ready denominator data for age adjustment: grab 2010 Census numbers
forageadj <- census2[ which(census2$Year == "2010" & census2$Sex =="2" & census2$Agenum >17 & 
                              census2$Agenum < 45),]

#Code race/eth variable (raceeth and raceeth_f)
forageadj$raceeth[forageadj$Race == 1 & forageadj$Origin == 0] <- 1
forageadj$raceeth[forageadj$Race == 2 & forageadj$Origin == 0] <- 2
forageadj$raceeth[forageadj$Race == 3 & forageadj$Origin == 0] <- 4
forageadj$raceeth[forageadj$Race == 4 & forageadj$Origin == 0] <- 5
forageadj$raceeth[forageadj$Origin == 1] <- 3
table(forageadj$raceeth)

forageadj$raceeth_f <- factor(forageadj$raceeth, levels = 1:5,
                              labels = c("NH White", "NH Black", "Hispanic",
                                         "NH American Indian", "NH Asian/PI"))
table(forageadj$raceeth_f, useNA = "always") #check

#Code 5 year age groups: 18-19, 20-24, 25-29, 30-34, 35-39, 40-44
summary(forageadj$Agenum)
forageadj$agegrp[forageadj$Agenum < 20] <- 1 #18-19, group 1
forageadj$agegrp[forageadj$Agenum > 19 & forageadj$Agenum < 25] <- 2 #20-24, grp 2
forageadj$agegrp[forageadj$Agenum > 24 & forageadj$Agenum < 30] <- 3 #25-29, grp 3
forageadj$agegrp[forageadj$Agenum > 29 & forageadj$Agenum < 35] <- 4 #30-34, grp 4
forageadj$agegrp[forageadj$Agenum > 34 & forageadj$Agenum < 40] <- 5 #35-39, grp 5
forageadj$agegrp[forageadj$Agenum > 39 & forageadj$Agenum < 45] <- 6 #40-44, grp 6
describe(forageadj$agegrp)

forageadj$agegrp_f <- factor(forageadj$agegrp, levels = 1:6,
                             labels = c("18-19 yrs", "20-24 yrs", "25-29 yrs", "30-34 yrs",
                                        "35-39 yrs", "40-44 yrs"))
table(forageadj$agegrp_f, useNA = "always") #check

#create counts by age and race of population in 2010 (1 year age groups)
agerace1 <- forageadj %>%
  group_by(raceeth_f, Agenum) %>%
  dplyr::summarise(Popsize=sum(Popsize)) #counts stratified by age and race (1 year groups)

racetotal1 <- agerace1 %>%
  group_by(raceeth_f) %>%
  dplyr::summarise(TotalPopRace=sum(Popsize)) #counts stratified by race only (1 year groups)

agetotal1 <- forageadj %>%
  group_by(Agenum) %>%
  dplyr::summarise(PopsizeAge=sum(Popsize)) #counts stratified by age only (1 year groups)

racetotal1 %>% dplyr::summarise(TotalTotal=sum(TotalPopRace)) 
#1765846 = pop size of females in NC between 18-44 years of age

#merge all estimates for calculating proportions (weights)
ageadjprop1 <- left_join(agerace1, racetotal1, by="raceeth_f") 
ageadjprop1 <- left_join(ageadjprop1, agetotal1, by="Agenum")
ageadjprop1$TotalTotal <- 1765846

#calculate proportions (weights) for standardizing
ageadjprop1$propxage <- ageadjprop1$PopsizeAge/ageadjprop1$TotalTotal #compare across race groups,
#standardized to NC age distribution in 2010

#rename variables for ease of later merging
ageadjprop1 <- rename(ageadjprop1, agerace2010=Popsize, race2010=TotalPopRace, age2010=PopsizeAge,
                      Total2010=TotalTotal)

#create denominator counts by race, age, and year and then average over 4 years
newdens <- census3 %>%
  group_by(raceeth_f, Year, Agenum) %>%
  dplyr::summarise(Popsize=sum(Popsize))

combo.newdens <- newdens %>%
  group_by(raceeth_f, Agenum) %>%
  dplyr::summarise(Popsize=sum(Popsize))

combo.newdens$avgPopsize <- combo.newdens$Popsize/4
combo.newdens <- rename(combo.newdens, Popsize4yr=Popsize)

#create numerator counts by race, age, and year and then average over 4 years
preadjnums <- fullfinal %>%
  group_by(raceeth_f, agey, fyear) %>%
  dplyr::summarise(hyst=n())

combo.preadjnums <- preadjnums %>%
  group_by(raceeth_f, agey) %>%
  dplyr::summarise(hyst4yr=sum(hyst))

combo.preadjnums$avghyst <- combo.preadjnums$hyst4yr/4

#merge 4 year denominators to numerator counts
new <- full_join(x=combo.newdens, y=combo.preadjnums, by=c("Agenum"="agey", "raceeth_f"))

#merge age adj proportions with 4yr avg counts
adjrates <- full_join(x=new, y=ageadjprop1, by = c("Agenum", "raceeth_f"))

#calculate crude rates for the whole state (double check)
cruderates.NC <- adjrates %>%
  dplyr::summarise(TotalHyst=sum(avghyst, na.rm=TRUE), TotalPop=sum(avgPopsize, na.rm = TRUE))
cruderates.NC$avgcrude <- (cruderates.NC$TotalHyst/cruderates.NC$TotalPop)*10000

#calculate age adjusted rates 
# standard pop = all NC females ages 18-44 in 2010
adjrates$rate1 <- adjrates$avghyst/adjrates$avgPopsize
adjrates$preageadj <- adjrates$rate1*adjrates$propxage
ageadjrates1 <- adjrates %>%
  group_by(raceeth_f) %>%
  dplyr::summarise(newrate=sum(preageadj, na.rm=TRUE))
ageadjrates1$ageadjrate <- ageadjrates1$newrate*10000

#DENOMINATOR ADJUSTED RATES
#read in hyst prev estimates (using NC pooled estimates, see code at end)
hyst <- read.csv("./Aim1and2/NCrace.4yrpooled.prev.csv", stringsAsFactors = TRUE)
hyst$X <- NULL
#reorder factors (order gets messy when writing and then reading in csv)
hyst$raceeth_f <- factor(hyst$raceeth_f, levels=c("NH White", "NH Black", "Hispanic",
                                                  "NH American Indian", "NH Asian/PI"))

#merge previous race estimates of hyst and popsize (cruderates)
#with race and year specific hyst prev estimates
racerates <- hyst %>%
  select(raceeth_f, prevpct.4yr) %>%
  full_join(cruderates.NC, hyst, by="raceeth_f")

#calc den adj rates
racerates$newden <- (racerates$TotalPop-(racerates$TotalPop*(racerates$prevpct.4yr/100)))
racerates$denrate <- (racerates$TotalHyst/racerates$newden)*10000

#AGE ADJUST RATES AND THEN DENOMINATOR ADJ
#create counts by age and race of population in 2010 (1 year age groups)
agerace1 <- forageadj %>%
  group_by(raceeth_f, Agenum) %>%
  dplyr::summarise(Popsize=sum(Popsize)) #counts stratified by age and race (1 year groups)

racetotal1 <- agerace1 %>%
  group_by(raceeth_f) %>%
  dplyr::summarise(TotalPopRace=sum(Popsize)) #counts stratified by race only (1 year groups)

agetotal1 <- forageadj %>%
  group_by(Agenum) %>%
  dplyr::summarise(PopsizeAge=sum(Popsize)) #counts stratified by age only (1 year groups)

racetotal1 %>% dplyr::summarise(TotalTotal=sum(TotalPopRace)) 
#1765846 = pop size of females in NC between 18-44 years of age
agetotal1$TotalTotal <- 1765846

#calculate proportions (weights) for standardizing
agetotal1$propxage <- agetotal1$PopsizeAge/agetotal1$TotalTotal #compare across race groups,
#standardized to NC age distribution in 2010

#merge all estimates for calculating proportions (weights)
ageadjprop1 <- left_join(agerace1, racetotal1, by="raceeth_f") 
ageadjprop1 <- left_join(ageadjprop1, agetotal1, by="Agenum")

#rename variables for ease of later merging
ageadjprop1 <- rename(ageadjprop1, agerace2010=Popsize, race2010=TotalPopRace, age2010=PopsizeAge,
                      Total2010=TotalTotal)

#create denominator counts by race, age, and year and then average over 4 years
newdens <- census3 %>%
  group_by(raceeth_f, Year, Agenum) %>%
  dplyr::summarise(Popsize=sum(Popsize))

combo.newdens <- newdens %>%
  group_by(raceeth_f, Agenum) %>%
  dplyr::summarise(Popsize=sum(Popsize))

combo.newdens$avgPopsize <- combo.newdens$Popsize/4
combo.newdens <- rename(combo.newdens, Popsize4yr=Popsize)

#create numerator counts by race, age, and year and then average over 4 years
preadjnums <- fullfinal %>%
  group_by(raceeth_f, agey, fyear) %>%
  dplyr::summarise(hyst=n())

combo.preadjnums <- preadjnums %>%
  group_by(raceeth_f, agey) %>%
  dplyr::summarise(hyst4yr=sum(hyst))

combo.preadjnums$avghyst <- combo.preadjnums$hyst4yr/4

#merge 4 year denominators to numerator counts
new <- full_join(x=combo.newdens, y=combo.preadjnums, by=c("Agenum"="agey", "raceeth_f"))

#merge age adj proportions with 4yr avg counts
adjratesl <- full_join(x=new, y=ageadjprop1, by = c("Agenum", "raceeth_f"))

#Direct standardizaton: standard pop = all NC females ages 18-44 in 2010
adjrates$adjhyst <- adjrates$avghyst*adjrates$propxage
adjrates$adjrate <- adjrates$adjhyst/adjrates$avgPopsize
ageadjrates1 <- adjrates %>%
  group_by(raceeth_f) %>%
  dplyr::summarise(newrate=sum(adjrate, na.rm=TRUE))
ageadjrates1$ageadjrate <- ageadjrates1$newrate*10000

#calculate total person time at risk by race/eth and join
ageadjrates1 <- full_join(ageadjrates1, combo.dens, by="raceeth_f")

ageadjrates1$adjhyst4yr <- ageadjrates1$avgPopsize*ageadjrates1$newrate #back calc adj hyst counts to check

#add hyst prev estimates
ageadjrates1 <- full_join(ageadjrates1, hyst, by="raceeth_f")

#calc den adj rates
ageadjrates1$newden <- (ageadjrates1$avgPopsize-
                                (ageadjrates1$avgPopsize*(ageadjrates1$prevpct.4yr/100)))
ageadjrates1$agedenrate <- (ageadjrates1$adjhyst4yr/ageadjrates1$newden)*10000

#create dataset for bootstrapping
forboot <- select(ageadjrates1, raceeth_f, adjhyst4yr, avgPopsize, ageadjrate, agedenrate, newden)

#.......................................................
# Pt 6: Bootstrap rate 95% CI ####
#.......................................................

#read in hyst prev estimates (.csv based on calculations done outside of R - used Lewis (2017)
# paper as guide for calculating 95% CI for pooled BRFSS prevalence estimates)
hystCI <- read.csv(".savefile/NCrace.4yrpooled.prev_95CI.csv", stringsAsFactors = TRUE)
#reorder factors (order gets messy when writing and then reading in csv, apparently)
hystCI$raceeth_f <- factor(hystCI$raceeth_f, levels=c("NH White", "NH Black", "Hispanic",
                                                  "NH American Indian", "NH Asian/PI"))

#add to bootstrap dataset
forboot <- full_join(forboot, hystCI, by="raceeth_f")


#NH WHITE: boot denominator
#simulate dataset using BRFSS data
set.seed(12345)
wh.rnorm <- as.data.frame(rnorm(n=2939243, mean=0.060833588, sd=0.005106829))
wh.rnorm$value <- wh.rnorm$`rnorm(n = 2939243, mean = 0.060833588, sd = 0.005106829)`
wh.rnorm$`rnorm(n = 2939243, mean = 0.060833588, sd = 0.005106829)` <- NULL


#this works for 10000 samples of possible denominator values
set.seed(12345)
boot.den.wh <- function(data, indices){
  d <- sample(data[indices, ], size = 10000, replace = TRUE)
  den <- 1087783
  prev <- d
  newden <- (den-(prev*den))
  rate <- (4618/newden)*10000
  return(rate)
}

z.wh <- boot(data=wh.rnorm, statistic=boot.den.wh, R=10000)
boot.ci(z.wh, type="perc") #percentile method of estimation


#hard coded option: take random sample (w/ replacement) of 10000 prev values
set.seed(12345)
wh.rnorm2 <- as.data.frame(sample(wh.rnorm$value, size=10000, replace = TRUE))
wh.rnorm2$prev <- wh.rnorm2$`sample(wh.rnorm$value, size = 10000, replace = TRUE)`
wh.rnorm2$`sample(wh.rnorm$value, size = 10000, replace = TRUE)` <- NULL
wh.rnorm2$newden <- (1087783-(wh.rnorm2$prev*1087783))
wh.rnorm2$rate <- ((4618/wh.rnorm2$newden)*10000)
wh.rate <- wh.rnorm2[order(wh.rnorm2$rate),]
ci2.5th <- wh.rate$rate[c(250)]
ci97.5th <- wh.rate$rate[c(9750)]


#NH BLACK: boot denominator
#simulate dataset using BRFSS data
set.seed(12345)
bl.rnorm <- as.data.frame(rnorm(n=1124324, mean=0.062036356, sd=0.008667414))
bl.rnorm$value <- bl.rnorm$`rnorm(n = 1124324, mean = 0.062036356, sd = 0.008667414)`
bl.rnorm$`rnorm(n = 1124324, mean = 0.062036356, sd = 0.008667414)` <- NULL
hist(bl.rnorm$value)

#this works for 10000 samples of possible denominator values
set.seed(12345)
boot.den.bl <- function(data, indices){
  d <- sample(data[indices, ], size = 10000, replace = TRUE)
  den <- 434652
  prev <- d
  newden <- (den-(prev*den))
  rate <- (2526/newden)*10000
  return(rate)
}

z.bl <- boot(data=bl.rnorm, statistic=boot.den.bl, R=10000)
boot.ci(z.bl, type="perc") #percentile method of estimation

#hard code check: take random sample (w/ replacement) of 10000 prev values
set.seed(12345)
bl.rnorm2 <- as.data.frame(sample(bl.rnorm$value, size=10000, replace = TRUE))
bl.rnorm2$prev <- bl.rnorm2$`sample(bl.rnorm$value, size = 10000, replace = TRUE)`
bl.rnorm2$`sample(bl.rnorm$value, size = 10000, replace = TRUE)` <- NULL
bl.rnorm2$newden <- (434652-(bl.rnorm2$prev*434652))
bl.rnorm2$rate <- ((2526/bl.rnorm2$newden)*10000)
bl.rate <- bl.rnorm2[order(bl.rnorm2$rate),]
ci2.5th <- bl.rate$rate[c(250)]
ci97.5th <- bl.rate$rate[c(9750)]


#NH AIAN: boot denominator
#simulate dataset using BRFSS data
set.seed(12345)
aian.rnorm <- as.data.frame(rnorm(n=51127, mean=0.09641432, sd=0.038156278))
aian.rnorm$value <- aian.rnorm$`rnorm(n = 51127, mean = 0.09641432, sd = 0.038156278)`
aian.rnorm$`rnorm(n = 51127, mean = 0.09641432, sd = 0.038156278)` <- NULL
hist(aian.rnorm$value)

#this works for 10000 samples of possible denominator values
set.seed(12345)
boot.den.ai <- function(data, indices){
  d <- sample(data[indices, ], size = 10000, replace = TRUE)
  den <- 23211
  prev <- d
  newden <- (den-(prev*den))
  rate <- (179/newden)*10000
  return(rate)
}

z.ai <- boot(data=ai.rnorm, statistic=boot.den.ai, R=10000)
boot.ci(z.ai, type="perc") #percentile method of estimation

#hard code check: take random sample (w/ replacement) of 10000 prev values
set.seed(12345)
aian.rnorm2 <- as.data.frame(sample(aian.rnorm$value, size=10000, replace = TRUE))
aian.rnorm2$prev <- aian.rnorm2$`sample(aian.rnorm$value, size = 10000, replace = TRUE)`
aian.rnorm2$`sample(aian.rnorm$value, size = 10000, replace = TRUE)` <- NULL
aian.rnorm2$newden <- (23211-(aian.rnorm2$prev*23211))
aian.rnorm2$rate <- ((179/aian.rnorm2$newden)*10000)
aian.rate <- aian.rnorm2[order(aian.rnorm2$rate),]
ci2.5th <- aian.rate$rate[c(250)]
ci97.5th <- aian.rate$rate[c(9750)]


#NH HISPANIC: boot denominator
#simulate dataset using BRFSS data
set.seed(12345)
hisp.rnorm <- as.data.frame(rnorm(n=524736, mean=0.05086678, sd=0.013330633))
hisp.rnorm$value <- hisp.rnorm$`rnorm(n = 524736, mean = 0.05086678, sd = 0.013330633)`
hisp.rnorm$`rnorm(n = 524736, mean = 0.05086678, sd = 0.013330633)` <- NULL
hisp.rnorm <- as.data.frame(hisp.rnorm[ which(hisp.rnorm$value >= 0),])
hisp.rnorm$value <- hisp.rnorm$`hisp.rnorm[which(hisp.rnorm$value >= 0), ]`
hisp.rnorm$`hisp.rnorm[which(hisp.rnorm$value >= 0), ]` <- NULL
hist(hisp.rnorm$value)

#this works for 10000 samples of possible denominator values
set.seed(12345)
boot.den.hisp <- function(data, indices){
  d <- sample(data[indices, ], size = 10000, replace = TRUE)
  den <- 176151
  prev <- d
  newden <- (den-(prev*den))
  rate <- (335/newden)*10000
  return(rate)
}

z.hisp <- boot(data=hisp.rnorm, statistic=boot.den.hisp, R=10000)
boot.ci(z.hisp, type="perc") #percentile method of estimation

#hard code to check: take random sample (w/ replacement) of 10000 prev values
set.seed(12345)
hisp.rnorm2 <- as.data.frame(sample(hisp.rnorm$value, size=10000, replace = TRUE))
hisp.rnorm2$prev <- hisp.rnorm2$`sample(hisp.rnorm$value, size = 10000, replace = TRUE)`
hisp.rnorm2$`sample(hisp.rnorm$value, size = 10000, replace = TRUE)` <- NULL
hisp.rnorm2$newden <- (176151-(hisp.rnorm2$prev*176151))
hisp.rnorm2$rate <- ((335/hisp.rnorm2$newden)*10000)
hisp.rate <- hisp.rnorm2[order(hisp.rnorm2$rate),]
ci2.5th <- hisp.rate$rate[c(250)]
ci97.5th <- hisp.rate$rate[c(9750)]


#NH AS/PI: boot denominator
#simulate 
set.seed(12345)
as.rnorm <- as.data.frame(rnorm(n=105119, mean=0.00811241, sd=0.006329597))
as.rnorm$value <- as.rnorm$`rnorm(n = 105119, mean = 0.00811241, sd = 0.006329597)`
as.rnorm$`rnorm(n = 105119, mean = 0.00811241, sd = 0.006329597)` <- NULL
as.rnorm <- as.data.frame(as.rnorm[ which(as.rnorm$value >= 0),])
as.rnorm$value <- as.rnorm$`as.rnorm[which(as.rnorm$value >= 0), ]`
as.rnorm$`as.rnorm[which(as.rnorm$value >= 0), ]` <- NULL
hist(as.rnorm$value)

#this works for 10000 samples of possible denominator values
set.seed(12345)
boot.den.as <- function(data, indices){
  d <- sample(data[indices, ], size = 10000, replace = TRUE)
  den <- 63496
  prev <- d
  newden <- (den-(prev*den))
  rate <- (51/newden)*10000
  return(rate)
}

z.as <- boot(data=as.rnorm, statistic=boot.den.as, R=2000)
boot.ci(z.as, type="perc") #percentile method of estimation

#hard code to check: take random sample (w/ replacement) of 10000 prev values
set.seed(12345)
as.rnorm2 <- as.data.frame(sample(as.rnorm$value, size=10000, replace = TRUE))
as.rnorm2$prev <- as.rnorm2$`sample(as.rnorm$value, size = 10000, replace = TRUE)`
as.rnorm2$`sample(as.rnorm$value, size = 10000, replace = TRUE)` <- NULL
as.rnorm2$newden <- (63496-(as.rnorm2$prev*63496))
as.rnorm2$rate <- ((51/as.rnorm2$newden)*10000)
as.rate <- as.rnorm2[order(as.rnorm2$rate),]
ci2.5th <- as.rate$rate[c(250)]
ci97.5th <- as.rate$rate[c(9750)]

#.......................................................
# Pt 7: Bootstrap rate difference 95% CI ####
#.......................................................

#difference estimation (B-W Rate)
#this works for 10000 samples of possible denominator values for B and W and then subtractions
#merge B and W datasets
bwdiff <- cbind(wh.rnorm2,bl.rnorm2)

bw1 <- function(data,i){
  d1 <- data[i,1]
  den1 <- 434652
  prev1 <- d1
  newden1 <- (den1-(prev1*den1))
  rate1 <- (2526/newden1)*10000
  d2 <- data[i,2]
  den2 <- 1087783
  prev2 <- d2
  newden2 <- (den2-(prev2*den2))
  rate2 <- (4618/newden2)*10000
  return(rate1-rate2)
}

set.seed(12345)
bw2 <- boot(data=bwdiff, statistic=bw1, R=2000)
boot.ci(bw2, type="perc") #percentile method of estimation

#AI & W difference
aiwdiff <- cbind(wh.rnorm2,ai.rnorm2)

aiw1 <- function(data,i){
  d1 <- data[i,2]
  den1 <- 23211
  prev1 <- d1
  newden1 <- (den1-(prev1*den1))
  rate1 <- (179/newden1)*10000
  d2 <- data[i,1]
  den2 <- 1087783
  prev2 <- d2
  newden2 <- (den2-(prev2*den2))
  rate2 <- (4618/newden2)*10000
  return(rate1-rate2)
}

set.seed(12345)
aiw2 <- boot(data=aiwdiff, statistic=aiw1, R=2000)
boot.ci(aiw2, type="perc") #percentile method of estimation


#Hisp & W difference
hispwdiff <- cbind(wh.rnorm2,hisp.rnorm2)

hispw1 <- function(data,i){
  d1 <- data[i,2]
  den1 <- 176151
  prev1 <- d1
  newden1 <- (den1-(prev1*den1))
  rate1 <- (335/newden1)*10000
  d2 <- data[i,1]
  den2 <- 1087783
  prev2 <- d2
  newden2 <- (den2-(prev2*den2))
  rate2 <- (4618/newden2)*10000
  return(rate1-rate2)
}

set.seed(12345)
hispw2 <- boot(data=hispwdiff, statistic=hispw1, R=2000)
boot.ci(hispw2, type="perc") #percentile method of estimation


#AS & W difference
aswdiff <- cbind(wh.rnorm2,as.rnorm2)

asw1 <- function(data,i){
  d1 <- data[i,2]
  den1 <- 63496
  prev1 <- d1
  newden1 <- (den1-(prev1*den1))
  rate1 <- (51/newden1)*10000
  d2 <- data[i,1]
  den2 <- 1087783
  prev2 <- d2
  newden2 <- (den2-(prev2*den2))
  rate2 <- (4618/newden2)*10000
  return(rate1-rate2)
}

set.seed(12345)
asw2 <- boot(data=aswdiff, statistic=asw1, R=2000)
boot.ci(asw2, type="perc") #percentile method of estimation
#.......................................................
# NOTES:



#...................................
# Pt 8: Comorbidity (CCI) calculations ####
#...................................
#create smaller dataset for ease of working
myvars <- names(fullfinal) %in% c("agey", "raceeth_f", "admitdx", "shepsid", "fyear", (paste("diag", 1:25, sep="")))
smallfull <- fullfinal[myvars]
smallfull$shepsid <- as.character(smallfull$shepsid)

#turn wide into long format
longfull <- smallfull %>%
  gather(key = codenum, value=code, admitdx, diag1:diag25) #ready data, 791154 obs

#drop missing values (double check b/c this results in loss of lots of observations)
describe(longfull$code)
longfull$code[longfull$code == ""] <- NA #reset blanks to missing 
nomissfull <- longfull[ which(!is.na(longfull$code)), ] #this results in loosing a lot but that is okay, obs=229792
describe(nomissfull$code)

#check that I didn't loose any women...
describe(longfull$shepsid) #30429 unique surgeries
describe(smallfull$shepsid) #30429 unique surgeries
describe(nomissfull$shepsid) #30428 unique surgeries, lost one b/c one surgery had all missing dx codes

#avg count of codes per patient
count <- count(fullfinal, c("shepsid"))
count2 <- count(nomissfull, c("shepsid"))
cat("average count of ICD dx codes per patient is: ", count2$n/count$n)
#7.550301 codes per patient

y <- count(nomissfull, vars = code)
kable(head(y[order(-y$n),], n=5), row.names=F, 
      caption = "Top 5 most popular ICD-9 dx codes")

#|vars  |     n|
#|:-----|-----:|
#|6262  | 25378| excessive bleeding
#|2189  | 16004| fibroids
#|6253  | 11000| dysmenorrhea
#|6259  | 10578| unspecified genital-related symptoms
#|V5869 | 7749 | long term use of meds


#use icd package - gets comorbidity score but not indicators for specific dx
nomissfull$icd9 <- nomissfull$code #rename for using icd9 package
charlson <- icd_charlson(nomissfull, visit_name = "shepsid", scoring_system = "charlson", 
                         return_df = TRUE, stringsAsFactors = getOption("stringsAsFactors"))

quan <- icd_charlson(nomissfull, visit_name = "shepsid", scoring_system = "quan",
                     return_df = TRUE, stringsAsFactors = getOption("stringsAsFactors"))
hist(quan$Charlson)
hist(charlson$Charlson)

names(icd9_map_quan_deyo)
names(icd9_map_quan_elix) 
names(icd9_map_ahrq) #preferring this one (same as quan_elix, except no arrhythmias)

#using Quan version
trial2 <- full_join(nomissfull, quan, by = "shepsid")

#these are clunky but work
fig.trial2 <- trial2 %>%
  group_by(shepsid, raceeth_f) %>%
  dplyr::summarise(score=sum(Charlson))
racetotal2 <- fig.trial2 %>%
  group_by(raceeth_f) %>%
  dplyr::summarise(Total=sum(score))

fig2 <- fig.trial2 %>%
  group_by(raceeth_f, score) %>%
  dplyr::summarise(num=n())

denominator1 <- trial2 %>%
  group_by(shepsid, raceeth_f) %>%
  dplyr::summarise(all=n())
denominator2 <- denominator1 %>%
  group_by(raceeth_f) %>%
  dplyr::summarise(all=n())
fig2.2 <- full_join(fig2, denominator2, by = "raceeth_f")
fig2.2$perc <- fig2.2$num/fig2.2$all

write.csv(fig2.2, ".savefile/CCS_Quan_allyearsXrace.csv")
#..........................................................
#NOTES: Used excel file to group scores into categories and 
#calculate mean scores


#..........................................................
# Pt 9: Sensitivity Analyses ####
#..........................................................
#misclassification of 35% AI/AN as White or Black 
#179.26 AI/AN procedures (original numerator estimate); 96 others misclassified
# 48 as White and 48 as Black

#misclassification of 18% AI/AN as White or Black
#179 AI/AN procedures (original numerator estimate); 39 others misclassified
# 19 as White and 20 as Black

#misclassification of 5% of Hispanics as White

#omission in denominator, by race/ethnicity (see appendix table)

#underreporting of hyst prev among non-White groups (see appendix table)

#write.csv(forboot, "./Aim1and2/useforsensitivity.csv")
#used excel to finish and create table



#Restrict to non-border counties and calculate race/eth specific rates
#see below
#this code creates 1-yr incidence rates by border or non border status

#CRUDE RATES
#merge border indicator
fullfinalb <- full_join(fullfinal, names, by = c("ptcnty" = "cnty"))

#create numerator for border and non-border counties then average over the 4 years
# to create the average population during the 4 year period
nums <- fullfinalb %>%
  group_by(border) %>%
  dplyr::summarise(total4yr=n())
nums$avghyst <- nums$total4yr/4

#create denominator for border and non-border counties then average over the 4 years
# to create the average population during the 4 year period
dens <- census3 %>%
  group_by(border) %>%
  dplyr::summarise(Popsize=sum(Popsize))
dens$avgPopsize <- dens$Popsize/4

#merge by cnty and calc rates
crude <- full_join(dens, nums, by = "border")

#calculate rates
crude$rate <- crude$avghyst/crude$avgPopsize
crude$totrateby10K <- crude$rate*10000
crude$totrateby10K

#AGE ADJUSTED RATES
# Calc age-adj cnty-level rates (5 year groups age adj) for border and non border cnties
#for age adj: create counts by age grp of population in 2010
age <- census4 %>%
  group_by(agegrp_f) %>%
  dplyr::summarise(Popsize=sum(Popsize)) %>%
  as.data.frame

#for age adj: calculate total population for the state
NCstate <- age %>%
  dplyr::summarise(AllNC=sum(Popsize))
NCstate
#1765846 is the pop of NC in 2010

#for age adj: calculate proportion for direct age standardization (prop by age group for the state)
age$propxage <- age$Popsize/1765846 #to be able to compare across all counties

ageadjprop <- age #rename data
ageadjprop$TotalPop2010 <- 1765846 #add for reference

#for age adj: rename variables for ease of later merging and duplicate dataset for each year
ageadjprop <- rename(ageadjprop, Popsize2010=Popsize)
#use ageadjprop as denominators

#for rates: create numerator counts by county, age (5 year groups), and year
#(agegrp and agegrp_f)
fullfinalb$agegrp[fullfinalb$agey < 20] <- 1 #18-19, group 1
fullfinalb$agegrp[fullfinalb$agey > 19 & fullfinalb$agey < 25] <- 2 #20-24, grp 2
fullfinalb$agegrp[fullfinalb$agey > 24 & fullfinalb$agey < 30] <- 3 #25-29, grp 3
fullfinalb$agegrp[fullfinalb$agey > 29 & fullfinalb$agey < 35] <- 4 #30-34, grp 4
fullfinalb$agegrp[fullfinalb$agey > 34 & fullfinalb$agey < 40] <- 5 #35-39, grp 5
fullfinalb$agegrp[fullfinalb$agey > 39 & fullfinalb$agey < 45] <- 6 #40-44, grp 6
describe(fullfinalb$agegrp)

fullfinalb$agegrp_f <- factor(fullfinalb$agegrp, levels = 1:6,
                              labels = c("18-19 yrs", "20-24 yrs", "25-29 yrs", "30-34 yrs",
                                         "35-39 yrs", "40-44 yrs"))
table(fullfinalb$agegrp_f, useNA = "always")

#for rates: create numerator counts by age and year and border status
agenums <- fullfinalb %>%
  group_by(border, agegrp_f) %>%
  dplyr::summarise(totalhyst=n()) %>%
  as.data.frame

#for rates: merge age adjustment to numerator counts
adjrates <- full_join(x=ageadjprop, y=agenums, by=c("agegrp_f")) 

#for rates: create denominator for each cnty
#this will be the total population by cnty and age group during the 4 year period
agedens <- census3 %>%
  group_by(border, agegrp_f) %>%
  dplyr::summarise(Popsize4yr=sum(Popsize)) %>%
  as.data.frame

#for rates: merge dens to age adjustment to numerator counts
adjrates <- full_join(x=adjrates, y=agedens, by=c("agegrp_f", "border")) 

#divide by 4 to get average 1 year rate estimates
adjrates$avgPopsize <- adjrates$Popsize4yr/4
adjrates$avghyst <- adjrates$totalhyst/4

#standard pop = all NC females ages 18-44 in 2010
adjrates$rate1 <- adjrates$avghyst/adjrates$avgPopsize #crude age specific rate by border status and age grp
adjrates$preageadj <- adjrates$rate1*adjrates$propxage #mult crude rate by age specific prop

ageadjrates <- adjrates %>%
  group_by(border) %>%
  dplyr::summarise(newrate=sum(preageadj, na.rm=TRUE)) %>%
  as.data.frame
ageadjrates$ageadjrate4yr <- ageadjrates$newrate*10000
ageadjrates <- select(ageadjrates, -newrate) #clean up a bit

#DENOMINATOR ADJUSTED RATES
#read in hyst prev estimates (using NC pooled estimates)
hyst <- read.csv("./Aim1and2/NCrace.4yrpooled.prev.csv", stringsAsFactors = TRUE)
hyst$X <- NULL
#reorder factors (order gets messy when writing and then reading in csv, apparently)
hyst$raceeth_f <- factor(hyst$raceeth_f, levels=c("NH White", "NH Black", "Hispanic",
                                                  "NH American Indian", "NH Asian/PI"))

#calculate statewide prevalences
hyst2 <- hyst %>%
  dplyr::summarise(Hyst=sum(Hyst), Pop=sum(Pop))

hyst2$prev.4yr <- hyst2$Hyst/hyst2$Pop #calc prev as percent
hyst2$prevpct.4yr <- hyst2$prev.4yr*100

#add state prev to crude rates
denomadj <- crude #create df
denomadj$prev <- 5.923164 #add prev to df

denomadj$newden <- (denomadj$avgPopsize-(denomadj$avgPopsize*(denomadj$prev/100)))
denomadj$denrate <- (denomadj$avghyst/denomadj$newden)*10000

# NOTES: denomadj df includes crude rates and denominator
# adjust rates

#DENOMINATOR + AGE ADJUSTED RATES
#start with adjrates df and add prev estimate
denomageadj <- adjrates
denomageadj$prev <- 5.923164 #add prev to df

#calculate new denominators for age-specific den adj rates
denomageadj$newden <- (denomageadj$avgPopsize-
                         (denomageadj$avgPopsize*(denomageadj$prev/100))) 

#Direct standardization: standard pop = all NC females ages 18-44 in 2010
denomageadj$rate1 <- denomageadj$avghyst/denomageadj$newden
denomageadj$preageadj <- denomageadj$rate1*denomageadj$propxage

agedenadj <- denomageadj %>%
  group_by(border) %>%
  dplyr::summarise(newrate=sum(preageadj, na.rm=TRUE)) %>%
  as.data.frame
agedenadj$agedenadjrate <- agedenadj$newrate*10000
agedenadj$newrate <- NULL #clean up a bit
#.......................................................
# NOTES: 


#.............................................
# Prologue: BRFSS Prevalence calculations ####
#..............................................
# Description: Create Hyst prevalence estimates

#download data from: https://www.cdc.gov/brfss/annual_data/annual_data.htm
#load data (note, smaller datasets saved as .csv for easy load - see below)
brfss16 <- sasxport.get("CDBRFS16.XPT")
brfss14 <- sasxport.get("CDBRFS14.XPT")
brfss12 <- sasxport.get("CDBRFS12.XPT")
brfss10 <- sasxport.get("CDBRFS10.XPT")
names(brfss16)
names(brfss14)
names(brfss12)
names(brfss10)

# Create smaller datasets & merge years

vars1 = c("x.state", "x.geostr", "x.psu", "x.ststr", "x.wt2", "age", "sex", "hadhyst2",
          "educa", "income2", "race2", "x.llcpwt")

vars2 = c("x.state", "x.geostr", "x.psu", "x.ststr", "x.wt2", "age", "sex", "hadhyst2",
          "educa", "income2", "race2", "x.finalwt")

vars4 = c("x.state", "x.geostr", "x.psu", "x.ststr", "sex", "hadhyst2",
          "educa", "income2", "race2", "x.llcpwt", "x.ageg5yr")

brfss12 <- brfss12[vars1]
brfss10 <- brfss10[vars2]
brfss14 <- brfss14[vars4]

#create year indicator
brfss12$year <- "2012"
brfss10$year <- "2010"
brfss14$year <- "2014"

#create csvs for easier load in
write.csv(brfss12, ".savefile/brfss12.csv")
write.csv(brfss10, ".savefile/brfss10.csv")
write.csv(brfss14, ".savefile/brfss14.csv")

#start here once data has been converted to .csv
brfss12 <- read.csv(".savefile/brfss12.csv", header = T, stringsAsFactors = F)
brfss10 <- read.csv(".savefile/brfss10.csv", header = T, stringsAsFactors = F)
brfss14 <- read.csv(".savefile/brfss14.csv", header = T, stringsAsFactors = F)

#create year variable for 2014 and rename race variable in 2014 and 2016 to match others
# rename 5 yr age group variable for 2016
brfss14$year <- 2014
brfss14 <- rename(brfss14, race2=x.race)
brfss16 <- rename(brfss16, race2=x.race)
brfss16 <- rename(brfss16, x.ageg5yr=x.ageg5yr.1)

#append datasets
allhyst <- dplyr::bind_rows(brfss12, brfss10, brfss14)
names(allhyst)

# Note: brfss12 and brfss14 x.llcpwt variable, brfss10 has
# x.finalwt instead
# brfss14 doesn't provide one year age groupings; only provides
# 5 year age groupings (x.ageg5yr)
# brfss14 had x.race; brfss10 and brfss12 had race2


# Recode Variables, 2010-16
# recode race into 6 categories (NH white, NH Black, Hisp, NH Asian/PI, NH AI/AN, NH Other)
# (raceeth); create factor too (raceeth_f)

str(allhyst$race2)
describe(allhyst$race2)
racehelper <- read.csv("./Aim1and2/Racehelper.csv", header=TRUE, stringsAsFactors = TRUE)
racecoded <- full_join(allhyst, racehelper, by= c("year", "race2"="oldcode"))

#check coding
table(racecoded$raceeth_f, racecoded$year, useNA = "always")
table(allhyst$race2, allhyst$year == 2010, useNA = "always")
table(allhyst$race2, allhyst$year == 2012, useNA = "always")
table(allhyst$race2, allhyst$year == 2014, useNA = "always")

#recode race factor
racecoded$raceeth_f <- NULL
racecoded$raceeth_f <- factor(racecoded$raceeth, levels=1:6,
                              labels = c("NH White", "NH Black", "Hispanic", "NH American Indian",
                                         "NH Asian/PI", "NH Other")) #create raceeth_f
table(racecoded$raceeth_f, useNA = "always")

# Age: recode missing
describe(racecoded$age)
racecoded$age[racecoded$age == 7 | racecoded$age == 9 | is.nan(racecoded$age)] <- NA
racecoded$agegroup[racecoded$age > 17 & racecoded$age < 25] <- 1
racecoded$agegroup[racecoded$age > 24 & racecoded$age < 30] <- 2
racecoded$agegroup[racecoded$age > 29 & racecoded$age < 35] <- 3
racecoded$agegroup[racecoded$age > 34 & racecoded$age < 40] <- 4
racecoded$agegroup[racecoded$age > 39 & racecoded$age < 45] <- 5
racecoded$agegroup[racecoded$age > 44] <- 6
racecoded$agegroup[is.na(racecoded$age)] <- NA
describe(racecoded$agegroup)

describe(racecoded$x.ageg5yr)
racecoded$agegroup[racecoded$x.ageg5yr == 1] <- 1
racecoded$agegroup[racecoded$x.ageg5yr == 2] <- 2
racecoded$agegroup[racecoded$x.ageg5yr == 3] <- 3
racecoded$agegroup[racecoded$x.ageg5yr == 4] <- 4
racecoded$agegroup[racecoded$x.ageg5yr == 5] <- 5
racecoded$agegroup[racecoded$x.ageg5yr > 5 & racecoded$x.ageg5yr < 14] <- 6
racecoded$agegroup[racecoded$x.ageg5yr == 14 & racecoded$year == 2014] <- NA
describe(racecoded$agegroup)

#check coding
table(racecoded$agegroup, racecoded$year, useNA = "always")
table(racecoded$x.ageg5yr, racecoded$year == 2014, useNA = "always")
#table(racecoded$x.ageg5yr, racecoded$year == 2016, useNA = "always")
table(racecoded$age, (racecoded$age > 17 & racecoded$age < 25 & racecoded$year == 2010), useNA = "always")
table(racecoded$age, (racecoded$age > 17 & racecoded$age < 25 & racecoded$year == 2012), useNA = "always")

# Had hysterectomy: recode missing as NA; make indicator (hyst) and factor (hyst_f)
table(racecoded$hadhyst2, useNA = "always")
racecoded$hadhyst2[racecoded$hadhyst2 == 7 |racecoded$hadhyst2 == 9 | 
                     is.na(racecoded$hadhyst2) | is.nan(racecoded$hadhyst2)] <- NA
racecoded$hyst[racecoded$hadhyst2 == 1] <- 1
racecoded$hyst[racecoded$hadhyst2 == 2] <- 0 
table(racecoded$hyst, useNA = "always")

racecoded$hyst_f <- factor(racecoded$hyst, levels = 0:1, 
                           labels = c("No hyst", "Had hyst"))
table(racecoded$hyst_f, useNA = "always")

write.csv(racecoded, "./Aim1and2/brfssallyearsrecoded.csv")
#Notes: Saved recoded dataset as (brfssallyearsrecoded.csv)

# Create design objects
#separate into years
fordesign.10 <- racecoded[which(racecoded$year==2010),]
fordesign.12 <- racecoded[which(racecoded$year==2012),]
fordesign.14 <- racecoded[which(racecoded$year==2014),]

#subset to nonmissing weights - 2010
sub10 <- c("agegroup", "raceeth", "raceeth_f", "hyst", "hyst_f", "sex", "x.state", "x.psu",
           "x.ststr", "x.geostr", "year", "x.finalwt")
pre10 <- fordesign.10[sub10]
post10 <- pre10[ which(!is.na(pre10$x.ststr) & !is.na(pre10$x.finalwt) & 
                         !is.na(pre10$x.psu)), ]

#create design object - 2010
design.10 <- post10 %>%
  as_survey_design(ids=x.psu,
                   weight=x.finalwt,
                   nest=T,
                   strata=x.ststr,
                   variables= c(agegroup, raceeth, raceeth_f, hyst, hyst_f, sex, x.state))

#subset to nonmissing weights - 2012
sub12 <- c("agegroup", "raceeth", "raceeth_f", "hyst", "hyst_f", "sex", "x.state", "x.psu",
           "x.ststr", "x.geostr", "year", "x.llcpwt")
pre12 <- fordesign.12[sub12]
post12 <- pre12[ which(!is.na(pre12$x.ststr) & !is.na(pre12$x.llcpwt) & 
                         !is.na(pre12$x.psu)), ]

#create design object - 2012
design.12 <- post12 %>%
  as_survey_design(ids=x.psu,
                   weight=x.llcpwt,
                   nest=T,
                   strata=x.ststr,
                   variables= c(agegroup, raceeth, raceeth_f, hyst, hyst_f, sex, x.state))

#subset to nonmissing weights - 2014
sub14 <- c("agegroup", "raceeth", "raceeth_f", "hyst", "hyst_f", "sex", "x.state", "x.psu",
           "x.ststr", "x.geostr", "year", "x.llcpwt")
pre14 <- fordesign.14[sub14]
post14 <- pre14[ which(!is.na(pre14$x.ststr) & !is.na(pre14$x.llcpwt) & 
                         !is.na(pre14$x.psu)), ]

#create design object - 2014
design.14 <- post14 %>%
  as_survey_design(ids=x.psu,
                   weight=x.llcpwt,
                   nest=T,
                   strata=x.ststr,
                   variables= c(agegroup, raceeth, raceeth_f, hyst, hyst_f, sex, x.state))


# Prevalence estimates - 2010 
# Age: recode missing
describe(fordesign.10$age)
fordesign.10$age[fordesign.10$age == 7 | fordesign.10$age == 9 | is.nan(fordesign.10$age)] <- NA
fordesign.10$agegroup[fordesign.10$age > 17 & fordesign.10$age < 25] <- 1
fordesign.10$agegroup[fordesign.10$age > 24 & fordesign.10$age < 30] <- 2
fordesign.10$agegroup[fordesign.10$age > 29 & fordesign.10$age < 35] <- 3
fordesign.10$agegroup[fordesign.10$age > 34 & fordesign.10$age < 40] <- 4
fordesign.10$agegroup[fordesign.10$age > 39 & fordesign.10$age < 45] <- 5
fordesign.10$agegroup[fordesign.10$age > 44] <- 6
fordesign.10$agegroup[is.na(fordesign.10$age)] <- NA
describe(fordesign.10$agegroup)

#look at racial distribution of hyst counts in NC
check.2010 <- fordesign.10 %>%
  filter(sex == 2 & (agegroup > 0 & agegroup < 6) & x.state == 37) %>%
  group_by(hyst_f, raceeth_f) %>%
  dplyr::summarise(total=n())

#look at racial distribution of hyst counts in NC by age group
check.2010 <- fordesign.10 %>%
  filter(sex == 2 & x.state == 37) %>%
  group_by(hyst_f, raceeth_f, agegroup) %>%
  dplyr::summarise(total=n())

#Overall hyst prevalence - NC
NC.2010 <- design.10 %>%
  filter(sex == 2 & (agegroup > 0 & agegroup < 6) & x.state == 37) %>%
  group_by(hyst_f) %>%
  summarize(prevalence = survey_mean(na.rm=T, vartype = c("ci", "var", "se", "cv")),
            N = survey_total(na.rm=T))
NC.2010 #overall prevalence (7.14% 95% CI: 5.60-8.68%)

#hyst prev by race - NC
NC.race.2010 <- design.10 %>%
  filter(sex == 2 & (agegroup > 0 & agegroup < 6) & x.state == 37) %>%
  group_by(raceeth_f, hyst_f) %>%
  summarize(prevalence = survey_mean(na.rm=T, vartype = c("ci", "var", "se", "cv")),
            N = survey_total(na.rm=T))

#NH White: 0.076307630 95% CI: 0.057755093-0.09486017
#NH Black 0.072442960 95%CI: 0.033312132-0.11157379
#Hispanic: 0.050433471 95% CI: 0.001727882-0.09913906
#NH American Indian: 0.143051389 95% CI: -0.041068801-0.32717158
#NH Asian/PI: 0.020844181 95% CI: -0.011311328-0.05299969

#hyst prev by race - US specific
us.2010 <- design.10 %>%
  filter(sex == 2 & (agegroup > 0 & agegroup < 6)) %>%
  group_by(raceeth_f, hyst_f) %>%
  summarise(prevalence = survey_mean(na.rm=T, vartype = "ci", level = 0.95),
            N = survey_total(na.rm=T, vartype = "ci"))
us.2010 <- us.2010[ which(us.2010$hyst_f == "Had hyst"),]

# Prevalence estimates - 2012 
# Age: recode missing
describe(fordesign.12$age)
fordesign.12$age[fordesign.12$age == 7 | fordesign.12$age == 9 | is.nan(fordesign.12$age)] <- NA
fordesign.12$agegroup[fordesign.12$age > 17 & fordesign.12$age < 25] <- 1
fordesign.12$agegroup[fordesign.12$age > 24 & fordesign.12$age < 30] <- 2
fordesign.12$agegroup[fordesign.12$age > 29 & fordesign.12$age < 35] <- 3
fordesign.12$agegroup[fordesign.12$age > 34 & fordesign.12$age < 40] <- 4
fordesign.12$agegroup[fordesign.12$age > 39 & fordesign.12$age < 45] <- 5
fordesign.12$agegroup[fordesign.12$age > 44] <- 6
fordesign.12$agegroup[is.na(fordesign.12$age)] <- NA
describe(fordesign.12$agegroup)

#look at racial distribution of hyst counts in NC
check.2012 <- fordesign.12 %>%
  filter(sex == 2 & (agegroup > 0 & agegroup < 6) & x.state == 37) %>%
  group_by(hyst_f, raceeth_f) %>%
  dplyr::summarise(total=n())

#look at racial distribution of hyst counts in NC by age group
check.2012 <- fordesign.12 %>%
  filter(sex == 2 & x.state == 37) %>%
  group_by(hyst_f, raceeth_f, agegroup) %>%
  dplyr::summarise(total=n())

#Overall hyst prevalence - NC
NC.2012 <- design.12 %>%
  filter(sex == 2 & (agegroup > 0 & agegroup < 6) & x.state == 37) %>%
  group_by(hyst_f) %>%
  summarize(prevalence = survey_mean(na.rm=T, vartype = c("ci", "var", "se", "cv")),
            N = survey_total(na.rm=T))
NC.2012 #overall prevalence (6.18% 95% CI: 4.80-7.56%)

#hyst prev by race - NC
NC.race.2012 <- design.12 %>%
  filter(sex == 2 & (agegroup > 0 & agegroup < 6) & x.state == 37) %>%
  group_by(raceeth_f, hyst_f) %>%
  summarize(prevalence = survey_mean(na.rm=T, vartype = c("ci", "var", "se", "cv")),
            N = survey_total(na.rm=T))

#NH White: 0.076307630 95% CI: 0.057755093-0.09486017  
#NH Black 0.072442960 95%CI:0.033312132-0.11157379
#Hispanic: 0.050433471 95% CI0.001727882-0.09913906 
#NH American Indian: 0.143051389 95% CI: -0.041068801-0.32717158 
#NH Asian/PI: 0.020844181 95% CI: -0.011311328-0.05299969

# Prevalence estimates - 2014
#look at racial distribution of hyst counts in NC
check.2014 <- fordesign.14 %>%
  filter(sex == 2 & (agegroup > 0 & agegroup < 6) & x.state == 37) %>%
  group_by(hyst_f, raceeth_f) %>%
  dplyr::summarise(total=n())

#look at racial distribution of hyst counts in NC by age group
check.2014 <- fordesign.14 %>%
  filter(sex == 2 & x.state == 37) %>%
  group_by(hyst_f, raceeth_f, agegroup) %>%
  dplyr::summarise(total=n())

#Overall hyst prevalence - NC
NC.2014 <- design.14 %>%
  filter(sex == 2 & (agegroup > 0 & agegroup < 6) & x.state == 37) %>%
  group_by(hyst_f) %>%
  summarize(prevalence = survey_mean(na.rm=T, vartype = c("ci", "var", "se", "cv")),
            N = survey_total(na.rm=T))
NC.2014 #overall prevalence (4.55% 95% CI: 3.15-5.95%)

#hyst prev by race - NC
NC.race.2014 <- design.14 %>%
  filter(sex == 2 & (agegroup > 0 & agegroup < 6) & x.state == 37) %>%
  group_by(raceeth_f, hyst_f) %>%
  summarize(prevalence = survey_mean(na.rm=T, vartype = c("ci", "var", "se", "cv")),
            N = survey_total(na.rm=T))

#NH White: 0.0443700 95% CI: 0.024968559-0.06377145
#NH Black 0.05877764 95%CI: 0.029288282-0.08826700
#Hispanic: 0.03603022 95% CI: 0.004706218-0.06735423
#NH American Indian: 0.00 95% CI: 
#NH Asian/PI: 0.00 95% CI: 



# Pooled Estimates (2010-14)
#add year and merge years
NC.race.2010$year <- 2010
NC.race.2012$year <- 2012
NC.race.2014$year <- 2014

#merge data sets together
bind1 <- bind_rows(NC.race.2010, NC.race.2012)
NC.race.allyears <- bind_rows(bind1, NC.race.2014)
table(NC.race.allyears$year)

#subset to only 5 groups
NC.race.allyears <- NC.race.allyears[ which(NC.race.allyears$raceeth_f != "NH Other"),]

#reorder factors (order gets messy when writing and then reading in csv, apparently)
#NC.race.allyears$raceeth_f <- factor(NC.race.allyears$raceeth_f, levels=c("NH White", "NH Black", "Hispanic",
#                                              "NH American Indian", "NH Asian/PI"))

#turn prev into percents
NC.race.allyears$prevpct <- NC.race.allyears$prevalence*100
NC.race.allyears$lowerpct <- NC.race.allyears$prevalence_low*100
NC.race.allyears$upperpct <- NC.race.allyears$prevalence_upp*100

pool <- NC.race.allyears %>%
  group_by(raceeth_f, hyst_f) %>%
  dplyr::summarise(Total=sum(N))

pool2 <- pool %>%
  group_by(raceeth_f) %>%
  dplyr::summarise(sum=sum(Total))

pool <- filter(pool, hyst_f=="Had hyst") 

pool3 <- full_join(pool, pool2, by=c("raceeth_f")) %>%
  dplyr::mutate(prev.4yr=Total/sum, prevpct.4yr=prev.4yr*100) %>%
  rename(Pop=sum, Hyst=Total) %>%
  select(-hyst_f)

#eventually used these data to do denominator adjustment
#after writing csv, used Lewis(2017) paper to calculated 95% for pooled estimates
