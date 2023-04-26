## === preparing data ===
# using meta-analyses for TBM in kids to get numbers
library(here)
library(glue)
library(data.table)
library(ggplot2)
library(ggthemes)
library(scales)
library(countrycode)
library(readxl)
id <- 'modelling/indata/'
dd <- 'modelling/data/'
xd <- 'modelling/plots/'
md <- 'metaanalyses/data/'
mp <- 'metaanalyses/plots/'

ssum <- function(x) sqrt(sum(x^2)) ## summaries
gh <- function(x) glue(here(x))
S <- function(...)sum(...,na.rm=TRUE)
## S(1,2,3)
## S(1,NA,3)


## ========== NOTE can skip to [*] if data interim data prep done

## --- WHO data
## notification data in children by country
## WHO load last year's estimates
## https://www.who.int/teams/global-tuberculosis-programme/data

## age-stratified incidence
IS <- fread(gh('{id}TB_burden_age_sex_2023-03-29.csv'))
IS[,se:=(hi-lo)/3.92]
IS <- IS[age_group %in% c('0-4','5-14'),
         .(inc.age=sum(best),inc.age.sd=ssum(se)),
         by=.(iso3,age_group)] #sum over sex
IS[,age_group:=gsub('-','_',age_group)]
ESA <- copy(IS)

## notifications
tb <- fread(gh('{id}TB_notifications_2023-03-29.csv'))
nnmz <- c('n.0.4','n.5.9','n.10.14','n.5.14','n.0.14')

TBN <- tb[year==max(year)]
TBN[,(nnmz):=.(S(newrel_m04,newrel_f04,newrel_sexunk04),
               S(newrel_m59,newrel_f59),
               S(newrel_m1014,newrel_f1014),
               S(newrel_m514,newrel_f514,newrel_sexunk514),
               S(newrel_m014,newrel_f014,newrel_sexunk014)
               ),
    by=iso3]

nnmz1 <- c('iso3',nnmz)
TBN <- TBN[,..nnmz1]

## merge with incidence
nnmz1 <- c('iso3','n.0.4','n.5.14')
TBN <- TBN[,..nnmz1]
names(TBN)[2:3] <- c('0_4','5_14')
TBN <- melt(TBN,id='iso3')
names(TBN)[2] <- 'age_group'
names(TBN)[3] <- 'notes'
TBN <- merge(TBN,ESA,by=c('iso3','age_group'))

## CDR by country
TBN[,cdr.mn:=(notes+1)/(inc.age+1)]
TBN[,cdr.vr:=(cdr.mn*inc.age.sd/(1+inc.age))^2]
TBN[,cdr.mu:=log(cdr.mn^2/sqrt(cdr.mn^2+cdr.vr))]
TBN[,cdr.sg:=sqrt(log(1+cdr.vr/cdr.mn^2))]
sg.m <- TBN[,.(fix=mean(cdr.sg)),by=age_group]
TBN <- merge(TBN,sg.m,by='age_group',all.x=TRUE)
TBN[cdr.sg==0,cdr.sg:=fix]
TBN[,fix:=NULL]

TBN[,cdr.ab:=cdr.mn*(1-cdr.mn)/(cdr.mn*inc.age.sd/(1+inc.age))^2-1] #v=k/U,dv=k.dU/U^2=v.dU/U
TBN[,c('cdr.a','cdr.b'):= .(cdr.mn * cdr.ab, (1-cdr.mn) * cdr.ab)]
mn <- TBN[is.finite(cdr.a+cdr.b),.(ab=mean(cdr.ab)),by=age_group]
mn2 <- TBN[is.finite(cdr.a+cdr.b),.(mn=sum(notes)/sum(inc.age)),by=age_group]
mn <- merge(mn,mn2,by='age_group')
mn[,c('a','b'):= .(mn * ab, (1-mn) * ab)]
TBN <- merge(TBN,mn,by='age_group',all.x=TRUE)
TBN[!is.finite(cdr.a+cdr.b),c('cdr.a','cdr.b'):=.(a,b)]
TBN[,c('a','b'):=NULL]
TBN[,c('cdr.ab','ab','mn'):=NULL]

save(TBN,file=gh('{dd}TBN.Rdata'))

## BCG by country
## from http://www.who.int/entity/immunization/monitoring_surveillance/data/coverage_estimates_series.xls?ua=1
BCG <- read_xls(gh('{id}coverage_estimates_series.xls'),sheet="BCG")
BCG <- as.data.table(BCG)[Vaccine=='BCG',.(iso3=ISO_code,bcg=`2019`)]

save(BCG,file=gh('{dd}BCG.Rdata'))


## FOI estimates by country
P <- fread(gh('{id}IHME-GBD_2019_DATA-043d7d7b-1.csv'))
## make key
key <- data.table(ihme=P[,unique(location_name)])
key[,iso3:=countrycode(ihme,origin = 'country.name',destination = 'iso3c')]
P <- merge(P,key,by.x = 'location_name',by.y = 'ihme')
P[,val.sd:=(upper-lower)/3.92]
P[iso3=='ZAF',val*1e5] #not actually percent
P <- P[,.(prev=1e5*sum(val),prev.sd=1e5*ssum(val.sd)),by=iso3]

save(P,file=gh('{dd}P.Rdata'))


## IHME populations 
POP <- fread(gh('{id}IHME_GBD_2019_POP_SYA_2019_Y2021M01D28.CSV'))
POP <- merge(POP,key,by.x='location_name',by.y='ihme',all.x=FALSE)
lagz <- c('<1 year','1','2','3','4','5','6','7','8','9','10','11','12','13','14')
lagg <- c('<1',rep('1-4',4),rep('5-9',5),rep('10-14',5))
akey <- data.table(age_group_name=lagz,acat=lagg)

POP <- POP[sex_name=='both']
POP <- POP[age_group_name %in% lagz]
POP <- merge(POP,akey,by='age_group_name',all.x = TRUE)
POP <- POP[,.(pop=sum(val)),by=.(iso3,acat)]

save(POP,file=gh('{dd}POP.Rdata'))


## ============================================ [*]
## --- join and output

## --- data from meta-analyses
## proportion of TBM in notifications by age
load(gh('{mp}BS.Rdata')) #
BS

## proportion of notifications in each age group
load(gh('{mp}AP.Rdata')) #
AP

## progression
load(gh('{mp}PD.Rdata')) #
PD

## age-dependent ORs
load(gh('{mp}cfrOR.Rdata'))
cfrOR

## HIV-related ORs
load(gh('{mp}HIVdOR.Rdata'))
load(gh('{mp}HIViOR.Rdata'))

## HIV in TB
load(file=gh('{md}hivintb.Rdata'))
hivintb <- hivintb[year==2019]
hivintb[,h.sd:=(h.hi-h.lo)/3.92]

## ---outputs prepared above
load(file=gh('{dd}TBN.Rdata'))
load(file=gh('{dd}P.Rdata'))
load(file=gh('{dd}POP.Rdata'))

## country list
length(TBN[,unique(iso3)]) #215
length(POP[,unique(iso3)]) #204
length(P[,unique(iso3)]) #204
isoz <- P[,unique(iso3)]
isoz <- intersect(isoz,TBN[,unique(iso3)])
length(isoz) #202

## master list
isoz <- sort(isoz)
cat(isoz,file=gh('{dd}isoz.txt'))

## notes & CDR
TBNW <- dcast(data=TBN,
              formula = iso3 ~ age_group,
              value.var = c('notes','cdr.mu','cdr.sg'))

## NOTE check beta version
setkey(TBNW,iso3)

## bcg
setkey(BCG,iso3)
BCG <- BCG[isoz]
bcg.mean <- mean(BCG$bcg,na.rm = TRUE)
BCG[is.na(bcg),bcg:=bcg.mean]

## pop
POPW <- dcast(POP,iso3~acat,value.var='pop')
POPW <- POPW[,.(iso3,`<1`,`1-4`,`5-9`,`10-14`)]
setkey(POPW,iso3)
sag <- unique(lagg)

## prev
setkey(P,iso3)

## HIV prev
setkey(hivintb,iso3)
hivintb <- hivintb[isoz,.(iso3,h.mid,h.sd)]
hivintb[,lmn:=log(h.mid/sqrt(1+h.sd^2/h.mid^2))]
hivintb[,lsg:=log(1+h.sd^2/h.mid^2)]
hivintb[!is.finite(lmn),lmn:=-5]
hivintb[!is.finite(lsg),lsg:=0.5]
hivintb[lsg<1e-3,lsg:=0.1]
## sdata$hiv_lmn <- hivintb$lmn
## sdata$hiv_lsg <- hivintb$lsg


## stan data
sdata <- list(
  N = length(isoz),
  CDR_mu = as.matrix(TBNW[isoz][,.(cdr.mu_0_4,cdr.mu_0_4,cdr.mu_5_14,cdr.mu_5_14)]),
  CDR_sg = as.matrix(TBNW[isoz][,.(cdr.sg_0_4,cdr.sg_0_4,cdr.sg_5_14,cdr.sg_5_14)]),
  notif_data = as.matrix(TBNW[isoz][,.(notes_0_4,notes_5_14)]),
  BCG_coverage = BCG,
  POPS = as.matrix(POPW[isoz][,..sag]),
  PREVM = P[isoz][,prev], PREVS = P[isoz][,prev.sd],
  hiv_lmn = hivintb$lmn,hiv_lsg = hivintb$lsg
)

## gamma: var=shape/rate^2, mean=shape/rate (sd=sqrt(shape)/rate)
## shape=(mean/sd)^2, rate = mean/sd^2
## sdata$PREVk <- (sdata$PREVM/sdata$PREVS)^2
## sdata$PREVr <- sdata$PREVM/(sdata$PREVS)^2

## summary(sdata$PREVk)
## summary(sdata$PREVr)

## LN version
sdata$PREVmu <- log(sdata$PREVM^2/sqrt(sdata$PREVM^2+sdata$PREVS^2))
sdata$PREVsg <- log(1+(sdata$PREVS/sdata$PREVM)^2)
summary(sdata$PREVmu)
summary(sdata$PREVsg)

sdata$PREVM <- sdata$PREVS <- NULL


## beta risk
sdata$mu_styblo <- 1.677891
sdata$sig_styblo <- 0.3714703

## proportion of note TBM
sdata$mu_notes <- BS$pred
sdata$sig_notes <- BS$se

## progression risk to TBM
## sdata$mu_prog <- BD$pred.l[1:4]
## sdata$sig_prog <- BD$pred.sdl[1:4]

## ## just cohorts:
plmn <- log(PD$prog.mn/sqrt(1+PD$prog.sd^2/PD$prog.mn^2))
plsg <- log(1+PD$prog.sd^2/PD$prog.mn^2)
## summary(rlnorm(1e4,plmn[1],plsg[1])) #check
sdata$mu_prog <- plmn
sdata$sig_prog <- plsg

## proportion of note in each age
sdata$mu_age <- AP$pred
sdata$sig_age <- AP$se

## BCG
## ## Marais
## sdata$bcgProtA <- 1.25
## sdata$bcgProtB <- 2.5
## mean(rbeta(1e4,1.25,2.5)) #HR

## Bourdin Trunz: 73 67-79
pab <- HEdtree::getAB(0.73,(0.79-0.67)^2/3.92^2)
## NOTE swap for HR
sdata$bcgProtA <- pab$b
sdata$bcgProtB <- pab$a
mean(rbeta(1e4,sdata$bcgProtA,sdata$bcgProtB)) #HR

## correct coverate
sdata$BCG_coverage <- sdata$BCG_coverage$bcg


## ## check beta parms
## summary(sdata$CDR_A)
## summary(sdata$CDR_B)

## ## u5
## bad <- (rowSums(sdata$CDR_A[,1:2]<0) + rowSums(sdata$CDR_B[,1:2]<0)) > 0
## av.a <- sum(sdata$CDR_A[!bad,1])/sum(!bad)
## av.b <- sum(sdata$CDR_B[!bad,1])/sum(!bad)
## sdata$CDR_A[bad,1:2] <- av.a
## sdata$CDR_B[bad,1:2] <- av.b

## ## o5
## bad <- (rowSums(sdata$CDR_A[,3:4]<0) + rowSums(sdata$CDR_B[,3:4]<0)) > 0
## av.a <- sum(sdata$CDR_A[!bad,3])/sum(!bad)
## av.b <- sum(sdata$CDR_B[!bad,3])/sum(!bad)
## sdata$CDR_A[bad,3:4] <- av.a
## sdata$CDR_B[bad,3:4] <- av.b

## notification safety
sdata$notif_data[sdata$notif_data[,1]==0,1] <- 1
sdata$notif_data[sdata$notif_data[,2]==0,2] <- 1


sdata$couple <- 1


## region key
regkey <- unique(tb[,.(iso3,g.whoregion=g_whoregion)])
setkey(regkey,iso3)
regkey <- regkey[isoz]

test <- model.matrix(iso3 ~ 0 + g.whoregion,data=regkey)
rkey <- as.matrix(test)
head(rkey,n=10)

for(nm in colnames(rkey)){
  nnm <- gsub("g\\.","",nm)
  sdata[[nnm]] <- which(rkey[,nm]==1)
  ## sdata[[nm]] <- NULL
}

RN <- rep(0,6)
for(i in 1:6) RN[i] <- sum(rkey[,i])

sdata$RN <- RN

names(sdata)

## Silvia
## 19·3% (95% CI 14·0–26·1)
bab <- HEdtree::getAB(0.193,(14.0-26.1)^2/392^2)

CFRtxA <- rep(bab$a,4)
CFRtxB <- rep(bab$b,4)

sdata$CFRtxA <- CFRtxA
sdata$CFRtxB <- CFRtxB


## age-dependent ORs
sdata$lgOR <- cfrOR$lgOR
sdata$lgORsd <- cfrOR$lgOR.sd


## sdata$mu_age
## sdata$sig_age

## sequelae from Silvia review
## 53·9 (42·6–64·9)
tmp <- HEdtree::getAB(0.539,(64.9-42.6)^2/392^2)
rtmp <- rbeta(1e4,tmp$a,tmp$b)
summary(rtmp)
## hist(rtmp)
sdata$seqA <- rep(tmp$a,4)
sdata$seqB <- rep(tmp$b,4)

## wide notif data
sdata$notif_dataW <- cbind(sdata$notif_data[,c(1,1)],sdata$notif_data[,c(2,2)])


## -- HIV ORs

## incidence
HIViOR[,S:=(pi.ub-pi.lb)/3.92]
mn <- HIViOR[,log(pred/sqrt(1+S^2/pred^2))]
sg <- HIViOR[,log(1+S^2/pred^2)]
## summary(rlnorm(1e4,mn,sg)) #check
sdata$mu_hivi <- mn
sdata$sig_hivi <- sg

## mortality
HIVdOR[,S:=(pi.ub-pi.lb)/3.92]
mn <- HIVdOR[,log(pred/sqrt(1+S^2/pred^2))]
sg <- HIVdOR[,log(1+S^2/pred^2)]
## summary(rlnorm(1e4,mn,sg)) #check
sdata$mu_hivd <- mn
sdata$sig_hivd <- sg


save(sdata,file=gh('{dd}sdata.Rdata'))

## TODO
## consider beta?
## tidy
