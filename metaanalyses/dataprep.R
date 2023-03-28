## looking at the data from countries in paediatric age groups
library(here)
library(ggplot2)
library(scales)
library(data.table)
library(glue)
gh <- function(x) glue(here(x))

wd <- 'metaanalyses/'      #working directory
dd <- 'metaanalyses/data/' #data directory
xd <- 'metaanalyses/plots/explore/' #exploratory plot directory

## === Brazil

## data
fn <- gh('{dd}BRA/BRA.csv')
if(!file.exists(fn)){
  ## read in
  hnk <- fread(gh('{dd}BRA/HnTBM.csv'))
  hnn <- fread(gh('{dd}BRA/HnTB.csv'))
  hpk <- fread(gh('{dd}BRA/HpTBM.csv'))
  hpn <- fread(gh('{dd}BRA/HpTB.csv'))
  names(hnk)[2] <- names(hnn)[2] <- names(hpk)[2] <- names(hpn)[2] <- '<1'
  ## reshape
  ## HIV-
  hnk <- melt(hnk[,.(year,`<1`,`1-4`,`5-9`,`10-14`)],id='year')
  names(hnk) <- c('year','age','TBM')
  hnn <- melt(hnn[,.(year,`<1`,`1-4`,`5-9`,`10-14`)],id='year')
  names(hnn) <- c('year','age','TB')
  hn <- merge(hnk,hnn,by=c('year','age'))
  hn[,hiv:='HIV-ve']
  ## HIV+
  hpk <- melt(hpk[,.(year,`<1`,`1-4`,`5-9`,`10-14`)],id='year')
  names(hpk) <- c('year','age','TBM')
  hpn <- melt(hpn[,.(year,`<1`,`1-4`,`5-9`,`10-14`)],id='year')
  names(hpn) <- c('year','age','TB')
  hp <- merge(hpk,hpn,by=c('year','age'))
  hp[,hiv:='HIV+ve']
  ## join
  BRA <- rbind(hn,hp)
  ## write out
  fwrite(BRA,file=fn)
} else {
    BRA <- fread(fn)
}
BRA$age <- factor(BRA$age,levels=c('<1','1-4','5-9','10-14'),ordered=TRUE)

## plot
ggplot(BRA,aes(age,TBM/TB)) +
    geom_boxplot(outlier.alpha = 0,notch=FALSE,alpha=0,col='grey')+
    geom_point(shape=1,aes(size=TB)) +
    scale_y_continuous(label=percent)+
    facet_wrap(~hiv)+
    theme_light()+
    ggtitle('BRA')+
    guides(shape="none",alpha="none",fill="none",col="none")+
    theme(legend.position = 'top')

ggsave(gh('{xd}BRA.png'),w=10,h=7)

## === USA

## data
fn <- gh('{dd}USA/USA.csv')
if(!file.exists(fn)){
  ## read in
  US <- fread(gh('{dd}USA/combined.csv'))
  US[,age:="<1"]
  US[grepl('1 ',Age),age:='1']
  US[grepl('2 ',Age),age:='2-4']
  US[grepl('5 ',Age),age:='5-9']
  US[grepl('10 ',Age),age:='10-14']
  unique(US[,.(age,Age)])
  US[,Age:=NULL]
  US <- melt(US,id=c('year','age'))
  US[,sex:='M']
  US[grepl('Females',variable),sex:='F']
  US[,type:='TB']
  US[grepl('gitis',variable),type:='TBM']
  US <- dcast(US,year+age+sex~type,value.var='value')
  ## write out
  fwrite(US,file=fn)
} else {
  USA <- fread(fn)
}
USA$age <- factor(USA$age,levels=c('<1','1','2-4','5-9','10-14'),ordered=TRUE)


## plot
ggplot(USA,aes(age,TBM/TB)) +
    geom_boxplot(outlier.alpha = 0,notch=FALSE,alpha=0,col='grey')+
    geom_point(shape=1,aes(size=TB)) +
    scale_y_continuous(label=percent)+
    facet_wrap(~sex)+
    theme_light()+
    ggtitle('USA')+
    guides(shape="none",alpha="none",fill="none",col="none")+
    theme(legend.position = 'top')

ggsave(gh('{xd}USA.png'),w=10,h=7)

## === ZAF
fn <- gh('{dd}ZAF/ZAF.csv')
if(!file.exists(fn)){
  ## read in
  ZA <- fread(gh('{dd}ZAF/cases.csv'))
  ZA[,sex:=ifelse(sex=='Female','F','M')]
  ZA[,hiv:=hiv2]
  ZA[,c('hiv2','propTBM'):=NULL]
  ZA[,TB:=TBcase]
  ZA[,TBcase:=NULL]
  ## write out
  fwrite(ZA,file=fn)
} else {
  ZAF <- fread(fn)
}
ZAF$age <- factor(ZAF$age,levels=c(paste0(0:14)),ordered=TRUE)


## plot
ggplot(ZAF,aes(age,TBM/TB)) +
    geom_boxplot(outlier.alpha = 0,notch=FALSE,alpha=0,col='grey')+
    geom_point(shape=1,aes(size=TB)) +
    scale_y_continuous(label=percent)+
    facet_grid(sex~hiv)+
    theme_light()+
    ggtitle('ZAF')+
    guides(shape="none",alpha="none",fill="none",col="none")+
    theme(legend.position = 'top')

ggsave(gh('{xd}ZAF.png'),w=10,h=10)


## === UKR

## data
fn <- gh('{dd}UKR/tbm_tables.csv')
if(!file.exists(fn)){
  ## read in
  UKR <- fread(fn)
  UKR[,age:="<1"]
  UKR[grepl('<2',AgeCat),age:='1']
  UKR[grepl('<5 ',AgeCat),age:='2-4']
  UKR[grepl('<10 ',AgeCat),age:='5-9']
  UKR[grepl('<15 ',AgeCat),age:='10-14']
  unique(UKR[,.(age,AgeCat)])
  UKR[,AgeCat:=NULL]
  UKR <- melt(UKR,id=c('year','age'))
  UKR[,sex:='M']
  UKR[grepl('Female',variable),sex:='F']
  UKR[,type:='TB']
  UKR[grepl('TBM',variable),type:='TBM']
  UKR[grepl('deaths',variable),type:='TBM deaths']
  UKR <- UKR[type != 'TBM deaths'] #NOTE ditch these for now
  UKR <- dcast(UKR,year+age+sex~type,value.var='value')
  ## write out
  fwrite(UKR,file=fn)
} else {
  UKR <- fread(fn)
}
UKR$age <- factor(UKR$age,levels=c('<1','1','2-4','5-9','10-14'),
                  ordered=TRUE)


## === comparison
## aggregate over ages
BRAs <- BRA
USAs <- USA
ZAFs <- ZAF
UKRs <- UKR

BRAs[,acat:=age]
BRAs[,age:=NULL]
BRAs[,sex:='unknown']
USAs[,acat:=age]
USAs[age=='1',acat:='1-4']
USAs[age=='2-4',acat:='1-4']
USAs <- USAs[,.(TB=sum(TB),TBM=sum(TBM)),by=.(year,sex,acat)]
USAs[,hiv:="HIV unknown"]

ZAFs[,acat:="<1"]
ZAFs[age %in% paste0(1:4),acat:="1-4"]
ZAFs[age %in% paste0(5:9),acat:="5-9"]
ZAFs[age %in% paste0(10:14),acat:="10-14"]
ZAFs[,age:=NULL]

UKRs[,acat:=age]
UKRs[age=='1',acat:='1-4']
UKRs[age=='2-4',acat:='1-4']
UKRs <- UKRs[,.(TB=sum(TB),TBM=sum(TBM)),by=.(year,sex,acat)]
UKRs[,hiv:="HIV-"] #NOTE assumption


## join
BRAs[,iso3:='BRA']
USAs[,iso3:='USA']
ZAFs[,iso3:='ZAF']
UKRs[,iso3:='UKR']


BRAs$acat <- factor(BRAs$acat,levels=c('<1','1-4','5-9','10-14'),ordered=TRUE)
USAs$acat <- factor(USAs$acat,levels=c('<1','1-4','5-9','10-14'),ordered=TRUE)
ZAFs$acat <- factor(ZAFs$acat,levels=c('<1','1-4','5-9','10-14'),ordered=TRUE)
UKRs$acat <- factor(UKRs$acat,levels=c('<1','1-4','5-9','10-14'),ordered=TRUE)


ALL <- rbindlist(list(BRAs,USAs,UKRs,ZAFs),use.names = TRUE)
fwrite(ALL,file=gh('{dd}TBMinTB_notECDC.csv'))


ggplot(ALL,aes(iso3,TBM/TB)) +
    geom_boxplot(outlier.alpha = 0,notch=FALSE,alpha=0,col='grey')+
    geom_point(shape=1,aes(size=TB)) +
    scale_y_continuous(label=percent)+
    facet_grid(acat~hiv+sex)+
    theme_light()+
    ggtitle('COMPARE')+
    guides(shape="none",alpha="none",fill="none",col="none")+
    theme(legend.position = 'top')


ggplot(ALL,aes(acat,TBM/TB,fill=iso3)) +
    geom_boxplot(outlier.alpha = 0,notch=FALSE,alpha=0.5)+
    ## geom_point(shape=1,aes(size=TB)) +
    scale_y_continuous(label=percent)+
    facet_grid(hiv~sex)+
    theme_light()+
    ggtitle('COMPARE')+
    theme(legend.position = 'top')

ALL2 <- ALL[hiv!='HIV+']
ALL2 <- ALL2[,.(TB=sum(TB),TBM=sum(TBM)),by=.(iso3,year,sex,acat)]


ggplot(ALL2,aes(acat,TBM/TB,fill=iso3)) +
    geom_boxplot(outlier.alpha = 0,notch=FALSE,alpha=0.5)+
    scale_y_continuous(label=percent)+
    facet_wrap(~sex)+
    theme_light()+
    ggtitle('COMPARE HIV- (NB assuming USA/UKR all HIV-ve, including ZAF HIV unknown)')+
    theme(legend.position = 'top')

ggsave(gh('{xd}compare.png'),w=12,h=7)


## === ECDC data
EU <- fread(gh('{dd}ECDC/TUBE_Case.csv'))

EU[,.N,by=.(CountryName)]
EU[,range(DateUsedForStatisticsYear)] #2010-2019

## summary
EU[,table(Outcome12Months)]
EU[,table(MajorSiteOfTB)] #376 TBM

elvls <- EU[,unique(AgeGroup)]
elvls <- elvls[c(5,1,2,4,3)]
EU$AgeGroup <- factor(EU$AgeGroup,levels=elvls,ordered=TRUE)

EUT <- EU[,.(cases=.N),by=.(AgeGroup,Gender,HIVStatus,MajorSiteOfTB)]
EUT[,sum(cases),by=HIVStatus]


ggplot(EUT,aes(AgeGroup,cases,fill=MajorSiteOfTB)) +
  facet_grid(HIVStatus~Gender)+
  geom_bar(stat='identity')+
  theme_light()

ggsave(gh('{xd}EUT.png'),w=12,h=7)

## aggregate over HIV, TBM vs not


## NOTE including minor site too
EU[,isTBM:=ifelse(MajorSiteOfTB=='MENING' | MinorSiteOfTB=='MENING',1,0)]
EUC <- EU[Gender!='UNK',
          .(TB=.N,TBM=sum(isTBM)),
          by=.(AgeGroup,sex=Gender,CountryName,
               year=DateUsedForStatisticsYear)]
EUC[,summary(TBM)]
EUC <- EU[Gender!='UNK' & HIVStatus!='POS',.(TB=.N,TBM=sum(isTBM)),by=.(AgeGroup,sex=Gender,CountryName)] #aggregate over years
EUC[,summary(TBM)]

EUC[,unique(CountryName)]
EUC[CountryName=='United Kingdom']
EU[CountryName=='United Kingdom'] #odd?

## NOTE exclusions following discussion
EUC <- EUC[!CountryName %in% c('United Kingdom','Finland')]

EUZC <- EUC[,.(TBM=sum(TBM)),by=.(CountryName)]
EUZC[,zerocountry:=ifelse(TBM==0,TRUE,FALSE)]

EUCT <- EUC[,.(TBM=sum(TBM),TB=sum(TB)),by=.(sex,AgeGroup)]
EUCT[,.(sum(TBM),sum(TB))] # 367 26829

EUC2 <- merge(EUC,EUZC[,.(CountryName,zerocountry)],by='CountryName')
EUCT2 <- EUC2[zerocountry!=TRUE,.(TBM=sum(TBM),TB=sum(TB)),by=.(sex,AgeGroup)] #only including those reporting TBM


## looking by country
EUCC <- EUC[,.(TBM=sum(TBM),TB=sum(TB)),by=.(CountryName)]
EUCC <- EUCC[order(TB)]
EUCC$CountryName <- factor(EUCC$CountryName,levels=EUCC$CountryName,ordered=TRUE)
EUCC[,notTBM:=TB-TBM]
EUCCM <- melt(EUCC[,.(CountryName,TBM,notTBM)],id='CountryName')
EUCCM[,type:=factor(variable,levels=c('notTBM','TBM'),ordered=TRUE)]
EUCC[,type:=NA]
## EUCCM <- merge(EUCCM,EUCC[,.(CountryName,TBtot=TB)],by='CountryName')

ggplot(EUCCM,aes(CountryName,value,fill=type))+
  geom_bar(stat='identity')+
  geom_text(data=EUCC,aes(CountryName,y=TB+200,label=TBM),col=2)+
  scale_y_continuous(label=comma)+
  coord_flip()+
  theme_light()+
  ylab('Total TB cases aged <15 years, 2010-2019')+
  xlab('Country')+
  ggtitle('(TBM count)')+
  theme(legend.position = 'top',plot.title = element_text(colour = "red"))
ggsave(file=gh('{xd}EUcountrybar.png'),w=6,h=7)


## --- Treat EU as single country:
agz <- c('<1','1-4','5-9','10-14')
EUY <- EU[Gender!='UNK' & HIVStatus!='POS',.(TB=.N,TBM=sum(isTBM)),
          by=.(AgeGroup,sex=Gender,year=DateUsedForStatisticsYear)] #aggregate over years

## age cats
EUY[,age:=agz[2]]
EUY[AgeGroup=='1<',age:=agz[1]]
EUY[AgeGroup=='5 - <10',age:=agz[3]]
EUY[AgeGroup=='10 - <15',age:=agz[4]]
EUY[,table(AgeGroup,age)]

## plot
ggplot(EUY,aes(age,TBM/TB)) +
  geom_boxplot(outlier.alpha = 0,notch=FALSE,alpha=0,col='grey')+
  geom_point(shape=1,aes(size=TB)) +
  scale_y_continuous(label=percent)+
  facet_grid(.~sex)+
  theme_light()+
  ggtitle('EU')+
  guides(shape="none",alpha="none",fill="none",col="none")+
  theme(legend.position = 'top')
ggsave(file=gh('{xd}EUY.png'),w=10,h=7)


## --- amend all data
ALL <- fread(gh('{dd}TBMinTB_notECDC.csv'))

ALL <- rbind(ALL,
             EUY[,.(year,TBM,TB,acat=age,
                    hiv="HIV-ve & unknown",sex,iso3="ECDC")])


fwrite(ALL,file=gh('{dd}TBMinTB.csv'))
ALL$acat <- factor(ALL$acat,levels=agz,ordered=TRUE)

ggplot(ALL,aes(iso3,TBM/TB)) +
    geom_boxplot(outlier.alpha = 0,notch=FALSE,alpha=0,col='grey')+
    geom_point(shape=1,aes(size=TB)) +
    scale_y_continuous(label=percent)+
    facet_grid(acat~hiv+sex)+
    theme_light()+
    ggtitle('COMPARE')+
    guides(shape="none",alpha="none",fill="none",col="none")+
    theme(legend.position = 'top')

ALL2 <- ALL[hiv!='HIV+']
ALL2 <- ALL2[,.(TB=sum(TB),TBM=sum(TBM)),by=.(iso3,year,sex,acat)]


ggplot(ALL2,aes(acat,TBM/TB,fill=iso3)) +
    geom_boxplot(outlier.alpha = 0,notch=FALSE,alpha=0.5)+
    scale_y_continuous(label=percent)+
    facet_wrap(~sex)+
    theme_light()+
    ggtitle('COMPARE HIV- (NB assuming USA/UKR all HIV-ve, including ZAF/EU HIV unknown)')+
    theme(legend.position = 'top')

ggsave(gh('{xd}compare_all.png'),w=12,h=7)



smy <- ALL[,.(TB=sum(TB),TBM=sum(TBM),sy=min(year),ey=max(year)),by=iso3]
smy <- smy[,.(iso3,years=paste0(sy," - ",ey),TBM,TB)]

fwrite(smy,file=gh('{dd}smy.csv'))


## looking at HIV data compared against estimates in:
## https://pubmed.ncbi.nlm.nih.gov/28807188/
## from: https://github.com/petedodd/4PM/blob/master/tables/countries_HIVpmFALSE_privFALSE.csv
Kids <- fread(gh('{dd}/other/countries_HIVpmFALSE_privFALSE.csv'))

## from: https://extranet.who.int/tme/generateCSV.asp?ds=estimates
Adz <- fread(gh('{dd}/other/TB_burden_countries_2023-03-28.csv'))

Adzr <- Adz[year==2015 & iso3 %in% Kids$iso3,.(year,iso3,e_tbhiv_prct)]

CF <- merge(Adzr,Kids,by='iso3')
md <- lm(data = CF[,.(kids=HIVinTB/1e3,all=e_tbhiv_prct/100)],all ~ kids+0)
txt <- round(coef(md),3)
CF[,hiva:=e_tbhiv_prct/100]
CF[,hivk:=HIVinTB/1e3]

GP <- ggplot(CF[,.(iso3,hivk,hiva)],
             aes(x=hiva,y=hivk,label=iso3))+
  geom_smooth(method='lm',formula = y~0+x)+
  geom_point()+
  ggrepel::geom_text_repel()+
  scale_x_continuous(label = percent)+xlab('HIV in TB (all ages)')+
  scale_y_continuous(label = percent)+ylab('HIV in TB (0-14 years)')+
  annotate(geom='text',label=glue('slope = {txt}'),col='blue',x=0.3,y=0.3,size=7)+
  ggtitle('2015: Child estimates from mortality paper vs WHO all TB estimates')

ggsave(GP,file=gh('{xd}KidHIV.png'),w=7,h=7)

confint(md) #1.6,1.9
-1e2*diff(rev(confint(md)))/txt #18%


tmp <- Adz[year>2015,.(iso3,year,e_tbhiv_prct, e_tbhiv_prct_lo, e_tbhiv_prct_hi)]
tmp[!is.na(e_tbhiv_prct),summary(-1e2*(e_tbhiv_prct_lo-e_tbhiv_prct_hi)/e_tbhiv_prct)] #84% >> regression

## NOTE neglect uncertainty from regression compared to HIV estimate
fac <- 1/txt
hivintb <- tmp[,.(iso3,year,
                  h.mid=fac*e_tbhiv_prct/1e2,
                  h.lo=fac*e_tbhiv_prct_lo/1e2,
                  h.hi=fac*e_tbhiv_prct_hi/1e2)]

save(hivintb,file=gh('{dd}hivintb.Rdata'))
