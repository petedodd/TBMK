library(here)
library(glue)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(scales)
library(data.table)
library(stringr)
library(rstan)
options(mc.cores = 4)

xd <- 'modelling/plots/'
dd <- 'modelling/data/'
id <- 'modelling/indata/'
zd <- 'modelling/stan/'

gh <- function(x) glue(here(x))
ssum <- function(x) sqrt(sum(x^2))
source(gh('modelling/utils.R'))

## load data for stan
load(file=gh('{dd}sdata.Rdata'))
## str(sdata)


## compile stan model
## mdl <- stan_model(file=gh('{zd}combine.stan'))
mdlH <- stan_model(file=gh('{zd}combineH.stan'))


## NOTE vectors and in [0,1]
## sample
set.seed(2345)
sdata$couple <- 1
smps <- sampling(mdlH, data=sdata, chains=4, iter=5000, cores=4)
## smps <- sampling(mdlH, data=sdata, chains=1, iter=1000, cores=1) #test

## global
A <- summary(smps,pars=c('global_mnotes','global_minc',
                         'global_deaths','global_morbs',
                         'global_fracDeathsUntreated',
                         'global_fracHIVInc','global_fracHIVDeaths'))$summary

A

save(A,file=gh('{xd}A.Rdata'))
write.csv(A,file=gh('{xd}A.csv'))

H <- summary(smps,pars=c('regional_fracHIVinc','regional_fracHIVmort',
                         'global_fracHIVincAge','global_fracHIVmortAge'))$summary

save(H,file=gh('{xd}H.Rdata'))


round(A[1:4,c('mean','2.5%','97.5%')],-2)
round(A[5:7,c('mean','2.5%','97.5%')]*100,1)

## ## NOTE checks
## ## foi,tbivec, and prog
## (tmp <- summary(smps,pars=c('foi'))$summary[,'mean'])
## summary(tmp*1e2) #range 0.2-0.4% py which is reasonable
## (tmp <- summary(smps,pars=c('BCGprot'))$summary[,'mean'])
## summary(tmp*1e2) #range 0.2-0.4% py which is reasonable
## (tmp2 <- summary(smps,pars=c('tbivec'))$summary[,'mean'])
## summary(tmp2*1e2) #also OK
## (tmp3 <- summary(smps,pars=c('prog'))$summary[,'mean'])
## summary(tmp3*1e2) #range 2-4% 
## sum(sdata$POPS)/1e9 #2 billion
## (tmp4 <- summary(smps,pars=c('propTBM'))$summary[,'mean'])
## summary(tmp4*1e2) # ~ 3% 

## regional
B1 <- summary(smps,pars=c('regional_AFR_mnotes','regional_AMR_mnotes',
                    'regional_EMR_mnotes','regional_EUR_mnotes',
                    'regional_SEA_mnotes','regional_WPR_mnotes'))$summary[,c('mean','sd')]
B2 <- summary(smps,pars=c('regional_AFR_minc','regional_AMR_minc',
                    'regional_EMR_minc','regional_EUR_minc',
                    'regional_SEA_minc','regional_WPR_minc'))$summary[,c('mean','sd')]
B3 <- summary(smps,pars=c('regional_AFR_deaths','regional_AMR_deaths',
                    'regional_EMR_deaths','regional_EUR_deaths',
                    'regional_SEA_deaths','regional_WPR_deaths'))$summary[,c('mean','sd')]
B4 <- summary(smps,pars=c('regional_AFR_morbs','regional_AMR_morbs',
                          'regional_EMR_morbs','regional_EUR_morbs',
                          'regional_SEA_morbs','regional_WPR_morbs'))$summary[,c('mean','sd')]


## age summaries
agep <- c('age_u1_minc',
          'age_1to4_minc',
          'age_5to9_minc',
          'age_10to14_minc',
          'age_u1_deaths',
          'age_1to4_deaths',
          'age_5to9_deaths',
          'age_10to14_deaths',
          'age_u1_morbs',
          'age_1to4_morbs',
          'age_5to9_morbs',
          'age_10to14_morbs'
          )
ageout <- summary(smps,pars=agep)$summary[,c('mean','sd')]

## global/regional table
A <- mkdt(A)
B1 <- mkdt(B1)
B2 <- mkdt(B2)
B3 <- mkdt(B3)
B4 <- mkdt(B4)
A <- rbindlist(list(A,B1,B2,B3,B4))
A[,c('type','region','quantity'):=tstrsplit(variable,split='_')]
A[is.na(quantity),quantity:=region]
A[type=='global',region:='global']

AW <- dcast(A,type+region~quantity,value.var = c('mean','sd'))
setcolorder(AW,c('type','region',
                 'mean_minc','sd_minc',
                 'mean_mnotes','sd_mnotes',
                 'mean_deaths','sd_deaths',
                 'mean_morbs','sd_morbs'))
AW

save(AW,file=gh('{xd}AW4.Rdata'))
fwrite(AW,file=gh('{xd}AW4.csv'))

## age split table
atab <- mkdt(ageout)
acts <- c('<1','1-4','5-9','10-14')

atab[,acat:=rep(acts,3)]
atab[,quantity:=rep(c('incidence','deaths','morbs'),each=4)]

AWA <- dcast(atab,acat~quantity,value.var = c('mean','sd'))
setcolorder(AWA,c('acat',
                 'mean_incidence','sd_incidence',
                 'mean_deaths','sd_deaths',
                 'mean_morbs','sd_morbs'))

AWA$acat <- factor(AWA$acat,levels=acts,ordered=TRUE)
setkey(AWA,acat)

AWA

save(AWA,file=gh('{xd}AWA3.Rdata'))
fwrite(AWA,file=gh('{xd}AWA3.csv'))



out <- summary(smps)$summary
nmz <- rownames(out)
out <- as.data.table(out)
out[,variable:=nmz]

out[Rhat>1.1,variable]
out

save(out,file=gh('{xd}out.Rdata'))


notes <- out[grepl('TBMnotes',variable)]
incs <- out[grepl('TBMI\\[',variable)]
incpcs <- out[grepl('TBMIpc',variable)]
cdr <- out[grepl('NoI',variable)]
tbmp <- out[grepl('propTBM',variable)]

## NOTE check ball park with means
incs[,sum(mean)]/1e3 #81K without notifications? about 31 with
summary(cdr)
summary(tbmp)


## ============ neater outputs
load(file=gh('{xd}out.Rdata'))
## out[,unique(variable)]
load(file=gh('{xd}A.Rdata'))
## A
load(file=gh('{xd}H.Rdata'))


## keep only those with a comma
outr <- out[grepl(",",variable)]

test <- head(outr$variable)

## test utils functions
getcno(test)
getAno(test)
getcn(test)
getac(test)


outr[,iso3:=getcn(variable)]
outr[,acat:=getac(variable)]


dropbrkts <- function(x) gsub("\\[(.*?)\\]","",x)

outr[,variable:=dropbrkts(variable)]
outr[,unique(variable)]

outr2 <- outr[variable %in% c('deaths','untreated','TBMI','TBMnotes','morbs')]
outr <- outr[variable %in% c('deaths','untreated','TBMI','TBMnotes')]


## f1 - inc by age, col= tx status
tmp <- outr[variable %in% c('untreated','TBMnotes'),.(value=sum(mean)),by=.(acat,variable)]
tmp[,var:=ifelse(variable=='untreated','untreated','treated')]
tmp$acat <- factor(tmp$acat,levels = acts,ordered=TRUE)

ggplot(tmp,aes(acat,value,fill=var))+
  geom_bar(stat='identity')+
  scale_fill_colorblind()+
  xlab('Age (years)')+
  ylab('TBM incidence')+
  scale_y_continuous(label=comma)+
  theme_classic()+ggpubr::grids()+
  theme(legend.position = 'top',legend.title = element_blank())

ggsave(gh('{xd}f1.png'),w=5,h=5)


## f2 - country bar chart
tmp <- outr[variable %in% c('TBMI'),.(value=sum(mean)),by=.(iso3)]
tmp <- tmp[rev(order(value))]
tmp <- tmp[1:10]
tmp$iso3 <- factor(tmp$iso3,levels=rev(tmp$iso3),ordered=TRUE)

ggplot(tmp,aes(iso3,value))+
  geom_bar(stat='identity')+
  coord_flip()+
  theme_classic()+ggpubr::grids()+
  scale_y_continuous(label=comma)+
  ylab('TBM incidence')+xlab('Country (ISO3 code)')

ggsave(gh('{xd}f2.png'),w=5,h=6)


## same as above but mortality
## f1 - inc by age, col= tx status
tmp <- outr[variable %in% c('untreated','deaths'),.(value=sum(mean)),by=.(acat,variable)]
tmp <- dcast(tmp,acat~variable,value.var = 'value')
tmp[,treated:=deaths-untreated]
tmp <- melt(tmp[,.(acat,treated,untreated)],id='acat')
tmp$acat <- factor(tmp$acat,levels = acts,ordered=TRUE)

ggplot(tmp,aes(acat,value,fill=variable))+
  geom_bar(stat='identity')+
  scale_fill_colorblind()+
  xlab('Age (years)')+
  ylab('TBM deaths')+
  scale_y_continuous(label=comma)+
  theme_classic()+ggpubr::grids()+
  theme(legend.position = 'top',legend.title = element_blank())

ggsave(gh('{xd}f3.png'),w=5,h=5)

## f2 - country bar chart
tmp <- outr[variable %in% c('deaths'),.(value=sum(mean)),by=.(iso3)]
tmp <- tmp[rev(order(value))]
tmp <- tmp[1:10]
tmp$iso3 <- factor(tmp$iso3,levels=rev(tmp$iso3),ordered=TRUE)

ggplot(tmp,aes(iso3,value))+
  geom_bar(stat='identity')+
  coord_flip()+
  theme_classic()+ggpubr::grids()+
  scale_y_continuous(label=comma)+
  ylab('TBM deaths')+xlab('Country (ISO3 code)')

ggsave(gh('{xd}f4.png'),w=5,h=6)


## table functions
## looking at tables
outr2 <- outr2[,.(mean,sd,variable,iso3,acat)]

save(outr2,file=gh('{xd}outr2.Rdata'))

load(file=gh('{xd}outr2.Rdata'))



tmp <- dcast(outr2,iso3 + acat ~ variable,value.var = c('mean','sd'))
tmp <- tmp[,.(iso3,acat,
              ## --- means
              ## inc
              mean_inc_tx=mean_TBMnotes,
              mean_inc_utx=mean_TBMI-mean_TBMnotes,
              mean_inc_tot=mean_TBMI,
              ## mort
              mean_mort_tx=mean_deaths-(mean_TBMI-mean_TBMnotes),
              mean_mort_utx=mean_TBMI-mean_TBMnotes,
              mean_mort_tot=mean_deaths,
              ## morb
              mean_morb_tx=mean_morbs,
              mean_morb_utx = 0,
              mean_morb_tot=mean_morbs,
              ## --- SDs
              ## inc
              sd_inc_tx=sd_TBMnotes,
              sd_inc_utx=sqrt(sd_TBMI^2+sd_TBMnotes^2),
              sd_inc_tot=sd_TBMI,
              ## mort
              sd_mort_tx=sqrt(sd_deaths^2+sd_TBMI^2+sd_TBMnotes^2),
              sd_mort_utx=sqrt(sd_TBMI^2+sd_TBMnotes^2),
              sd_mort_tot=sd_deaths,
              ## morb
              sd_morb_tx=sd_morbs,
              sd_morb_utx = 0,
              sd_morb_tot=sd_morbs
              )]

tmp <- melt(tmp,id=c('iso3','acat'))

tmp[,c('stat','qty','txstatus'):=tstrsplit(variable,split='_')]
tmp[,variable:=NULL]
tmp <- dcast(tmp,iso3+acat+qty+txstatus~stat,value.var = 'value')


tb <- fread(gh('{id}TB_notifications_2023-03-29.csv'))
rkey <- unique(tb[,.(iso3,g.whoregion=g_whoregion)])

tmp <- merge(tmp,rkey,by='iso3',all.x=TRUE)


## global
gtots <- tmp[,.(mid=sum(mean),s=ssum(sd)),by=.(acat,qty,txstatus)]
gtota <- gtots[,.(mid=sum(mid),s=ssum(s)),by=.(qty,txstatus)]
gtota[,acat:='all']
gtots <- rbind(gtots,gtota)
gtotsw <- dcast(gtots,qty + txstatus ~ acat,value.var = c('mid','s'))

gtotsw <- gtotsw[,.(`<1`=brkt(`mid_<1`,`s_<1`),
                    `1-4`=brkt(`mid_1-4`,`s_1-4`),
                    `5-9`=brkt(`mid_5-9`,`s_5-9`),
                    `10-14`=brkt(`mid_10-14`,`s_10-14`),
                    all=brkt(mid_all,s_all)),
                 by=.(qty,txstatus)]
gtotsw$qty <- factor(gtotsw$qty,levels=c('inc','morb','mort'),ordered = TRUE)
gtotsw$txstatus <- factor(gtotsw$txstatus,levels=c('tx','utx','tot'),ordered = TRUE)
setkey(gtotsw,qty,txstatus)
str(gtotsw)

fwrite(gtotsw,file=gh('{xd}Table1a.csv'))


## regional
rtots <- tmp[,.(mid=sum(mean),s=ssum(sd)),by=.(acat,qty,txstatus,g.whoregion)]
rtota <- rtots[,.(mid=sum(mid),s=ssum(s)),by=.(qty,txstatus,g.whoregion)]
rtota[,acat:='all']
rtota2 <- rtots[,.(mid=sum(mid),s=ssum(s)),by=.(qty,txstatus,acat)]
rtota2[,g.whoregion:='GLOBAL']
rtota3 <- rtots[,.(mid=sum(mid),s=ssum(s)),by=.(qty,txstatus)]
rtota3[,c('acat','g.whoregion'):=.('all','GLOBAL')]
rtots <- rbindlist(list(rtots,rtota,rtota2,rtota3),use.names = TRUE)
rtotsw <- dcast(rtots,qty + txstatus + g.whoregion ~ acat,value.var = c('mid','s'))

rtotsw <- rtotsw[,.(`<1`=brkt(`mid_<1`,`s_<1`),
                    `1-4`=brkt(`mid_1-4`,`s_1-4`),
                    `5-9`=brkt(`mid_5-9`,`s_5-9`),
                    `10-14`=brkt(`mid_10-14`,`s_10-14`),
                    all=brkt(mid_all,s_all)),
                 by=.(qty,txstatus,g.whoregion)]

rtotsw$qty <- factor(rtotsw$qty,levels=c('inc','morb','mort'),ordered = TRUE)
rtotsw$txstatus <- factor(rtotsw$txstatus,levels=c('tx','utx','tot'),ordered = TRUE)
rtotsw$g.whoregion <- factor(rtotsw$g.whoregion,levels=c('AFR','AMR','EMR','EUR','SEA','WPR','GLOBAL'),
                             ordered = TRUE)
setkey(rtotsw,qty,txstatus,g.whoregion)
rtotsw

fwrite(rtotsw,file=gh('{xd}Table1c.csv'))


## regional without age
atots <- tmp[,.(mid=sum(mean),s=ssum(sd)),by=.(g.whoregion,qty,txstatus)]
atota <- atots[,.(mid=sum(mid),s=ssum(s)),by=.(qty,txstatus)]
atota[,g.whoregion:='GLOBAL']
atots <- rbind(atots,atota)
atotsw <- dcast(atots,qty + txstatus ~ g.whoregion,value.var = c('mid','s'))

atotsw <- atotsw[,.(AFR=brkt(mid_AFR,s_AFR),
                    AMR=brkt(mid_AMR,s_AMR),
                    EMR=brkt(mid_EMR,s_EMR),
                    EUR=brkt(mid_EUR,s_EUR),
                    SEA=brkt(mid_SEA,s_SEA),
                    WPR=brkt(mid_WPR,s_WPR),
                    GLOBAL=brkt(mid_GLOBAL,s_GLOBAL)),
                 by=.(qty,txstatus)]
atotsw$qty <- factor(atotsw$qty,levels=c('inc','morb','mort'),ordered = TRUE)
atotsw$txstatus <- factor(atotsw$txstatus,levels=c('tx','utx','tot'),ordered = TRUE)
setkey(atotsw,qty,txstatus)
atotsw

fwrite(atotsw,file=gh('{xd}Table1b.csv'))


## bar plot
btmp <- copy(rtots[acat!='all' & txstatus!='tot'])
btmp$acat <- factor(btmp$acat,levels=c('<1','1-4','5-9','10-14'),ordered=TRUE)
btmp[,QTY:=qty]
btmp[qty=='inc',QTY:='TBM incidence']
btmp[qty=='mort',QTY:='TBM deaths']
btmp[qty=='morb',QTY:='TBM sequelae']
btmp$QTY <- factor(btmp$QTY,levels=c('TBM incidence','TBM sequelae','TBM deaths'),ordered = TRUE)
btmp[,TX:=txstatus]
btmp[txstatus=='utx',TX:='untreated for tuberculosis']
btmp[txstatus!='utx',TX:='treated for tuberculosis']
btmp$TX <- factor(btmp$TX,levels=c('untreated for tuberculosis','treated for tuberculosis'),
                  ordered = TRUE)
whoregkey <- data.table(g.whoregion=c("AFR", "AMR", "EMR", "EUR", "SEA", "WPR",'GLOBAL'),
                        Region=c('Africa','The Americas','Eastern Mediterranean',
                                 'Europe','South-East Asia','Western Pacific','GLOBAL'))
btmp <- merge(btmp,whoregkey,by='g.whoregion')
btmp$Region <- factor(btmp$Region,levels=whoregkey$Region,ordered = TRUE)

GG <- ggplot(btmp,
       aes(acat,mid,fill=TX)) +
  geom_bar(stat='identity')+
  facet_grid(QTY~Region)+
  scale_fill_colorblind(name=element_blank())+
  scale_y_continuous(label=comma)+
  ylab('Number in 2019')+
  xlab('Age category (years)')+
  theme_light()+
  theme(legend.position = 'top')

ggsave(GG,file=gh('{xd}MG.png'),w=12,h=7)

## === new version of MG with HIV
brktonly <- function(x) as.numeric(gsub("\\]","",gsub("^.+\\[","",x)))

HD <- mkdt(H)
HD[,qty:=ifelse(grepl('mort',HD$variable),'mort','inc')]
HD[,acat:=getac(variable)]
HD[is.na(acat),acat:=acts[brktonly(variable)]]
whorg <- c("AFR", "AMR", "EMR", "EUR", "SEA", "WPR")
HD[,g.whoregion:=whorg[getcno(variable)]]
HD[is.na(g.whoregion),g.whoregion:='GLOBAL']
HD[,variable:=NULL]
tmp <- HD[qty=='inc']
tmp[,qty:='morb']
HD <- rbind(HD,tmp)
H <- HD[,.(g.whoregion,acat,qty,value=mean,S=sd)]

btmp[Region=='GLOBAL'][,sum(mid),by=qty] #check

## global
htmp <- merge(btmp,H,by=c('qty','acat','g.whoregion'),all.x=TRUE)
## NOTE HIV uncertainty much smaller so neglect
htmp[,`HIV-infected`:=value*mid]
htmp[,`HIV-uninfected`:=(1-value)*mid]
htmp[,c('value','g.whoregion','mid'):=NULL]
htmp[,c('qty','s','txstatus','S'):=NULL]
htmp <- melt(htmp,id=c('QTY','TX','acat','Region'))
htmp[,combo:=paste0(variable,', ',TX)]
hlvs <- c(
  "HIV-infected, untreated for tuberculosis",
  "HIV-uninfected, untreated for tuberculosis",
  "HIV-infected, treated for tuberculosis",
  "HIV-uninfected, treated for tuberculosis"
)
htmp$combo <- factor(htmp$combo,levels=hlvs,ordered=TRUE)
htmp$acat <- factor(htmp$acat,levels=acts,ordered=TRUE)
str(htmp)

htmp[Region=='GLOBAL'][,sum(value),by=QTY] #check


GG <- ggplot(htmp,
             aes(acat,value,fill=combo)) +
  geom_bar(stat='identity')+
  facet_grid(QTY~Region)+
  scale_fill_colorblind(name=element_blank())+
  scale_y_continuous(label=comma)+
  ylab('Number in 2019')+
  xlab('Age category (years)')+
  theme_light()+
  theme(legend.position = 'top')
GG

ggsave(GG,file=gh('{xd}MGH.png'),w=12,h=7)

## Two piece version

GGR <- ggplot(htmp[Region!='GLOBAL'],
             aes(acat,value,fill=combo)) +
  geom_bar(stat='identity')+
  facet_grid(QTY~Region)+
  scale_fill_colorblind(name=element_blank())+
  scale_y_continuous(label=comma)+
  ylab('')+
  xlab('Age category (years)')+
  theme_light()+
  theme(legend.position = 'top',axis.title.x = element_text(hjust=0.4))
## GGR
GGG <- ggplot(htmp[Region=='GLOBAL'],
              aes(acat,value,fill=combo)) +
  geom_bar(stat='identity')+
  facet_grid(QTY~Region,switch='y')+
  scale_fill_colorblind(name=element_blank())+
  scale_y_continuous(label=comma)+
  ylab('Number in 2019')+
  xlab('')+
  theme_light()+
  theme(legend.position = 'top')
## GGG
w <- 1
GGG <- GGG + theme(legend.spacing.x = unit(w, 'cm'))
GGR <- GGR + theme(legend.spacing.x = unit(w, 'cm'))

GG <- ggarrange(GGG,GGR,
                widths=c(0.25,1),
                common.legend = TRUE)
ggsave(GG,file=gh('{xd}MGHS.png'),w=14,h=10)
ggsave(GG, file = gh("{xd}MGHS.pdf"), w = 14, h = 10)



## ================= supplemental tables ===========
htmp <- merge(btmp,H,by=c('qty','acat','g.whoregion'),all.x=TRUE)
## NOTE HIV uncertainty much smaller so neglect
htmp[,`HIV-infected`:=value*mid]
htmp[,`HIV-uninfected`:=(1-value)*mid]
htmp[,`Total`:=mid]
htmp[,his:= s * sqrt(value)]
htmp[,hus:= s * sqrt(1-value)]
htmp[,c('value','mid'):=NULL]
htmp[,c('qty','txstatus','S'):=NULL]
htmp[,g.whoregion:=NULL]

## make the no age but regional table:
noage <- htmp[,.(Total=sum(Total),
                 `HIV-infected`=sum(`HIV-infected`),`HIV-uninfected`=sum(`HIV-uninfected`),
                 s=ssum(s),his=ssum(his),hus=ssum(hus)
                 ),by=.(Region,QTY,TX)]

noagetab <- noage[,.(QTY,Region,TX,
                    `HIV-uninfected`=brkt(`HIV-uninfected`,hus),
                    `HIV-infected`=brkt(`HIV-infected`,his),
                    Total=brkt(Total,s)
                    )]

tmp <- melt(noagetab,id=c('QTY','Region','TX'))
noagetab <- dcast(tmp,QTY + TX + variable ~ Region,value.var='value')
names(noagetab)[names(noagetab)=='variable'] <- 'HIV'

## save:
fwrite(noagetab,file=gh('{xd}noagetab.csv'))

## make the age but no region table
noreg <- htmp[Region != 'GLOBAL',
              .(Total=sum(Total),
                `HIV-infected`=sum(`HIV-infected`),`HIV-uninfected`=sum(`HIV-uninfected`),
                s=ssum(s),his=ssum(his),hus=ssum(hus)
                ),by=.(AGE=acat,QTY,TX)]

noregtab <- noreg[,.(QTY,AGE,TX,
                     `HIV-uninfected`=brkt(`HIV-uninfected`,hus),
                     `HIV-infected`=brkt(`HIV-infected`,his),
                     Total=brkt(Total,s)
                     )]


tmp <- melt(noregtab,id=c('QTY','AGE','TX'))
noregtab <- dcast(tmp,QTY + TX + variable ~ AGE,value.var='value')

names(noregtab)[names(noregtab)=='variable'] <- 'HIV'

noregtab <- noregtab[order(QTY,TX,HIV),.(QTY,TX,HIV,`<1`,`1-4`,`5-9`,`10-14`)]

## save:
fwrite(noregtab,file=gh('{xd}noregtab.csv'))

