## meta-analyses for TBM in kids work
library(here)
library(metafor)
library(data.table)
library(ggplot2)
library(ggthemes)
library(scales)


## add CIs
## helper functions
getCI1 <- function(x) binom.test(x[1],x[2],p=.025)$conf.int #k,N
getCI <- function(k,N) t(apply(cbind(k,N),1,getCI1))
getCI(c(5,5,5),c(10,10,10))
## function to add binomial CIs
MLH <- function(k,N) {
  if(length(k)!=length(N)) stop('k and N have different lengths!')
  k <- as.integer(k); N <- as.integer(N)
  mid <- lo <- hi <- rep(NA,length(k))
  who <- which(!is.na(k) & !is.na(N) & N>0)
  if(any(k[who]<0)){ stop('k<0!')}
  if(any(k[who]>N[who])){ stop('k>N!')}
  HL <- getCI(k[who],N[who])
  mid[who] <- k[who]/N[who]
  lo[who] <- HL[,1]; hi[who] <- HL[,2]
  list(mid*1,lo*1,hi*1)
}
MLH(c(5,5,5),c(10,10,10))


## data
D <- fread(here('TBMinTB2.csv'))
D[,unique(hiv)]

DR <- D[!hiv %in% c('HIV+','HIV+ve')]

## 2 step: countrywise age x sex
## then aggregate these

unique(DR$sex)
(acts <- unique(DR$acat))
cnz <- unique(DR$iso3)
L <- list()

for(cn in cnz){
  DL <- DR[iso3==cn]
  sx <- unique(DL$sex)
  for(s in sx){
    for(a in acts){
      ## print(c(cn,sx,a))
      tmp <- DL[sex==s & acat==a]
      md <- rma.glmm(data=tmp,measure="PLO",xi=TBM,ni=TB)
      L[[paste0(cn,'_',s,'_',a)]] <- md
    }
  }
}

L[[1]]
names(L)

save(L,file=here('L.Rdata'))

load(file=here('L.Rdata'))



M <- list()
for(i in 1:length(L)){
  M[[i]] <- list(qty=names(L)[i],
                 lgte=coef(L[[i]]),
                 lgt.se=L[[i]]$se,
                 tau2=L[[i]]$tau2)
}
M <- rbindlist(M)
M[,c('iso3','sex','acat'):=tstrsplit(qty,'_')]
M[,qty:=NULL]

save(M,file=here('M.Rdata'))

load(file=here('M.Rdata'))

acts <- M[,unique(acat)]
M$acat <- factor(M$acat,levels=acts,ordered = TRUE)

ggplot(M,aes(x=acat,y=lgte,ymin=lgte-2*lgt.se,ymax=lgte+2*lgt.se,
             col=iso3,shape=sex))+
  geom_point(position=position_dodge(width=0.5)) +
  geom_pointrange(position=position_dodge(width=0.5))+
  scale_color_colorblind()+
  theme_classic() + ggpubr::grids()+
  xlab('Age category')+ylab('Logit RE MA estimate over years')

ggsave(here('plots/MAstep1.pdf'),w=6,h=5)


## ---  as above but not including sex:

## 2 step: countrywise age
## then aggregate these

L0 <- list()
for(cn in cnz){
  DL <- DR[iso3==cn,.(TBM=sum(TBM),TB=sum(TB)),by=.(year,acat,iso3)]
  for(a in acts){
      ## print(c(cn,sx,a))
    tmp <- DL[acat==a]
    md <- rma.glmm(data=tmp,measure="PLO",xi=TBM,ni=TB)
    L0[[paste0(cn,'_',a)]] <- md
  }
}

L0[[1]]
names(L0)

save(L0,file=here('L0.Rdata'))

load(file=here('L0.Rdata'))



M0 <- list()
for(i in 1:length(L0)){
  M0[[i]] <- list(qty=names(L0)[i],
                 lgte=coef(L0[[i]]),
                 lgt.se=L0[[i]]$se,
                 tau2=L0[[i]]$tau2)
}
M0 <- rbindlist(M0)

M0[,c('iso3','acat'):=tstrsplit(qty,'_')]
M0[,qty:=NULL]

save(M0,file=here('M0.Rdata'))

load(file=here('M0.Rdata'))

acts <- M0[,unique(acat)]
M0$acat <- factor(M0$acat,levels=acts,ordered = TRUE)

ggplot(M0,aes(x=acat,y=lgte,ymin=lgte-2*lgt.se,ymax=lgte+2*lgt.se,
             col=iso3))+
  geom_point(position=position_dodge(width=0.5)) +
  geom_pointrange(position=position_dodge(width=0.5))+
  scale_color_colorblind()+
  theme_classic() + ggpubr::grids()+
  xlab('Age category')+ylab('Logit RE MA estimate over years')+
  theme(legend.title=element_blank())

ggsave(here('plots/MAstep1nosex.pdf'),w=6,h=5)



## --- look at review papers also

## consider spline model as elsewhere
## spline model
R <- fread(here('TBMinTBreviews.csv'))
R <- R[Setting!='Spain'] #exclude because this data already in ECDC
tmp <- strsplit(R$acat,split='-')

for(i in 1:nrow(R)){
  R[i,age.lo:=as.numeric(tmp[[i]][1])]
  if(length(tmp[[i]])>1){
    R[i,age.hi:=as.numeric(tmp[[i]][2])]
  } else {
    R[i,age.hi:=NA_real_]
  }
}

R
R[is.na(age.hi),age.hi:=1.0]
R[,age.mid:=(age.lo+age.hi)/2]
R[,c('mid','lo','hi'):=MLH(TBM,TB)]

## sqrt v sqrt to look like quadratic:
## sqrt_p = C + B*sqrt_a + A*a
## assume uniform distribution of ages within range
## driver =
## (1/W)\int_L^H da (C+B.sqrt(a)+A.a)
## (1/W)(C.W+(2/3).B.(H^1.5-L^1.5)+(1/2).A.(H^2-L^2))
R[,sqrtT:=(2/3)*(age.hi^1.5-age.lo^1.5)/(age.hi-age.lo)]
R[,linT:=(1/2)*(age.hi^2-age.lo^2)/(age.hi-age.lo)] #not simplifying for clarity


A <- rma.glmm(measure = "PLO",
              xi=TBM,ni=TB,data=R,
              mods=~1+sqrtT,
              ## +linT,
              slab=Author)

summary(A)
forest(A,transf = transf.ilogit,refline=NA)
## 
cz <- coef(A)

## s1 <- function(x) strsplit(x,split = "_")[[1]][2]
## s2 <- function(x) unlist(lapply(x,s1))
## P[,Paper2:=s2(Paper)]

## predictions
agz <- seq(from=0.25,to=15,by=0.1)
PP <- predict(object = A,newmods = sqrt(agz))
PP <- as.data.table(PP)
PP[,age.mid:=agz]
PP[,c('age.lo','age.hi','Author'):=NA]
PP[,c('mid','lo','hi'):=.(exp(pred),exp(ci.lb),exp(ci.ub))]
PP[,c('lo2','hi2'):=.(exp(pi.lb),exp(pi.ub))]



GP <- ggplot(R,aes(age.mid,mid,
                   ymin=lo,ymax=hi,
                   xmin=age.lo,xmax=age.hi,
                   col=Author,shape=Author))+
  geom_pointrange()+geom_errorbarh(height=0)+
  geom_point()+scale_shape_manual(values=1:8)+
  scale_y_continuous(label=scales::percent)+
  theme_classic()+ggpubr::grids()+
  theme(legend.position = c(0.8,0.6))+
  xlab('Age (years)') + ylab('TBM proportion')

GP
GP+scale_y_sqrt() + scale_x_sqrt()



## including predz
GP <- ggplot(R,aes(age.mid,mid,
                   ymin=lo,ymax=hi,
                   xmin=age.lo,xmax=age.hi,
                   col=Author,shape=Author))+
  geom_ribbon(data=PP,aes(ymin=lo,ymax=hi),col=NA,alpha=0.3)+
  geom_line(data=PP,aes(y=lo2),lty=2,col='darkgrey')+
  geom_line(data=PP,aes(y=hi2),lty=2,col='darkgrey')+
  geom_line(data=PP,size=1,col='grey41')+
  geom_pointrange()+geom_errorbarh(height=0)+
  geom_point()+scale_shape_manual(values=1:8)+
  ## scale_y_continuous(label=scales::percent)+
  ## scale_x_continuous()+
  scale_y_sqrt(label=scales::percent)+
  scale_x_sqrt()+
  theme_classic()+ggpubr::grids()+
  theme(legend.position = c(0.8,0.8))+
  guides(color=guide_legend(ncol=2,byrow = TRUE),
         shape=guide_legend(ncol=2,byrow = TRUE))+
  xlab('Age (years, square root scale)') +
  ylab('TBM proportion (square root scale)')
GP


ggsave(GP,file=here('plots/TBMinTBreviewsSpline.pdf'),w=8,h=6)


## ---

## MA for relevant age cats
macts <- c('<1','1-4','5-9','10-14')
mlo <- c(0,1,5,10)
mhi <- c(1,5,10,15)
MM <- data.table(acat=macts,age.lo=mlo,age.hi=mhi)
MM[,sqrtT:=(2/3)*(age.hi^1.5-age.lo^1.5)/(age.hi-age.lo)]


MP <- predict(object = A,newmods = MM$sqrtT)
MP <- as.data.table(MP)
MP <- cbind(MP,MM)
MP[,c('mid','lo','hi'):=.(exp(pred),exp(ci.lb),exp(ci.ub))]
MP[,c('lo2','hi2'):=.(exp(pi.lb),exp(pi.ub))]

save(MP,file=here('MP.Rdata'))
load(file=here('MP.Rdata'))


B <- rbind(M[,.(source=paste0(iso3,": ",sex),acat,lgte,lgt.se,type='Surveillance')],
           MP[,.(source='meta-analysis',acat,lgte=pred,lgt.se=se,type='Review')])
B$acat <- factor(B$acat,levels=c('<1','1-4','5-9','10-14'),ordered = TRUE)
B[,sex:='Both']
B[grepl(': M',source),sex:='Male']
B[grepl(': F',source),sex:='Female']

psn <- position_dodge(width = 0.3)
GP <- ggplot(B,aes(acat,lgte,
                   ymin=lgte-lgt.se*1.96,ymax=lgte+lgt.se*1.96,
                   col=source,shape=type,lty=sex))+
  geom_errorbar(position=psn,width=0)+geom_point(position=psn,size=2)+
  xlab('Age category')+ylab('Proportion of TB that is TBM (logit scale)')+
  scale_y_continuous()+
  theme_bw()
GP

ggsave(GP,file=here('plots/TBMinTB_CFmaAs1_sex.pdf'),w=10,h=6)
ggsave(GP,file=here('plots/png/TBMinTB_CFmaAs1_sex.png'),w=10,h=6)





## no sex version

B <- rbind(M0[,.(source=iso3,acat,lgte,lgt.se,type='Surveillance')],
           MP[,.(source='meta-analysis',acat,lgte=pred,lgt.se=se,type='Review')])
B$acat <- factor(B$acat,levels=c('<1','1-4','5-9','10-14'),ordered = TRUE)


## meta-analysis
ma <- rma(data=B,yi=lgte,sei=lgt.se,mods = ~acat-1)

## prediction
BS <- as.data.table(predict(ma))
BS[,acat:=B$acat]
BS <- unique(BS)
BS[,type:=NA]
BS$acat <- factor(BS$acat,levels=c('<1','1-4','5-9','10-14'),ordered = TRUE)


## plot
psn <- position_dodge(width = 0.3)
GP <- ggplot(B,aes(## acat,lgte,
                   ## ymin=lgte-lgt.se*1.96,ymax=lgte+lgt.se*1.96,
                   col=source,shape=type))+
  geom_errorbar(aes(acat,ymin=lgte-lgt.se*1.96,ymax=lgte+lgt.se*1.96),position=psn,width=0)+
  geom_point(aes(acat,lgte),position=psn,size=2)+
  geom_errorbar(data=BS,aes(acat,ymin=ci.lb,ymax=ci.ub),
                  alpha=0.65,col=2,size=4,width=0)+
  geom_point(data=BS,aes(acat,pred),col=2,shape=5,size=4)+
  xlab('Age category')+ylab('Proportion of TB that is TBM (logit scale)')+
  scale_y_continuous()+
  theme_bw()
GP

ggsave(GP,file=here('plots/TBMinTB_CFmaAs1_nosex.pdf'),w=10,h=6)
ggsave(GP,file=here('plots/png/TBMinTB_CFmaAs1_nosex.png'),w=10,h=6)


ilogit <- function(x)1/(1+exp(-x))

BS[,c('prop','prop.lo','prop.hi','prop.lo2','prop.hi2'):=
      .(ilogit(pred),ilogit(ci.lb),ilogit(ci.ub),ilogit(pi.lb),ilogit(pi.ub))
   ]

BS

save(BS,file=here('plots/BS.Rdata'))
fwrite(BS,file=here('plots/BS.csv'))

load(file=here('plots/BS.Rdata'))

## ====================== HIV ======================
## BRA, ECDC?, ZAF
bra <- fread(here('BRA/BRA.csv'))
zaf <- fread(here('ZAF/ZAF.csv'))
eu <- fread(here('ECDC/TUBE_Case.csv'))

eu[,table(HIVStatus)] #68 positive NOTE don't use
bra[,sum(TB),by=hiv]  #1984
zaf[,sum(TB),by=hiv]

akey <- data.table(age=c(acts[1],rep(acts[2],4),rep(acts[3],5),rep(acts[4],5)),
                   agey=0:14)

## look at raw and then OR jj
H1 <- zaf[hiv!='HIV unknown']
H1 <- merge(H1,akey[,.(age=agey,acat=age)],by = 'age')
H1 <- H1[,.(TB=sum(TB),TBM=sum(TBM),iso3='ZAF'),by=.(year,age=acat,hiv)]
bra[,iso3:='BRA']

HD <- rbind(H1,bra)
HD$age <- factor(HD$age,levels=acts,ordered=TRUE)
HD[hiv=='HIV+ve',hiv:='HIV+']
HD[hiv=='HIV-ve',hiv:='HIV-']
HD$hiv <- factor(HD$hiv,levels=c('HIV-','HIV+'),ordered=TRUE)

summary(HD) #check

save(HD,file=here('plots/HD.Rdata'))

load(file=here('plots/HD.Rdata'))

GP <- ggplot(HD,aes(age,TBM/TB,size=TB,col=hiv))+
  geom_point(shape=1)+
  facet_wrap(~iso3)+
  scale_y_continuous(label=percent)+
  theme_light()+
  xlab('Age (years)')+ylab('TBM in TB')

ggsave(GP,file=here('plots/png/hiv_p_raw.png'),w=10,h=6)

HDW <- dcast(HD,iso3+year+age~hiv,value.var = c('TB','TBM'))
HDW[,RR:= (`TBM_HIV+`/`TB_HIV+`) / (`TBM_HIV-`/`TB_HIV-`)]
HDW[,TB:=`TB_HIV+`+`TB_HIV-`]

GP <- ggplot(HDW,aes(age,RR))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(size=TB),shape=1)+
  facet_wrap(~iso3)+
  theme_light()+
  geom_hline(yintercept = 1,col=2)+
  xlab('Age (years)')+ylab('RR for TBM given HIV')

ggsave(GP,file=here('plots/png/hiv_p_rawRR.png'),w=10,h=6)


## means?
write.csv(HDW[,(sum(`TBM_HIV+`)/sum(`TB_HIV+`)) / (sum(`TBM_HIV-`)/sum(`TB_HIV-`)),by=age],
          file=here('plots/HIV_RR_byage.csv'))
HDW[,(sum(`TBM_HIV+`)/sum(`TB_HIV+`)) / (sum(`TBM_HIV-`)/sum(`TB_HIV-`))] #1.18


## ====================== DEATHS ====================
Zm <- fread(here('ZAF/cases.csv'))
Zd <- fread(here('ZAF/deaths.csv'))

## --- HIV+ve first
ZD <- merge(Zd[hiv2!='HIV unknown',.(year,age,sex,hiv=hiv2,TBMdeaths)],
            Zm[hiv2!='HIV unknown',.(year,age,sex,hiv=hiv2,TBM)],
             by=c('year','age','sex','hiv'),
             all=TRUE)


ZD <- merge(ZD,akey[,.(age=agey,acat=age)],by = 'age')
summary(ZD)

ZD[is.na(TBMdeaths)]

ZD <- ZD[!is.na(TBMdeaths),.(TBMdeaths=sum(TBMdeaths),TBM=sum(TBM),iso3='ZAF'),by=.(year,age=acat,hiv)]
ZD$age <- factor(ZD$age,levels=acts,ordered=TRUE)
ZD$hiv <- factor(ZD$hiv,levels=c('HIV-','HIV+'),ordered=TRUE)


GP <- ggplot(ZD,aes(age,TBMdeaths/TBM,size=TBM,col=hiv))+
  geom_point(shape=1)+
  scale_y_continuous(label=percent)+
  theme_light()+
  xlab('Age (years)')+ylab('TBM deaths in notified TBM')

ggsave(GP,file=here('plots/png/hiv_d_raw.png'),w=6,h=6)


ZDW <- dcast(ZD,iso3+year+age~hiv,value.var = c('TBMdeaths','TBM'))
ZDW[,RR:= (`TBMdeaths_HIV+`/`TBM_HIV+`) / (`TBMdeaths_HIV-`/`TBM_HIV-`)]
ZDW[,TBM:=`TBM_HIV+`+`TBM_HIV-`]

summary(ZDW)

GP <- ggplot(ZDW[is.finite(RR)],aes(age,RR))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(size=TBM),shape=1)+
  theme_light()+
  geom_hline(yintercept = 1,col=2)+
  xlab('Age (years)')+ylab('RR for TBM given HIV')

ggsave(GP,file=here('plots/png/hiv_d_rawRR.png'),w=6,h=6)

## means?
write.csv(ZDW[,(sum(`TBMdeaths_HIV+`)/sum(`TBM_HIV+`)) / (sum(`TBMdeaths_HIV-`)/sum(`TBM_HIV-`)),by=age],
          file=here('plots/HIVd_RR_byage.csv'))
ZDW[,(sum(`TBMdeaths_HIV+`)/sum(`TBM_HIV+`)) / (sum(`TBMdeaths_HIV-`)/sum(`TBM_HIV-`))] #1.93


## ---- HIV-ve
Zbf <- merge(Zd[,.(year,age,sex,hiv2,TBMdeaths)],
             Zm[,.(year,age,sex,hiv2,TBM)],
             by=c('year','age','sex','hiv2'),
             all=TRUE)

GP <- ggplot(Zbf,aes(age,TBMdeaths/TBM,col=sex,size=TBM))+
  geom_point(shape=1)+facet_wrap(~hiv2)+theme_light()
GP
ggsave(GP,file=here('plots/png/deathsZAFraw.png'),w=10,h=6)

## aggregate
acts <- c('<1','1-4','5-9','10-14')
akey <- data.table(age=0:14,acat=c(acts[1],rep(acts[2],4),rep(acts[3],5),rep(acts[4],5)))
Zbf <- merge(Zbf,akey,by='age')
Zbf[is.na(TBMdeaths)]

Zb <- Zbf[!is.na(TBMdeaths),.(TBMdeaths=sum(TBMdeaths),TBM=sum(TBM)),by=.(year,hiv2,acat)]
Zb$acat <- factor(Zb$acat,levels=acts,ordered = TRUE)

## overall raw
Zb[,.(k=sum(TBMdeaths),N=sum(TBM)),by=.(hiv2)] #69  1329


GP <- ggplot(Zb,aes(acat,TBMdeaths/TBM,size=TBM))+
  geom_point(shape=1)+facet_wrap(~hiv2)+theme_light()
GP

ggsave(GP,file=here('plots/png/deathsZAFaggr.png'),w=10,h=6)

GP <- ggplot(Zb[hiv2!='HIV unknown'],aes(acat,TBMdeaths/TBM,size=TBM))+
  geom_point(shape=1)+facet_wrap(~hiv2)+theme_light()
GP

ggsave(GP,file=here('plots/png/deathsZAFaggr_noUnk.png'),w=10,h=6)

## meta-analyses
F <- list()
for(a in acts){
  for(hiv in c('HIV+','HIV-')){
    tmp <- Zb[hiv2==hiv & acat==a]
    md <- rma.glmm(data=tmp,measure="PLO",xi=TBMdeaths,ni=TBM)
    F[[paste0(a,',',hiv)]] <- md
  }
}

## extraction
CFR <- list()
for(i in 1:length(F)){
  CFR[[i]] <- list(qty=names(F)[i],
                 lgte=coef(F[[i]]),
                 lgt.se=F[[i]]$se,
                 tau2=F[[i]]$tau2)
}
CFR <- rbindlist(CFR)

CFR[,c('acat','hiv'):=tstrsplit(qty,',')]
CFR[,qty:=NULL]
CFR$acat <- factor(CFR$acat,levels=acts,ordered = TRUE)

save(CFR,file=here('CFR.Rdata'))

load(file=here('CFR.Rdata'))


CFR[,c('TBM','TBMdeaths'):=1.0]
CFR[,hiv2:=hiv]

GP <- ggplot(Zb[hiv2!='HIV unknown'],aes(acat,TBMdeaths/TBM,size=TBM))+
  geom_point(shape=1)+
  facet_wrap(~hiv2)+theme_light()+
  geom_pointrange(data=CFR,size=1,
                  aes(x=acat,y=transf.ilogit(lgte),
                      ymin=transf.ilogit(lgte-2*lgt.se),
                      ymax=transf.ilogit(lgte+2*lgt.se)))+
  scale_y_continuous(label=percent)+
  xlab('Age category') + ylab('Case fatality ratio')


GP

ggsave(GP,file=here('plots/png/deathsZAFaggrMA_noUnk.png'),w=10,h=6)


ggplot(M,aes(x=acat,y=lgte,ymin=lgte-2*lgt.se,ymax=lgte+2*lgt.se,
             col=iso3,shape=sex))+
  scale_color_colorblind()+
  theme_classic() + ggpubr::grids()+
  xlab('Age category')+ylab('Logit RE MA estimate over years')

ggsave(here('plots/MAstep1.pdf'),w=6,h=5)





## Silvia:
silv <- c(19.3,14,26.1)
diff(silv[-1])/silv[1]/3.92


## new version - use the overall raw CFR for all ages as reference
## 69/1329

## Use these to compute ORs?
CFR[,lgt.se/lgte]
tmp <- CFR[hiv=='HIV-']
tmp[,lgOR:=lgte-transf.logit(69/1329)]
tmp[,lgOR.sd:=lgt.se]
tmp[,OR:=exp(lgOR)]
tmp[,OR.lo:=exp(lgOR-1.96*lgOR.sd)]
tmp[,OR.hi:=exp(lgOR+1.96*lgOR.sd)]


cfrOR <- tmp[,.(acat,lgOR,lgOR.sd,OR,OR.lo,OR.hi)]

save(cfrOR,file = here('cfrOR.Rdata'))
fwrite(cfrOR,file = here('cfrOR.csv'))



## ## Use these to compute ORs?
## CFR[,lgt.se/lgte]
## tmp <- CFR[hiv=='HIV-']
## rf <- tmp[acat=='10-14'] #oldies
## tmp[,lgOR:=lgte-rf$lgte]
## tmp[acat=='10-14',lgt.se:=0]
## tmp[,lgOR.sd:=sqrt(rf$lgt.se^2+lgt.se^2)]
## tmp[,OR:=exp(lgOR)]
## tmp[,OR.lo:=exp(lgOR-1.96*lgOR.sd)]
## tmp[,OR.hi:=exp(lgOR+1.96*lgOR.sd)]
## tmp[acat=='10-14',c('lgOR.sd','OR.lo','OR.hi'):=.(0,1,1)]

## cfrOR <- tmp[,.(acat,lgOR,lgOR.sd,OR,OR.lo,OR.hi)]

## save(cfrOR,file = here('cfrOR.Rdata'))
## fwrite(cfrOR,file = here('cfrOR.csv'))
