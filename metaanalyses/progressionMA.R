## meta-analyses for TBM progression
library(here)
library(metafor)
library(data.table)
library(ggplot2)
library(ggthemes)
library(scales)
library(glue)
wd <- 'metaanalyses/'      #working directory
dd <- 'metaanalyses/data/' #data directory
td <- 'metaanalyses/data/tmp/' #tmp data directory
xd <- 'metaanalyses/plots/' #export/plot directory

hi <- function(x)quantile(x,0.975)
lo <- function(x)quantile(x,0.025)
ilogit <- function(x)1/(1+exp(-x))
gh <- function(x) glue(here(x))

## tidying paper names
s1 <- function(x) strsplit(x,split = "_")[[1]][2]
s2 <- function(x) unlist(lapply(x,s1))

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

## main data
P <- fread(gh('{dd}TBMprogressionReview.csv'))
P[,c('mid','lo','hi'):=MLH(TBM,infected)]
P[,age.mid:=(age.lo+age.hi)/2]

## exclude ML5 - due to overlap
P[,unique(Paper)]
P <- P[Paper!="ML5_Miller"]

GP <- ggplot(P,aes(age.mid,mid,
                   ymin=lo,ymax=hi,
                   xmin=age.lo,xmax=age.hi,
                   col=Paper,shape=Paper))+
  geom_pointrange()+geom_errorbarh(height=0)+
  geom_point()+
  scale_y_continuous(label=scales::percent)+
  theme_classic()+ggpubr::grids()+
  theme(legend.position = c(0.8,0.6))+
  xlab('Age (years)') + ylab('TBM risk')

GP

GP+scale_y_sqrt() + scale_x_sqrt()

## sqrt v sqrt to look like quadratic:
## sqrt_p = C + B*sqrt_a + A*a
## assume uniform distribution of ages within range
## driver =
## (1/W)\int_L^H da (C+B.sqrt(a)+A.a)
## (1/W)(C.W+(2/3).B.(H^1.5-L^1.5)+(1/2).A.(H^2-L^2))
P[,sqrtT:=(2/3)*(age.hi^1.5-age.lo^1.5)/(age.hi-age.lo)]
P[,linT:=(1/2)*(age.hi^2-age.lo^2)/(age.hi-age.lo)] #not simplifying for clarity


A <- rma.glmm(measure = "PLO",
              xi=TBM,ni=infected,data=P,
              mods=~1+sqrtT,
              ## +linT,
              slab=Paper)

summary(A)
GP+scale_y_sqrt() + scale_x_sqrt()
cz <- coef(A)

## tidy paper names
P[,Paper2:=s2(Paper)]
P[,Paper3:=Paper2]
P[Paper2=='Cam&Mil',Paper3:='Cammock']

## predictions
agz <- seq(from=0.25,to=15,by=0.1)
PP <- predict(object = A,newmods = sqrt(agz))
PP <- as.data.table(PP)
PP[,age.mid:=agz]
PP[,c('age.lo','age.hi','Paper','Paper2','Paper3'):=NA]
PP[,c('mid','lo','hi'):=.(exp(pred),exp(ci.lb),exp(ci.ub))]
PP[,c('lo2','hi2'):=.(exp(pi.lb),exp(pi.ub))]

brks <- c(0,1,5,10,15)

## graph
GP <- ggplot(P,aes(age.mid,mid,
                   ymin=lo,ymax=hi,
                   xmin=age.lo,xmax=age.hi,
                   col=Paper3,shape=Paper3))+
  geom_ribbon(data=PP,aes(ymin=lo,ymax=hi),col=NA,alpha=0.3)+
  geom_line(data=PP,aes(y=lo2),lty=2,col='darkgrey')+
  geom_line(data=PP,aes(y=hi2),lty=2,col='darkgrey')+
  geom_line(data=PP,size=1,col='grey41')+
  geom_pointrange()+geom_errorbarh(height=0)+
  geom_point()+
  scale_y_sqrt(label=scales::percent,limits=c(0,0.6))+
  scale_x_sqrt(breaks=brks,labels=brks)+
  theme_classic()+ggpubr::grids()+
  xlab('Age (years, square root scale)') +
  ylab('TBM risk (square root scale)')+
  theme(legend.position = c(0.7,0.8),legend.title = element_blank())+
  guides(color=guide_legend(ncol=2,byrow = TRUE),
         shape=guide_legend(ncol=2,byrow = TRUE))+
  scale_color_colorblind()

GP

## Note save for composition
Fig.A <- GP
save(Fig.A,file=gh('{xd}Fig.A.Rdata'))

ggsave(GP,file=gh('{xd}progression_withWallgren.png'),w=6,h=4)


## sqrt v sqrt to look like quadratic:
## sqrt_p = C + B*sqrt_a + A*a
## assume uniform distribution of ages within range
## driver =
## (1/W)\int_L^H da (C+B.sqrt(a)+A.a)
## (1/W)(C.W+(2/3).B.(H^1.5-L^1.5)+(1/2).A.(H^2-L^2))


PD <- data.table(acat=c('<1','1-4','5-9','10-14'),
                 L = c(0,1,5,10),
                 H = c(1,5,10,15))
PD[,W:=H-L]
PD[,prog:=(1/W)*(coef(A)[1]*W+(2/3)*coef(A)[2]*(H^1.5-L^1.5))]
PD[,prog:=transf.ilogit(prog)]

## uncertainty
N <- 5e3
PDL <- PD[rep(1:nrow(PD),each=N)]
PDL[,id:=rep(1:N,nrow(PD))]

## sample parms
pms <- MASS::mvrnorm(n=N,mu=coef(A),Sigma = vcov(A))
pms <- as.data.table(pms)
pms[,id:=1:N]

PDL <- merge(PDL,pms,by='id',all.x = TRUE)
PDL[,prog:=(1/W)*(intrcpt*W+(2/3)*sqrtT*(H^1.5-L^1.5))]
PDL[,prog:=transf.ilogit(prog)]

PD <- PDL[,.(prog.mn=mean(prog),prog.sd=sd(prog),prog.lo=lo(prog),prog.hi=hi(prog)),
    by=acat]
PD

save(PD,file=gh('{xd}PD.Rdata'))
fwrite(PD,file=gh('{xd}PD.csv'))


## === SA removing Wallgren
P <- P[Paper2!='Wallgren']


A <- rma.glmm(measure = "PLO",
              xi=TBM,ni=infected,data=P,
              mods=~1+sqrtT,
              slab=Paper)

## predictions
PP <- predict(object = A,newmods = sqrt(agz))
PP <- as.data.table(PP)
PP[,age.mid:=agz]
PP[,c('age.lo','age.hi','Paper','Paper2','Paper3'):=NA]
PP[,c('mid','lo','hi'):=.(exp(pred),exp(ci.lb),exp(ci.ub))]
PP[,c('lo2','hi2'):=.(exp(pi.lb),exp(pi.ub))]


## graph
GP <- ggplot(P,aes(age.mid,mid,
                   ymin=lo,ymax=hi,
                   xmin=age.lo,xmax=age.hi,
                   col=Paper3,shape=Paper3))+
  geom_ribbon(data=PP,aes(ymin=lo,ymax=hi),col=NA,alpha=0.3)+
  geom_line(data=PP,aes(y=lo2),lty=2,col='darkgrey')+
  geom_line(data=PP,aes(y=hi2),lty=2,col='darkgrey')+
  geom_line(data=PP,size=1,col='grey41')+
  geom_pointrange()+geom_errorbarh(height=0)+
  geom_point()+
  scale_y_sqrt(label=scales::percent,limits=c(0,0.6))+
  scale_x_sqrt(breaks=brks,labels=brks)+
  theme_classic()+ggpubr::grids()+
  xlab('Age (years, square root scale)') +
  ylab('TBM risk (square root scale)')+
  theme(legend.position = c(0.7,0.8),legend.title = element_blank())+
  guides(color=guide_legend(ncol=2,byrow = TRUE),
         shape=guide_legend(ncol=2,byrow = TRUE))+
  scale_color_colorblind()

## GP
ggsave(GP,file=gh('{xd}progression_noWallgren.png'),w=6,h=4)



PD <- data.table(acat=c('<1','1-4','5-9','10-14'),
                 L = c(0,1,5,10),
                 H = c(1,5,10,15))
PD[,W:=H-L]
PD[,prog:=(1/W)*(coef(A)[1]*W+(2/3)*coef(A)[2]*(H^1.5-L^1.5))]
PD[,prog:=transf.ilogit(prog)]

## uncertainty
N <- 5e3
PDL <- PD[rep(1:nrow(PD),each=N)]
PDL[,id:=rep(1:N,nrow(PD))]

## sample parms
pms <- MASS::mvrnorm(n=N,mu=coef(A),Sigma = vcov(A))
pms <- as.data.table(pms)
pms[,id:=1:N]

PDL <- merge(PDL,pms,by='id',all.x = TRUE)
PDL[,prog:=(1/W)*(intrcpt*W+(2/3)*sqrtT*(H^1.5-L^1.5))]
PDL[,prog:=transf.ilogit(prog)]

PD <- PDL[,.(prog.mn=mean(prog),prog.sd=sd(prog),prog.lo=lo(prog),prog.hi=hi(prog)),
    by=acat]
PD

save(PD,file=gh('{xd}PD_noWallgren.Rdata'))
fwrite(PD,file=gh('{xd}PD_noWallgren.csv'))
