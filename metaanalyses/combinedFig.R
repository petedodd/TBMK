## combining figures from other analyses
library(here)
library(ggpubr)
library(metafor)
library(data.table)
library(ggplot2)
library(ggthemes)
library(scales)
library(glue)
ilogit <- function(x)1/(1+exp(-x))
gh <- function(x) glue(here(x))
xd <- 'metaanalyses/plots/' #export/plot directory

## load figures
load(file=gh('{xd}Fig.A.Rdata'))
load(file=gh('{xd}Fig.B.Rdata'))
load(file=gh('{xd}Fig.C.Rdata'))

## modify
Fig.B <- Fig.B + theme_classic() + ggpubr::grids()
Fig.C <- Fig.C + theme_classic() + ggpubr::grids() +
  theme(strip.background = element_rect(fill = NA,colour = NA))



## combine & save
GP <- ggarrange(Fig.A,Fig.B,Fig.C,ncol=1,labels=c('A','B','C'))
ggsave(GP,file=gh('{xd}Fig.ALL.pdf'),w=7,h=15)
ggsave(GP,file=gh('{xd}Fig.ALL.png'),w=7,h=15)
