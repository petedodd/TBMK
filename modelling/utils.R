## some utilities for formatting & outputs



## converter from matrix to data.table
mkdt <- function(D){
  nmz <- rownames(D)
  D <- as.data.table(D)
  D[,variable:=nmz]
  return(D)
}

## functions for reformatting outputs
getcno <- function(x){
  a <- str_extract(x,"\\[(.*?),")
  a <- gsub('\\[','',a)
  a <- gsub(',','',a)
  as.integer(a)
}

getAno <- function(x){
  a <- str_extract(x,",(.*?)\\]")
  a <- gsub('\\]','',a)
  a <- gsub(',','',a)
  as.integer(a)
}

acts <- c('<1','1-4','5-9','10-14')
isoz <- scan(gh('{dd}isoz.txt'),what='char')

getcn <- function(x) isoz[getcno(x)]
getac <- function(x) acts[getAno(x)]



## rounding tables
## always round 0.5 up
round2 <- function(x, digits = 0) sign(x) * trunc(abs(x) * 10^digits + 0.5) / 10^digits
ft <- Vectorize(function(x){
  smallpos <- x > 0 & x < 0.01
  one2ten <- x >= 1 & x < 10
  zero2one <- x >= 0.1 & x < 1
  dg <- ifelse(abs(x) > 0.01 & abs(x) < 100, 2, 3)
  x2 <- signif(x, dg)
  trailing.0 <- x2 == round2(x) & one2ten == TRUE
  trailing0 <- x2 * 10 == round2(x * 10) & zero2one == TRUE & x2 < 1
  format(
    x2,
    digits = dg,
    nsmall = 0L,
    big.mark = " ",
    justify = 'right',
    drop0trailing = TRUE,
    scientific = FALSE
  )
})
brkt <- function(x,y) paste0(ft(pmax(0,x)),' (',
                             ft(pmax(0,x-1.96*y)),' to ',
                             ft(pmax(0,x+1.96*y)),')')
