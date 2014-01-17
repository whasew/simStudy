Design and Simulation of an gene expression experiment
======================================================

# Questions
- is gene expression different in differentiated cells?
- monocytes differentiates to macrophages
- drug concentrations induce gene expression?

# Design

## 5 cases
 case | differentiated  | drug concentration  
 ------------- | ------------- | -------------  
 1 | - | 0  
 2  | + | 0  
 3 | + | 2.5  
 4 | + | 5  
 5 | + | 7  

5 genes and duplicate measurement of each sample ID


```r
library(MASS)
library(nlme)
source("../src/myFun.R")
```

```
## Loading required package: survival
```

```
## Loading required package: splines
```

```
## Hmisc library by Frank E Harrell Jr
## 
## Type library(help='Hmisc'), ?Overview, or ?Hmisc.Overview') to see
## overall documentation.
## 
## NOTE:Hmisc no longer redefines [.factor to drop unused levels when
## subsetting.  To get the old behavior of Hmisc type
## dropUnusedLevels().
```

```
## Attaching package: 'Hmisc'
```

```
## Das folgende Objekt ist maskiert from 'package:survival':
## 
## untangle.specials
```

```
## Die folgenden Objekte sind maskiert from 'package:base':
## 
## format.pval, round.POSIXt, trunc.POSIXt, units
```

```
## Attaching package: 'plyr'
```

```
## Die folgenden Objekte sind maskiert from 'package:Hmisc':
## 
## is.discrete, summarize
```

```
## Loading required package: mvtnorm
```

```r
source("../src/helpers.R")
```



```r
mv <- c(HPRT = 1000, ABCD1 = 900, ABCD2 = 10, ABCD3 = 700)
conc <- c(0, 0, 2.5, 5, 7)
diff <- c(0, rep(1, length(conc) - 1))
hill <- 20
fold <- c(HPRT = 1.3, ABCD1 = 1.3, ABCD2 = 5, ABCD3 = 5)
ec50 <- c(HPRT = 10, ABCD1 = 7, ABCD2 = 5, ABCD3 = 7)
eps <- 0.2
Sigma <- diag(rep(0.1^2, length(mv)))
selGene <- c("HPRT", paste("ABCD", 1:3, sep = ""))
myF <- function(id, n) paste(id, n, sep = "")
yy <- function(x, ymin, fa, ec50) ymin * (1 + (fa - 1) * plogis(hill * 
    log(x/ec50)))
nRep <- 3
nId <- length(conc) * nRep
dd <- expand.grid(Case = seq(along = conc), Rep = 1:nRep, Obs = 1:2, Gene = selGene)
dd$FCase <- factor(myF("C", dd$Case), levels = myF("C", seq(along = conc)))
dd$ID <- paste("S", dd$Case, dd$Rep, sep = "")

dd <- transform(dd, ID = factor(ID, level = paste("S", rep(seq(along = conc), 
    each = nRep), rep(1:nRep, length(conc)), sep = "")), Conc = conc[Case], 
    mRNA = yy(conc[Case], mv[Gene], fold[Gene], ec50[Gene]))

rand <- with(dd, exp(mvrnorm(length(levels(ID)), rep(0, length(mv)), Sigma = Sigma)[ID]) * 
    mRNA)
dd$RNA <- rand
dd$obs <- exp(rnorm(nrow(dd), sd = eps)) * rand
head(dd)
```

```
##   Case Rep Obs Gene FCase  ID Conc mRNA    RNA    obs
## 1    1   1   1 HPRT    C1 S11  0.0 1000  951.7  954.6
## 2    2   1   1 HPRT    C2 S21  0.0 1000 1005.9  791.7
## 3    3   1   1 HPRT    C3 S31  2.5 1000  884.6 1181.8
## 4    4   1   1 HPRT    C4 S41  5.0 1000  927.6 1211.8
## 5    5   1   1 HPRT    C5 S51  7.0 1000  910.0  754.9
## 6    1   2   1 HPRT    C1 S12  0.0 1000 1062.2  980.8
```



```r
gd <- groupedData(obs ~ Conc | ID, dd)
plot(gd, scale = list(y = xyLogscale()), aspect = "fill", outer = ~Gene)
```

![plot of chunk gd](figure/gd1.png) 

```r
plot(gd, outer = ~Gene, aspect = "fill", scale = list(y = c(xyLogscale(c(1, 
    50000)), list(relation = "free"))))
```

![plot of chunk gd](figure/gd2.png) 



```r
gd <- groupedData(obs ~ FCase | ID/Obs, dd)
plot(gd, displayLevel = 1, collapseLevel = 1, , scale = list(x = xyLogscale(c(1, 
    10000))), aspect = "fill", inner = ~Gene)
```

![plot of chunk gdF](figure/gdF.png) 



# lme


```r
fms1 <- lme(log10(obs) ~ Gene * FCase, gd, random = ~1 | ID)
anova(fms1)
```

```
##             numDF denDF F-value p-value
## (Intercept)     1    90   37528  <.0001
## Gene            3    90    2969  <.0001
## FCase           4    10      20   1e-04
## Gene:FCase     12    90      21  <.0001
```

```r
summary(fms1)$tT
```

```
##                       Value Std.Error DF    t-value   p-value
## (Intercept)        3.004204   0.04256 90  70.588413 1.357e-80
## GeneABCD1         -0.049510   0.05062 90  -0.978062 3.307e-01
## GeneABCD2         -1.944737   0.05062 90 -38.417659 1.296e-57
## GeneABCD3         -0.193569   0.05062 90  -3.823899 2.418e-04
## FCaseC2           -0.012134   0.06019 10  -0.201608 8.443e-01
## FCaseC3           -0.029122   0.06019 10  -0.483843 6.389e-01
## FCaseC4            0.021335   0.06019 10   0.354466 7.304e-01
## FCaseC5           -0.023888   0.06019 10  -0.396896 6.998e-01
## GeneABCD1:FCaseC2  0.026404   0.07159 90   0.368824 7.131e-01
## GeneABCD2:FCaseC2 -0.047504   0.07159 90  -0.663574 5.087e-01
## GeneABCD3:FCaseC2  0.066664   0.07159 90   0.931213 3.542e-01
## GeneABCD1:FCaseC3  0.005019   0.07159 90   0.070109 9.443e-01
## GeneABCD2:FCaseC3 -0.131758   0.07159 90  -1.840489 6.899e-02
## GeneABCD3:FCaseC3  0.029189   0.07159 90   0.407734 6.844e-01
## GeneABCD1:FCaseC4 -0.007448   0.07159 90  -0.104040 9.174e-01
## GeneABCD2:FCaseC4  0.429220   0.07159 90   5.995625 4.141e-08
## GeneABCD3:FCaseC4 -0.000407   0.07159 90  -0.005686 9.955e-01
## GeneABCD1:FCaseC5  0.110350   0.07159 90   1.541448 1.267e-01
## GeneABCD2:FCaseC5  0.590636   0.07159 90   8.250398 1.240e-12
## GeneABCD3:FCaseC5  0.487153   0.07159 90   6.804873 1.089e-09
```

```r
gdInt <- intervals(fms1, levels = 0.68)
```


# barplot



```r
int0 <- gdInt$fix
int <- int0[grep("^Gene", rownames(int0)), ]
int
```

```
##                      lower      est.    upper
## GeneABCD1         -0.15008 -0.049510  0.05106
## GeneABCD2         -2.04530 -1.944737 -1.84417
## GeneABCD3         -0.29414 -0.193569 -0.09300
## GeneABCD1:FCaseC2 -0.11582  0.026404  0.16863
## GeneABCD2:FCaseC2 -0.18973 -0.047504  0.09472
## GeneABCD3:FCaseC2 -0.07556  0.066664  0.20889
## GeneABCD1:FCaseC3 -0.13720  0.005019  0.14724
## GeneABCD2:FCaseC3 -0.27398 -0.131758  0.01047
## GeneABCD3:FCaseC3 -0.11303  0.029189  0.17141
## GeneABCD1:FCaseC4 -0.14967 -0.007448  0.13478
## GeneABCD2:FCaseC4  0.28700  0.429220  0.57144
## GeneABCD3:FCaseC4 -0.14263 -0.000407  0.14182
## GeneABCD1:FCaseC5 -0.03187  0.110350  0.25257
## GeneABCD2:FCaseC5  0.44841  0.590636  0.73286
## GeneABCD3:FCaseC5  0.34493  0.487153  0.62938
```

```r
rownames(int) <- gsub("FCase", "", gsub("Gene", "", rownames(int)))
sel <- grep(":", rownames(int))
int <- int[sel, ]

tmp <- unPaste(rownames(int), sep = ":")
tmp
```

```
## [[1]]
##  [1] "ABCD1" "ABCD2" "ABCD3" "ABCD1" "ABCD2" "ABCD3" "ABCD1" "ABCD2"
##  [9] "ABCD3" "ABCD1" "ABCD2" "ABCD3"
## 
## [[2]]
##  [1] "C2" "C2" "C2" "C3" "C3" "C3" "C4" "C4" "C4" "C5" "C5" "C5"
```

```r
inta <- cbind(data.frame(int), Gene = tmp[[1]], FCase = tmp[[2]])
intal <- daply(inta, .(Gene), function(x) {
    y <- as.matrix(x[, 1:3])
    dimnames(y)[[1]] <- x$FCase
    y
})
# daply(inta,.(Gene),function(x) {y<-as.matrix(x[,1:3]);y})

intal
```

```
## , ,  = lower
## 
##        
## Gene          C2      C3      C4       C5
##   ABCD1 -0.11582 -0.1372 -0.1497 -0.03187
##   ABCD2 -0.18973 -0.2740  0.2870  0.44841
##   ABCD3 -0.07556 -0.1130 -0.1426  0.34493
## 
## , ,  = est.
## 
##        
## Gene          C2        C3        C4     C5
##   ABCD1  0.02640  0.005019 -0.007448 0.1104
##   ABCD2 -0.04750 -0.131758  0.429220 0.5906
##   ABCD3  0.06666  0.029189 -0.000407 0.4872
## 
## , ,  = upper
## 
##        
## Gene         C2      C3     C4     C5
##   ABCD1 0.16863 0.14724 0.1348 0.2526
##   ABCD2 0.09472 0.01047 0.5714 0.7329
##   ABCD3 0.20889 0.17141 0.1418 0.6294
```

```r
# barPlotsInt(intal,T)
x <- intal
ym <- x[, , 2]
yup <- x[, , 3]
ylow <- x[, , 1]
LOG <- F
if (!LOG) {
    ym <- 10^ym
    yup <- 10^yup
    ylow <- 10^ylow
}
# LAS=2;lay=NULL;ylim=NULL;col=0 barx <-
# barplot(ym,beside=T,col=col,las=LAS,ylim=ylim,ylab=ylab,legend.text=T,args.legend
# = list(x = leg,ncol=2,bty = 'n'),layout=lay)

barx <- barplot(ym, beside = T, col = 2:4, las = 1, legend.text = T, args.legend = list(x = "topleft", 
    ncol = 2, bty = "n"))
abline(h = 1)
error.bar(barx, ym, yup - ym, ym - ylow)
```

![plot of chunk barplot](figure/barplot.png) 




