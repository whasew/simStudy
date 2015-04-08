f2n<-function(f) as.numeric(levels(f))[f]

geomSeries <- function(base, max) {
    base^(0:floor(log(max, base)))
}

cumulationFactor<-function(tstar,k,n,tau) (1-exp(-n*k*tau))/(1-exp(-k*tau))

bateman<-function(Dose,time,fm,ka,CL,V,n=1,tau=24){
  K<-CL/V
  tstar<-time-(n-1)*tau
  cumulationFactor(tstar,k,n,tau)
  fm*ka*Dose/(ka-K)/V*(exp(-K*time)-exp(-ka*time))
}


xyLogscale <- function (lims = c(1e-04, 1000),labAdd=NULL) 
{
  labels <- c(1:9) * rep(10^c(-4:5), rep(9, 10))
  labels <- labels[labels >= lims[1] & labels <= lims[2]]
    if (!is.null(labAdd)) labels<-sort(c(labels,labAdd))
  at <- log10(labels)
  sel10 <- round(at - round(at), 6)==0
  if(!is.null(labAdd)) {
    att<-(at-log10(labAdd))
    selAdd <- round(att - round(att), 6)==0
  } else selAdd <-rep(FALSE,length(labels))
  
  sel<-labels%in%c(labels[sel10],labels[selAdd])
  if (sum(sel) > 1) 
    labels[!sel] <- ""
  list(log = TRUE, at = 10^at, labels = labels)
}
