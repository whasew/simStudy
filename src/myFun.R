######################################################################
## options
options(replace.assign=TRUE,width=70)
dev=c('pdf','png')
dev=c('tiff')
dev=c('tikz')
dev='CairoPNG'
dev='CairoPDF'

library(knitr)
library(reshape2)
library(languageR)
library(Hmisc)
library(nlme)
library(plyr)
library(ggplot2)
library(multcomp)
library(lattice)
if(F){
  lattice.options(default.args = list(page = function(n) { 
    panel.text(lab = sprintf("%s", date()), cex=.5,x = 0.95, y = 0.10, adj = 1) 
  })) 
}

super.sym<-super.sym.org<-trellis.par.get("superpose.symbol")
super.sym$col<-c("#0080ff" ,  "#ff00ff" ,  "darkgreen", "orange" ,   "#00ff00" ,  "brown",     "#ff0000"  )
trellis.par.set("superpose.symbol",super.sym)
##symbole()
##show.settings()
########################################################################
### functions
studDesign<-function(gd){
  table(getGroups(gd))
  table(getCovariate(gd), getGroups(gd))
}

panel.mydotplot <- function(x, y, subscripts, groups,ablv) { 
  dot.line <- trellis.par.get("dot.line")
  panel.abline(h = y, lwd = dot.line$lwd, lty = dot.line$lty, 
               col = dot.line$col)
  panel.abline(v=ablv,col=seq(along=ablv)+1)
  panel.superpose(x, y, subscripts, groups)
}

myPP<-function(gd,dL,cL,main,ablv,xlim=c(1E-1,2E0),inner=formula(~FA)){
  ppk<-list()
  ppk[[1]]<-plot(gd,displayLevel=dL,collapseLevel=cL,FUN=function(x) median(x),inner=inner,main=main,ablv=ablv)#,xlim=xlim)
  ppk[[1]]$panel<-panel.mydotplot
  ppk[[2]]<-plot(gd,displayLevel=dL,collapseLevel=cL,FUN=function(x) median(x),inner=inner,main=main,ablv=log10(ablv),xlim=xlim,scale=list(x=xyLogscale(xlim)))
  ppk[[2]]$panel<-panel.mydotplot
  ppk
}

ggpl<-function(gg,hline=0,exp=T,ylab="value",angle=NULL){
  #gg=ggplot(x,aes(x=CellType,y=estimate ,fill=Gene)
  if(exp) {
    gg$data[,chVec("estimate lwr upr")] <- 10^gg$data[,chVec("estimate lwr upr")] 
    #gg<-gg + scale_y_log10()
  }                    
  gg<- gg+ ylab(ylab)  +geom_bar(width=0.7,position=position_dodge(0.9),stat="identity")+geom_errorbar(aes(ymin=lwr,ymax=upr),width=.2,position=position_dodge(0.9)) + theme_bw()+scale_colour_grey(start = 0.2, end = 0.8,na.value = "grey50") +scale_fill_grey(start = 0.3, end = 0.8,  na.value = "grey50") 
  if(!is.null(angle)) gg<-gg+theme(axis.text.x=element_text(angle=45,hjust=1))
  if(hline==1) gg<-gg+geom_hline(yintercept=1) 
  gg
}
