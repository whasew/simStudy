helper.function <- function()
{
  return(1)
}


myMelt <- function(x,mvar,varn="variable",valn="value"){
melt(x,measure.var=mvar,variable.name=varn,value.name=valn) 
}


stripLocal <-function(...) strip.default(...,strip.names=c(T,T),style=1)

redFactor <- function(x){
 for(v in names(x)[sapply(x,data.class)=="character"]) x[,v] <-  factor(gsub(" ","",x[,v]))
  for(v in names(x)[sapply(x,data.class)=="factor"]) x[,v] <-  factor(gsub(" ","",x[,v]))
for(v in names(x)[sapply(x,data.class)=="ordered"]) x[,v] <-  ordered(gsub(" ","",x[,v]))
  x
}


lettern <- function(n,selCol=1:26){
  sel <- selCol+rep(c(0:n)*26,each=26)
  names(sel) <- c(LETTERS,paste(rep(LETTERS[1:(n)],each=26),LETTERS,sep=""))
  sel
}
plot2pp <- function(x,rel=NULL){
    tt <- paste("plot(",x,")",sep="")
    if(!is.null(rel)) tt <- paste("plot(",x,rel,")",sep="")
    ttlog <- paste("plot(",x,",scale=list(y=xyLogscale()))",sep="")
  list(
    eval(parse(text=tt)),
    eval(parse(text=ttlog))
       )
  }

preXlsx <- function(myFile,gen,fixClasses){
  colLetter <- c(LETTERS[1:7],gen)
  letterNo <- seq(along=LETTERS)
  names(letterNo) <- LETTERS
  colClass <- rep("numeric", length(gen));names(colClass) <- names(gen)
  list(file=myFile,colIndex=letterNo[names(letterNo)%in%colLetter],
  colClasses=c(fixClasses,colClass))
}


error.bar <- function(x,y,upper,lower=upper,length=0.03,...){
if(length(x)!=length(y)|length(y)!=length(lower)|length(lower)!=length(upper))
  stop("vectors must be same length")
arrows(x,y+upper,x,y-lower,angle=90,code=3,length=length,...)
}

bilde.pfad<-function(...){
  file.path(mypath,...)
}

chVec <- function(x,sep=" ") unlist(strsplit(x,sep))


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

pdfPrintpp<-function(pp,prefix="graphs/pp",color = TRUE,dev="pdf") {
  trellis.device(device = pdf,
                 file=paste(prefix,"%03d.pdf",sep=""),
                 onefile=FALSE,
                 color = color,
                 theme = lattice.getOption("default.theme"),
                 new = TRUE,
                 retain = FALSE)
  print(pp)
  graphics.off()
  invisible()
}

"show.info" <-
function(message,variable,DEBUG=TRUE){
  if(DEBUG){
    if(!missing(message )) cat("\n",message)
    if(!missing(variable)){
      cat(" value of ",substitute(variable),": \n",sep="")
      print(variable)
    }
  }
}


cor2cov <-
    function(cov = NULL, sd = NULL)
{
## cor2cov and cov2cor
## Re: [S] How to calculate cavriance from correlation matrix Mo 25.02.2002 20:22
## typo (%% modulus to  %*% matrix multiplication)corrected by
## Prof Brian D Ripley [ripley@stats.ox.ac.uk], M�inen Jussi [jussi.makinen@tapiola.fi]
	warning("UNDER CONSTRUCTION, MUST !!!")
	if(is.null(cov)) {
		## Re: [S] simulate correlation matrices, Di 26.02.2002 03:56
		set.seed(260202)
		x <- matrix(rnorm(5 * 5), ncol = 5)
		cov <- t(x) %*% x
		##    corMatrix <- cor(matrix( rnorm(50), 10,5))
		print(cov)
	}
	dd <- diag(cov)
	corMat <- length(unique(dd))
	if(dd[corMat] == 1) {
		sd2matrix <- outer(sd, sd, "*")
		cov <- sd2matrix * cov
		tmp <- list(covMatrix = cov)
	}
	else {
		sd <- sqrt(diag(cov))
		sd2matrix <- outer(sd, sd, "*")
		cov <- cov/sd2matrix
		tmp <- list(corMatrix = cov, SD = sd)
	}
        
	tmp
}

###################################################
### chunk number 3: symbole
###################################################
symbole<-function (PRINT = FALSE, ask = FALSE) 
{                                     
    oldpar <- par(c("mar", "oma", "mfrow", "ask"))
    par(ask = ask)                                
    par(mar = c(0, 0, 0, 0), mfrow = c(6, 1))     
    farbe <- colours()                            
    xx <- 0:25                                    
    plot(xx, xx/(xx + 1.9), type = "n", axes = FALSE, xlab = "", 
        ylim = c(0.2, 0.55), ylab = "", cex = 2)                 
    for (i in xx) {                                              
        points(i, 0.5, pch = i, cex = 4, bg = farbe[i])          
        if (i < 9)                                               
            points(i, 0.33, pch = 16, cex = 6, col = i)          
        text(i, 0.42, i, cex = 1.6)                              
        if (i < 7)                                               
            lines(c(i, i), c(0.2, 0.325), lty = i, lwd = 2)      
    }                                                            
    text(16, 0.3, "following colour names throught:", cex = 2)
    text(16, 0.22, "> colours()[<number>]           ", cex = 2)
    ppl <- 25
    if (ask)
        nplots <- ceiling(length(farbe)/ppl)
    else nplots <- par("mfrow")[1] - 2
    for (y in 0:nplots) {
        xx <- 1:ppl + y * ppl
        plot(xx, xx/(xx + 1.9), type = "n", axes = FALSE, xlab = "",
            ylim = c(0.2, 0.55), ylab = "", cex = 2)
        for (i in xx) {
            if (i <= length(farbe)) {
                points(i, 0.33, pch = 16, cex = 6, col = farbe[i])
                text(i, 0.42, i, cex = 1)
            }
        }
    }
    par(oldpar)
   invisible()
}


error.bar <- function(x,y,upper,lower=upper,length=0.03,...){
if(length(x)!=length(y)|length(y)!=length(lower)|length(lower)!=length(upper))
  stop("vectors must be same length")
arrows(x,y+upper,x,y-lower,angle=90,code=3,length=length,...)
}

