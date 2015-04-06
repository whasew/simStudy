
library(foreach)
library(ggplot2)
library(plyr)
source("src/stan_run.R")
source("src/matrixODE.R")
source("lib/helpers.R")
theme_set(theme_bw())

## use two cores
options(mc.cores=2)

logit <- function(p) log(p/(1-p))
inv_logit <- function(l) 1/(1+exp(-l))


## population medians, means for the log-normals (or the geometric mean)
MRT <- c(2.5, 10, 15)
V <- c(10, 20, 15)
fm1 <- 0.6 ## note: this is FIXED right now
fm2 <- 0.3 ## note: this is FIXED right now

truth <- c(MRT, V)
Ttruth <- log(truth)

## number of patients
J <- 30
time <- exp(seq(log(0.1),log(48),length=16))
##J <- 10
time <- exp(seq(log(0.1),log(48),length=11))
time <- exp(c(seq(log(0.1),log(24),length=11), seq(log(24.5), log(48), length=6)))

##J <- 2
##time <- seq(0,10,by=1)

## variance components
## random effects on CL and V
omega <- c(0.3, 0.2)
## residual error
sigma <- c(0.25, 0.1, 0.15)

## IV dose injected
IVdose <- 50
rate <- IVdose/(1.5*48) ## infusion which is longer than the overall time-span
rate <- IVdose/(0.5*48) ## infusion which is shorter than the overall time-span

## simulated weight distribution
weight_gmean <- 70
weight_sd95 <- 0.2


## utility functions to generate system matrices
systemODE_1cmt_m2 <- function(lMRT10, lMRT20, lMRT30, lfm1, lfm2) {
    k10 <- exp(-lMRT10)
    k12 <- exp(lfm1) * k10
    k13 <- exp(lfm2) * k10
    k20 <- exp(-lMRT20)
    k30 <- exp(-lMRT30)
    A <- rbind(
        c(-k10,    0,    0),
        c( k12, -k20,    0),
        c( k13,    0, -k30))
    A
}

## k   = CL/V
## MRT = V/CL
## CL  = V/MRT



## generate parameters per patient
pars <- data.frame(id=1:J
                  ,xi_lMRT10=rnorm(J, 0, 1)
                  ,xi_lV1=rnorm(J, 0, 1)
                  ,lMRT20=log(MRT[2])
                  ,lMRT30=log(MRT[3])
                  ,lV2=log(V[2])
                  ,lV3=log(V[3])
                  ,lfm1=log(fm1)
                  ,lfm2=log(fm2)
                  ,weight=rlnorm(J, log(weight_gmean), weight_sd95/1.96)
                   )
## add in the covariate effect of weight into the clearance and volume

pars <- mutate(pars
              ,lV1=log(V[1]) + xi_lV1 * omega[2] + 1.00 * log(weight/70)
              ,lCL=log(V[1]) - log(MRT[1]) + 0.75 * log(weight/70)
              ,lMRT10=lV1-lCL + xi_lMRT10 * omega[1]
               )

## lMRT10 = log(MRT10[1]) - 0.75 * log(weight/70) + eta_lMRT10

pars <- transform(pars
                 ,lk12=log(fm1) - lMRT10
                 ,lk13=log(fm2) - lMRT10
                 )

##qplot(log(weight), lCL, data=pars) + stat_smooth(se=FALSE, method=lm)
##qplot(log(weight), lV, data=pars) + stat_smooth(se=FALSE, method=lm)

summary(lm(lV1 ~ log(weight/70), data=pars))
summary(lm(lCL ~ log(weight/70), data=pars))

## generate Nonmem type data set
sim <- foreach(j=1:J, .combine=rbind) %do% {
    sol <- matrixODE(do.call(systemODE_1cmt_m2, pars[j,c("lMRT10", "lMRT20", "lMRT30", "lfm1", "lfm2")]))
    sol1 <- cmt_select(sol, 1)
    sol2 <- cmt_select(sol, 2)
    sol3 <- cmt_select(sol, 3)
    IVdose_j <- IVdose ##rlnorm(1, log(IVdose), log(1.1))
    ##x0 <- c(IVdose_j,0,0)
    ##b <- c(0, 0, 0)
    x0 <- c(0,0,0)
    b <- c(rate, 0, 0)
    b2 <- c(0, 0, 0)
    TinfEnd <- IVdose/rate
    time1 <- time[time <= TinfEnd + 1e-6]
    time2 <- time[time >  TinfEnd + 1e-6] - TinfEnd
    time2
    dd1 <- as.data.frame(t(sol(time1, x0, b)))
    dd1$time <- time1
    sol(TinfEnd, x0, b)
    dd2 <- as.data.frame(t(sol(time2, as.vector(sol(TinfEnd, x0, b)), b2)))
    dd2$time <- time2 + TinfEnd
    dd <- transform(rbind(dd1, dd2), V1=V1/exp(pars$lV1[j]), V2=V2/exp(pars$lV2[j]), V3=V3/exp(pars$lV3[j]))
    M <- melt(dd, id.vars="time")
    obs <- data.frame(id=j, sdv=M$value, time=M$time, amt=0, rate=0, evid=0, mdv=0, cmt=as.numeric(M$variable))
    obs$dv <- obs$sdv
    obs$dv[obs$cmt==1] <- obs$dv[obs$cmt==1] * rlnorm(nrow(dd), 0, sigma[1])
    obs$dv[obs$cmt==2] <- obs$dv[obs$cmt==2] * rlnorm(nrow(dd), 0, sigma[2])
    obs$dv[obs$cmt==3] <- obs$dv[obs$cmt==3] * rlnorm(nrow(dd), 0, sigma[3])
    amt <- data.frame(id=j, sdv=0, time=0, amt=IVdose_j, rate=b[1], evid=1, mdv=1, cmt=1, dv=0 )
    res <- arrange(rbind(amt, obs), time, cmt)
    cbind(res, weight=pars$weight[j])
}

## discard dv when system is still equal to 0 as we calculate in log-space
sim <- subset(sim, !(time==0 & evid==0))

head(sim)


## look at simulation
qplot(time, dv, data=subset(sim, mdv==0 & cmt==1), colour=factor(id), geom="line") + theme(legend.position="none")
qplot(time, dv, data=subset(sim, mdv==0 & cmt==2), colour=factor(id), geom="line")+ theme(legend.position="none")
qplot(time, dv, data=subset(sim, mdv==0 & cmt==3), colour=factor(id), geom="line")+ theme(legend.position="none")
qplot(time, dv, data=subset(sim, mdv==0 & cmt==1), colour=factor(id), geom=c("line", "point"), log="y") + theme(legend.position="none")
qplot(time, dv, data=subset(sim, mdv==0 & cmt==2), colour=factor(id), geom=c("line", "point"), log="y")+ theme(legend.position="none")



## define prior prior for Q parametrization; 4th parameter is
## logit(fm) which must be given an informative prior given
## identifiability issues; see at the end of the script for guidance
priorQ <- list(prior_ltheta_mean=Ttruth,
               prior_ltheta_sd95 =rep(log(2), length(Ttruth)),
               prior_sigma_sd95 = c(5, 5, 5),
               prior_sigma_gmean = c(0.2, 0.2, 0.2),
               random_effect=c(1, 0, 0, 1, 0, 0),
               prior_omega_sd95 = array(5, dim=2),
               prior_omega_gmean =array(0.25, dim=2)
              )

## we use the rate constants by default
prior <- priorQ
E <- sum(prior$random_effect)

## this init would be cheating, i.e. we would start of the true values
init <- function()
    list(ltheta=prior$prior_ltheta_mean,
             sigma=c(0.1, 0.1, 0.1),
             omega=array(0.5, dim=E),
             xi_eta=array(0, dim=c(J,E)))
             ##xi_eta=array(rnorm(J*E, 0, 0.5), dim=c(J,E)))


## first use the truth to let Stan calculate the simulated values

sim_params <- list(
    ltheta=Ttruth,
    xi_eta=array(cbind(pars$xi_lMRT10, pars$xi_lV1), dim=c(J,2)),
    sigma=sigma,
    omega=array(c(omega[1], omega[2]+1e-6), 2)
)


N <- nrow(sim)
pkMetaFix <- stan_run("SWlinearPkMeta2FmFixInf.stan",
                     data=c(as.list(sim), list(N=N, BC=1, base=as.matrix(log(sim$weight/70))), prior),
                     init=list(sim_params),
                     ##init=0,
                     algorithm="Fixed_param",chains=1,iter=1, warmup=0
                     )

stanRes <- extract(pkMetaFix, "ipre")$ipre[1,,]
head(sim)

cmp <- transform(sim, standv=exp(stanRes[cbind(1:N,sim$cmt)]))
cache("cmp")
qplot(standv, sdv, data=cmp)
qplot(time, standv - sdv, data=cmp)

qplot(standv, sdv, data=cmp) + facet_wrap(~cmt) + geom_abline()

qplot(time, log(standv)-log(sdv), data=cmp) + facet_wrap(~cmt) + geom_hline(yintercept=0)

qplot(time, standv - sdv, data=cmp) + facet_wrap(~cmt) + geom_hline(yintercept=0)

qplot(time, sdv, data=cmp, colour=factor(cmt), geom="line") +
    geom_line(aes(y=standv), linetype=2) + facet_wrap(~id)

qplot(time, sdv, data=cmp, colour=factor(cmt), geom="line") +
    geom_line(aes(y=standv), linetype=2) + facet_wrap(~id)


N <- nrow(sim)
pkMetaFit <- stan_run("SWlinearPkMeta2FmFixInf.stan",
                      data=c(as.list(sim), list(N=N, BC=1, base=as.matrix(log(sim$weight/70), ncol=1)), prior),
                      init=init,
                      ##init=0,
                      chains=2, iter=1000/2, warmup=500/2, parallel=TRUE,
                      control=list(stepsize=1e-2)
                      ##control=list(stepsize=1e-3, adapt_delta=0.95)
                      )
#cache("pkMetaFit")
print(pkMetaFit,
      pars=c("ltheta", "theta", "CL", "omega", "sigma"),
      digits_summary=3,
      prob=c(0.025,0.5,0.975))

post <- extract(pkMetaFit, pars=c("ltheta", "omega", "sigma"))

## check if estimates are in line with our simulated values (these are
## posterior p-values and should be not to close to 0 or 1)
colMeans(sweep(post$ltheta, 2, Ttruth) > 0)
##colMeans(sweep(post$omega, 2, omega) > 0)
mean(post$omega-omega > 0)
colMeans((post$sigma-sigma) > 0)

traceplot(pkMetaFit, pars="ltheta", window=100)
traceplot(pkMetaFit, pars=c("sigma", "omega"), window=100)


############################# BELOW HERE IS EXTRA CHECKING ##################################


## try out Ben's code to make evaluate_model function available in R;
## does not work for me as of now
if(FALSE) {
    library(devtools)
    devtools::source_url("https://raw.githubusercontent.com/stan-dev/rstan/develop/rstan/rstan/R/testify.R")
    testify(pkTwoFit)

    pp <- extract(pkTwoFit, c("xi_eta", "omega", "ltheta"))

    eta <- sweep(pp$xi_eta, c(1,3), pp$omega, "*")

    ltheta_est <- colMeans(pp$ltheta)
    eta_est <- t(apply(eta, c(2), colMeans))

    params <- matrix(ltheta_est, ncol=4, nrow=J, byrow=TRUE)
    params[,1] <- params[,1] + eta_est[,1]
    params[,2] <- params[,2] + eta_est[,2]

    M <- tapply(sim$id,sim$id,length)

    ldvFit <- evaluate_model(M, params, sim$time, sim$amt, sim$rate, sim$evid)

    res <- transform(sim, ipred=exp(ldvFit))

    ggplot(subset(res, id < 10), aes(time)) +
        geom_line(aes(y=sdv), linetype=2) +
            geom_point(aes(y=dv)) +
                geom_line(aes(y=ipred)) +
                    facet_wrap(~id)

}


## compare histograms of the posterior with the truth
library(reshape2)

fitVars <- paste("l", c("k", "k21", "kMean", "V"), sep="")
M <- transform(melt(post$ltheta), par=factor(Var2, levels=1:4, labels=fitVars))

ref <- data.frame(par=factor(1:4, labels=fitVars), value=log(alt))

ggplot(M, aes(value)) + geom_histogram(alpha=0.8) + facet_wrap(~par, scale="free") +
    geom_vline(data=ref, aes(xintercept=value))




N <- nrow(sim)
pkTwoFix <- stan_run("SWtwoCompMultiAnalyticFitL2.stan",
                     data=c(as.list(sim), list(N=N, BC=0, base=matrix(0, nrow=N, ncol=0)), prior),
                     init=init,
                     ##init=0,
                     algorithm="Fixed_param",chains=1,iter=1, warmup=0
                     )

simS <- subset(sim, !(time==0 & evid==0))

N <- nrow(simS)
pkMetaFix <- stan_run("SWlinearPkMultiAnalyticFitL2.stan",
                     data=c(as.list(simS), list(N=N, BC=1, base=as.matrix(log(simS$weight/70), ncol=1)), prior),
                     init=init,
                     ##init=0,
                     algorithm="Fixed_param",chains=1,iter=1, warmup=0
                     )


## get prior for logit(fm) which we expect at 60% and assume it is in
## the range of 50% to 70% -> adjust accordingly!

## formula for prior specification of normals by setting two quantiles
normPrior <- function(p1, q1, p2, q2) {
  z1 <- qnorm(p1)
  z2 <- qnorm(p2)
  mu <- (z2 * q1 - z1 * q2)/(z2-z1)
  sigma <- (q2-q1)/(z2-z1)
  c(mu=mu, sigma=sigma)
}

priorFm <- normPrior(0.025, logit(0.55), 0.975, logit(0.65))

priorFm

## prior in logit space
curve(dnorm(x, priorFm["mu"], priorFm["sigma"]), -2, 2)

## in 0-1 space
samp <- inv_logit(rnorm(5000, priorFm["mu"], priorFm["sigma"]))
hist(samp)


k10 <- 0.5
k20 <- 0.1
k12 <- 0.2

A <- matrix(NA, 2, 2)

A[1,1] <- -k10;
A[1,2] <- 0;
A[2,1] <- k12;
A[2,2] <- -k20;


Ai <- matrix(NA, 2, 2)

Ai[1,1] <- - 1. / k10 ;
Ai[1,2] <- 0.;
Ai[2,1] <- - k12 / ( k10 * k20 );
Ai[2,2] <- - 1. / k20 ;

print(A)
print(Ai)
print(A %*% Ai)
print(Ai %*% A)

ev[1] <- -k10;
ev[2] <- -k20;

detU <- (k20 - k10) / k12;
U <- matrix(NA, 2, 2)
U[1,1] <- detU;
U[2,1] <- 1.;
U[1,2] <- 0.;
U[2,2] <- 1.;

Ui <- matrix(NA, 2, 2)
Ui[1,1] <-  1. / detU;
Ui[1,2] <-  0.;
Ui[2,1] <- -1. / detU;
Ui[2,2] <- 1.;

print(U)
print(Ui)
print(U %*% Ui)
print(Ui %*% U)


sim_params <- list(
    ltheta=Ttruth,
    xi_eta=array((pars$lCL - log(CL)) / omega[1], dim=c(J,1)),
    sigma=sigma,
    omega=array(omega[1], 1)
)


N <- nrow(sim)
pkMetaFix <- stan_run("SWlinearPkMetaFmFix.stan",
                     data=c(as.list(sim), list(N=N, BC=1, base=as.matrix(log(sim$weight/70))), prior),
                     init=list(sim_params),
                     ##init=0,
                     algorithm="Fixed_param",chains=1,iter=1, warmup=0
                     )

stanRes <- extract(pkMetaFix, "ipre")$ipre[1,,]

head(sim)

cmp <- transform(sim, standv=exp(stanRes[cbind(1:N,sim$cmt)]))

qplot(standv, sdv, data=cmp)


## check the analytic solution

solA1 <- function(lCL, lV, lk20, Lfm) {
    k10 <- exp(lCL - lv)
    k12 <- inv_logit(lFm) * k10
    k20 <- exp(lk20)

    nldetU <- log( (k10 - k20) / k12 )

    function(time, lx0) {
        c(lx0[1] - k10 * time, log( exp(x0[1] - nldetU + log(exp(-k20 * time) - exp(-k10 * time) )) + exp(lx0[2] - k20 * time  ) )   )
    }
}

solA2 <- function(lk10, lk12, lk20) {
    k10 <- exp(lk10)
    k12 <- exp(lk12)
    k20 <- exp(lk20)

    nldetU <- log( (k10 - k20) / k12 )

    function(time, lx0) {
        cbind(cmt1=lx0[1] - k10 * time, cmt2=log( exp(lx0[1] - nldetU + log(exp(-k20 * time) - exp(-k10 * time) )) + exp(lx0[2] - k20 * time  ) )   )
    }
}


lk10 = -1.04678
lk12 = -1.55761
lk20 = -2.6691
nldetU = 0.290875


sf <- solA2(lk10, lk12, lk20)

dd <- subset(sim, id==10)


lx0 <- log(c(IVdose, 1e-15))
as <- sf(dd$time, lx0)

dd$ref <- exp(as[cbind(1:25,dd$cmt)])

dd


library(deSolve)
model <- function(time, y, p) {
  with(as.list(c(y,p)), {
    Dcent <- -k10 * cent
    Dmet <- k12 * cent - k20 * met
    list(c(Dcent, Dmet))
  })
}


y0 <- c(cent=IVdose, met=0)

out <- ode(y = y0,
           times = time,
           func = model,
           parms = c(k10=CL/V, k12=fm * CL/V, k20=k20)
           )

sa <- solA2(log(CL/V), log(fm * CL/V), log(k20))

outA <- exp(sa(time, c(log(IVdose), -30)))

dd <- cbind(as.data.frame(out), outA)
qplot(time, cent, data=dd, log="y", geom="line") + geom_line(aes(y=met), linetype=2)

qplot(time, cent-cmt1, data=dd, log="y", geom="line") + geom_line(aes(y=met-cmt2), linetype=2)

with(dd, cent-cmt1)
