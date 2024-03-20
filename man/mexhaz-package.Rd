\name{mexhaz-package}
\alias{mexhaz-package}
\docType{package}
\title{Mixed effect parametric excess hazard models}
\description{

Fit an (excess) hazard regression model using different shapes for the
baseline hazard (Weibull, piecewise constant, exponential of a B-spline
of degree 1 to 3, exponential of a restricted cubic spline), with the
possibility to include time-dependent effects of variable(s) and a
random effect defined at the cluster level. Follow-up time can be
entered in the right-censored or counting process input style. The
latter allows the modelling of survival data with delayed entries. The
time-dependent effect of a covariable is modelled by adding interaction
terms between the covariable and a function of time of the same class as
the one used for the baseline hazard (in particular, with the same knots
for piecewise constant hazards; and with the same degree and the same
knots for B-spline or restricted cubic spline functions). The random
effect is assumed to be normally distributed with mean 0 and standard
deviation sigma. The optimisation process uses adaptive Gaussian
quadrature to calculate the cluster-specific marginal likelihoods. The
logarithm of the full marginal likelihood, defined as the sum of the
logarithms of the cluster-specific marginal likelihoods, is then
maximised using an optimisation routine such as \code{nlm} (default) or
\code{optim}. Functions to compute and plot the predicted (excess)
hazard and (net) survival are provided. In the case of a random
intercept model, three types of predictions can be obtained : marginal
(population-averaged), cluster-specific posterior (cluster-averaged), or
conditional (for a given quantile of the ditribution of the random effect).

}

\author{
Hadrien Charvat, Aurelien Belot
}

\references{Charvat H, Remontet L, Bossard N, Roche L, Dejardin O,
  Rachet B, Launoy G, Belot A; CENSUR Working Survival Group. A
  multilevel excess hazard model to estimate net survival on
  hierarchical data allowing for non-linear and non-proportional effects
  of covariates. Stat Med 2016;35(18):3066-3084 (doi: 10.1002/sim.6881)

Charvat H, Belot A. An R package for fitting flexible hazard-based
regression models for overall and excess mortality with a random
effect. J Stat Softw 2021;98(14):1-36 (doi: 10.18637/jss.v098.i14)
} \keyword{package}
%\seealso{
%}

\examples{

data(simdatn1)

## Fit of a mixed-effect excess hazard model, with the baseline hazard
## described by a Weibull distribution (without covariables)

Mod_weib_mix <- mexhaz(formula=Surv(time=timesurv,
event=vstat)~1, data=simdatn1, base="weibull",
expected="popmrate", verbose=0, random="clust")


## Examples of syntax for various models (not run)

## Fit of a fixed-effect excess hazard model, with the baseline hazard
## described by a Weibull distribution and with effects of age (agecr),
## deprivation index (depindex) and sex (IsexH) using the optim
## procedure and the BFGS method (see help of optim).

# Mod_weib <- mexhaz(formula=Surv(time=timesurv,
# event=vstat)~agecr+depindex+IsexH, data=simdatn1, base="weibull",
# expected="popmrate", verbose=1000, fnoptim="optim",
# method="BFGS")


## Fit of a mixed-effect excess hazard model, with the baseline hazard
## described by a cubic B-spline with two knots at 1 and 5 year and with
## effects of age (agecr), deprivation index (depindex) and sex (IsexH)

# Mod_bs3_2mix <- mexhaz(formula=Surv(time=timesurv,
# event=vstat)~agecr+depindex+IsexH, data=simdatn1, base="exp.bs",
# degree=3, knots=c(1,5), expected="popmrate", random="clust",
# verbose=1000)


## Fit of a fixed-effect overall hazard model, with the baseline hazard
## described by a piecewise constant function with the following vector
## of knots (defining the endpoints of the intervals on which the hazard
## is constant): (1,3,5,8), and with effects of age (agecr), deprivation
## index (depindex) and sex (IsexH)

# Mod_pw <- mexhaz(formula=Surv(time=timesurv, event=vstat)~
# agecr+depindex+IsexH, data= simdatn1, base="pw.cst", knots=c(1,3,5,8),
# verbose=1000)


## Fit of a fixed-effect excess hazard model, with the baseline hazard
## described by a cubic B-spline with two knots at 1 and 5 year and with
## effects of age (agecr), deprivation index (depindex) and sex (IsexH)

# Mod_bs3_2 <- mexhaz(formula=Surv(time=timesurv,
# event=vstat)~agecr+depindex+IsexH, data=simdatn1, base="exp.bs",
# degree=3, knots=c(1,5), expected="popmrate", verbose=1000)


## Fit of a mixed-effect excess hazard model, with the baseline hazard
## described by a cubic B-spline with two knots at 1 and 5 year and with
## effects of age (agecr), deprivation index (depindex) and sex (IsexH)

# Mod_bs3_2mix <- mexhaz(formula=Surv(time=timesurv,
# event=vstat)~agecr+depindex+IsexH, data=simdatn1, base="exp.bs",
# degree=3, knots=c(1,5), expected="popmrate", random="clust",
# verbose=1000)


## Fit of a mixed-effect excess hazard model, with the baseline hazard
## described by a cubic B-spline with two knots at 1 and 5 year, with
## effects of age (agecr), deprivation index (depindex) and sex (IsexH)
## and with a time-dependent effect for age (agecr) and sex (IsexH).

# Mod_bs3_2mixnph <- mexhaz(formula=Surv(time=timesurv,
# event=vstat)~agecr+depindex+IsexH + nph(agecr+IsexH), data=simdatn1,
# base="exp.bs", degree=3, knots=c(1,5), expected="popmrate",
# random="clust", verbose=1000)

}