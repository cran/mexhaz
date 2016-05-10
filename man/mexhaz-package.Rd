\name{mexhaz-package}
\alias{mexhaz-package}
\docType{package}
\title{Mixed effect parametric excess hazard models}
\description{
Fit an (excess) hazard regression model using different shapes for the
baseline hazard (Weibull, piecewise constant and exponential of a
B-spline), with the possibility to include time-dependent and/or
non-linear effect(s) of variable(s) and a random effect defined at the
cluster level. The time-dependent effect of a covariate is modelled by
adding interaction terms between the covariate and a function of time of
the same class as the one used for the baseline hazard (in particular,
with the same knots for piecewise constant hazards; and with the same
degree and the same knots for B-spline functions). The random effect is
assumed to be normally distributed with mean 0 and standard deviation
sigma. The optimisation process uses adaptive Gaussian quadrature to
calculate the cluster-specific marginal likelihoods. The logarithm of
the full marginal likelihood, defined as the sum of the logarithms of
the cluster-specific marginal likelihoods, is then maximised using
optimisation routines such as \code{nlm} or \code{optim}. Functions to
compute and plot the predicted (excess) hazard and (net) survival are
provided.
}

\details{
\tabular{ll}{
Package: \tab mexhaz\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2016-03-24\cr
License: \tab GPL (>=2)\cr
}
}

\author{
Hadrien Charvat, Aurelien Belot
}

\references{Charvat H, Remontet L, Bossard N, Roche L, Dejardin O,
  Rachet B, Launoy G, Belot A; CENSUR Working Survival Group. A
  multilevel excess hazard model to estimate net survival on
  hierarchical data allowing for non-linear and non-proportional effects
  of covariates. Stat Med 2016. (doi: 10.1002/sim.6881) } \keyword{package}
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
# egree=3, knots=c(1,5), expected="popmrate", verbose=1000)


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
