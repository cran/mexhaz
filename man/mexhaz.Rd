\name{mexhaz}
\alias{mexhaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mexhaz function}
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
optimisation routine such as \code{nlm} or \code{optim}.

}

\usage{mexhaz(formula, data, expected = NULL, base = c("weibull",
"exp.bs", "pw.cst"), degree = 3, knots = NULL, bo.max = NULL,
n.gleg = 20, init = NULL, random = NULL, n.aghq = 10,
fnoptim = c("nlm", "optim"), verbose = 100, method = "Nelder-Mead",
iterlim=10000,print.level=1,...)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a formula object, with the response on the left of the \code{~}
operator, and the linear predictor on the right. The response must be of
the form \code{Surv(time, event)}. The linear predictor accepts a
special instruction \code{nph()} for specifying variables for which a
time-dependent effect should be modelled (if several variables are
modelled with time-dependent effects, separate these variables inside the
\code{nph()} instruction with a \code{+} sign).

In case \code{time} takes the value 0 for some observations, it is
assumed that these observations refer to events/censoring that occurred
on the first day of follow-up. Consequently, a value of 1/730.5 (half a
day) is substituted in order to make computations possible. However, it
should be stressed that this is just a convention and that it does not
make much sense if the time scale is not expressed in years. We therefore
advise the analyst to deal with 0 time values during the dataset preparation stage.
}
  \item{data}{
a \code{data.frame} containing the variables referred to in the \code{formula},
as well as in the \code{expected} and \code{random} arguments if these arguments are used.
}
  \item{expected}{
name of the variable (must be given in quotes) representing the
population (i.e., expected) hazard. By default, \code{expected=NULL}, which
means that the function estimates the overall hazard (and not the excess
hazard).
}
  \item{base}{
functional form that should be used to model the baseline
hazard. Selection can be made between the following options: \code{"weibull"}
for a Weibull hazard, \code{"exp.bs"} for a hazard described by the exponential
of a B-spline (only B-splines of degree 1, 2 or 3 are accepted),
\code{"pw.cst"} for a piecewise constant hazard. By default, \code{base="weibull"}.
}
  \item{degree}{
 if \code{base="exp.bs"}, \code{degree} represents the degree of the B-spline used. Only integer values between 1 and 3 are accepted, and 3 is the default.
}
\item{knots}{
if \code{base="exp.bs"}, \code{knots} is the vector of interior knots of
the B-spline. If \code{base="pw.cst"}, \code{knots} is the vector
defining the endpoints of the time intervals on which the hazard is
assumed to be constant. By default, \code{knots=NULL} (that is, it
produces a B-spline with no interior knots if \code{base="exp.bs"} or a
constant hazard over the whole follow-up period if
\code{base="pw.cst"}).
}
\item{bo.max}{
  
  if \code{base="exp.bs"}, computation of the B-spline basis requires
  that boundary knots be given. By default, these are set to
  \code{c(0,max(time))} . Provided that it is equal or greater than
  \code{max(time)} (where \code{time} is the time variable defined in
  the \code{Surv()} formula), the upper boundary knot can
  (theoretically) be set to any value, specified by \code{bo.max}. Using
  different values of \code{bo.max} will result in models with different
  estimated values of the parameters corresponding to the B-spline
  basis.  However, the resulting baseline hazard as well as the
  proportional effects of covariables will be almost identical (up to
  numerical approximations). By default, \code{bo.max=NULL} and the
  B-spline boundary knots are set to \code{c(0,max(time))}.
   
}
  \item{n.gleg}{
if \code{base="exp.bs"} and degree is equal to 2 or 3, the cumulative
hazard is computed via Gauss-Legendre quadrature and \code{n.gleg} is the number
of quadrature nodes to be used to compute the cumulative hazard. By
default, \code{n.gleg=20}.
}
  \item{init}{vector of initial values. By default \code{init=NULL}
    and the initial values are internally set to the following values:

    for the baseline hazard:

    if \code{base="weibull"}, the scale and shape parameters are set
    to 0.1;

    if \code{base="exp.bs"}, the parameters of the B-spline are all
    set to -1;

    if \code{base="pw.cst"}, the logarithm of the piecewise-constant
    hazards are set to -1;

     the parameters describing the effects of the covariates are all set
     to 0;

    the parameter representing the standard deviation of the random
    effect is set to 0.1.
}
  \item{random}{
name of the variable to be entered as a random effect (must be given
between quotes), representing the cluster membership. By default,
\code{random=NULL} which means that the function fits a fixed effects model.
}
  \item{n.aghq}{
number of quadrature points to be used for estimating the
cluster-specific marginal likelihoods by adaptive Gauss-Hermite
quadrature. By default, \code{n.aghq=10}.
}
  \item{fnoptim}{
name of the R optimisation procedure used to maximise the
    likelihood. Selection can be made between \code{"nlm"} (by default) and \code{"optim"}.
}
  \item{verbose}{
integer parameter representing the frequency at which the current state
of the optimisation process is displayed. Internally, an 'evaluation' is
defined as an estimation of the log-likelihood for a given vector of
parameters. This means that the number of evaluations is increased each
time the optimisation procedure updates the value of any of the
parameters to be estimated. If \code{verbose=n} (with \code{n} an integer),
the function will display the current values of the parameters, the
log-likelihood and the time elapsed every \code{n} evaluations. If
\code{verbose=0}, nothing is displayed.
}
  \item{method}{
if \code{fnoptim="optim"}, \code{method} represents the optimisation
    method to be used by \code{optim}. By
    default, \code{method="Nelder-Mead"}. This parameter is not used if \code{fnoptim="nlm"}.
}
  \item{iterlim}{
if \code{fnoptim="nlm"}, \code{iterlim} represents the maximum number of
    iterations before the \code{nlm} optimisation
    procedure is terminated. By default, \code{iterlim} is set to
    10000. This parameter is not used if \code{fnoptim="optim"} (in this
    case, the maximum number of iterations must be given as part of a
    list of control parameters via the \code{control} argument: see the
    help page of \code{optim} for further details).
  }
  \item{print.level}{
    this argument is only used if \code{fnoptim="nlm"}. It determines the
    level of printing during the optimisation process. The default value
    (for the \code{mexhaz} function) is set to '1' which means that
    details on the initial and final step of the optimisation procedure
    are printed (see the help page of \code{nlm} for further details).
    }
  \item{\dots}{
represents additional parameters directly passed to \code{nlm} or
\code{optim} to control the optimisation process.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  An object of class \code{mexhaz} containing the following elements:
  \item{dataset}{name of the dataset used to fit the model.}
  \item{call}{function call on which the model is based.}
  \item{formula}{formula part of the call.}
  \item{xlevels}{information concerning the levels of the categorical
    variables used in the model (used by \code{predMexhaz}).}
  \item{n.obs.tot}{total number of observations in the dataset.}
  \item{n.obs}{number of observations used to fit the model (after
    exclusion of missing values).}
  \item{n.events}{number of events (after exclusion of missing values).}
  \item{n.clust}{number of clusters.}
  \item{n.time.0}{number of observations for which the observed
    follow-up time was equal to 0.}
  \item{base}{function used to model the baseline hazard.}
  \item{max.time}{maximal observed time in the dataset.}
  \item{bounds}{vector of boundary values used to define the B-spline basis.}
  \item{degree}{degree of the B-spline used to model the logarithm of
    the baseline hazard.}
  \item{knots}{vector of interior knots used to define the B-spline
    basis.}
  \item{names.ph}{names of the covariables with a proportional effect.}
  \item{random}{name of the variable defining cluster membership (set to
    \code{NA} in the case of a purely fixed effects model).}
  \item{coefficients}{a vector containing the parameter estimates.}
  \item{std.errors}{a vector containing the standard errors of the
    parameter estimates.}
  \item{vcov}{the variance-covariance matrix of the estimated parameters.}
  \item{mu.hat}{a \code{data.frame} with the shrinkage estimates predicted
    for each cluster.}
  \item{n.par}{number of estimated parameters.}
  \item{n.gleg}{number of Gauss-Legendre quadrature points used to
    calculate the cumulative (excess) hazard (only relevant if a
    B-spline of degree 2 or 3 was used to model the logarithm of the
    baseline hazard).}
  \item{n.aghq}{number of adaptive Gauss-Hermite quadrature points used to
    calculate the cluster-specific marginal likelihoods (only relevant
    if a multi-level model is fitted).}
  \item{fnoptim}{name of the R optimisation procedure used to maximise the likelihood.}
  \item{method}{optimisation method used by \code{optim}.}
  \item{code}{code (integer) indicating the status of the optimisation
    process (this code has a different meaning for \code{nlm} and for
    \code{optim}).}
  \item{loglik}{value of the log-likelihood at the end of the
    optimisation procedure.}
  \item{iter}{number of iterations used in the
    optimisation process.}
  \item{eval}{number of evaluations used in the
    optimisation process.}
 \item{time.elapsed}{total time required to reach convergence.}
}
\references{Charvat H, Remontet L, Bossard N, Roche L, Dejardin O,
  Rachet B, Launoy G, Belot A; CENSUR Working Survival Group. A
  multilevel excess hazard model to estimate net survival on
  hierarchical data allowing for non-linear and non-proportional effects
  of covariates. Stat Med 2016. (doi: 10.1002/sim.6881)}
\author{
Hadrien Charvat, Aurelien Belot
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{print.mexhaz}}, \code{\link{summary.mexhaz}}, \code{\link{predMexhaz}}
}

\examples{

data(simdatn1)

## Fit of a mixed-effect excess hazard model, with the baseline hazard
## described by a Weibull distribution (without covariables)

Mod_weib_mix <- mexhaz(formula=Surv(time=timesurv,
event=vstat)~1, data=simdatn1, base="weibull",
expected="popmrate", verbose=0, random="clust")


## A more complex example (not run)

## Fit of a mixed-effect excess hazard model, with the baseline hazard
## described by a cubic B-spline with two knots at 1 and 5 year and with
## effects of age (agecr), deprivation index (depindex) and sex (IsexH)

# Mod_bs3_2mix_nph <- mexhaz(formula=Surv(time=timesurv,
# event=vstat)~agecr+depindex+IsexH+nph(agecr), data=simdatn1,
# base="exp.bs", degree=3, knots=c(1,5), expected="popmrate",
# random="clust", verbose=1000)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{models}
