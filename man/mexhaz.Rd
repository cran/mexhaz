\name{mexhaz}
\alias{mexhaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mexhaz function}
\description{

Fit an (excess) hazard regression model using different shapes for the
baseline hazard (Weibull, piecewise constant, exponential of a B-spline
of degree 1 to 3, exponential of a restricted cubic spline), with the
possibility to include time-dependent effects of variable(s) and a
random effect defined at the cluster level. The function accepts
right-censored and counting process input styles for the follow-up
time. The latter allows the modelling of survival data with delayed
entries. The time-dependent effect of a covariable is modelled by adding
interaction terms between the covariable and a function of time of the
same class as the one used for the baseline hazard (in particular, with
the same knots for piecewise constant hazards; and with the same degree
and the same knots for B-spline or restricted cubic spline
functions). The random effect is assumed to be normally distributed with
mean 0 and standard deviation sigma. The optimisation process uses
adaptive Gaussian quadrature to calculate the cluster-specific marginal
likelihoods. The logarithm of the full marginal likelihood, defined as
the sum of the logarithms of the cluster-specific marginal likelihoods,
is then maximised using an optimisation routine such as \code{nlm}
(default) or \code{optim}.

}

\usage{
mexhaz(formula, data, expected = NULL, base = c("weibull",
"exp.bs", "exp.ns", "pw.cst"), degree = 3, knots = NULL,
bound = NULL, n.gleg = 20, init = NULL, random = NULL,
n.aghq = 10, fnoptim = c("nlm", "optim"), verbose = 0,
method = "Nelder-Mead", iterlim = 10000, numHess = FALSE,
print.level = 1, exactGradHess = TRUE, gradtol =
ifelse(exactGradHess, 1e-8, 1e-6), testInit = TRUE, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a formula object, with the response on the left of the \code{~}
operator, and the linear predictor on the right. The response must be of
the form \code{Surv(time, event)} for right censored data or
\code{Surv(time, time2, event)} for counting process style data. The linear predictor accepts a
special instruction \code{nph()} for specifying variables for which a
time-dependent effect should be modelled (if several variables are
modelled with time-dependent effects, separate these variables inside the
\code{nph()} instruction with a \code{+} sign).

In case \code{time} takes the value 0 for some observations when data
are entered in the right censored style, it is
assumed that these observations refer to events (or censoring) that occurred
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

\item{base}{ functional form that should be used to model the baseline
hazard. Selection can be made between the following options:
\code{"weibull"} for a Weibull hazard, \code{"exp.bs"} for a hazard
described by the exponential of a B-spline (only B-splines of degree 1,
2 or 3 are accepted), \code{"exp.ns"} for a hazard described by the
exponential of a restricted cubic spline (also called 'natural
spline'), \code{"pw.cst"} for a piecewise constant hazard. By default,
\code{base="weibull"}.

For the Weibull hazard model, the cumulative hazard is given by the
following relationship:

H(t,x,z) = lambda*exp(x'b)*t^(rho*exp(z'g))

where lambda and rho are the parameters of the Weibull baseline hazard,
x represent variables with proportional effect (with corresponding
regression coefficients 'b') and z represent variables with
time-dependent effects (with corresponding regression coefficients
'g'). The \code{mexhaz()} function does not estimate directly lambda and
rho but their logarithms (in the output of the function, these are named
respectively 'logLambda' and 'logRho').

For the spline-based hazards, it should be noted that the B-spline
and restricted cubic spline bases created internally in the
\code{mexhaz()} function are identical to those obtained by the use of,
respectively, the \code{bs()} and \code{ns()} functions from the
\code{splines} package.}

\item{degree}{
  if \code{base="exp.bs"}, \code{degree} represents the
 degree of the B-spline used. Only integer values between 1 and 3 are
 accepted, and 3 is the default.
}

\item{knots}{
if \code{base="exp.bs"} or \code{"exp.ns"}, \code{knots} is the vector of interior knots of
the spline. If \code{base="pw.cst"}, \code{knots} is the vector
defining the endpoints of the time intervals on which the hazard is
assumed to be constant. By default, \code{knots=NULL} (that is, it
produces a B-spline with no interior knots if \code{base="exp.bs"}, a
linear B-spline with no interior knots if \code{base="exp.ns"}, or a
constant hazard over the whole follow-up period if
\code{base="pw.cst"}).
}
\item{bound}{
  a vector of two numerical values corresponding to the boundary knots
  of the spline functions. If \code{base="exp.bs"} or
  \code{base="exp.ns"}, computation of the B-spline basis requires that
  boundary knots be given. The \code{bound} argument allows the user to
  specify these boundary knots. If \code{base="exp.bs"}, the interval
  defined by the boundary knots must at least include the interval
  c(0,max(time)) (otherwise, there could be problems with
  ill-conditioned bases). If \code{base="exp.ns"}, the boundary knots
  correspond to the knots beyond which the spline is contrained to be
  linear (in that case, the boundary knots can be values contained in
  c(0,max(time))). By default, the
  boundary knots are set to c(0,max(time)).
}
  \item{n.gleg}{
if \code{base="exp.bs"} and degree is equal to 2 or 3, or if
  \code{base="exp.ns"}, the cumulative
hazard is computed via Gauss-Legendre quadrature and \code{n.gleg} is the number
of quadrature nodes to be used to compute the cumulative hazard. By
default, \code{n.gleg=20}.
}
  \item{init}{vector of initial values. By default \code{init=NULL}
    and the initial values are internally set to the following values:

    for the baseline hazard:

    if \code{exactGradHess=TRUE} (except for the excess hazard random
    effect model for which this argument is ignored), the intercept is
    set to 0.5*log(Number of events/Person-years of observation) and all
    other parameters set to 0. In case of failed convergence, and if the
    \code{testInit} argument is set to \code{TRUE} (default), several
    trials are run with an adaptation of the value of the intercept.

    if \code{exactGradHess=FALSE}, the following values are used:
    
    - if \code{base="weibull"}, the logarithm of the scale and shape parameters is set to 0.1;

    - if \code{base="exp.bs"}, the parameters of the B-spline are all
    set to -1;

    - if \code{base="exp.ns"}, the parameters of the restricted cubic spline are all
    set to -1;

    - if \code{base="pw.cst"}, the logarithm of the piecewise-constant
    hazards are set to -1;

    the parameters describing the effects of the covariables are all set
    to 0;

    the parameter representing the standard deviation of the random
    effect is set to 0.1. In case of failed convergence, if
    \code{exactGradHess=TRUE} and if \code{testInit=TRUE}, several
    trials are run with an adaptation of the value of this random effect
    parameter.

  }
  \item{random}{
name of the variable to be entered as a random effect (must be given
between quotes), representing the cluster membership. By default,
\code{random=NULL} which means that the function fits a fixed effects model.
}
  \item{n.aghq}{
number of quadrature points used for estimating the
cluster-specific marginal likelihoods by adaptive Gauss-Hermite
quadrature. By default, \code{n.aghq=10}.
}
  \item{fnoptim}{
name of the R optimisation procedure used to maximise the
    likelihood. Selection can be made between \code{"nlm"} (by default)
    and \code{"optim"}. Note: if \code{exactGradHess=TRUE}, this
    argument will be ignored (\code{fnoptim} will be set automatically to \code{"nlm"}).
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
\code{verbose=0} (default), nothing is displayed.
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
  \item{numHess}{
logical value allowing the user to choose between the Hessian returned
    by the optimization algorithm (default) or the Hessian estimated by
    the \code{hessian} function from the \code{numDeriv} package. The
    latter might be more accurate but its estimation is more
    time-consuming. We suggest to use the default Hessian estimation
    procedure during model building and estimate the \code{numDeriv}-based
    Hessian only on the final model. Note: if \code{exactGradHess=TRUE}, this argument is ignored.
}
  \item{print.level}{
    this argument is only used if \code{fnoptim="nlm"}. It determines the
    level of printing during the optimisation process. The default value
    (for the \code{mexhaz} function) is set to '1' which means that
    details on the initial and final step of the optimisation procedure
    are printed (see the help page of \code{nlm} for further details).
  }
  \item{exactGradHess}{ logical value allowing the user to decide
    whether maximisation of the likelihood should be based on the
    analytic gradient and Hessian computed internally (default,
    corresponding to \code{exactGradHess=TRUE}).  In that case,
    optimisation is performed with the \code{nlm} function. Note: even
    if set to \code{TRUE}, this argument is ignored when the user wants
    to fit an excess hazard model including a random effect because in
    that case, there is no simple way to obtain the analytic gradient
    and Hessian.}
  \item{gradtol}{
    this argument is only used if \code{fnoptim="nlm"}. It corresponds
    to the tolerance at which the scaled gradient is considered close
    enough to zero to terminate the algorithm. The default value depends
    on the value of the argument \code{exactGradHess}. 
  }
  \item{testInit}{ this argument is used only when
    \code{exactGradHess=TRUE} and when the model is not an excess hazard
    random effect model. It instructs the \code{mexhaz} function to try
    several vectors of initial values in case optimization was not
    successful with the default (or user-defined) initial
    values. Because optimization based on the analytical gradient and
    Hessian is usually fast, this simple and empirical procedure proves
    useful to increase the probability of convergence in cases when it
    is difficult to specify appropriate initial values. Only the initial
    values for the intercept and for the parameter corresponding to the
    random effect are modified.}
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
  \item{expected}{name of the variable corresponding to the expected
    hazard (takes the value \code{NA} for standard, i.e., 'non-excess' hazard models).}
  \item{xlevels}{information concerning the levels of the categorical
    variables used in the model (used by \code{predict.mexhaz}).}
  \item{n.obs.tot}{total number of observations in the dataset.}
  \item{n.obs}{number of observations used to fit the model (after
    exclusion of missing values).}
  \item{n.events}{number of events (after exclusion of missing values).}
  \item{n.clust}{number of clusters.}
  \item{n.time.0}{number of observations for which the observed
    follow-up time was equal to 0 (only for right censored type data).}
  \item{base}{function used to model the baseline hazard.}
  \item{max.time}{maximal observed time in the dataset.}
  \item{boundary.knots}{vector of boundary values used to define the B-spline (or
    natural spline) bases.}
  \item{degree}{degree of the B-spline used to model the logarithm of
    the baseline hazard.}
  \item{knots}{vector of interior knots used to define the B-spline (or
    natural spline) bases.}
  \item{names.ph}{names of the covariables with a proportional effect.}
  \item{random}{name of the variable defining cluster membership (set to
    \code{NA} in the case of a purely fixed effects model).}
  \item{init}{a vector containing the initial values of the parameters.}
  \item{coefficients}{a vector containing the parameter estimates.}
  \item{std.errors}{a vector containing the standard errors of the
    parameter estimates.}
  \item{vcov}{the variance-covariance matrix of the estimated parameters.}
  \item{gradient}{the gradient of the log-likelihood function evaluated at the estimated parameters.}
  \item{hessian}{the Hessian of the log-likelihood function evaluated at the estimated parameters.}
  \item{mu.hat}{a \code{data.frame} containing the estimated cluster-specific
    random effects (shrinkage estimators).}
  \item{var.mu.hat}{the covariance matrix of the
    cluster-specific shrinkage estimators.}
  \item{vcov.fix.mu.hat}{a \code{matrix} containing the covariances
    between the fixed effect and the
    cluster-specific shrinkage estimators. More specifically, the i-th line
    of the matrix represents the covariances between the shrinkage
    estimator of the i-th cluster and the fixed effect estimates. This
    matrix is used by the function \code{predict.mexhaz} to make
    cluster-specific predictions.}
  \item{n.par}{number of estimated parameters.}
  \item{n.gleg}{number of Gauss-Legendre quadrature points used to
    calculate the cumulative (excess) hazard (only relevant if a
    B-spline of degree 2 or 3 or a cubic restricted spline was used to
    model the logarithm of the baseline hazard).}
  \item{n.aghq}{number of adaptive Gauss-Hermite quadrature points used to
    calculate the cluster-specific marginal likelihoods (only relevant
    if a multi-level model is fitted).}
  \item{fnoptim}{name of the R optimisation procedure used to maximise the likelihood.}
  \item{method}{optimisation method used by \code{optim}.}
  \item{code}{code (integer) indicating the status of the optimisation
    process (this code has a different meaning for \code{nlm} and for
    \code{optim}: see their respective help page for details).}
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
  of covariates. Stat Med 2016;35(18):3066-3084 (doi: 10.1002/sim.6881)

Charvat H, Belot A. An R package for fitting flexible hazard-based
regression models for overall and excess mortality with a random
effect. J Stat Softw 2021;98(14):1-36 (doi: 10.18637/jss.v098.i14)
}
\author{
Hadrien Charvat, Aurelien Belot
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fixef.mexhaz}}, \code{\link{predict.mexhaz}}, \code{\link{print.mexhaz}},
\code{\link{ranef.mexhaz}}, \code{\link{summary.mexhaz}}, \code{\link{update.mexhaz}}, \code{\link{vcov.mexhaz}}
}

\examples{

data(simdatn1)

## Fit of a mixed-effect excess hazard model, with the baseline hazard
## described by a Weibull distribution (without covariables)

Mod_weib_mix <- mexhaz(formula=Surv(time=timesurv,
event=vstat)~1, data=simdatn1, base="weibull",
expected="popmrate", verbose=0, random="clust")


## More complex examples (not run)

## Fit of a mixed-effect excess hazard model, with the baseline hazard
## described by a cubic B-spline with two knots at 1 and 5 year and with
## effects of age (agecr), deprivation index (depindex) and sex (IsexH)

# Mod_bs3_2mix_nph <- mexhaz(formula=Surv(time=timesurv,
# event=vstat)~agecr+depindex+IsexH+nph(agecr), data=simdatn1,
# base="exp.bs", degree=3, knots=c(1,5), expected="popmrate",
# random="clust", verbose=1000)

## Fit of a mixed-effect excess hazard model, with the baseline hazard
## described by a restricted cubic spline with two knots at the
## tertiles of the distribution of follow-up times for individuals who
## presented the event and with effects of age (agecr) and sex (IsexH)

# Mod_ns3_2mix_nph <- mexhaz(formula=Surv(time=timesurv,
# event=vstat)~agecr+nph(agecr), data=simdatn1, base="exp.ns", degree=3,
# knots=quantile(simdatn1[simdatn1$vstat==1,]$timesurv, probs=c(1:2/3)),
# expected="popmrate", random="clust", verbose=1000)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{models}
