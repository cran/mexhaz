\name{update.mexhaz}
\alias{update.mexhaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Update of a mexhaz model}
\description{
Function allowing the user to update an existing mexhaz model. All the
arguments of the model can be updated. If the argument 'init' is not
provided, the function uses the estimated values of the existing model
as starting values for the corresponding parameters of the new model.
}

\usage{
\method{update}{mexhaz}(object, formula, data, expected = NULL, base = NULL, 
degree = 3, knots = NULL, bound = NULL, n.gleg = 20, init = NULL,
random = NULL, n.aghq = 10, fnoptim = c("nlm", "optim"), 
verbose = 0, method = "Nelder-Mead", iterlim = 10000, numHess = FALSE,
print.level = 1, exactGradHess = TRUE, gradtol = 
ifelse(exactGradHess, 1e-8, 1e-6), envir = parent.frame(), ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{mexhaz}, corresponding to a survival model fitted with the \code{mexhaz} function.
}
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

See \code{\link{update.formula}} for more details.}
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
exponential of a restricted cubic B-spline (also called 'natural
spline'), \code{"pw.cst"} for a piecewise constant hazard. By default,
\code{base=NULL}.  }

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
  If \code{base="exp.bs"} or \code{base="exp.ns"}, computation of the
  B-spline basis requires that boundary knots be given. The \code{bound}
  argument allows the user to specify these boundary knots. If
  \code{base="exp.bs"}, the interval defined by the boundary knots must
  at least include the interval c(0,max(time)) (otherwise, there
  could be problems with ill-conditioned bases). If
  \code{base="exp.ns"}, the boundary knots correspond to the knots
  beyond which the spline is contrained to be linear (in that case, the
  boundary knots can be values contained in c(0,max(time))). By
  default, the boundary knots are set to
  c(0,max(time)).
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
    other parameters set to 0. In case of failed convergence, several
    trials are run with an adaptation of the value of the intercept.

    if \code{exactGradHess=FALSE}, the following values are used:

    - if \code{base="weibull"}, the scale and shape parameters are set
    to 0.1;

    - if \code{base="exp.bs"}, the parameters of the B-spline are all
    set to -1;

    - if \code{base="exp.ns"}, the parameters of the restricted cubic B-spline are all
    set to -1;

    - if \code{base="pw.cst"}, the logarithm of the piecewise-constant
    hazards are set to -1;

     the parameters describing the effects of the covariates are all set
     to 0;

    the parameter representing the standard deviation of the random
    effect is set to 0.1. (if \code{exactGradHess=TRUE}, several trials are run with an adaptation of the value in case of failed convergence).
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
    likelihood. Selection can be made between \code{"nlm"} (by default) and \code{"optim"}. Note: if \code{exactGradHess=TRUE}, this
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
    on the value of the argument \code{exactGradHess}.}
  \item{envir}{environment in which the objects' names given as arguments to the updated
model are to be found.}
  \item{\dots}{
represents additional parameters directly passed to \code{nlm} or
\code{optim} to control the optimisation process.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  An object of class \code{mexhaz}. See \code{\link{mexhaz}} for more details.
}
\references{Charvat H, Remontet L, Bossard N, Roche L, Dejardin O,
  Rachet B, Launoy G, Belot A; CENSUR Working Survival Group. A
  multilevel excess hazard model to estimate net survival on
  hierarchical data allowing for non-linear and non-proportional effects
  of covariates. Stat Med 2016;35:3066-3084 (doi: 10.1002/sim.6881)}
\author{
Hadrien Charvat, Aurelien Belot
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mexhaz}}
}

\examples{

data(simdatn1)

## Fit of a mixed-effect excess hazard model, with the baseline hazard
## described by a Weibull distribution (without covariables)

Mod_weib <- mexhaz(formula=Surv(time=timesurv,
event=vstat)~1, data=simdatn1, base="weibull", verbose=0)

## Add an effect of gender
Mod_weib_2 <- update(Mod_weib, formula=~.+IsexH)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{models}
