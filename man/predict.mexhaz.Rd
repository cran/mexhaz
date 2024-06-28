\name{predict.mexhaz}
\alias{predict.mexhaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predictions based on a mexhaz model}
\description{Function for predicting the (excess) hazard and the
corresponding (net) survival from a model fitted with the \code{mexhaz}
function for a particular vector of covariates. If the survival model
was fitted with an expected hazard (excess hazard model), the estimates
obtained are excess hazard and net survival estimates. Corresponding
variance estimates are based on the Delta Method or Monte Carlo
simulation (based on the assumption of multivariate normality of the
model parameter estimates). This function allows the computation of the
hazard and the survival at one time point for several vectors of
covariates or for one vector of covariates at several time
points. When the model includes a random effect, three types of
predictions can be made: (i) marginal predictions (obtained by
integration over the random effect distribution), (ii) cluster-specific
posterior predictions for an existing cluster, or (iii) conditional
predictions for a given quantile of the random effect distribution (by
default, for the median value, that is, 0).}

\usage{
\method{predict}{mexhaz}(object, time.pts, data.val = data.frame(.NotUsed = NA),
marginal = FALSE, quant.rdm = 0.5, cluster = NULL, conf.int = c("delta", "simul", "none"),
level = 0.95, delta.type.h = c("log", "plain"), delta.type.s = c("log-log",
"log", "plain"), nb.sim = 10000, keep.sim = FALSE, include.gradient =
FALSE, dataset = NULL, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{mexhaz}, corresponding to a
    hazard-based regression model fitted with the \code{mexhaz} function.}
  \item{time.pts}{a vector of numerical values representing the time
    points at which predictions are requested. Time values greater than
    the maximum follow-up time on which the model estimation was based
    are discarded.}
  \item{data.val}{a \code{data.frame} containing the values of the
    covariates at which predictions should be calculated.}
  \item{marginal}{logical value controlling the type of predictions
    returned by the function when the model includes a random
    intercept. When \code{TRUE}, marginal predictions are computed. The
    marginal survival is obtained by integrating the predicted survival over the
    distribution of the random effect. The marginal hazard rate is
    obtained as the opposite of the marginal time
    derivative of the survival divided by the marginal survival. When
    \code{FALSE} (default value), cluster-specific posterior predictions
    or conditional predictions are calculated depending on the
    value of the \code{cluster} argument.}
  \item{quant.rdm}{numerical value (between 0 and 1) specifying the
    quantile of the random effect distribution that should be used when
    requesting conditional predictions. The default value is set to 0.5
    (corresponding to the median, that is a value of the random effect
    of 0). This argument is ignored if the model is a fixed effect
    model, when the \code{marginal} argument is set to \code{TRUE}, or
    the \code{cluster} argument is not \code{NULL}.}
  \item{cluster}{a single value corresponding to the name of the cluster for
    which posterior predictions should be calculated. These predictions
    are obtained by integrating over the cluster-specific posterior
    distribution of the random effect and thus require the original
    dataset. The dataset can either be provided as part of the
    \code{mexhaz} object given as argument or by specifying the name of
    the dataset in the \code{dataset} argument (see below). The \code{cluster} argument is not used if the model is a fixed effect
    model. The default value is \code{NULL}: this corresponds to
    marginal predictions (if \code{marginal} is set to \code{TRUE}, the preferred
    option), or to conditional predictions for a given quantile (by
    default, the median) of the distribution of the random
    effect (if \code{marginal} is set to \code{FALSE}).}
  \item{conf.int}{method to be used to compute confidence
    limits. Selection can be made between the following options:
    \code{"delta"} for the Delta Method (default value);
    \code{"simul"} for Monte Carlo simulations (can be time-consuming,
    especially for models using B-splines for the logarithm of the
    baseline hazard);
    \code{"none"} indicates absence of confidence limits estimation.}
  \item{level}{a number in (0,1) specifying the level of confidence for
  computing the confidence intervals of the hazard and the survival. The
  default value is set to \code{0.95}.}
    \item{delta.type.h}{type of confidence limits for the hazard when
      using the Delta Method. With the default value (\code{"log"}), the
      confidence limits are based on a Wald-type confidence interval for
      the logarithm of the hazard, otherwise they are based directly on
      a Wald-type CI for the hazard.}
    \item{delta.type.s}{type of confidence limits for the survival when
      using the Delta Method. With the default value (\code{"log-log"}),
      the confidence limits are based on a Wald-type confidence interval
      for the logarithm of the cumulative hazard; when the argument is
      set to \code{"log"}, they are based on a Wald-type CI for the
      cumulative hazard; otherwise they are based directly on a
      Wald-type CI for the survival.}
  \item{nb.sim}{integer value representing the number of simulations
  used to estimate the confidence limits for the (excess) hazard
  and the (net) survival. This argument is used only if
  \code{conf.int} is set to \code{"simul"}.}
  \item{keep.sim}{logical value determining if the simulated hazard and
    survival values should be returned (only used when
    \code{conf.int} is set to \code{"simul"}). These simulated values can be used by the
    \code{riskfunc} function to compute simulation-based confidence
    intervals for hazard ratios / risk ratios. The default value is set
    to \code{FALSE}.}
  \item{include.gradient}{logical value allowing the function to return
  the components of the gradient of the logarithm of the hazard and of
  the logarithm of the cumulative hazard for each prediction. This
  argument is used only if \code{conf.int} is set to \code{"delta"}. The
  default value is \code{FALSE}.}
  \item{dataset}{original dataset used to fit the \code{mexhaz} object
    given as argument to the function. This argument is only necessary
    if cluster-specific posterior predictions are requested (and if the
    dataset is not already provided in the \code{mexhaz} object). The
    default value is set to \code{NULL}.}
\item{\dots}{for potential additional parameters.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  An object of class \code{predMexhaz} that can be used by
  the function \code{plot.predMexhaz} to produce graphics of the (excess) hazard and
  the (net) survival. It contains the following elements:
  \item{call}{the \code{mexhaz} function call on which the model is based.}
\item{results}{a \code{data.frame} consisting of: the time
  points at which the (excess) hazard and the (net) survival have been
  calculated; the values of the covariates used to estimate the
  (excess) hazard and the (net) survival; the name of the cluster (when
  cluster-specific posterior predictions are requested), or the quantile
  and its value (when conditional predictions are requested); the
  (excess) hazard values with their confidence limits; and the (net) survival values with their confidence limits.}
\item{variances}{a \code{data.frame} consisting of two columns: the
  variance of the logarithm of the (excess) hazard and the variance of
  the (excess) cumulative hazard for each time points or each
  vector of covariates. The object \code{variances} is produced only when
  \code{conf.int} is set to \code{"delta"}.}
\item{grad.loghaz}{a \code{data.frame} consisting of the components of
  the gradient of the logarithm of the (excess) hazard for each time points or each
  vector of covariates. The number of columns corresponds to the number
  of model parameters. This object is produced only when
  \code{conf.int} is set to \code{"delta"} and \code{include.gradient}
  to \code{TRUE}.}
\item{grad.logcum}{a \code{data.frame} consisting of the components of
  the gradient of the logarithm of the (excess) cumulative hazard for each time points or each
  vector of covariates. The number of columns corresponds to the number
  of model parameters. This object is produced only when
  \code{conf.int} is set to \code{"delta"} and \code{include.gradient}
  to \code{TRUE}.}
\item{vcov}{a matrix corresponding to the covariance matrix used to
  compute the confidence intervals.}
\item{type}{this item can take the value
    \code{multitime} (computation of the hazard and the survival at
  at several time points for one vector of covariates) or
    \code{multiobs} (computation of the hazard and the survival at
  at one time point for several vectors of covariates). It is
  used by \code{plot.predMexhaz} and \code{points.predMexhaz}.}
\item{type.me}{the type of predictions produced in case of a model
  including a random intercept. Can take the values \code{marginal}
  (marginal predictions), \code{posterior} (cluster-specific posterior
  predictions), or \code{conditional} (conditional prediction for a
  given quantile of the distribution of the random intercept).}
\item{ci.method}{the method used to compute confidence limits.}
\item{level}{level of confidence used to compute confidence limits.}
\item{delta.type}{type of confidence limits for the hazard and the
  survival when using the Delta Method.}
\item{nb.sim}{number of simulations used to estimate the confidence
  limits when \code{ci.method} is set to \code{"simul"}.}
\item{sim.haz}{matrix containing the simulated hazards (each column
  representing a simulated vector of values); only returned when
  \code{keep.sim} is set to \code{TRUE}.}
\item{sim.surv}{matrix containing the simulated survival probabilities (each column
  representing a simulated vector of values); only returned when
  \code{keep.sim} is set to \code{TRUE}.}
}
\references{Charvat H, Remontet L, Bossard N, Roche L, Dejardin O,
  Rachet B, Launoy G, Belot A; CENSUR Working Survival Group. A
  multilevel excess hazard model to estimate net survival on
  hierarchical data allowing for non-linear and non-proportional effects
  of covariates. Stat Med 2016;35:3066-3084 (doi: 10.1002/sim.6881).

  Skrondal A, Rabe-Hesketh S. Prediction in multilevel generalized
  linear models. J R Stat Soc A Stat Soc 2009;172(3):659-687 (doi: 10.1111/j.1467-985X.2009.00587.x).
}
\author{
Hadrien Charvat, Aurelien Belot
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{print.predMexhaz}}, \code{\link{plot.predMexhaz}},
  \code{\link{points.predMexhaz}}, \code{\link{lines.predMexhaz}}
}

\examples{

data(simdatn1)

## Fit of a fixed-effect hazard model, with the baseline hazard
## described by a linear B-spline with two knots at 1 and 5 year and
## with effects of age (agecr), deprivation index (depindex) and sex
## (IsexH)

Mod_bs1_2 <- mexhaz(formula=Surv(time=timesurv,
event=vstat)~agecr+depindex+IsexH, data=simdatn1, base="exp.bs",
degree=1, knots=c(1,5), verbose=0)

## Prediction at several time points for one vector of covariates
Pred_Modbs1_2A <- predict(Mod_bs1_2, time.pts=seq(0.1,10,by=0.1),
data.val=data.frame(agecr=0,depindex=0.5,IsexH=1))

## Prediction for several vectors of covariates at one time point
Pred_Modbs1_2B <- predict(Mod_bs1_2, time.pts=10,
data.val=data.frame(agecr=c(-0.2,-0.1,0), depindex=c(0.5,0.5,0.5),
IsexH=c(1,1,1)))

## Prediction for all individuals of the study population
## at one time point

Pred_Modbs1_2C <- predict(Mod_bs1_2, time.pts=10,
data.val=simdatn1)


# Example of cluster-specific posterior prediction (not run)

## Fit of a mixed-effect excess hazard model, with the baseline hazard
## described by a cubic B-spline with two knots at 1 and 5 year

# Mod_bs3_2mix <- mexhaz(formula=Surv(time=timesurv,
# event=vstat)~agecr+IsexH, data=simdatn1, base="exp.bs", degree=3,
# knots=c(1,5), expected="popmrate", random="clust", verbose=1000)

## Posterior predictions at several time points for an individual
## in cluster 15 with a specific vector of covariates
# Pred_Modbs3_2A <- predict(Mod_bs3_2mix,
# time.pts=seq(0.1,10,by=0.1), data.val=data.frame(agecr=0.2, IsexH=1),
# cluster=15, dataset=simdatn1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{predict}
