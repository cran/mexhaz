\name{predict.mexhaz}
\alias{predict.mexhaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predictions based on a mexhaz model}
\description{
  Function for predicting the (excess) hazard and the corresponding
  (net) survival from a model fitted with the \code{mexhaz} function for
  a particular vector of covariables. If the survival model was fitted
  with an expected hazard (excess hazard model), the estimates obtained
  are excess hazard and net survival estimates. Corresponding variance
  estimates are based on the Delta Method or Monte Carlo simulation
  (based on the assumption of multivariate normality of the model
  parameter estimates). When the model includes a
  random effect, predictions can be made for a particular cluster (using
  the corresponding shrinkage estimate) or for the value 0 of the random
  effect. This function allows the computation of the hazard and the
  survival at one time point for several vectors of covariables or for
  one vector of covariables at several time points.
}

\usage{\method{predict}{mexhaz}(object, time.pts, data.val = data.frame(.NotUsed=NA),
cluster=NULL, conf.int=c("delta","simul","none"), level=0.95,
delta.type.h=c("log", "plain"), delta.type.s=c("log-log", "log",
"plain"), nb.sim = 10000, include.gradient=FALSE, ...)
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
    covariables at which predictions should be calculated.}
  \item{cluster}{a single number or character string corresponding to
    the name of the cluster for which predictions should be
    calculated. In that case, standard errors of prediction are obtained
    with the approximation of the conditional mean squared error of
    prediction (without the second-order bias correction term) described in Booth and
    Hobert. This argument is not used if the model is a fixed effect
    model. The default value is \code{NULL}: this corresponds to
    prediction for the value 0 of the random effect.}
  \item{conf.int}{method to be used to compute confidence
    limits. Selection can be made between the following options:
    \code{"delta"} for the Delta Method (default value);
    \code{"simul"} for Monte Carlo simulations (can be time-consuming,
    especially for models using B-splines for the logarithm of the
    baseline hazard);
    \code{"none"} indicates absence of confidence limits estimation.
}
  \item{level}{a number in (0,1) specifying the level of confidence for computing the confidence
    intervals of the hazard and the survival. By default,
    \code{level=0.95}.}
    \item{delta.type.h}{type of confidence limits for the hazard when
      using the Delta Method. If
      \code{delta.type.h="log"} (default value), the confidence limits are based on a Wald-type confidence interval for the
      logarithm of the hazard, otherwise they are based directly on a
      Wald-type CI for the hazard.}
    \item{delta.type.s}{type of confidence limits for the survival when
      using the Delta Method. If
      \code{delta.type.s="log-log"} (default value), the confidence limits are based on a Wald-type confidence interval for the
      logarithm of the cumulative hazard; if \code{delta.type.s="log"},
      they are based on a Wald-type CI for the
      cumulative hazard; otherwise they are based directly on a
      Wald-type CI for the survival.}
  \item{include.gradient}{logical value allowing the function to return
  the components of the gradient of the logarithm of the hazard and of
  the logarithm of the cumulative hazard for each prediction. This
  argument is used only if \code{conf.int="delta"}. The default value
  is \code{FALSE}.}
  \item{nb.sim}{integer value representing the number of simulations
  used to estimate the confidence limits for the (excess) hazard
  and the (net) survival. This argument is used only if
  \code{conf.int="simul"}.}
\item{\dots}{for potential additional parameters.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  An object of class \code{predMexhaz} that can be used by
  the functions \code{plot.predMexhaz} and \code{points.predMexhaz} to produce graphics of the (excess) hazard and
  the (net) survival. It contains the following elements:
  \item{call}{the \code{mexhaz} function call on which the model is based.}
\item{results}{a \code{data.frame} consisting of: the time
  points at which the (excess) hazard and the (net) survival have been
  calculated; the values of the covariables used to estimate the
  (excess) hazard and the (net) survival; the
  (excess) hazard values with their confidence limits; and the (net) survival values with their confidence limits.}
\item{variances}{a \code{data.frame} consisting of two columns: the
  variance of the logarithm of the (excess) hazard and the variance of
  the (excess) cumulative hazard for each time points or each
  vector of covariables. The object \code{variances} is produced only when
  \code{conf.int="delta"}.}
\item{grad.loghaz}{a \code{data.frame} consisting of the components of
  the gradient of the logarithm of the (excess) hazard for each time points or each
  vector of covariables. The number of columns corresponds to the number
  of model parameters. This object is produced only when
  \code{conf.int="delta"} and \code{include.gradient=TRUE}.}
\item{grad.logcum}{a \code{data.frame} consisting of the components of
  the gradient of the logarithm of the (excess) cumulative hazard for each time points or each
  vector of covariables. The number of columns corresponds to the number
  of model parameters. This object is produced only when
  \code{conf.int="delta"} and \code{include.gradient=TRUE}.}
\item{vcov}{a matrix corresponding to the covariance matrix used to
  compute the confidence intervals.}
\item{type}{the type of predictions produced. Can take the value
    \code{"multitime"} (computation of the hazard and the survival at
  at several time points for one vector of covariables) or
    \code{multiobs} (computation of the hazard and the survival at
  at one time point for several vectors of covariables). This value is
    used by \code{plot.predMexhaz} and \code{points.predMexhaz}.}
\item{ci.method}{the method used to compute confidence limits.}
\item{level}{level of confidence used to compute confidence limits.}
\item{delta.type}{type of confidence limits for the hazard and the
  survival when using the Delta Method.}
\item{nb.sim}{number of simulations used to estimate the confidence
  limits when \code{ci.method = "simul"}.}
}
\references{Charvat H, Remontet L, Bossard N, Roche L, Dejardin O,
  Rachet B, Launoy G, Belot A; CENSUR Working Survival Group. A
  multilevel excess hazard model to estimate net survival on
  hierarchical data allowing for non-linear and non-proportional effects
  of covariates. Stat Med 2016;35:3066-3084 (doi: 10.1002/sim.6881)

  Booth JG, Hobert JP. Standard errors of prediction in generalized
  linear mixed models. J Am Stat Assoc 1998;93:262-272 (doi: 10.2307/2669622).
}
\author{
Hadrien Charvat, Aurelien Belot
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{print.predMexhaz}}, \code{\link{plot.predMexhaz}}, \code{\link{points.predMexhaz}}
}

\examples{

data(simdatn1)

## Fit of a fixed-effect hazard model, with the baseline hazard
## described by a linear B-spline with two knots at 1 and 5 year and with
## effects of age (agecr), deprivation index (depindex) and sex (IsexH)

Mod_bs1_2 <- mexhaz(formula=Surv(time=timesurv,
event=vstat)~agecr+depindex+IsexH, data=simdatn1, base="exp.bs",
degree=1, knots=c(1,5), verbose=0)

## Prediction at several time points for one vector of covariates
Pred_Modbs1_2A <- predict.mexhaz(Mod_bs1_2, time.pts=seq(0.1,10,by=0.1),
data.val=data.frame(agecr=0,depindex=0.5,IsexH=1))

## Prediction for several vectors of covariates at one time point
Pred_Modbs1_2B <- predict.mexhaz(Mod_bs1_2, time.pts=10,
data.val=data.frame(agecr=c(-0.2,-0.1,0), depindex=c(0.5,0.5,0.5),
IsexH=c(1,1,1)))

## Prediction for all individuals of the study population at one time point
Pred_Modbs1_2C <- predict.mexhaz(Mod_bs1_2, time.pts=10,
data.val=simdatn1)


# Example of cluster-specific prediction (not run)

## Fit of a mixed-effect excess hazard model, with the baseline hazard
## described by a cubic B-spline with two knots at 1 and 5 year

# Mod_bs3_2mix <- mexhaz(formula=Surv(time=timesurv,
# event=vstat)~agecr+IsexH, data=simdatn1, base="exp.bs", degree=3,
# knots=c(1,5), expected="popmrate", random="clust", verbose=1000)

## Prediction at several time points for an individual in cluster 15
## with a specific vector of covariables
# Pred_Modbs3_2A <- predict.mexhaz(Mod_bs3_2mix,
# time.pts=seq(0.1,10,by=0.1), data.val=data.frame(agecr=0.2, IsexH=1),
# cluster=15)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{predict}
