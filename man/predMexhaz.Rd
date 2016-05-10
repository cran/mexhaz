\name{predMexhaz}
\alias{predMexhaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predictions based on a mexhaz model}
\description{
  Function for predicting the (excess) hazard and the corresponding
  (net) survival from
  a model fitted with the \code{mexhaz} function for a particular vector
  of covariates. If the survival model
  was fitted with an expected hazard, the estimates obtained are excess
  hazard and net survival estimates. When the model includes a random
  effect, the predicted values are obtained for the value 0 of the
  random effect. Confidence limits
  can be obtained by Monte-Carlo simulation (for all types of baseline
  hazard) and by the Delta Method (not available for Weibull hazard).
  This function allows the computation of the hazard and the survival at
  one time point for several vectors of covariables or for one vector of
  covariables at
  several time points.
}

\usage{predMexhaz(model, time.pts, data.val = data.frame(.NotUsed=NA),
conf.int=c("none","delta","simul"), nb.sim = 10000)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{an object of class \code{mexhaz}, corresponding to a
    survival model fitted with the \code{mexhaz} function.}
  \item{time.pts}{a vector of numerical values representing the time
    points at which predictions are requested. Time values greater than
    the maximum follow-up time on which the model estimation was based
    are discarded.}
  \item{data.val}{a \code{data.frame} containing the values of the
    covariables at which predictions should be calculated.}
  \item{conf.int}{method to be used to compute confidence
    limits. Selection can be made between the following options:

    \code{"none"} indicates absence of confidence limits estimation
    (default value);
    
    \code{"delta"} for the Delta Method (this option is not available for
    models using a Weibull hazard);
    
    \code{"simul"} for Monte Carlo simulations (can be time-consuming, especially for models using B-splines for the logarithm of the baseline hazard).}
  \item{nb.sim}{integer value representing the number of simulations
  used to estimate the 95\% confidence limits for the (excess) hazard and the (net) survival. This argument is used only if \code{conf.int="simul"}.}
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
  vector of covariables. These variances can be used to estimate the population
  (net) survival. The object \code{variances} is produced only when
  B-splines are used to model the logarithm of the
  baseline (excess) hazard with the option \code{conf.int="delta"}.}
\item{type}{the type of predictions produced. Can take the value
    \code{"multitime"} (computation of the hazard and the survival at
  at several time points for one vector of covariables) or
    \code{multiobs} (computation of the hazard and the survival at
  at one time point for several vectors of covariables). This value is
    used by \code{plot.predMexhaz} and \code{points.predMexhaz}.}
\item{ci.method}{the method used to compute confidence limits.}
\item{nb.sim}{number of simulations used to estimate the 95\% confidence
  limits (set to \code{NA} if confidence limits were not requested or
  calculated with the Delta Method).}
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
Pred_Modbs1_2A <- predMexhaz(Mod_bs1_2, time.pts=seq(0.1,10,by=0.1),
data.val=data.frame(agecr=0,depindex=0.5,IsexH=1), conf.int="delta")

## Prediction for several vectors of covariates at one time point
Pred_Modbs1_2B <- predMexhaz(Mod_bs1_2, time.pts=10,
data.val=data.frame(agecr=c(-0.2,-0.1,0), depindex=c(0.5,0.5,0.5),
IsexH=c(1,1,1)), conf.int="delta")

## Prediction for all individuals of the study population at one time point
Pred_Modbs1_2C <- predMexhaz(Mod_bs1_2, time.pts=10,
data.val=simdatn1, conf.int="delta")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{predict}
