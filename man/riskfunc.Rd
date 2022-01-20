\name{riskfunc}
\alias{riskfunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computation of hazard ratio and risk ratio estimates based on a mexhaz
  model}

\description{Function for computing hazard ratio and risk ratio (ratio
  of cumulative probabilities of failure) estimates from a model fitted
  with the \code{mexhaz} function. Corresponding confidence intervals
  are based on the Delta Method or Monte Carlo simulation (based on the
  assumption of multivariate normality of the model parameter
  estimates). This function allows the computation of estimates at one
  time point for several vectors of covariates or for one vector of
  covariates at several time points. When the model includes a random
  effect, two types of predictions can be made: (i) marginal predictions
  (obtained by integration over the random effect distribution) or (ii)
  conditional predictions either for a particular cluster (using the
  corresponding shrinkage estimate) or for the value 0 of the random
  effect.}

\usage{
riskfunc(object, time.pts, data, data.0, clust.name = NULL, 
marginal = TRUE, type = c("hr", "rr"), conf.int = c("delta",
"simul"), level = 0.95, nb.sim = 10000, seed = NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{mexhaz}, corresponding to a
    hazard-based regression model fitted with the \code{mexhaz} function.}
  \item{time.pts}{a vector of numerical values representing the time
    points at which predictions are requested. Time values greater than
    the maximum follow-up time on which the model estimation was based
    are discarded.}
  \item{data}{a \code{data.frame} containing the values of the
    covariates of the population for which hazard ratios or risk ratios
    are to be calculated.}
  \item{data.0}{a \code{data.frame} containing the values of the
    covariates of the reference population. Each row of \code{data.0} is
    used as the reference for the corresponding row of \code{data}.}
  \item{clust.name}{name of the variable in \code{data} (and
    \code{data.0}) containing information on cluster membership to be
    used to obtain conditional predictions from a mixed effect
    model. Standard errors of prediction are obtained with the
    approximation of the conditional mean squared error of prediction
    (without the second-order bias correction term) described in Booth
    and Hobert. This argument is not used if the model is a fixed effect
    model. The default value is \code{NULL}: this corresponds to
    marginal predictions (if \code{marginal} is set to \code{TRUE}, the
    preferred option), or to conditional predictions for the value 0 of
    the random effect (if \code{marginal} is set to \code{FALSE}).}
  \item{marginal}{logical value controlling the type of predictions
    returned by the function when the model includes a random
    intercept. When \code{TRUE}, marginal predictions are computed. The
    marginal survival is obtained by integrating the predicted survival
    over the distribution of the random effect. When \code{FALSE}
    (default value), conditional predictions depending on the value of
    the \code{cluster} argument are calculated.}
  \item{type}{argument specifying the type of predictions to be
  calculated. Selection can be made between \code{"hr"} (hazard ratio)
  and \code{"rr"} (risk ratio, i.e., ratio of cumulative failure
  probabilities).}
  \item{conf.int}{method to be used to compute confidence
    limits. Selection can be made between the following options:
    \code{"delta"} for the Delta Method (default value);
    \code{"simul"} for Monte Carlo simulations (can be time-consuming,
    especially for models using B-splines for the logarithm of the
    baseline hazard).}
  \item{level}{a number in (0,1) specifying the level of confidence for
    computing the confidence intervals of the hazard and the
    survival. By default, this argument is set to \code{0.95}.}
 \item{nb.sim}{integer value representing the number of simulations
   used to estimate the confidence limits for the (excess) hazard
   and the (net) survival. This argument is used only if
   \code{conf.int} is set to \code{"simul"}.}
 \item{seed}{argument allowing the user to set a random seed for
 simulations (only relevant when \code{conf.int} is set to
 \code{"simul"}). The default value is set to \code{NULL} in which case
 a random seed is automatically generated inside the function.}

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
    An object of class \code{resMexhaz} that can be used by the function
  \code{plot.resMexhaz} to produce graphics
  of the hazard ratio or risk ratio curve. It contains the following
  elements:
 \item{results}{a \code{data.frame} consisting of: the time points at
  which values have been calculated; the
  hazard ratio / risk ratio values with their confidence limits.}
\item{type}{type of results returned by the function. The value is
  used by \code{plot.resMexhaz} and \code{lines.resMexhaz}, and can take
  the values \code{"hr"} (hazard ratio) or \code{"rr"} (risk ratio).}
\item{multiobs}{value used by
  \code{plot.resMexhaz} and \code{lines.resMexhaz}, and set to
  \code{FALSE} when estimates are computed at several time
  points for one vector of covariates.}
\item{ci.method}{method used to compute confidence limits.}
\item{level}{level of confidence used to compute confidence limits.}
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
\code{\link{plot.resMexhaz}}, \code{\link{lines.resMexhaz}}
}

\examples{

data(simdatn1)

## Fit of a fixed-effect hazard model, with the baseline hazard
## described by a linear B-spline with two knots at 1 and 5 year and with
## effects of age (agecr), deprivation index (depindex) and sex (IsexH)

Mod_bs1_2 <- mexhaz(formula=Surv(time=timesurv,
event=vstat)~agecr+depindex+IsexH, data=simdatn1, base="exp.bs",
degree=1, knots=c(1,5), verbose=0)

## Risk ratio along time for agecr=0.2 compared to agecr=0.1

RR_Modbs1_2 <- riskfunc(Mod_bs1_2, time.pts=seq(0,10,le=101),
data=data.frame(agecr=0.2,depindex=0,IsexH=1),
data.0=data.frame(agecr=0.1,depindex=0,IsexH=1),type="rr",
conf.int="delta")

plot(RR_Modbs1_2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{predict}
