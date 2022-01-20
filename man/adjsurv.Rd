\name{adjsurv}
\alias{adjsurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computation of direct adjusted survival estimates based on a mexhaz
  model} 

\description{Function for computing direct adjusted survival estimates
from a model fitted with the \code{mexhaz}. It can be used to obtain
direct adjusted survival estimates for one or two populations. In the
latter case, survival difference estimates are also
computed. Corresponding variance estimates are based on the Delta Method
(based on the assumption of multivariate normality of the model
parameter estimates). When the model includes a random effect, two types
of predictions can be made: (i) marginal predictions (obtained by
integration over the random effect distribution) or (ii) conditional
predictions either for a particular cluster (using the corresponding
shrinkage estimate) or for the value 0 of the random effect.}

\usage{
adjsurv(object, time.pts, data, data.0 = NULL, weights = NULL,
clust.name = NULL, marginal = TRUE, level = 0.95)
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
    covariates of the population for which direct adjusted estimates
    are to be calculated.}
  \item{data.0}{an optional \code{data.frame} containing the values of
    the covariates of a second population for which direct adjusted
    estimates can also be calculated (and compared with those of the
    first population). The default value is set to \code{NULL}.}
  \item{clust.name}{name of the variable in \code{data} (and
    \code{data.0}) containing information on cluster membership to be
    used to obtain conditional predictions from a mixed effect
    model. Standard errors of prediction are obtained with the
    approximation of the conditional mean squared error of prediction
    (without the second-order bias correction term) described in Booth
    and Hobert. This argument is not used if the model is a fixed effect
    model. The default value is \code{NULL}: this corresponds to
    marginal predictions (if \code{marginal} is set to \code{TRUE}, the preferred
    option), or to conditional predictions for the value 0 of the random
    effect (if \code{marginal} is set to \code{FALSE}).}
  \item{weights}{optional argument specifying the weights to be
    associated with each row of \code{data} (and \code{data.0}). the
    default value is set to \code{NULL} which corresponds to attributing
    to each row of the dataset(s) a weight equal to one over the total
    number of rows.}
  \item{marginal}{logical value controlling the type of predictions
    returned by the function when the model includes a random
    intercept. When \code{TRUE}, marginal predictions are computed. The
    marginal survival is obtained by integrating the predicted survival
    over the distribution of the random effect. When \code{FALSE}
    (default value), conditional predictions depending on the value of
    the \code{cluster} argument are calculated.}
  \item{level}{a number in (0,1) specifying the level of confidence for
    computing the confidence intervals of the hazard and the
    survival. By default, \code{level=0.95}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  An object of class \code{resMexhaz} that can be used by the function
  \code{plot.resMexhaz} to produce graphics
  of the direct adjusted survival curve. It contains the following
  elements:
 \item{results}{a \code{data.frame} consisting of: the time points at
  which the direct adjusted survival values have been calculated; the
  direct ajusted survival values with their confidence limits for
  population \code{data}; the direct ajusted survival values with their
  confidence limits for population \code{data.0}; the direct adjusted
  survival difference estimates with their confidence limits.}
\item{type}{type of results returned by the function. The value is
  used by \code{plot.resMexhaz} and \code{lines.resMexhaz}, and set to
  \code{"as"} (adjusted survival).}
\item{multiobs}{value used by
  \code{plot.resMexhaz} and \code{lines.resMexhaz}, and set to
  \code{FALSE} (computation of the adjusted survival at several time
  points for one vector of covariates).}
\item{ci.method}{method used to compute confidence limits. Currently set
  to \code{"delta"} as only the Delta Method is implemented.}
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

## Direct adjusted survival for the simdatn1 population
DAS_Modbs1_2 <- adjsurv(Mod_bs1_2, time.pts=seq(1,10),
data=simdatn1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{predict}
