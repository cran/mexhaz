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
parameter estimates). When the model includes a random effect, three types of
predictions can be made: (i) marginal predictions (obtained by
integration over the random effect distribution), (ii) cluster-specific
posterior predictions for an existing cluster, or (iii) conditional
predictions for a given quantile of the random effect distribution (by
default, for the median value, that is, 0).}

\usage{
adjsurv(object, time.pts, data, data.0 = NULL, weights = NULL,
marginal = TRUE, quant.rdm = 0.5, cluster = NULL, quant.rdm.0 = 0.5,
cluster.0 = NULL, level = 0.95, dataset = NULL)
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
  \item{weights}{optional argument specifying the weights to be
    associated with each row of \code{data} (and \code{data.0}). the
    default value is set to \code{NULL} which corresponds to attributing
    to each row of the dataset(s) a weight equal to one over the total
    number of rows.}
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
  \item{quant.rdm.0}{random effect distribution quantile value to be used with \code{data.0} (see
    argument \code{quant.rdm} for details).}
  \item{cluster.0}{cluster value to be used with \code{data.0} (see
    argument \code{cluster} for details).}
  \item{level}{a number in (0,1) specifying the level of confidence for
    computing the confidence intervals of the hazard and the
    survival. By default, \code{level=0.95}.}
    \item{dataset}{original dataset used to fit the \code{mexhaz} object
    given as argument to the function. This argument is only necessary
    if cluster-specific posterior predictions are requested (and if the
    dataset is not already provided in the \code{mexhaz} object). The
    default value is set to \code{NULL}.}

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
