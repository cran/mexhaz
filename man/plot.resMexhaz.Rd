\name{plot.resMexhaz}
\alias{plot.resMexhaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot method for a resMexhaz object}
\description{Function for plotting various quantities calculated from a
  \code{mexhaz} model.}

\usage{
\method{plot}{resMexhaz}(x, conf.int = TRUE, lty.pe = "solid", lty.ci =
"blank", col.ci = "blue", alpha.col.ci = 0.25, ylim = NULL, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{resMexhaz}, corresponding to
    various predictions based on a survival
    model fitted with the \code{mexhaz} function. Predictions can be
    obtained for multiple times for one vector of covariables
    (\code{"multitime"}) or for several vectors of covariables at one
    time point (\code{"multiobs"}). The \code{plot()} function only
    applies to the \code{"multitime"} type of predictions.
  }
  \item{conf.int}{logical value allowing the user to decide whether to
    plot the confidence limits of the survival (or hazard).}
  \item{lty.pe}{type of line used for drawing the hazard/survival estimate.}
  \item{lty.ci}{type of line used for drawing the confidence
    limits.}
  \item{col.ci}{color used to fill in the polygon defined by the
    confidence limits.}
  \item{alpha.col.ci}{parameter used internally by the \code{rgb()}
    function to set the color transparency.}
  \item{ylim}{when set to \code{NULL} (default value), the range of the
    y-axis is based on the range of the confidence limits of the
    quantity to be plotted.}
  \item{...}{additional parameters that are directly passed to the
    \code{plot} function. These parameters will apply simultaneously to
    the point estimate and confidence limit curves.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{adjsurv}}, \code{\link{riskfunc}},
\code{\link{lines.resMexhaz}} }

\examples{

data(simdatn1)

## Fit of a fixed-effect hazard model, with the baseline hazard
## described by a linear B-spline with two knots at 1 and 5 year and with
## effects of age (agecr), deprivation index (depindex) and sex (IsexH)

Mod_bs2 <- mexhaz(formula=Surv(time=timesurv,
event=vstat)~agecr+depindex+IsexH, data=simdatn1, base="exp.bs",
degree=2, knots=c(1,5), verbose=0)

## Relative risk (ratio of cumulative incidence curves) for men versus
## women

RR_Modbs2 <- riskfunc(Mod_bs2, time.pts=seq(0,10, by=0.1),
data=data.frame(agecr=0.1, IsexH=1, depindex=0),
data.0=data.frame(agecr=0.1, IsexH=0, depindex=0), conf.int="delta",
type="rr")

plot(RR_Modbs2)
abline(h=1, lty="dashed")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
