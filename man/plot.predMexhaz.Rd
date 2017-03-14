\name{plot.predMexhaz}
\alias{plot.predMexhaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot method for a predMexhaz object}
\description{Function for plotting the predicted (excess) hazard or
  (net) survival based on a \code{predMexhaz} object.}

\usage{
\method{plot}{predMexhaz}(x, which=c("surv", "hazard"), conf.int=TRUE,
lty.pe="solid", lty.ci="dashed", ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{predMexhaz}, corresponding to
    (excess) hazard and (net) survival predictions based on a  survival
    model fitted with the \code{mexhaz} function. Predictions can be
    obtained for multiple times for one vector of covariables
    (\code{"multitime"}) or for several vectors of covariables at one
    time point (\code{"multiobs"}). The \code{plot()} function only
    applies to the \code{"multitime"} type of predictions.
  }
  \item{which}{type of curve to be plotted. Selection can be made
    between \code{"surv"} (default value) for the (net) survival curve
    and \code{"hazard"} for the (excess) hazard.}
  \item{conf.int}{logical value allowing the user to decide whether to
    plot the confidence limits of the survival (or hazard).}
  \item{lty.pe}{type of line used for drawing the hazard/survival estimate.}
  \item{lty.ci}{type of line used for drawing the confidence
    limits.}
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

\seealso{
\code{\link{predict.mexhaz}}, \code{\link{points.predMexhaz}}
}

\examples{

data(simdatn1)

## Fit of a fixed-effect hazard model, with the baseline hazard
## described by a linear B-spline with two knots at 1 and 5 year and with
## effects of age (agecr), deprivation index (depindex) and sex (IsexH)

Mod_bs2_2 <- mexhaz(formula=Surv(time=timesurv,
event=vstat)~agecr+depindex+IsexH, data=simdatn1, base="exp.bs",
degree=2, knots=c(1,5), verbose=0)

## Prediction at several time points for one vector of covariates
Pred_Modbs2_2A <- predict.mexhaz(Mod_bs2_2, time.pts=seq(0,10,by=0.1),
data.val=data.frame(agecr=0,depindex=0.5,IsexH=1))

plot(Pred_Modbs2_2A, which="hazard")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
