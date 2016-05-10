\name{points.predMexhaz}
\alias{points.predMexhaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Points method for a predMexhaz object}
\description{Function for adding to an already existing graphical window
  the predicted (excess) hazard or (net) survival based on a \code{predMexhaz} object.}

\usage{
\method{points}{predMexhaz}(x,which=c("surv","hazard"),conf.int=TRUE, ...)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{predMexhaz}, corresponding to
    (excess) hazard and (net) survival predictions based on a  survival
    model fitted with the \code{mexhaz} function. Predictions can be
    obtained for multiple times for one vector of covariates
    (\code{"multitime"}) or for several vectors of covariates at one
    time point (\code{"multiobs"}). The \code{points()} function only
    applies to the \code{"multitime"} type of predictions.
  }
  \item{which}{the type of curve to be plotted. Selection can be made
    between \code{"surv"} (default value) for the (net) survival curve
    and \code{"hazard"} for the (excess) hazard.}
  \item{conf.int}{logical values allowing the user to decide whether to
    plot the 95\% confidence limits of the survival (or hazard), if
    confidence limits are available in the \code{predMexhaz} object.}
  \item{...}{additional parameters that are directly passed to the \code{points} function.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{predMexhaz}}, \code{\link{plot.predMexhaz}}
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

plot(Pred_Modbs1_2A,which="hazard",col="red")
points(Pred_Modbs1_2A,which="hazard",conf.int=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
