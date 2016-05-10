\name{print.predMexhaz}
\alias{print.predMexhaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print method for a predMexhaz object}
\description{
Display the first lines of the \code{data.frame} containing the
predictions provided by the \code{predMexhaz} function.
}

\usage{
\method{print}{predMexhaz}(x, ...)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{predMexhaz}.}
  \item{\dots}{
represents additional parameters directly passed to \code{print}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{predMexhaz}}
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

print(Pred_Modbs1_2A)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
