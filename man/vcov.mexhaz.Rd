\name{vcov.mexhaz}
\alias{vcov.mexhaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Method for extracting the covariance matrix of the fixed effects}
\description{
  Display a the covariance matrix of the fixed effects of a mexhaz model.
}

\usage{
\method{vcov}{mexhaz}(object, ...)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{mexhaz}.}
  \item{\dots}{not used.}
}

\value{
The estimated covariance matrix of the fixed effects.
  }

\seealso{
 \code{\link{mexhaz}}}


\examples{

data(simdatn1)

## Fit of a mixed-effect excess hazard model, with the baseline hazard
## described by a Weibull distribution (without covariables)

Mod_weib_mix <- mexhaz(formula=Surv(time=timesurv,
event=vstat)~1, data=simdatn1, base="weibull",
expected="popmrate", verbose=0, random="clust")

vcov(Mod_weib_mix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{models}
