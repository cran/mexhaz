\name{ranef.mexhaz}
\alias{ranef.mexhaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Method for extracting random effects from a mexhaz object}
\description{
  Display a data frame containing the cluster-specific random effects
  with their standard errors.
}

\usage{
\method{ranef}{mexhaz}(x, ...)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{mexhaz}.}
  \item{\dots}{not used.}
}

\value{
A data frame with three columns containing the cluster names, the random
effect estimates, and their standard errors.
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

ranef(Mod_weib_mix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{models}
