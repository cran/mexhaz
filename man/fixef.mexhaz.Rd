\name{fixef.mexhaz}
\alias{fixef.mexhaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Method for extracting fixed effects from a mexhaz object}
\description{
  Display a vector containing the fixed effects of a mexhaz model.
}

\usage{
\method{fixef}{mexhaz}(x, ...)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{mexhaz}.}
  \item{\dots}{not used.}
}

\value{
A vector containing the fixed effect estimates.
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

fixef(Mod_weib_mix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{models}
