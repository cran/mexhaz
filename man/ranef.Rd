\name{ranef}
\alias{ranef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Method for extracting random effects}
\description{
 This is a generic function; currently, only the \code{mexhaz} class has
 a method implemented for this function.}

\usage{
ranef(x, ...)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{currently, an object of class \code{mexhaz}.}
  \item{\dots}{not used.}
}

\value{
see the documentation for \code{ranef.mexhaz}.}

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
