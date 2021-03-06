\name{summary.mexhaz}
\alias{summary.mexhaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary method for a mexhaz object}
\description{
 Produces a summary of a mexhaz object.}

\usage{
\method{summary}{mexhaz}(object, ...)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{mexhaz}.}
  \item{\dots}{
represents additional parameters directly passed to \code{summary}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mexhaz}}, \code{\link{print.mexhaz}}, \code{\link{print.summary.mexhaz}}
}

\examples{

data(simdatn1)

## Fit of a mixed-effect excess hazard model, with the baseline hazard
## described by a Weibull distribution (without covariables)

Mod_weib <- mexhaz(formula=Surv(time=timesurv,
event=vstat)~1, data=simdatn1, base="weibull", verbose=0)

summary(Mod_weib)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{models}
