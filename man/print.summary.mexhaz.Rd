\name{print.summary.mexhaz}
\alias{print.summary.mexhaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Printing method for a summary.mexhaz object}
\description{
 Display the model call, the values of the estimated model parameters,
 as well as the corresponding hazard ratios (only for proportional effects).}

\usage{
\method{print}{summary.mexhaz}(x, ...)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{summary.mexhaz}.}
  \item{\dots}{
represents additional parameters directly passed to \code{print}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mexhaz}}, \code{\link{print.mexhaz}}, \code{\link{summary.mexhaz}}
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
