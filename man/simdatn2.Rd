% File mexhaz/man/simdatn1.Rd
% Part of the mexhaz package for R

\name{simdatn2}
\alias{simdatn2}
\title{Simulated dataset}
\description{
  The \code{simdatn2} dataset has 4000 rows and 8 columns. The
  \code{IsexH} variable is simulated with a non-proportional effect.
}
\format{
    This dataset contains the following columns:
    \describe{
      
        \item{age}{Age at diagnosis (continuous).
        }
        \item{agecr}{Centred and rescaled age variable (age-70)/100.
        }
        \item{depindex}{Deprivation index (continuous).
        }
        \item{IsexH}{Sex (0 = Female, 1 = Male).
        }
        \item{clust}{ID number of the cluster.
          }
	  \item{vstat}{Vital status (0 = Alive, 1 = Dead).}
	  \item{timesurv}{Follow-up time (years).}
	  \item{popmrate}{Population (expected) mortality rate at the
        time of censoring.}
    }
}
%\details{}
%\source{
%}
%\examples{}
%\keyword{datasets}
