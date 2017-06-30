\name{pop2011}
\alias{pop2011}
\alias{popDE}
\docType{data}
\title{Berlin and German Population by Age Group, 2011}
\description{
  Population numbers from Berlin are available in the
  city district x age group (5-year intervals) matrix \code{pop2011}.
  The corresponding age distribution for whole Germany is stored in the
  vector \code{popDE}.
}
\usage{
## Berlin population by city district and age group, 2011
pop2011

## German population by age group, 2011
popDE
}
\format{
  \describe{
    \item{pop2011:}{a named, integer-valued
      12 (city districts) x 15 (age groups) matrix.}
    \item{popDE:}{a named integer vector of length 15 (age groups).}
  }
}
\source{
  \describe{
    \item{pop2011:}{numbers extracted from
      \url{https://www.statistik-berlin-brandenburg.de/webapi/opendatabase?id=BevBBBE}
      as of 2011-12-31 (before census), published by
      \emph{Amt fuer Statistik Berlin-Brandenburg}
      (Statistical Office of Berlin-Brandenburg)
      under the \sQuote{CC BY 3.0 DE} license
      (\url{https://creativecommons.org/licenses/by/3.0/de/}).}
    \item{popDE:}{numbers extracted from
      \url{https://www-genesis.destatis.de/genesis/online/link/tabellen/12411-0005}
      as of 2010-12-31, published by
      \emph{Statistisches Bundesamt}
      (Destatis, Federal Statistical Office of Germany)
      under the \sQuote{Data licence Germany - attribution - Version 2.0}
      (\url{https://www.govdata.de/dl-de/by-2-0}).}
  }
}
\author{Sebastian Meyer}
\keyword{datasets}
