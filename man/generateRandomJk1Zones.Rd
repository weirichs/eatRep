\name{generateRandomJk1Zones}
\alias{generateRandomJk1Zones}

\title{Generates random jackknife-1 zones based on sampling units in the data set.}

\description{
Function adds randomly generated jackknife-1 zones to the data.
}

\usage{
generateRandomJk1Zones (datL, unit, nZones, name = "randomCluster")
}
\arguments{
  \item{datL}{Data frame containing at least the primary sampling unit variable}
  \item{unit}{Variable name or column number of the primary sampling unit (i.e.
  student or class identifier}
  \item{nZones}{integer: number of jackknife zones. Note: The umber of jackknife
  zones must not exceed the number of distinct sampling units}
  \item{name}{New name of the jackknife-zone variable in the data set}
}
\value{
The original data with an additional column of the jackknife-zone variable
}

\author{
Sebastian Weirich
}
\examples{
data(lsa)

### We only consider year 2010
lsa10<- lsa[which(lsa[,"year"] == 2010),]
lsa10<- generateRandomJk1Zones(datL = lsa10, unit="idstud", nZones = 80)
}
