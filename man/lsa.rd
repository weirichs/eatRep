\name{lsa}

\docType{data}

\alias{lsa}

\title{Achievement data from two large-scale assessments of 2010 and 2015.}

\description{
This example data set contains fictional achievement scores of 11637 students from three countries
and two times of measurement in two domains (reading and listening comprehension) in the long format.
The data set contains nested multiple imputed plausible values of achievement scores as well as some
demographic variables. Illustrating trend analyses, data from two fictional time points (2010 and 2015)
are included.

The data set can be used for several illustration purposes. For example, if only multiple imputation
should be considered (without nesting), simply use only cases from the first nest (by subsetting). If
only one time of measurement should be considered (i.e., without any trend analyses), simply choose
only cases from 2010 or 2015. If only reading or listening should be considered, choose the desired
domain by subsetting according to the \code{domain} column.
}

\usage{data(lsa)}

\format{'data.frame':   77322 obs. of  25 variables
  \describe{
    \item{year}{Year of evaluation}
    \item{idstud}{individual student identification}
    \item{idclass}{class identifier}
    \item{wgt}{Total case weight}
    \item{L2wgt}{School weight (level 2 weight)}
    \item{L1wgt}{Student weight (level 1 weight)}
    \item{jkzone}{jackknifing zone (jk2) }
    \item{jkrep}{jackknife replicate}
    \item{imp}{Number of imputation}
    \item{nest}{Number of nest (for nested imputation only)}
    \item{country}{The country an examinee stems from}
    \item{sex}{student's sex}
    \item{ses}{student's socio-economical status}
    \item{mig}{student's migration background}
  	\item{domain}{The domain the corresponding score belongs to}
  	\item{score}{student's achievement score (corresponding to the domain reading or listening, and to the imputation 1, 2, or 3)}
  	\item{comp}{student's competence level}
  	\item{failMin}{dichotomous indicator whether the student fails to fulfill the minimal standard}
  	\item{passReg}{dichotomous indicator whether the student fulfills at least the regular standard}
  	\item{passOpt}{dichotomous indicator whether the student fulfills the optimal standard}
  	\item{leSore}{linking error of each student's achievement score}
  	\item{leComp}{linking error of each student's competence level}
  	\item{leFailMin}{linking error of each student's indicator of failing to fulfill the minimal standard}
  	\item{lePassReg}{linking error of each student's indicator of fulfilling the regular standard}
  	\item{lePassOpt}{linking error of each student's indicator of fulfilling the optimal standard}
 }
}

\source{Simulated data}

%\references{
%}

\keyword{datasets}


