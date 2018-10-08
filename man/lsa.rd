\name{lsa}

\docType{data}

\alias{lsa}

\title{Achievement data from large-scale assessment.}

\description{
This data set contains fictional achievement scores of 7518 students in two domains (reading and
listening comprehension) in the long format. The data set contains multiple imputed plausible
values of achievement scores as well as some demographic variables. Illustrating trend analyses,
data from two fictional time points (2010 and 2015) are included.
}

\usage{data(lsa)}

\format{'data.frame':   90216 obs. of  22 variables
  \describe{
    \item{year}{Year of evaluation}
    \item{idstud}{individual student identification}
    \item{wgt}{Individual student case weight}
    \item{jkzone}{jackknifing zone (jk2) }
    \item{jkrep}{jackknife replicate}
    \item{imp}{Number of imputation}
    \item{nest}{Number of nest (for nested imputation only)}
    \item{country}{The country an examinee stems from}
    \item{sex}{student's sex}
    \item{ses}{student's socio-economical status}
    \item{mig}{student's migration background}
  	\item{domain}{The domain the correponding score belongs to}
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


