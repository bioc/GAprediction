\name{predictGA}
\alias{predictGA}

\title{
Predict gestational age in days from conception
}
\description{
The function predictGA takes a matrix with Illumina HumanMethylation450 type DNA methylation data. Column names must designate CpG sites (i.e. 'cgXXXXXX', X=number) and row names samples IDs.
}
\usage{
predictGA(mldat, transp=TRUE, se=TRUE)
}

\arguments{
  \item{mldat}{
  A matrix containing DNA methylation beta values (0<=beta<=1)
  }
  \item{transp}{
  If TRUE (default), the transpose is automatically taken if the number of rows is greater than the number of columns.
  }
  \item{se}{
  If se=TRUE, the estimated coefficients are based on the prediction model with the lambda penalty term being allowed to vary up to one standard error within the minimum. If se=FALSE, the minimum lambda is assumed.
  }
}

\details{
The minimum lambda (se=FALSE) may result in slightly better predictions, however substantially more CpG sites are needed for estimation. Since the prediction difference is hardly noticeable se=TRUE is the default option.
}
\value{
The function returns estimated gestational age predictions, together with samples IDs as row names, in a data.frame object.
}
\references{

Jon Bohlin, Siri E. Haaberg, Per Magnus, et al. (2016). Prediction of gestational age based on genome-wide differentially methylated regions. Genome Biology (in review)

Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
  Regularization Paths for Generalized Linear Models via Coordinate
  Descent. Journal of Statistical Software, 33(1), 1-22. URL
  http://www.jstatsoft.org/v33/i01/.
}
\author{
Jon Bohlin
}
\note{
Requires quite a bit of memory due to the large DNA methylation matrix required for the prediction model.
}

\examples{
## Make a mock Illumina HumanMethylation450 type DNA methylation matrix
cpgs <- extractSites( type="se" )
allcpgs <- extractSites( type="all" )
numsamples <- 100
mlmatr <- matrix( NA, ncol=length( allcpgs ), nrow=numsamples )
mlmatr <- data.frame( mlmatr )
for( i in cpgs )
  mlmatr[,i] <- runif( numsamples, min=0, max=1 )
## Perform gestational age prediction
mypred <- predictGA( mlmatr )

}
