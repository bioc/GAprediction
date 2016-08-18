\name{UL.mod.cv}
\alias{UL.mod.cv}
\docType{data}
\title{
A glmnet-object trained to perform gestational age prediction.
}
\description{
The glmnet-object consists of a Lasso-regression model 'trained' to perform gestational age predictions. It is called by the wrapper function \code{\link{predictGA}}, which is more user-friendly.
}
\details{
The trained Lasso-model contains cross-validated estimates of the penalty term lambda that regulates the number of CpG sites needed for gestational age prediction. It is called by the glmnet-inherited predict function with a matrix of CpG betas (with values between 0 and 1) that conforms to the Illumina HumanMethylation450 platform. The gestational age estimates used to train the regression model were taken from the MoBa cohort and are based on ultrasound.
}
\source{
Magnus P, Irgens LM, Haug K, Nystad W, Skjaerven R, Stoltenberg C, MoBa Study Group. Cohort profile: the Norwegian mother and child cohort study (MoBa). International journal of epidemiology. 2006 Oct 1;35(5):1146-50.
}
\references{
Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1-22. URL http://www.jstatsoft.org/v33/i01/.
}
\examples{
## Extract all non-zero regression coefficients
temp <- as.matrix( coef( UL.mod.cv ) )
allNonZeroCoefs <- rownames( temp )[ temp[,1]!=0 ]
allNonZeroCoefs[ -1 ]
}
\keyword{datasets}
