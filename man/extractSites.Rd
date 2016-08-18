\name{extractSites}
\alias{extractSites}

\title{
Extract CpG sites for gestational age prediction
}
\description{
The function allows the user to extract CpG sites used for gestational age prediction with the function \code{\link{predictGA}}.
}
\usage{
extractSites(type="se")
}

\arguments{
  The argument type=c("se", "min", "all") can be used to specify which CpGs are to be extracted. "se" designates the CpGs needed by the predictGA function if the penalty term lambda is to be set to one standard error within the minimum, "min" specifies the minimum lambda, while "all" returns the complete sets of CpGs in the \code{\link{UL.mod.cv}} object.
  
  \item{type}{ - a string that can be "se" (default), "min" or "all", depending on which CpGs is wanted by the user.}
  
}
\details{
Use this function if predictGA fails due to missing predictor CpGs, or to see which CpGs are used by \code{\link{predictGA}} for gestational age prediction.
}
\value{
Returns a vector with the requested CpG sites.
}
\author{
Jon Bohlin
}

\seealso{

\code{\link{predictGA}}, \code{\link{UL.mod.cv}}
}
\examples{
CpGs <- extractSites( type="se" )
}
