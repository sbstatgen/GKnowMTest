\name{getDF-methods}
\docType{methods}
\alias{getDF-methods}
\alias{getDF,annotatedSNPset-method}
\alias{getDF}

\title{\code{getDF}}
\description{
 To extract \code{snp.df} slot of \code{\link{annotatedSNPset-class}}object.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "annotatedSNPset")}}{
getDF(x). x must be an object of \code{\link{annotatedSNPset-class}}
}
}}

\seealso{
\code{\link{anno.create}},\code{\link{annotatedSNPset}}
}

\examples{

  obfile <- system.file("sampleData", "obData.rda", package="GKnowMTest")
  load(obfile) # loads a mapping object ob1 
  res=getDF(ob1)
  head(res)
  
}
\keyword{methods}

