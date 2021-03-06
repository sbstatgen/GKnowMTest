\name{getMAT-methods}
\docType{methods}
\alias{getMAT-methods}
\alias{getMAT,annotatedSNPset-method}
\alias{getMAT}

\title{\code{getMAT}}
\description{
 It returns no.of equivalence class by no.of annotation matrix of an object of \code{\link{annotatedSNPset-class}}.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "annotatedSNPset")}}{
getMAT(x). x must be an object of \code{\link{annotatedSNPset-class}}
}
}}

\seealso{
\code{\link{anno.create}},\code{\link{annotatedSNPset}}
}

\examples{
  
  obfile <- system.file("sampleData", "obData.rda", package="GKnowMTest")
  load(obfile) # loads a mapping object ob1 
  res=getMAT(ob1)
  res[1:10,1:6]  

}
\keyword{methods}

