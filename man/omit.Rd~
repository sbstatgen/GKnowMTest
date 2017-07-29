\name{omit}
\alias{omit}

\title{
Omit SNPs from an Object of \code{\link{annotatedSNPset-class}}
}
\description{
This function removes input SNPs from an mapping object of \code{\link{annotatedSNPset-class}} and returns a new object.
}

\usage{
omit(ob,val)
}

\arguments{
  \item{ob}{Mapping object from which input SNPs to be removed}
  \item{val}{A vector of either SNPs(rsIDs) or indices of SNPs}
}

\value{
 It returns an object of annotatedSNPset-class.See \code{\link{annotatedSNPset-class}} for details.
}

\seealso{
\code{\link{anno.create}}
}

\examples{

obfile <- system.file("sampleData", "obData.rda", package="GKnowMTest")
load(obfile) # loads a mapping object ob1 

ob3<-omit(ob1,5:10)
ob3

}

