\name{anno.merge}
\alias{anno.merge}

\title{ Merges multiple annotatedSNPset objects
}

\description{
Function to merge multiple annotatedSNPset objects into a single object.
}

\usage{
anno.merge(anno.list)
}

\arguments{
  \item{anno.list}{
 A list of multiple \code{\link{annotatedSNPset}} objects. 
}

}

\details{
If several objects merged as a single object then all equivalence classes of input objects are merged to form a new equivalence class. \code{eq.map} matrix stores this mapping information i.e which input equivalence classes form the new equivalence classes. For a merged object,\code{snp.eq} is created by using this new equivalence class. The \code{eq.mat} lists are concatenated to hold one binary matrix from each merged object.
}

\value{
It returns an object of annotatedSNPset-class.See \code{\link{annotatedSNPset-class}} for details.
}

\seealso{
\code{\link{anno.create}}
}

\examples{

\dontrun{
  snpfile <- system.file("sampleData", "snpData.rda", package="GKnowMTest")
  pathfile <- system.file("sampleData", "pathData.rda", package="GKnowMTest")
  anmfile <- system.file("sampleData", "anmData.rda", package="GKnowMTest")

  load(snpfile) ## loads snpdf
  load(pathfile) ## loads a R list of gene symbols
  load(anmfile) ##loads annotation matrix

  snp<-rownames(snpdf)
  ob1=anno.create(snp,path.def=pathlist) ## Genelist as Pathlist

  ##creation of annotation matrix
  res=create.anno.mat(snp,base.path="/home/Datasets/Encode/",
      fl.suffix=".annot.wdist.wcoding")
  anm=res[[1]]

  ob2<-anno.create(snp,anno.mat=anm) #annotation matrix as pathlist
 }

obfile <- system.file("sampleData", "obData.rda", package="GKnowMTest")
load(obfile) ##loads two mapping objects ob1 and ob2.
ob3=anno.merge(list(ob1,ob2)) 
ob3

}

