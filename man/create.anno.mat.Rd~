\name{create.anno.mat}
\alias{create.anno.mat}

\title{
Creation of annotation matrix of a given a set of rsIDs.
}

\description{
This method creates a SNP level annotation matrix using 454 ENCODE functional annotations (Pickrell, 2016) that can be then passed to anno.create function to create a annotatedSNPset object.
}

\usage{
create.anno.mat(snp,base.path='annotation-',fl.suffix=".txt",ncols=454,
                start.anno=4,snp.col=3,header=TRUE,index=NULL)
}

\arguments{
  \item{snp}{
A character vector of rsIDs
}
  \item{base.path}{
The path of the files (physical location on the user's computer) where the ENCODE annotations are stored.Annotation file of \eqn{J^{th}}  chromosome must contain 'chrJ' term in its file name.If the location of an Encode file of chromosome 22 is '/home/Encode/anno-chr22.Encode.txt',then \code{base.path} must be specified as '/home/Encode/anno-'.
}  
  \item{fl.suffix}{
Suffix of the ENCODE annotations files.If an Encode file is 'anno-chr22.Encode.txt',
then \code{fl.suffix} must be specified as '.Encode.txt'.
}  
  \item{ncols}{
Total number of columns in ENCODE annotations files. 
} 
  \item{start.anno}{
Column number that indicates starting of annotation column in ENCODE annotations files.
}  
  \item{snp.col}{
Column Number or Column Name that indicates the rsIDs in ENCODE annotations files.
}  
  \item{header}{
Logical value to indicate column names of ENCODE annotation files.Default is TRUE.
}  
  \item{index}{
 An index file or NULL. See details.  
}   
  
}

\details{
To run this function, all the 24 ENCODE annotation files (Pickrell, 2016) for the different chromosomes should be downloaded an uncompressed before running this function or user can use his/her own ENCODE annotation file.If the 'annotaion-CHR.txt' files are not found at the \code{base.path} for any of the chromosomes required, the function terminates with an error. Initially input rsIDs of \eqn{i^{th}} chromosome are matched to SNPs of the annotation file of \eqn{i^{th}} chromosome. Matching entries from each of the chromosomes are extracted, merged and returned as an annotation matrix. If the code{index} is NULL, an index file((A list of choromosome wise SNPs of all ENCODE annotation files) is created and returned. For faster execution, it is recomended to create the index file for the first time and use the same in future.     
}

\value{
 It returns a R list containning the following items.
 \item{anno.mat}{Annotation matrix of input SNPs}
 \item{snp.count}{Chromosome wise summary matrix(total No. Input SNPs and Matched SNPs)}
 \item{miss.snp}{A list of chromosome wise input SNPs that are not  annotated (missing from file specified by base.path)}
 \item{miss.dbsnp}{A vector of input SNPs missing in dbSnp}
 \item{index}{R list of length 22.\code{index[[i]]} stores the all rsIDs of annotation file(specified by base.path) of \eqn{i^{th}}  chromosome}   
}

\seealso{
\code{\link{anno.merge}},\code{\link{anno.create}}
}

\references{
'https://github.com/joepickrell/fgwas/'
}

\examples{
\dontrun{

  ### Following block of codes may several minutes to complete ###
  snpfile <- system.file("sampleData", "snpData.rda", package="GKnowMTest")
  load(snpfile) ## loads snpdf
  res=create.anno.mat(rownames(snpdf),base.path="/home/Datasets/Encode/",
      fl.suffix=".annot.wdist.wcoding")
  anm=res[[1]] #anm is the annotation matrix
}

}

