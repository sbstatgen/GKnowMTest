\name{GKnowMTest-package}
\alias{GKnowMTest-package}
\alias{GKnowMTest}
\docType{package}
\title{
\packageTitle{GKnowMTest}
}
\description{
Tools for integrating prior knowledge into genomic studies. 
}
\details{
Currently the package provides two functions for integrating prior knowledge into SNP-level p-values (e.g. from a GWAS). One is\code{\link{anno.create}} and other is a statistical function \code{\link{prior.adjust}}.First function creates a Mapping Object of class \code{\link{annotatedSNPset}}.It supports two kinds of priors on SNPs, namely SNP-level functional annotations (e.g. from ENCODE) and gene-level annotations (e.g. pathways).SNP-level functional annotations can also be obtained using another function \code{\link{create.anno.mat}}. A mapping object needs to be first created using the \code{\link{anno.create}} or \code{\link{anno.merge}}functions before running the statistical routines.
The statistical function \code{\link{prior.adjust}} integrates prior knowledge into SNP-level p-values (e.g. from a GWAS).The adjusted p-values of the SNPs are calculated based on the prior probability of being associated to the disease. The prior association probabilies are calculated either by Method of Moments or by Regression, as specified by the user. It is incorporated as SNP groups such as functional annotations or pathway information as obtained from the mapping object. The adjusted p-values can be calculated by quantile adjustment or by simple or optimal weighting methods, as specified by the user.  
}

\references{
Pickrell, Joseph K. "Joint analysis of functional genomic data and genome-wide association studies of 18 human traits." The American Journal of Human Genetics 94.4 (2014): 559-573. 
}

\keyword{ SNP, GWAS, Prior Knowledge, Pathways, ENCODE, Annotation }


