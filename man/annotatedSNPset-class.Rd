\name{annotatedSNPset-class}
\docType{class}
\alias{annotatedSNPset}
\alias{annotatedSNPset-class}


\title{annotatedSNPset objects}

\description{
This class is for data objects of annotatedSNPset-class
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{annotatedSNPset(snp.df,snp.eq,eq.mat,eq.map,dim)}:
      Creates an object of annotatedSNPset class.
      \describe{
        \item{\code{snp.df}}{
          A data.frame with input rsIDs as rownames and both chromosome number and base pair position as columns.
        }
        \item{\code{snp.eq}}{
         This is a numeric vector indicating which snp mapped to which  equivalence class.
        }
        \item{\code{eq.mat}}{
         A list of binary matrix or matrices. Dimension of i^{th} matrix is no.of equivalence class by no.of annotations of i^{th} object.
        }
        \item{\code{eq.map}}{
         Either NULL or an integer matrix. It is NULL for a simple annotatedSNPset object. If more than one annotatedSNPset objects are merged into a composite object then this field holds a map matrix. Number of columns of the map matrix is the number of objects that are merged and the number of rows is the total number of new equivalence classes after merging.
        }
        \item{\code{dim}}{
         An R list of length 3. First element is a numeric vector of dimensions (number of SNPs \code{no_snps}, total number of equivalence classes \code{no_eq_class},total number of annotation \code{no_anno} and total number of merged objects \code{used_ob}. Other two elements are either NULL(for simple object) or numeric vectors of length equal to the number of objects used for merging i.e. \code{used_ob}. These two vectors hold the number of rows and number of columns respectively of the binary matrices in the list \code{eq.mat}.
        }  
      }
    }
  }
}


\section{Accessors}{
In the code snippets below, \code{x} is a annotatedSNPset object.
\describe{
\item{}{\code{getEQ(x)}:To extract 'snp.eq' slot}
\item{}{\code{getDF(x)}:To extract 'snp.df' slot}
\item{}{\code{getDIM(x)}:To extract first element 'dim' slot}
\item{}{\code{getMAT(x)}:To extract equivalence class by annotation matrix}
}
}


\examples{

  obfile <- system.file("sampleData", "obData.rda", package="GKnowMTest")
  load(obfile) # loads a mapping object ob1 
  ob1

  res=getMAT(ob1) # use of acessor function 'getMAT'
  res[1:10,1:6]
  
  res=getDF(ob1)
  head(res)
  
  res=getDIM(ob1)
  res
  
  res=getEQ(ob1)
  head(res,20)

}






