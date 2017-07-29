\name{prior.adjust}
\alias{prior.adjust}

\title{
Function for calculating adjusted p-values based on prior knowledge on SNPs.
}

\description{
This function calculates adjusted p-values based on some prior knowledge on  the likelihood of SNPs to be associated. Currently prior knowledge can only be specified as groups of SNP sharing certain properties such as functional annotations or sets of pathways that they map to. The prior is specified as an \code{\link{annotatedSNPset}} object. Either Z-scores or raw p-values are taken as input, and prior-adjusted p-values are returned that can be compared to the usual genome-wide (e.g. 5e-08) or other multiple testing cutoffs.
}

\usage{
prior.adjust(data.df,map.obj,prior.meth=c("MoM","Reg"),
adj.meth=c("quantile","opt.wt","simple.wt","quad.wt","pair.wt"),
control.prior=list(alpha=1,lambda=NULL,alt.dist="norm",use.cv=FALSE),
control.adj=list(level=0.05),verbose=TRUE,
control.locfdr=list(plot=FALSE),control.pair=list(npth.sc=-1))
}

\arguments{
  \item{data.df}{
A data frame containing a z-score column with column name "zscore" and/or p-value column with column name "pvalue".. If a z-score column exists, the p-value column is ignored (p-values are recalculated assuming these are 2-tailed Z-scores. The rownames attribute of the data frame must have rsID-s of the SNPs.
}
  \item{map.obj}{
An \code{\link{annotatedSNPset}} object that maps input SNPs to groups derived from prior knowledge. This object is usually obtained as an output of the \code{\link{anno.create}} function.
}
  \item{prior.meth}{
A character string specifying the method using which the prior probabilities are to be calculated. For Method of Moments "MoM" and for penalized logistic regession "Reg". Three options are .... See Details.
}
  \item{adj.meth}{
A character string specifying the method of calculation of adjusted p.values. See Details.
}
  \item{control.prior}{
A list containing the arguments to be passsed to prior calculation function specified by \code{prior.meth} . The arguments "alpha" and "lambda" and "fit" are required for the "Reg" option. The first two are directly passed to the \code{\link{glmnet}} function for penalized logistic regression. "alt." is a character string denoting the model to fit before running the penalized regression. It can be either "norm" (non-parametric mixture estimation as in package \code{\link{locfdr}}) or "npar" (using package \code{\link{locfdr}} to get the prior probabilities and non-parametically estimate the alternate cdf) or "mixnorm" (mixture normal density as in package \code{mixtools}).
}
  \item{control.adj}{
A list containing the arguments to be passed to  the p-value adjustment function specified by \code{adj.meth}. The nominal level of significance  \code{level} (a number between 0 and 1) is only required when "opt.wt" is selected. Default is 0.05.
}
  \item{verbose}{
A logical operator. If set to TRUE, messages indicating progress are printedand output is turned off if  set to FALSE. Default is TRUE.
}
  \item{control.locfdr}{
Further arguments to be passed to the fitting function \code{\link{locfdr}}. By default (if nothing is specified), default options of "locfdr" are used with \code{plot} set to FALSE. If locfdr fails with an error, parameters (e.g. degrees of freedom \code{df}) may need to be changed through this list.
}
\item{control.pair}{
Argument to be passed to the p.value weighting function \code{pair.wt}. By default (if nothing is specified), default option of -1 is used.
}
}

\details{

The function proceeds in two steps. The first step is deriving the "prior" probability of SNPs to be associated. The two-tailed sided z-.scores (or p-values) of SNPs are taken as input and converted to 1-tailed Z-scores. Also, the groups (equivalance classes) to which they SNPs map are derived from the input \code{map.obj}. Currently there are two alternative methods for obtaining the prior probabilities of each group. These are calculated either using the method-moments approach of Roeder and Wasserman (2009) or using a new (penalized-) regression approach (Biswas et. al., 2016). In the first approach, a simple mixture-normal model is used, variance of moment estimates are not accounted, so this approach may suffer when groups are small. In the second approach, posterior probabilities of the SNPs (obtained by fitting locfdr or mixnorm models) on the input Z-scores are regressed on the annotations using the \code{\link{glmnet}} package.

The second step is "adjustment" of p-values based on the prior probabilities or 'proportion of truths' in each group as derived in the first step. Currently, 3 alternative methods are supported, namely "quad.wt" (cubic weighting), "opt.wt" (Optimal-weighting), and "pair.wt" (pair-wise weighting). These methods are discussed in Biswas et. al. (2016). 'Optimal weighting' tries to maximize the prior subject to type-1 error constraint. 'Pair-wise weighting' is more like opt.wt but sequencially optimizes weight for a pair of SNP groups. 'Cubic weighting' derives weight as a cubic obtimization function by optimizing power.

}

\value{
The function returns a list containing a Data Frame of the names of SNPs as rownames, their two-tailed z-scores (if given as input), raw p-values,prior-adjusted p-values bonferroni adjusted pvalues, and prior-adjusted local-FDR values, prior 'probabilities of truth', alternate cdf and the weights assigned to the SNPs.
}

\examples{

snpfile <- system.file("sampleData", "snpData.rda", package="GKnowMTest")
obfile <- system.file("sampleData", "obData.rda", package="GKnowMTest")

load(snpfile) ## loads snp
load(obfile) ## loads mapping object

res.opt <- prior.adjust(snpdf,map.obj=ob1,prior.meth="Reg",adj.meth="pair.wt")
head(res.opt)

}

\references{
Roeder, Kathryn, and Larry Wasserman. "Genome-wide significance levels and weighted hypothesis testing." Statistical science: a review journal of the Institute of Mathematical Statistics 24.4 (2009): 398.

Biswas et. al. (2016, In preparation): "Pathway-guided Search Improves the Power to Discover SNPs from Genome-wide Association Studies."

Efron, Bradley, Brit Turnbull, and Balasubramanian Narasimhan. "locfdr: Computes local false discovery rates." R package (2011).
}

