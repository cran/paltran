\name{wapls}
\alias{wapls}

\title{weighted averaging - partial least square (WA-PLS) regression for paleoecology }
\description{This function computes with a given training set and environmental parameter a weighted averaging - partial least square (WA-PLS) transfer function as used in paleolimnology. For the calculation of the model predicting error 10 fold cross validation, bootstrap, ore Leave-on-out can bee chosen. }
\usage{
wapls(..., comp = 4, d.plot = TRUE, plot.comp = "RMSEP", env.trans = FALSE,
	spec.trans = FALSE, diagno = TRUE, seed = 1, run = 10,
	val = c("none", "10-cross", "loo", "boot"), scale =FALSE, 
	out = TRUE, drop.non.sig = FALSE, min.occ = 1)
}

\arguments{
  \item{\dots}{required x,y: a matrix or data frame of the species training set (x) and a vector or data frame of the related environmental parameter (y). optional: core samples (z) - vector or data frame of species data from a sediment core.}
  \item{comp}{number of components that will be calculated}
  \item{d.plot}{TRUE/FALSE: if TRUE diagnostic plots are given at the end of the analysis.}
  \item{plot.comp}{if "RMSEP" is chosen, the diagnostic plot for that component is given with the lowest RMESP }
  \item{env.trans}{should the environmental parameter bee transformed? "sqrt" for square root and "log10" for the logarithm to the basis 10 are possible choices, default is FALSE.}
  \item{spec.trans}{should the species data bee transformed? "sqrt" for square root and "log10" for the logarithm to the basis 10 are possible choices, default is FALSE.}
  \item{diagno}{should N2,number of non zero values bee calculated for the training set and test set? Default is TRUE }
  \item{seed}{set the seed for the random generator (using boot or 10-cross), default = 1}
  \item{run}{if "boot" or "10-cross" were chosen: number of cycles to run }
  \item{val}{validation method: one of "boot"(bootstrap), "loo"(Leave-on-out), or "10-cross"(10-fold cross validation) }
  \item{scale}{should the data scaled up to 100 percent? (Default is FALSE)}
  \item{out}{should the results printed on the console?}
  \item{drop.non.sig}{should a taxon that have non significant response to the environmental variable bee deleted? The calculation, if there is a significant relation between a taxa and the environmental variable of interest, is undertaken using a generalized additive model (GAM) and the package mgcv. As a GAM only works if a taxon occurred several times, only those taxa will be included that occurred more than 5 times (k=3). }
  \item{min.occ}{minimum occurrence: all taxa with less than min.occ will be deleted from the training set}
}
\details{
 The 10-fold cross validation is much more slower than the bootstrap or Leave-one-out, because 10 times more wapls-runs must bee performed than using e.g. bootstrap (within the same number of runs). The RMSEP of Leave one out is slightly different from C2. In this algorithm before each run of the loop the new training set (each time one sample is taken out) is controlled for zero species and removed (the same procedure as in wa, there C2 does the same). If that row is deleted from the algorithm, the results are equal for LOO. As C2 and R runs with different random numbers, the results of boottrap and 10 fold cross validation are only equal when using a high number of runs. 
}
\value{

  \item{species in train.set }{Number of non zero species in each sample of the training set}
  \item{N2 train.set  }{Hill's N2 of each sample of the training set}
\item{updated opt. }{updated optima (see reference)}  
\item{sample scores }{sample scores of the training set}  
\item{inferred train.set  }{inferred environmental parameter for the training set}
\item{performance  }{performance of the wa-pls-regression }
\item{inferred train.set.val  }{nferred environmental parameter for the training set using Leave-on-out }
\item{species in core.samples }{Number of non zero species in each sample of the core data set}
\item{n species core.samples in train.set }{How many species in the core samples are represented in the training set}
\item{N2 in core.samples  }{Hill's N2 of each sample of the core data}
\item{reconstruction_core.samples  }{reconstructed environmental parameter for the samples of the core }
\item{mean(reconstruction_core.samples).val  }{mean reconstructed environmental parameter for the samples of the core using "boot" or "loo" }
\item{sd(reconstruction_core.samples).val }{standard deviation of the reconstructed environmental parameter for the samples of the core using "boot" or "loo" }
\item{s1 (boot)}{component s1 of the bootstrap}
\item{s2 (boot)}{component s1 of the bootstrap}
\item{mean(inferred train.set).val }{mean inferred environmental variable for the training set using "boot"}
\item{sd(inferred train.set).val }{standard deviation of inferred environmental variable for the training set using "boot"}

}
\references{ter Braak, C.J.F. & Juggins, S. 1993. Weighted averaging partial least squares regression WA-PLS: an improved method for reconstructing environmental variables from species assemblages. Hydrobiologia 269:485-502. }
\author{Sven Adler}

\seealso{package analogue by G. Simpson  }
\examples{
data(train_set.MV)
data(train_env.MV)
data(dud.df)
try<-wapls(train_set.MV,train_env.MV,dud.df,val="boot")


}

\keyword{models }