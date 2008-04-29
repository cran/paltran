\name{wa}
\alias{wa}

\title{weighted averaging (WA) regression for paleolimnology}
\description{This function computes with a given training set and a given environmental parameter a weighted averaging transfer function as used in paleolimnology. For the calculation of the model predicting error 10 fold cross validation, bootstrap ore Leave-on-out can bee chosen. Inverse or classical deshrinking are supported.} 
\usage{
wa(..., d.plot = TRUE, env.trans = FALSE, spec.trans = FALSE, 
	diagno = TRUE, val = c("none", "n.cross", "loo", "boot"), 
	run = 10, scale =FALSE, seed = 1,out = TRUE, 
	desh.meth =c("class", "inverse"), 
	drop.non.sig = FALSE,min.occ = 1,nfold = 10)
}

\arguments{
  \item{\dots}{required x,y: a matrix or data frame of the species training set (x) and a vector or data frame of the related environmental parameter (y). optional: core samples (z) - vector or data frame of species data from a sediment core.}
  \item{d.plot}{TRUE/FALSE: if TRUE diagnostic plots are given at the end of the analysis. }
  \item{env.trans}{Should the environmental parameter bee transformed? Type "sqrt" for square root and "log10" for the logarithm to the basis 10 are possible choices, default is FALSE.}
  \item{spec.trans}{Should the species data bee transformed? "sqrt" for square root and "log10" for the logarithm to the basis 10 are possible choices, default is FALSE.}
  \item{diagno}{should N2,number of non zero values bee calculated for the training set and test set? Default is TRUE }
  \item{val}{validation method: one of "boot" (bootstrap), "loo" (Leave-on-out), or "n.cross" (n-fold cross validation)}
  \item{run}{ if "boot" or "n.cross" were chosen: number of cycles to run }
  \item{scale}{should the data scaled up to 100 percent? (Default is FALSE)}
  \item{seed}{set the seed for the random generator (using boot or 10-cross), default = 1}
  \item{out}{should the results printed on the console?}
  \item{desh.meth}{what kind of deshrinking method should bee used "class"(classical deshrinking), or "inverse" (inverse deshrinking), default is "inverse" }
  \item{drop.non.sig}{should a taxon that have non significant response to the environmental variable bee deleted? The calculation, if there is a significant relation between a taxa and the environmental variable of interest, is undertaken using a generalized additive model (GAM) and the package mgcv. As a GAM only works if a taxon occurred several times, only those taxa will be included that occurred more than 5 times (k=3). }
  \item{min.occ}{minimum occurrence: all taxa with less than min.occ will be deleted from the training set}
  \item{nfold}{what type of n-fold cross validation should be used?}
}

\value{
 
  \item{species in train.set }{Number of non zero species in each sample of the training set }
  \item{N2 train.set  }{Hill's N2 of each sample of the training set}
\item{species.optima }{wa-optima of each species}
\item{inferred train.set }{inferred environmental parameter for the training set}
\item{performance  }{performance of the wa-regression }
\item{species in core.samples  }{Number of none zero species in each sample of the core data set }
\item{n species core.samples in train.set }{How many species in the core samples are represented in the training set}
\item{N2 in core.samples  }{Hill's N2 of each sample of the core data  }
\item{reconstruction_core.samples }{reconstructed environmental parameter for the samples of the core }
\item{inferred train.set.val }{mean inferred environmental parameter for the training set using cross validation}
\item{mean(reconstruction_core.samples).val  }{reconstructed environmental parameter for the samples of the core using "boot" or "loo" }
\item{sd(reconstruction_core.samples).val  }{standard deviation of the reconstructed environmental parameter for the samples of the core using "boot" or "loo" }
\item{reconstruction_core.samples.val }{reconstructed environmental parameter for the samples of the core for each run of "boot" or "loo" }

}
\references{ter Braak, C.J.F. & van Dam, H. 1989. Inferring pH from diatoms: a comparison of old and new calibration methods. Hydrobiologia 178:209-23. }
\author{Sven Adler}

\seealso{package analogue }
\examples{
data(train_set.MV)
data(train_env.MV)
data(dud.df)
try<-wa(train_set.MV,train_env.MV)
try<-wa(train_set.MV,train_env.MV,desh.meth="class")
names(try)
try<-wa(train_set.MV,train_env.MV,dud.df,val="boot",run=10)
}
\keyword{models }

