\name{wapls}
\alias{wapls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ weighted averaging - partial least square (wa-pls) regression for paleoecology }
\description{
  This function computes with a given training set and environmental parameter a weighted averaging -
partial least square (WA-PLS) transfer function as used in paleolimnology. For error estimation a 10 fold cross validation procedure
is choosen: For large data sets (and a high number of components) the calculation time will be 25 seconds and more. The last column of the
output of the performance is needed for MW, for simple WA-PLS this can bee ignored.)
}
\usage{
wapls(..., n_comp = 4, d.plot = TRUE, plot.comp = "RMSEP", env.trans = FALSE, spec.trans = FALSE, diagno = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{x,y,z: required: species training set (x) as matrix and related environmental parameter (y). 
optional: test set(z) -  species data from a sediment core }
  \item{n_comp}{ number of components that should bee extract}
  \item{d.plot}{TRUE/FALSE: if TRUE diagnostic plots are given at the end of the analysis  }
  \item{plot.comp}{"RMSEP": the diagnostig plots are done automatical for that component with the lowest RMSEP}
  \item{env.trans}{should the environmental parameter bee transformed? "sqrt" for square root
and "log10" for logarithm to the basis 10 are possible, choices, default is FALSE}
  \item{spec.trans}{should the species data bee transformed? "sqrt" for square root
and "log10" for logarithm to the basis 10 are possible, choices, default is FALSE}
  \item{diagno}{should N2,number of non zero values bee calculated for the training set and test set? }
}

\value{
 
  \item{non_z_train }{Number of non zero species in each sample of the training set}
  \item{N2_train }{Hill's N2 of each sampel of the training set}
  \item{scores }{site scores of the training set samples}
  \item{inferred.env_train }{inferred environmental parameter for the training set}
  \item{non_z_test }{Number of non zero species in each sample of the test set}
  \item{N2_test}{Hill's N2 of each sampel of the test set}
  \item{reconstruction }{reconstructed environmental parameter for the samples of the core}
  \item{performance}{performance of the wa-pls regression}
  
}
\references{ter Braak, C.J.F. & Juggins, S. 1993. Weighted averaging partial least squares regression WA-PLS: 
an improved method for reconstructing environmental variables from species assemblages.
Hydrobiologia 269:485-502.}
\author{Sven Adler, sven.adler2@uni-rostock.de
University Rostock, 
Institute for Biosciences, 
General and Systematic Botany, 
Germany}


\seealso{wa, mwtraf, pom, package analogue (Simpson, 2008) for wa and MAT }
\examples{
data(age_dud)
data(dud.df)
data(train_set.MV)
data(train_env.MV)

fit1<-wapls(train_set.MV,train_env.MV,dud.df)
names(fit1)
palplot(fit1$reconstruction[,1],age_dud)
palplot(fit1$reconstruction[,1],age_dud,trans="log10")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{environment }
