\name{mwtraf}
\alias{mwtraf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Dynamic adjustment of training sets - moving windows reconstruction
}
\description{The moving window method identifies the best number of nearest neighbours (window size) for each fossil diatom assemblage and the 
best  transfer function based on the error statistic for each sample of the test set (e.g. sediment core) 
 
}
\usage{
mwtraf(..., method = "wapls", rplot = "TRUE", mwsize = c(40, 60, 80, 100, 120, 140, 160), mw_type = "dca", ncomp = 4, rmsep.incl = TRUE, env.trans = FALSE, spec.trans = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{ x,y,z: required: species training set (x) as matrix and related environmental parameter (y). 
 optional: test set(z) -  species data from a sediment core. For all data: row- and columnnames are required! Otherwize the 
 sample specific training set and the related environmenal parameter can not be found! }
  \item{method}{currently only wapls is running: choosing ncomp=1 a simple wa is done, for each sample specific training set.}
  \item{rplot}{TRUE/FALSE should the results plotted?}
  \item{mwsize}{vector of the number of samples within each window: default runs only with large data sets,
 e.g. combined TP data set from the EDDI data base}
  \item{mw_type}{How should the nearest neighborous be found? curently only a dca is running, chord distance
 and NMDS is in preperation.}
  \item{ncomp}{ number of components that should bee extract using wa-pls}
  \item{rmsep.incl}{should the RMSEP used for modelselection?}
  \item{env.trans}{ should the environmental parameter bee transformed? "sqrt" for square root
	and "log10" for logarithm to the basis 10 are possible choices, default is FALSE.}
  \item{spec.trans}{should the species data bee transformed? "sqrt" for square root
	and "log10" for logarithm to the basis 10 are possible choices, default is FALSE. }
}
\details{
WARNINGS: 
One run of a WA-PLS with a normal sample size takes about 20 seconds. Using mwtraf, for each sample 7 WA-PLS runs (default) 
are calculated, that takes 140 seconds. The reconstruction for a whole sedimet (80-100 samples) core can take 5 h and more! Please try first wih a smal test set, before 
running the whole reconstruction!
Data must be organised in the same way as running da cca or dca in package vegan! (see examples)
This approach needs large training sets, like the combinded TP data set from EDDI.
}
\value{
  
  \item{sample.performance }{perforamnce of the sample specific wa-pls model}
  \item{reconstruction }{reconstructed environmental parameter}
 
}
\references{ Huebener, T., Dressler,M., Schwarz,A.,Langner, K., Adler,S. 2007. Dynamic adjustment of training sets (`moving windows` reconstruction) by using
 transfer functions in paleolimnology -a new approache, J. o. Paleolimnology. DOI 10.2007/s10933-007-9145-7}
\author{Sven Adler, sven.adler2@uni-rostock.de, 
University Rostock, 
Institute for Biosciences, 
General and Systematic Botany, 
Germany}
\note{choosing mwsize=30 and ncomp =1 this approach equals the locally weighted WA regression (Battarby et al. 2006). 

 
}
\seealso{wa, wapls, wa package analogue}
\examples{

data(dud.df)
data(train_set.MV)
data(train_env.MV)

fit1<-mwtraf(train_set.MV,train_env.MV,dud.df[1:3,],mwsize = c(40, 60, 80))
fit1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ environment }
