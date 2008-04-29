\name{mw}
\alias{mw}

\title{Dynamic adjustment of training sets - moving windows reconstruction}
\description{The moving window method identifies for each single fossil diatom sample an optimal training set using DCA (CA, CCA) or simple distance measurement in combination with WA or WA-PLS error statistik. Downweighting of rare taxa can be chosen as non significant taxa can bee excluded }
\usage{
mw(train_set, train_env, core_data, method = c("wapls", "wa","mft"), 
	comp = 4, val = c("boot", "loo", "10-cross"), run = 100, 
	mwsize = c(20, 40, 60), dim = c(2, 3, 4), 
	mw.type = c("dca", "ca", "cca","sample"), 
	dist.m = "euclidean", rmsep.incl = TRUE, env.trans = FALSE,
	spec.trans = FALSE,
	rplot = TRUE, drop.non.sig = FALSE, min.occ = 1,scale=FALSE, dw=FALSE)
}

\arguments{
  \item{train_set}{required: matrix or data frame including species of the complete training set. rows = samples, columns = species, row and column names are required }
  \item{train_env}{required: environmental variable belongs to the training set  }
  \item{core_data}{required: species data from a core, those taxa that are not in the training set will be omitted. Minimum number is two samples.}
  \item{method}{type "wa" for weighted averaging regression or "wapls" for weighted averaging-partial least square regression. Which type of transfer function should be used to inferred the environmental variable to the core samples, default is "wapls"}
  \item{comp}{if wapls is used, how many components should be extract? Default is 4.  }
  \item{val}{validation method for the transfer function, on of "loo" for leave-on-out, "boot" for bootstrap, or "10-cross" for 10-fold cross validation, default is boot  }
  \item{run}{if "boot" or "10-cross" was chosen: how many cycles should be done? Should be low when running a new data set the first time, high values 1000 and more results in a large computing time}
  \item{mwsize}{vector of window size: how many nearest neighbours should be included? Default is 20,40,60.  }
  \item{dim}{how many dimensions should be used when the nearest neighbours are calculate using the sample scores of DCA ore CA, default is 2.}
  \item{mw.type}{type "dca" for DCA, "ca" for CA or "sample" for simple distance measurement. When choosing "dca" or "ca" the core samples are plotted in the training set samples using predict.cca or predict.decorana  (package vegan) and than the nearest training set samples to each single core sample are analysed. Using "sample" the ditances of the samples are analysed using the original species data instead of sample scores. Chosing "cca" a CCA is done and the scores of the first axis are used to analyse the nearest neighbours}
  \item{dist.m}{how to analyse the distance of the sample scores between training samples and core samples? All distances that are incorporated in vegdist (package vegan) are possible to use.}
  \item{rmsep.incl}{should the RMSEP be include in model selection or only R2.cross, mean(error).cross and max(error).cross}
  \item{env.trans}{should the environmental parameter bee transformed? "sqrt" for square root and "log10" for the logarithm to the basis 10 are possible choices, default is FALSE. }
  \item{spec.trans}{ should the species data bee transformed? "sqrt" for square root and "log10" for the logarithm to the basis 10 are possible choices, default is FALSE.}
  \item{rplot}{should a plot during the analysis be shown? Is set to bee FALSE if mw.type equals "sample" or "cca"}
  \item{drop.non.sig}{should a taxon that have non significant response to the environmental variable within the mw-traning set bee deleted? The calculation, if there is a significant relation between a taxa and the environmental variable of interest, is undertaken using a generalized additive model (GAM) and the package mgcv. As a GAM only works if a taxon occurred several times, only those taxa will be included that occurred more than 5 times (k=3). If the mwsize is too smal, it can happens, that no taxa have a significant response and the function stops}
  \item{min.occ}{minimum occurrence: all taxa with less than min.occ will be deleted from the training set}
  \item{scale}{should the data scaled up to 100 percent? (Default = FALSE)}
  \item{dw}{should rare taxa be downweighted? (see function downweight in the vegan package by J. Oksanen)}
}
\details{
Using mw, for each sample 3 WA-PLS runs (default) are calculated using
100 bootstrap runs for each. This takes time. The reconstruction for a whole
sedimet (80-100 samples) core can take several minutes. Please try first with
a small test set or with a low value for run (see examples), before running the whole reconstruction! At least the number of componentes (default is 4) can bee reduced to make the function faster. Data must be
organised in the same way as running cca or dca in package vegan. This approach
needs large training sets, like the combinded TP data set from EDDI. 
}

\value{
  \item{sample.performance }{gives all information for each core sample, which window size was used and the performance of the related transfer function }
  \item{reconstruction }{reconstructed values for the core sample}
  \item{mean(reconstruction).val }{mean values for the reconstruction for the core sample using bootstrap or 10-fold cross validation}
  \item{sd(reconstruction).val }{standard deviation of the reconstructed values for the core sample using bootstrap or 10-fold cross validation}
}
\references{Huebener, T., Dressler,M., Schwarz,A.,Langner, K., Adler,S. 2008. Dynamic adjustment of training sets (`moving windows` reconstruction) by using transfer functions in paleolimnology -a new approache, J. o. Paleolimnology 40: 79-95}
\author{ Sven Adler}



\seealso{wa,wapls, package analogue (G. Simpson and J Oksanen) and package vegan (J. Oksanen)}


\examples{
data(dud.df)
data(train_set.MV)
data(train_env.MV)
test<-dud.df[1:3,]
fit<-mw(train_set.MV,train_env.MV,test,mwsize = c(40,60),val="boot",run=5,comp=3)
names(fit)

fit<-mw(train_set.MV,train_env.MV,test,mwsize = c(40,60,80)
	,comp=3,method="wa",val="loo",)

fit<-mw(train_set.MV,train_env.MV,test,mwsize = c(40,60),run=5,
    	mw.type="sample",dist.m="bray",dw=TRUE)
}

\keyword{ models }

