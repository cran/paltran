
\name{pom}

\alias{pom}



\title{ Proportional Odds modelling for paleolimnology }

\description{
  Proportional odds modells are used to model the distribution of diatoms along environmental gradients. Than these
  models are applied to samples of a sediment core. This function is just a test version! The algorithm is slow, one run with the given data 
  takes about 5 to 10 minutes. 
}

\usage{

pom(..., d.plot = "TRUE")
}



\arguments{

  \item{\dots}{ x,y,z: required: species training set (x) as matrix and related environmental parameter (y). 
 optional: test set(z) -  species data from a sediment core }
 
 \item{d.plot}{TRUE/FALSE: if TRUE diagnostic plots are given at the end of the analysis }

}


\details{
 The relative abundances of a taxon within a data set is transformed to abundance classes (class 1 = <0 to <1% (sporadic), 
 class 2 = 1 < 5% (rarely), class 3 = 5 < 10% (regular), class 4 = 10 <30% (common), class 5 = 30 < 60% 
 ( frequent), class 6 = 60 - 100% (dominant). 
 In contrast to WA, WA-LS and MW the POM function use the abundance class specific optimum of a species,
 i.e. when a taxa occurred in a sediment core within the class 1, it might be interesting,
 where the optimum for this specific species abundance class to a related environmental factor in the training set is.
 This optimum can than be used to infer past environmental parameters. In general: 
 Not the overall species optimum is used to infer environmental parameter, but several different optima were
 estimate related to different abundance classes.  
 The pom() function uses mainly the function polr from the package MASS (Venables, W.N. Ripley, B.D. (2002).
 Modern applied statistics with S-Plus (third edition). Springer, New York, 501 pp.).
 WARNING: This version was build under R 2.5.1 - a test run with 2.7.0 leads to warning massages using the mgcv package!
 Use mgcv 1.3-31 to avoid the warning massage.
}
\value{
 
  \item{env.train }{environmental parameter of the train set}
  \item{env.test.pred }{parameter on which the prediction is done}
  \item{env.matrix.train}{density curve for the inferred environmental parameter for the train set}
  \item{inf.env.train }{inferred environmental parameter for the train set}
  \item{env.matrix.test }{density curve for the inferred environmental parameter for the test set}
  \item{reconstruction }{inferred environmental parameter for the test set}
  \item{inf.env.cross_train }{inferred environmental parameter for the train set using cross validation}
  \item{performance}{performance of the model (R2,RMSE,RMSEP,.....}
  
}
\references{Adler,S., Huebener T.,Anderson J.N., Lotter,A.F.,Modelling the distribution of diatoms along contrasting
environmental gradients in Europe using GLM, WA and Proportional Odds Models, submitted in Journal of European 
Phycology  }
\author{ Sven Adler, sven.adler2@uni-rostock.de, 
University Rostock, 
Institute for Biosciences, 
General and Systematic Botany, 
Germany}

\seealso{wa, wapls, mwtraf }
\examples{
data(dud.df)
data(age_dud)
data(train_set.MV)
data(train_env.MV)

fit1<-pom(train_set.MV,train_env.MV,dud.df)
fit1$performance
palplot(fit1$reconstruction)
palplot(fit1$reconstruction,age_dud)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ environment}
