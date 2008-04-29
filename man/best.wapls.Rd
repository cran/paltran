\name{best.wapls}


\alias{best.wapls}


\title{model selection of best WA-PLS component}



\description{
Required in MW transfer function - model selection of best WA-PLS component: Not only the RMSEP is used: A combination of
R2, min(error), max(error) and RMSEP is used for model election, as the RMSEP alone might be problematic. The component with the 
highest R2 (lowest min(error),...) gets 2 points, the second highest 1 point (column MW). That component with the highest
number of points (MW) will bee chosen for the reconstruction.  The RMSEP can bee excluded compleetly from the modelselection.
}



\usage{
best.wapls(error, nro = 4,rmsep.incl=TRUE)
}



\arguments{
  \item{error}{error statistic from a wa-pls.}
  \item{nro}{number of components, must be the same as in wapls(). }
  \item{rmsep.incl}{should the RMSEP included in the model selection?}
}

\value{Error statistic of a WA-PLS run ordered by XR2, XAveBias, XMaxBias ,  RMSEP. The component with the highest MW value is the best.
}

\references{Huebener, T., Dressler,M., Schwarz,A.,Langner, K., Adler,S. 2007. Dynamic adjustment of training sets (`moving windows` reconstruction) by using
 transfer functions in paleolimnology -a new approache,J. o. Paleolimnology. DOI 10.2007/s10933-007-9145-7.
}

\author{Sven Adler, sven.adler2@uni-rostock.de,
 University Rostock, 
 Institute for Biosciences, 
 General and Systematic Botany, 
 Germany 
}

\seealso{ mw, wapls, wa and package analogue }


\examples{
data(train_set.MV)
data(train_env.MV)
data(dud.df)
try1<-wapls(train_set.MV,train_env.MV,dud.df)
try1$performance
best.wapls(try1$performance)
}

\keyword{environment  }

