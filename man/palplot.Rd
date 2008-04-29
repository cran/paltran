\name{palplot}
\alias{palplot}
\title{plotting a simple plot for paleolimnology}
\description{Before exporting the results of wa, wapls or mw functions to C2 or TILIA simple stratigraphic plots can bee done with palplot including several variables to cheque the results. For each variable the colour and the plot style can be chosen, error bars and trend lines can be added.}
\usage{palplot(..., age = NA, error.plot = FALSE, error = 0, ptype = "b",
	error.ptype = "b", p.col = "black", 
	p.lty = 1, p.xlab = NULL, cex.main = 1,
	cex.axis = 1, dis = 0.15, p.main = NULL,
	y.lab = "sample", trend = 0, span = 0.1)}

\arguments{
  \item{...}{required: matrixes, vectors or data frames of species data and or reconstructed environmental parameters}
  \item{age}{if available: a vector including the age of the samples}
  \item{error.plot}{should error bars be included?}
  \item{error}{vector or matrix included the error}
  \item{ptype}{how should the variables be plotted? String or vector: "h" for bars, "p" for points, "l" for line, "b" for both line and points}
  \item{error.ptype}{how should the errors be plotted? String or vector: "b" for error bars, "p" for shaded regions}
  \item{p.col}{string or vector: colour of the plots}
  \item{p.lty}{if "l" is chosen for p.type, the line type can be chosen. string or vector }
  \item{p.xlab}{string or vector: labels of the x axis of each variable }
  \item{cex.main}{numeric or vector: size of the main of the plots for each variable }
  \item{cex.axis}{numeric or vector: size of x axis labels }
  \item{dis}{distance between the single plots}
  \item{p.main}{string or vector: main of the plots for each variable }
  \item{y.lab}{string: label of the y axis }
  \item{trend}{should a trend line be plotted? String or vector}
  \item{span}{what span should be used for the trend lines}
}
\author{ Sven Adler}

\seealso{ package analogue (G. Simpson) or package palaeo (S. Juggins) }
\examples{
data(age.dud)
data(dud.df)
data(train_set.MV)
data(train_env.MV)
try<-wa(train_set.MV,train_env.MV,dud.df,val="boot",run=10)     

spec<-c(23,24,217,218,223,228)
c<-c(rep(4,6),2,3)
t<-c(rep("h",6),"l","p")
x<-c(rep("abundance",6),"n","mug/l")
tr<-c(rep(0,7),1)
palplot(dud.df[,spec],try[[12]],10^try[[13]],trend=tr,
	cex.axis=0.8,p.col=c,ptype=t,p.main=c(names(dud.df[,spec]),
	"N2","inf. log(TP)"),dis=0.2,p.xlab=x)

palplot(try[[14]],ptype="l",error=1.96*try[[15]],
	error.ptype="p",error.plot=TRUE,p.col="blue",age=age.dud[,1])                        
}

\keyword{hplot}

