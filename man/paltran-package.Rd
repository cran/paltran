\name{paltran-package}
\alias{paltran-package}
\alias{paltran}
\docType{package}
\title{
WA, WA-PLS and MW for Paleolimnology }
\description{
This package contains functions for paleolimnology -
        WA-regression (see also package analogue!), WA-PLS and MW
        regression. The function palplot allows a first plot of the
        data including a trend analysis. Additional a transfer function based on proportional odds (POM) is added, but this function is just a
        first test version based on class data.
}

\details{
\tabular{ll}{
Package: \tab paltran\cr
Type: \tab Package\cr
Version: \tab 1.3-0\cr
Date: \tab 2010-02-06\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

}
\author{
Sven Adler <sven.adler2@uni-rostock.de>}
\references{

ter Braak, C.J.F. & van Dam, H. 1989. Inferring pH from diatoms: a comparison of old and new calibration methods. Hydrobiologia 178:209-23. 
ter Braak, C.J.F. & Juggins, S. 1993. Weighted averaging partial least squares regression WA-PLS: an improved method for reconstructing environmental variables from species assemblages. Hydrobiologia 269:485-502. 
Huebener, T., Dressler,M., Schwarz,A.,Langner, K., Adler,S. 2008. Dynamic adjustment of training sets (`moving windows` reconstruction) by using transfer functions in paleolimnology -a new approache, J. o. Paleolimnology 40, 79-95
}

\keyword{ package }
\seealso{
package analogue by G. Simpson and J. Oksanen
}
\examples{

data(age.dud)
data(dud.df)
data(train_set.MV)
data(train_env.MV)
try<-wa(train_set.MV,train_env.MV,dud.df,val="boot",run=10)     
# run should be set higher 100-1000

spec<-c(23,24,217,218,223,228)
c<-c(rep(4,6),2,3)
t<-c(rep("h",6),"l","p")
x<-c(rep("abundance",6),"n","mug/l")
tr<-c(rep(0,7),1)
palplot(dud.df[,spec],try[[12]],10^try[[13]],trend=tr,
	cex.axis=0.8,p.col=c,ptype=t,p.main=c(names(dud.df[,spec]),
	"N2","inf. log(TP)"),dis=0.2,p.xlab=x,age=age.dud[,1])

palplot(try[[14]],ptype="l",error=1.96*try[[15]],
	error.ptype="p",error.plot=TRUE,p.col="blue",age=age.dud[,1])  

}
