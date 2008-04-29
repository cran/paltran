\name{mtf}
\alias{mtf}

\title{ multinomial transfer function for Palaeolimnology}
\description{ Here the species abundances were assigned to different
abundance classes, for each class a species response curve is generated using GLM. These response
curves are used to infer the environmental parameter of interest (s. Details)}
\usage{
mtf(..., d.plot = TRUE, scale = TRUE, datatype = TRUE,
	n= 4,set.zero = FALSE, not.av = c("zero", "lower", "max"),
	val = "loo",class=c(1,5,10,30,60),out=TRUE)
}

\arguments{
  \item{...}{x,y,z: required: species training set (x) as matrix and related environmental parameter (y). optional: core data (z) - species data from a sediment core }
  \item{d.plot}{ if TRUE diagnostic plots are given at the end of the analysis}
  \item{scale}{should the data scaled up to 100 percent? (Default)}
  \item{datatype}{What kind of data input? TRUE for percentage data, otherwise no transformation to class date will bee done}
  \item{n}{minimum number of species occurrence within on class } 
  \item{set.zero}{If a species occurs within on class less than n time the values will be set as zero, if set.zero = TRUE. If set.zero = FALSE, the values will be allocate to the next lower class }
  \item{not.av}{If a species occurs within a abundance class that is not given in the training set, three possibilities are available. "zero" means that this species will be excluded, "lower" stands for to use the next lower class density. If "max" is chosen, the cumulative species response curve of all classes will be calculate and use}
  \item{val}{ validation method: So fare only Leave-one-out is performed }
  \item{class}{what classes should be used? Default is 0,1,5,10,30,60}
  \item{out}{should the result bee showen on the console?}
}
\details{
The relative abundances of a taxon within a data set is transformed to abundance classes (class 1:
0<x<=1,class 2: 1<x<=5, class 3: 5<x<=10, class 4 10<x<=30, class 5 30<x<=60, class 6: 60-100 percent abundance, number of classes and class borders can bee changed). Occurs a taxa in a sediment core sample within the 
class 1, it might give better reconstructions, when the optimum for this specific species abundance 
class to a related environmental factor in the training set is used instead of the optimum of the hole 
response. This optimum can than be used to infer past environmental parameters. In general: Not
the overall species optimum is used to infer environmental parameter, but several different optima
were estimate related to different abundance classes using GLM with binomial distribution. The
inferred value for the core sample is than calculated following the ML-Method: multiplying the class 
density functions of the species abundance class values from the core. The first version of this function
(paltran 1.0) uses mainly the function polr from the package MASS (Venables, W.N. Ripley, B.D. (2002). 
Instead, each class is modelled 
using GLM. If a species occurs only rarely (<n) in the training set tow options are possible. First, it can
removed (set.zero=TRUE), secondly the occurrence can add to the next lower class. On one hand this
is a misclassification, on the other hand, the data are not completely lost, and might improve the
model. If a species class from a core sample is not represent in the training set, the same options can
be chosen. The species can bee removed (not.av="zero"), or the next lower density function can bee
used. Thirdly the overall density function (product of all single class density functions) - assuming
a symmetric species distribution - can be chosen (not.av="max"). As validation method until now
only Leave One Out is available, as the algorithm is still very slow. With improving computing time,
other methods like bootstrap might be incorporated.

}
\value{
 
  
  \item{inferred train.set }{inferred environmental parameter for the training set}
  \item{performance}{performance of the pom-regression }
  \item{spec distribution}{species distribution curves for the single classes}
  \item{inferred train.set (loo) }{inferred environmental parameter for the training set using leave one out as cross validation method}
  \item{reconstruction_core.samples}{reconstructed environmental parameter for the samples of the core }
  \item{mean (reconstruction_core.samples) (loo)}{reconstructed environmental parameter for the samples of the core using leave one out}
  \item{sd (reconstruction_core.samples) (loo)}{standard deviation of the reconstructed environmental parameter for the samples of the core using "loo"}

}
\references{ in preperation }
\author{ Sven Adler }

\examples{
data(train_set.MV)
data(train_env.MV)
data(dud.df)
try<-mtf(train_set.MV,train_env.MV,dud.df[1:3,],not.av="max",val=FALSE)

#using "loo" takes several minutes for computing


}

\keyword{ models }
