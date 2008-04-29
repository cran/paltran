\name{rand.test}
\alias{rand.test}

\title{
randomisation test
}
\description{
This function computes a randomisation test to compare two different models
}
\usage{
rand.test(error1, error2, ran.nb = 1000, seed = 1)
}

\arguments{
  \item{error1}{error of model 1}
  \item{error2}{error of model 2}
  \item{ran.nb}{number of runs}
  \item{seed}{set the value of the function set.seed)
}
}
\details{

}
\value{

  \item{MSEP1 }{mean of error of model 1}
  \item{MSEP2 }{mean of error of model 1}
  \item{p }{if <0.05 both models are of significant different quality}
}
\references{
van der Voet, H. (1994) Comparing the predictive accuracy of models using a simple randomisation test.
Chemometrics and Intelligent Laboratory Systems 28: 165-180.
}
\author{
Sven Adler
}


\examples{
data(train_set.MV)
data(train_env.MV)
try<-wapls(train_set.MV,train_env.MV,val="boot")
e1<-try[[5]][,2]-train_env.MV
e2<-try[[5]][,3]-train_env.MV
rand.test(e1,e2)

}

\keyword{ models }
