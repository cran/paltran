\name{dud.df}
\alias{dud.df}
\docType{data}
\title{Diatom samples of a sediment core - lake Dudinghausen}
\description{Lake Dudinghausen is a recent eutrophic lake.  
The lake morphology of lake Dudinghausener (DUD) (53.54.30 N, 12.12.40 E)
is characterized by a nearly homogeneously shaped basin. The lake is almost completely 
surrounded in agricultural used watershed, and only a small part of the northern shore is forested. 
Large areas of cattail are located at the western and eastern edges. Species codes are following the EDDI code
http://craticula.ncl.ac.uk/Eddi/jsp/index.jsp
}
\usage{data(dud.df)}
\format{
 data.frame
}


\references{
 Dressler, M., Selig, U., Doerfler, W., Adler, S., Schubert, H., Huebener, T. 2006.
  Environmental changes and the migration period in Europe by the example of Lake Dudinghausen 
  northern Germany. Quat. Res. 66:25-37.
}
\examples{
data(age_dud)
data(dud.df)
data(train_set.MV)
data(train_env.MV)

fit1<-wa(train_set.MV,train_env.MV,dud.df)

}
\keyword{datasets}
