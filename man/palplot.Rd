\name{palplot}
\alias{palplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting a simple stratigraphic plot for paleoecology}
\description{
  Plotting a simple stratigraphic plot for paleoecology. To analyse the reconstructed environmental parameters a trend line 
can bee drawn. 
}
\usage{
palplot(..., col = "red", trend = TRUE, trendcol = "gray90", lty = 1, lwd = 1, k = FALSE, polyg = TRUE, colpoly = "gray90", trans = FALSE, meth = "loess", span = 0.2, conf.t = TRUE, xlab = "inferred env. factor", ylab = "sample", main = "reconstruction", cex = 1.5, pch = 20, type = "p")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{x,y: x (required) - inferred environmental parameter, y (if in hand) age of the samples}
  \item{col}{ color of points}
  \item{trend}{TRUE/FALSE - should a trend line bee drawn}
  \item{trendcol}{ color of the trend line}
  \item{lty}{ line type of the trend line,see ?plot for more information  }
  \item{lwd}{ line width of the trend line,see ?plot for more information  }
  \item{k}{if a GAM is used for the trend: number of knots }
  \item{polyg}{Should the trend shaded?}
  \item{colpoly}{color of the shaded trend }
  \item{trans}{ are the data tranformed? Here they can be backtransformed: "log10" for log10 transformed data,
 "log" for ln transformed data and "sqrt" for sqrt transformed data }
  \item{meth}{ method for analysing the trend: "loess" or "gam" }
  \item{span}{ span of the loess function }
  \item{conf.t}{should a confidence intervall bee calculated for the trend? }
  \item{xlab}{label x axis}
  \item{ylab}{label x axis  }
  \item{main}{label main}
  \item{cex}{size of symbols}
  \item{pch}{symbols, see ?plot for more information  }
  \item{type}{"p" for points, "l" for line, "b" for both}
}

\value{
  
  \item{env.inferred }{inferred environmental parameter}
  \item{env.trend }{fitted values of the estimated trend}
  \item{trend.conf.up }{fitted values of the estimated confidence intervall}
\item{trend.conf.un }{fitted values of the estimated confidence intervall}
\item{summary.trend }{summary of the trend}

}

\author{Sven Adler, sven.adler2@uni-rostock.de, 
University Rostock, 
Institute for Biosciences, 
General and Systematic Botany, 
Germany  }


\examples{
data(age_dud)
data(dud.df)
data(train_set.MV)
data(train_env.MV)

fit1<-wa(train_set.MV,train_env.MV,dud.df)
palplot(fit1$reconstruction,age_dud)
palplot(fit1$reconstruction,age_dud,trans="log10",main="TP for lake Dudinghausen", ylab="TP")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{environment }
