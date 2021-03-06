% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saveFigure.R
\name{saveFigure}
\alias{saveFigure}
\title{Save Figure}
\usage{
saveFigure(fileN)
}
\arguments{
\item{fileN}{File name witout extension. Default "default".}
}
\value{
Save a figure uder the current directory. If the file name exist an _k will be added and k will increment by one on each call of saveFigure.
}
\description{
Save Figure
}
\examples{
saveFigure()
saveFigure("populationStructure")
saveFigure("inbreedingCoefficient")
saveFigure("test")
}

