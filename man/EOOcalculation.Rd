% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{EOOcalculation}
\alias{EOOcalculation}
\title{Calculate Extent of Occurrence (EOO)}
\usage{
EOOcalculation(map)
}
\arguments{
\item{map}{An \code{sf} object representing a species' range.}
}
\value{
Numeric. EOO in km².
}
\description{
Uses a minimum convex polygon (MCP) to compute EOO area in km².
}
