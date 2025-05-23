% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{convert_to_decimal}
\alias{convert_to_decimal}
\title{Convert degrees, minutes, and seconds to decimal coordinates}
\usage{
convert_to_decimal(data)
}
\arguments{
\item{data}{A data frame with columns: \code{degreeLat}, \code{minLat}, \code{secLat}, \code{degreeLon}, \code{minLon}, \code{secLon}, \code{cardinalDirectionLat}, \code{cardinalDirectionLon}.}
}
\value{
A data frame with added \code{decimalLatitude} and \code{decimalLongitude} columns.
}
\description{
Converts occurrence coordinates stored in DMS format to decimal degrees.
}
