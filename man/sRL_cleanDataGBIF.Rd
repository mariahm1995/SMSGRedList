% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{sRL_cleanDataGBIF}
\alias{sRL_cleanDataGBIF}
\title{Clean GBIF occurrence records based on quality filters}
\usage{
sRL_cleanDataGBIF(
  flags,
  year_GBIF,
  uncertainty_GBIF,
  GBIF_xmin,
  GBIF_xmax,
  GBIF_ymin,
  GBIF_ymax,
  AltMIN,
  AltMAX,
  land
)
}
\arguments{
\item{flags}{A data frame of occurrence records to filter.}

\item{year_GBIF}{Minimum acceptable year (e.g., 1800).}

\item{uncertainty_GBIF}{Maximum acceptable uncertainty in km.}

\item{GBIF_xmin, GBIF_xmax, GBIF_ymin, GBIF_ymax}{Coordinate bounds for filtering.}

\item{AltMIN, AltMAX}{Elevation filtering limits (in meters).}

\item{land}{Logical. If TRUE, restricts to continental areas.}
}
\value{
A data frame with quality flags and a \code{Reason} column explaining rejections.
}
\description{
Applies temporal, spatial, uncertainty, elevation, and land/sea filters to flag invalid records.
}
