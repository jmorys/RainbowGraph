% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{simulate_expansion_of_fraction}
\alias{simulate_expansion_of_fraction}
\title{Simulate regeneration with each cell}
\usage{
simulate_expansion_of_fraction(
  graph,
  survivor_fraction = 0.3,
  expanding_fraction = 0.1
)
}
\arguments{
\item{graph}{undirected graph with edge weight parameter}

\item{survivor_fraction}{fraction of cells that survive}

\item{expanding_fraction}{fraction of cells among survivors that multiply and drive regeneration}
}
\value{
list containing vectors of surviving cells, multiplying cells, cells after expansion without non_multiplying cells and finally with non_multiplying cells.
The results of this function reference unique ids of cells, not colors.
}
\description{
Simulate regeneration with each cell
}
\examples{

}
