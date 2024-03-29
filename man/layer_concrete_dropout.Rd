% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{layer_concrete_dropout}
\alias{layer_concrete_dropout}
\title{Layer wrapper with concrete dropout}
\usage{
layer_concrete_dropout(
  object,
  layer,
  weight_regularizer = 1e-06,
  dropout_regularizer = 1e-06,
  init_min = 0.1,
  init_max = 0.02,
  is_mc_dropout = T,
  supports_masking = TRUE
)
}
\arguments{
\item{layer}{layer around which wrapper is built}

\item{weight_regularizer}{importance of weight regularizer in final loss - float}

\item{dropout_regularizer}{importance of dropout regularizer in final loss - float}

\item{init_min}{minimum value of dropout probability}

\item{init_max}{maximum value of dropout probability}

\item{is_mc_dropout}{boolean}

\item{supports_masking}{boolean}
}
\value{
custom layer
}
\description{
Layer wrapper with concrete dropout
}
