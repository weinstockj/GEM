% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{fit_model}
\alias{fit_model}
\title{Fit a GEM model}
\usage{
fit_model(
  x,
  y,
  sample_map,
  lr = 0.1,
  iters = 180,
  num_threads = 4L,
  size = NULL,
  verbose = FALSE
)
}
\arguments{
\item{x}{A numeric matrix of predictors}

\item{y}{A numeric response vector}

\item{sample_map}{A list of sample indices}

\item{lr}{Learning rate}

\item{iters}{Number of iterations}

\item{num_threads}{Number of threads to set}

\item{size}{Optional subsample size}
}
\value{
A list containing model parameters and fitted values
}
\description{
Fit a GEM model
}
