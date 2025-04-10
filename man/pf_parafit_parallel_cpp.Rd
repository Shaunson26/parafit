% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{pf_parafit_parallel_cpp}
\alias{pf_parafit_parallel_cpp}
\title{Run the parafit algorithm}
\usage{
pf_parafit_parallel_cpp(
  assoA,
  paraB,
  hostC,
  permutations,
  seed,
  verbose,
  print_n
)
}
\arguments{
\item{assoA}{association matrix}

\item{paraB}{parasite principal coordinates}

\item{hostC}{transposed host principal coordinates}

\item{permutations}{number of permutations}

\item{seed}{random seed value}

\item{verbose}{whether to print iteration number to the console}

\item{print_n}{the level to print iteration numbers by}
}
\value{
list of results, with elements global and links (the latter if test_links = TRUE).
}
\description{
Run the parafit algorithm to obtain the global test statistics and obtain a test
of significance using permutation. Optionally, obtain tests of individual associations.
}
