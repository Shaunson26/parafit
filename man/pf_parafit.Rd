% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pf_parafit.R
\name{pf_parafit}
\alias{pf_parafit}
\title{Run the parafit algorithm}
\usage{
pf_parafit(
  host_pcoa,
  parasite_pcoa,
  associations,
  permutations,
  test_links = FALSE,
  seed,
  parallel = F,
  cores,
  verbose = FALSE,
  .print_n = 100,
  use_r = FALSE
)
}
\arguments{
\item{host_pcoa}{host principal coordinates (samples in rows, axes in columns)}

\item{parasite_pcoa}{parasite principal coordinates (samples in rows, axes in columns)}

\item{associations}{association matrix (host in rows, parasite in columns)}

\item{permutations}{number of permutations to conduct for significance testing}

\item{test_links}{test individual associations in associations matrix}

\item{seed}{seed for randomisation}

\item{parallel}{test links with parallel computing}

\item{cores}{number of cores if using parallel computing}

\item{verbose}{whether to print iteration numbers to the console}

\item{.print_n}{what iteration values to print by}

\item{use_r}{use R version of parafit (slower)}
}
\value{
list of results, with elements global and links (the latter if test_links = TRUE).
}
\description{
Run the parafit algorithm to obtain the global test statistics and obtain a test
of significance using permutation. Optionally, obtain tests of individual associations.
}
