% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{groupThenSummarize}
\alias{groupThenSummarize}
\title{Customizable summary of data frame variables}
\usage{
groupThenSummarize(.data, var.group, .fns, ...)
}
\arguments{
\item{.data}{a data frame (extension)}

\item{var.group}{variables to group, either single name or gvars()}

\item{.fns}{named list of summary functions}

\item{...}{<\code{\link[dplyr:dplyr_tidy_select]{tidy-select}}> Selection of variables to summarize on}
}
\value{
a data frame (extension) of summary
}
\description{
Customizable summary of data frame variables
}
\examples{
# Single group variable
groupThenSummarize(mtcars, cyl, list(m=mean, s=sd), disp:wt)
# Multiple group variable
groupThenSummarize(mtcars, gvars(cyl,vs), list(m=mean, s=sd), disp:wt)
}
