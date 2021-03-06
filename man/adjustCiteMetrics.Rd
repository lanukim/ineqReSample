% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_RPackage.R
\name{adjustCiteMetrics}
\alias{adjustCiteMetrics}
\title{adjustCiteMetrics}
\usage{
adjustCiteMetrics(
  papers,
  cites,
  pubID = "paperID",
  pubYear = "publishedYear",
  citedID = "citedPaperID",
  citedYear = "citedPaperYear",
  citingYear = NULL,
  citationWindow = NULL,
  quantiles = c(0.2, 0.8),
  refYear = NULL,
  refPaperCount = NULL,
  refCiteCount = NULL,
  sims = 10,
  lowCitesInclusion = FALSE,
  lowCitesThreshold = 0.1,
  paperDupThreshold = 0.95,
  periods = NULL,
  verbose = FALSE
)
}
\arguments{
\item{papers}{required dataframe containing all citable papers studied (one row per paper). This dataframe must include two columns: the published paper IDs (\code{paperID}) and the publication year of the paper (\code{publishedYear}).}

\item{cites}{required dataframe containing all citations made to the papers listed in \code{papers} over a given set of years (one row per citation). This dataframe must include at least two columns: the IDs of the cited papers (\code{citedPaperID}), corresponding to the IDs given in \code{papers}, and the cited paper's publication year (\code{citedPaperYear}). Optionally, this dataframe may include the published year of paper sending the citation (\code{citingYear}) may be included optional.  \code{citingYear} must be included if users wish to filter citations made within \code{citationWindow} years of the cited papers' publication dates.}

\item{pubID}{character string, the column of \code{papers} containing publication IDs (default is `paperID').  IDs may be numeric, character, or factor.}

\item{pubYear}{character string, the column of \code{papers} containing publication years of each paper in numeric format (default is `publishedYear')}

\item{citedID}{character string, the column of \code{cites} containing the IDs of cited papers (default is `citedPaperID'). IDs may be numeric, character, or factor.}

\item{citedYear}{character string, the column of \code{cites} containing the years (as numeric data) in which papers were cited (default is `citedPaperYear')}

\item{citingYear}{character string, the column of \code{cites} containing the year (as numeric data) in which each citing paper was published (default = NULL)}

\item{citationWindow}{numeric scalar, the number of years covered by the citation window. For example, 2 means the function only analyzes citations made within 2 years after the cited paper is published. (default = NULL)}

\item{quantiles}{numeric vector in the interval \(0,1\] indicating the citation quantiles to be computed. Default = c(.20, .80), for the percentage of papers accounting for 20 percent and 80 percent of all citations.}

\item{refYear}{numeric scalar, the reference (or baseline) year (or other indicated period) among available publication years of papers in \code{papers}. Either \code{refYear} or \code{refPaperCount}/\code{refCiteCount} must be provided. (default = NULL)}

\item{refPaperCount}{a numeric scalar identifying the referenced publication count. Either \code{refYear} or \code{refPaperCount}/\code{refCiteCount} must be provided. (default = NULL)}

\item{refCiteCount}{a numeric scalar identifying the referenced citation count. Either \code{refYear} or \code{refPaperCount}/\code{refCiteCount} must be provided. (default = NULL)}

\item{sims}{optional; a numeric scalar identifying the number of times to run the simulation (default = 10)}

\item{lowCitesInclusion}{logical, TRUE: reports all results; FALSE: excludes any years from the results where the total number of citations is lower than the reference year in at least \code{lowCitesThreshold} percent of simulations (default=FALSE)}

\item{lowCitesThreshold}{numeric scalar in the interval \[0,1\] identifying threshold to define lowCites (default = 0.1)}

\item{paperDupThreshold}{numeric scalar in the interval \(0,1\] (default = 0.95)}

\item{periods}{optional; numeric vectors identifying the years of data to adjust. If NULL (default), the function automatically selects all years in \code{papers}.}

\item{verbose}{logical; returns detailed warnings related to lowCites (default=FALSE)}
}
\value{

}
\description{
Correct citation inequality measures for marginals bias via resampling
}
\details{
Kim, Adolph, West, and Stovel (2020) show that measures of inequality are subject to ``marginals bias,'' which impairs the ability to compare measures across time or contexts.  Consider the example of citation distributions, which describe the relative shares of incoming citations sent to each member of a population of published papers, e.g., for a specific field and year.  Inequality measures like the Gini coefficient (and even typically robust measures such as quantiles of the distribution) are only comparable across disciplinary fields or time periods if both the total number of papers potentially receiving citations and the total number of citations sent to this population remain constant.  If these ``marginals'' are rising, inequality may be understated, and if they are falling, inequality may be overstated.

This function corrects inequality measures for marginals bias by resampling the observed citations to have common marginals fixed to some reference level (the method suggested by Kim, Adolph, West, and Stovel, 2020).  Users may request either that inequality measures be adjusted to be comparable to those of a base period (via \code{refYear}), or to a fixed, user-provided count of papers and citations (vias \code{refPaperCount} and \code{refCiteCount}).  (The inputs to the function are named to correspond to the marginals in the citation example, but should have more general utility.)

Available corrected inequality measures include the percentage of ever cited papers, the Gini coefficient, the Herfindahl-Hirschman index, and (possibly user-defined) quantile measures. Results can be most easily accessed by applying the extractor function \code{citeIneq} to the object returned by this function.

Because results can vary across re-sampling iterations, the function averages adjusted measures across \code{sims} iterations. By default, this is 10 runs, but if adjusted values are unstable when re-running \code{adjustCiteMetrics}, users should increase \code{sims} until the results are consistent to the desired precision.  This is more likely to be necessary if the total number of reference papers or citations is small.

Note the resampling correction is only feasible for a target population if the number of citations sent to the target population is at least a bit larger than the reference number of total citations.  Otherwise, sampling without replacement will exhaust the pool of citations before a sufficient number of citations have been sampled to match the reference count.  When the total available citations are only slightly larger than the reference level, this can happen by chance in resampling, though typically the shortage of citations is quite small and may only affect a small fraction of iterations.  By default, \code{adjustCiteMetrics} will suppress results for any period in which the proportion of runs with a shortage of citations is greater than \code{lowCitesThreshold}.  If this poses a problem, consider using a smaller reference value of total citations.
}
\examples{
data(papers)
data(cites)
# May takes several minutes to run
result <- adjustCiteMetrics(papers = papers, 
                            cites = cites, 
                            pubID = "paperID", 
                            pubYear = "publishedYear", 
                            citedID = "citedID",
                            citedYear = "citedYear",
                            refYear = 1996)
citeIneq(result)
}
\references{
Lanu Kim, Christopher Adolph, Jevin West, and Katherine Stovel. 2020. ``The Influence of Changing Marginals on Measures of Inequality in Scholarly Citations: Evidence of Bias and a Resampling Correction.'' forthcoming in Sociological Science.

\url{http://faculty.washington.edu/cadolph/articles/kawsCitations.pdf}
}
\seealso{
\code{\link{citeIneq}}
}
