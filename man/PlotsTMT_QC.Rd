% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotsTMT_QC.R
\name{PlotsTMT_QC}
\alias{PlotsTMT_QC}
\title{Visualization for explanatory data analysis - TMT experiment}
\usage{
PlotsTMT_QC(
  data.peptide,
  data.summarization,
  ylimUp = FALSE,
  ylimDown = FALSE,
  x.axis.size = 10,
  y.axis.size = 10,
  text.size = 4,
  text.angle = 90,
  legend.size = 7,
  which.Protein = "all"
)
}
\arguments{
\item{data.peptide}{name of the data with peptide level, which can be the output of converter functions(\code{\link{PDtoMSstatsTMTFormat}}, \code{\link{MaxQtoMSstatsTMTFormat}}, \code{\link{SpectroMinetoMSstatsTMTFormat}}).}

\item{data.summarization}{name of the data with protein-level, which can be the output of \code{\link{proteinSummarization}} function.}

\item{ylimUp}{upper limit for y-axis in the log scale.
FALSE(Default) for Profile Plot and QC Plot uses the upper limit as rounded off maximum of log2(intensities) after normalization + 3..}

\item{ylimDown}{lower limit for y-axis in the log scale. FALSE(Default) for Profile Plot and QC Plot uses 0..}

\item{x.axis.size}{size of x-axis labeling for "Run" and "channel in Profile Plot and QC Plot.}

\item{y.axis.size}{size of y-axis labels. Default is 10.}

\item{text.size}{size of labels represented each condition at the top of Profile plot and QC plot. Default is 4.}

\item{text.angle}{angle of labels represented each condition at the top of Profile plot and QC plot. Default is 0.}

\item{legend.size}{size of legend above Profile plot. Default is 7.}

\item{which.Protein}{Protein list to draw plots. List can be names of Proteins or order numbers of Proteins.
Default is "all", which generates all plots for each protein. For QC plot, "allonly" will generate one QC plot with all proteins.}
}
\value{
QC plot
}
\description{
To illustrate the quantitative data and quality control of MS runs,
PlotsTMT_QC takes the quantitative data from converter functions (\code{\link{PDtoMSstatsTMTFormat}},
\code{\link{MaxQtoMSstatsTMTFormat}}, \code{\link{SpectroMinetoMSstatsTMTFormat}}) as input
and generate quality control plot, to evaluate the systematic bias between MS runs.
It is totally based on the function \code{\link{dataProcessPlotsTMT}}
}
\examples{
PlotsTMT_QC(mqpar25_MSstat, quant_mqpar25,
            which.Protein = "Q13868")
}
