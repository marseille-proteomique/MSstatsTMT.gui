#' @title Volcano plot
#'
#' @description
#' \code{plo_volc} Plot a volcano from the output data from the function \code{\link{groupComparionTMT}}
#'
#'
#' @import ggplot2
#' @import gghighlight
#' @import stats
#' @import stringr
#' @importFrom matlab linspace
#'
#' @param data The output from \code{\link{groupComparionTMT}}
#' @param lim_pv The p-value threshold
#' @param lim_dif A double, corresponding to the two log2 fold-change thresholds
#' @param correction A logical, if TRUE will use the p-v.adjust column from the output from \code{\link{groupComparionTMT}}
#' @param perf_corr A logical, if TRUE will correct the p-values using the method from p.adjust
#' @param your_corr A character vector, if perf_corr TRUE will use this method for correct the p-values;
#'                  it will print the mehtod used for correcting the p-values on the plot
#' @param comp A character vector, indicating the comparison you want to see the volcano plot. Has to match one of the element of
#'             'Label' column from the data
#' @param curve A logical, if TRUE will plot a curvature on the volcano plot
#' @param curvature A numeric, indicating the value of the curvature
#' @param tit The title of the volcano plot
#' @param ytit The title of the y axis
#'
#' @return A volcano plot
#'
#' @seealso  \code{\link{groupComparionTMT}} from MSstatsTMT package for more details
#'
#' @export
#'
#' @examples
#' library(MSstatsTMT)
#' quant.pd.msstats <- proteinSummarization(input.pd,
#'                                          method="msstats",
#'                                          global_norm=TRUE,
#'                                          reference_norm=TRUE)
#'
#' test.pairwise <- groupComparisonTMT(quant.pd.msstats, moderated = TRUE)
#'
#' plo_volc(test.pairwise)

plo_volc <- function(data, lim_pv = 0.05, lim_dif = c(-1,1),
                     correction = FALSE, perf_corr = FALSE, your_corr = "none",
                     comp = "0.125-Norm", curve = FALSE, curvature = 0.1, tit = "mqpar25 MSstatsTMT",
                     ytit = "-log10(p-value)"){
  data_ <- data[,c(1,2,3,6,7)]
  if(!correction){
    data_ <- data_[,-ncol(data_)]
  }
  else{
    data_ <- data_[,-ncol(data_)+1]
  }
  data_ <- data_[which(data_$Label == comp), -2]
  n <- nrow(data_)
  data_ <- na.omit(data_)
  n <- n - nrow(data_)

  colnames(data_)[3] <- "pv"
  if(perf_corr){
    data_$pv <- p.adjust(data_$pv, your_corr)
  }
  data_$pv <- -log10(data_$pv)

  data_$reg <- rep("no diff", nrow(data_))

  data_$reg[which(data_$log2FC <= lim_dif[1])] <- "diff neg"
  data_$reg[which(data_$log2FC >= lim_dif[2])] <- "diff pos"

  g <- ggplot(data_, aes(log2FC, pv)) + geom_point(aes(color = reg)) +
    scale_color_manual(values = c("no diff" = "black",
                                  "diff neg" = "red",
                                  "diff pos" = "blue")) +
    geom_vline(xintercept = lim_dif[1], linetype = "dashed") +
    geom_vline(xintercept = lim_dif[2], linetype = "dashed") +
    geom_hline(yintercept = -log10(lim_pv), linetype = "dashed")

  if(curve){
    volc <- data.frame(x=linspace(-range(data_$log2FC)[2] -1, range(data_$log2FC)[2] +1, 500), y=1:500)
    volc$y[which(volc$x <= lim_dif[1])] <- curvature/abs(volc$x[which(volc$x <= lim_dif[1])] - lim_dif[1]) + -log10(lim_pv)
    volc$y[which(volc$x >= lim_dif[2])] <- curvature/abs(volc$x[which(volc$x >= lim_dif[2])] - lim_dif[2]) + -log10(lim_pv)
    volc$y[which(volc$x > lim_dif[1] & volc$x < lim_dif[2])] <- NA # curvature, lim_dif cutoff, lim_pv_list[[n]] cutoff


    data_$z <- data_$log2FC
    data_$z[which(data_$log2FC <= lim_dif[1])] <- curvature/abs(data_$log2FC[which(data_$log2FC <= lim_dif[1])] - lim_dif[1]) + -log10(lim_pv)
    data_$z[which(data_$log2FC >= lim_dif[2])] <- curvature/abs(data_$log2FC[which(data_$log2FC >= lim_dif[2])] - lim_dif[2]) + -log10(lim_pv)
    data_$z[which(data_$log2FC > lim_dif[1] & data_$log2FC < lim_dif[2])] <- NA

    g <- g + geom_line(aes(x,y), volc) +
      gghighlight(pv >= data_$z,
                  label_key = Protein,
                  use_direct_label = TRUE,
                  max_highlight = 50)

    subt <- paste("p-value threshold =", lim_pv, ", p-value correction :", your_corr,
                  "\nlog2 fold change threshold =", lim_dif[1], "to", lim_dif[2],
                  "\ncurvature of", curvature,
                  "\nremoved", n, "proteins with missing data")
  }

  else{
    g <- g + gghighlight(pv >= -log10(lim_pv) & reg != "no diff",
                         label_key = Protein,
                         use_direct_label = TRUE,
                         max_highlight = 50)

    subt <- paste("p-value threshold =", lim_pv, ", p-value correction =", your_corr,
                  "\n log2 fold change threshold =", lim_dif[1], "to", lim_dif[2],
                  "\nremoved", n, "proteins with missing data")
  }

  g <- g + xlim(-max(abs(data_$log2FC)) - 0.5, max(abs(data_$log2FC)) + 0.5) +
    ylim(0, max(data_$pv) + 1) +
    labs(title = tit, x = paste("log2FC", comp), y = ytit,
         subtitle = subt) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(size = 9))

  return(g)
}
