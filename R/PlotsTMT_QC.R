#' Visualization for explanatory data analysis - TMT experiment
#'
#' To illustrate the quantitative data and quality control of MS runs,
#' PlotsTMT_QC takes the quantitative data from converter functions (\code{\link{PDtoMSstatsTMTFormat}},
#' \code{\link{MaxQtoMSstatsTMTFormat}}, \code{\link{SpectroMinetoMSstatsTMTFormat}}) as input
#' and generate quality control plot, to evaluate the systematic bias between MS runs.
#' It is totally based on the function \code{\link{dataProcessPlotsTMT}}
#'
#' @import ggplot2
#' @importFrom graphics axis image legend mtext par plot.new title plot
#' @importFrom grDevices dev.off hcl pdf
#' @importFrom dplyr mutate
#' @importFrom reshape2 dcast
#'
#' @param data.peptide name of the data with peptide level, which can be the output of converter functions(\code{\link{PDtoMSstatsTMTFormat}}, \code{\link{MaxQtoMSstatsTMTFormat}}, \code{\link{SpectroMinetoMSstatsTMTFormat}}).
#' @param data.summarization name of the data with protein-level, which can be the output of \code{\link{proteinSummarization}} function.
#' @param ylimUp upper limit for y-axis in the log scale.
#' FALSE(Default) for Profile Plot and QC Plot uses the upper limit as rounded off maximum of log2(intensities) after normalization + 3..
#' @param ylimDown lower limit for y-axis in the log scale. FALSE(Default) for Profile Plot and QC Plot uses 0..
#' @param x.axis.size size of x-axis labeling for "Run" and "channel in Profile Plot and QC Plot.
#' @param y.axis.size size of y-axis labels. Default is 10.
#' @param text.size size of labels represented each condition at the top of Profile plot and QC plot. Default is 4.
#' @param text.angle angle of labels represented each condition at the top of Profile plot and QC plot. Default is 0.
#' @param legend.size size of legend above Profile plot. Default is 7.
#' @param which.Protein Protein list to draw plots. List can be names of Proteins or order numbers of Proteins.
#' Default is "all", which generates all plots for each protein. For QC plot, "allonly" will generate one QC plot with all proteins.
#'
#' @return QC plot
#'
#' @export
#'
#' @examples
#' PlotsTMT_QC(mqpar25_MSstat, quant_mqpar25,
#'             which.Protein = "Q13868")

PlotsTMT_QC <- function(data.peptide,
                        data.summarization,
                        ylimUp = FALSE,
                        ylimDown = FALSE,
                        x.axis.size = 10,
                        y.axis.size = 10,
                        text.size = 4,
                        text.angle = 90,
                        legend.size = 7,
                        which.Protein = "all"){
  Condition = Run = xorder = Channel = NULL
  PeptideSequence = PSM = NULL
  groupAxis = cumGroupAxis = abundance = analysis = NULL

  datafeature <- data.peptide
  datarun <- data.summarization

  # conditions in feature data
  fea.conds <- as.character(unique(datafeature$Condition))
  # conditions in protein data
  run.conds <- as.character(unique(datarun$Condition))

  # only keep the overlapped conditions between feature data and protein data
  shared.conds <- intersect(fea.conds, run.conds)
  datafeature <- datafeature[datafeature$Condition %in% shared.conds,]
  datarun <- datarun[datarun$Condition %in% shared.conds,]

  # make sure condition is factor
  datafeature$Condition <- factor(datafeature$Condition)
  datarun$Condition <- factor(datarun$Condition)

  colnames(datafeature)[colnames(datafeature) == 'ProteinName'] <- 'Protein'
  datafeature$Protein <- factor(datafeature$Protein)
  datarun$Protein <- factor(datarun$Protein)

  ## feature level data : log2 transform
  datafeature$abundance <- log2(datafeature$Intensity)
  datafeature[!is.na(datafeature$Intensity) &
                datafeature$Intensity < 1, 'abundance'] <- 0


  ## QC plot (Quality control plot) ##
  ## ---------------------------------

  ## y-axis labeling
  yaxis.name <- 'Log2-intensities'


  ## assign upper or lower limit
  y.limup <- ceiling(max(datafeature$abundance, na.rm = TRUE) + 3)

  if (is.numeric(ylimUp)) {
    y.limup <- ylimUp
  }

  y.limdown <- 0
  if (is.numeric(ylimDown)) {
    y.limdown <- ylimDown
  }

  datafeature <- datafeature[with(datafeature, order(Run, Condition, Channel)), ]

  ## !! important: order of x-axis
  ## can be reorder by group and then channel, WITHIN Run
  ## first make new column for x-axis
  datafeature$group.channel <- paste(datafeature$Condition, datafeature$Channel, sep = "_")

  ## not sure better way for coding
  ## potentially change it.
  datafeature$xorder <- NA

  for (k in seq_along(unique(datafeature$Run))) {

    runid <- unique(datafeature$Run)[k]
    datafeature[datafeature$Run == runid, ]$xorder <- factor(datafeature[datafeature$Run == runid, ]$group.channel,
                                                             levels <- unique(datafeature[datafeature$Run == runid, ]$group.channel),
                                                             labels <- seq(1, length(unique(datafeature[datafeature$Run == runid, ]$group.channel))))
  }

  ## check
  ## unique(datafeature[datafeature$Run == 'PAMI-176_Mouse_K-T', c('Channel', 'Condition', 'Run', 'xorder','group.channel')])

  ## need to make data.frame with same variables for condition name
  datafeature$xorder <- as.numeric(datafeature$xorder)
  ## keep unique information for x-axis labeling. will be used in plotting
  tempGroupName <- unique(datafeature[, c("Condition", "xorder", "Run", "Channel")])

  ## count # per condition per Run
  #groupline <- unique(datafeature[, c('Condition', 'Run')])
  #groupline$groupAxis <- as.numeric(xtabs(~Condition+Run, tempGroupName))
  groupline <- tempGroupName %>% dplyr::group_by(Condition, Run) %>% dplyr::mutate(groupAxis = n())
  groupline <- groupline %>% dplyr::select(-xorder, -Channel)
  groupline <- groupline[!duplicated(groupline), ]

  ## make accumurated # as condition increase
  groupline <- groupline %>% dplyr::group_by(Run) %>% dplyr::mutate(cumGroupAxis = cumsum(groupAxis))

  groupline$cumGroupAxis <- groupline$cumGroupAxis + 0.5

  ## add coordinate for group id
  groupline$xorder <- groupline$cumGroupAxis - groupline$groupAxis / 2
  groupline$abundance <- y.limup - 0.5

  ## save all information, for labeling group in plot
  groupline.all <- groupline

  ## remove last condition for vertical line between groups
  groupline <- groupline[-which(groupline$Condition %in% levels(groupline$Condition)[nlevels(groupline$Condition)]), ]

  ## all protein
  if (which.Protein == 'all' | which.Protein == 'allonly') {

    ## for annotation of condition
    groupline.tmp <- data.frame(groupline,
                                "PSM" = unique(datafeature$PSM)[1],
                                "PeptideSequence" = unique(datafeature$PeptideSequence)[1])

    groupline.all.tmp <- data.frame(groupline.all,
                                    "PSM" = unique(datafeature$PSM)[1],
                                    "PeptideSequence" = unique(datafeature$PeptideSequence)[1])

    ## 1st plot for original plot
    ## for boxplot, x-axis, xorder should be factor
    datafeature$xorder <- factor(datafeature$xorder)

    ptemp <- ggplot(aes_string(x = 'xorder', y = 'abundance'), data = datafeature) +
      facet_grid(~Run) +
      geom_boxplot(aes_string(fill = 'Condition'), outlier.shape = 1, outlier.size = 1.5) +
      labs(title = 'All proteins',
           x = 'MS runs') +
      scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
      geom_vline(data = groupline.tmp,
                 aes(xintercept = cumGroupAxis),
                 colour = "grey", linetype = "longdash") +
      geom_text(data = groupline.all.tmp,
                aes(x = xorder, y = abundance, label = Condition),
                size = text.size,
                angle = text.angle,
                color = "black") +
      theme(
        panel.background = element_rect(fill = 'white', colour = "black"),
        legend.key = element_rect(fill = 'white', colour = 'white'),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'gray95'),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = y.axis.size, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
        axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
        title = element_text(size = x.axis.size + 8, vjust = 1.5),
        legend.position = "none")


    message("Drew the Quality Contol plot(boxplot) for all proteins.")
  }

  ## each protein
  ## choose Proteins or not
  if (which.Protein != 'allonly') {
    if (which.Protein != "all") {
      ## check which.Protein is name of Protein
      if (is.character(which.Protein)) {

        temp.name <- which.Protein

        ## message if name of Protein is wrong.
        if (length(setdiff(temp.name, unique(datafeature$Protein))) > 0) {
          dev.off()
          stop(paste0("Please check protein name. Data set does not have this protein. - ",
                      toString(temp.name)))
        }
      }

      ## check which.Protein is order number of Protein
      if (is.numeric(which.Protein)) {
        temp.name <- levels(datafeature$Protein)[which.Protein]

        ## message if name of Protein is wrong.
        if (length(levels(datafeature$Protein)) < max(which.Protein)) {
          dev.off()
          stop(paste0("Please check your ion of proteins. There are ",
                      length(levels(datafeature$Protein)), " proteins in this dataset."))
        }
      }

      ## use only assigned proteins
      datafeature <- datafeature[which(datafeature$Protein %in% temp.name), ]
      datafeature$Protein <- factor(datafeature$Protein)
    }

    for (i in seq_len(nlevels(datafeature$Protein))) {
      sub <- datafeature[datafeature$Protein == levels(datafeature$Protein)[i], ]
      sub <- sub[!is.na(sub$abundance), ]

      ## if all measurements are NA,
      if (nrow(sub) == sum(is.na(sub$abundance))) {
        message(paste("Can't the Quality Control plot for ", unique(sub$Protein),
                      "(", i, " of ", length(unique(datafeature$Protein)),
                      ") because all measurements are NAs."))
        next()
      }

      ## for annotation of condition
      groupline.tmp <- data.frame(groupline,
                                  "PSM" = unique(sub$PSM)[1],
                                  "PeptideSequence" = unique(sub$PeptideSequence)[1])

      groupline.all.tmp <- data.frame(groupline.all,
                                      "PSM" = unique(sub$PSM)[1],
                                      "PeptideSequence" = unique(sub$PeptideSequence)[1])

      ## 1st plot for original plot
      ## for boxplot, x-axis, xorder should be factor
      sub$xorder <- factor(sub$xorder)

      ptemp <- ggplot(aes_string(x = 'xorder', y = 'abundance'), data = sub) +
        facet_grid(~Run) +
        geom_boxplot(aes_string(fill = 'Condition'), outlier.shape = 1, outlier.size = 1.5) +
        labs(title = unique(sub$Protein),
             x = 'MS runs') +
        scale_y_continuous(yaxis.name, limits = c(y.limdown, y.limup)) +
        geom_vline(data = groupline.tmp,
                   aes(xintercept = cumGroupAxis),
                   colour = "grey", linetype = "longdash") +
        geom_text(data = groupline.all.tmp,
                  aes(x = xorder, y = abundance, label = Condition),
                  size = text.size,
                  angle = text.angle,
                  color = "black") +
        theme(
          panel.background = element_rect(fill = 'white', colour = "black"),
          legend.key = element_rect(fill = 'white', colour = 'white'),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = 'gray95'),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = y.axis.size, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
          axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
          title = element_text(size = x.axis.size + 8, vjust = 1.5),
          legend.position = "none")



      message(paste("Drew the Quality Contol plot(boxplot) for ", unique(sub$Protein),
                    "(", i, " of ", length(unique(datafeature$Protein)), ")"))

    } # end-loop
  }

  return(ptemp)

  # end QC plot
}
