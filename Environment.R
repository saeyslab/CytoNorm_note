### Packages ###################################################################
library(flowCore)
library(FlowSOM)
library(CytoNorm)
library(ggplot2) #Make plots
library(ggpubr) #Arrange plots
library(umap) #UMAP

### File organization ##########################################################
if (!exists("RDS")) dir.create("RDS")
if (!exists("Results")) dir.create("Results")

### Colors #####################################################################
batch_colors <- c("P1_C1" = "#900002",
                  "P1_C2" = "#FF3437",
                  "P2_C1" = "#0046A1",
                  "P2_C2" = "#6FCFFF",
                  "C1_goal" = "#7209b7",
                  "Goal distribution" = "#7209b7")
# MC_colors <- c("1" = "#9ef01a",
#                "2" = "#38b000",
#                "3" = "#006400")
# MC_colors <- c("1" = "#F8766D",
#                "2" = "#00BA38",
#                "3" = "#619CFF")
MC_colors <- c("1" = "#ffca3a",
               "2" = "#ff8600",
               "3" = "#8ac926")

### List all files per batch and define the cell type and state channels #######
files_P1_C1 <- list.files(path = "Data/Preprocessed/Panel1/Cohort1", 
                          pattern = "ID.*fcs", full.names = TRUE)
ff_P1 <- read.FCS(files_P1_C1[1])
files_P1_C2 <- list.files(path = "Data/Preprocessed/Panel1/Cohort2", 
                          pattern = "ID.*fcs", full.names = TRUE)

files_P2_C1 <- list.files(path = "Data/Preprocessed/Panel2/Cohort1", 
                          pattern = "ID.*fcs", full.names = TRUE)
ff_P2 <- read.FCS(files_P2_C1[1])

files_P2_C2 <- list.files(path = "Data/Preprocessed/Panel2/Cohort2", 
                          pattern = "ID.*fcs", full.names = TRUE)

cellTypeMarkers <- c("CD27", "CD38", "IgD", "IgM") 
cellTypeChannels <- GetChannels(object = ff_P1, markers = cellTypeMarkers)
P1_cellStateMarkers <- c("CD5", "CD268", "CD24", "PDL1", "PD1", "CD86", "KI67")
P1_cellStateChannels <- GetChannels(object = ff_P1, markers = P1_cellStateMarkers)
P2_cellStateMarkers <- c("CD40", "CD21", "HLA-DR", "ICOSL", "PD1", "CD86", "KI67")
P2_cellStateChannels <- GetChannels(object = ff_P2, markers = P2_cellStateMarkers)

### Functions ##################################################################
plotDensities <- function (input, channels, colors, model = NULL, transformList = NULL, 
                           show_goal = FALSE, suffix = c(original = "", normalized = "_norm")) 
{
  batch_names <- unique(gsub(suffix["original"], "", gsub(suffix["normalized"], 
                                                          "", names(input))))
  data <- list(original = list(), normalized = list())
  for (batch in batch_names) {
    for (type in c("original", "normalized")) {
      i <- paste0(batch, suffix[type])
      if (i %in% names(input)) {
        if (is.character(input[[i]])) {
          if (length(input[[i]] > 1)) {
            if (type == "original"){
              set.seed(2023)
              data[[type]][[batch]] <- FlowSOM::AggregateFlowFrames(fileNames = input[[i]], 
                                                                    cTotal = length(input[[i]]) * 10000)
            } else {
              data[[type]][[batch]] <- data[["original"]][[batch]]
              for (i_file in unique(flowCore::exprs(data[["original"]][[batch]])[,"File"])){
                ff <- read.FCS(input[[i]][i_file])
                i_IDs <- flowCore::exprs(data[["original"]][[batch]])[exprs(data[["original"]][[batch]])[,"File"] == i_file,"Original_ID2"] #Always Original_ID2???
                flowCore::exprs(data[[type]][[batch]])[flowCore::exprs(data[["original"]][[batch]])[,"File"] == i_file,1:(ncol(data[[type]][[batch]])-3)] <- flowCore::exprs(ff)[i_IDs, ]

                # i_IDs <- flowCore::exprs(data[["original"]][[batch]])[exprs(data[["original"]][[batch]])[,"File"] == i_file,"Original_ID2"] #Always Original_ID2???
                # flowCore::exprs(data[[type]][[batch]])[exprs(data[["original"]][[batch]])[,"File"] == i_file,1:(ncol(data[[type]][[batch]])-3)] <- flowCore::exprs(flowCore::read.FCS(input[[i]][i_file], which.lines = i_IDs))
                # 
              }
            }
          }
          else {
            data[[type]][[batch]] <- flowCore::read.FCS(input[[i]], 
                                                        truncate_max_range = FALSE)
          }
        }
        else if (methods::is(input[[i]], "flowFrame")) {
          data[[type]][[batch]] <- input[[i]]
        }
      }
      if (!is.null(transformList)) {
        data[[type]][[batch]] <- flowCore::transform(data[[type]][[batch]], 
                                                     transformList)
      }
    }
  }
  dfs <- list()
  for (type in c("original", "normalized")) {
    print(type)
    dfs[[type]] <- data.frame(do.call(rbind, lapply(data[[type]], 
                                                    function(x) flowCore::exprs(x))), check.names = FALSE)
    dfs[[type]]$Batch <- unlist(lapply(batch_names, function(i) rep(x = i, 
                                                                    times = nrow(data[[type]][[i]]))))
    if (!is.null(model) & type == "original") {
      mapped <- FlowSOM::NewData(model$fsom, as.matrix(dfs[[type]][model$fsom$map$colsUsed]))
      dfs[[type]][, "Cluster"] <- FlowSOM::GetMetaclusters(mapped)
      dfs[[type]] <- dfs[[type]][, c(channels, "File", 
                                     "Batch", "Cluster")]
    }
    else {
      dfs[[type]] <- dfs[[type]][, c(channels, "File", 
                                     "Batch")]
    }
  }
  if (!is.null(model)) {
    dfs[["normalized"]][, "Cluster"] <- dfs[["original"]][, 
                                                          "Cluster"]
  }
  if (show_goal) {
    if (is.null(model)) {
      stop("To show the goal, a model should be provided")
    }
    quantiles <- CytoNorm:::getCytoNormQuantiles(model)
  }
  plotlist <- list()
  for (channel in channels) {
    x_range <- c(min(sapply(dfs, function(df) min(df[[channel]], 
                                                  na.rm = TRUE))), max(sapply(dfs, function(df) max(df[[channel]], 
                                                                                                    na.rm = TRUE))))
    for (type in c("original", "normalized")) {
      df <- dfs[[type]]
      plotlist[[paste(channel, type)]] <- ggplot2::ggplot(df, ggplot2::aes(x = !!dplyr::sym(channel), 
                                            color = .data$Batch)) + ggplot2::stat_density(ggplot2::aes(group = paste(.data$Batch, 
                                                                                                                     .data$File)), geom = "line", position = "identity", 
                                                                                          alpha = 0.2) + ggplot2::stat_density(geom = "line", 
                                                                                                                               position = "identity", linewidth = 0.7) + ggplot2::scale_color_manual(values = colors) + 
        ggplot2::xlab(paste0(FlowSOM::GetMarkers(data[["original"]][[1]], 
                                                 channel), " <", channel, ">")) + ggplot2::ylab(type) + 
        ggplot2::theme_minimal() + ggplot2::xlim(x_range) + 
        ggplot2::ggtitle("")
      leg <- ggpubr::get_legend(plotlist[[paste(channel, type)]])
      plotlist[[paste(channel, type)]] <- plotlist[[paste(channel, type)]] +
        ggplot2::theme(legend.position = "none")
      if (type == "original" & !is.null(model)){
        plotlist[[paste(channel, type)]] <- plotlist[[paste(channel, type)]] +
          ggplot2::ggtitle(paste0("All cells"))
      }
      if (!is.null(model)) {
        for (cluster in levels(df$Cluster)){
          plotlist[[paste(channel, type, "cluster", cluster)]] <- ggplot2::ggplot(df[df$Cluster == cluster,], 
                                                                             ggplot2::aes(x = !!dplyr::sym(channel), color = .data$Batch))
          if (show_goal) {
            quantiles_df <- do.call(rbind, lapply(seq_along(quantiles), 
                                                  function(x) {
                                                    df_tmp <- data.frame(Value = 1/2 * (quantiles[[x]][-1, 
                                                                                                       channel] + quantiles[[x]][-nrow(quantiles[[x]]), 
                                                                                                                                 channel]), Density = 1/nrow(quantiles[[x]])/diff(quantiles[[x]][, 
                                                                                                                                                                                                 channel]), Cluster = x)
                                                    df_tmp$Density_smooth <- (stats::splinefun(df_tmp$Value, 
                                                                                               df_tmp$Density, method = "monoH.FC"))(df_tmp$Value)
                                                    df_tmp
                                                  }))
            quantiles_df$Batch <- "Goal distribution"
            plotlist[[paste(channel, type, "cluster", cluster)]] <- plotlist[[paste(channel, 
                                                                               type, "cluster", cluster)]] + ggplot2::geom_line(aes(x = Value, 
                                                                                                                               y = Density_smooth), data = quantiles_df[quantiles_df$Cluster == cluster,])
          }
          plotlist[[paste(channel, type, "cluster", cluster)]] <- plotlist[[paste(channel, 
                                                                             type, "cluster", cluster)]] + ggplot2::stat_density(ggplot2::aes(group = paste(.data$Batch, 
                                                                                                                                                       .data$File)), geom = "line", position = "identity", 
                                                                                                                            alpha = 0.2) + ggplot2::stat_density(geom = "line", 
                                                                                                                                                                 position = "identity", linewidth = 0.7) + 
            ggplot2::scale_color_manual(values = colors) + 
            ggplot2::xlab(paste0(FlowSOM::GetMarkers(data[["original"]][[1]], 
                                                     channel), " <", channel, ">")) + 
            ggplot2::ylab("") + 
            ggplot2::theme_minimal() + ggplot2::xlim(x_range) +
            ggplot2::ggtitle("")
          if (type == "original"){
            plotlist[[paste(channel, type, "cluster", cluster)]] <- plotlist[[paste(channel, type, "cluster", cluster)]] +
            ggplot2::ggtitle(paste0("Cluster ", cluster))
          }
         leg <- ggpubr::get_legend(plotlist[[paste(channel, type, "cluster", cluster)]]) 
          plotlist[[paste(channel, type, "cluster", cluster)]] <- plotlist[[paste(channel, type, "cluster", cluster)]] +
            ggplot2::theme(legend.position = "none")
               
 
        }
      }
    }
  }
  plotlist[["legend"]] <- leg
  return(plotlist)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
