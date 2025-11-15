#' Decision-curve analysis contour plot for survival outcomes
#'
#' @description
#' Produces a contour plot of net benefit from time-dependent decision-curve
#' analysis (DCA) for a Cox proportional hazards model. The x-axis represents
#' follow-up time and the y-axis represents the risk-threshold probability
#' used to define a positive prediction. At each time point the function
#' computes net benefit for the strategies "treat all", "treat none", and a
#' model-based strategy using \code{dcurves::dca()}. Optionally, a second
#' model can be supplied to visualize the difference in net benefit between
#' the two models.
#'
#' @param data A \code{data.frame} containing the survival outcome and all
#'   predictors used in \code{model} (and \code{model2}, if supplied).
#' @param time Numeric vector of follow-up times at which to evaluate net
#'   benefit and draw contours (x-axis).
#' @param threshold Numeric vector of risk-threshold probabilities at which to
#'   evaluate net benefit (y-axis).
#' @param model A fitted Cox proportional hazards model, typically a
#'   \code{\link[survival]{coxph}} object. This is the primary model whose
#'   net benefit surface is plotted.
#' @param option Character string specifying what to display in the contour
#'   plot. Allowed values are
#'   \itemize{
#'     \item \code{"Treat All"}: net benefit of the treat-all strategy.
#'     \item \code{"Treat None"}: net benefit of the treat-none strategy.
#'     \item \code{"Compare with All"}: difference in net benefit between the
#'       model-based strategy and the treat-all strategy.
#'     \item \code{"Compare two models"}: difference in net benefit between
#'       \code{model} and \code{model2}.
#'   }
#'   The default (\code{NULL}) is to plot the net benefit of the model-based
#'   strategy itself.
#' @param model2 Optional second \code{coxph} model to compare against
#'   \code{model} when \code{option = "Compare two models"}.
#'
#' @details
#' For each time point in \code{time}, the function first converts the model
#' predictions to time-dependent risks using
#' \code{\link[survival]{survfit}} and then calls
#' \code{\link[dcurves]{dca}} to obtain net benefit for the treat-all,
#' treat-none, and model-based strategies across the grid of threshold
#' probabilities. When there are not enough events to estimate survival
#' probability at a given time and threshold, the function replaces missing
#' net benefit values using a conservative backup calculation based on
#' true-positive and false-positive rates and reports this in the console.
#'
#' The resulting net benefit matrix is visualized with
#' \code{\link[plotly]{plot_ly}} as a contour plot, optionally showing
#' differences between strategies, which facilitates identifying clinically
#' meaningful time–threshold regions where a model provides incremental value.
#'
#' @return
#' A \code{\link[plotly]{plotly}} object showing a contour heatmap of net
#' benefit (or net benefit difference) over time and threshold.
#'
#' @seealso
#' \code{\link[dcurves]{dca}} for the underlying decision-curve analysis and
#' \code{\link[survival]{coxph}} for fitting Cox models.
#'
#' @examples
#' if (requireNamespace("dcurves", quietly = TRUE) &&
#'     requireNamespace("readr", quietly = TRUE) &&
#'     requireNamespace("labelled", quietly = TRUE) &&
#'     requireNamespace("dplyr", quietly = TRUE) &&
#'     requireNamespace("survival", quietly = TRUE) &&
#'     requireNamespace("plotly", quietly = TRUE)) {
#'
#'   library(dplyr)
#'   library(dcurves)
#'   library(survival)
#'   library(plotly)
#'
#'   df_time_to_cancer_dx <-
#'     readr::read_csv(
#'       file = paste0(
#'         "https://raw.githubusercontent.com/",
#'         "ddsjoberg/dca-tutorial/main/data/df_time_to_cancer_dx.csv"
#'       )
#'     ) %>%
#'     labelled::set_variable_labels(
#'       patientid        = "Patient ID",
#'       cancer           = "Cancer Diagnosis",
#'       ttcancer         = "Years to Diagnosis/Censor",
#'       risk_group       = "Risk Group",
#'       age              = "Patient Age",
#'       famhistory       = "Family History",
#'       marker           = "Marker",
#'       cancerpredmarker = "Prediction Model",
#'       cancer_cr        = "Cancer Diagnosis Status"
#'     )
#'
#'   coxmod <- coxph(
#'     Surv(ttcancer, cancer) ~ age + famhistory + marker,
#'     data = df_time_to_cancer_dx
#'   )
#'
#'   time_seq <- seq(
#'     min(df_time_to_cancer_dx$ttcancer),
#'     stats::quantile(df_time_to_cancer_dx$ttcancer, 0.9),
#'     length.out = 21
#'   )
#'   thres_seq <- seq(0, 0.5, 0.01)
#'
#'   # Net benefit surface of the model versus default strategies
#'   tempPlot <- dcaContour(
#'     data      = df_time_to_cancer_dx,
#'     time      = time_seq,
#'     threshold = thres_seq,
#'     model     = coxmod,
#'     option    = "Compare with All"
#'   )
#'
#'   # Difference in net benefit between two Cox models
#'   coxmod2 <- coxph(
#'     Surv(ttcancer, cancer) ~ age,
#'     data = df_time_to_cancer_dx
#'   )
#'
#'   tempPlot2 <- dcaContour(
#'     data      = df_time_to_cancer_dx,
#'     time      = time_seq,
#'     threshold = thres_seq,
#'     model     = coxmod,
#'     option    = "Compare two models",
#'     model2    = coxmod2
#'   )
#' }
#'
#' @import survival
#' @importFrom plotly plot_ly layout
#' @importFrom dcurves dca
#' @importFrom magrittr %>%
#'
#' @export
dcaContour <- function(data, time, threshold, model, option = NULL, model2 = NULL) {
  treatAllMat   <- matrix(NA, ncol = length(time), nrow = length(threshold))
  treatNoneMat  <- matrix(NA, ncol = length(time), nrow = length(threshold))
  treatModelMat <- matrix(NA, ncol = length(time), nrow = length(threshold))
  for (i in 1:length(time)) {
    data$tempVar <- 1 - summary(
      survfit(model, newdata = data),
      times = time[i]
    )$surv[1, ]

    eval(parse(text = paste0(
      "temp<-dca(",
      sub("(.*) ~.*", "\\1", as.character(model$formula))[2],
      "~ tempVar,data=data,
               time = time[i],
               thresholds = threshold)"
    )))
    nb <- temp$dca$net_benefit
    if (sum(is.na(nb[(1 + 2 * length(threshold)):(3 * length(threshold))])) > 0) {
      tp_backup <- rep(0, length(temp$dca$tp_rate))
      fp_backup <- temp$dca$pos_rate
      nb_backup <- tp_backup - fp_backup * rep(threshold / (1 - threshold), 3)
      nb_backup[nb_backup < 0] <- 0
      print(paste0(
        "At time ", round(time[i], 3),
        ", with threshold value ",
        threshold[is.na(nb[(1 + 2 * length(threshold)):(3 * length(threshold))])],
        ", there are not enough observations to calculate the survival ",
        "probability, and we impute NA value with ",
        nb_backup[is.na(nb[(1 + 2 * length(threshold)):(3 * length(threshold))])],
        "."
      ))
      nb[is.na(nb)] <- nb_backup[is.na(nb)]
      tp <- temp$dca$tp_rate
      tp[is.na(tp)] <- tp_backup[is.na(tp)]
      fp <- temp$dca$fp_rate
      fp[is.na(fp)] <- fp_backup[is.na(fp)]
    }
    treatAllMat[, i]   <- nb[1:length(threshold)]
    treatNoneMat[, i]  <- nb[(1 + length(threshold)):(2 * length(threshold))]
    treatModelMat[, i] <- nb[(1 + 2 * length(threshold)):(3 * length(threshold))]
  }

  if (!is.null(model2)) {
    treatModel2Mat <- matrix(NA, ncol = length(time), nrow = length(threshold))
    for (i in 1:length(time)) {
      data$tempVar <- 1 - summary(
        survfit(model2, newdata = data),
        times = time[i]
      )$surv[1, ]
      eval(parse(text = paste0(
        "temp<-dca(",
        sub("(.*) ~.*", "\\1", as.character(model$formula))[2],
        "~ tempVar,data=data,
               time = time[i],
               thresholds = threshold)"
      )))
      nb <- temp$dca$net_benefit
      if (sum(is.na(nb[(1 + 2 * length(threshold)):(3 * length(threshold))])) > 0) {
        tp_backup <- rep(0, length(temp$dca$tp_rate))
        fp_backup <- temp$dca$pos_rate
        nb_backup <- tp_backup - fp_backup * rep(threshold / (1 - threshold), 3)
        nb_backup[nb_backup < 0] <- 0
        print(paste0(
          "At time ", round(time[i], 3),
          ", with threshold value ",
          threshold[is.na(nb[(1 + 2 * length(threshold)):(3 * length(threshold))])],
          ", there are not enough observations to calculate the survival ",
          "probability, and we impute NA value with ",
          nb_backup[is.na(nb[(1 + 2 * length(threshold)):(3 * length(threshold))])],
          " for model 2."
        ))
        nb[is.na(nb)] <- nb_backup[is.na(nb)]
        tp <- temp$dca$tp_rate
        tp[is.na(tp)] <- tp_backup[is.na(tp)]
        fp <- temp$dca$fp_rate
        fp[is.na(fp)] <- fp_backup[is.na(fp)]
      }
      treatModel2Mat[, i] <- nb[(1 + 2 * length(threshold)):(3 * length(threshold))]
    }
  }

  outplot <- plot_ly(
    x = time, y = threshold, z = treatModelMat, type = "contour",
    hovertemplate = paste(
      "At time %{x:.2f} <br>with",
      "threshold %{y:.2f},<br>the net",
      "benefit is %{z:.2f}<extra></extra>"
    )
  )
  outplot <- outplot %>%
    layout(title = list(text = "Net benefit"))

  if (option == "Treat All") {
    outplot <- plot_ly(
      x = time, y = threshold, z = treatAllMat, type = "contour",
      hovertemplate = paste(
        "At time %{x:.2f} <br>with",
        "threshold %{y:.2f},<br>the net",
        "benefit is %{z:.2f}<extra></extra>"
      )
    )
    outplot <- outplot %>%
      layout(title = list(text = "Net benefit by treating all"))
  }
  if (option == "Treat None") {
    outplot <- plot_ly(
      x = time, y = threshold, z = treatNoneMat, type = "contour",
      hovertemplate = paste(
        "At time %{x:.2f} <br>with",
        "threshold %{y:.2f},<br>the net",
        "benefit is %{z:.2f}<extra></extra>"
      )
    )
    outplot <- outplot %>%
      layout(title = list(text = "Net benefit by treating none"))
  }
  if (option == "Compare with All") {
    outplot <- plot_ly(
      x = time, y = threshold, z = treatModelMat - treatAllMat, type = "contour",
      hovertemplate = paste(
        "At time %{x:.2f} <br>with",
        "threshold %{y:.2f},<br>the net",
        "benefit is %{z:.2f}<extra></extra>"
      )
    )
    outplot <- outplot %>%
      layout(title = list(text = "Net benefit compared with treating all"))
  }
  if (option == "Compare two models") {
    outplot <- plot_ly(
      x = time, y = threshold, z = treatModelMat - treatModel2Mat, type = "contour",
      hovertemplate = paste(
        "At time %{x:.2f} <br>with",
        "threshold %{y:.2f},<br>the difference in",
        "net benefit is %{z:.2f}<extra></extra>"
      )
    )
    outplot <- outplot %>%
      layout(title = list(text = "Net benefit difference of two models"))
  }

  outplot <- outplot %>%
    layout(
      xaxis = list(title = "Time"),
      yaxis = list(title = "Threshold")
    )
  outplot
}
