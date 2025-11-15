# File: R/CIndexContour.R

#' Internal imports for C-index / AUC contour functions
#'
#' @name cindex_auc_internal_imports
#' @keywords internal
#'
#' @importFrom stats as.formula quantile model.matrix pnorm predict sd
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom survival Surv coxph survfit
#' @importFrom pec cindex
#' @importFrom riskRegression selectCox Score
#' @importFrom randomForestSRC rfsrc predict.rfsrc
#' @importFrom prodlim prodlim
#' @importFrom plotly plot_ly layout subplot
#' @importFrom xgboost xgb.DMatrix
NULL


#' Time-dependent AUC contour (dispatcher for Cox / RSF / XGB)
#'
#' @description
#' High-level wrapper that computes a time-dependent AUC surface over follow-up
#' time (x-axis) and thresholds of a continuous covariate (y-axis) using a
#' fitted Cox model, Random Survival Forest, or XGBoost AFT model.
#'
#' \itemize{
#'   \item \code{coxph}: uses \code{\link{coxAUC}}
#'   \item \code{rfsrc}: uses \code{\link{rsfAUC}}
#'   \item \code{xgb.Booster}: uses \code{\link{xgbAUC}}
#' }
#'
#' @param data A data.frame containing variables used in the fitted model,
#'   including \code{timeVar}, \code{statusVar} and \code{contCov}.
#' @param model A fitted survival model of class one of \code{coxph}, \code{rfsrc},
#'   or \code{xgb.Booster}.
#' @param contCov Character string; name of the continuous covariate used for
#'   constructing thresholds.
#' @param statusVar Character string; event indicator (0/1).
#' @param timeVar Character string; follow-up time variable.
#' @param contCovName Optional y-axis label for the contour plot. Defaults
#'   to \code{contCov}.
#' @param eval_times Optional numeric vector of evaluation time points. If
#'   \code{NULL}, suitable times are computed internally by the model-specific
#'   methods.
#' @param nCovEval Integer; number of thresholds between the 5--95\% quantiles
#'   of \code{contCov}.
#' @param nCores Integer; number of CPU cores used by the underlying methods.
#' @param drawHistogram Logical; if \code{TRUE}, a marginal histogram of
#'   \code{contCov} is displayed to the right of the contour plot.
#'
#' @return
#' A \pkg{plotly} object containing either a single contour heatmap or a
#' subplot of the contour and the marginal histogram.
#'
#' @seealso \code{\link{coxAUC}}, \code{\link{rsfAUC}}, \code{\link{xgbAUC}}
#'
#' @export
AUCContour <- function(
    data, model, contCov,
    statusVar, timeVar,
    contCovName = NULL,
    eval_times = NULL,
    nCovEval = 30, nCores = 2,
    drawHistogram = TRUE
) {
  if (is.null(contCovName)) contCovName <- contCov

  if (inherits(model, "coxph")) {
    plt <- coxAUC(
      data         = data,
      model        = model,
      contCov      = contCov,
      statusVar    = statusVar,
      timeVar      = timeVar,
      contCovName  = contCovName,
      eval_times   = eval_times,
      nCovEval     = nCovEval,
      nCores       = nCores,
      drawHistogram = drawHistogram
    )
  } else if (any(grepl("rfsrc", class(model)))) {
    plt <- rsfAUC(
      data         = data,
      model        = model,
      contCov      = contCov,
      statusVar    = statusVar,
      timeVar      = timeVar,
      contCovName  = contCovName,
      eval_times   = eval_times,
      nCovEval     = nCovEval,
      nCores       = nCores,
      drawHistogram = drawHistogram
    )
  } else if (inherits(model, "xgb.Booster")) {
    plt <- xgbAUC(
      data         = data,
      model        = model,
      contCov      = contCov,
      statusVar    = statusVar,
      timeVar      = timeVar,
      contCovName  = contCovName,
      eval_times   = eval_times,
      nCovEval     = nCovEval,
      nCores       = nCores,
      drawHistogram = drawHistogram
    )
  } else {
    stop("Unsupported model class: ", paste(class(model), collapse = ", "))
  }

  attr(plt, "type") <- "AUC"
  plt
}


#' Time-dependent AUC contour for Cox proportional hazards models
#'
#' @description
#' Computes a time-dependent IPCW-adjusted AUC across follow-up time (x-axis)
#' and a sequence of thresholds of a continuous covariate (y-axis) for a fitted
#' Cox proportional hazards model.
#'
#' @param data A data.frame containing at least \code{timeVar}, \code{statusVar},
#'   and \code{contCov}.
#' @param model A fitted \code{\link[survival]{coxph}} model.
#' @param contCov Character, name of the continuous covariate.
#' @param statusVar Character, event indicator (0/1).
#' @param timeVar Character, follow-up time variable.
#' @param contCovName Optional y-axis label. Defaults to `contCov`.
#' @param nCovEval Integer, number of thresholds between 5–95\% quantiles.
#' @param nCores Integer, CPU cores for parallelization.
#' @param eval_times Optional numeric vector of evaluation times. If missing,
#'   extracted from \code{survival::survfit(model)}.
#' @param drawHistogram Logical; if TRUE attach histogram of \code{contCov}.
#'
#' @details
#' For each threshold value \eqn{c}, a binary variable \eqn{I(\mathrm{contCov} > c)}
#' is created by anchoring \code{contCov} at its 5–95\% quantiles and a Cox model
#' is refit. IPCW-adjusted time-dependent AUC is then computed via \pkg{riskRegression}.
#' \pkg{riskRegression}.
#'
#' @return
#' A Plotly contour heatmap (with or without histogram), showing AUC as a
#' function of time and the \code{contCov} threshold.
#'
#' @examples
#' \dontrun{
#' vet <- survival::veteran
#' fit <- survival::coxph(Surv(time, status) ~ age, data = vet)
#' p <- coxAUC(vet, fit, contCov = "age",
#'                statusVar = "status", timeVar = "time",
#'                nCovEval = 20, nCores = 2, drawHistogram = TRUE)
#' }
#'
#' @seealso \code{\link{AUCContour}}, \code{\link{rsfAUC}}, \code{\link{xgbAUC}}
#'
#' @export
coxAUC <- function(
    data, model, contCov,
    statusVar, timeVar,
    contCovName = NULL,
    eval_times = NULL,
    nCovEval = 30, nCores = 2,
    drawHistogram = TRUE
){

  if (is.null(contCovName)) contCovName <- contCov

  # eval times
  if (is.null(eval_times)) {
    sf <- survival::survfit(model, newdata = data)
    eval_times <- sf$time
  }
  eval_times <- sort(unique(eval_times))

  # cutoffs & anchors
  cont_vals <- seq(
    stats::quantile(data[[contCov]], 0.05, na.rm = TRUE),
    stats::quantile(data[[contCov]], 0.95, na.rm = TRUE),
    length.out = nCovEval
  )
  j  <- match(contCov, names(data))
  lo <- as.numeric(stats::quantile(data[[j]], 0.05, na.rm = TRUE))
  hi <- as.numeric(stats::quantile(data[[j]], 0.95, na.rm = TRUE))

  runner <- function(cutoff){
    newdata <- data
    newdata[[j]] <- ifelse(data[[j]] > cutoff, hi, lo)
    lp <- as.numeric(stats::predict(model, newdata = newdata, type = "lp"))

    sc <- riskRegression::Score(
      object       = list(risk = lp),
      formula      = as.formula(paste0("Surv(", timeVar, ", ", statusVar, ") ~ 1")),
      data         = data,
      times        = eval_times,
      metrics      = "auc",
      summary      = "risk",
      cens.model   = "marginal",
      split.method = "none",
      ROC          = FALSE
    )

    tab <- as.data.frame(sc$AUC$score)[, c("times","AUC")]
    tab <- tab[order(tab$times), , drop = FALSE]
    mt  <- match(eval_times, tab$times)
    out <- rep(NA_real_, length(eval_times))
    out[!is.na(mt)] <- tab$AUC[mt[!is.na(mt)]]
    out
  }

  # parallel
  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  Z_list <- foreach::foreach(
    cutoff = cont_vals,
    .packages = c("survival","riskRegression","prodlim")
  ) %dopar% runner(cutoff)

  Z <- do.call(rbind, Z_list)
  rownames(Z) <- round(cont_vals, 3)
  colnames(Z) <- paste0("t=", round(eval_times, 3))

  plt <- plotly::plot_ly(
    x = eval_times, y = cont_vals, z = Z, type = "contour",
    contours = list(coloring = "heatmap"),
    colorbar = list(title = "AUC"),
    hovertemplate = paste0(
      "Time=%{x:.3f}<br>", contCovName,
      " cutoff=%{y:.3f}<br>AUC=%{z:.3f}<extra></extra>"
    )
  ) |>
    plotly::layout(
      title = "Cox time-dependent AUC",
      xaxis = list(title = "Time"),
      yaxis = list(title = contCovName)
    )

  if (isTRUE(drawHistogram)) {
    histData <- data[[contCov]]
    hist_plot <- plotly::plot_ly(
      y = histData, type = "histogram", hoverinfo = "none"
    ) |>
      plotly::layout(
        xaxis = list(title = "Count"),
        yaxis = list(title = contCovName)
      )

    s <- plotly::subplot(
      plt, hist_plot,
      nrows = 1, widths = c(0.8, 0.2), margin = 0.01,
      shareY = TRUE, titleX = TRUE, titleY = TRUE
    ) |>
      plotly::layout(margin = list(t = 50, l = 5))

    return(s)
  } else {
    plt |>
      plotly::layout(margin = list(t = 50, l = 5))
  }
}


#' Time-dependent AUC contour for Random Survival Forests
#'
#' @description
#' Computes a matrix of time-dependent AUC values for a fitted
#' \pkg{randomForestSRC} survival forest by perturbing a continuous covariate
#' around a sequence of thresholds and evaluating IPCW-adjusted AUC via \pkg{riskRegression}.
#'
#' @param data A data.frame containing \code{timeVar}, \code{statusVar},
#'   and \code{contCov}.
#' @param model A fitted `rfsrc` survival forest.
#' @param contCov Character, continuous covariate name.
#' @param timeVar Character, follow-up time variable.
#' @param statusVar Character, event indicator (0/1).
#' @param contCovName Optional y-axis label. Defaults to `contCov`.
#' @param nCovEval Number of thresholds (5–95\% quantiles).
#' @param nCores CPU cores for parallel computing.
#' @param eval_times Optional vector of times to evaluate C-index.
#' @param drawHistogram Logical: show histogram of \code{contCov}.
#'
#' @return
#' A Plotly contour heatmap (with or without histogram), showing AUC as a
#' function of time and the \code{contCov} threshold.
#'
#' @examples
#' \dontrun{
#' vet <- survival::veteran
#' fit <- randomForestSRC::rfsrc(Surv(time, status) ~ age + karno, data = vet)
#' p <- rsfAUC(vet, fit, contCov = "age",
#'                statusVar = "status", timeVar = "time",
#'                nCovEval = 20, nCores = 2, drawHistogram = TRUE)
#' }
#'
#' @seealso \code{\link{AUCContour}}, \code{\link{coxAUC}}, \code{\link{xgbAUC}}
#'
#' @export
rsfAUC <- function(
    data, model, contCov,
    statusVar, timeVar,
    contCovName = NULL,
    eval_times = NULL,
    nCovEval = 30, nCores = 2,
    drawHistogram = TRUE
){

  if (is.null(contCovName)) contCovName <- contCov

  if (is.null(eval_times)) {
    tt <- data[[timeVar]]
    eval_times <- as.numeric(stats::quantile(
      tt[tt > 0], probs = seq(0.05, 0.95, length.out = 40), na.rm = TRUE
    ))
  }
  eval_times <- sort(unique(eval_times))

  cont_vals <- seq(
    stats::quantile(data[[contCov]], 0.05, na.rm = TRUE),
    stats::quantile(data[[contCov]], 0.95, na.rm = TRUE),
    length.out = nCovEval
  )
  sdj <- stats::sd(data[[contCov]], na.rm = TRUE)
  eps <- ifelse(is.finite(sdj) && sdj > 0, 1e-6 * sdj, 1e-6)

  j <- match(contCov, names(data))

  runner <- function(cutoff){
    newdata <- data
    newdata[[j]] <- ifelse(data[[j]] > cutoff, cutoff + eps, cutoff - eps)

    pr <- randomForestSRC::predict.rfsrc(
      model, newdata = newdata, na.action = "na.impute"
    )
    risk <- pr$chf[, ncol(pr$chf)]

    sc <- riskRegression::Score(
      object       = list(rsf = risk),
      formula      = as.formula(paste0("Surv(", timeVar, ", ", statusVar, ") ~ 1")),
      data         = data,
      times        = eval_times,
      metrics      = "auc",
      summary      = "risk",
      cens.model   = "marginal",
      split.method = "none",
      ROC          = TRUE
    )

    tab <- as.data.frame(sc$AUC$score)[, c("times","AUC")]
    tab <- tab[order(tab$times), , drop = FALSE]
    mt  <- match(eval_times, tab$times)
    out <- rep(NA_real_, length(eval_times))
    out[!is.na(mt)] <- tab$AUC[mt[!is.na(mt)]]
    out
  }

  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  Z_list <- foreach::foreach(
    cutoff = cont_vals,
    .packages = c("randomForestSRC","riskRegression","survival")
  ) %dopar% runner(cutoff)

  Z <- do.call(rbind, Z_list)
  rownames(Z) <- round(cont_vals, 3)
  colnames(Z) <- paste0("t=", round(eval_times, 2))

  plt <- plotly::plot_ly(
    x = eval_times, y = cont_vals, z = Z, type = "contour",
    contours = list(coloring = "heatmap"),
    colorbar = list(title = "AUC"),
    hovertemplate = paste0(
      "time=%{x:.2f}<br>", contCovName,
      " cutoff=%{y:.2f}<br>AUC=%{z:.3f}<extra></extra>"
    )
  ) |>
    plotly::layout(
      xaxis = list(title = "Time"),
      yaxis = list(title = contCovName),
      title = "RSF time-dependent AUC"
    )

  if (isTRUE(drawHistogram)) {
    histData <- data[[contCov]]
    hist_plot <- plotly::plot_ly(
      y = histData, type = "histogram", hoverinfo = "none"
    ) |>
      plotly::layout(
        xaxis = list(title = "Count"),
        yaxis = list(title = contCovName)
      )

    s <- plotly::subplot(
      plt, hist_plot,
      nrows = 1, widths = c(0.8, 0.2), margin = 0.01,
      shareY = TRUE, titleX = TRUE, titleY = TRUE
    ) |>
      plotly::layout(margin = list(t = 50, l = 5))

    return(s)
  } else {
    plt |>
      plotly::layout(margin = list(t = 50, l = 5))
  }
}


#' Time-dependent AUC contour for XGBoost AFT survival models
#'
#' @description
#' Computes IPCW-adjusted time-dependent AUC using \pkg{riskRegression} for a
#' fitted \pkg{xgboost} Accelerated Failure Time (AFT) survival model across
#' evaluation time and continuous covariate thresholds. Produces a Plotly
#' contour heatmap, optionally with a marginal histogram of \code{contCov}.
#'
#' @param data A data.frame containing \code{timeVar}, \code{statusVar},
#'   and \code{contCov}.
#' @param model A fitted `xgb.Booster` with AFT objective.
#' @param contCov Character, continuous covariate name.
#' @param statusVar Character, event indicator (0/1).
#' @param timeVar Character, follow-up time variable.
#' @param contCovName Optional y-axis label. Defaults to `contCov`.
#' @param nCovEval Number of thresholds evaluated.
#' @param eval_times Optional numeric time vector.
#' @param nCores Integer, CPU cores for parallelization.
#' @param drawHistogram Attach histogram of \code{contCov} if TRUE.
#'
#' @return
#' A Plotly contour heatmap (with or without histogram), showing AUC as a
#' function of time and the \code{contCov} threshold.
#'
#' @examples
#' \dontrun{
#' # xgb_fit <- xgboost::xgb.train(params = list(objective = "survival:aft"), ...)
#' p <- xgbAUC(vet, xgb_fit, contCov = "age",
#'                statusVar = "status", timeVar = "time",
#'                nCovEval = 20, nCores = 2, drawHistogram = TRUE)
#' }
#'
#' @seealso \code{\link{AUCContour}}, \code{\link{coxAUC}}, \code{\link{rsfAUC}}
#'
#' @export
xgbAUC <- function(
    data, model, contCov,
    statusVar, timeVar,
    contCovName = NULL,
    eval_times = NULL,
    nCovEval = 30, nCores = 2,
    drawHistogram = TRUE
){

  if (is.null(contCovName)) contCovName <- contCov

  if (is.null(eval_times)) {
    tt <- data[[timeVar]]
    eval_times <- as.numeric(stats::quantile(
      tt[tt > 0], probs = seq(0.05, 0.95, length.out = 40), na.rm = TRUE
    ))
  }
  eval_times <- sort(unique(eval_times))

  cont_vals <- seq(
    stats::quantile(data[[contCov]], 0.05, na.rm = TRUE),
    stats::quantile(data[[contCov]], 0.95, na.rm = TRUE),
    length.out = nCovEval
  )
  sdj <- stats::sd(data[[contCov]], na.rm = TRUE)
  eps <- ifelse(is.finite(sdj) && sdj > 0, 1e-6 * sdj, 1e-6)

  feat <- model$feature_names
  build_X <- function(df){
    X <- stats::model.matrix(
      ~ . - 1,
      data = df[, setdiff(names(df), c(statusVar, timeVar)), drop = FALSE]
    )
    if (!is.null(feat)) {
      miss <- setdiff(feat, colnames(X))
      if (length(miss)) {
        X <- cbind(
          X,
          matrix(0, nrow(X), length(miss), dimnames = list(NULL, miss))
        )
      }
      X <- X[, feat, drop = FALSE]
    }
    X
  }

  X0 <- build_X(data)
  j_col <- match(contCov, colnames(X0))

  runner <- function(cutoff){
    X_cut <- X0
    X_cut[, j_col] <- ifelse(
      data[[contCov]] > cutoff, cutoff + eps, cutoff - eps
    )
    mu <- as.numeric(predict(
      model, newdata = xgboost::xgb.DMatrix(X_cut), outputmargin = TRUE
    ))
    risk <- -mu

    sc <- riskRegression::Score(
      object       = list(xgb = risk),
      formula      = as.formula(paste0("Surv(", timeVar, ", ", statusVar, ") ~ 1")),
      data         = data,
      times        = eval_times,
      metrics      = "auc",
      summary      = "risk",
      cens.model   = "marginal",
      split.method = "none",
      ROC          = TRUE
    )

    tab <- as.data.frame(sc$AUC$score)[, c("times","AUC")]
    tab <- tab[order(tab$times), , drop = FALSE]
    mt  <- match(eval_times, tab$times)
    out <- rep(NA_real_, length(eval_times))
    out[!is.na(mt)] <- tab$AUC[mt[!is.na(mt)]]
    out
  }

  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  Z_list <- foreach::foreach(
    cutoff = cont_vals,
    .packages = c("xgboost","riskRegression","survival")
  ) %dopar% runner(cutoff)

  Z <- do.call(rbind, Z_list)
  rownames(Z) <- round(cont_vals, 3)
  colnames(Z) <- paste0("t=", round(eval_times, 2))

  plt <- plotly::plot_ly(
    x = eval_times, y = cont_vals, z = Z, type = "contour",
    contours = list(coloring = "heatmap"),
    colorbar  = list(title = "AUC"),
    hovertemplate = paste0(
      "time=%{x:.2f}<br>", contCovName,
      " cutoff=%{y:.2f}<br>AUC=%{z:.3f}<extra></extra>"
    )
  ) |>
    plotly::layout(
      xaxis = list(title = "Time"),
      yaxis = list(title = contCovName),
      title = "XGB-AFT time-dependent AUC"
    )

  if (isTRUE(drawHistogram)) {
    histData <- data[[contCov]]
    hist_plot <- plotly::plot_ly(
      y = histData, type = "histogram", hoverinfo = "none"
    ) |>
      plotly::layout(
        xaxis = list(title = "Count"),
        yaxis = list(title = contCovName)
      )

    s <- plotly::subplot(
      plt, hist_plot,
      nrows = 1, widths = c(0.8, 0.2), margin = 0.01,
      shareY = TRUE, titleX = TRUE, titleY = TRUE
    ) |>
      plotly::layout(margin = list(t = 50, l = 5))

    return(s)
  } else {
    plt |>
      plotly::layout(margin = list(t = 50, l = 5))
  }
}
