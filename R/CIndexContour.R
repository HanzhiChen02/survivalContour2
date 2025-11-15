# File: R/CIndexContour.R

#' Internal imports for C-index contour functions
#'
#' @name cindex_internal_imports
#' @keywords internal
#' @importFrom stats as.formula quantile model.matrix pnorm predict
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom survival Surv coxph
#' @importFrom pec cindex
#' @importFrom riskRegression selectCox
#' @importFrom randomForestSRC rfsrc
#' @importFrom prodlim prodlim
#' @importFrom plotly plot_ly layout subplot
#' @importFrom xgboost xgb.DMatrix
NULL


#' Time-dependent C-index contour (dispatcher for Cox / RSF / XGB)
#'
#' @description
#' Generic wrapper for computing a time-dependent C-index surface over
#' follow-up time (x-axis) and a continuous covariate threshold (y-axis).
#' It automatically dispatches to the appropriate backend based on the
#' class of \code{model}:
#'
#' \itemize{
#'   \item \code{coxph}: uses \code{\link{coxCIndex}}
#'   \item \code{rfsrc}: uses \code{\link{rsfCIndex}}
#'   \item \code{xgb.Booster}: uses \code{\link{xgbCIndex}}
#' }
#'
#' @param data A \code{data.frame} containing variables used in the fitted model.
#' @param model A fitted survival model: one of \code{coxph}, \code{rfsrc},
#'   or \code{xgb.Booster}.
#' @param contCov Character string; name of the continuous covariate used
#'   to define thresholds on the y-axis.
#' @param statusVar Character string; name of the event indicator (0/1).
#' @param timeVar Character string; name of the follow-up time variable.
#' @param contCovName Optional label for the y-axis. Defaults to \code{contCov}.
#' @param nCovEval Integer; number of thresholds between the 5\% and 95\%
#'   quantiles of \code{contCov}.
#' @param nCores Integer; number of CPU cores passed to the underlying
#'   C-index backends (for parallelization, where applicable).
#' @param eval_times Optional numeric vector of evaluation time points on
#'   the x-axis. If \code{NULL}, reasonable defaults are chosen within the
#'   backend functions.
#' @param drawHistogram Logical; if \code{TRUE}, a histogram of \code{contCov}
#'   is attached to the right of the contour plot.
#'
#' @return
#' A Plotly figure. Depending on \code{drawHistogram}, this is either
#' a single contour heatmap or a side-by-side subplot combining the
#' contour with a histogram of the continuous covariate.
#'
#' @seealso \code{\link{coxCIndex}}, \code{\link{rsfCIndex}}, \code{\link{xgbCIndex}}
#'
#' @export
CIndexContour <- function(
    data, model, contCov,
    statusVar = NULL, timeVar = NULL,
    contCovName = NULL,
    nCovEval = 30, nCores = 4,
    eval_times = NULL,
    drawHistogram = TRUE
) {
  if (is.null(contCovName)) contCovName <- contCov

  if (inherits(model, "coxph")) {
    return(
      coxCIndex(
        data = data, model = model, contCov = contCov,
        statusVar = statusVar, timeVar = timeVar,
        contCovName = contCovName, nCovEval = nCovEval,
        nCores = nCores, eval_times = eval_times,
        drawHistogram = drawHistogram
      )
    )
  } else if (any(grepl("rfsrc", class(model)))) {
    return(
      rsfCIndex(
        data = data, model = model, contCov = contCov,
        timeVar = timeVar, statusVar = statusVar,
        contCovName = contCovName, nCovEval = nCovEval,
        nCores = nCores, eval_times = eval_times,
        eps = 1e-6, drawHistogram = drawHistogram
      )
    )
  } else if (inherits(model, "xgb.Booster")) {
    return(
      xgbCIndex(
        data = data, model = model, contCov = contCov,
        statusVar = statusVar, timeVar = timeVar,
        contCovName = contCovName, nCovEval = nCovEval,
        eval_times = eval_times, sigma = 1.6, eps = 1e-6,
        drawHistogram = drawHistogram
      )
    )
  } else {
    stop("Unsupported model class: ", paste(class(model), collapse = ", "))
  }
}

#' Time-dependent C-index contour for Cox proportional hazards models
#'
#' @description
#' Computes a time-dependent IPCW-adjusted C-index across follow-up time (x-axis)
#' and a sequence of thresholds of a continuous covariate (y-axis).
#' Returns a Plotly contour plot, optionally with an attached histogram.
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
#' For each threshold value \eqn{c}, a binary variable \eqn{I(contCov > c)} is created
#' and a Cox model is refit. IPCW-adjusted C-index is computed via \pkg{pec} /
#' \pkg{riskRegression}.
#'
#' @return A Plotly contour heatmap (with or without histogram).
#'
#' @examples
#' \dontrun{
#' vet <- survival::veteran
#' fit <- survival::coxph(Surv(time, status) ~ age, data = vet)
#' p <- coxCIndex(vet, fit, contCov = "age",
#'                statusVar = "status", timeVar = "time")
#' }
#'
#' @seealso \code{\link{CIndexContour}}, \code{\link{rsfCIndex}}, \code{\link{xgbCIndex}}
#'
#' @export
coxCIndex <- function(
    data, model, contCov,
    statusVar = NULL, timeVar = NULL,
    contCovName = NULL,
    nCovEval = 30, nCores = 4,
    eval_times = NULL,
    drawHistogram = TRUE
) {
  if (is.null(timeVar) || is.null(statusVar)) {
    stop("For Cox C-index, 'timeVar' and 'statusVar' must be provided.")
  }
  if (is.null(contCovName)) contCovName <- contCov

  # parallel backend
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)

  # evaluation times
  if (is.null(eval_times)) {
    s0 <- survival::survfit(model, newdata = data)
    eval_times <- s0$time
  }
  eval_times <- sort(unique(eval_times))

  # thresholds
  cont_vals <- seq(
    quantile(data[[contCov]], 0.05, na.rm = TRUE),
    quantile(data[[contCov]], 0.95, na.rm = TRUE),
    length.out = nCovEval
  )

  # each threshold → one row of C-index over time
  Z <- foreach(i = seq_along(cont_vals), .combine = rbind) %dopar% {
    cutoff <- cont_vals[i]
    di <- data
    di$cont_binary <- ifelse(di[[contCov]] > cutoff, 1, 0)

    mod_bin <- survival::coxph(
      as.formula(paste0("Surv(", timeVar, ", ", statusVar, ") ~ cont_binary")),
      data = di, x = TRUE, y = TRUE
    )

    ci <- pec::cindex(
      object      = mod_bin,
      formula     = as.formula(paste0("Surv(", timeVar, ", ", statusVar, ") ~ cont_binary")),
      data        = di,
      eval.times  = eval_times,
      splitMethod = "none",
      cens.model  = "marginal"
    )
    as.numeric(ci$AppCindex[[1]])
  }

  rownames(Z) <- round(cont_vals, 3)
  colnames(Z) <- paste0("t=", round(eval_times, 2))

  # contour plot
  plt <- plot_ly(
    x = eval_times, y = cont_vals, z = Z, type = "contour",
    colorbar = list(title = "C-index", side = "right"),
    contours = list(coloring = "heatmap"),
    hovertemplate = paste(
      "At time %{x:.2f} <br>with ", contCovName,
      " cutoff %{y:.2f}<br>C-index %{z:.3f}<extra></extra>"
    )
  ) |>
    layout(
      title = list(
        text = paste("C-index Contour by Time and", contCovName),
        font = list(size = 20), x = 0.15
      ),
      xaxis = list(title = "Time", range = c(0, max(eval_times))),
      yaxis = list(title = contCovName, range = range(cont_vals)),
      margin = list(t = 50, l = 5)
    )

  if (isTRUE(drawHistogram)) {
    histData <- data[[contCov]]
    hist_plot <- plot_ly(y = histData, type = "histogram", hoverinfo = "none") |>
      layout(xaxis = list(title = "Count"),
             yaxis = list(title = contCovName))
    subplot(
      plt, hist_plot, nrows = 1, widths = c(0.8, 0.2), margin = 0.01,
      shareY = TRUE, titleX = TRUE, titleY = TRUE
    ) |>
      layout(margin = list(t = 50, l = 5))
  } else {
    plt
  }
}

#' Time-dependent C-index contour for Random Survival Forests
#'
#' @description
#' Computes a matrix of time-dependent C-index values for a fitted
#' \pkg{randomForestSRC} survival forest by perturbing a continuous covariate
#' around a sequence of thresholds and evaluating IPCW C-index via \pkg{pec}.
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
#' @param eps Numeric jitter added around threshold for stability.
#' @param drawHistogram Logical: show histogram of \code{contCov}.
#'
#' @return A Plotly contour plot (matrix of C-index over time × threshold).
#'
#' @examples
#' \dontrun{
#' # vet <- survival::veteran
#' # fit <- randomForestSRC::rfsrc(Surv(time, status) ~ age + karno, data = vet)
#' # rsfCIndex(vet, fit, contCov = "age",
#' #           timeVar = "time", statusVar = "status")
#' }
#'
#' @seealso \code{\link{CIndexContour}}, \code{\link{coxCIndex}}, \code{\link{xgbCIndex}}
#'
#' @export
rsfCIndex <- function(
    data, model, contCov,
    timeVar, statusVar,
    contCovName = NULL,
    nCovEval = 30, nCores = 4,
    eval_times = NULL, eps = 1e-6,
    drawHistogram = TRUE
){
  if (is.null(contCovName)) contCovName <- contCov

  # eval times
  if (is.null(eval_times)) {
    tt <- data[[timeVar]]
    tt <- tt[tt > 0 & is.finite(tt)]
    eval_times <- as.numeric(
      quantile(tt, probs = seq(0.05, 0.95, length.out = 40), na.rm = TRUE)
    )
  }
  eval_times <- sort(unique(eval_times))

  # thresholds
  cutoffs <- seq(
    quantile(data[[contCov]], 0.05, na.rm = TRUE),
    quantile(data[[contCov]], 0.95, na.rm = TRUE),
    length.out = nCovEval
  )

  # parallel
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)

  Z_list <- foreach(i = seq_along(cutoffs)) %dopar% {
    cutoff <- cutoffs[i]
    di <- data
    di[[contCov]] <- ifelse(di[[contCov]] > cutoff, cutoff + eps, cutoff - eps)

    ci <- pec::cindex(
      object      = model,
      formula     = as.formula(paste0("Surv(", timeVar, ",", statusVar, ") ~ 1")),
      data        = di,
      eval.times  = eval_times,
      splitMethod = "none",
      cens.model  = "marginal"
    )
    as.numeric(ci$AppCindex[[1]])
  }

  Z <- do.call(rbind, Z_list)
  rownames(Z) <- round(cutoffs, 3)
  colnames(Z) <- paste0("t=", round(eval_times, 2))

  plt <- plot_ly(
    x = eval_times, y = cutoffs, z = Z, type = "contour",
    contours = list(coloring = "heatmap"),
    colorbar = list(title = "C-index"),
    hovertemplate = paste0(
      "time=%{x:.2f}<br>", contCovName,
      " cutoff=%{y:.2f}<br>C=%{z:.3f}<extra></extra>"
    )
  ) |>
    layout(
      xaxis = list(title = "Time"),
      yaxis = list(title = contCovName),
      title = "RSF time-dependent C-index"
    )

  if (isTRUE(drawHistogram)) {
    histData <- data[[contCov]]
    hist_plot <- plot_ly(y = histData, type = "histogram", hoverinfo = "none") |>
      layout(xaxis = list(title = "Count"),
             yaxis = list(title = contCovName))
    subplot(
      plt, hist_plot, nrows = 1, widths = c(0.8, 0.2), margin = 0.01,
      shareY = TRUE, titleX = TRUE, titleY = TRUE
    ) |>
      layout(margin = list(t = 50, l = 5))
  } else {
    plt
  }
}

#' Time-dependent C-index contour for XGBoost AFT survival models
#'
#' @description
#' Computes IPCW-adjusted C-index using \pkg{pec} for a fitted
#' \pkg{xgboost} AFT survival model across evaluation time and continuous
#' covariate thresholds. Produces a Plotly contour heatmap, optionally
#' with histogram.
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
#' @param sigma Numeric AFT scale parameter used in survival approximation.
#' @param eps Numeric jitter for covariate dichotomization.
#' @param drawHistogram Attach histogram of \code{contCov} if TRUE.
#'
#' @return A Plotly contour heatmap (with or without histogram).
#'
#' @examples
#' \dontrun{
#' # xgb_fit <- xgboost::xgb.train(params = list(objective = "survival:aft"), ...)
#' # xgbCIndex(vet, xgb_fit, contCov = "age",
#' #           statusVar = "status", timeVar = "time")
#' }
#'
#' @seealso \code{\link{CIndexContour}}, \code{\link{coxCIndex}}, \code{\link{rsfCIndex}}
#'
#' @export
xgbCIndex <- function(
    data, model, contCov,
    statusVar, timeVar,
    contCovName = NULL,
    nCovEval = 30,
    eval_times  = NULL,
    sigma = 1.6, eps = 1e-6,
    drawHistogram = TRUE
){
  if (is.null(contCovName)) contCovName <- contCov

  # eval times
  if (is.null(eval_times)) {
    tt <- data[[timeVar]]
    eval_times <- as.numeric(
      quantile(tt[tt > 0],
               probs = seq(0.05, 0.95, length.out = 40),
               na.rm = TRUE)
    )
  }
  eval_times <- unique(sort(eval_times))

  # thresholds
  cont_vals <- seq(
    quantile(data[[contCov]], 0.05, na.rm = TRUE),
    quantile(data[[contCov]], 0.95, na.rm = TRUE),
    length.out = nCovEval
  )

  # build design matrix helper
  feat <- model$feature_names
  build_X <- function(df){
    X <- model.matrix(
      ~ . - 1,
      data = df[, setdiff(names(df), c(timeVar, statusVar)), drop = FALSE]
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

  # define predictSurvProb method for pec in the package namespace
  ns_env <- parent.env(environment())
  assign(
    "predictSurvProb.xgb.Booster",
    function(object, newdata, times, ...) {
      X_new <- build_X(newdata)
      mu    <- as.numeric(
        predict(object, newdata = xgboost::xgb.DMatrix(X_new),
                outputmargin = TRUE)
      )
      tt    <- pmax(as.numeric(times), .Machine$double.eps)
      Zs    <- outer(log(tt), mu, FUN = function(lt, m) (lt - m)/sigma)
      St    <- 1 - pnorm(t(Zs))               # n × K
      pmin(pmax(St, 1e-6), 1 - 1e-6)
    },
    envir = ns_env
  )

  # compute Z (rows = thresholds, cols = times)
  Z <- matrix(NA_real_, nrow = length(cont_vals), ncol = length(eval_times))
  for (i in seq_along(cont_vals)) {
    cutoff <- cont_vals[i]
    di <- data
    di[[contCov]] <- ifelse(di[[contCov]] > cutoff, cutoff + eps, cutoff - eps)

    ci_pec <- pec::cindex(
      object      = model,
      formula     = as.formula(paste0("Surv(", timeVar, ", ", statusVar, ") ~ 1")),
      data        = di,
      eval.times  = eval_times,
      splitMethod = "none",
      cens.model  = "marginal"
    )
    Z[i, ] <- as.numeric(ci_pec$AppCindex[[1]])
  }

  rownames(Z) <- round(cont_vals, 3)
  colnames(Z) <- paste0("t=", round(eval_times, 2))

  plt <- plot_ly(
    x = eval_times, y = cont_vals, z = Z, type = "contour",
    contours = list(coloring = "heatmap"),
    colorbar  = list(title = "C-index"),
    hovertemplate = paste0(
      "time=%{x:.2f}<br>", contCovName,
      " cutoff=%{y:.2f}<br>C=%{z:.3f}<extra></extra>"
    )
  ) |>
    layout(
      xaxis = list(title = "Time"),
      yaxis = list(title = contCovName),
      title = "XGB-AFT time-dependent C-index"
    )

  if (isTRUE(drawHistogram)) {
    histData <- data[[contCov]]
    hist_plot <- plot_ly(y = histData, type = "histogram", hoverinfo = "none") |>
      layout(xaxis = list(title = "Count"),
             yaxis = list(title = contCovName))
    subplot(
      plt, hist_plot, nrows = 1, widths = c(0.8, 0.2), margin = 0.01,
      shareY = TRUE, titleX = TRUE, titleY = TRUE
    ) |>
      layout(margin = list(t = 50, l = 5))
  } else {
    plt
  }
}

# silence R CMD check about foreach loop variable
utils::globalVariables("i")
