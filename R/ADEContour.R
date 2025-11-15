# File: R/ADEContour.R

#' Internal imports for ADE contour functions
#'
#' @name ade_internal_imports
#' @keywords internal
#'
#' @importFrom stats quantile approx predict model.matrix pnorm
#' @importFrom survival Surv coxph survfit basehaz
#' @importFrom randomForestSRC predict.rfsrc
#' @importFrom plotly plot_ly layout subplot
#' @importFrom xgboost xgb.DMatrix
NULL


#' Time-dependent Average Derivative Effect (ADE) contour
#'   (dispatcher for Cox / RSF / XGB)
#'
#' @description
#' Generic wrapper for computing a contour surface of the Average
#' Derivative Effect (ADE) on survival probability over follow-up
#' time (x-axis) and a continuous covariate (y-axis). It automatically
#' dispatches to the appropriate backend based on `model`.
#'
#' \itemize{
#'   \item \code{coxph}: uses \code{\link{coxADE}}
#'   \item \code{rfsrc}: uses \code{\link{rsfADE}}
#'   \item \code{xgb.Booster}: uses \code{\link{xgbADE}}
#' }
#'
#' @param data A \code{data.frame} containing variables used in the fitted model.
#' @param model A fitted survival model: one of \code{coxph}, \code{rfsrc},
#'   or \code{xgb.Booster}.
#' @param contCov Character; name of the continuous covariate of interest.
#' @param statusVar Character; name of the event indicator (0/1).
#' @param timeVar Character; name of the follow-up time variable.
#' @param contCovName Optional y-axis label for the covariate. Defaults to
#'   \code{contCov}.
#' @param nCovEval Integer; number of thresholds between the 5--95\% quantiles
#'   of \code{contCov}.
#' @param eval_times Optional numeric vector of time points at which to
#'   evaluate survival probabilities (used by Cox/XGB backends).
#' @param otherCov Optional single-row data frame or list specifying values
#'   of other covariates for \code{rsfADE}; if \code{NULL}, they are set to
#'   typical values (mean / mode).
#'
#' @return
#' A \pkg{plotly} object containing a contour heatmap of ADE
#' over time and \code{contCov}.
#'
#' @seealso \code{\link{coxADE}}, \code{\link{rsfADE}}, \code{\link{xgbADE}}
#'
#' @export
ADEContour <- function(
    data, model, contCov,
    statusVar, timeVar,
    contCovName = NULL,
    nCovEval = 30,
    eval_times = NULL,
    otherCov = NULL
) {
  if (is.null(contCovName)) contCovName <- contCov

  if (inherits(model, "coxph")) {
    coxADE(
      data        = data,
      model       = model,
      contCov     = contCov,
      contCovName = contCovName,
      statusVar   = statusVar,
      timeVar     = timeVar,
      eval_times  = eval_times,
      nCovEval    = nCovEval
    )
  } else if (any(grepl("rfsrc", class(model)))) {
    rsfADE(
      data        = data,
      model       = model,
      contCov     = contCov,
      contCovName = contCovName,
      timeVar     = timeVar,
      statusVar   = statusVar,
      nCovEval    = nCovEval,
      otherCov    = otherCov
    )
  } else if (inherits(model, "xgb.Booster")) {
    xgbADE(
      data        = data,
      model       = model,
      contCov     = contCov,
      contCovName = contCovName,
      statusVar   = statusVar,
      timeVar     = timeVar,
      nCovEval    = nCovEval,
      eval_times  = eval_times
    )
  } else {
    stop("Unsupported model class: ", paste(class(model), collapse = ", "))
  }
}


# ---------- shared plotting helper ----------

#' @keywords internal
contourPart_ADE <- function(x, y, z, contCovName) {
  fig <- plot_ly(
    x = x, y = y, z = z, type = "contour",
    colorbar = list(
      title = "Delta Survival Probability",
      titleside = "right"
    ),
    hovertemplate = paste(
      "At time %{x:.2f}<br>",
      contCovName, " = %{y:.2f}<br>",
      "Delta S (per +1 in ", contCovName, ") = %{z:.4f}",
      "<extra></extra>"
    )
  ) |>
    layout(
      title = list(
        text = "Contour of Average Derivative Effect (ADE) on Survival Probability",
        font = list(size = 20),
        x    = 0.15
      ),
      xaxis = list(
        title = "Time",
        font  = list(size = 24),
        range = c(0, max(x, na.rm = TRUE))
      ),
      yaxis = list(
        title = contCovName,
        font  = list(size = 24),
        range = range(y, na.rm = TRUE)
      )
    )

  fig
}


# ---------- Cox backend ----------

#' Time-dependent ADE for Cox proportional hazards models
#'
#' @inheritParams ADEContour
#' @param eval_times Optional numeric vector of time points for evaluation.
#'
#' @return A Plotly contour heatmap of ADE.
#'
#' @export
coxADE <- function(
    data, model, contCov,
    contCovName = NULL,
    statusVar, timeVar,
    eval_times = NULL,
    nCovEval = 30
) {
  if (is.null(contCovName)) contCovName <- contCov

  ## 1) time
  if (is.null(eval_times)) {
    sf <- survival::survfit(model, newdata = data)
    eval_times <- sort(unique(sf$time))
  } else {
    eval_times <- sort(unique(as.numeric(eval_times)))
  }

  ## 2) covariate grid
  cov_vals <- seq(
    as.numeric(stats::quantile(data[[contCov]], 0.05, na.rm = TRUE)),
    as.numeric(stats::quantile(data[[contCov]], 0.95, na.rm = TRUE)),
    length.out = nCovEval
  )

  ## 3) anchor profile for other covariates
  anchor <- data[1, , drop = TRUE]
  for (nm in names(anchor)) {
    if (nm %in% c(timeVar, statusVar, contCov)) next
    v <- data[[nm]]
    if (is.numeric(v) || is.integer(v)) {
      anchor[[nm]] <- mean(v, na.rm = TRUE)
    } else if (is.factor(v)) {
      anchor[[nm]] <- factor(names(which.max(table(v))), levels = levels(v))
    } else {
      tab <- table(v)
      anchor[[nm]] <- as.character(names(which.max(tab)))
    }
  }
  newData0 <- as.data.frame(anchor)[rep(1, nCovEval), , drop = FALSE]
  newData0[[contCov]] <- cov_vals
  newData1 <- newData0
  newData1[[contCov]] <- cov_vals + 1  # Delta= +1

  ## 4) baseline cumulative hazard H0(t)
  bh <- survival::basehaz(model, centered = FALSE)
  H0 <- stats::approx(
    x = bh$time, y = bh$hazard, xout = eval_times,
    method = "linear", rule = 2, ties = "ordered"
  )$y
  H0 <- pmax(H0, 0)

  ## 5) linear predictors eta(x)
  lp0 <- as.numeric(stats::predict(model, newdata = newData0, type = "lp"))
  lp1 <- as.numeric(stats::predict(model, newdata = newData1, type = "lp"))

  ## 6) survival probabilities S(t | x) = exp(-H0(t) * exp(eta))
  make_S <- function(lp_vec, H0_vec) {
    n <- length(lp_vec); k <- length(H0_vec)
    eta <- matrix(lp_vec, nrow = n, ncol = k)
    H   <- matrix(H0_vec, nrow = n, ncol = k, byrow = TRUE)
    S   <- exp(-H * exp(eta))
    pmin(pmax(S, 1e-12), 1 - 1e-12)
  }
  S0 <- make_S(lp0, H0)
  S1 <- make_S(lp1, H0)

  ## 7) ADE surface
  Z <- S1 - S0
  rownames(Z) <- round(cov_vals, 3)
  colnames(Z) <- paste0("t=", round(eval_times, 3))

  ## 8) plot
  contourPart_ADE(eval_times, cov_vals, Z, contCovName)
}


# ---------- RSF backend ----------

#' Time-dependent ADE for Random Survival Forests
#'
#' @inheritParams ADEContour
#'
#' @return A Plotly contour heatmap of ADE.
#'
#' @export
rsfADE <- function(
    data, model, contCov,
    contCovName = NULL,
    timeVar = NULL, statusVar = NULL,
    nCovEval = 30,
    otherCov = NULL
) {
  if (is.null(contCovName)) contCovName <- contCov

  xvars <- model$xvar.names
  xdf   <- data[, xvars, drop = FALSE]

  if (!contCov %in% xvars)
    stop("contCov not found in model$xvar.names")

  # anchor for other covariates
  if (is.null(otherCov)) {
    newData <- xdf[1, , drop = FALSE]
    for (nm in names(newData)) {
      v <- xdf[[nm]]
      if (is.numeric(v) || is.integer(v)) {
        newData[[nm]] <- mean(v, na.rm = TRUE)
      } else if (is.factor(v)) {
        newData[[nm]] <- factor(names(which.max(table(v))), levels = levels(v))
      } else {
        newData[[nm]] <- names(which.max(table(v)))
      }
    }
  } else {
    newData <- as.data.frame(otherCov, stringsAsFactors = FALSE)[, xvars, drop = FALSE]
    miss <- setdiff(xvars, names(newData))
    for (m in miss) newData[[m]] <- xdf[[m]][1]
    newData <- newData[, xvars, drop = FALSE]
  }

  cov_vals <- seq(
    as.numeric(stats::quantile(xdf[[contCov]], 0.05, na.rm = TRUE)),
    as.numeric(stats::quantile(xdf[[contCov]], 0.95, na.rm = TRUE)),
    length.out = nCovEval
  )

  newData0 <- newData[rep(1, nCovEval), , drop = FALSE]
  newData0[[contCov]] <- cov_vals

  newData1 <- newData0
  newData1[[contCov]] <- cov_vals + 1  # Delta= +1

  # RSF survival probabilities
  pred0 <- randomForestSRC::predict.rfsrc(model, newdata = newData0)
  pred1 <- randomForestSRC::predict.rfsrc(model, newdata = newData1)

  S0 <- pmin(pmax(as.matrix(pred0$survival), 1e-12), 1 - 1e-12)
  S1 <- pmin(pmax(as.matrix(pred1$survival), 1e-12), 1 - 1e-12)
  tt <- as.numeric(pred0$time.interest)

  # ADE: S(x+1) - S(x)
  Z <- S1 - S0

  contourPart_ADE(tt, cov_vals, Z, contCovName)
}


# ---------- XGBoost backend ----------

#' Time-dependent ADE for XGBoost AFT survival models
#'
#' @inheritParams ADEContour
#'
#' @return A Plotly contour heatmap of ADE.
#'
#' @export
xgbADE <- function(
    data, model, contCov,
    contCovName = NULL,
    statusVar, timeVar,
    nCovEval = 30,
    eval_times = NULL
) {
  if (is.null(contCovName)) contCovName <- contCov

  # 0) evaluation times
  if (is.null(eval_times)) {
    tt <- data[[timeVar]]
    eval_times <- as.numeric(stats::quantile(
      tt[tt > 0],
      probs = seq(0.05, 0.95, length.out = 40),
      na.rm = TRUE
    ))
  }
  eval_times <- unique(sort(eval_times))

  # 1) covariate grid
  cov_vals <- seq(
    as.numeric(stats::quantile(data[[contCov]], 0.05, na.rm = TRUE)),
    as.numeric(stats::quantile(data[[contCov]], 0.95, na.rm = TRUE)),
    length.out = nCovEval
  )

  # 2) anchor profile
  anchor_row <- data[1, , drop = FALSE]
  for (nm in names(anchor_row)) {
    if (nm %in% c(timeVar, statusVar)) next
    v <- data[[nm]]
    if (is.numeric(v) || is.integer(v)) {
      anchor_row[[nm]] <- mean(v, na.rm = TRUE)
    } else if (is.factor(v)) {
      anchor_row[[nm]] <- factor(names(which.max(table(v))), levels = levels(v))
    } else {
      tab <- table(v)
      anchor_row[[nm]] <- as.character(names(which.max(tab)))
    }
  }

  newData0 <- anchor_row[rep(1, nCovEval), , drop = FALSE]
  newData0[[contCov]] <- cov_vals
  newData1 <- newData0
  newData1[[contCov]] <- cov_vals + 1  # Delta= +1

  # 3) design matrix
  feat <- model$feature_names
  build_X <- function(df) {
    X <- stats::model.matrix(
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
  X0 <- build_X(newData0)
  X1 <- build_X(newData1)

  # 4) sigma (AFT scale parameter)
  get_sigma <- function(m) {
    s <- tryCatch(m$params$aft_loss_distribution_scale,
                  error = function(e) NULL)
    if (is.null(s)) {
      s <- tryCatch(m$call$params$aft_loss_distribution_scale,
                    error = function(e) NULL)
    }
    if (is.null(s)) s <- 1.6
    as.numeric(s)
  }
  sigma <- get_sigma(model)

  # 5) survival S(t) = 1 - Phi((log t - mu)/sigma)
  surv_from_mu <- function(mu_vec, times, sigma) {
    n <- length(mu_vec); K <- length(times)
    lt <- log(pmax(times, .Machine$double.eps))
    Z  <- (matrix(lt, n, K, byrow = TRUE) -
             matrix(mu_vec, n, K)) / sigma
    stats::pnorm(Z, lower.tail = FALSE)
  }

  # 6) S0/S1
  mu0 <- as.numeric(predict(
    model, newdata = xgboost::xgb.DMatrix(X0), outputmargin = TRUE
  ))
  mu1 <- as.numeric(predict(
    model, newdata = xgboost::xgb.DMatrix(X1), outputmargin = TRUE
  ))

  S0 <- surv_from_mu(mu0, eval_times, sigma)
  S1 <- surv_from_mu(mu1, eval_times, sigma)

  # 7) ADE matrix
  Z <- S1 - S0
  rownames(Z) <- round(cov_vals, 3)
  colnames(Z) <- paste0("t=", round(eval_times, 2))

  # 8) plot
  contourPart_ADE(eval_times, cov_vals, Z, contCovName)
}
