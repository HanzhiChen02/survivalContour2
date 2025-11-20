# File: R/effectContour.R

#' Internal imports for effect-size contour functions
#'
#' @name effect_internal_imports
#' @keywords internal
#'
#' @importFrom stats quantile loess model.matrix pnorm qnorm
#' @importFrom plotly plot_ly layout
#' @importFrom randomForestSRC predict.rfsrc
#' @importFrom xgboost xgb.DMatrix
#' @importFrom riskRegression predictRisk
NULL


# ---------- Shared helpers ----------------------------------------------------

# Contour Plot for Estimated Hazard Ratios
contourPart <- function(x, y, z, contCovName) {
  fig <- plot_ly(
    x = x, y = y, z = z, type = "contour",
    colorbar = list(title = "Hazard Ratio", titleside = "right"),
    hovertemplate = paste(
      "At time %{x:.2f} <br>with ", contCovName,
      " = %{y:.2f},<br>HR is %{z:.2f}<extra></extra>"
    )
  ) |>
    layout(
      title = list(
        text = "Contour Plot of the Estimated Hazard Ratio",
        font = list(size = 20), x = 0.15
      ),
      xaxis = list(
        title = "Time", font = list(size = 24), range = c(0, max(x, na.rm = TRUE))
      ),
      yaxis = list(
        title = contCovName, font = list(size = 24), range = range(y, na.rm = TRUE)
      )
    )
  fig
}

# Compute Column-wise Derivative over Time
colDeri <- function(mat, time) {
  time1 <- time[-length(time)]
  time2 <- time[-1]
  timeDiff <- time2 - time1

  mat1 <- mat[, -1, drop = FALSE]
  mat2 <- mat[, -ncol(mat), drop = FALSE]
  matDiff <- mat1 - mat2

  matDiff / matrix(timeDiff, nrow = nrow(mat), ncol = length(timeDiff), byrow = TRUE)
}

#' Compute Row-wise Derivative over Covariate
#'
#' @param mat A matrix of cumulative hazard.
#' @param cov A numeric vector of covariate values.
#' @return A matrix of row-wise derivatives.
rowDeri <- function(mat, cov) {
  cov1 <- cov[-length(cov)]
  cov2 <- cov[-1]
  covDiff <- cov2 - cov1

  mat1 <- mat[-1, , drop = FALSE]
  mat2 <- mat[-nrow(mat), , drop = FALSE]
  matDiff <- mat1 - mat2

  matDiff / matrix(covDiff, nrow = nrow(matDiff), ncol = ncol(matDiff), byrow = FALSE)
}


# ---------- Dispatcher --------------------------------------------------------

#' Unified effect-size contour visualization for survival models
#'
#' @description
#' Computes a time-dependent effect-size surface (typically log-hazard ratio)
#' over follow-up time (x-axis) and a continuous covariate (y-axis) for a range
#' of survival models. Dispatches automatically based on the fitted model:
#' \itemize{
#'   \item \code{flexsurvreg} → \code{paraEffectContour}
#'   \item \code{pycox} models → \code{pycoxEffectContour}
#'   \item \code{rfsrc} → \code{rfsrcEffectContour}
#'   \item \code{xgb.Booster} → \code{xgbEffectContour}
#'   \item \code{CSC} (Fine–Gray) → \code{FGEffectContour}
#' }
#'
#' @param data Data frame used to fit the model.
#' @param model Fitted survival model (one of the classes above).
#' @param contCov Character; name of the continuous covariate of interest.
#' @param contCovName Optional pretty label for the y-axis; defaults to
#'   \code{contCov}.
#' @param timeVar Character; name of time variable (for RSF / XGB backends).
#' @param statusVar Character; name of event indicator (for RSF / XGB backends).
#' @param nCovEval Integer; number of covariate grid points between the 5–95%
#'   quantiles of \code{contCov}.
#' @param otherCov Optional covariate profile for certain backends.
#' @param smooth Logical; whether to 2D smooth the effect surface (when
#'   supported by backend).
#' @param spanx,spany Numeric; LOESS spans in x and y directions for smoothing.
#' @param sigma Optional AFT scale (for XGBoost AFT backend).
#'
#' @return A Plotly contour plot of effect size over time and covariate.
#'
#' @export
effectContour <- function(
    data, model, contCov, contCovName = NULL,
    timeVar, statusVar,
    nCovEval = 30, otherCov = NULL,
    smooth = FALSE, spanx = 0.75, spany = 0.75,
    sigma = NULL
) {
  if (is.null(contCovName)) contCovName <- contCov

  if (inherits(model, "flexsurvreg")) {
    # Parametric survival model
    Plot <- paraEffectContour(
      data        = data,
      model       = model,
      contCov     = contCov,
      contCovName = contCovName,
      nCovEval    = nCovEval,
      otherCov    = otherCov,
      smooth      = smooth,
      spanx       = spanx,
      spany       = spany
    )

  } else if (any(grepl("pycox", class(model)))) {
    # Deep learning survival model
    Plot <- pycoxEffectContour(
      data        = data,
      trainModel  = model,
      contCov     = contCov,
      contCovName = contCovName,
      nCovEval    = nCovEval,
      otherCov    = otherCov,
      smooth      = smooth,
      spanx       = spanx,
      spany       = spany
    )

  } else if (any(grepl("rfsrc", class(model)))) {
    # Random survival forest model
    Plot <- rfsrcEffectContour(
      data        = data,
      model       = model,
      contCov     = contCov,
      contCovName = contCovName,
      timeVar     = timeVar,
      statusVar   = statusVar,
      nCovEval    = nCovEval,
      otherCov    = otherCov,
      smooth      = smooth,
      spanx       = spanx,
      spany       = spany
    )

  } else if (inherits(model, "xgb.Booster")) {
    # XGBoost AFT survival model
    Plot <- xgbEffectContour(
      data        = data,
      model       = model,
      contCov     = contCov,
      contCovName = contCovName,
      statusVar   = statusVar,
      timeVar     = timeVar,
      nCovEval    = nCovEval,
      sigma       = sigma
    )

  } else if (inherits(model, "CSC")) {
    # Fine & Gray model
    Plot <- FGEffectContour(
      data        = data,
      trainModel  = model,
      contCov     = contCov,
      contCovName = contCovName,
      nCovEval    = nCovEval,
      otherCov    = otherCov,
      smooth      = smooth,
      spanx       = spanx,
      spany       = spany
    )

  } else {
    stop("The model provided cannot be handled by effectContour.")
  }

  Plot
}


# ---------- Fine & Gray backend ----------------------------------------------

FGEffectContour <- function(
    data, trainModel, contCov, contCovName = NULL,
    nCovEval = 30, otherCov = NULL,
    smooth = FALSE, spanx = NA, spany = NA
) {
  if (is.null(otherCov)) {
    newData <- data.frame(data[1, , drop = FALSE])
    for (i in seq_along(data)) {
      if (is.numeric(data[[i]]) || is.integer(data[[i]])) {
        newData[[i]] <- mean(data[[i]], na.rm = TRUE)
      } else {
        newData[[i]] <- as.character(na.omit(data[[i]]))[which.max(table(data[[i]]))]
      }
    }
  } else {
    newData <- otherCov
  }

  cont_range <- stats::quantile(data[[contCov]], probs = c(0.025, 0.975), na.rm = TRUE)
  cont_seq <- seq(from = cont_range[1], to = cont_range[2], length.out = nCovEval)

  newData0 <- newData[rep(1, nCovEval), ]
  newData0[[contCov]] <- cont_seq

  ## time grid
  ctime <- trainModel$response[trainModel$response[, 2] == trainModel$cause, 1]
  ctime <- sort(unique(ctime))
  if (ctime[1] != 0) ctime <- c(0, ctime)

  ## cumulative incidence
  pred0 <- 1 - riskRegression::predictRisk(
    object = trainModel, newdata = newData0, times = ctime
  )

  cHaz <- -log(pred0)
  haz  <- colDeri(cHaz, ctime)
  loghaz <- log(haz)
  loghr  <- rowDeri(cHaz, cont_seq)

  if (is.null(contCovName)) contCovName <- contCov

  if (isTRUE(smooth)) {
    loghr <- apply(loghr, 2, function(vec)
      stats::loess(vec ~ I(1:length(vec)), span = spany)$fitted)
    loghr <- apply(loghr, 1, function(vec)
      stats::loess(vec ~ I(1:length(vec)), span = spanx)$fitted)
    loghr <- t(loghr)
  }

  contourPart(ctime, cont_seq, loghr, contCovName)
}


# ---------- Utility helpers (kept) ---------------------------------------

# Match column types of df2 to df1
match_types <- function(df1, df2) {
  df2[] <- lapply(seq_along(df1), function(i) {
    if (is.integer(df1[[i]])) {
      as.integer(df2[[i]])
    } else if (is.numeric(df1[[i]])) {
      as.numeric(df2[[i]])
    } else if (is.factor(df1[[i]])) {
      as.factor(df2[[i]])
    } else {
      df2[[i]]
    }
  })
  df2
}

predictPhreg <- function(object,newdata,times=NULL,individual.time=FALSE,tminus=FALSE,se=TRUE,robust=FALSE,conf.type="log",conf.int=0.95,km=FALSE,...)
{# {{{ default is all time-points from the object

  ### take baseline and strata from object# {{{
  strata <- object$strata[object$jumps]
  nstrata <- object$nstrata
  jumptimes <- object$cumhaz[,1]
  chaz <- object$cumhaz[,2]
  if (se) {
    if (!robust) {
      se.chaz <- object$se.cumhaz[,2]
      #varbeta <- object$ihessian
      varbeta <- object$II
      Pt <- apply(object$E/c(object$S0),2,cumsumstrata,strata,nstrata)
    } else {
      if (is.null(object$opt) | is.null(object$coef)) fixbeta<- 1 else fixbeta <- 0
      IsdM <- squareintHdM(object,ft=NULL,fixbeta=fixbeta,...)
      ###
      se.chaz <-   IsdM$varInt[object$jumps]^.5
      covv <- IsdM$covv[object$jumps,,drop=FALSE]
      varbeta <- IsdM$vbeta
      Pt <- IsdM$Ht[object$jumps,,drop=FALSE]
    }
  } # }}}


  ### setting up newdata with factors and strata
  desX <- readPhreg(object,newdata)
  X <- desX$X
  strataNew <- desX$strata

  if (is.null(times)) times <- sort(unique(c(object$exit)))
  if (individual.time & is.null(times)) times <- c(object$exit)
  if (individual.time & length(times)==1) times <- rep(times,length(object$exit))

  se.cumhaz <- NULL
  if (!individual.time) {
    pcumhaz <- surv <- matrix(0,nrow(X),length(times))
    if (se) se.cumhaz <- matrix(0,nrow(X),length(times))
  } else {
    pcumhaz <- surv <- matrix(0,nrow(X),1)
    if (se) se.cumhaz <- matrix(0,nrow(X),1)
  }
  hazt <- length(times)

  for (j in unique(strataNew)) {
    where <- sindex.prodlim(c(0,jumptimes[strata==j]),times,strict=tminus)
    plhazt <- hazt <- c(0,chaz[strata==j])
    if (km) { plhazt <- c(1,exp(cumsum(log(1-diff(hazt)))));  plhazt[is.na(hazt)] <- 0 }
    if (se) se.hazt <- c(0,se.chaz[strata==j])
    Xs <- X[strataNew==j,,drop=FALSE]
    ###	offs <- object$offsets[object$strata==j]
    if (object$p==0) RR <- rep(1,nrow(Xs)) else RR <- c(exp( Xs %*% coef(object)))
    if (se)  { # {{{ based on Hazard's
      if (object$p>0) {
        Ps <- Pt[strata==j,,drop=FALSE]
        Ps <- rbind(0,Ps)[where,,drop=FALSE]
        #  print(Xs); print(varbeta); print(dim(Ps)); print((Xs %*% varbeta))
        Xbeta <- Xs %*% varbeta
        seXbeta <- rowSums(Xbeta*Xs)^.5
        cov2 <- cov1 <- Xbeta %*% t(Ps*hazt[where])
        if (robust)	{
          covvs <- covv[strata==j,,drop=FALSE]
          covvs <- rbind(0,covvs)[where,,drop=FALSE]
          covv1 <- Xs %*% t((covvs*hazt[where]))
          cov1 <- cov1-covv1
        }
      } else cov1 <- 0
    }# }}}
    haztw <- hazt[where]
    if (se) se.haztw <- se.hazt[where]
    if (is.null(object$propodds)) {
      plhaztw <- plhazt[where]
      if (!individual.time) pcumhaz[strataNew==j,]  <- RR%o%haztw else pcumhaz[strataNew==j,] <- RR*haztw[strataNew==j]
      if (!km) {
        if (!individual.time) surv[strataNew==j,]  <- exp(- RR%o%haztw)
        else surv[strataNew==j,]  <- exp(-RR*haztw[strataNew==j])
      } else {
        if (!individual.time) surv[strataNew==j,]  <- exp( RR %o% log(plhaztw))
        else surv[strataNew==j,]  <- plhaztw[strataNew==j]^RR
      }
    } else {
      if (!individual.time) surv[strataNew==j,]  <- 1/(1+RR%o%haztw)
      else surv[strataNew==j,]  <- 1/(1+RR*haztw[strataNew==j])
    }
    if (se) {# {{{
      if (object$p>0)  {
        if (!individual.time) se.cumhaz[strataNew==j,]  <-
            ((RR %o% se.haztw)^2+(c(RR*seXbeta) %o% haztw)^2-2*RR^2*cov1)^.5
        else se.cumhaz[strataNew==j,]  <- RR* (se.haztw^2+(c(seXbeta)*haztw)^2-2*diag(cov1))^.5
      } else {
        if (!individual.time) se.cumhaz[strataNew==j,]  <- RR %o% (se.haztw)
        else se.cumhaz[strataNew==j,]  <- RR* se.haztw[strataNew==j]
      }
    }# }}}
  }


  zval <- qnorm(1 - (1 - conf.int)/2, 0, 1)
  std.err <- se.cumhaz
  cisurv  <- list()
  cisurv$upper <- NULL
  cisurv$lower <- NULL

  ### different conf-types for surv
  if (se) {# {{{
    if (conf.type == "plain") {# {{{
      temp1 <- surv + zval * std.err * surv
      temp2 <- surv - zval * std.err * surv
      cisurv <- list(upper = pmin(temp1, 1), lower = pmax(temp2,
                                                          0), conf.type = "plain", conf.int = conf.int)
    }
    if (conf.type == "log") {
      xx <- ifelse(surv == 0, 1, surv)
      temp1 <- ifelse(surv == 0, NA, exp(log(xx) + zval * std.err))
      temp2 <- ifelse(surv == 0, NA, exp(log(xx) - zval * std.err))
      cisurv <- list(upper = pmin(temp1, 1), lower = temp2,
                     conf.type = "log", conf.int = conf.int)
    }
    if (conf.type == "log-log") {
      who <- (surv == 0 | surv == 1)
      temp3 <- ifelse(surv == 0, NA, 1)
      xx <- ifelse(who, 0.1, surv)
      temp1 <- exp(-exp(log(-log(xx)) + zval * std.err/log(xx)))
      temp1 <- ifelse(who, temp3, temp1)
      temp2 <- exp(-exp(log(-log(xx)) - zval * std.err/log(xx)))
      temp2 <- ifelse(who, temp3, temp2)
      cisurv <- list(upper = temp1, lower = temp2,
                     conf.type = "log-log", conf.int = conf.int)
    }# }}}
  }# }}}

  if (object$p>0) RR <-  exp(X %*% coef(object)) else RR <- rep(1,nrow(X))

  ### non-cox setting
  if (!is.null(object$propodds)) pcumhaz <- -log(surv)

  out <- list(surv=surv,times=times,
              upper=cisurv$upper,lower=cisurv$lower,cumhaz=pcumhaz,se.cumhaz=se.cumhaz,strata=strataNew,X=X, RR=RR)
  if (length(class(object))==2 && substr(class(object)[2],1,3)=="cif") {
    out <- c(out,list(cif=1-out$surv,cif.lower=1-out$upper, cif.upper=1-out$lower))
  }
  class(out) <- c("predictphreg")
  if (length(class(object))==2) class(out) <- c("predictphreg",class(object)[2])
  return(out)
}# }}}

# ---------- PyCox backend -----------------------------------------------------

pycoxEffectContour <- function(
    data, trainModel, contCov, contCovName = NULL,
    nCovEval = 30, otherCov = NULL,
    smooth = FALSE, spanx = NA, spany = NA
) {
  if (is.null(contCovName)) contCovName <- contCov

  cont_range <- stats::quantile(data[[contCov]], probs = c(0.025, 0.975), na.rm = TRUE)
  cont_seq <- seq(from = cont_range[1], to = cont_range[2], length.out = nCovEval)

  newData0 <- data[rep(1, nCovEval), , drop = FALSE]
  newData0[[contCov]] <- cont_seq

  if (!is.null(otherCov)) {
    for (col in otherCov) {
      if (is.numeric(data[[col]]) || is.integer(data[[col]])) {
        newData0[[col]] <- mean(data[[col]], na.rm = TRUE)
      } else {
        newData0[[col]] <- names(which.max(table(data[[col]])))
      }
    }
  }

  surv0 <- predict(trainModel, newdata = newData0, type = "survival")
  time_seq <- as.numeric(colnames(surv0))
  surv_mat <- as.matrix(surv0)

  cHaz <- -log(surv_mat)
  haz  <- colDeri(cHaz, time_seq)
  loghaz <- log(haz)
  loghr  <- rowDeri(cHaz, cont_seq)

  if (isTRUE(smooth)) {
    loghr <- apply(loghr, 2, function(vec)
      stats::loess(vec ~ I(1:length(vec)), span = spany)$fitted)
    loghr <- apply(loghr, 1, function(vec)
      stats::loess(vec ~ I(1:length(vec)), span = spanx)$fitted)
    loghr <- t(loghr)
  }

  contourPart(
    x = time_seq[-length(time_seq)],
    y = cont_seq[-length(cont_seq)],
    z = loghr,
    contCovName = contCovName
  )
}


# ---------- RSF backend -------------------------------------------------------

rfsrcEffectContour <- function(
    data, model, contCov, contCovName = NULL,
    timeVar = NULL, statusVar = NULL,
    nCovEval = 30, otherCov = NULL,
    smooth = FALSE, spanx = 0.75, spany = 0.75,
    mode = "cumhaz"
) {
  mode <- match.arg(mode)
  if (is.null(contCovName)) contCovName <- contCov

  xvars <- model$xvar.names
  xdf   <- data[, xvars, drop = FALSE]

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
  newData2 <- newData[rep(1, nCovEval), , drop = FALSE]
  newData2[[contCov]] <- cov_vals

  predObj <- randomForestSRC::predict.rfsrc(model, newdata = newData2)
  S  <- pmin(pmax(as.matrix(predObj$survival), 1e-12), 1 - 1e-12)
  tt <- as.numeric(predObj$time.interest)

  cHaz <- -log(S)
  haz  <- colDeri(cHaz, tt)
  time_axis <- tt[-1]

  if (mode == "cumhaz") {
    dCH_dc <- rowDeri(cHaz, cov_vals)
    Z <- dCH_dc[, -ncol(dCH_dc), drop = FALSE]
    cov_axis <- cov_vals[-1]
  } else {
    loghaz  <- log(pmax(haz, 1e-12))
    dloghdc <- rowDeri(loghaz, cov_vals)
    Z <- dloghdc
    cov_axis <- (cov_vals[-1] + cov_vals[-length(cov_vals)]) / 2
  }

  if (isTRUE(smooth)) {
    Z <- apply(Z, 2, function(v)
      stats::loess(v ~ I(seq_along(v)), span = spany)$fitted)
    Z <- apply(Z, 1, function(v)
      stats::loess(v ~ I(seq_along(v)), span = spanx)$fitted)
    Z <- t(Z)
  }

  contourPart(time_axis, cov_axis, Z, contCovName)
}


# ---------- XGBoost backend ---------------------------------------------------

xgbEffectContour <- function(
    data, model, contCov, contCovName = NULL,
    statusVar, timeVar, nCovEval = 30, eval_times = NULL,
    sigma = 1.6, covStep = NULL, agg = c("median", "mean")
) {
  if (is.null(contCovName)) contCovName <- contCov
  agg <- match.arg(agg)

  # eval times
  if (is.null(eval_times)) {
    tt <- data[[timeVar]]
    eval_times <- as.numeric(stats::quantile(
      tt[tt > 0],
      probs = seq(0.05, 0.95, length.out = 40),
      na.rm = TRUE
    ))
  }
  eval_times <- unique(sort(eval_times))

  # covariate grid
  cov_vals <- seq(
    stats::quantile(data[[contCov]], 0.05, na.rm = TRUE),
    stats::quantile(data[[contCov]], 0.95, na.rm = TRUE),
    length.out = nCovEval
  )
  rng <- range(data[[contCov]], na.rm = TRUE)

  feat <- model$feature_names
  build_X <- function(df) {
    X <- stats::model.matrix(
      ~ . - 1,
      data = df[, setdiff(names(df), c(timeVar, statusVar)), drop = FALSE]
    )
    if (!is.null(feat)) {
      miss <- setdiff(feat, colnames(X))
      if (length(miss)) {
        X <- cbind(X, matrix(0, nrow(X), length(miss), dimnames = list(NULL, miss)))
      }
      X <- X[, feat, drop = FALSE]
    }
    X
  }
  X0 <- build_X(data)
  j_col <- match(contCov, colnames(X0))
  if (is.na(j_col)) {
    stop(sprintf("Column '%s' not found in design matrix.", contCov))
  }

  # log-hazard(t | mu, sigma) for log-normal AFT
  loghaz_from_mu <- function(mu_vec, times, sigma) {
    n <- length(mu_vec); K <- length(times)
    lt <- log(pmax(times, .Machine$double.eps))
    Z  <- (matrix(lt, n, K, byrow = TRUE) - matrix(mu_vec, n, K)) / sigma
    stats::dnorm(Z, log = TRUE) -
      log(sigma) -
      matrix(lt, n, K, byrow = TRUE) -
      stats::pnorm(Z, lower.tail = FALSE, log.p = TRUE)
  }

  if (is.null(covStep)) {
    covStep <- max(0.2 * stats::sd(data[[contCov]], na.rm = TRUE), 1.0)
  }

  Z <- matrix(NA_real_, nrow = length(cov_vals), ncol = length(eval_times))

  for (i in seq_along(cov_vals)) {
    v  <- cov_vals[i]
    hi <- min(v + covStep, rng[2])
    lo <- max(v - covStep, rng[1])
    grp <- data[[contCov]] > v

    X_plus  <- X0
    X_minus <- X0
    X_plus[,  j_col] <- ifelse(grp, hi, lo)
    X_minus[, j_col] <- ifelse(grp, lo, hi)

    mu_plus <- as.numeric(predict(
      model, newdata = xgboost::xgb.DMatrix(X_plus), outputmargin = TRUE
    ))
    mu_minus <- as.numeric(predict(
      model, newdata = xgboost::xgb.DMatrix(X_minus), outputmargin = TRUE
    ))

    loghaz_plus  <- loghaz_from_mu(mu_plus,  eval_times, sigma)
    loghaz_minus <- loghaz_from_mu(mu_minus, eval_times, sigma)

    logHR_mat <- loghaz_plus - loghaz_minus

    Z[i, ] <- if (agg == "median") {
      apply(logHR_mat, 2, median, na.rm = TRUE)
    } else {
      colMeans(logHR_mat, na.rm = TRUE)
    }
  }

  rownames(Z) <- round(cov_vals, 3)
  colnames(Z) <- paste0("t=", round(eval_times, 2))

  contourPart(eval_times, cov_vals, Z, contCovName)
}
