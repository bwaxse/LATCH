#####################################################################
# TITLE: Long COVID Phenotyping
# PURPOSE:
# - Simulate a semi-supervised classification setting
# - Generate surrogate (u099.flag) and true outcome (Y)
# - Train XGBoost surrogate model under multiple settings
# - Combine surrogate and EHR features in LOOCV with adaptive LASSO
# - Evaluate prediction using ROC/AUC
#####################################################################

#### =============================== ####
####      0. SETUP AND LIBRARIES     ####
#### =============================== ####
rm(list = ls())
library(xgboost)
library(pROC)
library(dplyr)
library(glmnet)
library(glmpath)

#### =============================== ####
#### 1. HELPER FUNCTIONS (AUXILIARY) ####
#### =============================== ####

# Logistic and inverse logit functions
g.logit <- function(xx) { exp(xx) / (1 + exp(xx)) }
logit   <- function(xx) { log(xx / (1 - xx)) }

# Utility matrix replicator
VTM <- function(vc, dm) {
  matrix(vc, ncol = length(vc), nrow = dm, byrow = TRUE)
}

# Rank-based cumulative sum helper
sum.I <- function(yy, FUN, Yi, Vi = NULL) {
  if (FUN == "<" | FUN == ">=") { yy <- -yy; Yi <- -Yi }
  pos <- rank(c(yy, Yi), ties.method = 'f')[1:length(yy)] - rank(yy, ties.method = 'f')
  if (substring(FUN, 2, 2) == "=") pos <- length(Yi) - pos
  if (!is.null(Vi)) {
    tmpind <- if (substring(FUN, 2, 2) == "=") order(-Yi) else order(Yi)
    Vi <- apply(as.matrix(Vi)[tmpind, , drop = FALSE], 2, cumsum)
    return(rbind(0, Vi)[pos + 1, ])
  } else return(pos)
}

# Empirical CDF estimate
S.FUN <- function(yy, Yi, Di, yes.smooth = FALSE) {
  (sum.I(yy, "<", Yi, Vi = Di) + sum.I(yy, "<=", Yi, Vi = Di)) / sum(Di) / 2
}

# ROC and AUC estimation from predicted probabilities
ROC.Est.FUN <- function(Di, yyi, yy0, fpr0 = NULL, wgti = NULL, yes.smooth = FALSE) {
  if (is.null(wgti)) wgti <- rep(1, length(Di))
  yyi <- as.matrix(yyi); pp <- ncol(yyi)
  mu0 <- sum(wgti * (1 - Di)) / sum(wgti); mu1 <- 1 - mu0
  out.yy <- out.pp <- out.AUC <- out.TPR <- out.FPR <- out.PPV <- out.NPV <- NULL
  
  for (k in 1:pp) {
    yy <- yy0
    if (!is.null(fpr0)) {
      tpr.all <- S.FUN(yyi[, k], Yi = yyi[, k], Di * wgti)
      fpr.all <- S.FUN(yyi[, k], Yi = yyi[, k], (1 - Di) * wgti)
      TPR <- approx(c(0, fpr.all, 1), c(0, tpr.all, 1), fpr0, method = "linear", rule = 2)$y
      TPR <- c(S.FUN(yy0, Yi = yyi[, k], Di * wgti), TPR)
      yy <- c(yy, Sinv.FUN(fpr0, Yi = yyi[, k], (1 - Di) * wgti))
      FPR <- S.FUN(yy, Yi = yyi[, k], (1 - Di) * wgti)
    } else {
      TPR <- S.FUN(yy, Yi = yyi[, k], Di * wgti)
      FPR <- S.FUN(yy, Yi = yyi[, k], (1 - Di) * wgti)
    }
    out.yy <- cbind(out.yy, yy)
    out.pp <- cbind(out.pp, S.FUN(yy, Yi = yyi[, k], wgti))
    out.TPR <- cbind(out.TPR, TPR); out.FPR <- cbind(out.FPR, FPR)
    PPV <- 1 / (1 + FPR * mu0 / (TPR * mu1)); NPV <- 1 / (1 + (1 - TPR) * mu1 / ((1 - FPR) * mu0))
    out.PPV <- cbind(out.PPV, PPV); out.NPV <- cbind(out.NPV, NPV)
    AUC <- sum(S.FUN(yyi[, k], Yi = yyi[, k], Di * wgti) * (1 - Di) * wgti) / sum((1 - Di) * wgti)
    out.AUC <- c(out.AUC, AUC)
  }
  c(out.AUC, out.yy, out.pp, out.FPR, out.TPR, out.PPV, out.NPV)
}

# Adaptive LASSO with BIC selection
Est.ALASSO.GLM <- function(data, Wi = NULL, rtn = "EST", nopen.ind = NULL, regularize = TRUE) {
  data <- as.matrix(data); y <- data[, 1]; x <- data[, -1, drop = FALSE]
  nn <- length(y); if (is.null(Wi)) Wi <- rep(1, nn); pp <- ncol(x)
  
  if (regularize) {
    lam.ridge <- pp / nn
    bini <- as.vector(coef(glmnet(x, y, weights = Wi, alpha = 0, standardize = FALSE,
                                  lambda = lam.ridge, family = "binomial")))
    w.b <- 1 / abs(bini[-1]); x.t <- x / VTM(w.b, nrow(x))
    tmpfit <- glmpath(x.t, y, nopenalty.subset = nopen.ind, family = binomial,
                      weight = Wi, standardize = FALSE, min.lambda = 0, lambda2 = lam.ridge)
    lam.all <- seq(min(tmpfit$lambda), max(tmpfit$lambda), length = 500)
    b.all <- predict(tmpfit, s = lam.all, type = "coefficients", mode = "lambda")
    b.all <- b.all / VTM(c(1, w.b), nrow(b.all))
    df.all <- apply(b.all[, -1, drop = FALSE] != 0, 1, sum)
    BIC.lam <- -2 * apply(predict(tmpfit, newx = x.t, newy = y, s = lam.all,
                                  type = "loglik", mode = "lambda"), 2, sum) +
      min(sum(y)^0.1, log(sum(y))) * df.all
    m.opt <- which.min(BIC.lam)
    bhat <- b.all[m.opt, ]; lamhat <- lam.all[m.opt]
  } else {
    bhat <- bini <- coef(glm(y ~ x, family = binomial, weights = Wi))
    lamhat <- 0; lam.all <- BIC.lam <- b.all <- NULL
  }
  out <- c("b" = bhat, "bini" = bini, "lamhat" = lamhat,
           "lam.all" = lam.all, "BIC.lam" = BIC.lam, "b.all" = b.all)
  if (rtn == "EST") return(out)
}

#### =============================== ####
#### 2. DATA SIMULATION              ####
#### =============================== ####

generate_data <- function(n, p, seed = 123) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- paste0("V", 1:p)
  X_df <- as.data.frame(X)
  X_df$period12 <- sample(0:1, n, replace = TRUE)
  X_df$inpat <- sample(0:1, n, replace = TRUE)
  X_df$U099_Count <- rpois(n, lambda = 3)
  
  logit_u099 <- 0.5 * X_df$V1 - 0.3 * X_df$V2 + 0.4 * X_df$period12 - 0.4 * X_df$inpat
  X_df$u099.flag <- rbinom(n, 1, g.logit(logit_u099))
  
  logit_Y <- 2 * X_df$u099.flag + 2 * X_df$V3 + 2 * X_df$U099_Count
  X_df$Y <- rbinom(n, 1, g.logit(logit_Y))
  ind.miss <- which(rbinom(n, 1, 0.6) == 1)
  X_df$Y[ind.miss] <- NA
  return(X_df)
}

n <- 1000; p <- 50
dat <- generate_data(n, p)

#### =============================== ####
#### 3. XGBOOST TRAINING             ####
#### =============================== ####

params <- list(
  booster = "gbtree", objective = "binary:logistic",
  eta = 0.1, max_depth = 6, eval_metric = "auc"
)

xgb.model <- list(); xgb.pred0 <- NULL
for (hospital.setting in c("allpat", "inpat", "outpat")) {
  for (period.setting in c("cotrain.ind", "12", "3")) {
    feature.sel.keep <- paste0("V", 1:30)
    if (period.setting == "cotrain.ind") feature.sel.keep <- c("period12", feature.sel.keep)
    if (hospital.setting == "allpat") feature.sel.keep <- c("inpat", feature.sel.keep)
    
    idx <- which(
      (period.setting == "cotrain.ind" |
         (period.setting == "12" & dat$period12 == 1) |
         (period.setting == "3" & dat$period12 == 0)) &
        (hospital.setting == "allpat" |
           (hospital.setting == "inpat" & dat$inpat == 1) |
           (hospital.setting == "outpat" & dat$inpat == 0))
    )
    
    dtrain <- xgb.DMatrix(data = data.matrix(dat[idx, feature.sel.keep]),
                          label = dat$u099.flag[idx])
    dtest <- xgb.DMatrix(data = data.matrix(dat[, feature.sel.keep]))
    
    model <- xgb.train(params = params, data = dtrain, nrounds = 50, verbose = 0)
    model_name <- paste0(period.setting, "_", hospital.setting)
    xgb.model[[model_name]] <- model
    xgb.pred0 <- cbind(xgb.pred0, predict(model, dtest))
  }
}

colnames(xgb.pred0) <- names(xgb.model)
xgb.pred <- data.frame(
  u099.flag = dat$u099.flag,
  period12 = dat$period12,
  inpat = dat$inpat,
  xgb.pred0
)

#### =============================== ####
#### 4. LOOCV WITH ADAPTIVE LASSO    ####
#### =============================== ####

# Prepare labeled data for LOOCV
ind.miss <- which(!is.na(dat$Y))
Y <- dat$Y[ind.miss]
dat$model.score <- xgb.pred$cotrain.ind_outpat
S <- dat[ind.miss, c("u099.flag", "period12", "inpat", "U099_Count", "model.score")]

ssl.pred.leave1out.fun <- function(Y, S, alg = "alasso") {
  data <- cbind(Y, S)
  ssl.pred <- rep(NA, length(Y))
  
  if (alg == "alasso") {
    for (k in 1:length(Y)) {
      dat.t <- data[-k, ]
      dat.v <- data[k, ]
      beta.t <- Est.ALASSO.GLM(dat.t, regularize = TRUE)
      beta.t <- beta.t[1:ncol(data)]
      x.v <- c(1, as.numeric(dat.v[-1]))
      ssl.pred[k] <- g.logit(sum(x.v * beta.t))
    }
  }
  return(ssl.pred)
}

preds <- ssl.pred.leave1out.fun(Y, S, alg = "alasso")

# Evaluate ROC
ROC.Est.FUN(Y, preds, yy0 = 0.5)
