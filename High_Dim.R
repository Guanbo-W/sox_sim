library(sox)
library(glmnet)

n <- 400
p <- 20
# No. of simulation replications
n.sim <- 2

set.seed(20230218)

# generate main effects
# x_1, ..., x_p
cov.mat <- matrix(NA, nrow = p, ncol = p)
auto.reg <- 0.3

for (i in 1:p) {
  for (j in i:p) {
    cov.mat[i, j] <- cov.mat[j, i] <- auto.reg^abs(i - j)
  }
}

cov.chol <- chol(cov.mat)

p.full <- p + p * (p - 1)/2

# generate interaction terms
# x_1 * x_2, x_1 * x_3, x_1 * x_p, x_2 * x_3, ..., x_p-1 * x_p
two.way.intx <- as.matrix(expand.grid(seq(p), seq(p)))
two.way.intx <- two.way.intx[, 2:1]
two.way.intx <- two.way.intx[two.way.intx[, 1] < two.way.intx[, 2], ]

# generate true beta
beta.true <- rep(0, p.full)

## nonzero interaction terms
interactions <- list(
  c(1, 2),
  c(1, 3),
  c(1, 7),
  c(1, 8),
  c(1, 9),
  c(4, 5),
  c(4, 6),
  c(7, 8),
  c(7, 9)
)
nz.main.id <- sort(unique(unlist(interactions)))
beta.true[nz.main.id] <- 0.4

nz.intx.id <- lapply(interactions, function(x) {
  which(colSums(apply(two.way.intx, MARGIN = 1, FUN = "==", x)) == 2)
})
nz.intx.id <- unlist(nz.intx.id) + p
beta.true[nz.intx.id] <- 0.3

# sox
grp <- matrix(0L, nrow = p.full, ncol = p.full)
for (i in (p+1):p.full) {
  grp[i, two.way.intx[i - p, ]] <- 1L
}

grp.var <- diag(p.full)
eta.g <- rep(1, p.full)

lam.seq <- 10^seq(-0.5, -2, -0.05)
lam.seq.db <- 10^seq(0, -4, -0.1)

precision <- recall <- precision.db <- recall.db <- numeric(n.sim)
precision.glmnet <- recall.glmnet <- precision.glmnet.db <- recall.glmnet.db <- numeric(n.sim)

nonzero.esti <- nonzero.esti.db <- nonzero.esti.glmnet <- nonzero.esti.glmnet.db <- vector(mode = "list", length = n.sim)

sox.esti <- sox.db.esti <- glmnet.esti <- glmnet.db.esti <- matrix(NA, nrow = p.full, ncol = n.sim)

for (s in 1:n.sim) {
  
  x <- matrix(NA, nrow = n, ncol = p.full)
  x[, seq(p)] <- mvnfast::rmvn(n, rep(0, p),
                               sigma = cov.chol,
                               ncores = 8, isChol = T)
  
  for (i in seq(nrow(two.way.intx))) {
    x[, p + i] <- x[, two.way.intx[i, 1]] * x[, two.way.intx[i, 2]]
  }
  
  # simulate survival data
  T.max <- 100
  dat <- coxed::sim.survdata(T = T.max, X = x, censor = 0.3, beta = beta.true)
  
  cv <- sox_cv(x = x,
               ID = seq(n),
               time = rep(0, n),
               time2 = dat$data$y,
               event = dat$data$failed,
               lambda = lam.seq,
               group = grp,
               group_variable = grp.var,
               penalty_weights = eta.g,
               nfolds = 10,
               tol = 1e-4,
               maxit = 1e3,
               verbose = FALSE)
  plot_sox_cv(cv)
  
  plot_sox_sp(sox_obj = cv,
              plot_min = TRUE,
              plot_1se = TRUE,
              lwd = 1)
  
  nz.true <- c(nz.main.id, nz.intx.id)
  esti <- cv$sox.fit$estimates[, paste(cv$lambda.1se)]
  sox.esti[, s] <- esti
  
  nz.esti <- which(esti != 0)
  tp <- length(intersect(nz.true, nz.esti))
  fn <- length(setdiff(nz.true, nz.esti))
  fp <- length(setdiff(nz.esti, nz.true))
  precision[s] <- tp/(tp+fp)
  recall[s] <- tp/(tp+fn)
  
  nonzero.esti[[s]] <- nz.esti
  
  cat("precision:", precision[s], "\n")
  cat("recall:", recall[s], "\n")
  cat("nonzero:", nz.esti, "\n")
  
  
  ### debias with inverse-norm group weights
  eta.g.db <- 1/sqrt(pmax(abs(esti), .Machine$double.eps))
  
  cv.db <- sox_cv(x = x,
                  ID = seq(n),
                  time = rep(0, n),
                  time2 = dat$data$y,
                  event = dat$data$failed,
                  lambda = lam.seq.db,
                  group = grp,
                  group_variable = grp.var,
                  penalty_weights = eta.g.db,
                  nfolds = 10,
                  tol = 1e-4,
                  maxit = 1e3,
                  verbose = FALSE)
  plot_sox_cv(cv.db)
  
  plot_sox_sp(sox_obj = cv.db,
              plot_min = TRUE,
              plot_1se = TRUE,
              lwd = 1)
  
  esti.db <- cv.db$sox.fit$estimates[, paste(cv.db$lambda.1se)]
  sox.db.esti[, s] <- esti.db
  
  nz.esti.db <- which(esti.db != 0)
  tp.db <- length(intersect(nz.true, nz.esti.db))
  fn.db <- length(setdiff(nz.true, nz.esti.db))
  fp.db <- length(setdiff(nz.esti.db, nz.true))
  precision.db[s] <- tp.db/(tp.db+fp.db)
  recall.db[s] <- tp.db/(tp.db+fn.db)
  
  nonzero.esti.db[[s]] <- nz.esti.db
  
  cat("precision with debias:", precision.db[s], "\n")
  cat("recall with debias:", recall.db[s], "\n")
  cat("nonzero with debias:", nz.esti.db, "\n")
  
  # glmnet cox model
  message("glmnet")
  yy <- dat$data[, c("y", "failed")]
  yy <- sapply(yy, as.numeric)
  yy <- as.matrix(yy)
  colnames(yy) <- c("time", "status")
  
  # without debiasing
  glmnet.fit <- cv.glmnet(x = x,
                          y = yy,
                          family = "cox",
                          nfolds = 10,
                          type.measure = "deviance")
  plot(glmnet.fit)
  
  nz.esti.glmnet <- coef(glmnet.fit, s = "lambda.1se")@i + 1L
  glmnet.esti[, s] <- as.vector(coef(glmnet.fit, s = "lambda.1se"))
  
  tp.glmnet <- length(intersect(nz.true, nz.esti.glmnet))
  fn.glmnet <- length(setdiff(nz.true, nz.esti.glmnet))
  fp.glmnet <- length(setdiff(nz.esti.glmnet, nz.true))
  precision.glmnet[s] <- tp.glmnet/(tp.glmnet+fp.glmnet)
  recall.glmnet[s] <- tp.glmnet/(tp.glmnet+fn.glmnet)
  
  nonzero.esti.glmnet[[s]] <- nz.esti.glmnet
  
  cat("precision glmnet:", precision.glmnet[s], "\n")
  cat("recall glmnet:", recall.glmnet[s], "\n")
  cat("nonzero glmnet:", nz.esti.glmnet, "\n")
  
  # with debiasing
  message("glmnet with debias")
  
  penalty.factor <- 1/sqrt(pmax(abs(coef(glmnet.fit, s = "lambda.1se")), .Machine$double.eps))
  
  glmnet.fit.db <- cv.glmnet(x = x,
                             y = yy,
                             family = "cox",
                             penalty.factor = penalty.factor,
                             nfolds = 10,
                             type.measure = "deviance")
  plot(glmnet.fit.db)
  
  nz.esti.glmnet.db <- coef(glmnet.fit.db, s = "lambda.1se")@i + 1L
  glmnet.db.esti[, s] <- as.vector(coef(glmnet.fit.db, s = "lambda.1se"))
  
  tp.glmnet.db <- length(intersect(nz.true, nz.esti.glmnet.db))
  fn.glmnet.db <- length(setdiff(nz.true, nz.esti.glmnet.db))
  fp.glmnet.db <- length(setdiff(nz.esti.glmnet.db, nz.true))
  precision.glmnet.db[s] <- tp.glmnet.db/(tp.glmnet.db+fp.glmnet.db)
  recall.glmnet.db[s] <- tp.glmnet.db/(tp.glmnet.db+fn.glmnet.db)
  
  nonzero.esti.glmnet.db[[s]] <- nz.esti.glmnet.db
  
  cat("precision glmnet with debias:", precision.glmnet.db[s], "\n")
  cat("recall glmnet with debias:", recall.glmnet.db[s], "\n")
  cat("nonzero glmnet with debias:", nz.esti.glmnet.db, "\n")
}

## Results
# precision and recall
mean(precision, na.rm = TRUE); mean(recall, na.rm = TRUE)
sd(precision, na.rm = TRUE)/sqrt(20); sd(recall, na.rm = TRUE)/sqrt(20)

mean(precision.db, na.rm = TRUE); mean(recall.db, na.rm = TRUE)
sd(precision.db, na.rm = TRUE)/sqrt(20); sd(recall.db, na.rm = TRUE)/sqrt(20)

mean(precision.glmnet, na.rm = TRUE); mean(recall.glmnet, na.rm = TRUE)
sd(precision.glmnet, na.rm = TRUE)/sqrt(20); sd(recall.glmnet, na.rm = TRUE)/sqrt(20)

mean(precision.glmnet.db, na.rm = TRUE); mean(recall.glmnet.db, na.rm = TRUE)
sd(precision.glmnet.db, na.rm = TRUE)/sqrt(20); sd(recall.glmnet.db, na.rm = TRUE)/sqrt(20)

# MSE
mean(apply(sox.esti - beta.true, MARGIN = 2, FUN = function(x) {sum(x^2)}))
mean(apply(sox.db.esti - beta.true, MARGIN = 2, FUN = function(x) {sum(x^2)}))
mean(apply(glmnet.esti - beta.true, MARGIN = 2, FUN = function(x) {sum(x^2)}))
mean(apply(glmnet.db.esti - beta.true, MARGIN = 2, FUN = function(x) {sum(x^2)}))

# Violations
violations.glmnet <- violations.glmnet.db <- 0
for (i in 1:n.sim) {
  nz.tmp <- nonzero.esti.glmnet[[i]]
  intx.tmp <- nz.tmp[nz.tmp > p]
  main.tmp <- nz.tmp[nz.tmp <= p]
  
  for (j in intx.tmp) {
    mains <- two.way.intx[j - p, ]
    if (sum(mains %in% main.tmp) < 2) {
      violations.glmnet <- violations.glmnet + 1
    }
  }
  
  nz.tmp.db <- nonzero.esti.glmnet.db[[i]]
  intx.tmp.db <- nz.tmp.db[nz.tmp.db > p]
  main.tmp.db <- nz.tmp.db[nz.tmp.db <= p]
  
  for (j in intx.tmp.db) {
    mains.db <- two.way.intx[j - p, ]
    if (sum(mains.db %in% main.tmp.db) < 2) {
      violations.glmnet.db <- violations.glmnet.db + 1
    }
  }
}
violations.glmnet; violations.glmnet.db
