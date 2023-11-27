setwd("~/Desktop/Research/Collab/GuanboWang/SVS_TD_Cox/revision")

library(sox)
library(glmnet)
library(SGL)

n <- 800
p <- 20

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


# library(survival)
# cox.test <- coxph(Surv(dat$data$y, dat$data$failed) ~ x)
# cox.test$coefficients

if (p == 40) {
  lam.seq <- 10^seq(-0.5, -1.5, -0.05)
} else if (p <= 30) {
  lam.seq <- 10^seq(-0.5, -2, -0.05)
  lam.seq.db <- 10^seq(0, -4, -0.1)
}

# No. of simulation replications
n.sim <- 20

data.list <- vector(mode = "list", length = n.sim)

t.sox <- t.sox.db <- t.glmnet <- t.glmnet.db <- vector(mode = "list", length = n.sim)

cvs.sox <- cvs.sox.db <- cvs.glmnet <- cvs.glmnet.db <- cvs.sgl <- vector(mode = "list", length = n.sim)

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
  
  data.list[[s]] <- dat
  
  message("sox")
  
  cv <- sox_cv(x = x,
               ID = seq(n),
               time = rep(0, n),
               time2 = dat$data$y,
               event = dat$data$failed,
               penalty = "overlapping",
               lambda = lam.seq,
               group = grp,
               group_variable = grp.var,
               penalty_weights = eta.g,
               nfolds = 10,
               tol = 1e-4,
               maxit = 1e3,
               verbose = FALSE)
  cvs.sox[[s]] <- cv
  
  plot_sox_cv(cv)
  plot_sox_sp(sox_obj = cv,
              plot_min = TRUE,
              plot_1se = TRUE,
              lwd = 1)  
  
  
  ### debias with inverse-norm group weights
  message("sox with debias")
  
  eta.g.db <- 1/sqrt(pmax(abs(cv$sox.fit$estimates[, paste(cv$lambda.min)]),
                          .Machine$double.eps))
  
  cv.db <- sox_cv(x = x,
                  ID = seq(n),
                  time = rep(0, n),
                  time2 = dat$data$y,
                  event = dat$data$failed,
                  penalty = "overlapping",
                  lambda = lam.seq.db,
                  group = grp,
                  group_variable = grp.var,
                  penalty_weights = eta.g.db,
                  nfolds = 10,
                  tol = 1e-4,
                  maxit = 1e3,
                  verbose = FALSE)
  cvs.sox.db[[s]] <- cv.db
  
  plot_sox_cv(cv.db)
  plot_sox_sp(sox_obj = cv.db,
              plot_min = TRUE,
              plot_1se = TRUE,
              lwd = 1)
  
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
  cvs.glmnet[[s]] <- glmnet.fit
  
  plot(glmnet.fit)
  
  # with debiasing
  message("glmnet with debias")
  
  penalty.factor <- 1/sqrt(pmax(abs(coef(glmnet.fit, s = "lambda.min")), .Machine$double.eps))
  
  glmnet.fit.db <- cv.glmnet(x = x,
                             y = yy,
                             family = "cox",
                             penalty.factor = penalty.factor,
                             nfolds = 10,
                             type.measure = "deviance")
  cvs.glmnet.db[[s]] <- glmnet.fit.db
  
  plot(glmnet.fit.db)
  
  # sgl
  message("SGL")
  data.sgl <- list()
  data.sgl$x <- x
  data.sgl$time <- as.numeric(dat$data$y)
  data.sgl$status <- as.integer(dat$data$failed)
  
  # groups: (1, 2, 1 x 2), (3, 4, 3 x 4), ..., (p-1, p, p-1 x p)
  n.grp.sgl <- p.full - p/2
  intx.in.grps.id <- which(two.way.intx[, 2] - two.way.intx[, 1] == 1 & two.way.intx[, 2]%%2 == 0)
  main.in.grps <- two.way.intx[intx.in.grps.id, ]
  grp.sgl <- integer(p.full)
  grp.sgl[1:p] <- ceiling((1:p)/2)
  grp.sgl[intx.in.grps.id + p] <- ceiling(1:(p/2))
  sgl.grp.id <- p/2
  for (j in 1:p.full) {
    if (grp.sgl[j] == 0) {
      sgl.grp.id <- sgl.grp.id + 1
      grp.sgl[j] <- sgl.grp.id
    }
  }
  if (sgl.grp.id != p.full - p/2 * 2) cat("Wrong grp.sgl!")
  
  cv.sgl <- cvSGL(data.sgl,
                  index = grp.sgl,
                  type = "cox")
  cvs.sgl[[s]] <- cv.sgl
  
  # plot(cv.sgl)
  # abline(v = log(cv.sgl$lambdas[which.min(cv.sgl$lldiff)]), lty = 2)
  
}

file.name <- paste0("/Users/YLIAN/Desktop/Research/Collab/GuanboWang/SVS_TD_Cox/revision/n", n, "p", p, ".Rdata")
# save.image(file.name)

# non zeros
nz.true <- c(nz.main.id, nz.intx.id)

nz.sox <- lapply(cvs.sox,
                 FUN = function(x) {
                   which(x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.min)] != 0)
                 }
)
nz.sox.db <- lapply(cvs.sox.db,
                    FUN = function(x) {
                      which(x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.min)] != 0)
                    }
)
nz.glmnet <- lapply(cvs.glmnet,
                    FUN = function(x){
                      coef(x, s = "lambda.min")@i + 1L
                    }
)
nz.glmnet.db <- lapply(cvs.glmnet.db,
                       FUN = function(x){
                         coef(x, s = "lambda.min")@i + 1L
                       }
)
nz.sgl <- lapply(cvs.sgl,
                 FUN = function(x){
                   which(x$fit$beta[, which.min(x$lldiff)] != 0)
                 }
)

# joint detection rate
mean(sapply(nz.sox, FUN = function(x) {
  length(setdiff(nz.true, x)) == 0
}, simplify = TRUE))
mean(sapply(nz.sox.db, FUN = function(x) {
  length(setdiff(nz.true, x)) == 0
}, simplify = TRUE))
mean(sapply(nz.glmnet, FUN = function(x) {
  length(setdiff(nz.true, x)) == 0
}, simplify = TRUE))
mean(sapply(nz.glmnet.db, FUN = function(x) {
  length(setdiff(nz.true, x)) == 0
}, simplify = TRUE))
mean(sapply(nz.sgl, FUN = function(x) {
  length(setdiff(nz.true, x)) == 0
}, simplify = TRUE))

# missing rate
nz.length <- length(nz.true)
mean(sapply(nz.sox, FUN = function(x) {
  length(setdiff(nz.true, x))/nz.length
}))
mean(sapply(nz.sox.db, FUN = function(x) {
  length(setdiff(nz.true, x))/nz.length
}))
mean(sapply(nz.glmnet, FUN = function(x) {
  length(setdiff(nz.true, x))/nz.length
}))
mean(sapply(nz.glmnet.db, FUN = function(x) {
  length(setdiff(nz.true, x))/nz.length
}))
mean(sapply(nz.sgl, FUN = function(x) {
  length(setdiff(nz.true, x))/nz.length
}))

# false alarm rate
z.length <- p * (p - 1)/2 - nz.length
mean(sapply(nz.sox, FUN = function(x) {
  length(setdiff(x, nz.true))/z.length
}))
mean(sapply(nz.sox.db, FUN = function(x) {
  length(setdiff(x, nz.true))/z.length
}))
mean(sapply(nz.glmnet, FUN = function(x) {
  length(setdiff(x, nz.true))/z.length
}))
mean(sapply(nz.glmnet.db, FUN = function(x) {
  length(setdiff(x, nz.true))/z.length
}))
mean(sapply(nz.sgl, FUN = function(x) {
  length(setdiff(x, nz.true))/z.length
}))


# R1S
r1s.sox <- r1s.sox.db <- r1s.glmnet <- r1s.glmnet.db <- r1s.sgl <- 0
for (i in 1:n.sim) {
  nz.tmp <- nz.sox[[i]]
  intx.tmp <- nz.tmp[nz.tmp > p]
  main.tmp <- nz.tmp[nz.tmp <= p]
  
  for (j in intx.tmp) {
    mains <- two.way.intx[j - p, ]
    if (sum(mains %in% main.tmp) < 2) {
      r1s.sox <- r1s.sox + 1
    }
  }
  
  nz.tmp.db <- nz.sox.db[[i]]
  intx.tmp.db <- nz.tmp.db[nz.tmp.db > p]
  main.tmp.db <- nz.tmp.db[nz.tmp.db <= p]
  
  for (j in intx.tmp.db) {
    mains.db <- two.way.intx[j - p, ]
    if (sum(mains.db %in% main.tmp.db) < 2) {
      r1s.sox.db <- r1s.sox.db + 1
    }
  }
  
  nz.tmp <- nz.glmnet[[i]]
  intx.tmp <- nz.tmp[nz.tmp > p]
  main.tmp <- nz.tmp[nz.tmp <= p]
  
  for (j in intx.tmp) {
    mains <- two.way.intx[j - p, ]
    if (sum(mains %in% main.tmp) < 2) {
      r1s.glmnet <- r1s.glmnet + 1
    }
  }
  
  nz.tmp.db <- nz.glmnet.db[[i]]
  intx.tmp.db <- nz.tmp.db[nz.tmp.db > p]
  main.tmp.db <- nz.tmp.db[nz.tmp.db <= p]
  
  for (j in intx.tmp.db) {
    mains.db <- two.way.intx[j - p, ]
    if (sum(mains.db %in% main.tmp.db) < 2) {
      r1s.glmnet.db <- r1s.glmnet.db + 1
    }
  }
  
  nz.tmp.sgl <- nz.sgl[[i]]
  intx.tmp.sgl <- nz.tmp.sgl[nz.tmp.sgl > p]
  main.tmp.sgl <- nz.tmp.sgl[nz.tmp.sgl <= p]
  
  for (j in intx.tmp.sgl) {
    mains.sgl <- two.way.intx[j - p, ]
    if (sum(mains.sgl %in% main.tmp.sgl) < 2) {
      r1s.sgl <- r1s.sgl + 1
    }
  }
}
1 - (r1s.sox/n.sim/(p * (p - 1) /2))
1 - (r1s.sox.db/n.sim/(p * (p - 1) /2))
1 - (r1s.glmnet/n.sim/(p * (p - 1) /2))
1 - (r1s.glmnet.db/n.sim/(p * (p - 1) /2))
1 - (r1s.sgl/n.sim/(p * (p - 1) /2))

# RCI
rci.sox <- rci.sox.db <- rci.glmnet <- rci.glmnet.db <- rci.sgl <- numeric(n.sim)
for (s in 1:n.sim) {
  data <- data.list[[s]]
  x <- as.matrix(data$xdata)
  time = data$data$y
  event = data$data$failed
  
  x.sox <- as.matrix(x[, nz.sox[[s]]])
  x.sox.db <- x[, nz.sox.db[[s]]]
  x.glmnet <- x[, nz.glmnet[[s]]]
  x.glmnet.db <- x[, nz.glmnet.db[[s]]]
  x.sgl <- x[, nz.sgl[[s]]]
  
  cox.sox <- survival::coxph(survival::Surv(time, event) ~ x.sox)
  rci.sox[s] <- cox.sox$concordance["concordance"]
  
  cox.sox.db <- survival::coxph(survival::Surv(time, event) ~ x.sox.db)
  rci.sox.db[s] <- cox.sox.db$concordance["concordance"]
  
  cox.glmnet <- survival::coxph(survival::Surv(time, event) ~ x.glmnet)
  rci.glmnet[s] <- cox.glmnet$concordance["concordance"]
  
  cox.glmnet.db <- survival::coxph(survival::Surv(time, event) ~ x.glmnet.db)
  rci.glmnet.db[s] <- cox.glmnet.db$concordance["concordance"]
  
  cox.sgl <- if (length(nz.sgl[[s]]) == 0) {
    survival::coxph(survival::Surv(time, event) ~ 1)
  } else {
    survival::coxph(survival::Surv(time, event) ~ x.sgl)
  }
  rci.sgl[s] <- cox.sgl$concordance["concordance"]
}

mean(rci.sox)
mean(rci.sox.db)
mean(rci.glmnet)
mean(rci.glmnet.db)
mean(rci.sgl)

# MSE
mean(sapply(cvs.sox,
            FUN = function(x) {
              mean((x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.min)] - beta.true)^2)
            }))
mean(sapply(cvs.sox.db,
            FUN = function(x) {
              mean((x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.min)] - beta.true)^2)
            }))
mean(sapply(cvs.glmnet,
            FUN = function(x) {
              mean((coef(x, s = "lambda.min") - beta.true)^2)
            }))
mean(sapply(cvs.glmnet.db,
            FUN = function(x) {
              mean((coef(x, s = "lambda.min") - beta.true)^2)
            }))
betas.sgl <- sapply(cvs.sgl,
                    FUN = function(x) {
                      beta.tmp <- c(x$fit$beta[, which.min(x$lldiff)])
                      if (length(beta.tmp) == 0) beta.tmp <- rep(0, p.full)
                      return(beta.tmp)
                      # mean((x$fit$beta[, which.min(x$lldiff)] - beta.true)^2)
                    }, simplify = FALSE)

mean(sapply(betas.sgl,
            FUN = function(x) {
              mean((x - beta.true)^2)
            }))


# CV error
cve.sox <- cve.sox.db <- cve.glmnet <- cve.glmnet.db <- cve.sgl <- numeric(n.sim)
for (s in 1:n.sim) {
  cve.sox[s] <- cvs.sox[[s]]$cvm[which(cvs.sox[[s]]$lambdas == cvs.sox[[s]]$lambda.min)]
  cve.sox.db[s] <- cvs.sox.db[[s]]$cvm[which(cvs.sox.db[[s]]$lambdas == cvs.sox.db[[s]]$lambda.min)]
  cve.glmnet[s] <- cvs.glmnet[[s]]$cvm[which(cvs.glmnet[[s]]$lambda == cvs.glmnet[[s]]$lambda.min)]
  cve.glmnet.db[s] <- cvs.glmnet.db[[s]]$cvm[which(cvs.glmnet.db[[s]]$lambda == cvs.glmnet.db[[s]]$lambda.min)]
}

# data <- train.data <- vector(mode = "list", length = n.sim)
cve.sgl.full <- matrix(nrow = n.sim, ncol = 10)
for (s in 1:n.sim) {
  data <- data.list[[s]]
  beta <- cvs.sgl[[s]]$fit$beta[, which.min(cvs.sgl[[s]]$lldiff)]
  if (is.matrix(beta)) {
    if (ncol(beta) == 0)  beta <- rep(0, p.full)
  } 
  foldid <- cvs.sgl[[s]]$foldid
  
  time <- data$data$y
  event <- data$data$failed
  x <- data$xdata
  for (f in 1:10) {
    train.id <- which(foldid != f)
    time.train <- time[train.id]
    event.train <- event[train.id]
    x.train <- x[train.id, ]
    
    dev1 <- coxnet.deviance(y = Surv(time = time, event = event),
                            x = x,
                            beta = beta)
    dev2 <- coxnet.deviance(y = Surv(time = time.train, event = event.train),
                            x = x.train,
                            beta = beta)
    cve.sgl.full[s, f] <- (dev1 - dev2)/(sum(event) - sum(event.train))
  }
}

mean(cve.sox)
mean(cve.sox.db)
mean(cve.glmnet)
mean(cve.glmnet.db)
mean(cve.sgl.full)


for (s in 1:n.sim) {
  plot_sox_sp(cvs.sox[[s]], plot_min = TRUE)
  # plot(cvs.glmnet[[s]])
  # plot(cvs.sgl[[s]])
  # plot_sox_sp(cvs.sox[[s]])
  # plot(cvs.glmnet[[s]]$glmnet.fit)
  # matplot(x = cvs.sgl[[s]]$lambdas, y = t(cvs.sgl[[s]]$fit$beta), type = "l")
}
