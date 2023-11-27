# sox simulation: generate data with time-dependent covariates for the high-dimensional case

# set m=4: each subject has 4 time points
# generate one variable for one subject: 
# first generate 2 numbers~N(0,1), 
# repeat each number either 2 or 3 times, concatenate these repeated values together, 
# resulting in a single vector with a length between 4 and 6
# take the first 4 values,
# so that each subject has 4 values corresponds to time 0, 1, 2, 3
gen_con=function(m){
  X=rnorm(2)
  XX=NULL
  for (i in 1:length(X)) {
    if (length(XX)<m){
      X.rep=rep(X[i],round(runif(1,2,3),0))
      XX=c(XX,X.rep)
    }
  }
  return(XX[1:m])
}

# generate p covariates for one subject
gen_X=function(m,p){
  X=matrix(NA,m,p)
  for (j in seq(p)) {
    X[,j]=gen_con(m)
  }
  return(X)
}
# generate p covariates for all subject
gen_X_n=function(m,n,p){
  Xn=NULL
  for (i in 1:n) {
    X=gen_X(m,p)
    Xn=rbind(Xn,X)
  }
  return(Xn)
}

n <- 400
p <- 20
m <- 4
# generate interactions
######### should be the same as your way to generate interaction
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

nz.true <- c(nz.main.id, nz.intx.id)

## if you want to check the model fit: fit.original = coxph(Surv(Start, Stop, Event) ~ . ,data[,-c(1,3)])
## the values from the generated data on line 50 corresponds to the sox function:
## time = Start,
## time2 = Stop,
## event = Event,

# all the metrics see https://github.com/Guanbo-W/sox_sim/blob/main/Low_Dim.R line 148-end



# sox
library(sox)
grp <- matrix(0L, nrow = p.full, ncol = p.full)
for (i in (p+1):p.full) {
  grp[i, two.way.intx[i - p, ]] <- 1L
}

grp.var <- diag(p.full)
eta.g <- rep(1, p.full)

nlam <- 30
lam.seq.db <- 10^seq(-1, -4, length.out = nlam)
# if (p == 40) {
lam.seq <- 10^seq(-0.5, -2, length.out = nlam)
# } else if (p <= 30) {
# lam.seq <- 10^seq(-0.5, -2, -0.05)
# }

# No. of simulation replications
n.sim <- 20

# precision <- recall <- precision.db <- recall.db <- numeric(n.sim)
# precision.glmnet <- recall.glmnet <- precision.glmnet.db <- recall.glmnet.db <- numeric(n.sim)


# sox.esti <- sox.db.esti <- glmnet.esti <- glmnet.db.esti <- matrix(NA, nrow = p.full, ncol = n.sim)

data.list <- vector(mode = "list", length = n.sim)

t.sox <- t.sox.db <- t.glmnet <- t.glmnet.db <- vector(mode = "list", length = n.sim)

cvs.sox <- cvs.sox.db <- cvs.glmnet <- cvs.glmnet.db <- vector(mode = "list", length = n.sim)

for (s in 1:n.sim) {
  cat("Replication: ", s, "\n")
  # covariates= the covariate matrix (n * (p+(p,2)) )
  x <- matrix(NA, nrow = n * m, ncol = p.full)
  x[, seq(p)] <- gen_X_n(m, n, p)
  
  for (i in seq(nrow(two.way.intx))) {
    x[, p + i] <- x[, two.way.intx[i, 1]] * x[, two.way.intx[i, 2]]
  }
  ######### the end of generating interactions
  
  data <- PermAlgo::permalgorithm(n, m, x, 
                                  betas = beta.true,
                                  groupByD = FALSE)
  
  data.list[[s]] <- data
  
  x <- as.matrix(data[, -(1:5)])
  
  message("sox")
  t0 <- proc.time()
  cv <- sox_cv(x = x,
               ID = data$Id,
               time = data$Start,
               time2 = data$Stop,
               event = data$Event,
               penalty = "overlapping",
               lambda = lam.seq,
               group = grp,
               group_variable = grp.var,
               penalty_weights = eta.g,
               nfolds = 10,
               tol = 1e-4,
               maxit = 1e3,
               verbose = FALSE)
  t.sox[[s]] <- proc.time() - t0
  cvs.sox[[s]] <- cv
  
  plot_sox_cv(cv)
  plot_sox_sp(sox_obj = cv,
              plot_min = TRUE,
              plot_1se = TRUE,
              lwd = 1)
  
  esti <- cv$sox.fit$estimates[, paste(cv$lambda.1se)]
  # sox.esti[, s] <- esti
  # 
  # nz.esti <- which(esti != 0)
  # tp <- length(intersect(nz.true, nz.esti))
  # fn <- length(setdiff(nz.true, nz.esti))
  # fp <- length(setdiff(nz.esti, nz.true))
  # precision[s] <- tp/(tp+fp)
  # recall[s] <- tp/(tp+fn)
  # 
  # nonzero.sox[[s]] <- nz.esti
  # 
  # cat("precision:", precision[s], "\n")
  # cat("recall:", recall[s], "\n")
  # cat("nonzero:", nz.esti, "\n")
  
  message("sox with debias")
  ### debias with inverse-norm group weights
  eta.g.db <- 1/sqrt(pmax(abs(esti), .Machine$double.eps))
  
  t0 <- proc.time()
  cv.db <- sox_cv(x = x,
                  ID = data$Id,
                  time = data$Start,
                  time2 = data$Stop,
                  event = data$Event,
                  penalty = "overlapping",
                  lambda = lam.seq.db,
                  group = grp,
                  group_variable = grp.var,
                  penalty_weights = eta.g.db,
                  nfolds = 10,
                  tol = 1e-4,
                  maxit = 1e3,
                  verbose = FALSE)
  t.sox.db[[s]] <- proc.time() - t0
  cvs.sox.db[[s]] <- cv.db
  
  plot_sox_cv(cv.db)
  plot_sox_sp(sox_obj = cv.db,
              plot_min = TRUE,
              plot_1se = TRUE,
              lwd = 1)
  
  # esti.db <- cv.db$sox.fit$estimates[, paste(cv.db$lambda.1se)]
  # sox.db.esti[, s] <- esti.db
  # 
  # nz.esti.db <- which(esti.db != 0)
  # tp.db <- length(intersect(nz.true, nz.esti.db))
  # fn.db <- length(setdiff(nz.true, nz.esti.db))
  # fp.db <- length(setdiff(nz.esti.db, nz.true))
  # precision.db[s] <- tp.db/(tp.db+fp.db)
  # recall.db[s] <- tp.db/(tp.db+fn.db)
  # 
  # nonzero.sox.db[[s]] <- nz.esti.db
  # 
  # cat("precision with debias:", precision.db[s], "\n")
  # cat("recall with debias:", recall.db[s], "\n")
  # cat("nonzero with debias:", nz.esti.db, "\n")
  
  # glmnet cox model
  message("glmnet")
  yy <- Surv(data$Start, data$Stop, data$Event)
  
  # without debiasing
  t0 <- proc.time()
  glmnet.fit <- cv.glmnet(x = x,
                          y = yy,
                          family = "cox",
                          nfolds = 10,
                          type.measure = "deviance")
  t.glmnet[[s]] <- proc.time() - t0
  cvs.glmnet[[s]] <- glmnet.fit
  
  plot(glmnet.fit)
  
  # nz.esti.glmnet <- coef(glmnet.fit, s = "lambda.1se")@i + 1L
  # glmnet.esti[, s] <- as.vector(coef(glmnet.fit, s = "lambda.1se"))
  # 
  # tp.glmnet <- length(intersect(nz.true, nz.esti.glmnet))
  # fn.glmnet <- length(setdiff(nz.true, nz.esti.glmnet))
  # fp.glmnet <- length(setdiff(nz.esti.glmnet, nz.true))
  # precision.glmnet[s] <- tp.glmnet/(tp.glmnet+fp.glmnet)
  # recall.glmnet[s] <- tp.glmnet/(tp.glmnet+fn.glmnet)
  # 
  # nonzero.glmnet[[s]] <- nz.esti.glmnet
  # 
  # cat("precision glmnet:", precision.glmnet[s], "\n")
  # cat("recall glmnet:", recall.glmnet[s], "\n")
  # cat("nonzero glmnet:", nz.esti.glmnet, "\n")
  
  # with debiasing
  message("glmnet with debias")
  
  penalty.factor <- 1/sqrt(pmax(abs(coef(glmnet.fit, s = "lambda.1se")), .Machine$double.eps))
  
  t0 <- proc.time()
  glmnet.fit.db <- cv.glmnet(x = x,
                             y = yy,
                             family = "cox",
                             penalty.factor = penalty.factor,
                             nfolds = 10,
                             type.measure = "deviance")
  t.glmnet.db[[s]] <- proc.time() - t0
  cvs.glmnet.db[[s]] <- glmnet.fit.db
  
  plot(glmnet.fit.db)
  
  # nz.esti.glmnet.db <- coef(glmnet.fit.db, s = "lambda.1se")@i + 1L
  # glmnet.db.esti[, s] <- as.vector(coef(glmnet.fit.db, s = "lambda.1se"))
  # 
  # tp.glmnet.db <- length(intersect(nz.true, nz.esti.glmnet.db))
  # fn.glmnet.db <- length(setdiff(nz.true, nz.esti.glmnet.db))
  # fp.glmnet.db <- length(setdiff(nz.esti.glmnet.db, nz.true))
  # precision.glmnet.db[s] <- tp.glmnet.db/(tp.glmnet.db+fp.glmnet.db)
  # recall.glmnet.db[s] <- tp.glmnet.db/(tp.glmnet.db+fn.glmnet.db)
  # 
  # nonzero.glmnet.db[[s]] <- nz.esti.glmnet.db
  # 
  # cat("precision glmnet with debias:", precision.glmnet.db[s], "\n")
  # cat("recall glmnet with debias:", recall.glmnet.db[s], "\n")
  # cat("nonzero glmnet with debias:", nz.esti.glmnet.db, "\n")
}

# # mean(precision, na.rm = TRUE); mean(recall, na.rm = TRUE)
# # sd(precision, na.rm = TRUE)/sqrt(20); sd(recall, na.rm = TRUE)/sqrt(20)
# # 
# # mean(precision.db, na.rm = TRUE); mean(recall.db, na.rm = TRUE)
# # sd(precision.db, na.rm = TRUE)/sqrt(20); sd(recall.db, na.rm = TRUE)/sqrt(20)
# # 
# # mean(precision.glmnet, na.rm = TRUE); mean(recall.glmnet, na.rm = TRUE)
# # sd(precision.glmnet, na.rm = TRUE)/sqrt(20); sd(recall.glmnet, na.rm = TRUE)/sqrt(20)
# # 
# # mean(precision.glmnet.db, na.rm = TRUE); mean(recall.glmnet.db, na.rm = TRUE)
# # sd(precision.glmnet.db, na.rm = TRUE)/sqrt(20); sd(recall.glmnet.db, na.rm = TRUE)/sqrt(20)
# 
# # MSE
# mean(apply(sox.esti - beta.true, MARGIN = 2, FUN = function(x) {sum(x^2)}))
# mean(apply(sox.db.esti - beta.true, MARGIN = 2, FUN = function(x) {sum(x^2)}))
# mean(apply(glmnet.esti - beta.true, MARGIN = 2, FUN = function(x) {sum(x^2)}))
# mean(apply(glmnet.db.esti - beta.true, MARGIN = 2, FUN = function(x) {sum(x^2)}))
# 
# violations.glmnet <- violations.glmnet.db <- 0
# for (i in 1:n.sim) {
#   nz.tmp <- nonzero.glmnet[[i]]
#   intx.tmp <- nz.tmp[nz.tmp > p]
#   main.tmp <- nz.tmp[nz.tmp <= p]
#   
#   for (j in intx.tmp) {
#     mains <- two.way.intx[j - p, ]
#     if (sum(mains %in% main.tmp) < 2) {
#       violations.glmnet <- violations.glmnet + 1
#     }
#   }
#   
#   nz.tmp.db <- nonzero.glmnet.db[[i]]
#   intx.tmp.db <- nz.tmp.db[nz.tmp.db > p]
#   main.tmp.db <- nz.tmp.db[nz.tmp.db <= p]
#   
#   for (j in intx.tmp.db) {
#     mains.db <- two.way.intx[j - p, ]
#     if (sum(mains.db %in% main.tmp.db) < 2) {
#       violations.glmnet.db <- violations.glmnet.db + 1
#     }
#   }
# }
# violations.glmnet; violations.glmnet.db
# 
# beepr::beep(8)
# 
file.name <- paste0("/Users/YLIAN/Desktop/Research/Collab/GuanboWang/SVS_TD_Cox/revision/n", n, "p", p, ".Rdata")
# save.image(file.name)

# non zeros
nz.true <- c(nz.main.id, nz.intx.id)

nz.sox <- lapply(cvs.sox,
                 FUN = function(x) {
                   which(x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] != 0)
                 }
)
nz.sox.db <- lapply(cvs.sox.db,
                    FUN = function(x) {
                      which(x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] != 0)
                    }
)
nz.glmnet <- lapply(cvs.glmnet,
                    FUN = function(x){
                      coef(x, s = "lambda.1se")@i + 1L
                    }
)
nz.glmnet.db <- lapply(cvs.glmnet.db,
                       FUN = function(x){
                         coef(x, s = "lambda.1se")@i + 1L
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



# R1S
r1s.sox <- r1s.sox.db <- r1s.glmnet <- r1s.glmnet.db <- 0
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
}
1 - (r1s.sox/n.sim/(p * (p - 1) /2))
1 - (r1s.sox.db/n.sim/(p * (p - 1) /2))
1 - (r1s.glmnet/n.sim/(p * (p - 1) /2))
1 - (r1s.glmnet.db/n.sim/(p * (p - 1) /2))

# RCI
rci.sox <- rci.sox.db <- rci.glmnet <- rci.glmnet.db <- numeric(n.sim)
for (s in 1:n.sim) {
  data <- data.list[[s]]
  x <- as.matrix(data[, -(1:5)])
  time = data$Start
  time2 = data$Stop
  event = data$Event
  
  x.sox <- as.matrix(x[, nz.sox[[s]]])
  x.sox.db <- x[, nz.sox.db[[s]]]
  x.glmnet <- x[, nz.glmnet[[s]]]
  x.glmnet.db <- x[, nz.glmnet.db[[s]]]
  
  cox.sox <- survival::coxph(survival::Surv(time, time2, event) ~ x.sox)
  rci.sox[s] <- cox.sox$concordance["concordance"]
  
  cox.sox.db <- survival::coxph(survival::Surv(time, time2, event) ~ x.sox.db)
  rci.sox.db[s] <- cox.sox.db$concordance["concordance"]
  
  cox.glmnet <- survival::coxph(survival::Surv(time, time2, event) ~ x.glmnet)
  rci.glmnet[s] <- cox.glmnet$concordance["concordance"]
  
  cox.glmnet.db <- survival::coxph(survival::Surv(time, time2, event) ~ x.glmnet.db)
  rci.glmnet.db[s] <- cox.glmnet.db$concordance["concordance"]
}

mean(rci.sox)
mean(rci.sox.db)
mean(rci.glmnet)
mean(rci.glmnet.db)

# MSE
mean(sapply(cvs.sox,
            FUN = function(x) {
              mean((x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] - beta.true)^2)
            }))
mean(sapply(cvs.sox.db,
            FUN = function(x) {
              mean((x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] - beta.true)^2)
            }))
mean(sapply(cvs.glmnet,
            FUN = function(x) {
              mean((coef(x, s = "lambda.1se") - beta.true)^2)
            }))
mean(sapply(cvs.glmnet.db,
            FUN = function(x) {
              mean((coef(x, s = "lambda.1se") - beta.true)^2)
            }))

# CV error
cve.sox <- cve.sox.db <- cve.glmnet <- cve.glmnet.db <- numeric(n.sim)
for (s in 1:n.sim) {
  cve.sox[s] <- cvs.sox[[s]]$cvm[which(cvs.sox[[s]]$lambdas == cvs.sox[[s]]$lambda.1se)]
  cve.sox.db[s] <- cvs.sox.db[[s]]$cvm[which(cvs.sox.db[[s]]$lambdas == cvs.sox.db[[s]]$lambda.1se)]
  cve.glmnet[s] <- cvs.glmnet[[s]]$cvm[which(cvs.glmnet[[s]]$lambda == cvs.glmnet[[s]]$lambda.1se)]
  cve.glmnet.db[s] <- cvs.glmnet.db[[s]]$cvm[which(cvs.glmnet.db[[s]]$lambda == cvs.glmnet.db[[s]]$lambda.1se)]
}
mean(cve.sox)
mean(cve.sox.db)
mean(cve.glmnet)
mean(cve.glmnet.db)

# Timing
t.mean.sox <- 0
t.mean.sox.db <- 0
for (s in 1:n.sim) {
  t.mean.sox <- t.mean.sox + t.sox[[s]][3]
  t.mean.sox.db <- t.mean.sox.db + t.sox.db[[s]][3]
}
t.mean.sox/n.sim
t.mean.sox.db/n.sim



