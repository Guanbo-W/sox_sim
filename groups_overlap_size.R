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
p <- 25
m <- 4

beta.true <- c(rep(0, 13), 0.4, 0.4, 0.4, 0, 0, 0.4, 0, 0,0.4, 0, 0, 0)
nz.true <- which(beta.true != 0.0)

# seven different groups
## two groups of size 10, overlaps are of sizes 2, 5, 8
grp282 <- list(1, 2, 3, 4:13, 6:15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)
grp555 <- list(1, 2, 3, 4:13, 9:18, 19, 20, 21, 22, 23, 24, 25)
grp828 <- list(1, 2, 3, 4:13, 12:21, 22, 23, 24, 25)
## two groups of sizes 7, 10, 13, overlap of size 5
grp252 <- list(1, 2, 3, 4, 5, 6, 7:13, 9:15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)
grp858 <- list(1:13, 9:21, 22, 23, 24, 25)
## two groups of sizes 7, overlap of size 5, slide to the left to change sparsity level
grp252l <- list(1, 2, 3, 4, 5, 6:12, 8:14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)
grp252ll <- list(1, 2, 3, 4, 5:11, 7:13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)

graph282 <- graph_builder(grp282)
graph555 <- graph_builder(grp555)
graph828 <- graph_builder(grp828)
graph252 <- graph_builder(grp252)
graph858 <- graph_builder(grp858)
graph252l <- graph_builder(grp252l)
graph252ll <- graph_builder(grp252ll)

# sox
library(sox)

set.seed(20231125)

nlam <- 30
lam.seq <- 10^seq(0, -2, length.out = nlam)

# No. of simulation replications
n.sim <- 20

data.list <- vector(mode = "list", length = n.sim)
t282 <- t555 <- t828 <- t252 <- t858 <- t252l <- t252ll <- vector(mode = "list", length = n.sim)
cvs282 <- cvs555 <- cvs828 <- cvs252 <- cvs858 <- cvs252l <- cvs252ll <- vector(mode = "list", length = n.sim)

for (s in 1:n.sim) {
  cat("Replication: ", s, "\n")
  x <- matrix(NA, nrow = n * m, ncol = p)
  x[, seq(p)] <- gen_X_n(m, n, p)
  
  ######### the end of generating interactions
  
  data <- PermAlgo::permalgorithm(n, m, x, 
                                  betas = beta.true,
                                  groupByD = FALSE)
  
  data.list[[s]] <- data
  
  x <- as.matrix(data[, -(1:5)])
  
  message("Group 282.")
  t0 <- proc.time()
  cv282 <- sox_cv(x = x,
                  ID = data$Id,
                  time = data$Start,
                  time2 = data$Stop,
                  event = data$Event,
                  penalty = "overlapping",
                  lambda = lam.seq,
                  group = graph282$groups,
                  group_variable = graph282$groups_var,
                  penalty_weights = graph282$eta_g,
                  nfolds = 10,
                  tol = 1e-4,
                  maxit = 1e3,
                  verbose = FALSE)
  t282[[s]] <- proc.time() - t0
  cvs282[[s]] <- cv282
  plot_sox_cv(cv282)
  plot_sox_sp(sox_obj = cv282,
              plot_min = TRUE,
              plot_1se = TRUE,
              lwd = 1)
  
  message("Group 555.")
  t0 <- proc.time()
  cv555 <- sox_cv(x = x,
                  ID = data$Id,
                  time = data$Start,
                  time2 = data$Stop,
                  event = data$Event,
                  penalty = "overlapping",
                  lambda = lam.seq,
                  group = graph555$groups,
                  group_variable = graph555$groups_var,
                  penalty_weights = graph555$eta_g,
                  nfolds = 10,
                  tol = 1e-4,
                  maxit = 1e3,
                  verbose = FALSE)
  t555[[s]] <- proc.time() - t0
  cvs555[[s]] <- cv555
  plot_sox_cv(cv555)
  plot_sox_sp(sox_obj = cv555,
              plot_min = TRUE,
              plot_1se = TRUE,
              lwd = 1)
  
  message("Group 828.")
  t0 <- proc.time()
  cv828 <- sox_cv(x = x,
                  ID = data$Id,
                  time = data$Start,
                  time2 = data$Stop,
                  event = data$Event,
                  penalty = "overlapping",
                  lambda = lam.seq,
                  group = graph828$groups,
                  group_variable = graph828$groups_var,
                  penalty_weights = graph828$eta_g,
                  nfolds = 10,
                  tol = 1e-4,
                  maxit = 1e3,
                  verbose = FALSE)
  t828[[s]] <- proc.time() - t0
  cvs828[[s]] <- cv828
  plot_sox_cv(cv828)
  plot_sox_sp(sox_obj = cv828,
              plot_min = TRUE,
              plot_1se = TRUE,
              lwd = 1)
  
  message("Group 252.")
  t0 <- proc.time()
  cv252 <- sox_cv(x = x,
                  ID = data$Id,
                  time = data$Start,
                  time2 = data$Stop,
                  event = data$Event,
                  penalty = "overlapping",
                  lambda = lam.seq,
                  group = graph252$groups,
                  group_variable = graph252$groups_var,
                  penalty_weights = graph252$eta_g,
                  nfolds = 10,
                  tol = 1e-4,
                  maxit = 1e3,
                  verbose = FALSE)
  t252[[s]] <- proc.time() - t0
  cvs252[[s]] <- cv252
  plot_sox_cv(cv252)
  plot_sox_sp(sox_obj = cv252,
              plot_min = TRUE,
              plot_1se = TRUE,
              lwd = 1)
  
  message("Group 858.")
  t0 <- proc.time()
  cv858 <- sox_cv(x = x,
                  ID = data$Id,
                  time = data$Start,
                  time2 = data$Stop,
                  event = data$Event,
                  penalty = "overlapping",
                  lambda = lam.seq,
                  group = graph858$groups,
                  group_variable = graph858$groups_var,
                  penalty_weights = graph858$eta_g,
                  nfolds = 10,
                  tol = 1e-4,
                  maxit = 1e3,
                  verbose = FALSE)
  t858[[s]] <- proc.time() - t0
  cvs858[[s]] <- cv858
  plot_sox_cv(cv858)
  plot_sox_sp(sox_obj = cv858,
              plot_min = TRUE,
              plot_1se = TRUE,
              lwd = 1)
  
  message("Group 252l.")
  t0 <- proc.time()
  cv252l <- sox_cv(x = x,
                   ID = data$Id,
                   time = data$Start,
                   time2 = data$Stop,
                   event = data$Event,
                   penalty = "overlapping",
                   lambda = lam.seq,
                   group = graph252l$groups,
                   group_variable = graph252l$groups_var,
                   penalty_weights = graph252l$eta_g,
                   nfolds = 10,
                   tol = 1e-4,
                   maxit = 1e3,
                   verbose = FALSE)
  t252l[[s]] <- proc.time() - t0
  cvs252l[[s]] <- cv252l
  plot_sox_cv(cv252l)
  plot_sox_sp(sox_obj = cv252l,
              plot_min = TRUE,
              plot_1se = TRUE,
              lwd = 1)
  
  message("Group 252ll.")
  t0 <- proc.time()
  cv252ll <- sox_cv(x = x,
                    ID = data$Id,
                    time = data$Start,
                    time2 = data$Stop,
                    event = data$Event,
                    penalty = "overlapping",
                    lambda = lam.seq,
                    group = graph252ll$groups,
                    group_variable = graph252ll$groups_var,
                    penalty_weights = graph252ll$eta_g,
                    nfolds = 10,
                    tol = 1e-4,
                    maxit = 1e3,
                    verbose = FALSE)
  t252ll[[s]] <- proc.time() - t0
  cvs252ll[[s]] <- cv252ll
  plot_sox_cv(cv252ll)
  plot_sox_sp(sox_obj = cv252ll,
              plot_min = TRUE,
              plot_1se = TRUE,
              lwd = 1)
}


# beepr::beep(8)

# non zeros
nz282 <- lapply(cvs282,
                FUN = function(x) {
                  which(x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] != 0)
                }
)
nz555 <- lapply(cvs555,
                FUN = function(x) {
                  which(x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] != 0)
                }
)
nz828 <- lapply(cvs828,
                FUN = function(x) {
                  which(x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] != 0)
                }
)
nz252 <- lapply(cvs252,
                FUN = function(x) {
                  which(x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] != 0)
                }
)
nz858 <- lapply(cvs858,
                FUN = function(x) {
                  which(x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] != 0)
                }
)
nz252l <- lapply(cvs252l,
                 FUN = function(x) {
                   which(x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] != 0)
                 }
)
nz252ll <- lapply(cvs252ll,
                  FUN = function(x) {
                    which(x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] != 0)
                  }
)

# joint detection rate
mean(sapply(nz282, FUN = function(x) {
  length(setdiff(nz.true, x)) == 0
}, simplify = TRUE))
mean(sapply(nz555, FUN = function(x) {
  length(setdiff(nz.true, x)) == 0
}, simplify = TRUE))
mean(sapply(nz828, FUN = function(x) {
  length(setdiff(nz.true, x)) == 0
}, simplify = TRUE))
mean(sapply(nz252, FUN = function(x) {
  length(setdiff(nz.true, x)) == 0
}, simplify = TRUE))
mean(sapply(nz858, FUN = function(x) {
  length(setdiff(nz.true, x)) == 0
}, simplify = TRUE))
mean(sapply(nz252l, FUN = function(x) {
  length(setdiff(nz.true, x)) == 0
}, simplify = TRUE))
mean(sapply(nz252ll, FUN = function(x) {
  length(setdiff(nz.true, x)) == 0
}, simplify = TRUE))

# missing rate
nz.length <- length(nz.true)
mean(sapply(nz282, FUN = function(x) {
  length(setdiff(nz.true, x))/nz.length
}))
mean(sapply(nz555, FUN = function(x) {
  length(setdiff(nz.true, x))/nz.length
}))
mean(sapply(nz828, FUN = function(x) {
  length(setdiff(nz.true, x))/nz.length
}))
mean(sapply(nz252, FUN = function(x) {
  length(setdiff(nz.true, x))/nz.length
}))
mean(sapply(nz858, FUN = function(x) {
  length(setdiff(nz.true, x))/nz.length
}))
mean(sapply(nz252l, FUN = function(x) {
  length(setdiff(nz.true, x))/nz.length
}))
mean(sapply(nz252ll, FUN = function(x) {
  length(setdiff(nz.true, x))/nz.length
}))

# false alarm rate
z.length <- p - length(nz.true)
mean(sapply(nz282, FUN = function(x) {
  length(setdiff(x, nz.true))/z.length
}))
mean(sapply(nz555, FUN = function(x) {
  length(setdiff(x, nz.true))/z.length
}))
mean(sapply(nz828, FUN = function(x) {
  length(setdiff(x, nz.true))/z.length
}))
mean(sapply(nz252, FUN = function(x) {
  length(setdiff(x, nz.true))/z.length
}))
mean(sapply(nz858, FUN = function(x) {
  length(setdiff(x, nz.true))/z.length
}))
mean(sapply(nz252l, FUN = function(x) {
  length(setdiff(x, nz.true))/z.length
}))
mean(sapply(nz252ll, FUN = function(x) {
  length(setdiff(x, nz.true))/z.length
}))

# # R1S
# r1s.sox <- r1s.sox.db <- r1s.glmnet <- r1s.glmnet.db <- 0
# for (i in 1:n.sim) {
#   nz.tmp <- nz.sox[[i]]
#   intx.tmp <- nz.tmp[nz.tmp > p]
#   main.tmp <- nz.tmp[nz.tmp <= p]
#   
#   for (j in intx.tmp) {
#     mains <- two.way.intx[j - p, ]
#     if (sum(mains %in% main.tmp) < 2) {
#       r1s.sox <- r1s.sox + 1
#     }
#   }
#   
#   nz.tmp.db <- nz.sox.db[[i]]
#   intx.tmp.db <- nz.tmp.db[nz.tmp.db > p]
#   main.tmp.db <- nz.tmp.db[nz.tmp.db <= p]
#   
#   for (j in intx.tmp.db) {
#     mains.db <- two.way.intx[j - p, ]
#     if (sum(mains.db %in% main.tmp.db) < 2) {
#       r1s.sox.db <- r1s.sox.db + 1
#     }
#   }
#   
#   nz.tmp <- nz.glmnet[[i]]
#   intx.tmp <- nz.tmp[nz.tmp > p]
#   main.tmp <- nz.tmp[nz.tmp <= p]
#   
#   for (j in intx.tmp) {
#     mains <- two.way.intx[j - p, ]
#     if (sum(mains %in% main.tmp) < 2) {
#       r1s.glmnet <- r1s.glmnet + 1
#     }
#   }
#   
#   nz.tmp.db <- nz.glmnet.db[[i]]
#   intx.tmp.db <- nz.tmp.db[nz.tmp.db > p]
#   main.tmp.db <- nz.tmp.db[nz.tmp.db <= p]
#   
#   for (j in intx.tmp.db) {
#     mains.db <- two.way.intx[j - p, ]
#     if (sum(mains.db %in% main.tmp.db) < 2) {
#       r1s.glmnet.db <- r1s.glmnet.db + 1
#     }
#   }
# }
# 1 - (r1s.sox/n.sim/(p * (p - 1) /2))
# 1 - (r1s.sox.db/n.sim/(p * (p - 1) /2))
# 1 - (r1s.glmnet/n.sim/(p * (p - 1) /2))
# 1 - (r1s.glmnet.db/n.sim/(p * (p - 1) /2))
# 
# RCI
rci282 <- rci555 <- rci828 <- rci252 <- rci858 <- rci252l <- rci252ll <- numeric(n.sim)
for (s in 1:n.sim) {
  data <- data.list[[s]]
  x <- as.matrix(data[, -(1:5)])
  time = data$Start
  time2 = data$Stop
  event = data$Event

  x282 <- as.matrix(x[, nz282[[s]]])
  x555 <- x[, nz555[[s]]]
  x828 <- x[, nz828[[s]]]
  x252 <- x[, nz252[[s]]]
  x858 <- x[, nz858[[s]]]
  x252l <- x[, nz252l[[s]]]
  x252ll <- x[, nz252ll[[s]]]
  
  cox282 <- survival::coxph(survival::Surv(time, time2, event) ~ x282)
  rci282[s] <- cox282$concordance["concordance"]

  cox555 <- survival::coxph(survival::Surv(time, time2, event) ~ x555)
  rci555[s] <- cox555$concordance["concordance"]

  cox828 <- survival::coxph(survival::Surv(time, time2, event) ~ x828)
  rci828[s] <- cox828$concordance["concordance"]

  cox252 <- survival::coxph(survival::Surv(time, time2, event) ~ x252)
  rci252[s] <- cox252$concordance["concordance"]
  
  cox858 <- survival::coxph(survival::Surv(time, time2, event) ~ x858)
  rci858[s] <- cox858$concordance["concordance"]
  
  cox252l <- survival::coxph(survival::Surv(time, time2, event) ~ x252l)
  rci252l[s] <- cox252l$concordance["concordance"]
  
  cox252ll <- survival::coxph(survival::Surv(time, time2, event) ~ x252ll)
  rci252ll[s] <- cox252ll$concordance["concordance"]
}

mean(rci282)
mean(rci555)
mean(rci828)
mean(rci252)
mean(rci858)
mean(rci252l)
mean(rci252ll)

# MSE
mean(sapply(cvs282,
            FUN = function(x) {
              mean((x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] - beta.true)^2)
            }))
mean(sapply(cvs555,
            FUN = function(x) {
              mean((x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] - beta.true)^2)
            }))
mean(sapply(cvs828,
            FUN = function(x) {
              mean((x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] - beta.true)^2)
            }))
mean(sapply(cvs252,
            FUN = function(x) {
              mean((x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] - beta.true)^2)
            }))
mean(sapply(cvs858,
            FUN = function(x) {
              mean((x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] - beta.true)^2)
            }))
mean(sapply(cvs252l,
            FUN = function(x) {
              mean((x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] - beta.true)^2)
            }))
mean(sapply(cvs252ll,
            FUN = function(x) {
              mean((x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.1se)] - beta.true)^2)
            }))

# CV error
cve282 <- cve555 <- cve828 <- cve252 <- cve858 <- cve252l <- cve252ll <- numeric(n.sim)
for (s in 1:n.sim) {
  cve282[s] <- cvs282[[s]]$cvm[which(cvs282[[s]]$lambdas == cvs282[[s]]$lambda.1se)]
  cve555[s] <- cvs555[[s]]$cvm[which(cvs555[[s]]$lambdas == cvs555[[s]]$lambda.1se)]
  cve828[s] <- cvs828[[s]]$cvm[which(cvs828[[s]]$lambdas == cvs828[[s]]$lambda.1se)]
  cve252[s] <- cvs252[[s]]$cvm[which(cvs252[[s]]$lambdas == cvs252[[s]]$lambda.1se)]
  cve858[s] <- cvs858[[s]]$cvm[which(cvs858[[s]]$lambdas == cvs858[[s]]$lambda.1se)]
  cve252l[s] <- cvs252l[[s]]$cvm[which(cvs252l[[s]]$lambdas == cvs252l[[s]]$lambda.1se)]
  cve252ll[s] <- cvs252ll[[s]]$cvm[which(cvs252ll[[s]]$lambdas == cvs252ll[[s]]$lambda.1se)]
}
mean(cve282)
mean(cve555)
mean(cve828)
mean(cve252)
mean(cve858)
mean(cve252l)
mean(cve252ll)

# 
# # Timing
# t.mean.sox <- 0
# t.mean.sox.db <- 0
# for (s in 1:n.sim) {
#   t.mean.sox <- t.mean.sox + t.sox[[s]][3]
#   t.mean.sox.db <- t.mean.sox.db + t.sox.db[[s]][3]
# }
# t.mean.sox/n.sim
# t.mean.sox.db/n.sim
# 
# 
# 
