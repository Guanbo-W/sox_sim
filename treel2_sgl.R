library(sox)
library(SGL)
library(grpCox)
set.seed(20230218)

n <- 100
m <- 10
p <- 200

p.m <- p/m

g <- 1

# beta true
beta.true.m <- rep(0, p.m)
beta.true.m[1:5] <- (1:5)/10

beta.true <- rep(0, p)
beta.true[1:(p.m * g)] <- rep(beta.true.m, g)

nz.true <- which(beta.true != 0)

# sox
tree.list <- list(1:p)
tree.list <- c(tree.list, split(1:p, rep(1:m, each = p.m)))
tree.list <- c(tree.list, as.list(1:p))
names(tree.list) <- NULL

tree.pars <- tree_builder(tree.list)
tree.pars$eta_g[2 + ((1:m) - 1) * (p.m + 1)] <- tree.pars$eta_g[2 + ((1:m) - 1) * (p.m + 1)] * sqrt(20)

grp.grpcox <- tree.list[-1]

# No. of simulation replications
n.sim <- 20

cvs.sox <- cvs.sgl <- cvs.grp <- vector(mode = "list", length = n.sim)

for (s in 1:n.sim) {
  
  x <- matrix(rnorm(n * p), nrow = n)
  x <- scale(x)
  
  # simulate survival data
  T.max <- 100
  dat <- coxed::sim.survdata(T = T.max, X = x[, 1:(p.m * g)], censor = 0.3, beta = beta.true[1:(p.m * g)])
  
  data.sgl <- list()
  data.sgl$x <- x
  data.sgl$time <- as.numeric(dat$data$y)
  data.sgl$status <- as.integer(dat$data$failed)
  group.sgl <- rep(1:m, each = p.m)
  
  cv.sgl <- cvSGL(data.sgl,
                  index = group.sgl,
                  type = "cox",
                  standardize = FALSE,
                  alpha = 0.5)
  plot(cv.sgl)
  abline(v = log(cv.sgl$lambdas[which.min(cv.sgl$lldiff)]), lty = 2)
  
  cvs.sgl[[s]] <- cv.sgl
  
  cv.sox <- sox_cv(x = x,
                   ID = seq(n),
                   time = rep(0, n),
                   time2 = as.numeric(dat$data$y),
                   event = dat$data$failed,
                   penalty = "nested",
                   lambda = cv.sgl$lambdas,
                   group = tree.pars$groups,
                   own_variable = tree.pars$own_variables,
                   no_own_variable = tree.pars$N_own_variables,
                   penalty_weights = tree.pars$eta_g,
                   nfolds = 10,
                   foldid = cv.sgl$foldid,
                   tol = 1e-6,
                   maxit = 1e3,
                   verbose = FALSE)
  plot_sox_cv(cv.sox)
  
  cvs.sox[[s]] <- cv.sox
  
  y.grpcox <- data.frame(time = dat$data$y, status = dat$data$failed)
  cv.grp <- cv.grpCoxOverlap(X0 = x,
                             y = y.grpcox,
                             group = grp.grpcox,
                             penalty = "glasso",
                             # lambda = cv.sgl$lambdas,
                             standardize = TRUE)
  plot.llCV(cv.grp)
  
  cvs.grp[[s]] <- cv.grp
}


nz.sgl <- nz.sox <- nz.grp <- vector(mode = "list", length = n.sim)
for (s in 1:n.sim) {
  cv.sgl <- cvs.sgl[[s]]
  cv.sox <- cvs.sox[[s]]
  cv.grp <- cvs.grp[[s]]
  
  nz.sgl[[s]] <- which(cv.sgl$fit$beta[, which.min(cv.sgl$lldiff)] != 0)
  nz.sox[[s]] <- which(cv.sox$sox.fit$estimates[, which(cv.sox$lambdas == cv.sox$lambda.min)] != 0)
  nz.grp[[s]] <- which(cv.grp$aBetaOri[, which(cv.grp$lambda == cv.grp$lambda.max)] != 0)
}


# setdiff(nz.sgl, nz.sox)
# setdiff(nz.sox, nz.sgl)

# joint detection rate
mean(sapply(nz.sox, FUN = function(x) {
  length(setdiff(nz.true, x)) == 0
}, simplify = TRUE))
mean(sapply(nz.sgl, FUN = function(x) {
  length(setdiff(nz.true, x)) == 0
}, simplify = TRUE))
mean(sapply(nz.grp, FUN = function(x) {
  length(setdiff(nz.true, x)) == 0
}, simplify = TRUE))

# missing rate
nz.length <- length(nz.true)
mean(sapply(nz.sox, FUN = function(x) {
  length(setdiff(nz.true, x))/nz.length
}))
mean(sapply(nz.sgl, FUN = function(x) {
  length(setdiff(nz.true, x))/nz.length
}))
mean(sapply(nz.grp, FUN = function(x) {
  length(setdiff(nz.true, x))/nz.length
}))


# false alarm rate
z.length <- p - nz.length
mean(sapply(nz.sox, FUN = function(x) {
  length(setdiff(x, nz.true))/z.length
}))
mean(sapply(nz.sgl, FUN = function(x) {
  length(setdiff(x, nz.true))/z.length
}))
mean(sapply(nz.grp, FUN = function(x) {
  length(setdiff(x, nz.true))/z.length
}))

# MSE
mean(sapply(cvs.sox,
            FUN = function(x) {
              mean((x$sox.fit$estimates[, which(x$sox.fit$lambdas == x$lambda.min)] - beta.true)^2)
            }))
mean(sapply(cvs.sgl,
            FUN = function(x) {
              mean((x$fit$beta[, which.min(x$lldiff)] - beta.true)^2)
            }))
mean(sapply(cvs.grp,
            FUN = function(x) {
              mean((x$aBetaOri[, which(x$lambda == x$lambda.max)] - beta.true)^2)
            }))


# save.image("/Users/YLIAN/Desktop/Research/Collab/GuanboWang/SVS_TD_Cox/revision/n100_p200_m10_g1.Rdata")