# generate B
gen_con=function(m){
  X=rnorm(m/5)
  XX=NULL
  for (i in 1:length(X)) {
    if (length(XX)<m){
    X.rep=rep(X[i],round(runif(1,5,10),0))
    XX=c(XX,X.rep)
    }
  }
  return(XX[1:m])
}

# generate A and C
gen_cat=function(m){
  X=sample.int(3, m/5,replace = TRUE)
  XX=NULL
  for (i in 1:length(X)) {
    if (length(XX)<m){
      X.rep=rep(X[i],round(runif(1,5,10),0))
      XX=c(XX,X.rep)
    }
  }
  return(XX[1:m])
}

# generate covariate for one subject
gen_X=function(m){
  A=gen_cat(m);B=gen_con(m);C=gen_cat(m)
  A1=ifelse(A==1,1,0);A2=ifelse(A==2,1,0)
  C1=ifelse(C==1,1,0);C2=ifelse(C==2,1,0)
  A1B=A1*B;A2B=A2*B
  C1B=C1*B;C2B=C2*B
  return(as.matrix(cbind(A1,A2,C1,C2,B,A1B,A2B,C1B,C2B)))
}

# generate covariate for all subject
gen_X_n=function(m,n){
  Xn=NULL
  for (i in 1:n) {
    X=gen_X(m)
    Xn=rbind(Xn,X)
  }
  return(Xn)
}

# loss function and its first derivative
fit.cox=function(beta,data){coxph(Surv(Start, Stop, Event) ~ .,data,init = beta, control=list('iter.max'=0,timefix = FALSE))}
l=function(beta,data){-fit.cox(beta,data)$loglik[1]/n}
l.d=function(beta,data){-fit.cox(beta,data)$first/n}

# Proximal gradient descent with backtracking line search
SurvTDSelect=function(data, grp, grpV, eta_g, regul, betas, t, alpha, epsilon, lam1){
  graph = list('eta_g'= eta_g,'groups' = grp,'groups_var' = grpV)
  n=length(unique(data$Id));data=data[,-c(1,3)]
  U=matrix(betas-l.d(betas,data)*t, byrow=FALSE)
  V=spams.proximalGraph(U,graph=graph,lambda1=lam1,regul="graph")
  while(sum(abs(betas-V))>epsilon){
    while( l(V,data) > l(betas,data)+l.d(betas,data) %*% (V-betas)+sum((V-betas)^2)/(2*t) ){
      t = alpha * t
    }
    betas=V
    U=matrix(betas-l.d(betas,data)*t, byrow=FALSE)
    V=spams.proximalGraph(U,graph=graph,lambda1=lam1,regul="graph")
  }
  return(V)
}

# cross-validation
SurvTDSelect_CV = function(data,grp, grpV, eta_g, regul, betas, t, alpha, epsilon, lam1, K){
  error=NULL
  graph = list('eta_g'= eta_g,'groups' = grp,'groups_var' = grpV)
  n = length(unique(data$Id))
  # sample splitting
  cv.id = NULL
  for (j in 1: K){cv.id[data$Id %in% seq(n/K*(j-1)+1,n/K*j,1)] = j}
  for (i in 1:K){
    test.data = data[cv.id == i,]
    train.data = data[cv.id != i,]
    est.train = SurvTDSelect(train.data, grp, grpV, eta_g, regul, betas, t, alpha, epsilon, lam1)
    error[i] = coxnet.deviance(y=Surv(data$Start, data$Stop, data$Event),
                               x=as.matrix(data[,-c(1:5)]),beta=est.train)-
      coxnet.deviance(y=Surv(train.data$Start, train.data$Stop, train.data$Event),
                      x=as.matrix(train.data[,-c(1:5)]),beta=est.train)
    error[i] = error[i]/length(which(test.data$Event==1))
  }
  return (error)
}
