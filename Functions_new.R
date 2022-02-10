
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

gen_X=function(m){
  A=gen_cat(m);B=gen_con(m);C=gen_cat(m)
  A1=ifelse(A==1,1,0);A2=ifelse(A==2,1,0)
  C1=ifelse(C==1,1,0);C2=ifelse(C==2,1,0)
  A1B=A1*B;A2B=A2*B
  C1B=C1*B;C2B=C2*B
  return(as.matrix(cbind(A1,A2,C1,C2,B,A1B,A2B,C1B,C2B)))
}

gen_X_n=function(m,n){
  Xn=NULL
  for (i in 1:n) {
    X=gen_X(m)
    Xn=rbind(Xn,X)
  }
  return(Xn)
}


# data <- permalgorithm(n, m, Xmat=gen_X_n(m,n),
#                       XmatNames=c("A1","A2","B","C1","C2","A1B","A2B","C1B","C2B"),
#                       #eventRandom = eventRandom, censorRandom=censorRandom,
#                       betas=c(log(3), log(3), rep(0,7)),groupByD=FALSE )
# coxph(Surv(Start, Stop, Event) ~ A1+A2+B+C1+C2+A1B+A2B+C1B+C2B,data)
# #lam1=2;regul='graph';betas=rep(0,5)  
# n=50;m=50
fit.cox=function(beta,data){coxph(Surv(Start, Stop, Event) ~ .,data,init = beta, control=list('iter.max'=0,timefix = FALSE))}
l=function(beta,data){-fit.cox(beta,data)$loglik[1]/n}
l.d=function(beta,data){-fit.cox(beta,data)$first/n}

# dataa <- permalgorithm(n, m, Xmat=gen_X_n(m,n)[,1],
#                       XmatNames=c("A1"),
#                       #eventRandom = eventRandom, censorRandom=censorRandom,
#                       betas=c(log(3)),groupByD=FALSE )
# coxph(Surv(Start, Stop, Event) ~ A1,dataa)
# x=as.matrix(seq(0.8,1.5,0.01),byrow=T)
# y=apply(x, 1, FUN = function(x) l(x,dataa[,-c(1,3)]))
# z=apply(x, 1, FUN = function(x) -l.d(x,dataa[,-c(1,3)]))
# plot(x,y);abline(v=1.1966)
# plot(x,z);abline(v=1.1966,h=0)
# https://stats.stackexchange.com/questions/187281/how-to-compute-partial-log-likelihood-function-in-cox-proportional-hazards-model 


SurvTDSelect=function(data, grp, grpV, eta_g, regul, betas, t, alpha, epsilon, lam1){
  graph = list('eta_g'= eta_g,'groups' = grp,'groups_var' = grpV)
  #lambda1=lam1/t
  n=length(unique(data$Id));data=data[,-c(1,3)]
  U=matrix(betas-l.d(betas,data)*t, byrow=FALSE)
  V=spams.proximalGraph(U,graph=graph,lambda1=lam1,regul="graph")
  while(sum(abs(betas-V))>epsilon){
    while( l(V,data) > l(betas,data)+l.d(betas,data) %*% (V-betas)+sum((V-betas)^2)/(2*t) ){
      t = alpha * t
      #U=matrix(betas-l.d(betas,data)*t, byrow=FALSE)
      #V=spams.proximalGraph(U,graph=graph,lambda1=lam1/t,regul="graph")
    }
    betas=V
    U=matrix(betas-l.d(betas,data)*t, byrow=FALSE)
    V=spams.proximalGraph(U,graph=graph,lambda1=lam1,regul="graph")
  }
  return(V)
}



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
    #error[i] = -2*n*(-l(est.train,data[,-c(1,3)])+l(est.train,train.data[,-c(1,3)]))
    #error[i] = 2*n*(l(est.train,test.data[,-c(1,3)]))
    # fit.original.test = coxph(Surv(Start, Stop, Event) ~ . ,test.data[,-c(1,3)])
    # error[i] = 2*(fit.original.test$loglik[2]+l(est.train,test.data[,-c(1,3)])*n)
    # # https://github.com/cran/glmnet/blob/master/R/coxnet.deviance.R line 170
    # https://arxiv.org/pdf/1905.10432.pdf 
    error[i] = coxnet.deviance(y=Surv(data$Start, data$Stop, data$Event),
                               x=as.matrix(data[,-c(1:5)]),beta=est.train)-
      coxnet.deviance(y=Surv(train.data$Start, train.data$Stop, train.data$Event),
                      x=as.matrix(train.data[,-c(1:5)]),beta=est.train)
    error[i] = error[i]/length(which(test.data$Event==1))
  }
  return (error)
}
