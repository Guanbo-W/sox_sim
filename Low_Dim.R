# load libraries
library(Rcpp)
library(doParallel)
library(doMC)
library(survival)
library(glmnet)
library(spams)
library(PermAlgo)


source("Low_Dim_Functions.R")

# for high performance computing
args = commandArgs(trailingOnly=TRUE)
start=as.numeric(args[1])
end=as.numeric(args[2])

SSS=NULL
for (i in start:end){
  set.seed(i)
  n=500;m=50
  covariates=gen_X_n(m,n)
  # generate outcomes
  data <- permalgorithm(n, m, covariates, 
                        XmatNames=c("A1","A2","C1","C2","B","A1B","A2B","C1B","C2B"),
                        #change according to scenario 1/2
                        betas=c(rep(log(3),2), rep(0,2), log(4), rep(log(3),2),rep(0,2)),groupByD=FALSE )
  fit.original = coxph(Surv(Start, Stop, Event) ~ . ,data[,-c(1,3)])

  # glmnet
  registerDoMC(cores = 50)
  A=cv.glmnet(as.matrix(data[,-c(1:5)]), Surv(data$Start, data$Stop, data$Event),family = "cox", nfolds = 10, parallel = TRUE)
  B=coef(glmnet(as.matrix(data[,-c(1:5)]), Surv(data$Start, data$Stop, data$Event),family = "cox", lambda = A$lambda.1se))
  B.min=coef(glmnet(as.matrix(data[,-c(1:5)]), Surv(data$Start, data$Stop, data$Event),family = "cox", lambda = A$lambda.min))
  fit.glm.1se=glmnet(as.matrix(data[,-c(1:5)]), Surv(data$Start, data$Stop, data$Event),family = "cox", lambda = A$lambda.1se)
  fit.glm.min=glmnet(as.matrix(data[,-c(1:5)]), Surv(data$Start, data$Stop, data$Event),family = "cox", lambda = A$lambda.min)
  
  # specify hyper-parameters reqiured for our method
  # number of covariate
  p=9
  # initial values
  betas=rep(0,p)
  # step size; convergence threshold; step shrinking multiplier
  t = 1;epsilon = 10^-5;alpha = 0.8
  # numnber of groups
  g=5
  # weights of each group
  eta_g = as.vector(rep(1,5),mode='double')
  # grouping structure, see http://thoth.inrialpes.fr/people/mairal/spams/doc-R/html/doc_spams006.html#sec27 
  grp =as(matrix(as.vector(c(0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0,
                             1, 1, 0, 0, 0,
                             0, 0, 0, 0, 0,
                             0, 1, 0, 1, 0),mode='logical'),ncol = 5,byrow = T),'CsparseMatrix')
  grpV =as(matrix(as.vector(c(1, 0, 0, 0, 0, #A1
                              1, 0, 0, 0, 0, #A2
                              0, 0, 0, 1, 0, #C1
                              0, 0, 0, 1, 0, #C2
                              0, 1, 0, 0, 0, #B
                              0, 0, 1, 0, 0, #A1B
                              0, 0, 1, 0, 0, #A2B
                              0, 0, 0, 0, 1, #C1B
                              0, 0, 0, 0, 1  #C2B
  ),mode='logical'),ncol = 5,byrow = T),'CsparseMatrix')
  
  # specify the lambda sequence/range
  lam1.max = 0.5
  lambda.min = 0.0001
  S_max = log(lam1.max)
  S_min = log(lambda.min)
  seq_S = seq(S_min, S_max,length.out=50)
  seq_lambda = exp(seq_S)
  length_lambda =length(seq_lambda)

  
  # parallel computing
  # nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))
  # cl = makeCluster(nodeslist, type = "PSOCK") 
  # registerDoParallel(cl)
  cl = makePSOCKcluster(50)
  registerDoParallel(cl)
  clusterEvalQ(cl, {
    library(survival)
    library(spams)
    library(glmnet)
  })
  clusterExport(cl, ls())
  
  coefs = NULL
  errors = NULL
  coefs=foreach(lam1 = seq_lambda, .combine=cbind) %dopar% {
    SurvTDSelect(data, grp, grpV, eta_g, regul, betas, t, alpha, epsilon, lam1)
  }
  
  devout=paste(".../coefs","_",start,"_",end,".csv", sep="")
  row.names(coefs)=c("A1","A2","C1","C2","B","A1B","A2B","C1B","C2B")
  write.csv(coefs, devout, row.names = T)

  #plot
  # coefs.col = matrix(t(coefs),ncol=1)
  # D1 = data.frame(coefs.col, var=as.character(rep(1:p, each= length_lambda)),lambda=rep(seq_lambda,p))
  # ggplot(D1, aes(x = lambda, y = coefs.col, color = var, group = var)) +
  #   geom_line(position=position_dodge(width=0.01)) + labs(color = 'var')+ xlab(expression(lambda))+ylab("Coefficients")+
  #   geom_dl(aes(label = var), method = list("first.points"), cex = 0.02)
  #dev.off()
  
  # CV error plot
  K = 10
  # error_on_each_fold=NULL
  # for (lam1 in seq_lambda) {
  #   temp=SurvTDSelect_CV(data,grp, grpV, eta_g, regul, betas, t, alpha, epsilon, lam1, K)
  #   error_on_each_fold=cbind(error_on_each_fold,temp)
  # }
  
  # cross-validation
  error_on_each_fold = foreach(lam1 = seq_lambda, .combine=cbind) %dopar% {
    SurvTDSelect_CV(data,grp, grpV, eta_g, regul, betas, t, alpha, epsilon, lam1, K)
  }
  stopCluster(cl)
  
  # 1se rule
  errors = error_uppers = error_lowers = NULL
  errors = apply(error_on_each_fold, 2, mean)
  error_uppers = errors + apply(error_on_each_fold, 2, function(x) sd(x)/sqrt(K)  )
  error_lowers = errors - apply(error_on_each_fold, 2, function(x) sd(x)/sqrt(K) )

  # D2 = data.frame(errors, error_uppers, error_lowers, lambda=seq_lambda)
  
  lambda.min = seq_lambda[errors == min(errors)]
  min_error.plus_one_se = error_uppers[errors == min(errors)]
  lambda.1se = max(seq_lambda[errors <= min_error.plus_one_se])
  
  
  # plot
  # D2 = data.frame(log(errors), error_uppers, error_lowers, lambda=log(seq_lambda))
  # ggplot(D2, aes(x=log(seq_lambda), y=errors, colour="red")) + theme(legend.position = "none")+
  #   geom_point()+ xlab(expression(log(lambda))) + ylab("Mean (1SE) Neg Ave Log Lik")+
  #   geom_errorbar(aes(ymin=error_lowers, ymax=error_uppers), width=.01, colour="grey", position =position_dodge(0.1))+
  #   geom_vline(xintercept = log(lambda.min), linetype="dotted",
  #              color = "blue", size=0.5) + geom_text(aes(log(lambda.min),error_lowers[log(seq_lambda)==log(lambda.min)],label="lambda.min"),size=3, colour="blue") +
  #   geom_vline(xintercept = log(lambda.1se), linetype="dotted",
  #              color = "darkgreen", size=0.5) + geom_text(aes(log(lambda.1se),error_lowers[log(seq_lambda)==log(lambda.min)],label="lambda.1se"),size=3, colour="darkgreen")

  
  COEF = coefs[,seq_lambda==lambda.1se]
  COEF.min= coefs[,seq_lambda==lambda.min]
  
  # testing measurement
  # 1) mean cross-validated error (log partial likelihood deviance)
  LPLD.1se=errors[seq_lambda == lambda.1se]
  LPLD.min=errors[seq_lambda == lambda.min]
  LPLD.1se.glm=A$cvm[A$lambda.1se==A$lambda]
  LPLD.min.glm=A$cvm[A$lambda.min==A$lambda]
  
  # 2) MSE
  truebetas=c(rep(log(3),2), rep(0,2), log(4), rep(log(3),2),rep(0,2))
  MSE.full = mean((truebetas-fit.original$coefficients)^2)
  MSE.1se = mean((truebetas-COEF)^2)
  MSE.min=mean((truebetas-COEF.min)^2)
  MSE.1se.glm=mean((truebetas-B)^2)
  MSE.min.glm=mean((truebetas-B.min)^2)
  #MSE.all = rbind(MSE.all, MSE)
  
  # 3) joint detection rate
  JDR.1se = as.numeric(sum(COEF[c(which (truebetas != 0))]!=0)==length(which (truebetas != 0)))
  JDR.min = as.numeric(sum(COEF.min[c(which (truebetas != 0))]!=0)==length(which (truebetas != 0)))
  JDR.1se.glm = as.numeric(sum(B[c(which (truebetas != 0))]!=0)==length(which (truebetas != 0)))
  JDR.min.glm = as.numeric(sum(B.min[c(which (truebetas != 0))]!=0)==length(which (truebetas != 0)))
  
  # 4) missing rate
  MR.1se = length(which(truebetas != 0 & COEF == 0))/length(which (truebetas != 0) )
  MR.min = length(which(truebetas != 0 & COEF.min == 0))/length(which (truebetas != 0) )
  MR.1se.glm = length(which(truebetas != 0 & B == 0))/length(which (truebetas != 0) )
  MR.min.glm = length(which(truebetas != 0 & B.min == 0))/length(which (truebetas != 0) )
  
  # 5) false alarm rate
  FAR.1se = length(which(truebetas == 0 & COEF != 0))/length(which (truebetas == 0) )
  FAR.min = length(which(truebetas == 0 & COEF.min != 0))/length(which (truebetas == 0) )
  FAR.1se.glm = length(which(truebetas == 0 & B != 0))/length(which (truebetas == 0) )
  FAR.min.glm = length(which(truebetas == 0 & B.min != 0))/length(which (truebetas == 0) )
  
  # 6) consistency of categorical variable selection
  # 1 if the dummy variables for each categorical variable are selected collectively
  all=COEF!=0
  A=all[1]==all[2];C=all[3]==all[4];AB=all[6]==all[7];BC=all[8]==all[9]
  CCVS.1se=(A+C+AB+BC)/4
  all=COEF.min!=0
  A=all[1]==all[2];C=all[3]==all[4];AB=all[6]==all[7];BC=all[8]==all[9]
  CCVS.min=(A+C+AB+BC)/4
  all=B!=0
  A=all[1]==all[2];C=all[3]==all[4];AB=all[6]==all[7];BC=all[8]==all[9]
  CCVS.1se.glm=(A+C+AB+BC)/4
  all=B.min!=0
  A=all[1]==all[2];C=all[3]==all[4];AB=all[6]==all[7];BC=all[8]==all[9]
  CCVS.min.glm=(A+C+AB+BC)/4
  
  # 7) consistency of strong heredity/conditional selection
  # 1 if the two strong heredity selection rules are respected
  all=COEF!=0
  CHS.AB1=ifelse( F %in% c(all[c(1,2,5,6,7)]==rep(T,5) ),0,1);CHS.AB2=ifelse( F %in% c(all[c(1,2,5,6,7)]==c(F,F,T,F,F) ),0,1);CHS.AB3=ifelse( F %in% c(all[c(1,2,5,6,7)]==c(T,T,F,F,F) ),0,1);CHS.AB4=ifelse( F %in% c(all[c(1,2,5,6,7)]==c(T,T,T,F,F) ),0,1);CHS.AB5=ifelse( F %in% c(all[c(1,2,5,6,7)]==rep(F,5) ),0,1)
  CHS.BC1=ifelse( F %in% c(all[c(3,4,5,8,9)]==rep(T,5) ),0,1);CHS.BC2=ifelse( F %in% c(all[c(3,4,5,8,9)]==c(F,F,T,F,F) ),0,1);CHS.BC3=ifelse( F %in% c(all[c(3,4,5,8,9)]==c(T,T,F,F,F) ),0,1);CHS.BC4=ifelse( F %in% c(all[c(3,4,5,8,9)]==c(T,T,T,F,F) ),0,1);CHS.BC5=ifelse( F %in% c(all[c(3,4,5,8,9)]==rep(F,5) ),0,1)
  CHS.1se=(CHS.AB1+CHS.AB2+CHS.AB3+CHS.AB4+CHS.AB5+CHS.BC1+CHS.BC2+CHS.BC3+CHS.BC4+CHS.BC5)/2
  all=COEF.min!=0
  CHS.AB1=ifelse( F %in% c(all[c(1,2,5,6,7)]==rep(T,5) ),0,1);CHS.AB2=ifelse( F %in% c(all[c(1,2,5,6,7)]==c(F,F,T,F,F) ),0,1);CHS.AB3=ifelse( F %in% c(all[c(1,2,5,6,7)]==c(T,T,F,F,F) ),0,1);CHS.AB4=ifelse( F %in% c(all[c(1,2,5,6,7)]==c(T,T,T,F,F) ),0,1);CHS.AB5=ifelse( F %in% c(all[c(1,2,5,6,7)]==rep(F,5) ),0,1)
  CHS.BC1=ifelse( F %in% c(all[c(3,4,5,8,9)]==rep(T,5) ),0,1);CHS.BC2=ifelse( F %in% c(all[c(3,4,5,8,9)]==c(F,F,T,F,F) ),0,1);CHS.BC3=ifelse( F %in% c(all[c(3,4,5,8,9)]==c(T,T,F,F,F) ),0,1);CHS.BC4=ifelse( F %in% c(all[c(3,4,5,8,9)]==c(T,T,T,F,F) ),0,1);CHS.BC5=ifelse( F %in% c(all[c(3,4,5,8,9)]==rep(F,5) ),0,1)
  CHS.min=(CHS.AB1+CHS.AB2+CHS.AB3+CHS.AB4+CHS.AB5+CHS.BC1+CHS.BC2+CHS.BC3+CHS.BC4+CHS.BC5)/2
  all=B!=0
  CHS.AB1=ifelse( F %in% c(all[c(1,2,5,6,7)]==rep(T,5) ),0,1);CHS.AB2=ifelse( F %in% c(all[c(1,2,5,6,7)]==c(F,F,T,F,F) ),0,1);CHS.AB3=ifelse( F %in% c(all[c(1,2,5,6,7)]==c(T,T,F,F,F) ),0,1);CHS.AB4=ifelse( F %in% c(all[c(1,2,5,6,7)]==c(T,T,T,F,F) ),0,1);CHS.AB5=ifelse( F %in% c(all[c(1,2,5,6,7)]==rep(F,5) ),0,1)
  CHS.BC1=ifelse( F %in% c(all[c(3,4,5,8,9)]==rep(T,5) ),0,1);CHS.BC2=ifelse( F %in% c(all[c(3,4,5,8,9)]==c(F,F,T,F,F) ),0,1);CHS.BC3=ifelse( F %in% c(all[c(3,4,5,8,9)]==c(T,T,F,F,F) ),0,1);CHS.BC4=ifelse( F %in% c(all[c(3,4,5,8,9)]==c(T,T,T,F,F) ),0,1);CHS.BC5=ifelse( F %in% c(all[c(3,4,5,8,9)]==rep(F,5) ),0,1)
  CHS.1se.glm=(CHS.AB1+CHS.AB2+CHS.AB3+CHS.AB4+CHS.AB5+CHS.BC1+CHS.BC2+CHS.BC3+CHS.BC4+CHS.BC5)/2
  all=B.min!=0
  CHS.AB1=ifelse( F %in% c(all[c(1,2,5,6,7)]==rep(T,5) ),0,1);CHS.AB2=ifelse( F %in% c(all[c(1,2,5,6,7)]==c(F,F,T,F,F) ),0,1);CHS.AB3=ifelse( F %in% c(all[c(1,2,5,6,7)]==c(T,T,F,F,F) ),0,1);CHS.AB4=ifelse( F %in% c(all[c(1,2,5,6,7)]==c(T,T,T,F,F) ),0,1);CHS.AB5=ifelse( F %in% c(all[c(1,2,5,6,7)]==rep(F,5) ),0,1)
  CHS.BC1=ifelse( F %in% c(all[c(3,4,5,8,9)]==rep(T,5) ),0,1);CHS.BC2=ifelse( F %in% c(all[c(3,4,5,8,9)]==c(F,F,T,F,F) ),0,1);CHS.BC3=ifelse( F %in% c(all[c(3,4,5,8,9)]==c(T,T,F,F,F) ),0,1);CHS.BC4=ifelse( F %in% c(all[c(3,4,5,8,9)]==c(T,T,T,F,F) ),0,1);CHS.BC5=ifelse( F %in% c(all[c(3,4,5,8,9)]==rep(F,5) ),0,1)
  CHS.min.glm=(CHS.AB1+CHS.AB2+CHS.AB3+CHS.AB4+CHS.AB5+CHS.BC1+CHS.BC2+CHS.BC3+CHS.BC4+CHS.BC5)/2
  
  # 8) concordance index
  CI.full = concordance(fit.original)$concordance
  VARIABLES = c("A1","A2","C1","C2","B","A1B","A2B","C1B","C2B")
  if ( sum(COEF) == 0 ){
    fit.select = coxph(Surv(Start, Stop, Event) ~ 1,data[,-c(1,3)])
    CI.1se = concordance(fit.select)$concordance
  } else {
    variables = VARIABLES[which (COEF != 0)]
    f = as.formula(paste("Surv(Start, Stop, Event)",paste(variables,collapse="+"),sep="~"))
    fit.select = coxph(f, data)
    CI.1se = concordance(fit.select)$concordance
  }
  if ( sum(COEF.min) == 0 ){
    fit.select = coxph(Surv(Start, Stop, Event) ~ 1,data[,-c(1,3)])
    CI.min = concordance(fit.select)$concordance
  } else {
    variables = VARIABLES[which (COEF.min != 0)]
    f = as.formula(paste("Surv(Start, Stop, Event)",paste(variables,collapse="+"),sep="~"))
    fit.select = coxph(f, data)
    CI.min = concordance(fit.select)$concordance
  }
  if ( sum(B) == 0 ){
    fit.select = coxph(Surv(Start, Stop, Event) ~ 1,data[,-c(1,3)])
    CI.1se.glm = concordance(fit.select)$concordance
  } else {
    variables = VARIABLES[which (B != 0)]
    f = as.formula(paste("Surv(Start, Stop, Event)",paste(variables,collapse="+"),sep="~"))
    fit.select = coxph(f, data)
    CI.1se.glm = concordance(fit.select)$concordance
  }
  if ( sum(B.min) == 0 ){
    fit.select = coxph(Surv(Start, Stop, Event) ~ 1,data[,-c(1,3)])
    CI.min.glm = concordance(fit.select)$concordance
  } else {
    variables = VARIABLES[which (B.min != 0)]
    f = as.formula(paste("Surv(Start, Stop, Event)",paste(variables,collapse="+"),sep="~"))
    fit.select = coxph(f, data)
    CI.min.glm = concordance(fit.select)$concordance
  }
  
  
  us.1se=c(COEF,JDR.1se, MR.1se, FAR.1se, CCVS.1se, CHS.1se, CI.1se, LPLD.1se, MSE.1se)
  us.min=c(COEF.min,JDR.min, MR.min, FAR.min, CCVS.min, CHS.min, CI.min, LPLD.min, MSE.min)
  glm.1se=c(as.vector(B),JDR.1se.glm, MR.1se.glm, FAR.1se.glm, CCVS.1se.glm, CHS.1se.glm, CI.1se.glm, LPLD.1se.glm, MSE.1se.glm)
  glm.min=c(as.vector(B.min),JDR.min.glm, MR.min.glm, FAR.min.glm, CCVS.min.glm, CHS.min.glm, CI.min.glm, LPLD.min.glm, MSE.min.glm)
  results = cbind(us.1se,us.min,glm.1se,glm.min)
  row.names(results)=c("A1","A2","C1","C2","B","A1B","A2B","C1B","C2B","JDR", "MR", "FAR", "CCVS", "CHS", "CI", "LPLD.1se", "MSE")
  
  SSS=cbind(SSS,results)
  
  devout=paste(".../result","_",start,"_",end,".csv", sep="")
  write.csv(SSS, devout, row.names = T)
  
}



