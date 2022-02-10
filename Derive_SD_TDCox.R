
library(ggm)
# unit function when the unit rule is select all variables in F
unit=function(F, V){
  # F is a vector, V is a vector
  b=powerset(setdiff(V, F), nonempty = FALSE)
  Result=lapply(b, FUN=function(x) union(F,x))
  return(Result)}

# select all variables in F or none (categorical variable selection)
cat=function(F, V){
  # F is a vector, V is a vector
  b=powerset(setdiff(V, F), nonempty = FALSE)
  Result=append(b, lapply(b, FUN=function(x) union(F,x)))
  return(Result)}

sort.list=function(A){
  return(lapply(A, FUN=function(x) sort(x)))
}

#change character to list
CCL=function(A){
  return(lapply(strsplit(A," "), FUN=function(x) as.numeric(x)))
}

#intersect.list
IL=function(A,B){
  A=sort.list(A)
  AA=sapply(A , paste, collapse = " ")
  B=sort.list(B)
  BB=sapply(B , paste, collapse = " ")
  return(as.list(AA[which(AA %in% BB)]))
}

#power set minus a set
set.difference=function(P,A){
  A=sort.list(A)
  P=sort.list(P)
  AA=sapply(A , paste, collapse = " ")
  PP=sapply(P , paste, collapse = " ")
  return(as.list(PP[-which(PP %in% AA)]))
}

#if then operation
ifthen=function(A,B,V){
  P=powerset(V, nonempty = FALSE)
  C=set.difference(P,A)
  D=IL(A,B)
  Result=append(C,D)
  Result=unique(sapply(Result , paste, collapse = " "))
  return(Result)
}



#################################################
V=c("A1","A2","A1B","A2B","B","C1","C2","C1B","C2B")
V=seq(9)
start_time <- Sys.time()
a1=as.list(cat(c(1,2),V))
a2=as.list(cat(c(3,4),V))
a3=as.list(cat(c(6,7),V))
a4=as.list(cat(c(8,9),V))
a5=as.list(ifthen(unit(c(3,4),V),unit(c(1,2,5),V),V))
a6=as.list(ifthen(unit(c(8,9),V),unit(c(5,6,7),V),V))
M=sapply(IL(IL(IL(IL(IL(a1,a2),a3),a4),a5),a6), paste, collapse = " ")
end_time <- Sys.time()
end_time - start_time
length(M)
