
library(quadprog)
library(mixtools)
library(dplyr)

#Parameters of optimisation methods
optim_trace <- 0
optim_maxiter <- 20000
EM_maxiter <- optim_maxiter
opt_reltol <- 10^(-13)
opt_factr=100
eps_opt=10^(-13)


## Implemented methods and helper functions
prod2 <- function(a1,b1,a2,b2){
  #scalar product of densities of two normal distributions
  #a1,a2 - means; b1,b2 - standard deviations
  sigma2=b1^2+b2^2
  return(1/(sqrt(2*pi*sigma2))*exp(-(a1-a2)^2/2/sigma2))}



find_weights <- function(mu,sigma,data,h,scale=TRUE){
  # Finds optimal weights for given means and standard deviations of a mixture
  # True density is approximated by kernel smoothing of data 
  # mu, sigma - means and standard deviations of mixture components
  # data - sample; h - kernel smoothing bandwidth
  # scale - if scaling of equations is done for numerical stability
  M=length(mu)
  N=length(data)
  F <- rep(0,M)
  A <- matrix(0,nrow=M,ncol=M)
  B <- matrix(0,nrow=M,ncol=N)
  for(i in 1:M){
    B[i,]=prod2(mu[i],sigma[i],data,h)
    F[i]=mean(B[i,])
    A[i,]=prod2(mu[i],sigma[i],mu,sigma)
  }
  Condition=rbind(rep(1,M),diag(1,M,M))
  #first condition: the sum of weigths is 1
  #the rest are non-negativity constraings
  bCond=c(1,rep(0,M))
  Dmat=(A+t(A))/2 #to force symmetricity
  if(scale){
    #scaling equations for numerical stability
    Dnorm=norm(Dmat,"2")
    answer=solve.QP(Dmat/Dnorm,F/Dnorm,t(Condition),bCond,meq=1) 
  } else {
    answer=solve.QP(Dmat,F,t(Condition),bCond,meq=1) 
  }
  return(list(answer$solution,A,B))
}


L2_distance=function(data,h,initial_mu,initial_sig,initial_v=NULL,dist_power=2){
  # 
  # Weights are found by solving a quadratic optimization problem
  # Kernel density estimate with bandwidth h is used to approximate the true density
  # initial_mu contains starting values for means
  # initial_sig contains starting values for standard deviations
  N=length(data)
  M=length(initial_mu)
  Pnorm2=0
  for(i in 1:N){
      Pnorm2=Pnorm2+sum(prod2(data[i],h,data,h))/N**2
  }
  distance<- function(mu_sig,m1=data,h1=h,norm=Pnorm2,exponent=dist_power/2){
    #finds the distance between the kernel density estimate
    #and the best approximation with the mixture of normal densities with given means and standard deviations
    #optimal parameters are found by using the default method of the optim command
    M=length(mu_sig)/2
    N=length(m1)
    mu=mu_sig[1:M]
    sigma=mu_sig[(M+1):(2*M)]
    tmp <- find_weights(mu,sigma,m1,h1)#optimal weights
    w<-tmp[[1]]
    A=tmp[[2]]
    B=tmp[[3]]
    return(abs(t(w)%*%A%*%w-2*t(w)%*%B%*%rep(1/N,N)+norm)^exponent)
  }
  result <- try(optim(c(initial_mu,sqrt(initial_sig^2+h^2)),distance,control=list(trace=optim_trace,maxit=optim_maxiter,reltol=opt_reltol)),TRUE) # reltol 16.05.2023
  if(class(result)!="try-error"){
    mu_sig=result$par
    mu=mu_sig[1:M]
    sigma=mu_sig[(M+1):(2*M)]
    tmp <- find_weights(mu,sigma,data,h)
    v<-tmp[[1]]
    iterations=result$counts[1]
    error=FALSE
    value_of_objective_fn=result$value
  } else {
    mu=rep(0,M)
    sigma=rep(0,M)
    v<-rep(0,M)
    iterations=0
    error=TRUE
    value_of_objective_fn=-1
  }
  return(list(mu=mu,v=v,sigma=sigma,iterations=iterations,error=error,value_of_objective_fn=value_of_objective_fn))
}

L2_distanceBFGS=function(data,h,initial_mu,initial_sig,initial_v=NULL,dist_power=2){
  #finds the distance between the kernel density estimate
  #and the best approximation with the mixture of normal densities with given means and standard deviations
  #optimal parameters are found by using the BFGS method for constrained optimisation (standard deviations have to be non-negative)
  N=length(data)
  M=length(initial_mu)
  Pnorm2=0
  for(i in 1:N){
      Pnorm2=Pnorm2+sum(prod2(data[i],h,data,h))/N**2
  }
  distance<- function(mu_sig,m1=data,h1=h,norm=Pnorm2,exponent=dist_power/2){
    M=length(mu_sig)/2
    N=length(m1)
    mu=mu_sig[1:M]
    sigma=mu_sig[(M+1):(2*M)]
    tmp <- find_weights(mu,sigma,m1,h1)
    w<-tmp[[1]]
    A=tmp[[2]]
    B=tmp[[3]]
    return(abs(t(w)%*%A%*%w-2*t(w)%*%B%*%rep(1/N,N)+norm)^exponent)
  }
  result <- try(optim(c(initial_mu,sqrt(initial_sig^2+h^2)),distance,method="L-BFGS-B",control=list(trace=optim_trace,maxit=optim_maxiter,factr=opt_factr),lower=c(rep(-Inf,length(initial_mu)),rep(0,length(initial_sig)))),TRUE) # algl?hendid muudetud; reltol kriteerium 17.05; kitsendused 19.V.23
  if(class(result)!="try-error"){
    mu_sig=result$par
    mu=mu_sig[1:M]
    sigma=mu_sig[(M+1):(2*M)]
    tmp <- find_weights(mu,sigma,data,h)
    v<-tmp[[1]]
    iterations=result$counts[1]
    error=FALSE
    value_of_objective_fn=result$value
  } else {
    mu=rep(0,M)
    sigma=rep(0,M)
    v<-rep(0,M)
    iterations=0
    error=TRUE
    value_of_objective_fn=-1
  }
  return(list(mu=mu,v=v,sigma=sigma,iterations=iterations,error=error,value_of_objective_fn=value_of_objective_fn))
}


L2_distance_corrected=function(data,h,initial_mu,initial_sig,initial_v=NULL,dist_power=2){
  #finds the distance between the kernel density estimate
  #and the best approximation with the mixture of normal densities with given means and bias corrected standard deviations
  #optimal parameters are found by using the BFGS method for constrained optimisation (standard deviations have to be non-negative)
   N=length(data)
   M=length(initial_mu)
   Pnorm2=0
   for(i in 1:N){
       Pnorm2=Pnorm2+sum(prod2(data[i],h,data,h))/N**2
   }
   distance<-function(mu_sig,m1=data,h1=h,norm=Pnorm2,exponent=dist_power/2){
     M=length(mu_sig)/2
     N=length(m1)
     mu=mu_sig[1:M]
     sigma=mu_sig[(M+1):(2*M)]
     tmp <- find_weights(mu,sigma,m1,h1)
     w<-tmp[[1]]
     A=tmp[[2]]
     B=tmp[[3]]
return(abs(t(w)%*%A%*%w-2*t(w)%*%B%*%rep(1/N,N)+norm)^exponent)
   }
   result <-try(optim(c(initial_mu,sqrt(initial_sig^2+h^2)),distance,method="L-BFGS-B",control=list(trace=optim_trace,maxit=optim_maxiter,factr=opt_factr),lower=c(rep(-Inf,length(initial_mu)),rep(h,length(initial_sig)))),TRUE) # factr 16.05.2022
   if(class(result)!="try-error"){
     mu_sig=result$par
     mu=mu_sig[1:M]
     sigma1=mu_sig[(M+1):(2*M)]
     sigma=sqrt(sigma1^2-h^2)
     tmp <- find_weights(mu,sigma1,data,h)
     v<-tmp[[1]]
     iterations=result$counts[1]
     error=FALSE
     value_of_objective_fn=result$value
   } else {
     mu=rep(0,M)
     sigma=rep(0,M)
     v<-rep(0,M)
     iterations=0
     error=TRUE
     value_of_objective_fn=-1
   }
return(list(mu=mu,v=v,sigma=sigma,iterations=iterations,error=error,value_of_objective_fn=value_of_objective_fn))
}



EM=function(data,h,initial_mu,initial_sig,initial_v,dist_power=2){
  #EM algorithm
  #maximal number of iterations should be given in the clobal variable EM_maxiter
  #Uses only given initial values (maxrestarts=0)
  #dist_power and h are not used, but included for compatibility with other methods
  # initial_v specifies initial weights for  v1-v_{m-1}
 result=try(normalmixEM(data,c(initial_v,1-sum(initial_v)),initial_mu,initial_sig,maxit=EM_maxiter,epsilon=eps_opt,maxrestarts = 0),TRUE)
  if(class(result)!="try-error"){
    mu=result$mu
    sigma=result$sigma
    v<-result$lambda
    iterations=length(result$all.loglik)-1
    error=FALSE
    value_of_objective_fn=result$loglik
  } else {
    M=length(initial_mu)
    mu=rep(0,M)
    sigma=rep(0,M)
    v<-rep(0,M)
    iterations=0
    error=TRUE
    value_of_objective_fn=-1
  }
  return(list(mu=mu,v=v,sigma=sigma,iterations=iterations,error=error,value_of_objective_fn=value_of_objective_fn))
}


pseudo_likelihood=function(data,h,initial_mu,initial_sig,initial_v,dist_power=2){
  #Maximizes pseudo likelihood function
  #initial_mu - initial values for means
  #initial_v initial values for first m-1 components
  #h - kernel bandwidth parameter
  m=data
  N=length(m)
  M=length(initial_mu)
  likelihood<- function(mu_sig,m1=m,h1=h){
    N=length(m1)
    M=length(mu_sig)%/%2
    mu=mu_sig[1:M]
    sigma=mu_sig[(M+1):(2*M)]
    w=find_weights(mu,sigma,m1,h1)[[1]]
    density=w[1]*dnorm(m1,mu[1],sigma[1])
    for(i in 2:M){
      density=density+w[i]*dnorm(m1,mu[i],sigma[i])
    }
    
    return(-sum(log(density)))
  }
  result <- try(optim(c(initial_mu,initial_sig),likelihood,control=list(trace=optim_trace,maxit=optim_maxiter,reltol=opt_reltol)))
   if(class(result)!="try-error"){
      params=result$par
      mu=params[1:M]
      sigma=params[(M+1):(2*M)]
      v=find_weights(mu,sigma,m,h)[[1]]
      iterations=result$counts[1]
      error=FALSE
      value_of_objective_fn=result$value
   } else {
      mu=rep(0,M)
      sigma=rep(0,M)
      v<-rep(0,M)
      iterations=0
       error=TRUE
       value_of_objective_fn=-1
   }
  return(list(mu=mu,v=v,sigma=sigma,iterations=iterations,error=error,value_of_objective_fn=-value_of_objective_fn))
}



pseudo_likelihood_corrected=function(data,h,initial_mu,initial_sig,initial_v,dist_power=2){
  #uses bias corrected variances when finding weights
  m=data
  N=length(m)
  likelihood<- function(mu_sig,m1=m,h1=h){
    N=length(m1)
    M=length(mu_sig)%/%2
    mu=mu_sig[1:M]
    sigma=mu_sig[(M+1):(2*M)]
    w=find_weights(mu,sqrt(sigma^2+h1^2),m1,h1)[[1]]
    density=w[1]*dnorm(m1,mu[1],sigma[1])
    for(i in 2:M){
      density=density+w[i]*dnorm(m1,mu[i],sigma[i])
    }
  return(-sum(log(density)))
  }
  result <- try(optim(c(initial_mu,initial_sig),likelihood,control=list(trace=optim_trace,maxit=optim_maxiter,reltol=opt_reltol)))
  # if(class(result)=="try-error"){
  #   result <- try(optim(c(initial_mu,initial_v,initial_sig),distance,method="BFGS"))
  # }
   if(class(result)!="try-error"){
      M=length(initial_mu)
      params=result$par
      mu=params[1:M]
      sigma=params[(M+1):(2*M)]
      v=find_weights(mu,sqrt(sigma^2+h^2),m,h)[[1]]   # Parandus, 18.10.2022
      iterations=result$counts[1]
      error=FALSE
      value_of_objective_fn=result$value
   } else {
      mu=rep(0,M)
      sigma=rep(0,M)
      v<-rep(0,M)
      iterations=0
       error=TRUE
       value_of_objective_fn=-1
   }
  return(list(mu=mu,v=v,sigma=sigma,iterations=iterations,error=error,value_of_objective_fn=-value_of_objective_fn))
}



gen_data=function(n,mu,sigma,v,seed=NULL){
  #data generatsion from given mixture distribution
  if(!is.null(seed)){
    set.seed(seed)
  }
  size=n #sample size
  component_number=sample.int(length(mu),size,prob=v,replace=TRUE)
  data=rep(NA,size)
  n_i=table(component_number)
  for(i in 1:length(mu)){
    data[component_number==i]=rnorm(n_i[i],mean=mu[i],sd=sigma[i])
  }
  return(data)  
}



compute_h1=function(n,sigma){
  #suggested rule
  return(max(sigma)/n^(0.2))
}

numerical_experiments=function(sample_sizes,repeats,methods,parameter_sets,dist_power=2,output_basename="experiment1",seeds=NULL,initial_values=2,gen_iv=gen_iv1, h_fn=compute_h1){
  #function for comparing methods by numerical experiments
  #sample_sizes - which sample sizes to consider
  #repeats - how many samples for the same parameter set and same sample size to generate
  #methods - list of methods to compare
  #parameter_sets - sets of parameters as a matrix or data frame, first M-1 weights, then M means and at the end M standard deviations
  #initial_values - how many different initial value sets to consider
  # gen_iv - function for generating/selecting initial values for parameters
  n_sets=ncol(parameter_sets)
  total_samples=n_sets*repeats
  if(length(seeds)!=total_samples){
    print("Generating new seeds")
    seeds=sample.int(10^9,total_samples)
  }
  total_size=length(sample_sizes)*repeats*length(methods)*n_sets*initial_values
  M=(nrow(parameter_sets)+1)/3
  column_names=c("size","sample_number","error","set_number","iv_set_number",paste0("v",1:M),paste0("mu",1:M),paste0("sigma",1:M),"time","iterations","value_of_objective_fn","seed")
  veerge=length(column_names)
  results=data.frame(matrix(NA,nrow=total_size,ncol=veerge))
  colnames(results)=column_names
  results$method=rep(" ",total_size)
  column_names=colnames(results)
  mu_location=which(column_names=="mu1")
  v_location=which(column_names=="v1")
  sigma_location=which(column_names=="sigma1")
  initial_value_sets=data.frame(matrix(NA,nrow=repeats*initial_values*length(sample_sizes)*n_sets,ncol=3*M+4))
  colnames(initial_value_sets)=c("size","set_number","sample_number","iv_set_number",paste0("v",1:M),paste0("mu",1:M),paste0("sigma",1:M))
  row=1
  iv_row=1
  method_names=names(methods)
  tmp=1:M
  for(n in sample_sizes){
    for(j in 1:n_sets){
      v_mu_sig=parameter_sets[,j]
      v=v_mu_sig[1:(M-1)]
      v=c(v,1-sum(v))
      mu=v_mu_sig[M-1+tmp]
      sigma=v_mu_sig[2*M-1+tmp]
      h=h_fn(n,sigma)
      for(j1 in 1:repeats){
        seed=seeds[(j-1)*repeats+j1]
        data=gen_data(n,mu,sigma,v,seed)
        for(iv_number in 1:initial_values){
              av=gen_iv(mu,sigma,v,data,iv_number)
              initial_value_sets[iv_row,]=c(n,j,j1,iv_number,av$v,av$mu,av$sigma)
              iv_row=iv_row+1 
            for(j2 in 1:length(methods)){
              method=methods[[j2]]
              time=system.time(estimates<-method(data,h,av$mu,av$sigma,av$v[-M],dist_power))[1]
              results$time[row]=time
              results$size[row]=n
              results$sample_number[row]=j1
              results$error[row]=estimates$error
              results$set_number[row]=j
              results$iv_set_number[row]=iv_number
              results$iterations[row]=estimates$iterations
              results[row,mu_location+tmp-1]=estimates$mu
              results[row,v_location+tmp-1]=estimates$v
              results[row,sigma_location+tmp-1]=estimates$sigma
              results$method[row]=method_names[j2]
              results$value_of_objective_fn[row]=estimates$value_of_objective_fn
              results$seed[row]=seed
              print(paste0(row,"/",total_size))
              row=row+1
            }
        }
      }
    }
  write.csv(filter(results,size==n),file=paste0(output_basename,"_",n,".csv"),row.names=FALSE)
  write.csv(filter(initial_value_sets,size==n),file=paste0(output_basename,"_initial_values_",n,".csv"),row.names=FALSE)
  }
  saveRDS(object = seeds,file=paste0(output_basename,"_seeds.rds"))
  return(list(results=results,initial_value_sets=initial_value_sets,seeds=seeds))
}




#example usage
param2set=matrix(c(0.75,-1,0,1,1,0.5,-1,0,1,2,0.75,-1,0,2,1,0.5,1,0,1,1,0.5,0,0,1,2,0.5,-1.5,0,1,1,0.5,-5,0,1,2),ncol=7,nrow=5)

#
m1init<-matrix(0,nrow=8,ncol=6)
m1init[1,]<-c(0.5,0.5,-5,0,1,2)
m1init[2,]<-c(0.5,0.5,-4,0.75,1.25,1.75)
m1init[3,]<-c(0.5,0.5,-3,1,1.5,1.5)
m1init[4,]<-c(0.5,0.5,-6,1,1,2)
m1init[5,]<-c(0.6,0.4,-4,0.75,1.25,1.75)
m1init[6,]<-c(0.7,0.3,-3,0.5,1,1.75)
m1init[7,]<-c(0.4,0.6,-6,-1,1,2.25)
m1init[8,]<-c(0.4,0.6,-5,2,1.25,1.75)

j=4
gen_iv=function(mu,sigma,v,data,j){
  v_initial=m1init[j,1:2]
  mu_initial=m1init[j,3:4]
  sigma_initial=m1init[j,5:6]
  return(list(v=v_initial,mu=mu_initial,sigma=sigma_initial))
}


methods=list(EM=EM,dude=L2_distanceBFGS,dudec=L2_distance_corrected,pseudo=pseudo_likelihood,pseudoc=pseudo_likelihood_corrected)

true_params=param2set[,6:7]


first_try=numerical_experiments(sample_sizes=c(100,200),repeats=3,methods=methods,parameter_sets=true_params,dist_power=2,output_basename="test1",seeds=NULL,initial_values=8,gen_iv=gen_iv,h_fn=compute_h1)


