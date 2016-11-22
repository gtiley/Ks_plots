##################################################################
# George P. Tiley - gtiley@ufl.edu
# 20151107
# Code borrowed extensively from source code by Sally Otto (unpublished)
# and mixtools R package (Tatiana Benaglia, Didier Chauveau, David R. Hunter, Derek Young (2009). mixtools: An R Package for Analyzing Finite Mixture Models. Journal of Statistical Software, 32(6), 1-29. URL http://www.jstatsoft.org/v32/i06/.)
# 20160330
##################################################################
##################################################################
#model 1 = exp + (k-1) normals
#model 2 = exp + (k-1) gammas
#model 3 = k normals
#model 4 = k gammas
##################################################################

##########################################################################################
# initMix() - propose starting parameters for expectation maximization algorithm
#
#
##########################################################################################
initMix <- function(x, lambda = NULL, alpha = NULL, beta = NULL, k = NULL, model=NULL)
{
	n <- length(x)
    if (is.null(k)) 
    {
    	stop("Missing number of components ... exiting!")
    }
    else
    {   
    	if (model == 1) 
    	{  
    		if (k == 1) #fit a single exponential
	    	{
    			x.bar=mean(x)
				alpha = 1
				beta = x.bar
				lambda = 1
#			cat ("starting params:", lambda, "\t", alpha, "\t", beta, "\n")
	    	}
		    else
    		{
	    		if (is.null(lambda))
    			{
	    			lambda = runif(k)
		    	    lambda = lambda/sum(lambda)
    				x.sort=sort(x)
					ind=floor(n*cumsum(lambda))
					x.part=list()
					x.part[[1]]=x.sort[1:(ind[1]+1)]
					for(j in 2:k)
					{
						x.part[[j]]=x.sort[ind[j-1]:ind[j]]
					}
					x.bar=sapply(x.part,mean)
					x.sd=sapply(x.part,sd)
	    			if(is.null(alpha))
				    {
						alpha=x.bar
						alpha[1] = 1
#						cat ("Checking starting alphas:\t", alpha, "\n")
				    }
					if(is.null(beta))
				    {
						beta=x.sd
						beta[1] = 1/x.bar[1]
#						cat ("Checking starting betas:\t", beta, "\n")
		    		}
			    }
			    else
			    {
    				x.sort=sort(x)
					ind=floor(n*cumsum(lambda))
					x.part=list()
					x.part[[1]]=x.sort[1:(ind[1]+1)]
					for(j in 2:k)
					{
						x.part[[j]]=x.sort[ind[j-1]:ind[j]]
					}
					x.bar=sapply(x.part,mean)
					x.sd=sapply(x.part,sd)
	    			if(is.null(alpha))
				    {
						alpha=x.bar
						alpha[1] = 1
#						cat ("Checking starting alphas:\t", alpha, "\n")
				    }
					if(is.null(beta))
				    {
						beta=x.sd
						beta[1] = x.bar[1]
#						cat ("Checking starting betas:\t", beta, "\n")
		    		}
		    	}
			}
		}
		if (model == 2)
		{
			if (k == 1) #fit a single exponential
	    	{
    			x.bar=mean(x)
				alpha = 1
				beta = x.bar
				lambda = 1
#			cat ("starting params:", lambda, "\t", alpha, "\t", beta, "\n")
	    	}
		    else
    		{
	    		if (is.null(lambda))
    			{
	    			lambda = runif(k)
		    	    lambda = lambda/sum(lambda)
    				x.sort=sort(x)
					ind=floor(n*cumsum(lambda))
					x.part=list()
					x.part[[1]]=x.sort[1:(ind[1]+1)]
					for(j in 2:k)
					{
						x.part[[j]]=x.sort[ind[j-1]:ind[j]]
					}
					x.bar=sapply(x.part,mean)
					x2.bar=sapply(lapply(x.part,"^",2),mean)
					x.var=sapply(x.part,var)
	    			if(is.null(alpha))
				    {
						alpha=x.bar^2/(x2.bar-x.bar^2)
						alpha[1] = 1
#						cat ("Checking starting alphas:\t", alpha, "\n")
				    }
					if(is.null(beta))
				    {
						beta=(x2.bar-x.bar^2)/x.bar
						beta[1] = x.bar[1]
#						cat ("Checking starting betas:\t", beta, "\n")
						
		    		}
			    }
			    else
			    {
			    	x.sort=sort(x)
					ind=floor(n*cumsum(lambda))
					x.part=list()
					x.part[[1]]=x.sort[1:(ind[1]+1)]
					for(j in 2:k)
					{
						x.part[[j]]=x.sort[ind[j-1]:ind[j]]
					}
					x.bar=sapply(x.part,mean)
					x2.bar=sapply(lapply(x.part,"^",2),mean)
					x.var=sapply(x.part,var)
	    			if(is.null(alpha))
				    {
						alpha=x.bar^2/(x2.bar-x.bar^2)
						alpha[1] = 1
#						cat ("Checking starting alphas:\t", alpha, "\n")
				    }
					if(is.null(beta))
				    {
						beta=(x2.bar-x.bar^2)/x.bar
						beta[1] = x.bar[1]
#						cat ("Checking starting betas:\t", beta, "\n")
		    		}
		    	}
			}
		}
		
		if (model == 3)
		{
			if (k == 1) #fit a single exponential
	    	{
    			x.bar=mean(x)
    			x.sd=sd(x)
				alpha = x.bar
				beta = x.sd
				lambda = 1
#				cat ("starting params:", lambda, "\t", alpha, "\t", beta, "\n")
	    	}
	    	else
	    	{
				if (is.null(lambda))
    			{
	    			lambda = runif(k)
		   	    	lambda = lambda/sum(lambda)
    				x.sort=sort(x)
					ind=floor(n*cumsum(lambda))
					x.part=list()
					x.part[[1]]=x.sort[1:(ind[1]+1)]
					for(j in 2:k)
					{
						x.part[[j]]=x.sort[ind[j-1]:ind[j]]
					}
					x.bar=sapply(x.part,mean)
					x.sd=sapply(x.part,sd)
    				if(is.null(alpha))
				    {
						alpha=x.bar
#						cat ("Checking starting alphas:\t", alpha, "\n")
				    }
					if(is.null(beta))
				    {
						beta=x.sd
#						cat ("Checking starting betas:\t", beta, "\n")
	    			}
			    }		
				else
				{
					x.sort=sort(x)
					ind=floor(n*cumsum(lambda))
					x.part=list()
					x.part[[1]]=x.sort[1:(ind[1]+1)]
					for(j in 2:k)
					{
						x.part[[j]]=x.sort[ind[j-1]:ind[j]]
					}
					x.bar=sapply(x.part,mean)
					x.sd=sapply(x.part,sd)
		    		if(is.null(alpha))
					{
						alpha=x.bar
#						cat ("Checking starting alphas:\t", alpha, "\n")
					}
					if(is.null(beta))
					{
						beta=x.sd
#						cat ("Checking starting betas:\t", beta, "\n")
		    		}
			   	}	
			}
		}	
		if (model == 4)
		{
			if (k==1)
			{
				x.bar = mean(x)
				x2.bar = mean(x * x)
				x.var = var(x)
				alpha = x.bar^2/(x2.bar-x.bar^2)
				beta = (x2.bar-x.bar^2)/x.bar
				lambda = 1
#				cat ("starting params:", lambda, "\t", alpha, "\t", beta, "\n")
			}
			else
			{
		    	if (is.null(lambda))
    			{
	    			lambda = runif(k)
		    	    lambda = lambda/sum(lambda)
    				x.sort=sort(x)
					ind=floor(n*cumsum(lambda))
					x.part=list()
					x.part[[1]]=x.sort[1:(ind[1]+1)]
					for(j in 2:k)
					{
						x.part[[j]]=x.sort[ind[j-1]:ind[j]]
					}
					x.bar=sapply(x.part,mean)
					x2.bar=sapply(lapply(x.part,"^",2),mean)
					x.var=sapply(x.part,var)
		    		if(is.null(alpha))
				    {
						alpha=x.bar^2/(x2.bar-x.bar^2)
#						cat ("Checking starting alphas:\t", alpha, "\n")
				    }
					if(is.null(beta))
				    {
						beta=(x2.bar-x.bar^2)/x.bar
#						cat ("Checking starting betas:\t", beta, "\n")
			    	}
				}
				else
				{
				   	x.sort=sort(x)
					ind=floor(n*cumsum(lambda))
					x.part=list()
					x.part[[1]]=x.sort[1:(ind[1]+1)]
					for(j in 2:k)
					{
						x.part[[j]]=x.sort[ind[j-1]:ind[j]]
					}
					x.bar=sapply(x.part,mean)
					x2.bar=sapply(lapply(x.part,"^",2),mean)
					x.var=sapply(x.part,var)
	    			if(is.null(alpha))
					{
						alpha=x.bar^2/(x2.bar-x.bar^2)
						alpha[1] = 1
#						cat ("Checking starting alphas:\t", alpha, "\n")
					}
					if(is.null(beta))
					{
						beta=(x2.bar-x.bar^2)/x.bar
#						cat ("Checking starting betas:\t", beta, "\n")
			    	}
		    	}
			}	
		}
	}
	list(lambda=lambda, alpha=alpha, beta=beta, k=k)
}

##########################################################################################
# mixEM() - optimize model parameters by expectation maximization
#
#
##########################################################################################
mixEM <- function (x, lambda = NULL, alpha = NULL, beta = NULL, k = NULL, model = 1, nstarts = 100, epsilon = 1e-08, maxit = 1000, maxrestarts=20, verb=0)
{
	bestlnL = -1000000000
	w = 0
	for (w in 1:nstarts)
	{
	
		check = 0
	
		x <- as.vector(x)
		lambda = NULL
		alpha = NULL
		beta = NULL
		tmp <- initMix(x = x, lambda = lambda, alpha=alpha, beta=beta, k=k, model=model)
		lambda <- tmp$lambda 
		alpha <- tmp$alpha
		beta <- tmp$beta
		theta <- c(alpha, beta)
		k <- tmp$k 	
		iter <- 0
		mr <- 0
		diff <- epsilon+1
		n <- length(x)
#		cat ("starting params:", lambda, "\t", alpha, "\t", beta, "\n")
		
############################## 
	# Cases for all k > 0
	# The first component is an expoential distribution while the rest are normals, lognormals, or gammas, depending on the model selected
	# flexible numbers of k
############################## 

		if (model == 1)
		{
			dens <- NULL
			dens <- function(lambda, theta, k)
			{
				alpha=theta[1:k]
				beta=theta[(k+1):(2*k)]
				psum = 0
				temp<-NULL			
				for(j in 1:k)
				{
					if (j == 1)
					{
						temp=cbind(temp,dexp(x,rate = 1/beta[j])) 
					}
					if (j > 1)
					{
						temp=cbind(temp,dnorm(x, mean = alpha[j], sd = beta[j])) 
					} 
					if (j < k)
					{
						psum = psum + lambda[j]
					}
					if (j == k)
					{
						lambda[j] = (1 - psum)
					}
				}
				temp=t(lambda*t(temp))
				temp			
			}	
		}
		if (model == 2)
		{	
			dens <- NULL
			dens <- function(lambda, theta, k)
			{
				alpha=theta[1:k]
				beta=theta[(k+1):(2*k)]
				psum = 0
				temp<-NULL			
				for(j in 1:k)
				{
					if (j == 1)
					{
						temp=cbind(temp,dexp(x,rate = 1/beta[j])) 
					}
					if (j > 1)
					{
						temp=cbind(temp,dgamma(x, shape = alpha[j], scale = 1/beta[j])) 
					} 
					if (j < k)
					{
						psum = psum + lambda[j]
					}
					if (j == k)
					{
						lambda[j] = (1 - psum)
					}
				}
				temp=t(lambda*t(temp))
				temp			
			}
		}
		if (model == 3)
		{
			dens <- NULL
			dens <- function(lambda, theta, k)
			{
				alpha=theta[1:k]
				beta=theta[(k+1):(2*k)]
				psum = 0
				temp<-NULL			
				for(j in 1:k)
				{
					temp=cbind(temp,dnorm(x, mean = alpha[j], sd = beta[j])) 
					if (j < k)
					{
						psum = psum + lambda[j]
					}
					if (j == k)
					{
						lambda[j] = (1 - psum)
					}
				}
				temp=t(lambda*t(temp))
				temp			
			}	
		}
		if (model == 4)
		{	
			dens <- NULL
			dens <- function(lambda, theta, k)
			{
				alpha=theta[1:k]
				beta=theta[(k+1):(2*k)]
				psum = 0
				temp<-NULL			
				for(j in 1:k)
				{
					temp=cbind(temp,dgamma(x, shape = alpha[j], scale = 1/beta[j])) 
					if (j < k)
					{
						psum = psum + lambda[j]
					}
					if (j == k)
					{
						lambda[j] = (1 - psum)
					}
				}
				temp=t(lambda*t(temp))
				temp			
			}
		}
	
		old.obs.ll <- sum(log(apply(dens(lambda, theta, k),1,sum)))
		ll <- old.obs.ll	
			
		model.ll <- function(theta,z,lambda,k) 
		{
			-sum(z*log(dens(lambda,theta,k)))
		}
		
		while(diff > epsilon && iter < maxit && check == 0 && is.na(diff)==FALSE)
		{
#update probabilities from previous optimization		
			dens1=dens(lambda,theta,k)
			z=dens1/apply(dens1,1,sum)
			lambda.hat=apply(z,2,mean)
				
#nlm -> newtonian minimization with full-blown 2nd order partials; set Hessian = TRUE for quasi-newton
			out=try(suppressWarnings(nlm(model.ll,p=theta,lambda=lambda.hat,k=k,z=z)),silent=FALSE)
			if(class(out)=="try-error")
			{
			
				if(mr==maxrestarts) 
				{
					check = 1
				}
				
				mr <- mr+1
				lambda = NULL
				alpha = NULL
				beta = NULL
				tmp <- initMix(x = x, lambda = lambda, alpha=alpha, beta=beta, k=k, model=model)
	 			lambda <- tmp$lambda 
				alpha <- tmp$alpha
 				beta <- tmp$beta
	  			theta <- c(alpha,beta)
	 			k <- tmp$k 
		  		iter <- 0
 				diff <- epsilon+1
				old.obs.ll <- sum(log(apply(dens(lambda, theta, k),1,sum)))
   				ll <- old.obs.ll
			} 
			else
  			{
				theta.hat=out$estimate
				alpha.hat=theta.hat[1:k]
				beta.hat=theta.hat[(k+1):(2*k)]
				new.obs.ll <- sum(log(apply(dens(lambda.hat, theta.hat, k),1,sum)))
				diff <- new.obs.ll-old.obs.ll
				old.obs.ll <- new.obs.ll
				ll <- c(ll,old.obs.ll)
				lambda=lambda.hat
				theta=theta.hat
				alpha=alpha.hat
				beta=beta.hat
				iter=iter+1
#				cat (diff,"\n")


##############################################################
#Additional logic block to prevent halting
				if (is.nan(new.obs.ll) == TRUE)
				{
			        lambda = NULL
					alpha = NULL
					beta = NULL
               		tmp <- initMix(x = x, lambda = lambda, alpha=alpha, beta=beta, k=k, model=model)
                    lambda <- tmp$lambda
   	                alpha <- tmp$alpha
       	            beta <- tmp$beta
           	        theta <- c(alpha,beta)
              	    k <- tmp$k
                    iter <- 0
                    diff <- epsilon+1
                    old.obs.ll <- sum(log(apply(dens(lambda, theta, k),1,sum)))
                    ll <- old.obs.ll
				}
##############################################################
        		if (verb == 1) 
	      		{
    	       		cat("iteration =", iter, " log-lik diff =", diff, " log-lik =", new.obs.ll, "\n")
    	   	   		cat("#=====Parameter Estimates=====#\nlambda:\t", lambda, "\nalpha:\t", alpha, "\nbeta:\t", beta, "\n#============================#\n")
      	  		}	   	
			}
    	}
# print optimized values/ save to object
		if (iter == maxit) 
		   {
    	   	cat("WARNING! NOT CONVERGENT!", "\n")
	    }
	    
	    if (check == 0)
		{
#	    	cat("number of iterations=", iter, "\n")
			theta=rbind(alpha,beta)
			rownames(theta)=c("alpha","beta")
			colnames(theta)=c(paste("comp", ".", 1:k, sep = ""))
	    	a=list(x=x, lambda = lambda, parameters = theta, loglik = new.obs.ll, posterior = z, all.loglik=ll, ft="mixEM")
			class(a) = "mixEM"
    		if (new.obs.ll > bestlnL && is.finite(new.obs.ll) == TRUE)
			{
   				bestlnL = new.obs.ll
			    result = list(x=x, lambda = lambda, parameters = theta, loglik = new.obs.ll, posterior = z, all.loglik=ll, ft="mixEM")
		   	}		
		}	
	}
	result
}

##########################################################################################
# fitMix() - calculate lnL for mixture models with fixed parameters
#
#
##########################################################################################	
fitMix <- function (x, lambda = NULL, alpha = NULL, beta = NULL, k = NULL, model = NULL)
{
	if (is.null(k) | is.null(lambda) | is.null(alpha) | is.null(beta) | is.null(model)) 
    {
    	stop("Missing one or more necessary parameters ... exiting!")
    }
    else
    {
    	theta <- c(alpha, beta)
    	if (model == 1)
		{
			dens <- NULL
			dens <- function(lambda, theta, k)
			{
				alpha=theta[1:k]
				beta=theta[(k+1):(2*k)]
				psum = 0
				temp<-NULL			
				for(j in 1:k)
				{
					if (j == 1)
					{
						temp=cbind(temp,dexp(x,rate = 1/beta[j])) 
					}
					if (j > 1)
					{
						temp=cbind(temp,dnorm(x, mean = alpha[j], sd = beta[j])) 
					} 
					if (j < k)
					{
						psum = psum + lambda[j]
					}
					if (j == k)
					{
						lambda[j] = (1 - psum)
					}
				}
				temp=t(lambda*t(temp))
				temp			
			}	
		}
		if (model == 2)
		{	
			dens <- NULL
			dens <- function(lambda, theta, k)
			{
				alpha=theta[1:k]
				beta=theta[(k+1):(2*k)]
				psum = 0
				temp<-NULL			
				for(j in 1:k)
				{
					if (j == 1)
					{
						temp=cbind(temp,dexp(x,rate = 1/beta[j])) 
					}
					if (j > 1)
					{
						temp=cbind(temp,dgamma(x, shape = alpha[j], scale = 1/beta[j])) 
					} 
					if (j < k)
					{
						psum = psum + lambda[j]
					}
					if (j == k)
					{
						lambda[j] = (1 - psum)
					}
				}
				temp=t(lambda*t(temp))
				temp			
			}
		}
		if (model == 3)
		{
			dens <- NULL
			dens <- function(lambda, theta, k)
			{
				alpha=theta[1:k]
				beta=theta[(k+1):(2*k)]
				psum = 0
				temp<-NULL			
				for(j in 1:k)
				{
					temp=cbind(temp,dnorm(x, mean = alpha[j], sd = beta[j])) 
					if (j < k)
					{
						psum = psum + lambda[j]
					}
					if (j == k)
					{
						lambda[j] = (1 - psum)
					}
				}
				temp=t(lambda*t(temp))
				temp			
			}	
		}
		if (model == 4)
		{	
			dens <- NULL
			dens <- function(lambda, theta, k)
			{
				alpha=theta[1:k]
				beta=theta[(k+1):(2*k)]
				psum = 0
				temp<-NULL			
				for(j in 1:k)
				{
					temp=cbind(temp,dgamma(x, shape = alpha[j], scale = 1/beta[j])) 
					if (j < k)
					{
						psum = psum + lambda[j]
					}
					if (j == k)
					{
						lambda[j] = (1 - psum)
					}
				}
				temp=t(lambda*t(temp))
				temp			
			}
		}
		
		obs.ll <- sum(log(apply(dens(lambda, theta, k),1,sum)))
		dens1=dens(lambda,theta,k)
		z=dens1/apply(dens1,1,sum)
	    theta=rbind(alpha,beta)
		rownames(theta)=c("alpha","beta")
		colnames(theta)=c(paste("comp", ".", 1:k, sep = ""))
	    a=list(x=x, lambda = lambda, parameters = theta, loglik = obs.ll, posterior = z, ft="mixEM")
		class(a) = "mixEM"
	   	a
    }
}

##########################################################################################
# fixComponents() - calculate lnL for mixture models with fixed lambda, but free alphas/betas
#
#
##########################################################################################	
fitComponents <- function (x, lambda = NULL, alpha = NULL, beta = NULL, k = NULL, model = 1, nstarts = 100, epsilon = 1e-08, maxit = 1000, maxrestarts=20, verb=0)
{
	if (is.null(k) | is.null(lambda) | is.null(model)) 
    {
    	stop("The fixed component proportions must be defined with lambda=c(p1,p2,...,pn)\n")
 	}   
 	else
 	{
 		bestlnL = -1000000000
 		w=0;
		for (w in 1:nstarts)
		{
		
			check = 0
			
			x <- as.vector(x)
			lambda = lambda
			alpha = NULL
			beta = NULL
			tmp <- initMix(x = x, lambda = lambda, alpha=alpha, beta=beta, k=k, model=model)
			alpha <- tmp$alpha
			beta <- tmp$beta
			theta <- c(alpha, beta)
			k <- tmp$k 	
			iter <- 0
			mr <- 0
			diff <- epsilon+1
			n <- length(x)
 	
 			if (model == 1)
			{
				dens <- NULL
				dens <- function(lambda, theta, k)
				{
					alpha=theta[1:k]
					beta=theta[(k+1):(2*k)]
					psum = 0
					temp<-NULL			
					for(j in 1:k)
					{
						if (j == 1)
						{
							temp=cbind(temp,dexp(x,rate = 1/beta[j])) 
						}
						if (j > 1)
						{
							temp=cbind(temp,dnorm(x, mean = alpha[j], sd = beta[j])) 
						} 
						if (j < k)
						{
							psum = psum + lambda[j]
						}
						if (j == k)
						{
							lambda[j] = (1 - psum)
						}
					}
					temp=t(lambda*t(temp))
					temp			
				}	
			}
			if (model == 2)
			{	
				dens <- NULL
				dens <- function(lambda, theta, k)
				{
					alpha=theta[1:k]
					beta=theta[(k+1):(2*k)]
					psum = 0
					temp<-NULL			
					for(j in 1:k)
					{
						if (j == 1)
						{
							temp=cbind(temp,dexp(x,rate = 1/beta[j])) 
						}
						if (j > 1)
						{
							temp=cbind(temp,dgamma(x, shape = alpha[j], scale = 1/beta[j])) 
						} 
						if (j < k)
						{
							psum = psum + lambda[j]
						}
						if (j == k)
						{
							lambda[j] = (1 - psum)
						}
					}
					temp=t(lambda*t(temp))
					temp			
				}
			}
			if (model == 3)
			{
				dens <- NULL
				dens <- function(lambda, theta, k)
				{
					alpha=theta[1:k]
					beta=theta[(k+1):(2*k)]
					psum = 0
					temp<-NULL			
					for(j in 1:k)
					{
						temp=cbind(temp,dnorm(x, mean = alpha[j], sd = beta[j])) 
						if (j < k)
						{
							psum = psum + lambda[j]
						}
						if (j == k)
						{
							lambda[j] = (1 - psum)
						}
					}
					temp=t(lambda*t(temp))
					temp			
				}	
			}
			if (model == 4)
			{	
				dens <- NULL
				dens <- function(lambda, theta, k)
				{
					alpha=theta[1:k]
					beta=theta[(k+1):(2*k)]
					psum = 0
					temp<-NULL			
					for(j in 1:k)
					{
						temp=cbind(temp,dgamma(x, shape = alpha[j], scale = 1/beta[j])) 
						if (j < k)
						{
							psum = psum + lambda[j]
						}
						if (j == k)
						{
							lambda[j] = (1 - psum)
						}
					}
					temp=t(lambda*t(temp))
					temp			
				}
			}

			old.obs.ll <- sum(log(apply(dens(lambda, theta, k),1,sum)))
			ll <- old.obs.ll	
			
			model.ll <- function(theta,z,lambda,k) 
			{
				-sum(z*log(dens(lambda,theta,k)))
			}
	
			while(diff > epsilon && iter < maxit && check == 0 && is.na(diff)==FALSE)
			{
#update probabilities from previous optimization		
				dens1=dens(lambda,theta,k)
				z=dens1/apply(dens1,1,sum)
############
#Don't need to update lambda for this optimization, fixed in this function
#				lambda.hat=apply(z,2,mean)
				
#nlm -> newtonian minimization with full-blown 2nd order partials; set Hessian = TRUE for quasi-newton
				out=try(suppressWarnings(nlm(model.ll,p=theta,lambda=lambda,k=k,z=z)),silent=FALSE)
				if(class(out)=="try-error")
				{
				
					if(mr==maxrestarts) 
					{
						check = 1
					}

					mr <- mr+1
	#				tmp <- initMix(x = x, lambda = lambda, alpha=alpha, beta=beta, k=k, model=model)
	# 				lambda <- tmp$lambda 
	#				alpha <- tmp$alpha
	# 				beta <- tmp$beta
	#				alpha = sapply(alpha)
	#				beta = sapply(beta)
					for (i in 1:k)
					{
						alpha[i] = alpha[i] + 0.05
						beta[i] = beta[i] + 0.05
						if (model == 1 | model == 2 | model == 3)
						{
							alpha[1] = 1
						}
					}
	  				theta <- c(alpha,beta)
# 					k <- tmp$k 
	  				iter <- 0
	 				diff <- epsilon+1
					old.obs.ll <- sum(log(apply(dens(lambda, theta, k),1,sum)))
   					ll <- old.obs.ll
				} 
	  			else
		  		{
					theta.hat=out$estimate
					alpha.hat=theta.hat[1:k]
					beta.hat=theta.hat[(k+1):(2*k)]
					new.obs.ll <- sum(log(apply(dens(lambda, theta.hat, k),1,sum)))
					diff <- new.obs.ll-old.obs.ll
					old.obs.ll <- new.obs.ll
					ll <- c(ll,old.obs.ll)
#				lambda=lambda.hat
					theta=theta.hat
					alpha=alpha.hat
					beta=beta.hat
					iter=iter+1
        			if (verb == 1) 
		       		{
    		       		cat("iteration =", iter, " log-lik diff =", diff, " log-lik =", new.obs.ll, "\n")
    		       		cat("#=====Parameter Estimates=====#\nlambda:\t", lambda, "\nalpha:\t", alpha, "\nbeta:\t", beta, "\n#============================#\n")
	       		  	}	   	
				}
    		}
# print optimized values/ save to object
			if (iter == maxit) 
		    {
    		   	cat("WARNING! NOT CONVERGENT!", "\n")
		    }
		    
		    if (check == 0)
		    {
			    cat("number of iterations=", iter, "\n")
				theta=rbind(alpha,beta)
				rownames(theta)=c("alpha","beta")
				colnames(theta)=c(paste("comp", ".", 1:k, sep = ""))
	    		a=list(x=x, lambda = lambda, parameters = theta, loglik = new.obs.ll, posterior = z, all.loglik=ll, ft="mixEM")
				class(a) = "mixEM"
			    if (new.obs.ll > bestlnL && is.finite(new.obs.ll) == TRUE)
				{
    				bestlnL = new.obs.ll
				    result = list(x=x, lambda = lambda, parameters = theta, loglik = new.obs.ll, posterior = z, all.loglik=ll, ft="mixEM")
		   		}
		   	}	
		}
		result
	}	
}