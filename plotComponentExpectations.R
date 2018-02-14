plotComponents <- function (x, lambda, params, xlim = NULL, ylim = NULL, xlab= NULL, ylab = NULL, main = NULL, model)
{
	if (is.null(x) || is.null(lambda) || is.null(params) || is.null(model))
	{
		stop ("Missing some required stuff: data, mixing proportions, parameters, model specification...exiting\n")
	}
	else
	{
		n = length(x)
		k = length(lambda)
		alpha <- c()
		beta <- c()
		for (i in 1:k)
		{
			alpha[i] = params[2*i - 1]
			beta[i] = params[2*i]
		}
		
		if (is.null(ylab))
		{
			ylab = ""
		}
		if (is.null(xlab))
		{
			xlab = ""
		}
		if (is.null(main))
		{
			main = ""
		}
		if (is.null(xlim))
		{
			xmax <- max(x)
			xlim = c(0,xmax)
		}
		if (is.null(ylim))
		{
			Denmax = 0
			for (i in 1:k)
			{
				if (model == 1)
				{
					if (i == 1)
					{
						theoreticalData <- rexp(round(lambda[i] * n), rate = 1/beta[i])
						theoreticalDensity <- dexp(theoreticalData, rate = 1/beta[i])
						for (w in 1:length(theoreticalDensity))
						{
							if (theoreticalDensity[w] > Denmax)
							{
								Denmax = theoreticalDensity[w]
							}
						}
					}
					else
					{
						theoreticalData <- rnorm(round(lambda[i] * n), mean = alpha[i], sd = beta[i])
						theoreticalDensity <- dnorm(theoreticalData, mean = alpha[i], sd = beta[i])
						for (w in 1:length(theoreticalDensity))
						{
							if (theoreticalDensity[w] > Denmax)
							{
								Denmax = theoreticalDensity[w]
							}
						}
					}
				}
				if (model == 2)
				{
					if (i == 1)
					{
						theoreticalData <- rexp(round(lambda[i] * n), rate = 1/beta[i])
						theoreticalDensity <- dexp(theoreticalData, rate = 1/beta[i])
						for (w in 1:length(theoreticalDensity))
						{
							if (theoreticalDensity[w] > Denmax)
							{
								Denmax = theoreticalDensity[w]
							}
						}
					}
					else
					{
						theoreticalData <- rgamma(round(lambda[i] * n), shape = alpha[i], scale = 1/beta[i])
						theoreticalDensity <- dgamma(theoreticalData, shape = alpha[i], scale = 1/beta[i])
						for (w in 1:length(theoreticalDensity))
						{
							if (theoreticalDensity[w] > Denmax)
							{
								Denmax = theoreticalDensity[w]
							}
						}
					}
				}
				if (model ==3)
				{
					if (i == 1)
					{
						theoreticalData <- rexp(round(lambda[i] * n), rate = 1/beta[i])
						theoreticalDensity <- dexp(theoreticalData, rate = 1/beta[i])
						for (w in 1:length(theoreticalDensity))
						{
							if (theoreticalDensity[w] > Denmax)
							{
								Denmax = theoreticalDensity[w]
							}
						}
					}
					else
					{
						theoreticalData <- rlnorm(round(lambda[i] * n), meanlog = alpha[i], sdlog = beta[i])
						theoreticalDensity <- dlnorm(theoreticalData, meanlog = alpha[i], sdlog = beta[i])
						for (w in 1:length(theoreticalDensity))
						{
							if (theoreticalDensity[w] > Denmax)
							{
								Denmax = theoreticalDensity[w]
							}
						}
					}				
				}
				if (model == 4)
				{
					theoreticalData <- rnorm(round(lambda[i] * n), mean = alpha[i], sd = beta[i])
					theoreticalDensity <- dnorm(theoreticalData, mean = alpha[i], sd = beta[i])
					for (w in 1:length(theoreticalDensity))
					{
						if (theoreticalDensity[w] > Denmax)
						{
							Denmax = theoreticalDensity[w]
						}
					}
				}
				if (model == 5)
				{
					theoreticalData <- rgamma(round(lambda[i] * n), shape = alpha[i], scale = 1/beta[i])
					theoreticalDensity <- dgamma(theoreticalData, shape = alpha[i], scale = 1/beta[i])
					for (w in 1:length(theoreticalDensity))
					{
						if (theoreticalDensity[w] > Denmax)
						{
							Denmax = theoreticalDensity[w]
						}
					}
				}	
				if (model == 6)
				{
					theoreticalData <- rlnorm(round(lambda[i] * n), meanlog = alpha[i], sdlog = beta[i])
					theoreticalDensity <- dlnorm(theoreticalData, meanlog = alpha[i], sdlog = beta[i])
					for (w in 1:length(theoreticalDensity))
					{
						if (theoreticalDensity[w] > Denmax)
						{
							Denmax = theoreticalDensity[w]
						}
					}
				}		
			}
			ylim = c(0,Denmax)	
		}
			
### Now actually begin the plotting		
		hist (x, breaks=100, col="royalblue2", xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,prob=TRUE)

		for (i in 1:k)
		{
			if (model == 1)
			{
				if (i == 1)
				{
					theoreticalData <- rexp(round(lambda[i] * n), rate = 1/beta[i])
					theoreticalDensity <- dexp(theoreticalData, rate = 1/beta[i])
					lines(smooth.spline(theoreticalData, theoreticalDensity))
					segments(beta[i],0,beta[i],dexp(beta[i],rate=1/beta[i]),col="firebrick",lwd=2)
				}
				else
				{
					theoreticalData <- rnorm(round(lambda[i] * n), mean = alpha[i], sd = beta[i])
					theoreticalDensity <- dnorm(theoreticalData, mean = alpha[i], sd = beta[i])
					lines(smooth.spline(theoreticalData, theoreticalDensity))
					segments(alpha[i],0,alpha[i],dnorm(alpha[i],mean=alpha[i], sd=beta[i]),col="firebrick",lwd=2)
				}
			}
			if (model == 2)
			{
				if (i == 1)
				{
					theoreticalData <- rexp(round(lambda[i] * n), rate = 1/beta[i])
					theoreticalDensity <- dexp(theoreticalData, rate = 1/beta[i])
					lines(smooth.spline(theoreticalData, theoreticalDensity))
					segments(beta[i],0,beta[i],dexp(beta[i],rate=1/beta[i]),col="firebrick",lwd=2)
				}
				else
				{
					theoreticalData <- rgamma(round(lambda[i] * n), shape = alpha[i], scale = 1/beta[i])
					theoreticalDensity <- dgamma(theoreticalData, shape = alpha[i], scale = 1/beta[i])
					lines(smooth.spline(theoreticalData, theoreticalDensity))
					segments(alpha[i]/beta[i],0,alpha[i]/beta[i],dgamma(alpha[i]/beta[i], shape = alpha[i], scale = 1/beta[i]),col="firebrick",lwd=2)
				}
			}
			if (model == 3)
			{
				if (i == 1)
				{
					theoreticalData <- rexp(round(lambda[i] * n), rate = 1/beta[i])
					theoreticalDensity <- dexp(theoreticalData, rate = 1/beta[i])
					lines(smooth.spline(theoreticalData, theoreticalDensity))
					segments(beta[i],0,beta[i],dexp(beta[i],rate=1/beta[i]),col="firebrick",lwd=2)
				}
				else
				{
					theoreticalData <- rlnorm(round(lambda[i] * n), meanlog = alpha[i], sdlog = beta[i])
					theoreticalDensity <- dlnorm(theoreticalData, meanlog = alpha[i], sdlog = beta[i])
					lines(smooth.spline(theoreticalData, theoreticalDensity))
					segments(exp(alpha[i]),0,exp(alpha[i]),dlnorm(alpha[i], meanlog = alpha[i], sdlog = beta[i]),col="firebrick",lwd=2)
				}
			}
			
##############
# Models 4-6 only have mixtures of one type of distribution
##############
			if (model == 4)
			{
				theoreticalData <- rnorm(round(lambda[i] * n), mean = alpha[i], sd = beta[i])
				theoreticalDensity <- dnorm(theoreticalData, mean = alpha[i], sd = beta[i])
				lines(smooth.spline(theoreticalData, theoreticalDensity))
				segments(alpha[i],0,alpha[i],dnorm(alpha[i],mean=alpha[i], sd=beta[i]),col="firebrick",lwd=2)
			}
			if (model == 5)
			{
				theoreticalData <- rgamma(round(lambda[i] * n), shape = alpha[i], scale = 1/beta[i])
				theoreticalDensity <- dgamma(theoreticalData, shape = alpha[i], scale = beta[i])
				lines(smooth.spline(theoreticalData, theoreticalDensity))
				segments(alpha[i] * beta[i],0,alpha[i] * beta[i],dgamma(alpha[i] * beta[i], shape = alpha[i], scale = beta[i]),col="firebrick",lwd=2, lty=2)
			}
			if (model == 6)
			{
				theoreticalData <- rlnorm(round(lambda[i] * n), meanlog = alpha[i], sdlog = beta[i])
				theoreticalDensity <- dlnorm(theoreticalData, meanlog = alpha[i], sdlog = beta[i])
				lines(smooth.spline(theoreticalData, theoreticalDensity))
				segments(exp(alpha[i]),0,exp(alpha[i]),dlnorm(alpha[i], meanlog = alpha[i], sdlog = beta[i]),col="firebrick",lwd=2)
			}			
		}
	}
}