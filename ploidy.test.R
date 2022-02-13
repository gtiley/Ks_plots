ploidy.test <- function (x, maxPloidy = 6, model = 4, nstarts = 100, outPrefix = NULL)
{
	n <- length(x)
	thisK <- 1
	resList <- c()
	bicList <- c()
	if (is.null(outPrefix))
	{
		outPrefix <- "ploidyResults"
	}
	if (maxPloidy > 8)
	{
		stop ("Currently not supporting ploidy levels greater than 8")
	}
	while(thisK <= (maxPloidy - 1))
	{
		if (thisK == 1)
		{
			obs.res.one <- onlyVariances(x,lambda=c(1),alpha=c(0.5),k=thisK,model=model,verb=0)
			bic.one <- -2 * obs.res.one$loglik + (thisK * 2 - 1) * log(n)
			resList[[thisK]] <- obs.res.one
			bicList[thisK] <- bic.one
		}
		if (thisK == 2)
		{
			obs.res.two <- onlyVariances(x,lambda=c(0.5,0.5),alpha=c(1/3,2/3),k=thisK,model=model,verb=0)
			bic.two <- -2 * obs.res.two$loglik + (thisK * 2 - 1) * log(n)
			resList[[thisK]] <- obs.res.two
			bicList[thisK] <- bic.two
		}
		if (thisK == 3)
		{
			obs.res.three <- onlyVariances(x,lambda=c(0.25,0.5,0.25),alpha=c(0.25,0.5,0.75),k=thisK,model=model,verb=0)
			bic.three <- -2 * obs.res.three$loglik + (thisK * 2 - 1) * log(n)
			resList[[thisK]] <- obs.res.three
			bicList[thisK] <- bic.three
		}
		if (thisK == 4)
		{
			obs.res.four <- onlyVariances(x,lambda=c(1/4,1/4,1/4,1/4),alpha=c(1/5,2/5,3/5,4/5),k=thisK,model=model,verb=0)
			bic.four <- -2 * obs.res.four$loglik + (thisK * 2 - 1) * log(n)
			resList[[thisK]] <- obs.res.four
			bicList[thisK] <- bic.four
		}
		if (thisK == 5)
		{
			obs.res.five <- onlyVariances(x,lambda=c(1/6,1/6,2/6,1/6,1/6),alpha=c(1/6,2/6,3/6,4/6,5/6),k=thisK,model=model,verb=0)
			bic.five <- -2 * obs.res.five$loglik + (thisK * 2 - 1) * log(n)
			resList[[thisK]] <- obs.res.five
			bicList[thisK] <- bic.five
		}
		if (thisK == 6)
		{
			obs.res.six <- onlyVariances(x,lambda=c(1/6,1/6,1/6,1/6,1/6,1/6),alpha=c(1/7,2/7,3/7,4/7,5/7,6/7),k=thisK,model=model,verb=0)
			bic.six <- -2 * obs.res.six$loglik + (thisK * 2 - 1) * log(n)
			resList[[thisK]] <- obs.res.six
			bicList[thisK] <- bic.six
		}
		if (thisK == 7)
		{
			obs.res.seven <- onlyVariances(x,lambda=c(1/8,1/8,1/8,2/8,1/8,1/8,1/8),alpha=c(1/8,2/8,3/8,4/8,5/8,6/8,7/8),k=thisK,model=model,verb=0)
			bic.seven <- -2 * obs.res.seven$loglik + (thisK * 2 - 1) * log(n)
			resList[[thisK]] <- obs.res.seven
			bicList[thisK] <- bic.seven
		}
		thisK = thisK + 1
	}
			
	minBIC <- 0
	expDeltaList <- c()
	expDeltaSum <- 0
	modelWeights <- c()
	modelWinner <- 0
		
	for (i in 1:(maxPloidy-1))
	{
		if (bicList[i] < minBIC)
		{
			minBIC <- bicList[i]
		}
	}
	
	for (i in 1:(maxPloidy-1))
	{
		expDeltaList[i] <- exp(minBIC - bicList[i])
		expDeltaSum <- expDeltaSum + expDeltaList[i]
	}
	
	for (i in 1:(maxPloidy-1))
	{
		modelWeights[i] <- expDeltaList[i]/expDeltaSum
		if (modelWeights[i] > 0.5)
		{
			modelWinner <- i
			
			
			sink(paste(outPrefix,".modelSelection.txt",sep=""))
			cat("Model Choice ", model, "\n")
			cat("Selected Model Output for ", i, " components\n")
			cat("Ploidy:", i+1, "\n")
			cat("lnL:\t", resList[[i]]$loglik, "\n")
			cat("mixing proportions:\t", resList[[i]]$lambda, "\n")
			cat("parameter estimates\n", resList[[i]]$parameters, "\n")
			cat("Matrix of posterior probabilities\n")
			
			for (l in 1:n)
			{
				for (m in 1:i)
				{
					if (m > 1)
					{
						cat("\t")
					}
					cat(resList[[i]]$posterior[l,m])
					if (m == i)
					{
						cat("\n")
					}
				}
			}
			sink()
			pdf(paste(outPrefix,".component.pdf",sep=""))
            plotComponents(x, lambda=resList[[i]]$lambda, params=resList[[i]]$parameters, model=model)
            dev.off()
		}
	}
	
	if (modelWinner == 0)
	{
		sink(paste(outPrefix,".modelSelection.txt",sep=""))
		cat("Ploidy:", 2, "\n")
		cat("Warning: Ploidy could not be derived from AB spectrum because multiple models were equally bad. Defaulting to 2 but inspection of data recommended!\n")
		sink()
	}
}		