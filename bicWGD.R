bic.test.wgd <- function (x, startK = 1, maxK = 2, model = 1, nstarts = 100, outPrefix = NULL)
{
	n <- length(x)
	thisK <- startK
	check = 0
	if (is.null(outPrefix))
	{
		outPrefix <- "wgdMix"
	}
	while(thisK < maxK && check == 0)
	{
		if (thisK == startK)
		{
			obs.res.null <- mixEM(x, k=thisK, model=model,nstarts=nstarts)
			obs.res.alt <- mixEM(x, k=(thisK + 1), model=model,nstarts=nstarts)
		}
		if (thisK > startK)
		{
			obs.res.null <- obs.res.alt
			obs.res.alt <- mixEM(x, k=(thisK + 1), model=model,nstarts=nstarts)
		}
		nullBIC <- -2 * obs.res.null$loglik + (thisK * 2 - 1) * log(n)
		altBIC <- -2 * obs.res.alt$loglik + ((thisK+1) * 2 - 1) * log(n)
		
		deltaBIC <- (nullBIC - altBIC)
		cat ("deltaBIC:\t",deltaBIC,"\nnull_lnL:\t",obs.res.null$loglik,"\nalt_lnL:\t",obs.res.alt$loglik,"\n")
		if (deltaBIC > 3.2)
		{
			if (thisK + 1 == maxK)
			{
				sink(paste(outPrefix,".modelSelection.txt",sep=""))
				cat("Model Choice ", model, "\n")
				cat("Selected Model Output for ", maxK, " components\n")
				cat("lnL:\t", obs.res.alt$loglik, "\n")
				cat("mixing proportions:\t", obs.res.alt$lambda, "\n")
				cat("parameter estimates\n", obs.res.alt$parameters, "\n")
				cat("Matrix of posterior probabilities\n")
			
				for (l in 1:n)
				{
					for (m in 1:maxK)
					{
						if (m > 1)
						{
							cat("\t")
						}
						cat(obs.res.alt$posterior[l,m])
						if (m == maxK)
						{
							cat("\n")
						}
					}
				}
				sink()
				pdf(paste(outPrefix,".component.pdf",sep=""))
                                plotComponents(x, lambda=obs.res.alt$lambda, params=obs.res.alt$parameters, model=model)
                                dev.off()
				check = 1
			}
			else
			{
				thisK = thisK + 1
			}
		}
		if (deltaBIC <= 3.2)
		{
			sink(paste(outPrefix,".modelSelection.txt", sep=""))
			cat("Model Choice ", model, "\n")
			cat("Selected Model Output for ", thisK, " components\n")
			cat("lnL:\t", obs.res.null$loglik, "\n")
			cat("mixing proportions:\t", obs.res.null$lambda, "\n")
			cat("parameter estimates\n", obs.res.null$parameters, "\n")
			cat("Matrix of posterior probabilities\n")
			
			for (l in 1:n)
			{
				for (m in 1:thisK)
				{
					if (m > 1)
					{
						cat("\t")
					}
					cat(obs.res.null$posterior[l,m])
					if (m == thisK)
					{
						cat("\n")
					}
				}
			}
			sink()
			pdf(paste(outPrefix,".component.pdf",sep=""))
                        plotComponents(x, lambda=obs.res.null$lambda, params=obs.res.null$parameters, model=model)
                        dev.off()
			check = 1
		} 
	}
}		