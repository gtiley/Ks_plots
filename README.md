## Ks_plots
R functions for Ks plot analyses with mixture models

### Updated 13 February 2022

Some bits of code were adapted to create the new functions `ploidy.test()` and `onlyVariances()`. These are used to evaluate the ploidy level of an individual based on distributions of allele balance from genotyping data. The ideas were similar enough to leverage the same mixture model optimization routines.  However, for ploidy, there are clear expectations for the mixing proportions and means. These are fixed on only the variances are optimized.

Any model with a weight greater than 0.5 is returned as the best ploidy level for that individual. If the results are ambiguous, with no model weight > 0.5, a ploidy level of 2 is returned by default and a warning provided in the output.

Only ploidy levels up to 8 are supported, and perhaps it is dubious to go even beyond 6. 


### Functionality prior to 13 February 2022

There are 6 model types available to use for mixture model analyses:
* model 1 = exp + (k-1) normals
* model 2 = exp + (k-1) gammas
* model 3 = exp + (k-1) lognormals
* model 4 = k normals
* model 5 = k gammas
* model 6 = k lognormals

Once the R source code is loaded, perhaps the function of most interest is
mixEM <- function (x, lambda = NULL, alpha = NULL, beta = NULL, k = NULL, model = 1, nstarts = 100, epsilon = 1e-08, maxit = 1000, maxrestarts=20, verb=0)
which optimizes model parameters by expectation maximization.
In this case, all model parameters are unknown except the number of components and x is simply a vector of numeric data. Thus, if I wanted to maximize the likelihood function for model that uses a mixture of 1 exponential distribution and 2 lognormal distributions, I could go
mixEM(x, k=3, model = 3, nstarts = 100, epsilon = 1e-08, maxit = 1000, maxrestarts=20, verb=0)

Usually we want to select the best mixture model from a range of models that use different numbers of components though. There are multiple ways to go about this, but using the difference in the BIC scores as an approximation for Bayes factors is popular. Thus, the model selection for a range of k components can be automated with the function
bic.test.wgd <- function (x, startK = 1, maxK = 2, model = 1, nstarts = 100, outPrefix = NULL)
Thus, if you wanted to pick the best mixture model of an exponential distribution and somewhere between 1 and 4 lognormal distributions, we can go
bic.test.wgd(x, startK = 2, maxK = 5, model = 3, nstarts = 100, outPrefix = myResults)
The output prefix will allow you to automatically print the results out to file, including a plot of your data with the theoretical mixtures of distributions on top.