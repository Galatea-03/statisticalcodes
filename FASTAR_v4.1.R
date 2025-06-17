#   FASTAR v.4.1 (based on MixSIR model [Moore and Semmens 2008])

####======  Contact: Aaron Galloway (aaron.galloway@gmail.com)=============================================####
####======  Adapted from MixSIR by Aaron Galloway, Eric Ward, Gordon Holtgrieve, and Mike Brett============####
####======  This and other files are all available for download on my EcologyBox site:=====================####
####======  http://conserver.iugo-cafe.org/user/gway ======================================================####   

#   Last updated 11-Nov-2014 by G. Holtgrieve and A. Galloway prior to PLoSONE revision submission
#   This version will run with any number of producers and consumers, and will accept unbalanced group n
#   Install Jags on machine prior to running. It is not an R "package": get it here: http://sourceforge.net/projects/mcmc-jags/
#   NOTE: variables may be faty acids (FA) or stable isotopes (SI) or possibly both
#   NOTE: FASTAR assumes FA/SI "modification/fractionation" by consumer is accounted for and included in producer file
#   .... therefore, 'producer' signatures are actually values of consumers fed pure diets of those producers in feeding trials
#   .... in 'producers' file, sources are rows, FA means and FA sds are columns (see example files)
#   .... in 'consumer.master' file, rows all consumers. Separate "groups" identified in the Group column (see example)
#   .... in the case of Freshwater Bio (Cladocera) example, each row is it's own group (with one result file and figure per row)

# 0. Clear the workspace 
rm(list = ls())

######################################################################################################################################
#### This is the user input section.  

dir <- "~/Documents/FASTAR_EcoBox_Nov2014"                                        # define your working directory
consumer.file <- "SampleFAs.csv"      # FAs of Samples file name    
producer.file <- "FASources.csv"          # FAs of Sources file name     

modelname = "mixSIR"                          # underlying model presently used is "mixSIR"
mcmc.chainLength <- as.integer(100000)        # post-burn 
mcmc.burn <- as.integer(50000) 
mcmc.thin = 50
mcmc.chains = 3 

#### NOTE: After setting the directory and files above, save changes to the file and 'source' the entire file 
######################################################################################################################################

# 1. Load libraries and all dependencies
library(R2jags)
library(gtools)
library(gdata)
library(ggplot2)

# 2. Set working directory
setwd("D:/R/FASTAR")

# 3. Write model files to working folder
model = cat("

model {
            # prior 
            for(i in 1:num.prey) {
            alpha[i] <- 1;
            }
            p[1:num.prey] ~ ddirch(alpha[]);  # these are weights for means
            # bookkeeping
            for(i in 1:num.prey) {
            p2[i] <- p[i]*p[i]; # these are weights for variances
            }
            
            # for each fa variable and population, calculate the predicted mixtures
            for(fa in 1:num.fa) {
            mix.mu[fa] <- inprod(u[,fa],p[]);
            mix.var[fa] <- inprod(sigma2[,fa],p2[]);
            mix.totalVar[fa] <- mix.var[fa];
            mix.prcsn[fa] <- 1/(mix.totalVar[fa]);
            }
            
            # This section does the likelihood / posterior, N data points
            for(i in 1:N) {
            for(fa in 1:num.fa) {
            X[i,fa] ~ dnorm(mix.mu[fa], mix.prcsn[fa]);
            }
            }
            
            }  
            
            ", file = "mixSIRmodel.txt")
              
# 3. Define the producers and master consumer lists from the files read in at the top:
producers <- read.csv(producer.file, stringsAsFactors=FALSE, header = T)
consumers.master <- read.csv(consumer.file, stringsAsFactors=FALSE, header = T)

# 4. Select 'consumer'(s) that the analysis will be run on from the consumer.master file. 
# A group-compiled figure and csv results file will be generated for each unique group.
consumer.groups <- unique(consumers.master$Group)

# 5. Iteratively run each of the individual consumers within a group, saving each p and then 
# reports the summary of all of the individual p's as one result and plot? 
for (k in consumer.groups){

    # 5.a. Defining data:
    index <- which(consumers.master$Group==k)
    X_group <- as.matrix(consumers.master[index,-c(1,2)])        # drop first 2 columns (individual and group ID info)
    num.fa <- dim(X_group)[2]                                    # number of fa / variables
    num.prey <- dim(producers)[1]                                # of producers 
    Sources <- producers[1:num.prey,1]                           # vector of names for the different producer sources
    N_group <- dim(X_group)[1]                                   # this is the number of consumers / data points 
    meansAndSDs = as.matrix(producers[,-1])                      # fractionation already accounted for in 'producer' file
    u <- meansAndSDs[,1:num.fa]
    sigma2 <- (meansAndSDs[,(num.fa+1):dim(meansAndSDs)[2]])^2
    
    # 5.b. Print summary of the data to the console
    cat("\n", "Number of consumer groups: ", length(consumer.groups), "\n",
        as.character(consumer.groups), "\n",
        "Number of consumers: ", N_group, "\n",
        "Numer of sources: ", num.prey, "\n",
        as.character(levels(Sources)), "\n",
        "Number of tracers: ", num.fa, "\n", "\n", "\n")
    
    # 5.c. Parameters that may be set by the user for MCMC estimation
    mcmc.nSavedDraws <- (mcmc.chainLength-mcmc.burn)/mcmc.thin*mcmc.chains
    runif(1)
    
    # 5.d. Run model. 
    jags.data = list("u", "sigma2", "N", "num.prey", "num.fa",  "X")
    jags.params = c("p", "mix.mu", "mix.var", "mix.prcsn")
    model.loc = ("mixSIRmodel.txt")
    
    # 5.e. Run model for each data point individually and combine into a mixed distribution
    #establish list to store results
    indivP <- lapply(seq_len(num.prey), matrix, data = NA, nrow = mcmc.nSavedDraws, ncol = N_group)
    names(indivP) <- paste("p", 1:num.prey)
    
    # 5.f. Loop through consumers running JAGS for each individually
    for (i in 1:N_group){
      X <- t(as.matrix(X_group[i,]))
      N <- dim(X)[1]
      jags.model = jags(jags.data, inits = NULL, parameters.to.save = jags.params, model.file = model.loc, n.chains = mcmc.chains, 
                        n.burnin = mcmc.burn, n.thin = mcmc.thin, n.iter = mcmc.chainLength, DIC = TRUE)  #runs the Bayesian
      attach.jags(jags.model) 
      
      # Store results from each model run in a list where # of elements is equal to number of sources (num.prey)
      for (ii in seq(1,num.prey, 1))  indivP[[ii]][,i] <- p[,ii]  
      
      detach.jags()  
    }
    
    # 5.g. Plot the individual and combined posterior densities for each source
    # the dimensions of plots are a starting point. User may want to adjust for examples with >5 sources (e.g., cladocera file)
    pdf(file=paste("sourceDensities_", k,".pdf", sep=""), width=2*num.prey, height=0.75*num.prey)  
    op <- par(mfrow=c(1,num.prey), mai =c(0.33, 0.6, 0.25, 0.1))            
    cols <- adjustcolor(rainbow(N_group),alpha.f=0.5)
    
    for (j in 1:num.prey){
        dens <- vector("list", N_group)
        for (i in 1:N_group){
          dens[[i]] <- density(indivP[[j]][,i])
          if(i==1) { y <- dens[[i]]$y} else if(i!=1) {y <- c(y,dens[[i]]$y)}
        }
        plot(dens[[1]], main="", ylim=c(0, max(y)), xlim=c(0,1), xlab=Sources[j], cex=2, type="n")  # establish empty plot
        for (i in 1:N_group) polygon(dens[[i]], col=cols[i], lwd=0.25)  # add polygon of posterior for each consumer
        # Plot the combined density for the group (e.g., contents of results file) as a heavy line
        par(new=T); plot(density(indivP[[j]]), axes=F, xlim=c(0,1), xlab="", font.main=1, main=Sources[j], lwd=1) 
    }
    dev.off()
    par(op)
    
    # 5.h. Store the percentiles into data summary table called Results, and save as a csv file  
    ## This currently saves percentiles after pooling the individuals in a given group
    Percentiles = seq(0, 1, 0.01)
    Results = array(NA, dim = c(num.prey, length(Percentiles)))
    rownames(Results) <- Sources; colnames(Results) <- Percentiles
     for(i in 1:num.prey) Results[i,] = quantile(indivP[[i]], Percentiles)
    
    # 5.i. Send results file to a results folder in the working directory
    write.csv(Results, paste("percentiles_", k,".csv", sep = ""))
}