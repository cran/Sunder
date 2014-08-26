\name{MCMCCV} 
\alias{MCMCCV}
\title{Inference and model selection for analysis of geographical
  genetic variation}
\description{MCMC inference and model selection by cross-validation for the analysis of geographical
  genetic variation}
\usage{MCMCCV(gen, D_G, D_E,
                nit, thinning, theta.max, theta.init,
                run, ud, n.validation.set,print.pct)
}



\arguments{
  \item{gen}{
An array with dimensions (n,l,a) n: number of geographical locations, l:
  number of loci, a: max number of alleles
}
  \item{D_G}{A matrix of geograpical distances}
  \item{D_E}{A matrix of environmental distances}
  \item{nit}{Number of iterations}
  \item{thinning}{Thinning of MCMC iterations}
  \item{theta.max}{Upper bounds for the vector of parameters in theta. Note that
  in theta, the parameters are assumed to be in this order:
  (alpha,beta_G, beta_E, gamma, delta)}
  \item{theta.init}{Initial state of theta}
  \item{run}{A vector of booleans of length 3 stating which sub-model is
  investigated  among G+E,G,E (in this order). For example, the default value
  \code{run=c(TRUE,FALSE,FALSE)} means that only one MCMC run will be performed to
  estimate paremeters in the G+E model while
  e.g. \code{run=c(TRUE,FALSE,TRUE)} means that MCMC runs will be
  performed for models G+E and E. If \code{n.validation.set >0}
  likelihoods under these two models will be returned for model selection.}
  \item{ud}{A vector  of booleans  of length 5 stating which entries in
    theta=(alpha,beta_E,beta_G,gamma,delta) 
    will be updated in the MCMC iterations. By default all parameters in
  theta will be updated. If one entry is not updated, the value used
  along the MCMC simulation for this parameter is the initial value.
}
  \item{n.validation.set}{The number of pairs (sites x locus) used as
    validation set}
  \item{print.pct}{A boolean stating whether Fortran prints percentage
  of computation achieved along each MCMC run. }
}
\value{A list with a component named mod.lik containing likelihoods on
  the   validation set for the various models compared.}
\author{Filippo Botta, Gilles Guillot}
\keyword{MCMC inference}
\keyword{Model selection}

\examples{
\dontrun{
data(toydata,package='Sunder')
     
     #### Computing options
     nit <- 10^2
     run  <- c(1,1,1)
     thinning  <- 1 # max(nit/10^3,1);
     ud   <- c(0,1,1,0,0) 
     theta.init <- c(1,2,1,1,0.01)
     n.validation.set  <- dim(gen)[1]*dim(gen)[2]/10 
     theta.max  <- c(10,10*max(D_G),10*max(D_E),1,0.01)
     
     plot  <- TRUE
     trace <- FALSE
     
     #### Call Sunder ####
     output <- MCMCCV(gen,D_G,D_E,
                                nit,thinning,
                                theta.max,
                                theta.init,
                                run,ud,n.validation.set)
     
     mod.lik <- output$mod.lik
     tvt <- output$theta
     
     
     ## plotting outputs
     upd=matrix(nrow=sum(run), ncol=length(theta.init), data=1)
     upd[2,3]=upd[3,2]=0
     
     plot(as.vector(D_G),as.vector(cor(t(gen[,,1]))),
          bg=colorRampPalette(c("blue", "red"))(dim(D_E)[1]^2)[order(order(as.vector(D_E)))],
          pch=21,
          xlab='Geographic distance',
          ylab='Empirical covariance genotypes')
     
     
     kol=c(4,2,3) 
     xseq=seq(thinning,nit,thinning)
     ylab=c(expression(paste(alpha)),
            expression(paste(beta[D])),
            expression(paste(beta[E])),
            expression(paste(gamma)),
            expression(paste(delta)))
     
     par(mfrow=c(sum(run),length(theta.init)))
     for (j in 1:sum(run))
     {
       for(k in 1:length(theta.init))
       {
         if (sum(upd[,k]==1)>0)
         {
           if(upd[j, k]==1)
           {
             if(exists("theta"))
               ylim=c(min(tvt[k,,j],theta[k]),max(tvt[k,,j],theta[k])) else
                 ylim=c(min(tvt[k,,j]),max(tvt[k,,j]))
             plot(0, type="n",xlab="",ylab="", xlim=c(0,nit), ylim=ylim)
             lines(xseq,tvt[k,,j],col=kol[j],xlab="",ylab="")
             if(exists("theta")) abline(h=theta[k],lty=2)
             title(xlab="iterations")        
             mtext(ylab[k], side=2, line=2.3,las=1)} else plot.new()
         }
       }
     }

print(mod.lik)
print(paste('The model achieving the highest likelihood on the validation set is:',
            names(mod.lik)[order(mod.lik,decreasing=TRUE)[1]]))
theta.GE <- apply(output$theta[,,1],1,mean)
print('Posterior mean theta under model G+E:')
print(theta.GE)

theta.G <- apply(output$theta[,,2],1,mean)
theta.G[3] <- NA
print('Posterior mean theta under model G:')
print(theta.G)

theta.E <- apply(output$theta[,,3],1,mean)
theta.E[2] <- NA
print('Posterior mean theta under model E:')
print(theta.E)

}
}


