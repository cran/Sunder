\name{SunderInference}
\alias{SunderInference}
\title{Inference and model selection for analysis of geographical
  genetic variation}
\description{MCMC inference and model selection by cross-validation for the analysis of geographical
  genetic variation}
\usage{SunderInference(gen, D_G, D_E, nit, thinning, theta.max, theta.init, run, ud, set)}


\arguments{
  \item{gen}{
An array with dimensions (n,l,a) n: number of geographical locations, l:
  number of loci, a: max number of alleles
}
  \item{D_G}{A matrix of geograpical distances}
  \item{D_E}{A matrix of environmental distances}
  \item{nit}{Number of iterations}
  \item{thinning}{Thinning of MCMC iterations}
  \item{theta.max}{Upper bounds for the parameters in theta}
  \item{theta.init}{Initial state of theta}
  \item{run}{A vector of 0/1 of length 3 stating which sub-model is
  investigated  among G+E,G,E (in this order)}
  \item{ud}{A vector of length 5 stating which entries in
    theta=(alpha,beta_E,beta_G,gamma,delta) 
    will be updated in the MCMC iterations. 
}
  \item{set}{The number of pairs (sites x locus) used as validation set}
}
\value{A list with a component named mod.lik containing likelihoods on
  the   validation set for the various models compared.}
\author{Filippo Botta}
\keyword{MCMC inference}
\keyword{Model selection}

\examples{
\dontrun{
data(toydata,package='Sunder')

#### Computing options
nit <- 10^2
run  <- c(1,1,1)
thinning  <- max(nit/10^3,1);
ud   <- c(0,1,1,0,0) 
theta.init <- c(1,2,1,1,0.01)
set  <- dim(gen)[1]*dim(gen)[2]/10 
theta.max  <- c(10,10*max(D_G),10*max(D_E),1,0.01)

plot  <- TRUE
trace <- FALSE

#### Call Sunder ####
output <- SunderInference(gen,D_G,D_E,
                           nit,thinning,
                           theta.max,
                           theta.init,
                           run,ud,set)

mod.lik <- output$mod.lik
tvt <- output$tvt


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
}
}

