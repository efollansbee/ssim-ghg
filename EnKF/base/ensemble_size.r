time.stamp <- "Time-stamp: <aj:/Users/andy/Desktop/ssim-ghg/EnKF/base/ensemble_size.r - 08 Jul 2025 (Tue) 14:53:53 MDT>"
cat(sprintf("[Script info] %s\n",time.stamp))

# This code applies the EnKF measurement update with a varying number
# of ensemble members, plus an anakytical KF result for
# comparison. Results are stored in a data frame "experiments". The
# rows are experimental conditions, the (named) columns are input or
# output quantities.

# Add local library for R packages
.libPaths(new=c(sprintf("%s/shared/lib/R-4.3/x86_64-pc-linux-gnu",Sys.getenv("HOME")),.libPaths())) 

make.pdfs <- TRUE # not for notebooks, but for CLI

source("../tools/enkf.r")
source("../tools/progress.bar.r")
source("../tools/find.indir.r")
source("../tools/time.r")
source("../tools/load.ncdf4.r")
options(warn=2) # error out on warnings (probably they are mistakes)
indir <- find.indir()

# Since these load() statements can take some time, we can supply
# them pre-loaded from memory or load them here.
if(!exists("H.orig")) {
    t0 <- proc.time()[3]
    cat("Loading Jacobians...")
    load(file.path(indir,"jacobians/trunc_full_jacob_030624_with_dimnames_sib4_4x5_mask.rda"))
    H <- jacob*(12/44) # Andrew reports units conversion needed
    rm(jacob)
    # We'll be subsetting H later; preserve the original matrix for
    # later use.
    H.orig <- H
    cat(sprintf('%.1fs\n',proc.time()[3]-t0))
}

if(!exists("obs_catalog")) {
    t0 <- proc.time()[3]
    cat("Loading obs_catalog...")
    load(file.path(indir,"obs/obs_catalog_042424_unit_pulse_hour_timestamp_witherrors_withdates.rda"))
    # change time zone to "UTC" for obs_catalog times (is by default
    # "GMT")
    attributes(obs_catalog$DATE)$tzone <- "UTC"
    cat(sprintf('%.1fs\n',proc.time()[3]-t0))
}

nreg <- 22
nmons <- 24
nparms <- nreg*nmons

# probability (0-1) that obs will be randomly selected for
# assimilation
obs.prob <- 0.005


load(file.path(indir,"misc/truth_array.rda"))

# Andrew uses a tm() function to limit the range of values in an
# array. This function does the same but is very explicit. Note that
# the which() function temporarily casts x from a multidimensional
# array to a vector, and the x[lx] assignment does the same. The
# dimensions attribute of x is unaffacted, so the array size is
# unaffected. 
trim <- function(x,limits) {
    if(length(limits)!=2) {
        stop("Expecting 2-vector for limits argument")
    }
    lx <- which(x < limits[1])
    if(length(lx)>0) {
        x[lx] <- limits[1]
    }
    lx <- which(x > limits[2])
    if(length(lx)>0) {
        x[lx] <- limits[2]
    }
    return(x)
}

truth_condition <- -1 * trim(truth_array[1:24,2:23,1,1],c(-2,2))
dim(truth_condition) <- c(nreg*nmons,1) 
# Finally, add back the 1.0 that Andrew removed
truth_condition <- truth_condition + 1.0

# generate obs
nobs <- dim(H)[1]
# "d" suffix means its the diagonal (of a diagonal matrix)
Szd.actual <- rep(0.3^2,nobs) # variance in ppm^2
#Szd.assumed <- rep((0.1)^2,nobs) # variance in ppm^2
Szd.assumed <- Szd.actual
# Note that supplying a Szd argument to the simulate_observed function
# will result in perturbations being added to the observations.  This
# would be the place to add biases to obs.
obs <- simulate_observed(H=H,
                         x=matrix(truth_condition,nrow=nreg*nmons),
                         Szd=Szd.actual)
dim(obs) <- c(nobs,1)

time.edges <- seq(ISOdatetime(2014,9,1,0,0,0,tz="UTC"),
                  ISOdatetime(2016,12,1,0,0,0,tz="UTC"),
                  by="1 month")

# Alternative to this probability-based approach, can
# sample(1:nobs,n.selected) to get exactly n.selected obs.
lx.selected <- which(runif(n=nobs)<=obs.prob)
n.selected <- length(lx.selected)

# Set up a data frame with all the experimental conditions we're going
# to run.
experiments <- NULL
iexp <- 0
for (nmembers in unique(sort(c(0,529,seq(1000,10000,by=1000),seq(1500,4500,by=1000),n.selected,20000,50000)))) {
    iexp <- iexp + 1
    experiments <- rbind(experiments,
                         data.frame(nmembers=nmembers,method="ensemble",
                                    chi2.obs=NA,chi2.state=NA,
                                    RMSE=NA,ndofs.patil=NA,
                                    elapsed.time=NA))
    if(nmembers==0) {
        experiments$method[iexp] <- "KF"
    }
}

nexperiments <- dim(experiments)[1]

pctage <- 100*n.selected/nobs
cat(sprintf("%d obs selected out of %d available (%.1f%%)\n",
            n.selected,nobs,pctage))

posterior.dofs <- FALSE

these.obs <- obs[lx.selected]
these.Szd <- Szd.assumed[lx.selected]

# diag() can form a diagonal matrix from a vector, or extract the
# diagonal of a matrix. Behavior depends on its argument.
x.prior <- matrix(1,nrow=nparms,ncol=1)
#Sx.prior <- diag(5,nparms) 
Sx.prior <- diag(1.4,nparms) 


# We need a state vector and an ensemble with the total number of
# parameters, 22x24, to multiply the Jacobian H. 
#t0 <- proc.time()[3]
y.prior <- simulate_observed(H=H[lx.selected,],
                             x=x.prior)
#cat(sprintf('y.prior computed: %.1fs\n',proc.time()[3]-t0))



these.obs <- obs[lx.selected]
these.Szd <- Szd.assumed[lx.selected]


for (iexp in 1:nexperiments) {

    nmembers <- experiments$nmembers[iexp]
    method <- experiments$method[iexp]
    
    state <- list()
    state$x.prior <- x.prior
    state$Sx.prior <- Sx.prior
    # Next lines are just allocating space
    state$x.post <- matrix(NA,nrow=nparms,ncol=1)
    state$Sx.post <- matrix(NA,nrow=nparms,ncol=nparms)

    
    t0 <- proc.time()[3]
    if(nmembers==0) { # special value to signify KF solution

        post <- kf_meas_update(x=state$x.prior,
                               Sx=state$Sx.prior,
                               H=H[lx.selected,],
                               z=these.obs,
                               Sz=diag(these.Szd))
        state$x.post <- post$x
        state$Sx.post <- post$Sx

    } else {
        
        state$dx.prior <- generate_ensemble(Sx=state$Sx.prior,nmemb=nmembers) # prior deviations
        state$dx.post <- matrix(NA,nrow=nmembers,ncol=nparms)

        dy.prior <- t(simulate_observed(H=H[lx.selected,],
                                        x=t(state$dx.prior)))
        
        post <- enkf_meas_update_loc(x=state$x.prior,
                                     dx=state$dx.prior,
                                     obs=these.obs,
                                     Szd=these.Szd,
                                     y=y.prior,dy=dy.prior,
                                     localization_mask=NULL)

        state$x.post <- post$x
        state$dx.post <- post$dx
        # Compute an estimate of Sx.post from posterior ensemble deviations
        state$Sx.post <- cov(state$dx.post)  
        
    }
    dt <- proc.time()[3]-t0

    ndofs <- ndofs.patil(state$Sx.post)
    # Except we can't have any more DOFs than ensemble members
    if(nmembers>0) { # special value to signify KF solution
        ndofs <- min(ndofs,nmembers)
    }

    # Notice the 1/nparms below. This is assuming that the no. DOFs is
    # nparms.
    chi2.state <- (1/nparms) * t(state$x.post - truth_condition) %*% solve(state$Sx.post) %*% (state$x.post - truth_condition)
    chi2.obs <- (1/n.selected) * t(obs[lx.selected] - y.prior) %*% solve(H[lx.selected,] %*% state$Sx.prior %*% t(H[lx.selected,]) + diag(Szd.assumed[lx.selected])) %*% (obs[lx.selected] - y.prior)

    rmse <- compute.rmse(state$x.post - truth_condition)
    cat(sprintf("[%s] chi2 means: state %.2f, obs %.2f on %d (%d) DOFs, RMSE %.2f (%d members) in %.1fs\n",
                method,chi2.state,chi2.obs,nparms,ndofs,rmse,nmembers,dt))

    experiments$chi2.state[iexp] <- chi2.state
    experiments$chi2.obs[iexp] <- chi2.obs
    experiments$chi2.obs[iexp] <- chi2.obs
    experiments$RMSE[iexp] <- rmse
    experiments$ndofs.patil[iexp] <- ndofs
    experiments$elapsed.time[iexp] <- dt
}

if(make.pdfs) {
    pdf("ensemble_size.pdf",width=10,height=3)
    layout(matrix(1:3,nrow=1))
    par(las=1,bg='white')
}

plot(experiments$nmembers[3:nexperiments],experiments$chi2.state[3:nexperiments],log='xy',
     main=expression(paste(text=chi[state]^2)),
     las=1,ylab='',xlab='ensemble size',bg='skyblue',pch=21)
abline(h=experiments$chi2.state[1],col='red')
abline(h=1.0)
legend(x='topleft',legend=c('ensemble','KF','target'),
       lty=c(NA,1,1),pch=c(21,NA,NA),pt.bg='skyblue',
       bty='n',
       col=c('black','red','black'))
       
plot(experiments$nmembers[3:nexperiments],experiments$RMSE[3:nexperiments],log='xy',main="RMSE",las=1,ylab='',xlab='ensemble size',bg='skyblue',pch=21)
abline(h=experiments$RMSE[1],col='red')
legend(x='topleft',legend=c('ensemble','KF'),
       lty=c(NA,1),pch=c(21,NA),pt.bg='skyblue',bty='n',
       col=c('black','red'))
       

plot(experiments$nmembers[3:nexperiments],experiments$chi2.obs[3:nexperiments],log='xy',
     main=expression(paste(text=chi[forecast]^2)),
     las=1,ylab='',xlab='ensemble size',bg='skyblue',pch=21)
abline(h=experiments$chi2.obs[1],col='red')
abline(h=1.0)
legend(x='topleft',legend=c('ensemble','KF','target'),
       lty=c(NA,1,1),pch=c(21,NA,NA),pt.bg='skyblue',bty='n',
       col=c('black','red','black'))
       

if(make.pdfs) {
    dev.close()
}
