time.stamp <- "Time-stamp: <aj:/Users/andy/Desktop/ssim-ghg/EnKF/base/ensemble_size.r - 07 Jul 2025 (Mon) 12:48:58 MDT>"
cat(sprintf("[Script info] %s\n",time.stamp))

# This code applies the EnKF measurement update with a varying number
# of ensemble members.

# Add local library for R packages
.libPaths(new=c(sprintf("%s/shared/lib/R-4.3/x86_64-pc-linux-gnu",Sys.getenv("HOME")),.libPaths())) 

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

pctage <- 100*n.selected/nobs
cat(sprintf("%d obs selected out of %d available (%.1f%%)\n",
            n.selected,nobs,pctage))

posterior.dofs <- FALSE

if(posterior.dofs) {
  ndofs.kf <- ndofs.patil(state.kf$Sx.post)
} else {
  ndofs.kf <- nmons*nreg
}

# assimilate
these.obs <- obs[lx.selected]
these.Szd <- Szd.assumed[lx.selected]



# Perform the analytical KF inversion as a baseline
state.kf <- list()
state.kf$x.prior <- matrix(1,nrow=nparms,ncol=1)
state.kf$x.post <- matrix(NA,nrow=nparms,ncol=1)

# We need a state vector and an ensemble with the total number of
# parameters, 22x24, to multiply the Jacobian H. 
t0 <- proc.time()[3]
y.prior <- simulate_observed(H=H[lx.selected,],
                             x=state.kf$x.prior)
cat(sprintf('y.prior computed: %.1fs\n',proc.time()[3]-t0))

# diag() can form a diagonal matrix from a vector, or extract the
# diagonal of a matrix. Behavior depends on its argument.
Sx.prior <- diag(1.3,nparms) 
state.kf$x.prior[,1] <- 1
state.kf$Sx.prior <- Sx.prior

these.obs <- obs[lx.selected]
these.Szd <- Szd.assumed[lx.selected]

t0 <- proc.time()[3]
post <- kf_meas_update(x=state.kf$x.prior,
                       Sx=state.kf$Sx.prior,
                       H=H[lx.selected,],
                       z=these.obs,
                       Sz=diag(these.Szd))

cat(sprintf('KF update: %.1fs\n',proc.time()[3]-t0))

state.kf$x.post <- post$x
state.kf$Sx.post <- post$Sx
y.kf.post <- simulate_observed(x=state.kf$x.post,H=H[lx.selected,])
y.kf.prior <- y.prior


chi2.state.kf <- (1/ndofs.kf) * t(state.kf$x.post - truth_condition) %*% solve(state.kf$Sx.post) %*% (state.kf$x.post - truth_condition)
chi2.prior.kf <- (1/ndofs.kf) * t(state.kf$x.prior - truth_condition) %*% solve(state.kf$Sx.prior) %*% (state.kf$x.prior - truth_condition)
chi2.innov.kf <- (1/ndofs.kf) * t(state.kf$x.post - state.kf$x.prior) %*% solve(state.kf$Sx.prior) %*% (state.kf$x.post - state.kf$x.prior)

chi2.obs.kf <- (1/n.selected) * t(obs[lx.selected] - y.kf.prior) %*% solve(H[lx.selected,] %*% state.kf$Sx.prior %*% t(H[lx.selected,]) + diag(Szd.assumed[lx.selected])) %*% (obs[lx.selected] - y.kf.prior)

cat(sprintf(" [KF] chi2 means: state %.2f, prior %.2f, innov %.2f, obs %.2f on %d (%d) DOFs, RMSE %.2f\n",
            chi2.state.kf,chi2.prior.kf,chi2.innov.kf,chi2.obs.kf,ndofs.kf,ndofs.patil(state.kf$Sx.post),compute.rmse(state.kf$x.post - truth_condition)))


for (nmemb in c(seq(1000,10000,by=1000),20000,100000)) {
    
    state.ens <- list()
    state.ens$x.prior <- matrix(1,nrow=nparms,ncol=1)
    state.ens$dx.prior <- matrix(NA,nrow=nmemb,ncol=nparms)
    state.ens$x.post <- matrix(NA,nrow=nparms,ncol=1)
    state.ens$dx.post <- matrix(NA,nrow=nmemb,ncol=nparms)

    state.ens$x.prior[,1] <- 1
    state.ens$dx.prior <- generate_ensemble(Sx=Sx.prior,nmemb=nmemb) # prior deviations

    

#    t0 <- proc.time()[3]
#    y.prior <- simulate_observed(H=H[lx.selected,],
#                                 x=state.ens$x.prior)
    y.ens.prior <- y.prior
    dy.prior <- t(simulate_observed(H=H[lx.selected,],
                                    x=t(state.ens$dx.prior)))
#    cat(sprintf('y,dy.prior for %d members computed: %.1fs\n',nmemb,proc.time()[3]-t0))

    
#    t0 <- proc.time()[3]
    post <- enkf_meas_update_loc(x=state.ens$x.prior,
                                 dx=state.ens$dx.prior,
                                 obs=these.obs,
                                 Szd=these.Szd,
                                 y=y.prior,dy=dy.prior,
                                 localization_mask=NULL)
    
#    cat(sprintf('EnKF computed: %.1fs\n',proc.time()[3]-t0))

    state.ens$x.post <- post$x
    state.ens$dx.post <- post$dx

    # Compute sample covariance matrices (528x528)
    state.ens$Sx.post <- matrix(0,nrow=nmons*nreg,ncol=nmons*nreg)
    state.ens$Sx.prior <- matrix(0,nrow=nmons*nreg,ncol=nmons*nreg)
    state.ens$Sx.post <- cov(state.ens$dx.post)  
    state.ens$Sx.prior <- cov(state.ens$dx.prior)



    if(posterior.dofs) {
        ndofs.ens <- ndofs.patil(state.ens$Sx)
    } else {
        ndofs.ens <- nmons*nreg
    }

    # Except we can't have any more DOFs than ensemble members
    ndofs.ens <- min(ndofs.ens,nmemb)

    chi2.state.ens <- (1/ndofs.ens) * t(state.ens$x.post - truth_condition) %*% solve(state.ens$Sx.post) %*% (state.ens$x.post - truth_condition)
    chi2.prior.ens <- (1/ndofs.ens) * t(state.ens$x.prior - truth_condition) %*% solve(state.ens$Sx.prior) %*% (state.ens$x.prior - truth_condition)
    chi2.innov.ens <- (1/ndofs.ens) * t(state.ens$x.post - state.ens$x.prior) %*% solve(state.ens$Sx.prior) %*% (state.ens$x.post - state.ens$x.prior)
    
    chi2.obs.ens <- (1/n.selected) * t(obs[lx.selected] - y.ens.prior) %*% solve(H[lx.selected,] %*% state.ens$Sx.prior %*% t(H[lx.selected,]) + diag(Szd.assumed[lx.selected])) %*% (obs[lx.selected] - y.ens.prior)
    
    cat(sprintf("[Ens] chi2 means: state %.2f, prior %.2f, innov %.2f, obs %.2f on %d (%d) DOFs, RMSE %.2f (%d members)\n",
                chi2.state.ens,chi2.prior.ens,chi2.innov.ens,chi2.obs.ens,ndofs.ens,ndofs.patil(state.ens$Sx.post),compute.rmse(state.ens$x.post - truth_condition),nmemb))
    
    if(FALSE) {
        plot.x.timeseries(ests=list(Truth=list(x=truth_condition),
                                    Ensemble=list(x=state.ens$x.post,Sx=state.ens$Sx),
                                    Analytical=list(x=state.kf$x.post,Sx=state.kf$Sx)))
        #                  pdf.name="common.x.pdf")

        plot.flux.timeseries(ests=list(Truth=list(x=truth_condition),
                                       Ensemble=list(x=state.ens$x.post,Sx=state.ens$Sx),
                                       Analytical=list(x=state.kf$x.post,Sx=state.kf$Sx)))
        #                  pdf.name="common.flux.pdf")


        plot.is.timeseries(xs=list(Truth=truth_condition,
                                   Ensemble=state.ens$x.post,
                                   Analytical=state.kf$x.post,
                                   Prior=x.prior),
                           dataset_names=c("co2_mlo_surface-insitu_1_allvalid",
                                           "co2_brw_surface-insitu_1_allvalid",
                                           "co2_smo_surface-insitu_1_allvalid",
                                           "co2_spo_surface-insitu_1_allvalid",
                                           "co2_lef_tower-insitu_1_allvalid-396magl"),
                           H=H.orig)
        #                   pdf.name='common.obs.pdf')
    }
}

