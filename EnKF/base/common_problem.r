time.stamp <- "Time-stamp: <aj:/Users/andy/Desktop/ssim-ghg/EnKF/base/common_problem.r - 09 Jul 2025 (Wed) 16:48:38 MDT>"
cat(sprintf("[Script info] %s\n",time.stamp))

# This code applies the EnKF measurement update to a truth condition
# generated from scaling factors derived from OCO-2 v10 MIP models,
# without any time propagation. This means that the state vector is
# exactly that of the toy problem, but in an ensemble context.
#
# state.ens contains all info for the ensemble solution
#
# state.kf contains a the state for a control case using analytical KF

# R quick syntax
#  <- is historically the assignment operator (can also use = nowadays)
#  t() is the matrix transpose operator
#  %*% is matrix multiplication
#  solve() computes the matrix inverse

# Add local library for R packages
.libPaths(new=c(sprintf("%s/shared/lib/R-4.3/x86_64-pc-linux-gnu",Sys.getenv("HOME")),.libPaths())) 

make.pdfs <- TRUE # not for notebooks, but for CLI

source("../tools/enkf.r")
source("../tools/progress.bar.r")
source("../tools/find.indir.r")
source("../tools/time.r")
source("../tools/load.ncdf4.r")
source("../tools/normality.test.r")
options(warn=2) # error out on warnings (probably they are mistakes)
indir <- find.indir()

# Since these load() statements can take some time, we use
# them pre-loaded from memory or if needed, load them here.
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
nmemb <- 5000
# probability (0-1) that obs will be randomly selected for
# assimilation
obs.prob <- 0.005


load(file.path(indir,"misc/truth_array.rda"))

# Andrew uses a tm() function to limit the range of values in an
# array. This function does the same but is more explicit. Note that
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
# "d" suffix means its the vector of the diagonal (of a diagonal
# matrix)
Szd.actual <- rep(0.3^2,nobs) # variance in ppm^2
Szd.assumed <- Szd.actual
#Szd.assumed <- rep((0.1)^2,nobs) # variance in ppm^2

# Note that supplying a Szd argument to the simulate_observed function
# will result in (zero-mean) perturbations being added to the
# observations.  This would be the place to add biases to obs.
obs <- simulate_observed(H=H,
                         x=matrix(truth_condition,nrow=nreg*nmons),
                         Szd=Szd.actual)
dim(obs) <- c(nobs,1)

# Note on the potentially-confusing use of the term prior. In
# the (En)KF, you have a prior state at time k which is conditional on
# all the obs already assimilated, from times 1 to k-1. This use of
# the term prior would perhaps better be called "unoptimized"; we use
# it here to represent the simulated values you'd get from the
# unscaled fluxes from our prior models (SiB4 and Landschutzer).

state.ens <- list()
state.ens$x.prior <- matrix(NA,nrow=nparms,ncol=1)
state.ens$dx.prior <- matrix(NA,nrow=nmemb,ncol=nparms)
state.ens$x.post <- matrix(NA,nrow=nparms,ncol=1)
state.ens$dx.post <- matrix(NA,nrow=nmemb,ncol=nparms)

state.kf <- list()
state.kf$x.prior <- matrix(NA,nrow=nparms,ncol=1)
state.kf$Sx.prior <- matrix(NA,nrow=nmemb,ncol=nparms)
state.kf$x.post <- matrix(NA,nrow=nparms,ncol=1)
state.kf$Sx.post <- matrix(NA,nrow=nmemb,ncol=nparms)

# diag() can form a diagonal matrix from a vector, or extract the
# diagonal of a matrix. Behavior depends on its argument.
Sx.prior <- diag(1.4,nparms)
state.kf$x.prior[,1] <- 1
state.ens$x.prior[,1] <- 1
state.kf$Sx.prior <- Sx.prior
state.ens$dx.prior <- generate_ensemble(Sx=Sx.prior,nmemb=nmemb) # prior deviations

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

if(n.selected > 0) {

    # assimilate
    these.obs <- obs[lx.selected]
    these.Szd <- Szd.assumed[lx.selected]

    # We need a state vector and an ensemble with the total number of
    # parameters, 22x24, to multiply the Jacobian H. 
    y.prior <- simulate_observed(H=H[lx.selected,],
                                 x=state.ens$x.prior)

    dy.prior <- t(simulate_observed(H=H[lx.selected,],
                                    x=t(state.ens$dx.prior)))


    # R stores matrices in column-major format. That means if a matrix
    # "foo" is 2x3 (nrow x ncolumn) and looks like:
    #   > foo
    #          [,1]   [,2]   [,3]
    #   [1,] "r1c1" "r1c2" "r1c3"
    #   [2,] "r2c1" "r2c2" "r2c3"
    #   > dim(foo)
    #   [1] 2 3
    
    # Then if we assign dims 6x1 to it, we will get:
    #   > dim(foo) <- c(6,1)
    #   > foo
    #          [,1]
    #   [1,] "r1c1"
    #   [2,] "r2c1"
    #   [3,] "r1c2"
    #   [4,] "r2c2"
    #   [5,] "r1c3"
    #   [6,] "r2c3"
    #
    # The r1c1, r2c1, r1c2, ... ordering is how column-major data are
    # stored in memory: the rows change quickest, the columns slowest.
    
    # The original matrix can be retrieved by just assigning 2x3
    # dimensions back to foo.  We will change dimensions on matrices
    # knowing that this storage order is preserved:
    #
    #   > dim(foo) <- c(2,3)
    #   > foo
    #          [,1]   [,2]   [,3]
    #   [1,] "r1c1" "r1c2" "r1c3"
    #   [2,] "r2c1" "r2c2" "r2c3"
    
    post <- enkf_meas_update_loc(x=state.ens$x.prior,
                                 dx=state.ens$dx.prior,
                                 obs=these.obs,
                                 Szd=these.Szd,
                                 y=y.prior,dy=dy.prior,
                                 localization_mask=NULL)
    
    state.ens$x.post <- post$x
    state.ens$dx.post <- post$dx

    post <- kf_meas_update(x=state.kf$x.prior,
                           Sx=state.kf$Sx.prior,
                           H=H[lx.selected,],
                           z=these.obs,
                           Sz=diag(these.Szd))
    
    state.kf$x.post <- post$x
    state.kf$Sx.post <- post$Sx
    
} # if n.selected > 0

# Compute sample covariance matrices (528x528)
state.ens$Sx.post <- matrix(0,nrow=nmons*nreg,ncol=nmons*nreg)
state.ens$Sx.prior <- matrix(0,nrow=nmons*nreg,ncol=nmons*nreg)
state.ens$Sx.post <- cov(state.ens$dx.post)  
state.ens$Sx.prior <- cov(state.ens$dx.prior)


posterior.dofs <- FALSE

if(posterior.dofs) {
  ndofs.ens <- ndofs.patil(state.ens$Sx)
  ndofs.kf <- ndofs.patil(state.kf$Sx.post)
} else {
  ndofs.ens <- nmons*nreg
  ndofs.kf <- nmons*nreg
}

# Except we can't have any more DOFs than ensemble members
ndofs.ens <- min(ndofs.ens,nmemb)

# This is where we finally use lx.selected, to make posterior
# simulated values for the measurements we assimilated.
#obs.enkf.post <- simulate_observed(x=state.ens$x.post.finals,H=H[lx.selected.all,])
obs.ens.prior <- y.prior
obs.ens.post <- simulate_observed(x=state.ens$x.post,H=H[lx.selected,])
obs.kf.prior <- y.prior
obs.kf.post <- simulate_observed(x=state.kf$x.post,H=H[lx.selected,])

chi2.state.ens <- (1/ndofs.ens) * t(state.ens$x.post - truth_condition) %*% solve(state.ens$Sx.post) %*% (state.ens$x.post - truth_condition)
chi2.obs.ens <- (1/n.selected) * t(obs[lx.selected] - obs.ens.prior) %*% solve(H[lx.selected,] %*% state.ens$Sx.prior %*% t(H[lx.selected,]) + diag(Szd.assumed[lx.selected])) %*% (obs[lx.selected] - obs.ens.prior)

cat(sprintf("[Ens] chi2 means: state %.2f, obs %.2f on %d (%d) DOFs, RMSE %.2f (%d members)\n",
            chi2.state.ens,chi2.obs.ens,ndofs.ens,ndofs.patil(state.ens$Sx.post),compute.rmse(state.ens$x.post - truth_condition),nmemb))

chi2.state.kf <- (1/ndofs.kf) * t(state.kf$x.post - truth_condition) %*% solve(state.kf$Sx.post) %*% (state.kf$x.post - truth_condition)
resid.kf <- obs[lx.selected] - obs.kf.prior
chi2.obs.kf <- (1/n.selected) * t(resid.kf) %*% solve(H[lx.selected,] %*% state.kf$Sx.prior %*% t(H[lx.selected,]) + diag(Szd.assumed[lx.selected])) %*% (resid.kf)

cat(sprintf(" [KF] chi2 means: state %.2f, obs %.2f on %d (%d) DOFs, RMSE %.2f\n",
            chi2.state.kf,chi2.obs.kf,ndofs.kf,ndofs.patil(state.kf$Sx.post),compute.rmse(state.kf$x.post - truth_condition)))

if(make.pdfs) {
    pdf("common.diagnostics.pdf",width=11,height=8)
}


resid.ens <- obs[lx.selected] - obs.ens.prior
denom <-  chol(solve(H[lx.selected,] %*% state.ens$Sx.prior %*% t(H[lx.selected,])+ diag(Szd.assumed[lx.selected]),tol=1e-8))
resid.z.ens <-  denom %*% resid.ens
cat(sprintf("\n\n[Ens] stats on normalized measurement forecast residuals:\n"))
normality.test(resid.z.ens,known.mean=0,known.sd=1)

denom <-  chol(solve(H[lx.selected,] %*% state.kf$Sx.prior %*% t(H[lx.selected,])+ diag(Szd.assumed[lx.selected]),tol=1e-8))
resid.z.kf <-  denom %*% resid.kf
cat(sprintf("\n\n[KF] stats on normalized measurement forecast residuals:\n"))
normality.test(resid.z.kf,known.mean=0,known.sd=1)

if(make.pdfs) {
    dev.off()
    x.pdf.name='common.x.pdf'
    flux.pdf.name='common.flux.pdf'
    is.pdf.name='common.obs.pdf'
} else {
    x.pdf.name=NULL
    flux.pdf.name=NULL
    is.pdf.name=NULL
}

plot.x.timeseries(ests=list(Truth=list(x=truth_condition),
                            Ensemble=list(x=state.ens$x.post,Sx=state.ens$Sx.post),
                            Analytical=list(x=state.kf$x.post,Sx=state.kf$Sx.post)),
                  pdf.name=x.pdf.name)

plot.flux.timeseries(ests=list(Truth=list(x=truth_condition),
                               Ensemble=list(x=state.ens$x.post,Sx=state.ens$Sx.post),
                               Analytical=list(x=state.kf$x.post,Sx=state.kf$Sx.post)),
                     pdf.name=flux.pdf.name)

plot.is.timeseries(xs=list(Truth=truth_condition,
                           Ensemble=state.ens$x.post,
                           Analytical=state.kf$x.post,
                           Prior=state.kf$x.prior),
                   dataset_names=c("co2_mlo_surface-insitu_1_allvalid",
                                   "co2_brw_surface-insitu_1_allvalid",
                                   "co2_smo_surface-insitu_1_allvalid",
                                   "co2_spo_surface-insitu_1_allvalid",
                                   "co2_lef_tower-insitu_1_allvalid-396magl"),
                   H=H.orig,
                   pdf.name=is.pdf.name)


