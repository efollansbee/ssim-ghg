# Time-stamp: <aj:/Users/andy/Desktop/ssim-ghg/EnKF/L-curve/lcurve.r - 05 Feb 2025 (Wed) 19:18:38 MST>

source("../tools/enkf.r")
source("../tools/progress.bar.r")
source("../tools/find.indir.r")
indir <- find.indir() # potentially replace with trivial yaml access

# explicitly include 1.0 since for this artificial 
# case, it is the perfect configuration
alphas <-unique(sort(c(1,10^seq(-3,3,length.out=40))))

nalphas <- length(alphas)

# allocate storage for results
var.resids <- matrix(NA,nrow=1,ncol=nalphas)
var.flux <- matrix(NA,nrow=1,ncol=nalphas)
                                      
# Load sensitivity matrices (Jacobians)
if(!exists("H.orig")) {
    t0 <- proc.time()[3]
    cat("Loading Jacobians...")
    H.orig <- list()
    load(file.path(indir,"jacobians/trunc_full_jacob_030624_with_dimnames_sib4_4x5_mask.rda"))
    load(file.path(indir,"jacobians/jacob_bgd_021624.rda"))
    H.orig$H <- jacob*(12/44) # Andrew reports units conversion needed
    H.orig$H_fixed <- jacob_bgd[,c(2,3)]
    rm(jacob,jacob_bgd)
    cat(sprintf('%.1fs\n',proc.time()[3]-t0))
}

nobs <- dim(H.orig$H)[1]

Szd.actual <- rep(0.5,nobs) # variance in ppm^2
Szd.assumed <- Szd.actual

nparms <- 22*24 # 22 regions, 24 months

# Real covariance of unknowns
Sx <- diag(rep(1,nparms))

# The generate_ensemble() function just chooses nmemb samples from 
# the multivariate normal distribution with covariance Sx and mean 0.
# We repurpose this to choose a random truth condition here.
truth_condition <- generate_ensemble(Sx=Sx,nmemb=1)
dim(truth_condition) <- c(nparms,1)


# Perturbed observations (because Szd is supplied)
obs <- simulate_observed(H=H.orig$H, x=truth_condition, H_fixed=H.orig$H_fixed, Szd=Szd.actual)
dim(obs) <- c(nobs,1)

# Restrict to nobs randomly sampled subset of measurements. Could use
# obs_catalog or row.names of H to do more systematically-chosen
# subsets.
#
# The L-curve is typical of /under/-constrained problems, so we
# restrict to just a small number of observations.
nobs <- 200
lx <- sample(x=1:length(obs),size=nobs)
obs <- obs[lx,]
Szd.assumed <- Szd.assumed[lx]
Szd.actual <- Szd.actual[lx]
H <- H.orig$H[lx,]
H_fixed <- H.orig$H_fixed[lx,]

obs_fixed <- apply(H_fixed,c(1),sum)

pb <- progress.bar.start(nalphas)

ialpha <- 0
for (alpha in alphas) {
  ialpha <- ialpha + 1

  Sx.prior <- alpha*Sx 
  x.prior <- rep(0,nparms)
  
  
  # Kalman filter measurement update
  kf <- kf_meas_update(x=x.prior,Sx=Sx.prior,H=H,z=obs-obs_fixed,
                       Sz=diag(Szd.assumed))
  
  dobs <- matrix(simulate_observed(H=H, x=kf$x,H_fixed=H_fixed) - obs,nrow=nobs,ncol=1)
  dx <- matrix(kf$x - x.prior,nrow=nparms,ncol=1)
  
  var.resids[ialpha] <- var(as.vector(dobs))
  var.flux[ialpha] <- var(as.vector(kf$x))

  pb <- progress.bar.print(pb,ialpha)
} # alpha

progress.bar.end(pb)

# save(chi2.resids,chi2.flux,var.flux,alphas,file="lcurve.rda")

options(jupyter.plot_scale=1,repr.plot.height=7,repr.plot.width=14)
layout(matrix(1:2,nrow=1))

lx <- c(which.min(alphas),which(alphas==1),which.max(alphas))
plot(x=var.resids, y=var.flux,xlab="variance of the obs resids",ylab="variance of optimized parameters",main="Linear scale")
points(x=var.resids[lx], y=var.flux[lx],pch=20,col='red')
text(x=var.resids[lx], y=var.flux[lx],labels=alphas[lx],pos=4,xpd=NA)

plot(x=var.resids, y=var.flux,xlab="variance of the obs resids",ylab="variance of optimized parameters",log='xy',main="Log scale")
text(x=var.resids[lx], y=var.flux[lx],labels=alphas[lx],pos=4,xpd=NA)
points(x=var.resids[lx], y=var.flux[lx],pch=20,col='red')
