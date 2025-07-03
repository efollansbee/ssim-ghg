#-- SETUP
setup = function()
{
       
###############################################
#-- Required Libraries
###############################################

require(ncdf4,warn.conflicts = FALSE)
require(plyr,warn.conflicts = FALSE)
require(dplyr,warn.conflicts = FALSE)
require(parallel,warn.conflicts = FALSE)
require(ggplot2,warn.conflicts = FALSE)
require(abind,warn.conflicts = FALSE)
require(Matrix,warn.conflicts = FALSE)
require(lattice,warn.conflicts = FALSE)
require(memuse,warn.conflicts = FALSE)
require(EnvStats,warn.conflicts = FALSE)
require(gridExtra,warn.conflicts = FALSE)
require(mvtnorm,warn.conflicts = FALSE)
require(plotly,warn.conflicts = FALSE)

    
###############################################
#-- Load Code
##############################################
source(file.path(Rcode_dir,"plot_concentrations.R"))
source(file.path(Rcode_dir,"inversion_032024.R"))
source(file.path(Rcode_dir,"write_inversion_2_netcdf_032024.R"))
source(file.path(Rcode_dir,"generate_transcom_flux_ensemble_from_inversion.R"))

########################
#--  Detect Cores
########################
print(paste("Num CPUs:",detectCores(),"cores"))
memuse::Sys.meminfo()
}
    

#-- Create the true flux field from prior and true state vector

pull_true_transcom_flux = function(prior_flux_file,state_true)
{
  con = nc_open(prior_flux_file)
  NEE = ncvar_get(con,"NEE")
  nc_close(con)

  NEE_1x1 = aaply(NEE,3,.fun=function(x){expand_5x4_2_1x1(x)}) %>%
           aperm(c(2,3,1))

  zz = aaply(NEE_1x1,3,.fun=function(x){grid2transcom(x,file_location=data_dir)})
  yy = matrix(state_true,nrow=24,byrow=FALSE)
  #qq = zz[,-1]*yy
  #transcom_fluxes_real = apply(qq,2,sum)  *12/44   *30.5*3600*24*1e3*0.5 * 1e-15 # ~ PgC/yr adjustment to prior
  
  #-- These are prior fluxes
  monthly_prior = zz[,-1] *12/44   *30.5*3600*24*1e3 * 1e-15 
  annual_2yr_prior = apply(monthly_prior,2,sum) * 0.5  # ~ PgC/yr adjustment to prior
  
  #-- These are prior fluxes times true state
  qq = zz[,-1]*yy *12/44   *30.5*3600*24*1e3 * 1e-15 
  transcom_fluxes_real = apply(qq,2,sum) * 0.5  # ~ PgC/yr adjustment to prior
  ret = list(monthly=qq,annual_2yr=transcom_fluxes_real,monthly_prior=monthly_prior,annual_2yr_prior=annual_2yr_prior)
  return(ret)
}

generate_observations = function(H,H_bgd,state_vector,err_obs=NULL){
  obs = H %*% (1+state_vector) #+ apply(H_bgd[,c(2,3)],1,sum)
  if(!is.null(err_obs)){obs = obs + rnorm(err_obs)}
  return(obs)
}


grid2transcom = function(mat,model.grid.x=1,model.grid.y=1,transcom=TRUE,
                         file_location=data_dir)
{
  tfile = paste(file_location,"transcom/TRANSCOM_mask_GEOS_Chem_",model.grid.y,"x",model.grid.x,".nc",sep="")
  afile = paste(file_location,"areas/area_",model.grid.x,"x",model.grid.y,".nc",sep="")

  grid_areas = load.ncdf(afile,vars=c("grid_cell_area"))$grid.cell.area
  
  transfil = load.ncdf(tfile)
  if(model.grid.x==1 & model.grid.y==1){
    regions = transfil$mask64
  }else{
      for(i in 1:23){
        if(i==1){ret = transfil[[1]]}else{
          transfil[[i]][transfil[[i]]==1] = i
          ret = ret + transfil[[i]]}
      }
  regions=ret
  }

  #-- weird 3.63216561952993e-40 value in transcom files
  regions[regions > 0 & regions < 1e-6] = 0
  
 if(transcom)
  {
    transcom.agg  = aggregate(as.vector(mat * grid_areas),list(as.vector(regions)),sum)
    transcom.agg.vect = transcom.agg[,2] #as.vector(unlist(transcom.agg))
    
    return(transcom.agg.vect)
  }else{
    return(mat * grid_areas)
  }
}



transcom2grid = function(vect=rep(2,23),model.grid.x=1,model.grid.y=1,transcom=TRUE,
                         file_location=data_dir)
{
  tfile = paste(file_location,"transcom/TRANSCOM_mask_GEOS_Chem_",model.grid.y,"x",model.grid.x,".nc",sep="")
  afile = paste(file_location,"areas/area_",model.grid.x,"x",model.grid.y,".nc",sep="")
  
  if(!file.exists(tfile)){stop(paste("file:",tfile," doesn't exist"))}
  if(!file.exists(afile)){stop(paste("file:",afile," doesn't exist"))}
  
  grid_areas = load.ncdf(afile,vars=c("grid_cell_area"))$grid.cell.area

  transfil = load.ncdf(tfile)
  if(model.grid.x==1 & model.grid.y==1){
    regions = transfil$mask64
  }else{
    for(i in 1:23){
      if(i==1){ret = transfil[[1]]}else{
        transfil[[i]][transfil[[i]]==1] = i
        ret = ret + transfil[[i]]}
    }
    regions=ret
  }
  
  regions[regions < 0.001] = 0
  
  for(i in 1:22)
  {
    regions[regions == i] = vect[i]
  }
  
 return(regions)
}

#-- checking effect of truncating big and small scalings
tm = function(x,low,high){
  x[x < low] = low;x[x > high] = high
  return(x)
}

expand_5x4_2_1x1 = function(x){
  v = as.vector(x)
  m = matrix(as.vector(matrix(rep(v,5),nrow=5,byrow=TRUE)),nrow=360)
  m2 = t(m)
  v2 = as.vector(m2)
  m3 = matrix(as.vector(matrix(rep(v2,4),nrow=4,byrow=TRUE)),nrow=184)
  m4 = t(m3)
  return(m4[,3:182])
}

#-- EXTRA CODE
#gridded_5x4_state = transcom2grid(as.vector(inv_object$x_hat),model.grid.x = 5,model.grid.y = 4,
#                                  file_location = tr_dir)
#
#gridded_5x4_mean_flux = aaply(NEE,c(3),.fun=function(x){x*gridded_5x4_state}) %>%
#  aperm(c(2,3,1))
#
#gridded_5x4_state_samples = aaply(samps,1,.fun=function(x){transcom2grid(x,model.grid.x = 5,model.grid.y = 4,file_location = tr_dir)}) %>%
#  aperm(c(2,3,1))
#
#gridded_5x4_flux_samples = array(NA,dim=c(dim(gridded_5x4_mean_flux),sample_number))
#
#for(i in 1:sample_number){
#  tmp_state = gridded_5x4_state_samples[,,i]
#  gridded_5x4_flux_samples[,,,i] = aaply(gridded_5x4_mean_flux,3,.fun=function(x){x*tmp_state}) %>%
#    aperm(c(2,3,1))
#}

# Time-stamp: <andyj.cmdl.noaa.gov:/ct/tools/R/load.ncdf4.r - 07 Nov 2011 (Mon) 15:26:41 MST>

interpret.udunits.time <- function(vals,unitstring,tz="UTC") {
  retval <- list()
  retval$is.time <- FALSE
  if(length(grep("^DAYS since", unitstring,ignore.case=TRUE))==1) {
    retval$is.time <- TRUE
    # length of string is a bad parsing heuristic, but works
    # for everything so far.
    #     21 chars if "days since 1900-01-01"
    #     19 chars if "days since 1900-1-1"
    retval$vals <- as.POSIXct(substr(unitstring,11,nchar(unitstring)),
                              tz=tz) + vals*86400
    retval$tz <- tz # UTC tzone is a presumption
  }
  if(length(grep("^SECONDS since", unitstring,ignore.case=TRUE))==1) {
    retval$is.time <- TRUE
    retval$vals <- as.POSIXct(substr(unitstring,15,nchar(unitstring)),
                              tz=tz) + vals
    retval$tz <- tz # UTC tzone is a presumption
  }
  if(length(grep("^HOURS since", unitstring,ignore.case=TRUE))==1) {
    retval$is.time <- TRUE
    retval$vals <- as.POSIXct(substr(unitstring,12,nchar(unitstring)),
                              tz=tz) + vals*3600
    retval$tz <- tz # UTC tzone is a presumption
  }
  if(length(grep("^decimal date", unitstring,ignore.case=TRUE))==1) {
    retval$is.time <- TRUE
    retval$vals <- decimal.to.POSIX(vals)
    retval$tz <- tz # UTC tzone is a presumption
  }
  return(retval)
}


load.ncdf <- function (ncname, lowercase = TRUE, dims = TRUE, attribs = NULL, 
                       vars = NULL, verbose = FALSE) 
{
  if (verbose) {
    cat(sprintf("[load.ncdf]  Called with ncname = %s,\n", 
                ncname))
    cat(sprintf("[load.ncdf]  and lowercase = %s,\n", lowercase))
    cat(sprintf("[load.ncdf]  and dims = %s,\n", dims))
    if (!is.null(attribs)) {
      cat(sprintf("[load.ncdf]  and attribs = %s,\n", attribs))
    }
    if (!is.null(vars)) {
      cat(sprintf("[load.ncdf]  and vars = %s.\n", vars))
    }
  }
  library(ncdf4)
  retval <- list()
  nc <- nc_open(ncname)
  i <- 0
  if (is.null(vars)) {
    vars <- names(nc$var)
  }
  for (var in vars) {
    if (!var %in% names(nc$var)) {
      warning(sprintf("variable \"%s\" not in file \"%s\".", 
                      var, ncname))
      next
    }
    i <- i + 1
    retval[[i]] <- ncvar_get(nc, var)
    vunits <- ncatt_get(nc, var, attname = "units")
    if (vunits$hasatt) {
      attributes(retval[[i]])$units <- vunits$value
      time.checked <- interpret.udunits.time(vals = retval[[i]], 
                                             unitstring = vunits$value)
      if (time.checked$is.time) {
        retval[[i]] <- time.checked$vals
        attributes(retval[[i]])$tzone <- time.checked$tz
      }
    }
    if (lowercase) {
      names(retval)[i] <- gsub("_", ".", var)
    }
    else {
      names(retval)[i] <- var
    }
  }
  if (dims) {
    retval$dim <- nc$dim
    for (d in names(nc$dim)) {
      ddot <- d
      if (lowercase) {
        ddot <- gsub("_", ".", ddot)
      }
      retval[[ddot]] <- nc$dim[[d]]$vals
      unitstring <- nc$dim[[d]]$units
      time.checked <- interpret.udunits.time(nc$dim[[d]]$vals, 
                                             nc$dim[[d]]$units)
      if (time.checked$is.time) {
        retval[[ddot]] <- time.checked$vals
        attributes(retval[[ddot]])$tzone <- time.checked$tz
      }
    }
  }
  for (att in attribs) {
    target.name <- att
    if (lowercase) {
      target.name <- gsub("_", ".", target.name)
    }
    if (target.name %in% names(attributes(retval))) {
      target.name <- paste(target.name, "2", sep = "")
    }
    attributes(retval)[[target.name]] <- NULL
    att.val <- ncatt_get(nc, 0, att)
    if (att.val$hasatt) {
      attributes(retval)[[target.name]] <- att.val$value
    }
    else {
      warning(sprintf("[load.ncdf] Cannot find attribute \"%s\" in file \"%s\".", 
                      att, nc$filename))
    }
  }
  nc_close(nc)
  if (lowercase) {
    names(retval) <- tolower(names(retval))
  }
  return(retval)
}

my.col = function (nl) 
{
  rgb(c(seq(0, 1, length = floor(nl/2)), rep(1, ceiling(nl/2))), 
      c(abs(sin((0:(nl - 1)) * pi/(nl - 1)))), c(c(rep(1, ceiling(nl/2)), 
                                                   seq(1, 0, length = floor(nl/2)))))
}

#-- Plotting of monthly flux averages for 9/2014 - 8/2016 for chosen Transcom Regions
#- need to make general so next line can take vector of region #s

plot_transcom_flux_by_month = function(ret){
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #transcom_region_plot_monthly = 1
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  #-- Plot transcom regions 
  for(transcom_region_plot_monthly in 1:22){
    subset = ret$prior_tr_monthly[,,transcom_region_plot_monthly]
    prior_df = cbind(FLUX=as.vector(subset),KIND=rep("Prior",length(subset)),MONTH=rep(1:24,(dim(subset)[2])))
    
    subset = ret$post_tr_monthly[,,transcom_region_plot_monthly]
    post_df = cbind(FLUX=as.vector(subset),KIND=rep("Post",length(subset)),MONTH=rep(1:24,(dim(subset)[2])))
    
    combined_df  = as.data.frame(rbind(post_df,prior_df))
    
    combined_df$KIND = factor(combined_df$KIND)
    combined_df$FLUX = as.numeric(combined_df$FLUX)
    combined_df$MONTH = factor(combined_df$MONTH,levels=1:24)
    require(ggplot2)
    
    options(repr.plot.width=20, repr.plot.height=8)
    
    #-- This plot is monthly avg regional flux (in units of PgC/yr), plotted for 2 years
    
    subset_truth = ret$transcom_fluxes_real_monthly_avg[,transcom_region_plot_monthly]
    new_data <- data.frame(MONTH =1:24, FLUX=c(subset_truth), KIND=rep("Truth",24))
    
    g = ggplot(combined_df, aes(x=MONTH, y= FLUX, fill=KIND)) +
      geom_boxplot(width=0.5) +   # outlier.shape = NA
      scale_fill_manual(values=c("red", "blue","green")) +
      ylab("PgC/year") +
      geom_point(data=new_data, aes(x=MONTH, y=FLUX, fill=KIND), color="black",bg="green", size=5, pch=21) +
      labs(title=transcom_names[transcom_region_plot_monthly]) + theme(axis.text=element_text(size=12),                                                               axis.title=element_text(size=14,face="bold"),title=element_text(size=16),legend.text=element_text(size=14)) +
      scale_x_discrete(breaks=c("1","2","3","4","5","6","7","8","9","10","11","12",
                               "13","14","15","16","17","18","19","20","21","22","23","24"),
                       labels=dts) + theme(axis.text.x = element_text(angle=70,size=15))
    
    
    
    plot(g)
    
  }
  
  
  #-- Plot global sum 
  subset = ret$prior_tr_monthly_global
  prior_df = cbind(FLUX=as.vector(subset),KIND=rep("Prior",length(subset)),MONTH=rep(1:24,(dim(subset)[2])))
  
  subset = ret$post_tr_monthly_global
  post_df = cbind(FLUX=as.vector(subset),KIND=rep("Post",length(subset)),MONTH=rep(1:24,(dim(subset)[2])))
  
  combined_df  = as.data.frame(rbind(post_df,prior_df))
  
  combined_df$KIND = factor(combined_df$KIND)
  combined_df$FLUX = as.numeric(combined_df$FLUX)
  combined_df$MONTH = factor(combined_df$MONTH,levels=1:24)
  require(ggplot2)
  
  options(repr.plot.width=20, repr.plot.height=8)
  
  #-- This plot is monthly avg regional flux (in units of PgC/yr), plotted for 2 years
  
  subset_truth = ret$transcom_fluxes_real_monthly_avg_global
  new_data <- data.frame(MONTH =1:24, FLUX=c(subset_truth), KIND=rep("Truth",24))
  
  g = ggplot(combined_df, aes(x=MONTH, y= FLUX, fill=KIND)) +
    geom_boxplot(width=0.5) +   # outlier.shape = NA
    scale_fill_manual(values=c("red", "blue","green")) +
    ylab("PgC/year") +
    geom_point(data=new_data, aes(x=MONTH, y=FLUX, fill=KIND), color="black",bg="green", size=5, pch=21) +
    labs(title="Global CO2 Flux") + theme(axis.text=element_text(size=12),
                                          axis.title=element_text(size=14,face="bold"),title=element_text(size=16),legend.text=element_text(size=14))
  
  
  
  plot(g)
}

#-- Use this to control plot size
plot_timeseries_flux_bytranscom = function(ret,include_prior_ocn_land=FALSE,include_fossil_biofires=FALSE)
{

  mat_table = data.frame(NUMBER=1:22,NAME=transcom_names)
  
  plt_df = ret$full_df
  plt_df$REGION = mat_table$NAME[plt_df$REGION]
 
  
  #-- Plot global sum 
  subset = ret$prior_tr_monthly_global
  subset = apply(subset,c(2),sum)  * 0.5   # 0.5 is to take average over 24 months ~ PgC/yr 
  prior_df = cbind(FLUX=as.vector(subset),KIND=rep("Prior",length(subset)))
  
  subset = ret$post_tr_monthly_global
  subset = apply(subset,c(2),sum)  * 0.5   # 0.5 is to take average over 24 months ~ PgC/yr 
  post_df = cbind(FLUX=as.vector(subset),KIND=rep("Post",length(subset)))  
  
  combined_df  = as.data.frame(rbind(post_df,prior_df))
  
  combined_df$KIND = factor(combined_df$KIND)
  combined_df$FLUX = as.numeric(combined_df$FLUX)
  combined_df = cbind(combined_df,REGION=rep("Global",dim(combined_df)[1]))
  combined_df$REGION = as.factor(combined_df$REGION)
  
  plt_df = rbind(plt_df,combined_df)
  
  smp_mean = aggregate(plt_df$FLUX,list(plt_df$KIND,plt_df$REGION),mean)
  
  ###############################################
  #--   FOR REGIONAL/TRANSCOM ESTIMATES     -- -#

  options(repr.plot.width=18, repr.plot.height=8)
  
  g = ggplot(plt_df[plt_df$REGION != "Global",], aes(x=REGION, y= FLUX, fill=KIND)) +
    geom_boxplot(width=0.5) +   # outlier.shape = NA
    ylab("PgC/year") +
    scale_fill_manual(values=c("red", "blue","green"))
  
  new_data <- data.frame(REGION = c(transcom_names,"Global"), FLUX=c(ret$transcom_fluxes_real_annual_avg,
                                                                     sum(ret$transcom_fluxes_real_annual_avg)), KIND=rep("Truth",23))
  
  h = g + geom_point(data=new_data[new_data$REGION != "Global",], aes(x=REGION, y=FLUX, fill=KIND), color="black",bg="green", size=5, pch=21) +
    ylab("PgC/year") +
    theme(axis.text.x = element_text(angle=70,size=15))
  
  plot(h)
  ###############################################
  
  
  ###############################################
  #--   FOR GLOBAL ESTIMATES     -- -#  

  #options(repr.plot.width=4, repr.plot.height=8)
  
  g2 = ggplot(plt_df[plt_df$REGION == "Global",], aes(x=REGION, y= FLUX, fill=KIND)) +
    geom_boxplot(width=0.5) +   # outlier.shape = NA
    ylab("PgC/month") +
    scale_fill_manual(values=c("red", "blue","green"))
  
  
  new_data <- data.frame(REGION = c(transcom_names,"Global"), FLUX=c(ret$transcom_fluxes_real_annual_avg,
                                                                     sum(ret$transcom_fluxes_real_annual_avg)), KIND=rep("Truth",23))
  
  h2 = g2 + geom_point(data=new_data[new_data$REGION == "Global",], aes(x=REGION, y=FLUX, fill=KIND), color="black",bg="green", size=5, pch=21) +
    ylab("PgC/month") +
    theme(axis.text.x = element_text(angle=70,size=15))
  
  plot(h2)
  ###############################################    
  
  
  
  prnt_data = new_data[,c(1,2)]
  names(prnt_data)[2] = "TRUTH"
  prnt_data$POST = NA
  prnt_data$PRIOR = NA  
  
  for(i in 1:dim(prnt_data)[1]){
    prnt_data$POST[i] = smp_mean[smp_mean$Group.1== "Post" & smp_mean$Group.2==prnt_data[i,1] , 3]
    prnt_data$PRIOR[i] = smp_mean[smp_mean$Group.1== "Prior" & smp_mean$Group.2==prnt_data[i,1] , 3]
  }
  
  prnt_data$TRUTH = round(prnt_data$TRUTH,2)
  prnt_data$PRIOR = round(prnt_data$PRIOR,2)
  prnt_data$POST = round(prnt_data$POST,2)
  
  print("Sample Means for Regions:")
  print(prnt_data)
  
}

#-- FORM DATA FRAME OF PRIOR/POSTERIOR USING SAMPLES FROM NETCDF FILES, FOR PLOTTING
#-- this needs transparent parapply wrapper on the monthly grid2transcom() to speed up (e.g. line 2 below)
organize_data_for_plotting = function(prior_flux_netcdf="/Users/aschuh/test_output_prior/gridded_fluxes.nc4",
                                      posterior_flux_netcdf="/Users/aschuh/test_output_posterior/gridded_fluxes.nc4",
                                      max.cores=4)
{
  cores_used = min(max.cores,detectCores() - 2)
  cl <- makeCluster(cores_used, type = "FORK")
   print(paste("you have",detectCores(),"cores to work with, working with ",cores_used))
  d = load.ncdf(posterior_flux_netcdf)
  post_tr_monthly <- parApply(cl,d$flux.samples,c(3,4),FUN=function(x){grid2transcom(x,file_location=data_dir)}) %>% aperm(c(2,3,1)) 
  #post_tr_monthly = aaply(d$flux.samples,c(3,4),.fun=function(x){grid2transcom(x)})
  post_tr_monthly = post_tr_monthly[,,2:23]  #drop transcom=0
  post_tr_monthly_global = apply(post_tr_monthly,c(1,2),sum)
  post_tr_annual = apply(post_tr_monthly,c(2,3),sum)*30.5*3600*24*1e3*0.5 * 1e-15# ~ gC/yr
  post_df = cbind(FLUX=as.vector(t(post_tr_annual)),KIND=rep("Post",length(post_tr_annual)),REGION=rep(1:22,dim(post_tr_annual)[1]))
  
  e = load.ncdf(prior_flux_netcdf)
  prior_tr_monthly <- parApply(cl,e$flux.samples,c(3,4),FUN=function(x){grid2transcom(x,file_location=data_dir)}) %>% aperm(c(2,3,1))
  #prior_tr_monthly = aaply(e$flux.samples,c(3,4),.fun=function(x){grid2transcom(x)})
  prior_tr_monthly = prior_tr_monthly[,,2:23] 
  prior_tr_monthly_global = apply(prior_tr_monthly,c(1,2),sum)
  prior_tr_annual = apply(prior_tr_monthly,c(2,3),sum)*30.5*3600*24*1e3*0.5 * 1e-15 # ~ ~ PgC/yr adjustment to prior, 0.5 is 2 yr -> 1yr
  prior_df = cbind(FLUX=as.vector(t(prior_tr_annual)),KIND=rep("Prior",length(prior_tr_annual)),REGION=rep(1:22,dim(prior_tr_annual)[1]))
  
  full_df  = as.data.frame(rbind(post_df,prior_df))
  
  stopCluster(cl)
  
  full_df$KIND = factor(full_df$KIND)
  full_df$FLUX = as.numeric(full_df$FLUX)
  full_df$REGION = factor(full_df$REGION,levels=1:22)
  
  
  #-- This calculate the "true" fluxes from original prior fluxes
  transcom_fluxes_real = pull_true_transcom_flux(prior_flux_file=file.path(data_dir,"priors/prior_SiB4.nc"),state_true=state_vector_true)
  transcom_fluxes_real_annual_avg = transcom_fluxes_real$annual_2yr
  transcom_fluxes_real_monthly_avg = transcom_fluxes_real$monthly
  transcom_fluxes_real_monthly_avg_global = apply(transcom_fluxes_real$monthly,c(1),sum)
  
  output = list(prior_tr_monthly=prior_tr_monthly,prior_tr_monthly_global=prior_tr_monthly_global,prior_tr_annual=prior_tr_annual,
                post_tr_monthly=post_tr_monthly,post_tr_monthly_global=post_tr_monthly_global,post_tr_annual=post_tr_annual,
                full_df=full_df,transcom_fluxes_real_annual_avg=transcom_fluxes_real_annual_avg,
                transcom_fluxes_real_monthly_avg=transcom_fluxes_real_monthly_avg,
                transcom_fluxes_real_monthly_avg_global=transcom_fluxes_real_monthly_avg_global)
  return(output)
}

ar_covariance <- function(n, rho, variance = 1) {
  variance * toeplitz(rho ^ (0 : (n - 1)))
}

plot_inversion_correlations = function(org_data)
{
  dts = c("Sep 2014","Oct 2014","Nov 2014","Dec 2014","Jan 2015","Feb 2015","Mar 2015","Apr 2015",
          "May 2015","Jun 2015","Jul 2015","Aug 2015","Sep 2015","Oct 2015","Nov 2015","Dec 2015",
          "Jan 2016","Feb 2016","Mar 2016","Apr 2016",
          "May 2016","Jun 2016","Jul 2016","Aug 2016")
  
  post_cor_flux = cor(org_data$post_tr_annual)
  prior_cor_flux = cor(org_data$prior_tr_annual)
  dimnames(post_cor_flux) = list(transcom_names,transcom_names)
  dimnames(prior_cor_flux) = list(transcom_names,transcom_names)
  
  post_time_cor_flux = cov2cor(cor(t(org_data$post_tr_monthly_global)))
  prior_time_cor_flux = cov2cor(cor(t(org_data$prior_tr_monthly_global)))
  dimnames(prior_time_cor_flux) = list(dts,dts)
  dimnames(post_time_cor_flux) = list(dts,dts)
  
  diag(post_cor_flux) = NA
  diag(prior_cor_flux) = NA
  diag(post_time_cor_flux) = NA
  diag(prior_time_cor_flux) = NA

  rng1 = max(abs(c(as.vector(prior_cor_flux),c(as.vector(post_cor_flux)))),na.rm=TRUE)
  rng2 = max(abs(c(as.vector(prior_time_cor_flux),c(as.vector(post_time_cor_flux)))),na.rm=TRUE)
  
  options(repr.plot.width = 20, repr.plot.height = 20)
  
  p1 = levelplot(prior_cor_flux,col.regions=my.col(40),
                 at=seq(-rng1,rng1,length=40),main="Prior Correlation in Avg Annual Flux between Transcom Region ",scales=list(x=list(rot=60)),
                 xlab="",ylab="")
  
  p2 = levelplot(post_cor_flux,col.regions=my.col(40),at=seq(-rng1,rng1,length=40),main="Posterior Correlation in Avg Annual Flux between Transcom Region ",scales=list(x=list(rot=60)),
                 xlab="",ylab="")
  
  p3 = levelplot(prior_time_cor_flux,col.regions=my.col(40),at=seq(-rng2,rng2,length=40),main="Prior Correlation in time in global flux",scales=list(x=list(rot=60)),
                 xlab="",ylab="")
  
  p4 = levelplot(post_time_cor_flux,col.regions=my.col(40),at=seq(-rng2,rng2,length=40),main="Posterior Correlation in time in global flux ",scales=list(x=list(rot=60)),
                 xlab="",ylab="")
  
  marrangeGrob(list(p1,p2,p3,p4),nrow=2,ncol=2)
}


plot_inversion_correlations_by_transcom = function(org_data)
{
  dts = c("Sep 2014","Oct 2014","Nov 2014","Dec 2014","Jan 2015","Feb 2015","Mar 2015","Apr 2015",
          "May 2015","Jun 2015","Jul 2015","Aug 2015","Sep 2015","Oct 2015","Nov 2015","Dec 2015",
          "Jan 2016","Feb 2016","Mar 2016","Apr 2016",
          "May 2016","Jun 2016","Jul 2016","Aug 2016")
  
  plt = function(i)
  {
  transcom_subset_post = t(org_data$post_tr_monthly[,,i])
  transcom_subset_prior = t(org_data$prior_tr_monthly[,,i])
  
  post_cor_flux = cor(transcom_subset_post)
  prior_cor_flux = cor(transcom_subset_prior)
  dimnames(post_cor_flux) = list(dts,dts)
  dimnames(prior_cor_flux) = list(dts,dts)
  
  diag(post_cor_flux) = NA
  diag(prior_cor_flux) = NA
  
  rng1 = max(abs(c(as.vector(prior_cor_flux),c(as.vector(post_cor_flux)))),na.rm=TRUE)
  
  options(repr.plot.width = 20, repr.plot.height = 20)
  
  p1 = levelplot(prior_cor_flux,col.regions=my.col(20),at=seq(-rng1,rng1,length=20),main=paste("Prior Correlation in Month to Month flux for",transcom_names[i]),scales=list(x=list(rot=60)),
                 xlab="",ylab="")
  
  p2 = levelplot(post_cor_flux,col.regions=my.col(20),at=seq(-rng1,rng1,length=20),main=paste("Posterior Correlation in Month to Month flux for",transcom_names[i]),scales=list(x=list(rot=60)),
                 xlab="",ylab="")
  
  marrangeGrob(list(p1,p2),nrow=1,ncol=2)
  }
  
  lapply(1:22,FUN=plt)

}

dts = c("Sep 2014","Oct 2014","Nov 2014","Dec 2014","Jan 2015","Feb 2015","Mar 2015","Apr 2015",
        "May 2015","Jun 2015","Jul 2015","Aug 2015","Sep 2015","Oct 2015","Nov 2015","Dec 2015",
        "Jan 2016","Feb 2016","Mar 2016","Apr 2016",
        "May 2016","Jun 2016","Jul 2016","Aug 2016")

#-- Transcom region labels for regions 1:22
transcom_names = c("North American Boreal    ", "North American Temperate ",
                   "South American Tropical  ", "South American Temperate ",
                   "Northern Africa          ", "Southern Africa          ",
                   "Eurasia Boreal           ", "Eurasia Temperate        ",
                   "Tropical Asia            ", "Australia                ",
                   "Europe                   ", "North Pacific Temperate  ",
                   "West Pacific Tropical    ", "East Pacific Tropical    ",
                   "South Pacific Temperate  ", "Northern Ocean           ",
                   "North Atlantic Temperate ", "Atlantic Tropical        ",
                   "South Atlantic Temperate ", "Southern Ocean           ",
                   "Indian Tropical          ", "South Indian Temperate   ")

parse.obspack_id <- function(x) {
  info <- data.frame(t(matrix(as.vector(unlist(strsplit(x,'~'))),nrow=3)),stringsAsFactors=FALSE)
  names(info) <- c("obspack_name","dataset_name","obspack_number")
  this.tempval <-        data.frame(t(matrix(as.vector(unlist(strsplit(info$dataset_name,"_"))),nrow=5)),stringsAsFactors=FALSE)
  info$species <- this.tempval[,1]
  info$site <- this.tempval[,2]
  info$project <- this.tempval[,3]
  info$labcode <- this.tempval[,4]
  info$selection <- this.tempval[,5]
  return(info)
} 


plot_Jacobian_cols_observations = function(transcom_region=1:3,month=1:3,
                                           plot_sum=FALSE,
                                           site_strings=c("brw_surface-insitu","mlo_surface-flask","smo_surface-insitu",
                                                          "co2_spo_surface-flask",
                                                          "lef","wkt","wbi","nwr","hun","cgo","cpt"),
                                          plot_fossil_instead=FALSE,
                                          plot_fires_instead=FALSE)
{
  month_string = c("2014-09","2014-10","2014-11","2014-12","2015-01","2015-02","2015-03","2015-04","2015-05","2015-06","2015-07","2015-08",
                   "2015-09","2015-10","2015-11","2015-12","2016-01","2016-02","2016-03","2016-04","2016-05","2016-06","2016-07","2016-08")
  cnt = 1
  state_ind = vector()
  for(i in 1:length(transcom_region)){
    for(j in 1:length(month)){
      if(transcom_region[i] <= 11){state_ind[cnt] = paste("nee_regionRegion",pad(transcom_region[i],width=2,fill="0"),"_month",month_string[month[j]],sep="")}
      if(transcom_region[i] > 11){state_ind[cnt] = paste("ocean_regionRegion",pad(transcom_region[i],width=2,fill="0"),"_month",month_string[month[j]],sep="")} 
      cnt = cnt + 1
    }
  }
  #print(state_ind)
  if(!plot_fossil_instead & !plot_fires_instead){
    if(plot_sum)
    {
      df = data.frame(ID=obs_catalog$ID,
                  DATE=obs_catalog$DATE,
                  VALUE = apply(H[,state_ind],1,sum))
    }else{
      df = data.frame(ID=obs_catalog$ID,
                      DATE=obs_catalog$DATE,
                      VALUE = H[,state_ind])      
    }
  }
  if(!plot_fossil_instead & plot_fires_instead){
      df = data.frame(ID=obs_catalog$ID,
                  DATE=obs_catalog$DATE,
                  VALUE = H_bgd[,2])
  }
  if(plot_fossil_instead & !plot_fires_instead){
      df = data.frame(ID=obs_catalog$ID,
                  DATE=obs_catalog$DATE,
                  VALUE = H_bgd[,3])
  }
  if(plot_fossil_instead & plot_fires_instead){
      df = data.frame(ID=obs_catalog$ID,
                  DATE=obs_catalog$DATE,
                  VALUE = H_bgd[,2] + H_bgd[,3])
  }
                  #VALUE=c(ret2$posterior$inputs$obs + pr,
                  #        ret2$prior$outputs$modeled_obs + pr,
                  #        ret2$posterior$outputs$modeled_obs + pr),
                  #TYPE=c(rep("OBS",dim(obs_catalog)[1]),
                  #       rep("PRIOR",dim(obs_catalog)[1]),
                  #       rep("POSTERIOR",dim(obs_catalog)[1])))
  
  
  #df$TYPE =factor(df$TYPE)
  #df$VALUE = as.numeric(df$VALUE)
  require(reshape2)
  
  for(i in 1:length(site_strings))
  {
    ind = 1:length(df$ID) %in% grep(site_strings[i],df$ID)
  
    df2 = df[ind,]
    df3 <- melt(df2[,-1] ,  id.vars = 'DATE', variable.name = 'series')
    
    
    options(repr.plot.width=20, repr.plot.height=8)
    

    g =  ggplot(df3, aes(x = DATE, y = value)) + 
      #scale_color_manual(values = c("blue","red","black")) +
      #OLDgeom_point(size = 0.5) + 
      geom_point(size = 0.5,aes(colour=series)) +
      #geom_smooth(method = "lm", formula = y ~ poly(x, 12), se = FALSE) +
      labs(title=paste(site_strings[i],"Time series"),x ="Date", y = "CO2") + 
      theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
            title=element_text(size=16),legend.text=element_text(size=14)) 
       
      if(dim(df2)[2]>5){ 
        print("w/ more than 5 contributions, suppressing legend")
        g = g + theme(legend.position = "none")
       }
    plot(g)
  }
}

plot_Jacobian_rows_fluxes = function(row_number){

  load(file.path(data_dir,"obs/obs_catalog_042424_unit_pulse_hour_timestamp_witherrors_withdates.rda")) 

   print(obs_catalog[row_number,c("TYPE","ID","LON","LAT","DATE")])
    
   month_string = c("2014-09","2014-10","2014-11","2014-12","2015-01","2015-02","2015-03","2015-04","2015-05","2015-06","2015-07","2015-08",
                 "2015-09","2015-10","2015-11","2015-12","2016-01","2016-02","2016-03","2016-04","2016-05","2016-06","2016-07","2016-08")
  
  for(i in 1:22)
  {
    transcom_regs_dims = sapply(dimnames(H)[[2]],FUN=function(x){gsub("regionRegion","",strsplit(x,"_")[[1]][2])})
    transcom_ind = transcom_regs_dims == pad(i,width=2,fill="0")
    
    df2 = H[row_number,transcom_ind]
    ord = order(names(df2))
    df2 = df2[ord]
    
    options(repr.plot.width=20, repr.plot.height=8)

#-- Transcom region labels for regions 1:22
transcom_names = c("North American Boreal    ", "North American Temperate ",
                   "South American Tropical  ", "South American Temperate ",
                   "Northern Africa          ", "Southern Africa          ",
                   "Eurasia Boreal           ", "Eurasia Temperate        ",
                   "Tropical Asia            ", "Australia                ",
                   "Europe                   ", "North Pacific Temperate  ",
                   "West Pacific Tropical    ", "East Pacific Tropical    ",
                   "South Pacific Temperate  ", "Northern Ocean           ",
                   "North Atlantic Temperate ", "Atlantic Tropical        ",
                   "South Atlantic Temperate ", "Southern Ocean           ",
                   "Indian Tropical          ", "South Indian Temperate   ")
      
    plot(1:24,as.numeric(df2),xaxt="n",ylab="sensitivity (ppm CO2) per flux region/month",
        xlab="Time",main=paste("Sensitivity of observation to a monthly flux from",gsub(" ","",transcom_names[i]),"at the given month"))
    grid()
    axis(side=1,at=1:24,labels=month_string)
  }
}


plot_base_pulse_flux = function(month=c(1),transcom_region=c(1)){
  
  old_opt = options()
  options(repr.plot.width=12, repr.plot.height=8)
  #-- Transcom region labels for regions 1:22
  transcom_names = c("North American Boreal    ", "North American Temperate ",
                     "South American Tropical  ", "South American Temperate ",
                     "Northern Africa          ", "Southern Africa          ",
                     "Eurasia Boreal           ", "Eurasia Temperate        ",
                     "Tropical Asia            ", "Australia                ",
                     "Europe                   ", "North Pacific Temperate  ",
                     "West Pacific Tropical    ", "East Pacific Tropical    ",
                     "South Pacific Temperate  ", "Northern Ocean           ",
                     "North Atlantic Temperate ", "Atlantic Tropical        ",
                     "South Atlantic Temperate ", "Southern Ocean           ",
                     "Indian Tropical          ", "South Indian Temperate   ")
  
  
  dts = c("2014-09","2014-10","2014-11","2014-12","2015-01","2015-02","2015-03","2015-04","2015-05","2015-06","2015-07","2015-08",
                            "2015-09","2015-10","2015-11","2015-12","2016-01","2016-02","2016-03","2016-04","2016-05","2016-06","2016-07","2016-08")
  require(ncdf4)
  
  for(i in 1:length(month)){
    for(j in 1:length(transcom_region)){
     f = nc_open(file.path(data_dir,"priors/sib4_NEE_pulses.nc4"))
     flx = ncvar_get(f,"NEE",start=c(1,1,transcom_region[j],month[i]),count=c(-1,-1,1,1))
     #-- I'm not positive but I've had to scale CO2 output by 12/44 and it might be needed here too
     #-- haven't identified where the issue is...
     flx = flx * 12/44
     nc_close(f)
     library(maps)
     w = map("world",plot = FALSE)
     grd = expand.grid(longitude=seq(-180,177.5,by=5),latitude=c(-89,seq(-86,86,by=4),89))
     grd$z = as.vector(flx)
     require(lattice)
     plot(levelplot(z ~ longitude + latitude,data=grd,
                 xlab="",ylab="",
                 #col.regions=ferret.palette("broad"),
                 col.regions=ferret_light_centered_palette_63,
                 at=seq(-max(abs(grd$z)),max(abs(grd$z)),length=63),
                 aspect="fill",useRaster=TRUE,
                 scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
                 main=paste("Avg NEE Flux (CO2), umolC/m2/sec for ",transcom_names[transcom_region[j]]," ",dts[month[i]],sep=""),
                 panel = function(x,y,z,...)
                 {
                   #llines(mm$x,mm$y,col="black",lty=1,lwd=3)
                   #llines(w$x,w$y,col="bl",lty=1,lwd=2)
                   panel.levelplot(x,y,z,...)
                   llines(w$x,w$y,col="black",lty=1,lwd=2)
                 }))
    }
  }
}


plot_flux_maps_annual_prior_post = function(prior_file_nc="/Users/aschuh/temp/ssim-ghg-output/gridded_fluxes_prior.nc4",
                                            posterior_file_nc="/Users/aschuh/temp/ssim-ghg-output/gridded_fluxes_posterior.nc4"){

fil_prior = nc_open(prior_file_nc)
prior_flux_samples = ncvar_get(fil_prior,"flux_samples")
prior_flux_mean_sampled = apply(prior_flux_samples,c(1,2,3),mean)
prior_flux_mean_analytical = ncvar_get(fil_prior,"flux_mean")
prior_flux_sd_sampled = apply(prior_flux_samples,c(1,2,3),sd)
prior_flux_annual_flux_mean_analytical = apply(prior_flux_mean_analytical,c(1,2),sum)*1/2
prior_flux_annual_flux_mean_sampled = apply(prior_flux_samples,c(1,2),sum)*1/2
prior_flux_annual_flux_sd_sampled = apply(apply(prior_flux_samples,c(1,2,4),sum)*1/2,c(1,2),sd)

rng_mn_prior = range(prior_flux_annual_flux_mean_analytical)
rng_sd_prior = range(prior_flux_annual_flux_sd_sampled)

fil_posterior = nc_open(posterior_file_nc)
posterior_flux_samples = ncvar_get(fil_posterior,"flux_samples")
posterior_flux_mean_sampled = apply(posterior_flux_samples,c(1,2,3),mean)
posterior_flux_mean_analytical = ncvar_get(fil_posterior,"flux_mean")
posterior_flux_sd_sampled = apply(posterior_flux_samples,c(1,2,3),sd)
posterior_flux_annual_flux_mean_analytical = apply(posterior_flux_mean_analytical,c(1,2),sum)*1/2
posterior_flux_annual_flux_mean_sampled = apply(posterior_flux_samples,c(1,2),sum)*1/2
posterior_flux_annual_flux_sd_sampled = apply(apply(posterior_flux_samples,c(1,2,4),sum)*1/2,c(1,2),sd)

prior_flux_annual_flux_reduction_error_sampled = 1 - (posterior_flux_annual_flux_sd_sampled/prior_flux_annual_flux_sd_sampled)
prior_flux_annual_flux_reduction_error_sampled[is.na(prior_flux_annual_flux_reduction_error_sampled)] = 0

rng_mn_posterior = range(posterior_flux_annual_flux_mean_analytical)
rng_sd_posterior = range(posterior_flux_annual_flux_sd_sampled)
rng_sd_reduction = range(prior_flux_annual_flux_reduction_error_sampled)

library(maps)
w = map("world",plot=FALSE)

grd = expand.grid(longitude=fil_prior$dim$longitude$vals,latitude=fil_prior$dim$latitude$vals)

#-- units need to go from kgCO2/m2/sec to gC/m2/yr, 
units_scaling = 1e3*12/44*3600*24*30.5
grd$prior.mean = as.vector(prior_flux_annual_flux_mean_analytical)*units_scaling
grd$posterior.mean = as.vector(posterior_flux_annual_flux_mean_analytical)*units_scaling
grd$prior.sd = as.vector(prior_flux_annual_flux_sd_sampled)*units_scaling
grd$posterior.sd = as.vector(posterior_flux_annual_flux_sd_sampled)*units_scaling
grd$reduction.sd = as.vector(prior_flux_annual_flux_reduction_error_sampled)

rng_mn_prior = rng_mn_prior * units_scaling
rng_sd_prior = rng_sd_prior * units_scaling
rng_mn_posterior = rng_mn_posterior * units_scaling
rng_sd_posterior = rng_sd_posterior * units_scaling
rng_sd_reduction = rng_sd_reduction 

plt1 = levelplot(prior.mean ~ longitude + latitude,data=grd,col.regions=my.col(50),
                 at=seq(-max(abs(c(rng_mn_prior,rng_mn_posterior))),max(abs(c(rng_mn_prior,rng_mn_posterior))),length=50),
                 main=c("Annual Prior Mean Flux (gC/m2/yr)"),xlab="",ylab="",aspect="iso",useRaster=TRUE,
                 panel = function(..., at, region,contour = FALSE, labels = NULL) {
                         panel.levelplot(..., at = at, contour = contour,labels = labels)
                         llines(w$x,w$y,col="black")}
)

plt2 = levelplot(posterior.mean ~ longitude + latitude,data=grd,col.regions=my.col(50),cuts=50,
                 at=seq(-max(abs(c(rng_mn_prior,rng_mn_posterior))),max(abs(c(rng_mn_prior,rng_mn_posterior))),length=50),
                 main=c("Annual Posterior Mean Flux (gC/m2/yr)"),xlab="",ylab="",aspect="iso",useRaster=TRUE,
                 panel = function(..., at, region,contour = FALSE, labels = NULL) {
                   panel.levelplot(..., at = at, contour = contour,labels = labels)
                   llines(w$x,w$y,col="black")})

plt3 = levelplot(prior.sd ~ longitude + latitude,data=grd,col.regions=my.col(100)[51:100],cuts=50,
                 at=seq(0,max(abs(c(rng_sd_prior,rng_sd_posterior))),length=50),
                 main=c("Annual Prior Flux Standard Deviation (gC/m2/yr)"),xlab="",ylab="",aspect="iso",useRaster=TRUE,
                 panel = function(..., at, region,contour = FALSE, labels = NULL) {
                   panel.levelplot(..., at = at, contour = contour,labels = labels)
                   llines(w$x,w$y,col="black")})

plt4 = levelplot(posterior.sd ~ longitude + latitude,data=grd,col.regions=my.col(100)[51:100],cuts=50,
                 at=seq(0,max(abs(c(rng_sd_prior,rng_sd_posterior))),length=50),
                 main=c("Annual Posterior Flux Standard Deviation (gC/m2/yr)"),xlab="",ylab="",aspect="iso",useRaster=TRUE,
                 panel = function(..., at, region,contour = FALSE, labels = NULL) {
                   panel.levelplot(..., at = at, contour = contour,labels = labels)
                   llines(w$x,w$y,col="black")})

plt5 = levelplot(reduction.sd ~ longitude + latitude,data=grd,col.regions=my.col(50),cuts=50,
                 at=seq(-max(abs(c(rng_sd_reduction))),max(abs(c(rng_sd_reduction))),length=50),
                 main=c("Annual Flux Reduction in Standard Deviation ( 1 - (posterior sd)/(prior sd) )"), xlab="",ylab="",aspect="iso",useRaster=TRUE,
                 panel = function(..., at, region,contour = FALSE, labels = NULL) {
                   panel.levelplot(..., at = at, contour = contour,labels = labels)
                   llines(w$x,w$y,col="black")})

options(jupyter.plot_scale=1)

print(plt1)
print(plt2)
print(plt3)
print(plt4)
print(plt5)

}

plot_flux_maps_annual_prior_post_truth=function (inv_object = ret2, true_state = state_vector_true, 
                                                 prior_mean_ncdf = file.path(data_dir, "priors/prior_SiB4.nc"), 
                                                 center_prior_on_zero = TRUE) 
{
  print("creating gridded fluxes....")
  con = nc_open(prior_mean_ncdf)
  longitude_prior = con$dim$longitude$vals
  latitude_prior = con$dim$latitude$vals
  NEE = ncvar_get(con, "NEE")
  nc_close(con)
  NEE_1x1 = aaply(NEE, 3, .fun = function(x) {
    expand_5x4_2_1x1(x)
  }) %>% aperm(c(2, 3, 1))
  
  NEE_transcom = aaply(NEE_1x1, 3, .fun = function(x) {
    grid2transcom(x)
  })
  
  tr_dir = file.path(data_dir, "/transcom/", sep = "")
  x_prior_matrix = matrix(inv_object$prior$x_hat, nrow = 24, 
                          byrow = FALSE)
  x_hat_matrix = matrix(inv_object$posterior$x_hat, nrow = 24, 
                        byrow = FALSE)
  true_state_matrix = matrix(true_state, nrow = 24, byrow = FALSE)
  
  #-- might need NEE_transcom[,2:23], not sure
  prior_flux_unc = diag(as.vector(NEE_transcom[,1:22])) %*% inv_object$prior$Sx %*%diag(as.vector(NEE_transcom[,1:22]))
  post_flux_unc = diag(as.vector(NEE_transcom[,1:22])) %*% inv_object$posterior$Sx  %*% diag(as.vector(NEE_transcom[,1:22]))
  
  A = cbind(diag(1,nrow=22),diag(1,nrow=22),diag(1,nrow=22),diag(1,nrow=22),diag(1,nrow=22),diag(1,nrow=22),
            diag(1,nrow=22),diag(1,nrow=22),diag(1,nrow=22),diag(1,nrow=22),diag(1,nrow=22),diag(1,nrow=22),
            diag(1,nrow=22),diag(1,nrow=22),diag(1,nrow=22),diag(1,nrow=22),diag(1,nrow=22),diag(1,nrow=22),
            diag(1,nrow=22),diag(1,nrow=22),diag(1,nrow=22),diag(1,nrow=22),diag(1,nrow=22),diag(1,nrow=22))
  
  annual_avg_prior_flux_cov = 0.5*A%*%prior_flux_unc%*%t(0.5*A) 
  annual_avg_post_flux_cov = 0.5*A%*%post_flux_unc%*%t(0.5*A)
  
  #-- these are plotted
  gridded_1x1_prior_sd_flux_annual = transcom2grid(sqrt(diag(annual_avg_prior_flux_cov))) * 1e-15
  gridded_1x1_post_sd_flux_annual =  transcom2grid(sqrt(diag(annual_avg_post_flux_cov)))* 1e-15
  gridded_1x1_post_sd_reduction_annual = 1 - (gridded_1x1_post_sd_flux_annual/gridded_1x1_prior_sd_flux_annual)
  ######
  
  if (!center_prior_on_zero) {
    x_prior_matrix = x_prior_matrix + 1
    x_hat_matrix = x_hat_matrix + 1
    true_state_matrix = true_state_matrix + 1
  }
  gridded_1x1_prior_state = aaply(x_prior_matrix, 1, .fun = function(x) {
    transcom2grid(x, model.grid.x = 1, model.grid.y = 1, 
                  file_location = data_dir)
  }) %>% aperm(c(2, 3, 1))
  gridded_1x1_posterior_state = aaply(x_hat_matrix, 1, .fun = function(x) {
    transcom2grid(x, model.grid.x = 1, model.grid.y = 1, 
                  file_location = data_dir)
  }) %>% aperm(c(2, 3, 1))
  gridded_1x1_true_state = aaply(true_state_matrix, 1, .fun = function(x) {
    transcom2grid(x, model.grid.x = 1, model.grid.y = 1, 
                  file_location = data_dir)
  }) %>% aperm(c(2, 3, 1))
  gridded_1x1_prior_mean_flux = NEE_1x1 * gridded_1x1_prior_state
  gridded_1x1_posterior_mean_flux = NEE_1x1 * gridded_1x1_posterior_state
  gridded_1x1_truth = NEE_1x1 * gridded_1x1_true_state
  
  #-- these are plotted
  gridded_1x1_prior_mean_flux_annual = apply(gridded_1x1_prior_mean_flux, 
                                             c(1, 2), sum)/2
  gridded_1x1_posterior_mean_flux_annual = apply(gridded_1x1_posterior_mean_flux, 
                                                 c(1, 2), sum)/2
  gridded_1x1_true_mean_flux_annual = apply(gridded_1x1_truth, 
                                            c(1, 2), sum)/2
  ######
  
  library(maps)
  w = map("world", plot = FALSE)
  grd = expand.grid(longitude = seq(-179.5, 179.5, by = 1), 
                    latitude = seq(-89.5, 89.5, by = 1))
  units_scaling = 1000 * 12/44 * 3600 * 24 * 30.5
  grd$prior.mean = as.vector(gridded_1x1_prior_mean_flux_annual) * 
    units_scaling
  grd$posterior.mean = as.vector(gridded_1x1_posterior_mean_flux_annual) * 
    units_scaling
  grd$truth = as.vector(gridded_1x1_true_mean_flux_annual) * 
    units_scaling
  grd$difference = grd$posterior.mean - grd$truth
  
  grd$prior_sd = as.vector(gridded_1x1_prior_sd_flux_annual)*units_scaling
  grd$post_sd = as.vector(gridded_1x1_post_sd_flux_annual)*units_scaling  
  grd$reduction_sd = as.vector(gridded_1x1_post_sd_reduction_annual)
  
  rng_mn = range(c(grd$prior.mean, grd$posterior.mean, grd$truth, 
                   grd$difference))

  rng_sd = range(c(grd$prior_sd, grd$post_sd))
  
  plt1 = levelplot(prior.mean ~ longitude + latitude, data = grd, 
                   col.regions = my.col(50), at = seq(-max(abs(c(rng_mn))), 
                                                      max(abs(c(rng_mn))), length = 50), main = c("Annual Prior Mean Flux (gC/m2/yr)"), 
                   xlab = "", ylab = "", aspect = "iso", useRaster = TRUE, 
                   panel = function(..., at, region, contour = FALSE, labels = NULL) {
                     panel.levelplot(..., at = at, contour = contour, 
                                     labels = labels)
                     llines(w$x, w$y, col = "black")
                   })
  plt2 = levelplot(posterior.mean ~ longitude + latitude, data = grd, 
                   col.regions = my.col(50), cuts = 50, at = seq(-max(abs(c(rng_mn))), 
                                                                 max(abs(c(rng_mn))), length = 50), main = c("Annual Posterior Mean Flux (gC/m2/yr)"), 
                   xlab = "", ylab = "", aspect = "iso", useRaster = TRUE, 
                   panel = function(..., at, region, contour = FALSE, labels = NULL) {
                     panel.levelplot(..., at = at, contour = contour, 
                                     labels = labels)
                     llines(w$x, w$y, col = "black")
                   })
  plt3 = levelplot(truth ~ longitude + latitude, data = grd, 
                   col.regions = my.col(50), cuts = 50, at = seq(-max(abs(c(rng_mn))), 
                                                                 max(abs(c(rng_mn))), length = 50), main = c("True Flux (gC/m2/yr)"), 
                   xlab = "", ylab = "", aspect = "iso", useRaster = TRUE, 
                   panel = function(..., at, region, contour = FALSE, labels = NULL) {
                     panel.levelplot(..., at = at, contour = contour, 
                                     labels = labels)
                     llines(w$x, w$y, col = "black")
                   })
  plt4 = levelplot(difference ~ longitude + latitude, data = grd, 
                   col.regions = my.col(50), cuts = 50, at = seq(-max(abs(c(rng_mn))), 
                                                                 max(abs(c(rng_mn))), length = 50), main = c("Annual Posterior Mean Flux - TRUTH (gC/m2/yr)"), 
                   xlab = "", ylab = "", aspect = "iso", useRaster = TRUE, 
                   panel = function(..., at, region, contour = FALSE, labels = NULL) {
                     panel.levelplot(..., at = at, contour = contour, 
                                     labels = labels)
                     llines(w$x, w$y, col = "black")
                   })
  plt5 = levelplot(prior_sd ~ longitude + latitude, data = grd, 
                   col.regions = my.col(50)[26:50], at = seq(0, 
                                                                 max(abs(c(rng_sd))), length = 25), main = c("Transcom Average* Annual Prior Mean Flux Standard Deviation (PgC/region/yr)"), 
                   xlab = "", ylab = "", aspect = "iso", useRaster = TRUE, 
                   panel = function(..., at, region, contour = FALSE, labels = NULL) {
                     panel.levelplot(..., at = at, contour = contour, 
                                     labels = labels)
                     llines(w$x, w$y, col = "black")
                   })
  plt6 = levelplot(post_sd ~ longitude + latitude, data = grd, 
                   col.regions = my.col(50)[26:50],  at = seq(0, 
                                                                 max(abs(c(rng_sd))), length = 25), main = c("Transcom Average* Annual Posterior Mean Flux Standard Deviation (PgC/region/yr)"), 
                   xlab = "", ylab = "", aspect = "iso", useRaster = TRUE, 
                   panel = function(..., at, region, contour = FALSE, labels = NULL) {
                     panel.levelplot(..., at = at, contour = contour, 
                                     labels = labels)
                     llines(w$x, w$y, col = "black")
                   })
  plt7 = levelplot(reduction_sd ~ longitude + latitude, data = grd, 
                   col.regions = rev(my.col(50)),  at = seq(-1,1, length = 50), main = c("Transcom Average* Flux Standard Deviation Reduction by Inversion"), 
                   xlab = "", ylab = "", aspect = "iso", useRaster = TRUE, 
                   panel = function(..., at, region, contour = FALSE, labels = NULL) {
                     panel.levelplot(..., at = at, contour = contour, 
                                     labels = labels)
                     llines(w$x, w$y, col = "black")
                   })
  
  options(jupyter.plot_scale = 1)
  print(plt1)
  print(plt2)
  print(plt3)
  print(plt4)
  print(plt5)
  print(plt6)
  print(plt7)
}



pad = function (x, width, fill = " ", left = TRUE) 
{
  xneg = FALSE
  nc = nchar(x)
  if (width > nc) {
    if (left) {
      str = paste(paste(rep(fill, width - nc), collapse = ""), 
                  x, sep = "")
    }
    else {
      str = paste(x, paste(rep(fill, width - nc), collapse = ""), 
                  sep = "")
    }
  }
  else {
    str = x
  }
  return(str)
}

ferret_light_centered_palette_63 = c("#00FFFFFF", "#03F4FFFF", "#05E9FFFF", "#08DEFFFF", "#0BD3FFFF",
                                     "#0EC8FFFF", "#10BDFFFF", "#13B2FFFF", "#16A7FFFF",
                                     "#199CFFFF", "#1B91FFFF", "#1E86FFFF", "#217BFFFF", 
                                     "#2470FFFF", "#2665FFFF", "#295AFFFF", "#2C50FFFF", "#2F45FFFF",
                                     "#313AFFFF", "#3A3AFFFF", "#4A4AFFFF", "#5A5AFFFF", "#6B6BFFFF",
                                     "#7B7BFFFF", "#8C8CFFFF", "#9C9CFFFF", "#ADADFFFF",
                                     "#BDBDFFFF", "#CECEFFFF", "#DEDEFFFF", "#EFEFFFFF", "#FFFFFFFF", 
                                     "#FFEFEFFF", "#FFDEDEFF", "#FFCECEFF", "#FFBDBDFF",
                                     "#FFADADFF", "#FF9C9CFF", "#FF8C8CFF", "#FF7B7BFF", "#FF6B6BFF",
                                     "#FF5A5AFF", "#FF4A4AFF", "#FF3A3AFF", "#FF3A31FF",
                                     "#FF452FFF", "#FF502CFF", "#FF5A29FF", "#FF6526FF", "#FF7024FF", 
                                     "#FF7B21FF", "#FF861EFF", "#FF911BFF", "#FF9C19FF",
                                     "#FFA716FF", "#FFB213FF", "#FFBD10FF", "#FFC80EFF", "#FFD30BFF",
                                     "#FFDE08FF", "#FFE905FF", "#FFF403FF", "#FFFF00FF")
