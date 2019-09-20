#Function to simulate species distribution and abundance with respect to spatiotemporal and environmental covariates
# Brodie et al.(2019).Ecography. Trade-offs in covariate selection for species distribution models: a methodological comparison. doi: 10.1111/ecog.04707

#Set working directory to the repository 
setwd('~/PatterOrProcess_Simulation/')

SimulateWorld <- function(){
  require(scales)
  require(raster)
  require(virtualspecies)
  
  #----Generate grid of locations through time----
  x_tot <- seq(1, 20, 1)
  y_tot <- seq(1, 20, 1)
  expand <- expand.grid(x = x_tot, y = y_tot)
  grid <- as.data.frame(cbind(x = rep(expand$x,20),y = rep(expand$y,20),year = rep(1:20,each=400)))
  years <- seq(1,20,1)
  
  #Generate Depth Covariate
  depth <- raster(ncol=20,nrow=20)
  ex <- extent(0.5,20.5,0.5,20.5)
  extent(depth) <- ex
  xy <- coordinates(depth)
  depth[] <- xy[,1]
  vals <- seq(-170,-80, 0.2255)
  depth[] <- rev(vals)
  depth <- flip(t(depth),1)
  depth <- raster::calc(depth, fun = function(x) jitter(x,amount=65)) #Add random noise
  
  #Generate Longitude Covariate
  lon <- raster(ncol=20,nrow=20)
  ex <- extent(0.5,20.5,0.5,20.5)
  extent(lon) <- ex
  xy <- coordinates(lon)
  lon[] <- xy[,1]
  vals <- seq(1,20)
  lon <-raster::calc(lon, fun = function(x) jitter(x,amount=1))
  lon[] <- rescale(lon@data@values, to = c(1,20))
  if (cor(lon@data@values,depth@data@values)>0.6) stop('lon~depth correlation is too high. Re-run.') 
  
  #Generate Latitude Covariate
  lat <- raster(ncol=20,nrow=20)
  ex <- extent(0.5,20.5,0.5,20.5)
  extent(lat) <- ex
  xy <- coordinates(lat)
  lat[] <- xy[,1]
  vals <- seq(1,20)
  vals <- vals[1:400]
  lat <- flip(t(lat),2)
  lat <-raster::calc(lat, fun = function(x) jitter(x,amount=1))
  lat[] <- rescale(lat@data@values, to = c(1,20))
  
  #Create Habitat Suitability for Latitude and Longitude
  latlon_stack <- stack( lat,lon)
  names(latlon_stack) <- c("lat","lon")
  latlon_parameters <- virtualspecies::formatFunctions(
    lat = c(fun="dnorm",mean=10,sd=20),
    lon = c(fun="dnorm",mean=10,sd=20))
  latlon_suitability <- virtualspecies::generateSpFromFun(latlon_stack,parameters=latlon_parameters, rescale = FALSE)
  # plot(latlon_suitability$suitab.raster) #Plot habitat suitability
  # virtualspecies::plotResponse(latlon_suitability) #Plot response curves
  
  for (y in years){
    print(paste0("Creating environmental simulation for Year ",y))
    
    #Generate Spatiotemporal Covariate
    high_years <- c(1,2,3,8,9,10,11) #Decide which years for each mode
    medium_years <- c(4,5,6,7,12,13,18,19,20)
    low_years <- c(14,15,16,17)
    if (y %in% high_years){
      spatiotemp <- raster("Shape_High.png")
    }
    if (y %in% medium_years){
      spatiotemp <- raster("Shape_Medium.png")
    }
    if (y %in% low_years){
      spatiotemp <- raster("Shape_Low.png")
    }
    ex <- extent(0.5,20.5,0.5,20.5)
    extent(spatiotemp) <- ex
    spatiotemp <- resample(spatiotemp,lat)
    spatiotemp <-raster::calc(spatiotemp, fun = function(x) jitter(x,amount=50))
    spatiotemp[] <- rescale(spatiotemp@data@values, to = c(0,1))
    
    #Generate Temperature Covariate
    temp_plain <- raster(ncol=20,nrow=20)
    ex <- extent(0.5,20.5,0.5,20.5)
    extent(temp_plain) <- ex
    xy <- coordinates(temp_plain)
    temp_plain[] <- xy[,2]
    min <- 0.1053*y + 1.8947 #temp increased by 2 degrees over 20 years
    max <- 0.1053*y + 5.8947
    vals <- seq(min,max,0.01)
    vals <- vals[1:400]
    temp_plain[] <- vals
    temp <- raster::calc (temp_plain, fun = function(x) jitter(x,amount=3.5))
    if (cor(temp@data@values,lat@data@values)<c(-0.6)) stop('lat~temp correlation >0.6')
    
    #Create Habitat Suitability for Temperature and Depth
    envir_stack <- stack(temp, depth)
    names(envir_stack) <- c('temp','depth')
    parameters <- virtualspecies::formatFunctions(depth = c(fun='logisticFun', beta=100, alpha=100),temp = c(fun="dnorm",mean=6,sd=12))
    envirosuitability <- virtualspecies::generateSpFromFun(envir_stack,parameters=parameters, rescale = FALSE)
    # plot(envirosuitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(envirosuitability) #plot response curves
    
    #----Combining all habitat suitability rasters----
    suitability_combo <- ((envirosuitability$suitab.raster*2)  + spatiotemp*1 + (latlon_suitability$suitab.raster*1)) / 4 
    suitability_combo[] <- rescale(suitability_combo@data@values, to=c(0,1))
    # plot(suitability_combo) #Plot final habitat suitability layer
    
    #----Convert suitability to Presence-Absence----
    suitability_PA <- convertToPA(suitability_combo, PA.method = "probability", beta = 0.5,
                                  alpha = -0.05, species.prevalence = NULL, plot = FALSE)
    
    #----Extract suitability for each survey location----
    print("Extracting suitability")
    for (i in 1:400){
      start_index <- min(which(grid$year==y))
      ii = (i + start_index) -1
      s <- raster::extract(suitability_combo,grid[ii,1:2]) 
      pa <- raster::extract(suitability_PA$pa.raster,grid[ii,1:2])
      t <- raster::extract(temp,grid[ii,1:2]) 
      d <- raster::extract(depth,grid[ii,1:2])
      h <- raster::extract(spatiotemp,grid[ii,1:2])
      grid$suitability[ii] <-  s
      grid$pres[ii] <-  pa
      grid$temp[ii] <-  t
      grid$depth[ii] <-  d
      grid$spatiotemp[ii] <-  h
    }
  }
  summary(grid)
  
  #----Create abundance as a function of the environment----
  grid$abundance <- ifelse(grid$pres==1,rlnorm(1,2,0.1)*grid$suitability,0)
  
  return(grid)
}


#----EXAMPLE----
# data <- SimulateWorld()
# #plot time-series of total catch
# plot(aggregate(abundance~year,data,FUN="sum"),type="l",ylab="Abundance")

