require(NicheMapR)
require(tidyr)
require(dplyr)
require(raster)
require(Rcpp)
library(ncdf4) 
library(rgdal)
require(GeoLight)
require(ggplot2)

############## Model parameters and functions ###########

require(Rcpp)
dir <- "" # set working directory

# 1. Load model functions
load(paste0(dir,"/Te_model.RData")) # Operative temperature model function
load(paste0(dir,"/Tb_t.RData")) # Transitent body temperature model
sourceCpp(paste0(dir,"/Cpp_functions.cpp")) # Probabilistic behavioral thermoregulation model (in Cpp)

# 2. Load the Leaf Area Index data and geographic coordinates of the tropical forest biome
load(paste0(dir,"/LAI/LAI_matrix_july_tropical_forest_HD.RData"))
load(paste0(dir,"/xy.values_HD.RData"))

# 3. Load the function compute microclimates (requires NicheMapR; Kearney and Porter 2017 Ecography) 
save(paste0(dir, "/compute_microclimates.RData"))

############## Run NicheMapR with simulated increments in air temperature ###########

M = 5.7         # Body mass (g)
a = 0.88        # skin absorbance
height = 0.01   # height from ground (cm)
month = 7       # month for simulations

Tpref_mean = 27  # prefered body temperature (ºC)
Tpref_sd = 1    # body temperature standard deviation
max_iter = 10   # number of iterations

deltaT <- seq(0, 4, length.out = 10) # simulated temperature increase (between 0 and 4 ºC)

meanTb_matrix_df <- SDTb_matrix_df <- Te_shade_matrix <- Te_sun_matrix <- array(NA, dim=c(1440, length(deltaT)))
meanTb_matrix_list <- SDTb_matrix_list <- Te_shade_matrix_list <- Te_sun_matrix_list <- list()
for(i in 1:nrow(xy.values)){
  t0 <- Sys.time()
  LAI <- LAI_matrix_july[i]
  
  for(j in 1:length(deltaT)){
    ## Run NicheMapR at each location (grid cell)
    microclimate <- compute_microclimates(run.gads = 0, lat = xy.values[i,2], long = xy.values[i,1], month = month, height = height, 
                                          LAI = LAI, delta_T = deltaT[j])
    
    microclimate_sun <- microclimate$microclimate_sun_min
    microclimate_shade <- microclimate$microclimate_shade_min
    
    S_sun = microclimate_sun$S
    Ta_sun = microclimate_sun$Ta
    Tsk_sun = microclimate_sun$Tsk
    Tg_sun = microclimate_sun$Tg
    V_sun = microclimate_sun$V
    
    S_shade = microclimate_shade$S
    Ta_shade = microclimate_shade$Ta
    Tsk_shade = microclimate_shade$Tsk
    Tg_shade = microclimate_shade$Tg
    V_shade = microclimate_shade$V
    
    Te_sun <- Te_model(M=M, a=a, ground=T, S=S_sun, Tsk=Tsk_sun, Ta=Ta_sun, Tg=Tg_sun, v=V_sun)
    Te_shade <- Te_model(M=M, a=a, ground=T, S=S_shade, Tsk=Tsk_shade, Ta=Ta_shade, Tg=Tg_shade, v=V_shade)

    ## Run thermoregulation model
    output <- behav_therm(M=M, a=a, ground=1, S_sun=S_sun, Ta_sun=Ta_sun, Tsk_sun=Tsk_sun, Tg_sun=Tg_sun, V_sun=V_sun,
                          S_shade=S_shade, Ta_shade=Ta_shade, Tsk_shade=Tsk_shade, Tg_shade=Tg_shade, V_shade=V_shade,
                          Tpref_mean=Tpref_mean, Tpref_sd=Tpref_sd, sigma_sun=1, sigma_shade=1, max_iter=max_iter)

    meanTb_matrix_df[,j] <- output[,1]
    SDTb_matrix_df[,j] <- output[,3]
    Te_shade_matrix[,j] <- Te_shade
    Te_sun_matrix[,j] <- Te_sun
  }
  
  meanTb_matrix_list[[i]] <- meanTb_matrix_df
  SDTb_matrix_list[[i]] <- SDTb_matrix_df
  Te_shade_matrix_list[[i]] <- Te_shade_matrix
  Te_sun_matrix_list[[i]] <- Te_sun_matrix
  
  t1 <- Sys.time()
  elapsed <- t1 - t0
  remaining <- (nrow(xy.values) - i) * elapsed

  print(paste(round(100 * i / nrow(xy.values),2), "%", "remaining:", round(remaining / 60 , 2), "min"))
}

############## Run NicheMapR across shade levels ###########

M = 5.7           # Body mass
a = 0.88         # skin absorbance
height = 0.01
month = 7

Tpref_mean = 27
Tpref_sd = 1
max_iter = 10

shade_reduction <- seq(0.5, 1, length.out = 10)

meanTb_matrix_df <- SDTb_matrix_df <- Te_shade_matrix <- array(NA, dim=c(1440, length(deltaT)))
meanTb_matrix_list <- SDTb_matrix_list <- Te_shade_matrix_list <- list()
for(i in 1:nrow(xy.values)){
  t0 <- Sys.time()
  LAI <- LAI_matrix_july[i] 

  for(j in 1:length(shade_reduction)){
    ## Run NicheMapR
    
    LAI_reduction <- LAI * shade_reduction[j]
  
    microclimate <- compute_microclimates(run.gads = 0, lat = xy.values[i,2], long = xy.values[i,1], month = month, height = height, 
                                          LAI = LAI_reduction, delta_T = 0)
    
    microclimate_sun <- microclimate$microclimate_sun_min
    microclimate_shade <- microclimate$microclimate_shade_min
    
    S_sun = microclimate_sun$S
    Ta_sun = microclimate_sun$Ta
    Tsk_sun = microclimate_sun$Tsk
    Tg_sun = microclimate_sun$Tg
    V_sun = microclimate_sun$V
    
    S_shade = microclimate_shade$S
    Ta_shade = microclimate_shade$Ta
    Tsk_shade = microclimate_shade$Tsk
    Tg_shade = microclimate_shade$Tg
    V_shade = microclimate_shade$V
    
    Te_sun <- Te_model(M=M, a=a, ground=T, S=S_sun, Tsk=Tsk_sun, Ta=Ta_sun, Tg=Tg_sun, v=V_sun)
    Te_shade <- Te_model(M=M, a=a, ground=T, S=S_shade, Tsk=Tsk_shade, Ta=Ta_shade, Tg=Tg_shade, v=V_shade)
    
    ## Run thermoregulation model
    output <- behav_therm(M=M, a=a, ground=1, S_sun=S_sun, Ta_sun=Ta_sun, Tsk_sun=Tsk_sun, Tg_sun=Tg_sun, V_sun=V_sun,
                          S_shade=S_shade, Ta_shade=Ta_shade, Tsk_shade=Tsk_shade, Tg_shade=Tg_shade, V_shade=V_shade,
                          Tpref_mean=Tpref_mean, Tpref_sd=Tpref_sd, sigma_sun=1, sigma_shade=1, max_iter=max_iter)
    
    meanTb_matrix_df[,j] <- output[,1]
    SDTb_matrix_df[,j] <- output[,3]
    Te_shade_matrix[,j] <- Te_shade
  }
  
  meanTb_matrix_list[[i]] <- meanTb_matrix_df
  SDTb_matrix_list[[i]] <- SDTb_matrix_df
  Te_shade_matrix_list[[i]] <- Te_shade_matrix
  
  t1 <- Sys.time()
  elapsed <- t1 - t0
  remaining <- (nrow(xy.values) - i) * elapsed
  
  print(paste(round(100 * i / nrow(xy.values),2), "%", "remaining:", round(remaining / 60 , 2), "min"))
}

############## Run NicheMapR across shade availabilities ###########

M = 5.7          # Body mass
a = 0.88         # skin absorbance
height = 0.01
month = 7

Tpref_mean = 27
Tpref_sd = 1
max_iter = 10

shade_availability <- seq(0, 1, length.out = 30)

meanTb_matrix_df <- SDTb_matrix_df <- Te_shade_matrix <- array(NA, dim=c(1440, length(shade_availability)))
meanTb_matrix_list <- SDTb_matrix_list <- Te_shade_matrix_list <- list()
for(i in 1:nrow(xy.values)){
  t0 <- Sys.time()
  
  ## Run NicheMapR
  LAI <- LAI_matrix_july[i]
  microclimate <- compute_microclimates(run.gads = 0, lat = xy.values[i,2], long = xy.values[i,1], month = month, height = height, 
                                        LAI = LAI, delta_T = 0)
  
  microclimate_sun <- microclimate$microclimate_sun_min
  microclimate_shade <- microclimate$microclimate_shade_min
  
  S_sun = microclimate_sun$S
  Ta_sun = microclimate_sun$Ta
  Tsk_sun = microclimate_sun$Tsk
  Tg_sun = microclimate_sun$Tg
  V_sun = microclimate_sun$V
  
  S_shade = microclimate_shade$S
  Ta_shade = microclimate_shade$Ta
  Tsk_shade = microclimate_shade$Tsk
  Tg_shade = microclimate_shade$Tg
  V_shade = microclimate_shade$V
  
  Te_sun <- Te_model(M=M, a=a, ground=T, S=S_sun, Tsk=Tsk_sun, Ta=Ta_sun, Tg=Tg_sun, v=V_sun)
  Te_shade <- Te_model(M=M, a=a, ground=T, S=S_shade, Tsk=Tsk_shade, Ta=Ta_shade, Tg=Tg_shade, v=V_shade)
  
  for(j in 1:length(shade_availability)){
    sigma_shade = shade_availability[j]
    
    ## Run thermoregulation model
    output <- behav_therm(M=M, a=a, ground=1, S_sun=S_sun, Ta_sun=Ta_sun, Tsk_sun=Tsk_sun, Tg_sun=Tg_sun, V_sun=V_sun,
                          S_shade=S_shade, Ta_shade=Ta_shade, Tsk_shade=Tsk_shade, Tg_shade=Tg_shade, V_shade=V_shade,
                          Tpref_mean=Tpref_mean, Tpref_sd=Tpref_sd, sigma_sun=1, sigma_shade=sigma_shade, max_iter=max_iter)
    
    
    meanTb_matrix_df[,j] <- output[,1]
    SDTb_matrix_df[,j] <- output[,3]
    Te_shade_matrix[,j] <- Te_shade
  }
  
  meanTb_matrix_list[[i]] <- meanTb_matrix_df
  SDTb_matrix_list[[i]] <- SDTb_matrix_df
  Te_shade_matrix_list[[i]] <- Te_shade_matrix
  
  t1 <- Sys.time()
  elapsed <- t1 - t0
  remaining <- (nrow(xy.values) - i) * elapsed
  
  print(paste(round(100 * i / nrow(xy.values),2), "%", "remaining:", round(remaining / 60 , 2), "min"))
}



























