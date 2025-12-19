# remove everything from workspace
rm(list = ls(all = TRUE)) 

# collects garbage memory
gc(reset = TRUE) 

## libraries
library(odin.dust)

## Human population parameters
# total number of CMP workers
N <- 1860000*0.23 ## 23% of total foreign workforce source: ?

# number of countries where the workers have been recruited
n_country <- 3

# list of countries (only considered 3 countries with most contribution)
country_list <- c("India", "Bangladesh", "China") 

# Distribution of nationalities(source: COVID-19 survey by Clapham et al.)
n_country_percent <- c(0.50, 0.35, 0.15)

# Population size per country
n_country_population <- n_country_percent*N

# Termination rate (as the age range for employment is 18-60 years)
termination_rate <- 1/(42*365)

# Recruitment rate (assumed same to keep the total employee constant)
# We don't have any data on this
recruitment_rate <- termination_rate 

## Mosquito parameters

# Daily death rate
daily_death_rate <- 1/10

# Recruitment rate for mosquitoes 
# Assumed to be same as death rate to keep the vector population constant
vector_recruitment <- 1/10

# Extrinsic incubation period
ex_incubation_period <- 1/10


## Dengue Transmission Parameters

# Number of serotypes (DENV-1,2,3,4)
n_serotypes <- 4

# Bite rate
b <- 15

# Transmission probability (human to vector)
beta_h_to_v <- c(0.20, 0.20, 0.20, 0.20)

beta_h_to_v_sec <- c(0.20,0.20,0.20,0.20)

# Transmission probability (vector to human)
beta_v_to_h <- c(0.166, 0.176, 0.152, 0.148)


# Infectious period (primary)
gamma_1 <- 0.25

# Infectious period (secondary)
gamma_2 <- 0.25

# Primary infection symptom probability
rho_1 <- 0.2

# Secondary infection symptom probability
rho_2 <- rho_1*2

# Cross-protection duration
phi <- 1/365

# Hospitalization rate for secondary infection
hospital_2 <- 0.30

# Hospitalization rate for primary infection
hospital_1 <- hospital_2/4

## FOI sampling

## mid_age calculation
mid_age <- mean(18:60) + 0.5 ## sensitivity analysis here

## min and max interval for FOI
## these ranges are taken form literature

foi_estimate_range = data.frame(
  countries = c("India", "Bangladesh", "China"),
  lower = c(0.0785, 0.015, 0.0003),
  upper = c(0.0872, 0.017, 0.0164)
)


# Function to calculate level of exposure 
seropositives <- function(foi, age=mid_age){
  no_exp <- exp(-4*foi*age)
  exactly_one_exp <- 4*(1 - exp(-foi*age))*exp(-3*foi*age)
  more_than_one_exp <- 1 - no_exp - exactly_one_exp
  
  return(c(no_exposure = no_exp, 
           exactly_one_exposure = exactly_one_exp, 
           more_than_one_exposure = more_than_one_exp))
}


## we will now calculate exposure given foi and mean age

get_exposure_matrix <- function(foi_estimate_range){
  
  ## sample foi from the given ranges for all 3 countries
  foi_sample <- runif(
    nrow(foi_estimate_range),
    min = foi_estimate_range$lower,
    max = foi_estimate_range$upper
  )
  
## Now calculate exposure for each of the countries
  exposure_matrix <- sapply(foi_sample/4, seropositives)
  exposure_matrix <- t(exposure_matrix)
  
  return(exposure_matrix)
}

target_foi <- 0.01
#factor (the ration of human to mosquitoes)
vector_factor <- 0.01617

#0.045 for FOI = 0.1
#0.0288 for FOI = 0.05
# 0.01617 for FOI = 0.01


# Susceptible vectors

S_v0 <- vector_factor * N

E_v0 <- c(0,0,0,0)            

I_v0 <- c(5, 1, 1e-5, 1e-6)


##  number of years to simulate
run_year <- 500

### functions to calculate annual numbers
calc_yearly<-function(mat){
  tapply(mat,(seq_along(mat)-1) %/%365,sum)
}

calc_yearly_wo_sum <- function(mat,time_frequency){
  mat[seq(time_frequency,length(mat),by=time_frequency)]
}



out_model <- function(par){
  
  # path to model
  # path_to_model <- here::here("~/internships/code/model_baseline_Keya.R")
  path_to_model <- here::here("model_baseline_Keya.R")
  dengue_model <- odin.dust::odin_dust(path_to_model)
  
  # create an instance
  model <- dengue_model$new(par, time = 1, n_particles = 1, ode_control = list(max_steps = 10000000, 
                                                                               step_size_min = 1e-14,
                                                                               debug_record_step_times = TRUE))
  ## create array of time in day (required for solving ode)
  
  t <- seq(1,(run_year*365), by = 1)
  
  ## Simulate the model
  out <- model$simulate(t)
  
  # indices
  index_S_all <- model$info()$index$S_all
  index_I <- model$info()$index$I
  index_A <- model$info()$index$A
  index_C <- model$info()$index$C
  index_S <- model$info()$index$S
  index_Iij <- model$info()$index$I_ij
  index_R <- model$info()$index$R
  
  index_Sm <- model$info()$index$S_v
  index_Em <- model$info()$index$E_v
  index_Im <- model$info()$index$I_v
  index_foi <- model$info()$index$foi_vectors
  index_total_infection <- model$info()$index$total_infection
  index_total_susceptible <- model$info()$index$total_susceptible
  index_pop_country <- model$info()$index$N_country
  index_seropositive <- model$info()$index$seropositive
  index_pri_inf <- model$info()$index$pri_inf
  index_sec_inf <- model$info()$index$sec_inf
  
  
  ## Now calculate annual numbers according your need
  ## As an example I add only FOI 
  ## NOte: Can Add other variables 

  ## calculate annual foi
  foi_annual <- calc_yearly(array(out[index_foi,,],dim = c(n_serotypes,length(t))))

  ## remove output just to clear memory
  rm(list = c("out"))
  
  
  return(list(foi_annual = foi_annual))
  
  
}

n_sample <- 100
## create space for output for accumulating them for sampling

sample_foi_annual <- array(NA, dim = c(n_sample,n_serotypes, run_year))



get_sample_output <- function() {
  

  
## for reproducibility
set.seed(543793)
  
## iterated loops for sampling
for (i in 1:n_sample){
  

  
  # get exposure matrix sampling uniformly from the given range
  exposure_matrix <- get_exposure_matrix(
                                foi_estimate_range = foi_estimate_range)
  
  ## now calculate the subsequent quantities depends on exposure
  exposure_0_array <- unname(exposure_matrix[, "no_exposure"])
  exposure_exact1_array <- unname(exposure_matrix[, "exactly_one_exposure"])
  exposure_1plus_array <- unname(exposure_matrix[, "more_than_one_exposure"])
  
  
  ## now introduce exposure with different serotype in different countries
  ## IF we assume all are same then all these values will be replaced by 0.25 right?
  ## 
  frac_sero_india <- c(0.28, 0.26, 0.32, 0.14)
  frac_sero_bangladesh <- c(0.11, 0.30, 0.57, 0.02)
  frac_sero_china <- c(0, 0.88, 0.12, 0)
  
  ## combine all
  frac_sero_mat <- t(array(c(frac_sero_india, 
                             frac_sero_bangladesh,
                             frac_sero_china), dim = c(4,3)))
  
  # Primary symptomatic infections assumed to be zero
  I0 <- array(0, dim = c(n_country, n_serotypes))
  
  ## no individual with exposure with exactly one serotype, is in cross-protected compartments
  C0 = array(0, dim = c(n_country, n_serotypes))
  ## all the individuals exposed with exactly one serotype is in S0
  ## First calculate the number of people in each countries with exactly one exposure 
  # i.e., array(rep(n_country_population*exposure_exact1_array,n_serotypes), dim = c(n_country,n_serotypes) )
  ## and then split them up for different serotype by  multiplying with "frac_sero_mat"
  
  S0 = frac_sero_mat*array(rep(n_country_population*exposure_exact1_array,n_serotypes),
                           dim = c(n_country,n_serotypes) )

  # Secondary symptomatic infections
  I_ij0 <- array(0, dim = c(n_country,n_serotypes))
  R0 <- n_country_population*exposure_1plus_array
  
  # Susceptible population
  S0_all =  n_country_population - rowSums(I0 + C0 + S0) - rowSums(I_ij0) - R0
  
  
  par <- list(
    n_country = n_country,
    n_serotypes = n_serotypes,
    b = b,
    beta_h_to_v = beta_h_to_v,
    beta_h_to_v_sec = beta_h_to_v_sec,
    beta_v_to_h = beta_v_to_h,
    rho_1 = rho_1,
    rho_2 = rho_2,
    phi = phi,
    hospital_1 = hospital_1,
    hospital_2 = hospital_2,
    termination_rate = termination_rate,
    recruitment_rate = recruitment_rate,
    daily_death_rate = daily_death_rate,
    vector_recruitment = vector_recruitment,
    ex_incubation_period = ex_incubation_period,
    S0_all = S0_all,
    I0 = I0,
    C0 = C0,
    S0 = S0,
    I_ij0 = I_ij0,
    R0 = R0,
    S0_all = S0_all,
    S_v0 = S_v0,
    E_v0 = E_v0,
    I_v0 = I_v0,
    gamma_1 = gamma_1,
    gamma_2 = gamma_2,
    frac_sero_mat = frac_sero_mat,
    exposure_0_array = exposure_0_array,
    exposure_exact1_array = exposure_exact1_array,
    exposure_1plus_array = exposure_1plus_array

  )
  
  ## call the out_model function repeatedly and accumulate the outputs
  xx <- out_model(par)
  
  ## now stack them like this:
  sample_foi_annual[i,,] <- xx$foi_annual
  
}
  
  return(list(sample_foi_annual=sample_foi_annual
             ))
  
  
}


model_sample_output <- get_sample_output()


### Save the file in .rds and add some info of the model run so that
### you can recognize which model run it when you will analyze the data file for plotting etc.
### Create another folder called "model_output" (for example) to save the data file
### In the name for example I put diff_sero_exposure and target_foi so that I can recognize which datafile refers to which scenario
saveRDS(list(
             n_sample = n_sample,
             n_country = n_country,
             n_sero = n_serotypes,
             n_year = run_year,
             output = model_sample_output), 
        file = here::here("model_output", paste0("diff_sero_exposure","_target_foi_" , target_foi, ".rds" )))



