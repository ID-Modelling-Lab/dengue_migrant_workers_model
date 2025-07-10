## libraries
library(dust)
# library(mcstate)
library(odin.dust)
# library(readxl)
library(ggplot2)
# library(ggpubr)
library(dplyr)

## Human population parameters


# total number of CMP workers
N <- 1860000*0.23

# number of countries where the workers have been recruited
n_country <- 3

# list of countries
country_list <- c("India", "Bangladesh", "China") 

# Distribution of nationalities
n_country_percent <- c(0.50, 0.35, 0.15)

# Population size per country
n_country_population <- n_country_percent*N

# Termination rate
termination_rate <- 1/(42*365)

# Recruitment rate
recruitment_rate <- termination_rate


## Mosquito parameters

# Daily death rate
daily_death_rate <- 1/10

# Recruitment rate for mosquitoes
vector_recruitment <- daily_death_rate

# Extrinsic incubation period
ex_incubation_period <- 1/10


## Dengue Transmission Parameters

# Number of serotypes (DENV-1,2,3,4)
n_serotypes <- 4

# Bite rate
b <- 15

# Transmission probability (human to vector)
beta_h_to_v <- c(0.20, 0.20, 0.20, 0.20)

# Transmission probability (vector to human)
beta_v_to_h <- c(0.166, 0.176, 0.152, 0.144)


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
hospital_1 <- 0.30/4



## Sero-prevalence data

foi_estimate = data.frame(
  countries = c("India", "Bangladesh", "China"),
  mean = c(0.08278, 0.045, 0.018),
  lower = c(0.0785, 0.026, 0.002),
  upper = c(0.0872, 0.071, 0.040)
)

# Constructing Simple Catalytic Model

mid_age <- mean(18:60) + 0.5

seropositives <- function(mean_foi, age = mid_age){
  no_exp <- exp(-4*mean_foi*age)
  exactly_one_exp <- 4*(1-exp(-mean_foi*age))*exp(-3*mean_foi*age)
  more_than_one_exp <- 1 - no_exp - exactly_one_exp
  
  return(c(no_exposure = no_exp, exactly_one_exposure = exactly_one_exp, more_than_one_exposure = more_than_one_exp))
}


# Sample FOI from the interval uniformly or with other suitable distribution

#sampled_foi_india <- runif(10, min=0.0785, max=0.0872)
#sampled_foi_bangladesh <- runif(10, min=0.026, max=0.071)
#sampled_foi_china <- runif(10, min=0.002, max=0.040)

# Exposure matrix

exposure_matrix <- sapply(foi_estimate$mean/4, seropositives)
exposure_matrix <- t(exposure_matrix)
rownames(exposure_matrix) <- foi_estimate$countries


# Levels of exposure
exposure_0_array <- unname(exposure_matrix[, "no_exposure"])
exposure_exact1_array <- unname(exposure_matrix[, "exactly_one_exposure"])
exposure_1plus_array <- unname(exposure_matrix[, "more_than_one_exposure"])

## Initial conditions 

# Primary symptomatic infections
I0 <- array(0, dim = c(n_country, n_serotypes))
frac_serotype <- c(0.25, 0.25, 0.25, 0.25)

C0 = 0*t(frac_serotype*t(array(rep(n_country_population*exposure_exact1_array,n_serotypes), dim = c(n_country,n_serotypes) )))
S0 = 1*t(frac_serotype*t(array(rep(n_country_population*exposure_exact1_array,n_serotypes), dim = c(n_country,n_serotypes) )))

# Secondary symptomatic infections
I_ij0 <- array(0, dim=c(n_country,n_serotypes))
R0 <- n_country_population*exposure_1plus_array

# Susceptible population
S0_all =  n_country_population - rowSums(I0 + C0 + S0) - rowSums(I_ij0) - R0

#factor
factor <- 0.02095

# Susceptible vectors
S_v0 <- factor*N

# Exposed vectors
E_v0 <- c(0,0,0,0)

# Infected vectors
I_v0 <- c(1e-1,2e-2,1e-4,1e-5)


# Total vector population

#N_v <- S_v0 + sum(E_v0) + sum(I_v0)

par <- list(
  n_country = n_country,
  n_serotypes = n_serotypes,
  b = b,
  beta_h_to_v = beta_h_to_v,
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
  exposure_0_array = exposure_0_array,
  exposure_exact1_array = exposure_exact1_array,
  exposure_1plus_array = exposure_1plus_array,
  n_country_population = n_country_population
)


path_to_model <- here::here("C:/Users/keyas/OneDrive/Documents/internships/code/model_baseline_Keya.R")

dengue_model <- odin.dust::odin_dust(path_to_model)

# Instance of the model
model <- dengue_model$new(par, time = 1, n_particles = 1, ode_control = list(max_steps = 10000000, 
                                                                             step_size_min = 1e-14,
                                                                             debug_record_step_times = TRUE))
## no of year we want to run the simulation
run_year <- 500

## create array of time in day (required for solving ode)
t <- seq(1,(run_year*365), by = 1)

## Simulate the model
out <- model$simulate(t)
#plot(out[index_I[1],,])

start = Sys.time()

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

#+(n_country_population[]*exposure_1plus_array[]))

### Annual infection for four different serotypes (Combined and Per country)
target_annual_foi <- 0.1

foi <- array(out[index_foi,,],dim=c(n_serotypes,length(t)))

#function to calculate annual foi
calc_yearly<-function(mat){
  tapply(mat,(seq_along(mat)-1) %/%365,sum)
}

calc_yearly_wo_sum <- function(mat,time_frequency){
  mat[seq(time_frequency,length(mat),by=time_frequency)]
}

## calculate annual foi
foi_annual <- calc_yearly(colSums(foi))

## Plot annual foi
plot(foi_annual, type="l",ylim=c(0,0.5))+
  lines(rep(mean(tail(foi_annual, 20)),length(foi_annual)),type="l",col="red")+
  lines(rep(target_annual_foi,length(foi_annual)),type="l",col="blue")

#print(mean(tail(foi_annual,20)))

#### Using the total infection calculations (Lines 31-39 on baseline file)

## calculate the percentage of susceptible population getting infected annually

infection <- array(out[index_total_infection,,],dim=c(n_country,n_serotypes,length(t)))


## summing over serotype and get country wise infection
infection_country <- apply(infection,c(1,3),sum)

## now summing over country and make it annual infection
total_infection_annual <- calc_yearly(colSums(infection_country))

## annual susceptible
susceptible <- array(out[index_total_susceptible,,],dim=c(length(t)))

susceptible_annual <- calc_yearly_wo_sum(mat=susceptible,time_frequency=365)

## Plot the percentage of susceptible getting dengue infection annually 
### Error: Going above a 100% for avg FOI=0.1

plot(100*total_infection_annual/susceptible_annual,type="l",ylim=c(0,110)) +
  lines(rep(mean(100*total_infection_annual/susceptible_annual),length(100*total_infection_annual/susceptible_annual)), col="red")+
  lines(rep(100*target_annual_foi,length(foi_annual)),type="l",col="blue")


### Further checks on Infection

## serotype wise

# calculate serotype wise
infection_serotype <- apply(infection, c(2,3), sum)

# make it annual
infection_serotype_annual <- t(apply(infection_serotype, 1, calc_yearly))


line_color <- c("red", "blue", "green", "black")

tail_years <- 20

matplot(1:tail_years,
        t(infection_serotype_annual)[(dim(infection_serotype_annual)[2] - tail_years + 1):dim(infection_serotype_annual)[2],],
        type = "l",
        lty  = 1,
        lwd = 3,
        ylim = c(0,max(infection_serotype_annual)),
        col  = line_color,
        xlab = "Year",
        ylab = "Incidence of infection (serotype-wise)",
        xaxt = "n")

axis(1, at = 1:tail_years, labels = 1:tail_years)
legend("topright",
       legend = c("D1", "D2", "D3", "D4"),
       col    = line_color,
       lty    = 1,
       lwd = 3)



## country wise

# calculate country wise
infection_country <- apply(infection , c(1,3), sum)

# make it annual
infection_country_annual <- t(apply(infection_country, 1, calc_yearly))


line_color <- c("darkred", "darkblue", "darkgreen")

tail_years <- 20

matplot(1:tail_years,
        t(infection_country_annual)[(dim(infection_country_annual)[2] - tail_years + 1):dim(infection_country_annual)[2],],
        type = "l",
        lty  = 1,
        lwd = 3,
        ylim = c(0,max(infection_country_annual)),
        col  = line_color,
        xlab = "Year",
        ylab = "Incidence of infection (Nationality-wise)",
        xaxt = "n")

axis(1, at = 1:tail_years, labels = 1:tail_years)
legend("topright",
       legend = c("India", "Bangladesh", "China" ),
       col    = line_color,
       lty    = 1,
       lwd = 3)


### incidence per 100000 population for each nationalities 

## ERROR?: (shows that the number of infections are the same regardless of nationality)

pop_country <- array(out[index_pop_country,,], dim = c(n_country,length(t)))

pop_country_annual <- t(apply(pop_country,1, calc_yearly_wo_sum, time_frequency = 365))

infection_country_annual_p_100k <- 100000*(infection_country_annual/pop_country_annual)


matplot(
  1:tail_years,
  t(infection_country_annual_p_100k)[(dim(infection_country_annual_p_100k)[2] - tail_years + 1):dim(infection_country_annual_p_100k)[2],],
  type = "l",
  lty  = 1,
  lwd = 3,
  ylim = c(0, max(infection_country_annual_p_100k)),
  col = line_color,
  xlab = "Year",
  ylab = "Incidence of infection (Nationality-wise) per 100k population",
  xaxt = "n"
)


axis(1, at = 1:tail_years, labels = 1:tail_years)
legend("topright",
       legend = c("India", "Bangladesh", "China" ),
       col    = line_color,
       lty    = 1,
       lwd = 3)

## Seroprevalence

seropositive_country <- array(out[index_seropositive,,], dim = c(n_country,length(t)))
print(seropositive_country[1:3,1:10])

seropositive_total <- apply(seropositive_country, c(2), sum) 


seropositive_total_annual <- calc_yearly_wo_sum(seropositive_total, 365) #apply(seropositive_total,1, calc_yearly_wo_sum, time_frequency = 365)

seroprevalence <- 100*seropositive_total_annual/rowSums(pop_country_annual)

seroprevalence_country <- 100*apply(seropositive_country,1, calc_yearly_wo_sum, time_frequency = 365)/(pop_country_annual)

matplot((seroprevalence_country), type = "l")

plot(seroprevalence, type = "l")
