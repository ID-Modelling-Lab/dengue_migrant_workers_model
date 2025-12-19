## libraries
library(dust)
# library(mcstate)
library(odin.dust)
# library(readxl)
library(ggplot2)
# library(ggpubr)
library(dplyr)
library(here)


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
#+ 0.000001


## Mosquito parameters

# Daily death rate
daily_death_rate <- 1/10

# Recruitment rate for mosquitoes
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

beta_h_to_v_sec <-c(0.20,0.20,0.20,0.20)

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
hospital_1 <- 0.30/4

## Sero-prevalence data

## FOI sampling
## could add for loop from here
## min and max interval for FOI
foi_india <- c(0.0785, 0.0872)
foi_bangladesh <- c(0.015, 0.017)
foi_china <- c(0.0003, 0.0164)

## number of samples generated
samples <- 250

## Getting samples for FOI
sample_foi_india <- round(runif(samples, min(foi_india), max(foi_india)),4)
sample_foi_bangladesh <- round(runif(samples, min(foi_bangladesh), max(foi_bangladesh)),4)
sample_foi_china <- round(runif(samples, min(foi_china), max(foi_china)),4)

## Combining samples by country into a data frame to be accessed later for seropostive/exposure matrix calculation
country_sample <- data.frame(
  India = sample_foi_india,
  Bangladesh = sample_foi_bangladesh,
  China = sample_foi_china
)

## assumption 2: unequal serotypes (how did you get these numbers?)
## If you assumed those, fine, we can do sensitivity analysis on this
## If equal serotype then, all these numbers will be 0.25 right?

frac_sero_india <- c(0.28, 0.26, 0.32, 0.14)
frac_sero_bangladesh <- c(0.11, 0.30, 0.57, 0.02)
frac_sero_china <- c(0, 0.88, 0.12, 0)

# frac_sero_mat <- t(array(c(0.28, 0.11, 0, 0.26, 0.30, 0.88, 0.32, 0.57, 0.12, 0.14, 0.02, 0), dim=c(3,4)))
## can write it in terms of variables you defined above, so no need type them again
frac_sero_mat <- t(array(c(frac_sero_india, frac_sero_bangladesh,frac_sero_china), dim = c(4,3)))


frac_sero_list <- list(
  India = frac_sero_india,
  Bangladesh = frac_sero_bangladesh,
  China = frac_sero_china
)

## mid_age calculation
mid_age <- mean(18:60) + 0.5 ## sensitivity analysis here

## unequal serotype distribution
seropositive_diff <- function(mean_foi, fractional_serotype, age = mid_age){
  
  lambda_i = 4*mean_foi*fractional_serotype
  lambda_tot = sum(lambda_i)
  
  no_exp <- exp(-lambda_tot*age)
  exactly_one_exp <- sum((1-exp(-lambda_i*age))*exp(-(lambda_tot - lambda_i)*age))
  more_than_one_exp <- 1- no_exp - exactly_one_exp
  
  c(no_exposure = no_exp, 
    exactly_one_exposure = exactly_one_exp,
    more_than_one_exposure = more_than_one_exp)
}


## looping over samples + exposure matrix
india_foi <- country_sample$India
bangladesh_foi <- country_sample$Bangladesh
china_foi <- country_sample$China


levels <- c("no_exposure", "exactly_one_exposure", "more_than_one_exposure")
countries <- colnames(country_sample)


exposure_mat <- array(
  0,
  dim=c(length(countries), length(levels), samples),
  dimnames = list(countries, levels, 1:samples)
)

for (s in 1:samples){
  for (ctry in countries){
    exposure_mat[ctry, ,s] <- seropositive_diff(
      mean_foi = country_sample[s, ctry],
      fractional_serotype = frac_sero_list[[ctry]]
    )
  }
}


## Exposure matrix
exposure_0_array <- unname(exposure_mat[,"no_exposure",])
exposure_exact1_array <- unname(exposure_mat[, "exactly_one_exposure",])
exposure_1plus_array <- unname(exposure_mat[, "more_than_one_exposure",])


print(exposure_0_array)

#factor
vector_factor <- 0.01617

#0.045 for FOI = 0.1
#0.0288 for FOI = 0.05
# 0.01617 for FOI = 0.01


# Susceptible vectors

S_v0 <- vector_factor * N

E_v0 <- c(0,0,0,0)            

I_v0 <- c(5, 1, 1e-5, 1e-6)

#function to calculate annual foi
calc_yearly<-function(mat){
  tapply(mat,(seq_along(mat)-1) %/%365,sum)
}

calc_yearly_wo_sum <- function(mat,time_frequency){
  mat[seq(time_frequency,length(mat),by=time_frequency)]
}

## set up foi average
foi_average <- vector("double", length = length(samples))

## common pars

common_par <- list(
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
  frac_serotype = frac_serotype,
  gamma_1 = gamma_1,
  gamma_2 = gamma_2
)

# path to model
path_to_model <- here::here("~/internships/code/model_baseline_Keya.R")
dengue_model <- odin.dust::odin_dust(path_to_model)

## main loop to set initial conditions + run code
for (s in seq_len(samples)){
  exp0_s <- exposure_0_array[, s]
  exp1_s <- exposure_exact1_array[, s]
  exp1plus_s <- exposure_1plus_array[, s]
  
  S0 <- t(frac_sero_mat* t(array(rep(n_country_population * exp1_s, n_serotypes),
                                  dim = c(n_country, n_serotypes))))

  C0   <- 0 * S0
  I0   <- array(0, dim = c(n_country, n_serotypes))
  I_ij0 <- array(0, dim = c(n_country, n_serotypes))
  R0   <- n_country_population * exp1plus_s
  S0_all <- n_country_population - rowSums(I0 + C0 + S0) - rowSums(I_ij0) - R0
  
  #print(exposure_mat)
  par <- c(common_par, list(
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
    exposure_0_array = exp0_s,
    exposure_exact1_array = exp1_s,
    exposure_1plus_array = exp1plus_s,
    n_country_population = n_country_population
  ))
  
  model <- dengue_model$new(par, time = 1, n_particles = 1, ode_control = list(max_steps = 10000000, 
                                                                               step_size_min = 1e-14,
                                                                               debug_record_step_times = TRUE))
  ## create array of time in day (required for solving ode)
  run_year <- 500
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
  index_pri_inf <- model$info()$index$pri_inf
  index_sec_inf <- model$info()$index$sec_inf
  
  
  target_annual_foi <- 0.05
  foi <- array(out[index_foi,,],dim=c(n_serotypes,length(t)))
  
  ## calculate annual foi
  foi_annual <- calc_yearly(colSums(foi))
  
  ## calculate mean foi
  foi_average[s] = mean(tail(foi_annual,20))
}

## Plot annual foi
plot(foi_average, type="l",ylim=c(max(foi_average)+0.05,min(foi_average)-0.05),xlab="Number of iterations", ylab="Annual FOI",main="Average FOI")+
  lines(rep(mean(foi_average),length(foi_average)),type="l",col="red")+
  lines(rep(0.1,length(foi_average)),type="l",col="blue")

print(mean(foi_average))

mean_foi <- mean(foi_average)
sd_foi <- sd(foi_average)
se <- sd_foi/(sqrt(length(foi_average)))


## Use quantile to calculate the 95% simulation interval
lower <- mean_foi - 1.96*se
upper <- mean_foi + 1.96*se


plot_foi_average <- data.frame(iterations = seq_along(foi_average),
                               foi = foi_average)

foi_plot <- ggplot(plot_foi_average, aes(x = iterations, y = foi)) +
  geom_line(color = "black") +
  geom_hline(yintercept = mean_foi, color = "red", linewidth = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2, fill = "darkgreen") +
  ggtitle("Iterations for average FOI") +
  xlab("Iterations") +
  ylab("FOI")

foi_plot + theme_bw()







#### Using the total infection calculations (Lines 31-39 on baseline file)

## calculate the percentage of susceptible population getting infected annually

infection <- array(out[index_total_infection,,],dim=c(n_country,n_serotypes,length(t)))
#print(infection[1:3,1:4,1:5])

## summing over serotype and get country wise infection
infection_country <- apply(infection,c(1,3),sum)
#dim(infection_country)

## now summing over country and make it annual infection
total_infection_annual <- calc_yearly(colSums(infection_country))
#plot(1:500,total_infection_annual,type="l",lwd=1.1,col="blue")

# ## annual susceptible
# susceptible <- array(out[index_total_susceptible,,],dim=c(n_country,length(t)))
# susceptible_annual <- calc_yearly_wo_sum(colSums(susceptible),time_frequency = 365)
# #plot(1:500,susceptible_annual,type="l",lwd=1.1,col="red")
# 
# ## Plot the percentage of susceptible getting dengue infection annually 
# 
# plot((total_infection_annual)/susceptible_annual*100,type="l",ylim=c(0,110),xlab="Time (Years)",ylab="Proportion of infected population",main="Proportion of infected population (FOI=0.1)") +
#   lines(rep(mean(100*total_infection_annual/susceptible_annual),length(100*total_infection_annual/susceptible_annual)), col="red")+
#   lines(rep(100*target_annual_foi,length(foi_annual)),type="l",col="blue")
# 
# print(mean(100*total_infection_annual/susceptible_annual))

## Plot on total infection

plot(total_infection_annual, type="l",xlab="Time",ylab="Infections per year",main="Total infection")

### Further checks on Infection

## serotype wise

# calculate serotype wise
infection_serotype <- apply(infection, c(2,3), sum)

# make it annual
infection_serotype_annual <- t(apply(infection_serotype, 1, calc_yearly))


line_color <- c("red", "blue", "green", "black")

tail_years <- 50

sero_wise_inf = t(infection_serotype_annual)[(dim(infection_serotype_annual)[2] - tail_years + 1):dim(infection_serotype_annual)[2],]

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
matplot(1:tail_years,
        sero_wise_inf,
        type = "l",
        lty  = 1,
        lwd = 1.5,
        ylim = c(0,max(sero_wise_inf)),
        col  = c("red", "blue", "green", "black"),
        xlab = "Year",
        ylab = "Incidence of infection",
        main = "Infection (serotype-wise)",
        xaxt = "n")


axis(1,at=1:tail_years,labels=1:tail_years)
legend("topright",
       inset=c(-0.25,0),
       legend = c("D1", "D2", "D3", "D4"),
       col    = c("red", "blue", "green", "black"),
       lty    = 1,
       lwd = 3,
       xpd=TRUE)

## Primary infection 
primary_infection <- array(out[index_pri_inf,,],dim=c(n_country,n_serotypes,length(t)))
primary_inf <- apply(primary_infection, 3, sum) 
annual_pri_inf <- calc_yearly(primary_inf)


## Primary infection (serotype-wise)

pri_inf_sero <- apply(primary_infection,c(2,3),sum)

annual_pri_inf_sero <- t(apply(pri_inf_sero,1,calc_yearly))


line_color <- c("red", "blue", "green", "black")

tail_years <- 50

plot(1:tail_years,
     annual_pri_inf[451:500],
     type="l",
     lty=1,
     lwd=1.5,
     ylim=c(0,5000),
     xlab = "Year",
     ylab = "Incidence of infection",
     main = "Primary Infection",
     xaxt = "n")

axis(1,at=1:tail_years,label=1:tail_years)

pri_sero_wise_inf = t(annual_pri_inf_sero)[(dim(annual_pri_inf_sero)[2] - tail_years + 1):dim(annual_pri_inf_sero)[2],]

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
matplot(1:tail_years,
        pri_sero_wise_inf,
        type = "l",
        lty  = 1,
        lwd = 1.5,
        col  = c("red", "blue", "green", "black"),
        xlab = "Year",
        ylim =c(0,1500),
        ylab = "Incidence of infection",
        main = "Primary Infection 
        (serotype-wise)",
        xaxt = "n")


axis(1,at=1:tail_years,labels=1:tail_years)
legend("topright",
       inset=c(-0.25,0),
       legend = c("D1", "D2", "D3", "D4"),
       col    = c("red", "blue", "green", "black"),
       lty    = 1,
       lwd = 3,
       xpd=TRUE)

## Total secondary infection
secondary_infection <- array(out[index_sec_inf,,],dim=c(n_country,n_serotypes,length(t)))
secondary_inf <- apply(secondary_infection, 3, sum) 
annual_sec_inf <- calc_yearly(secondary_inf)

plot(1:tail_years,
     annual_sec_inf[451:500],
     type="l",
     lty=1,
     lwd=1.5,
     ylim=c(0,5000),
     xlab = "Year",
     ylab = "Incidence of infection",
     main = "Secondary Infection",
     xaxt = "n")

axis(1,at=1:tail_years,labels=1:tail_years)

## Secondary infection (serotype-wise)


sec_inf_sero <- apply(secondary_infection,c(2,3),sum)

annual_sec_inf_sero <- t(apply(sec_inf_sero,1,calc_yearly))

line_color <- c("red", "blue", "green", "black")

tail_years <- 50

sec_sero_wise_inf = t(annual_sec_inf_sero)[(dim(annual_sec_inf_sero)[2] - tail_years + 1):dim(annual_sec_inf_sero)[2],]

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
matplot(1:tail_years,
        sec_sero_wise_inf,
        type = "l",
        lty  = 1,
        lwd = 1.5,
        ylim = c(0,3200),
        col  = c("red", "blue", "green", "black"),
        xlab = "Year",
        ylab = "Incidence of infection",
        main = "Secondary Infection 
        (serotype-wise)",
        xaxt = "n")
axis(1,at=1:tail_years,labels=1:tail_years)
legend("topright",
       inset=c(-0.25,0),
       legend = c("D1", "D2", "D3", "D4"),
       col    = c("red", "blue", "green", "black"),
       lty    = 1,
       lwd = 3,
       xpd=TRUE)


## country wise

# calculate country wise
infection_country <- apply(infection , c(1,3), sum)

# make it annual
infection_country_annual <- t(apply(infection_country, 1, calc_yearly))
dim(infection_country_annual)

line_color <- c("orange", "lightseagreen", "slateblue")

tail_years <- 50

country_wise_inf <-  t(infection_country_annual)[(dim(infection_country_annual)[2] - tail_years + 1):dim(infection_country_annual)[2],]

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
matplot(1:tail_years,
        country_wise_inf,
        type = "l",
        lty  = 1,
        lwd = 1.5,
        ylim = c(0,max(country_wise_inf)),
        col  = line_color,
        xlab = "Year",
        ylab = "Incidence of infection (Nationality-wise)",
        main = "Total infection (Nationality-wise)",
        xaxt = "n")

axis(1,at=1:tail_years,labels=1:tail_years)
legend("bottomright",
       legend = c("India", "Bangladesh", "China" ),
       col    = line_color,
       lty    = 1,
       lwd = 3,
       inset=c(-0.38,0))



### incidence per 100000 population for each nationalities 

pop_country <- array(out[index_pop_country,,], dim = c(n_country,length(t)))
dim(pop_country)

# plot(1:(50*365), pop_country[3,1:(50*365)], type = "l", lty = 1, lwd = 1.5,col="red")

pop_country_annual <- t(apply(pop_country,1, time_frequency = 365, calc_yearly_wo_sum))
# 
# plot(1:500, pop_country_annual[1,1:500], type="l", lty=1, lwd=1.5, col="green")

infection_country_annual_p_100k <- 100000*(infection_country_annual/pop_country_annual)
#print(infection_country_annual_p_100k[1:3,1:20])

tail_country_inf_per_100k <- t(infection_country_annual_p_100k)[(dim(infection_country_annual_p_100k)[2] - tail_years + 1):dim(infection_country_annual_p_100k)[2],1:3]

## This graph works, but hard to see the changes since they overlap significantly (Is there any way to edit the formatting? or should I do 3 different graphs?)

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
matplot(
  1:tail_years,
  tail_country_inf_per_100k  ,
  type = "l",
  lty  = 1,
  lwd = 2,
  ylim = c(0, max(tail_country_inf_per_100k)),
  xlim = c(0, tail_years),
  col = line_color,
  xlab = "Year",
  ylab = "Incidence of infection per 100k population",
  main = "Incidence of infection per 
  100k population (Nationality-wise)",
  xaxt = "n"
)

axis(1,at=1:tail_years,labels=1:tail_years)
legend("topright",
       legend = c("India", "Bangladesh", "China" ),
       col    = c("orange","lightseagreen","slateblue"),
       lty    = 1.5,
       lwd = 3,
       inset=c(-0.38,0))


## Percentage proportion of each country infection incidence

infection_country_proportion <- (infection_country_annual/pop_country_annual)*100

plot(1:tail_years,infection_country_proportion[1,451:run_year], type = "l", col="orange",lty="solid",lwd=1.5,ylim=c(0,5),xlab="Time (Years)",ylab="Percentage of infection", main="Percentage of total infection (nationality-wise)")
lines(1:tail_years,infection_country_proportion[2,451:run_year], col="lightseagreen",lty="solid",lwd=1.5)
lines(1:tail_years,infection_country_proportion[3,451:run_year], col="slateblue",lty="solid",lwd=1.5)

legend("topright",legend=c("India","Bangladesh","China"),col = c("orange","lightseagreen","slateblue"),lty=1.5,lwd=3)

## Seroprevalence

seropositive_country <- array(out[index_seropositive,,], dim = c(n_country,length(t)))
#print(seropositive_country[1:3,1:5])

seropositive_total <- apply(seropositive_country, 2, sum)

seropositive_total_annual <- calc_yearly_wo_sum(seropositive_total, 365) 

#apply(seropositive_total,1, calc_yearly_wo_sum, time_frequency = 365)

seroprevalence <- 100*seropositive_total_annual/colSums(pop_country_annual)

seroprevalence_country <- t(apply(seropositive_country,1, calc_yearly_wo_sum, time_frequency = 365))/pop_country_annual


plot(1:run_year, seroprevalence[1:run_year], type="l",lwd=1.1)+
  lines(rep(mean(seroprevalence),length(seroprevalence)),col="blue")


plot(1:run_year,seroprevalence_country[1,1:run_year],ylim=c(0.8,1.2),type="l",lwd=1.1,col="darkred")

lines(1:run_year,seroprevalence_country[2,1:run_year],type="l", lwd=1.1, col="darkblue")
lines(1:run_year,seroprevalence_country[3,1:run_year],type="l",lwd=1.1,col="darkgreen")

legend("bottomright",legend=c("India","Bangladesh", "China"),col=c("darkred","darkblue","darkgreen"),lwd=1)

