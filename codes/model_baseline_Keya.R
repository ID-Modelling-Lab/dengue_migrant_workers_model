
# total population
N <- sum(S_all[]) + sum(I[,]) + sum(C[,]) + sum(S[,]) + sum(I_ij[,]) + sum(R[])

# infected human population
## Loop over serotypes
infected_hum[] <- sum(I[,i]) + sum(I_ij[,i])
output(infected_hum[]) <- TRUE


# population each nationality
N_country[] <- S_all[i] + sum(I[i,]) + sum(C[i,]) + sum(S[i,]) + sum(I_ij[i,]) + R[i]
output(N_country[]) <- TRUE


# foi of infected vectors
foi_vectors[] <- b*beta_h_to_v[i]*I_v[i]/N
output(foi_vectors[]) <- TRUE

# foi of infected vectors (secondary)
foi_vectors_sec[] <- b*beta_h_to_v_sec[i]*I_v[i]/N
output(foi_vectors_sec[]) <- TRUE

# foi of infected humans
foi_hum[] <- b*beta_v_to_h[i]*infected_hum[i]/N
# output(foi_hum[]) <- TRUE

# vector recruitment (w/seasonality)
recruitment_rate_season <- vector_recruitment #vector_recruitment*(1+0.22*cos(2*3.141592653589793*t)/365 - 2.93)
output(recruitment_rate_season) <- TRUE


# secondary susceptible
S_sec[] <- sum(S[i,])

# Country specific symptomatic infection

# cases from primary infection
age_symp_inf_pri[] <- rho_1*sum(foi_vectors[])*S_all[i]
age_sero_symp_inf_sec[,] <- foi_vectors[j]*(S_sec[i]-S[i,j])

# cases from secondary infection
age_symp_inf_sec[] <- rho_2*sum(age_sero_symp_inf_sec[i,])

# Total infection
pri_inf[,] <- foi_vectors[j] * S_all[i]
sec_inf[,] <- foi_vectors[j] * (S_sec[i] - S[i,j])
total_infection[,] <- pri_inf[i,j] + sec_inf[i,j]
output(total_infection[,]) <- TRUE

# Total susceptible
total_susceptible[] <- S_all[i] + S_sec[i] 

output(total_susceptible[]) <- TRUE

# Total seropositive
seropositive[] <- sum(C[i,]) + sum(S[i,]) + R[i]
output(seropositive[]) <- TRUE



## Human Equations
#i = country, j = serotype
#Susceptible to all serotypes
deriv(S_all[1:n_country]) <- exposure_0_array[i]*recruitment_rate*N_country[i] - sum(foi_vectors[])*S_all[i] - termination_rate*S_all[i]

#Infected with serotype i
deriv(I[1:n_country,1:n_serotypes]) <- foi_vectors[j]*S_all[i] - gamma_1*I[i,j] - termination_rate*I[i,j] 

#Cross protected individuals
deriv(C[1:n_country, 1:n_serotypes]) <- gamma_1*I[i,j] - phi*C[i,j] - termination_rate*C[i,j]

#Susceptible for secondary infection
deriv(S[1:n_country, 1:n_serotypes]) <- exposure_exact1_array[i]*recruitment_rate*N_country[i]*frac_serotype[j] + phi*C[i,j] - (sum(foi_vectors_sec[])-foi_vectors_sec[j])*S[i,j] - termination_rate*S[i,j]

#Secondary infection
deriv(I_ij[1:n_country, 1:n_serotypes]) <- foi_vectors_sec[j]*(S_sec[i] - S[i,j]) - gamma_2*I_ij[i,j] - termination_rate*I_ij[i,j]

prod_term[] <- exposure_1plus_array[i]*N_country[i]
dim(prod_term) <- n_country 

#Recovery
deriv(R[1:n_country]) <- recruitment_rate*prod_term[i] + gamma_2*sum(I_ij[i,])  - termination_rate*R[i]


## Vector equations

deriv(S_v) <- recruitment_rate_season*(S_v + sum(E_v[]) + sum(I_v[])) - S_v*sum(foi_hum[]) - daily_death_rate*S_v

deriv(E_v[1:n_serotypes]) <- S_v*(foi_hum[i]) - ex_incubation_period*E_v[i] - daily_death_rate*E_v[i]

deriv(I_v[1:n_serotypes]) <- ex_incubation_period*E_v[i] - daily_death_rate*I_v[i]

## Parameters

## time dependent
recruitment_rate <- user()
termination_rate <-  user()
b <-  user()

# parameters

rho_1 <- user() 
rho_2 <- user()
ex_incubation_period <- user()

# disease progression related
phi <- user()
gamma_1 <- user()
gamma_2 <- user()

# demography related
vector_recruitment <- user()
daily_death_rate <- user()
#recruitment_rate_season <- user()

# age-group, serotype and transmission rates and age rate
n_country <- user() 
n_serotypes <- user() 
beta_h_to_v[] <- user()
beta_v_to_h[] <- user()
beta_h_to_v_sec[] <- user() ## new
# n_country_population[] <- user()
exposure_1plus_array[] <- user()
exposure_exact1_array[] <- user()
exposure_0_array[] <- user()
frac_serotype[] <- user()
dim(frac_serotype) <- n_serotypes


# initial conditions
S0_all[] <- user()
I0[,] <- user()
C0[,] <- user()
S0[,] <- user()
I_ij0[,] <- user()
R0[] <- user()
S_v0 <- user()
E_v0[] <- user()
I_v0[] <- user()


## parse initial conditions
initial(S_all[]) <- S0_all[i]
initial(I[,]) <- I0[i,j]
initial(C[,]) <- C0[i,j]
initial(S[,]) <- S0[i,j]
initial(I_ij[,]) <- I_ij0[i,j]
initial(R[]) <- R0[i]
initial(S_v) <- S_v0 
initial(E_v[]) <- E_v0[i]
initial(I_v[]) <- I_v0[i]


# DECLARE ALL THE DIMENSIONS OF VECTORS USED AS USER_DEFINED PARAMETERS

# initial conditions
dim(S0_all) <- n_country
dim(I0) <- c(n_country, n_serotypes)
dim(C0) <- c(n_country, n_serotypes)
dim(S0) <- c(n_country, n_serotypes)
dim(I_ij0) <- c(n_country, n_serotypes)
dim(R0) <- n_country
dim(E_v0) <- n_serotypes
dim(I_v0) <- n_serotypes


# state variables
dim(S_all) <- n_country
dim(I) <- c(n_country, n_serotypes)
dim(C) <- c(n_country, n_serotypes)
dim(S) <- c(n_country, n_serotypes)
dim(I_ij) <- c(n_country, n_serotypes)

# dim(H) <- n_age
dim(R) <- n_country
dim(E_v) <- n_serotypes
dim(I_v) <- n_serotypes
dim(N_country) <- n_country

# parameters
dim(beta_h_to_v) <- n_serotypes
dim(beta_h_to_v_sec) <- n_serotypes ## new
dim(beta_v_to_h) <- n_serotypes
dim(exposure_1plus_array) <- n_country
dim(exposure_exact1_array) <- n_country
dim(exposure_0_array) <- n_country
# dim(n_country_population) <- n_country
dim(seropositive) <- n_country

# extra 
dim(infected_hum) <- n_serotypes
dim(foi_vectors) <- n_serotypes
dim(foi_vectors_sec) <- n_serotypes
dim(foi_hum) <- n_serotypes

dim(pri_inf) <- c(n_country, n_serotypes)
dim(sec_inf) <- c(n_country, n_serotypes)
dim(total_infection) <- c(n_country, n_serotypes)

dim(total_susceptible) <- n_country

dim(age_symp_inf_pri) <- n_country
dim(age_symp_inf_sec) <- n_country
dim(age_sero_symp_inf_sec) <- c(n_country, n_serotypes)
dim(S_sec) <- n_country