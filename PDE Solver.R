##Author -- Nathan Wycoff
##Description -- An implementation of the finite difference method to 
##integrate the PDE's found in On the Spatial Spread of the Grey Squirrel in Britain by Okubo et al.

##Problem Parameters
time_interval <- c(0,10)
dist_interval <- c(0,50)

#Solver Parameters
delta_t <- 0.001
delta_x <- 1

#Initialization Variables
x_size <- (dist_interval[2] - dist_interval[1]) / delta_x#How many x (rows) do we have?
t_size <- (time_interval[2] - time_interval[1]) / delta_t#How many t (cols) do we have?
S <- matrix(rep(0, x_size * t_size), nrow = x_size)#Initialize a matrix of zeros

#Enforce boundary and initial conditions
S[2:49,1] <- 1#Enforce Initial Conditions, point mass at 25
S[1,1] <- 0# Boundary Condition on one end
S[50,1] <- 0# Boundary Condition on the other

#Do iterations
for (j in 1:(t_size-1)) {
  for (i in 2:(x_size-1)) {
    S[i,j+1] <- S[i,j] + (delta_t / delta_x^2) * (S[i+1,j] - 2 * S[i,j] + S[i-1,j])
  }
}

image(S)
