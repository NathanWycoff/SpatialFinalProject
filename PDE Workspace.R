## Author -- Nathan Wycoff
## Description -- An implementation of the finite difference method to 
## integrate the PDE's found in "On the Spatial Spread of the Grey Squirrel in Britain" by Okubo et al.
## in 2D.

require(distr)

##Problem Parameters
#Params in the PDE
a1 <- 0.82
b1 <- 10
D1 <- 17.9
a2 <- 0.62
b2 <- 0.75
D2 <- 17.9

#Get c_i from gamma_i
c1 <- 10
c2 <- 0.01

#Params for boundary conditions, Red squirrel begins at a point mass in the middle, grey begins uniformly
#distributed.
grey_density <- 1
red_density <- 0.39

time_interval <- c(0,8)
dist_interval_i <- c(0,24)# X interval is same in both directions.
dist_interval_j <- c(0,16)# X interval is same in both directions.

##Solver Parameters
delta_t <- 0.01
delta_x <- 1

##Initialization Variables
x_size_i <- (dist_interval_i[2] - dist_interval_i[1]) / delta_x#How many x (rows) do we have?
x_size_j <- (dist_interval_j[2] - dist_interval_j[1]) / delta_x#How many x (rows) do we have?
t_size <- (time_interval[2] - time_interval[1]) / delta_t#How many t (cols) do we have?
G <- array(0, dim=c(x_size_i, x_size_j, t_size))#Initialize a matrix of zeros for Grey Squirels
R <- array(red_density, dim=c(x_size_i, x_size_j, t_size))#Initialize a matrix of zeros for Red Squirels

##Enforce boundary and initial conditions
G[floor(x_size_i / 2), floor(x_size_j / 2),1] <- grey_density#Enforce Initial Conditions, point mass at 25

#Do iterations
for (t in 1:(t_size-1)) {
  for (i in 2:(x_size_i-1)) {
    for (j in 2:(x_size_j-1)) {
      #Grey Squirrel Update
      diffusion_term <- D1 * (G[i+1,j,t] - 2 * G[i,j,t] + G[i-1,j,t] + G[i,j+1,t] - 2 * G[i,j,t] + G[i,j-1,t])
      interaction_term <- delta_t * a1 * G[i,j,t] * (1 - b1 * G[i,j,t] - c1 * R[i,j,t])
      G[i,j,t+1] <- G[i,j,t] + (delta_t / delta_x^2) * diffusion_term + interaction_term
      
      #Red Squirrel Update
      #Grey Squirrel Update
      diffusion_term <- D2 * (R[i+1,j,t] - 2 * R[i,j,t] + R[i-1,j,t] + R[i,j+1,t] - 2 * R[i,j,t] + R[i,j-1,t])
      interaction_term <- delta_t * a2 * R[i,j,t] * (1 - b2 * R[i,j,t] - c2 * G[i,j,t])
      R[i,j,t+1] <- R[i,j,t] + (delta_t / delta_x^2) * diffusion_term + interaction_term
    }
  }
}

#Show result
image(G[,,t_size])
image(R[,,t_size])

