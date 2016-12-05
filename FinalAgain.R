##Author -- Nathan Wycoff
##Description -- Attempt to use PDE as prior information

#Requires spBayes
require(spBayes)

##Problem Parameters
time_interval <- c(0,10)
dist_interval <- c(0,50)
D <- 10#Diffusion Coefficient

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
    S[i,j+1] <- S[i,j] + D*(delta_t / delta_x^2) * (S[i+1,j] - 2 * S[i,j] + S[i-1,j])
  }
}

image(S)

######
# Now we evaluate several GP based methods on noise corrupted data
######
#Params
sigma <- 0.1#Variance of noise corruption
iters <- 100

results <- c()

#for (ayy in 1:iters) {
#Get some train and test data.
n_test <- 20
n_train <- 20
n_dynamic <- 200

train <- unique(data.frame(t(sapply(1:n_train, function(i) {
  x <<- sample(1:x_size,1)
  t <<- sample(1:t_size,1)
  return(c(x * delta_x + min(dist_interval),t * delta_t + min(dist_interval), S[x,t], S[x, t] + rnorm(1,0,sigma)))
}))))
colnames(train) <- c('x','t','Truth','Observed')

dynamic <- unique(data.frame(t(sapply(1:n_dynamic, function(i) {
  x <<- sample(1:x_size,1)
  t <<- sample(1:t_size,1)
  return(c(x * delta_x + min(dist_interval),t * delta_t + min(dist_interval), S[x,t], S[x, t] + rnorm(1,0,sigma)))
}))))
colnames(dynamic) <- c('x','t','Truth','Observed')

test <- unique(data.frame(t(sapply(1:n_test, function(i) {
  x <- sample(1:x_size,1)
  t <- sample(1:t_size,1)
  return(c(x * delta_x + min(dist_interval),t * delta_t + min(dist_interval), S[x,t], S[x, t] + rnorm(1,0,sigma)))
}))))
colnames(test) <- c('x','t','Truth','Observed')

train_coords <- as.matrix(train[,1:2])
dynamic_coords <- as.matrix(dynamic[,1:2])
test_coords <- as.matrix(test[,1:2])

#Baselin Mean Absolute error is the null model
baseline_mae <- mean(abs(test$Observed - mean(train$Observed)))

#Predictions from the Dynamical System
dynamical_mae <- mean(abs(test$Observed - test$Truth))

####
# Regular GP, fit using spLM
####
##Fit the model
#Prepare for the model
starting <- list(sigma.sq = 0.1, tau.sq = 0.1, phi = 3)
tuning <- list(sigma.sq = 1, tau.sq = 1, phi = 1)
priors <- list(sigma.sq.ig = c(0.1,0.1), tau.sq.ig = c(0.1,0.1), phi.unif = c(0.5,3.5))

#Actual model fitting
standard_fit <- spLM(Observed ~ 1, data = test, coords = train_coords, starting = starting, priors = priors,
                     tuning = tuning, cov.model= 'gaussian', n.samples = 1000)

##Evaluate the Model
#Get the posterior mean on each test point
standard_pred <- spPredict(standard_fit, test_coords, pred.covars = as.matrix(rep(1,dim(test_coords)[1])))
standard_pred <- rowMeans((standard_pred$p.y.predictive.samples))

#Calculate mean absolute error
standard_mae <- mean(abs(standard_pred - test$Observed))

####
# GP fit using spLM, but with priors derived from the differential equation.
####

###
#STEP 1: Fit a GP on the PDE solution, using vague priors with no nugget effect
###
##Fit the model
#Prepare for the model
starting <- list(sigma.sq = 1, phi = 1)
tuning <- list(sigma.sq = 1, phi = 1)
priors <- list(sigma.sq.ig = c(0.1,0.1),  phi.unif = c(0.5,3.5)) 

#Actual Model Fitting
precondition_fit <- spLM(Truth ~ 1, data = dynamic, coords = dynamic_coords, starting = starting, priors = priors,
                     tuning = tuning, cov.model= 'gaussian', n.samples = 1000)

##Moment Match our posterior on sigma square to an IG
y <- precondition_fit$p.theta.samples[,'sigma.sq']
alpha <- 2 + mean(y)^2 / var(y)
beta <- mean(y) * (1 + mean(y)^2/var(y))

print(c(alpha, beta))

###
#STEP 2: Use the new prior on sig.sq
###
starting <- list(sigma.sq = 0.1, tau.sq = 0.1, phi = 3)
tuning <- list(sigma.sq = 1, tau.sq = 1, phi = 1)
priors <- list(sigma.sq.ig = c(alpha, beta), tau.sq.ig = c(0.1,0.1), phi.unif = c(0.5,3.5))

conditioned_fit <- spLM(Truth ~ 1, data = train, coords = train_coords, starting = starting, priors = priors,
                         tuning = tuning, cov.model= 'gaussian', n.samples = 1000)

##Evaluate the Model
#Get the posterior mean on each test point
conditioned_pred <- spPredict(conditioned_fit, test_coords, pred.covars = as.matrix(rep(1,dim(test_coords)[1])))
conditioned_pred <- rowMeans((conditioned_pred$p.y.predictive.samples))

#Calculate mean absolute error
conditioned_mae <- mean(abs(conditioned_pred - test$Observed))

print(c(baseline_mae, dynamical_mae, standard_mae, conditioned_mae))

#results <- c(results, baseline_mae, dynamical_mae, standard_mae, conditioned_mae)
#}

#sapply(1:4, function(x) mean(results[seq(x, length(results), 4)]))

