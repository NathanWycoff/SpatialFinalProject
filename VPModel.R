##Author -- Nathan Wycoff
##Description -- Comparison of methods to estimate a noise corrupted model with perfect knowledge
## of the true model. Attempt to incorperate dynamical system information in a GP prior.

##TODO:
#  -- I don't think we currently have a burn in
#  -- Do a model spBayes(Observed ~ Truth,...)
#  -- Do a similar model, but with informed priors similar to the preconditioned model.
#  -- Do the nonparameteric kernel function model

#Requires spBayes
require(spBayes)

##Problem Parameters
time_interval <- c(0,10)
dist_interval <- c(0,50)
D_true <- 10#Actual Diffusion Coefficient
D_used <- 30#Diffusion coefficient used in modeling

#Solver Parameters
delta_t <- 0.001
delta_x <- 1


#######
## Integrate THE TRUE MODEL
#######
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
    S[i,j+1] <- S[i,j] + D_true*(delta_t / delta_x^2) * (S[i+1,j] - 2 * S[i,j] + S[i-1,j])
  }
}

S_true <- S


#######
## Integrate THE USED MODEL
#######
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
    S[i,j+1] <- S[i,j] + D_used*(delta_t / delta_x^2) * (S[i+1,j] - 2 * S[i,j] + S[i-1,j])
  }
}

S_used <- S

par(mfrow=c(2,1))
image(S_used, 'main' = 'USED')
image(S_true, main = 'TRUE')
par(mfrow=c(1,1))

######
# Now we evaluate several GP based methods on noise corrupted data
######
#Params
sigma <- 0.1#Variance of noise corruption
iters <- 100

#Get some train and test data.
n_samples <- 100
n_train <- 50
results <- c()

for (ayy in 1:iters) {
  x_samples <- c()
  t_samples <- c()
  while (length(x_samples) < n_samples) {
    if (n_samples > dim(S)[1] * dim(S)[2]) {
      print("Impossible sample without replacement")
      break
    }
    x <- sample(1:x_size,1)
    t <- sample(1:t_size,1)
    
    #Check if we already have that one
    if (x %in% x_samples) {
      taken_ts <- t_samples[which(x_samples==x)]
      if (t %in% taken_ts){
        next
      }
    }
    
    #Add it to our collection
    x_samples <- c(x_samples, x)
    t_samples <- c(t_samples, t)
  }
  sample_df <- data.frame(x_samples * delta_x + min(dist_interval),
                          t_samples * delta_t + min(dist_interval),
                          sapply(1:n_samples, function(i) S_used[x_samples[i],t_samples[i]]),
                          sapply(1:n_samples, function(i) S_true[x_samples[i],t_samples[i]] + rnorm(1,0,sigma)))
  colnames(sample_df) <- c('x','t','Estimated','Observed')
  train <- sample_df[1:n_train,]
  test <- sample_df[(n_train+1):dim(sample_df)[1],]
  
  train_coords <- as.matrix(train[,1:2])
  test_coords <- as.matrix(test[,1:2])
  
  #if (length(unique(rbind(train_coords,test_coords))) != length(rbind(train_coords,test_coords))) {
  #  print("We got em")
  #  break
  #}
  
  #Baselin Mean Absolute error is the null model
  baseline_mae <- mean(abs(test$Observed - mean(train$Observed)))
  
  #Predictions from the Dynamical System
  dynamical_mae <- mean(abs(test$Observed - test$Estimated))
  
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
  starting <- list(sigma.sq = 1, phi = 3)
  tuning <- list(sigma.sq = 1, phi = 1)
  priors <- list(sigma.sq.ig = c(0.1,0.1),  phi.unif = c(0.5,3.5)) 
  
  #Actual Model Fitting
  precondition_fit <- spLM(Estimated ~ 1, data = rbind(train,test), coords = rbind(train_coords,test_coords), starting = starting, priors = priors,
                           tuning = tuning, cov.model= 'gaussian', n.samples = 1000)
  
  ##Moment Match our posterior on sigma square to an IG
  y <- precondition_fit$p.theta.samples[,'sigma.sq']
  m <- mean(y)
  v <- var(y)
  alpha <- 2 + m^2 / v
  beta <- m * (1 + m^2/v)
  
  print(c(alpha, beta))
  
  ###
  #STEP 2: Use the new prior on sig.sq
  ###
  starting <- list(sigma.sq = 0.1, tau.sq = 0.1, phi = 3)
  tuning <- list(sigma.sq = 1, tau.sq = 1, phi = 1)
  priors <- list(sigma.sq.ig = c(alpha, 1/beta), tau.sq.ig = c(0.1,0.1), phi.unif = c(0.5,3.5))
  
  conditioned_fit <- spLM(Observed ~ 1, data = train, coords = train_coords, starting = starting, priors = priors,
                          tuning = tuning, cov.model= 'gaussian', n.samples = 1000)
  
  ##Evaluate the Model
  #Get the posterior mean on each test point
  conditioned_pred <- spPredict(conditioned_fit, test_coords, pred.covars = as.matrix(rep(1,dim(test_coords)[1])))
  conditioned_pred <- rowMeans((conditioned_pred$p.y.predictive.samples))
  
  #Calculate mean absolute error
  conditioned_mae <- mean(abs(conditioned_pred - test$Observed))
  
  ####
  # GP fit using spLM on residuals of Dynamical System model
  ####
  starting <- list(beta=c(1,1), sigma.sq = 0.1, tau.sq = 0.1, phi = 3)
  tuning <- list(sigma.sq = 1, tau.sq = 1, phi = 1)
  priors <- list('beta.flat', sigma.sq.ig = c(0.1,0.1), tau.sq.ig = c(0.1,0.1), phi.unif = c(0.5,3.5))
  
  #Actual model fitting
  resid_fit <- spLM(Observed ~ 1 + Estimated, data = test, coords = train_coords, starting = starting, priors = priors,
                       tuning = tuning, cov.model= 'gaussian', n.samples = 1000)
  
  ##Evaluate the Model
  #Get the posterior mean on each test point
  resid_pred <- spPredict(resid_fit, test_coords, pred.covars = matrix(c(rep(1,dim(test_coords)[1]), test$Estimated), ncol=2))
  resid_pred <- rowMeans((resid_pred$p.y.predictive.samples))
  
  #Calculate mean absolute error
  resid_mae <- mean(abs(resid_pred - test$Observed))
  
  ####
  # GP fit using spLM on residuals of the dynamical system AND with priors derived from the differential equation.
  ####
  # Use the prior from above on sig.sq again
  starting <- list(beta=c(1,1), sigma.sq = 0.1, tau.sq = 0.1, phi = 3)
  tuning <- list(sigma.sq = 1, tau.sq = 1, phi = 1)
  priors <- list(sigma.sq.ig = c(alpha, 1/beta), tau.sq.ig = c(0.1,0.1), phi.unif = c(0.5,3.5))
  
  #Actual model fitting
  resid_cond_fit <- spLM(Observed ~ 1 + Estimated, data = test, coords = train_coords, starting = starting, priors = priors,
                    tuning = tuning, cov.model= 'gaussian', n.samples = 1000)
  
  ##Evaluate the Model
  #Get the posterior mean on each test point
  resid_cond_pred <- spPredict(resid_cond_fit, test_coords, pred.covars = matrix(c(rep(1,dim(test_coords)[1]), test$Estimated), ncol=2))
  resid_cond_pred <- rowMeans((resid_cond_pred$p.y.predictive.samples))
  
  #Calculate mean absolute error
  resid_cond_mae <- mean(abs(resid_cond_pred - test$Observed))
  
  #print(c(baseline_mae, dynamical_mae, standard_mae, conditioned_mae))
  
  results <- c(results, c(baseline_mae, dynamical_mae, standard_mae, conditioned_mae, resid_mae, resid_cond_mae))
}

print(c(mean(results[seq(1,length(results), 6)]), mean(results[seq(2,length(results), 6)]),mean(results[seq(3,length(results), 6)]), mean(results[seq(4,length(results), 6)]), mean(results[seq(5,length(results), 6)]), mean(results[seq(6,length(results), 6)])))

system("zenity --title=\"R script info\" --text=\"Simulations Complete.\" --info") # for GTK dialog