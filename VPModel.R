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

par(mfrow=c(1,2))
image(S_used, 'main' = 'USED')
image(S_true, main = 'TRUE')
par(mfrow=c(1,1))

#######
## Use hyper GP to estimate test data using my own special sauce.
######

#Do GP regression (Point estimates only), using a kernel function
# nonparametrically fit to solutions of a dynamical system
#
# X -- Locations of training data
# y -- Response at locations
# XX -- Locations to predict at
# solx -- solution of dynamical system at training points
# solxx -- solution of dynamical system at testing points
nathans_gp_reg <- function(y, X, XX, solx, solxx) {
  #####
  ## Do hyper GP
  #####
  
  XXX <- rbind(X,XX)
  
  ##Get Covariance
  D <- dist(XXX)#Create a distance matrix
  bins <- hist(D, plot = FALSE)$breaks#create bins
  
  #Go over it all and get the covariances.
  n <- dim(XXX)[1]
  covs <- c()
  lens <- c()
  hyperdistance <- c()
  for (bi in 1:(length(bins)-1)) {
    y1 <- c()
    y2 <- c()
    b <- bins[bi+1]
    b1 <- bins[bi]
    for (j in 1:n) {
      for (i in 1:j) {
        if (i == j) {
          next
        }
        ind <- n*(i-1) - i *(i-1)/2 + j - i
        if (D[ind] < b && D[ind] >= b1) {
          y1 <- c(y1, c(solx, solxx)[i])
          y2 <- c(y2, c(solx, solxx)[j])
        }
      }
    }
    if (length(y1) > 1) {
      covs <- c(covs, cov(y1,y2))
      lens <- c(lens, length(y1))
      hyperdistance <- c(hyperdistance, b1)
    }
  }
  
  #Use weights to get multiple points for bin in the histo
  weights <- floor(lens / min(lens))#Weigh cov estimates approximately by their number of constituents.
  hyperdata <- rep(covs,weights)
  hyperdistance <- rep(hyperdistance,weights)
  
  ##Get rid of places with only one datum
  #hyperdata <- hyperdata[which(rep(weights > 1, weights))]
  #hyperdistance <- hyperdistance[which(rep(weights > 1, weights))]
  
  #Fit GP to these covars
  length_scale <- (bins[2] - bins[1])*1000#Make the length scale 1 bin
  gpi <- newGP(matrix(hyperdistance), hyperdata, d=length_scale, g=0.1, dK=FALSE)
  ga <- garg(list(mle=TRUE, min =0.01, max=15, start = 5), hyperdata)  #Make max bigger?
  mle<-mleGP(gpi, param="g", tmin=ga$min, tmax=ga$max, ab=ga$ab)
  grid <- matrix(seq(min(hyperdistance), max(hyperdistance)*5, length.out = 250))
  p <- predGP(gpi, grid, lite=TRUE)$mean
  
  #Visualize nonparam kernel
  plot(jitter(hyperdistance,amount=1), jitter(hyperdata, amount = 0.01))
  points(grid, p, type = 'l')
  
  #Store the predicted values and spit out each step.
  np_k <- function(x1,x2) {
    #print("hehehoho")
    #print(x1)
    h <- norm(x1 - x2, '2')
    return(p[min(which(grid > h))])
  }
  
  ##Create Covariance Matrix
  n <- length(y)
  K <- matrix(rep(0, n*n), nrow = n)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i,j] <- np_k(X[i,], X[j,])
    }
  }
  
  ##Do prediction
  sigma2 <- 10
  
  #Some precacluations
  L <- chol(K + diag(sigma2, n))
  alpha <- solve(t(L), solve(L, y))
  
  #Iterate over each data point and spit out our predictions
  f <- c()
  for (i in 1:dim(XX)[1]) {
    x_new <- XX[i,]
    k_new <- sapply(1:dim(X)[1], function(j) np_k(x_new, X[j,]))
    f <- c(f, t(k_new) %*% alpha)
  }
  
  return(f)
}


######
# Now we evaluate several GP based methods on noise corrupted data
######
#Params
sigma <- 0.1#Variance of noise corruption
iters <- 100

#Get some train and test data.
n_samples <- 100
hyper_results <- c()
for (n_samples in seq(10,150,by=10)) {
  n_train <- floor(n_samples/2)
  results <- c()
  
  for (ayy in 1:iters) {
    tryCatch({
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
      baseline_mae <- mean((test$Observed - mean(train$Observed))^2)
      
      #Predictions from the Dynamical System
      dynamical_mae <- mean((test$Observed - test$Estimated)^2)
      
      ####
      # Regular GP, fit using spLM
      ####
      ##Fit the model
      #Prepare for the model
      starting <- list(sigma.sq = 0.1, tau.sq = 0.1, phi = 3)
      tuning <- list(sigma.sq = 1, tau.sq = 1, phi = 1)
      priors <- list(sigma.sq.ig = c(0.1,0.1), tau.sq.ig = c(0.1,0.1), phi.unif = c(0.15,3.5))
      
      #Actual model fitting
      standard_fit <- spLM(Observed ~ 1, data = test, coords = train_coords, starting = starting, priors = priors,
                           tuning = tuning, cov.model= 'gaussian', n.samples = 1000)
      
      ##Evaluate the Model
      #Get the posterior mean on each test point
      standard_pred <- spPredict(standard_fit, test_coords, pred.covars = as.matrix(rep(1,dim(test_coords)[1])))
      standard_pred <- rowMeans((standard_pred$p.y.predictive.samples))
      
      #Calculate mean absolute error
      standard_mae <- mean((standard_pred - test$Observed)^2)
      
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
      priors <- list(sigma.sq.ig = c(0.1,0.1),  phi.unif = c(0.1,10)) 
      
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
      
      ##Store our phi estimates
      phi_hyper <- quantile(precondition_fit$p.theta.samples[,'phi'], c(0, 1))
      
      ###
      #STEP 2: Use the new prior on sig.sq
      ###
      starting <- list(sigma.sq = 0.1, tau.sq = 0.1, phi = 3)
      tuning <- list(sigma.sq = 1, tau.sq = 1, phi = 1)
      priors <- list(sigma.sq.ig = c(alpha, 1/beta), tau.sq.ig = c(0.1,0.1), phi.unif = phi_hyper)
      
      conditioned_fit <- spLM(Observed ~ 1, data = train, coords = train_coords, starting = starting, priors = priors,
                              tuning = tuning, cov.model= 'gaussian', n.samples = 1000)
      
      ##Evaluate the Model
      #Get the posterior mean on each test point
      conditioned_pred <- spPredict(conditioned_fit, test_coords, pred.covars = as.matrix(rep(1,dim(test_coords)[1])))
      conditioned_pred <- rowMeans((conditioned_pred$p.y.predictive.samples))
      
      #Calculate mean absolute error
      conditioned_mae <- mean((conditioned_pred - test$Observed)^2)
      
      ####
      # GP fit using spLM on residuals of Dynamical System model
      ####
      starting <- list(beta=c(1,1), sigma.sq = 0.1, tau.sq = 0.1, phi = 3)
      tuning <- list(sigma.sq = 1, tau.sq = 1, phi = 1)
      priors <- list('beta.flat', sigma.sq.ig = c(0.1,0.1), tau.sq.ig = c(0.1,0.1), phi.unif = c(0.15,3.5))
      
      #Actual model fitting
      resid_fit <- spLM(Observed ~ 1 + Estimated, data = test, coords = train_coords, starting = starting, priors = priors,
                        tuning = tuning, cov.model= 'gaussian', n.samples = 1000)
      
      ##Evaluate the Model
      #Get the posterior mean on each test point
      resid_pred <- spPredict(resid_fit, test_coords, pred.covars = matrix(c(rep(1,dim(test_coords)[1]), test$Estimated), ncol=2))
      resid_pred <- rowMeans((resid_pred$p.y.predictive.samples))
      
      #Calculate mean absolute error
      resid_mae <- mean((resid_pred - test$Observed)^2)
      
      ####
      # GP fit using spLM on residuals of Dynamical System model
      ####
      X <- cbind(train$x, train$t)
      XX <- cbind(test$x, test$t)
      y <- train$Observed
      solx <- train$Estimated
      solxx <- test$Estimated
      
      pred <- nathans_gp_reg(y,X,XX,solx,solxx)
      
      nate_mae <- mean((pred - test$Observed)^2)
      
      results <- rbind(results, c(baseline_mae, dynamical_mae, standard_mae, conditioned_mae, resid_mae, nate_mae))
    }, error = function(e) {})
  }
  results <- as.data.frame(results)
  results <- results[complete.cases(results),]
  hyper_results <- rbind(hyper_results, c(colMeans(results), n_samples))
}

hyper_results <- as.data.frame(hyper_results)
colnames(hyper_results) <- c("Null", "DiffEQ","StandardGP", "ConditionedGP", "GPonDiffEQ", "NPGP","N")
rownames(hyper_results) <- 1:dim(hyper_results)[1]

##Plot results
#Get the height for the plot
height <- c(min(hyper_results), max(hyper_results[,1:(dim(hyper_results)[2]-1)]))
height[2] <- height[2] + 1/2 * (height[2] - height[1])

#Plot lines of each method
plot(hyper_results$N, hyper_results[,1], type='l', col = 'red', 
     ylim = height,
     xlab = 'Sample Size (train + test)', ylab = 'Out of Sample MSE',
     main = 'Comparison of Methods for Model with Parameter Error')

points(hyper_results$N, hyper_results[,2], type='l', col = 'magenta')
points(hyper_results$N, hyper_results[,3], type='l', col = 'green')
points(hyper_results$N, hyper_results[,4], type='l', col = 'purple')
points(hyper_results$N, hyper_results[,5], type='l', col = 'blue')
points(hyper_results$N, hyper_results[,6], type='l', col = 'orange')

#Make the legend
legend(x = 50, y = height[2], legend = colnames(hyper_results)[-dim(hyper_results)[2]], cex = 0.70,
       lty=c(1,1), lwd=c(2.5,2.5),col=c('red','magenta','green','purple','blue','orange'))

#print(c(mean(results[seq(1,length(results), 5)]), mean(results[seq(2,length(results), 5)]),mean(results[seq(3,length(results), 5)]), mean(results[seq(4,length(results), 5)], na.rm = TRUE), mean(results[seq(5,length(results), 5)])))

system("zenity --title=\"R script info\" --text=\"Simulations Complete.\" --info") # for GTK dialog