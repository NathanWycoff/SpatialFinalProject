## Spatial Analysis of the interaction of the grey and red squirrel in easter UK.
##
##Author -- Nathan Wycoff

require(h5)
require(kernlab)

######
## Data Manip
#####

###Read in data
data <- h5file('/home/nathan/Documents/CurrentTerm/Space/Final Project/reynolds_data.h5', 'r')

grey <- as.array(data['grey_squirrel'][,,])
red <- as.array(data['red_squirrel'][,,])

h5close(data)

###Translate data to form suitable for GP fitting.
#Translate ndarray to n-4 array of x,y,t and response.
grey_df <- c()
for (t in 1:dim(grey)[1]) {
  for (i in 1:dim(grey)[2]) {
    for (j in 1:dim(grey)[3]) {
      grey_df <- rbind(grey_df, c(t,i,j,grey[t,i,j]))
    }
  }
}

red_df <- c()
for (t in 1:dim(red)[1]) {
  for (i in 1:dim(red)[2]) {
    for (j in 1:dim(red)[3]) {
      red_df <- rbind(red_df, c(t,i,j,red[t,i,j]))
    }
  }
}

#Make into df and rename cols
grey_df <- data.frame(grey_df)
colnames(grey_df) <- c('Time', 'Easting', 'Northing', 'Presence')
red_df <- data.frame(red_df)
colnames(red_df) <- c('Time', 'Easting', 'Northing', 'Presence')

#Visualize Squirrel Data
par(mfrow=c(2,2))
image(grey[1,,], main = 'Initial Red Squirrel Values')
image(grey[9,,], main = 'Final Red Squirrel Values')
image(red[1,,], main = 'Initial Grey Squirrel Values')
image(red[9,,], main = 'Final Grey Squirrel Values')

##Train test split
#Train data, all but the lat time point
train_ind <- which(grey_df$Time < 9)
grey_df_train <- grey_df[train_ind,]
red_df_train <- grey_df[train_ind,]

#Test data, only the last time point.
test_ind <- which(grey_df$Time == 9)
grey_df_test <- grey_df[test_ind,]
red_df_test <- grey_df[test_ind,]

#######
## Dynamical System Evaluation for models that require it
#######

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

time_interval <- c(0,9)
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
#image(G[,,t_size])
#image(R[,,t_size])




#######
## Kernel Evaluations
#######


### GP Only
#For grey squirrel
##Package data so it can be consumed by my custom function
X <- grey_df_train[,1:3]
XX <- grey_df_test[,1:3]
y <- grey_df_train$Presence
solx <- c()
for (i in 1:dim(X)[1]) {
  #Just to fix 0 indexing
  time_ind <- ifelse(X[i,1]==1, (X[i,1]-1)/delta_t+1 ,(X[i,1]-1)/delta_t)
  solx <- c(solx, G[X[i,3], X[i,2], time_ind])
}
solxx <- c()
for (i in 1:dim(XX)[1]) {
  #Just to fix 0 indexing
  time_ind <- ifelse(XX[i,1]==1, (XX[i,1]-1)/delta_t+1 ,(XX[i,1]-1)/delta_t)
  solxx <- c(solxx, G[XX[i,3], XX[i,2], time_ind])
}

#Do prediction for GP linear model
fit_grey_dyn <- glm(grey_df_train$Presence ~ solx, family = binomial)
mean((predict(fit_grey_dyn, newdata = data.frame(solx=solxx)) > 0.5) == grey_df_test$Presence)

#For Red squirrel
##Package data so it can be consumed by my custom function
X <- red_df_train[,1:3]
XX <- red_df_test[,1:3]
y <- red_df_train$Presence
solx <- c()
for (i in 1:dim(X)[1]) {
  #Just to fix 0 indexing
  time_ind <- ifelse(X[i,1]==1, (X[i,1]-1)/delta_t+1 ,(X[i,1]-1)/delta_t)
  solx <- c(solx, R[X[i,3], X[i,2], time_ind])
}
solxx <- c()
for (i in 1:dim(XX)[1]) {
  #Just to fix 0 indexing
  time_ind <- ifelse(XX[i,1]==1, (XX[i,1]-1)/delta_t+1 ,(XX[i,1]-1)/delta_t)
  solxx <- c(solxx, R[XX[i,3], XX[i,2], time_ind])
}


#Do prediction for GP linear model
fit_red_dyn <- glm(red_df_train$Presence ~ solx, family = binomial)
mean((predict(fit_red_dyn, newdata = data.frame(solx=solxx)) > 0.5) == red_df_test$Presence)

fit_grey <- gausspr(I(as.factor(Presence)) ~ Time + Easting + Northing, data = grey_df_train, type = 'classification', kernel = 'rbfdot', kpar = list(sigma=1))
mean(predict(fit_grey, grey_df_test) == grey_df_test$Presence)

#For red squirrel
fit_red <- gausspr(I(as.factor(Presence)) ~ Time + Easting + Northing, data = red_df_train, type = 'classification', kernel = 'rbfdot', kpar = list(sigma=1))
mean(predict(fit_red, red_df_test) == red_df_test$Presence)

### DS with GP resids
#For grey squirrel
fit_grey <- gausspr(I(as.factor(Presence)) ~ Time + Easting + Northing, data = grey_df_train, type = 'classification', kernel = 'rbfdot', kpar = list(sigma=1))
mean(predict(fit_grey, grey_df_test) == grey_df_test$Presence)

#For red squirrel
fit_red <- gausspr(I(as.factor(Presence)) ~ Time + Easting + Northing, data = red_df_train, type = 'classification', kernel = 'rbfdot', kpar = list(sigma=1))
mean(predict(fit_red, red_df_test) == red_df_test$Presence)


###Nathan's Special Kernel
nathans_gp_class <- function(y, X, XX, solx, solxx) {
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
    covs <- c(covs, cov(y1,y2))
    lens <- c(lens, length(y1))
  }
  
  #Use weights to get multiple points for bin in the histo
  weights <- floor(lens / min(lens))#Weigh cov estimates approximately by their number of constituents.
  hyperdata <- rep(covs,weights)
  hyperdistance <- rep(bins[2:length(bins)],weights)
  
  #Fit GP to these covars
  length_scale <- (bins[2] - bins[1])*1000#Make the length scale 1 bin
  gpi <- newGP(matrix(hyperdistance), hyperdata, d=length_scale, g=0.1, dK=FALSE)
  ga <- garg(list(mle=TRUE, min =0.01, max=15, start = 5), hyperdata)  #Make max bigger?
  mle<-mleGP(gpi, param="g", tmin=ga$min, tmax=ga$max, ab=ga$ab)
  grid <- matrix(seq(min(hyperdistance), max(hyperdistance)))
  p <- predGP(gpi, grid, lite=TRUE)$mean
  
  #Visualize nonparam kernel
  plot(jitter(hyperdistance,amount=1), jitter(hyperdata, amount = 0.1))
  points(grid, p, type = 'l')
  
  #Store the predicted values and spit out each step.
  np_k <- function(x1,x2) {
    #print("hehehoho")
    #print(x1)
    h <- norm(x1 - x2, '2')
    return(p[min(which(grid > h))])
  }
  
  #Our covar mat
  K <- matrix(rep(0, n*n), nrow = n)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i,j] <- np_k(X[i], X[j])
    }
  }
  
  
  ##Begin Laplace Approx Algorithm from Rasmussen and Williams pg 46
  get_mode <- function(K, y) {
    n <- length(y)
    #Initializations
    f <- rnorm(n,0,0.1)
    
    #Newton iteration
    tol <- 0.1
    while (TRUE) {
      
      ##Update our stuff
      W <- diag(pnorm(train$y * as.numeric(f)) * (1-pnorm(y * as.numeric(f))))
      L <- chol(diag(n) + sqrt(W) %*% K %*% sqrt(W))
      b <- W %*% f + (y+1)/2 - pnorm(y*f)
      a <- b - sqrt(W) %*% solve(t(L), solve(L, sqrt(W) %*% K %*% b))
      f_new <- K %*% a
      
      ##Check for convergence
      if (norm(f_new - f, 'I') < tol) {
        break
      }
      f <- f_new
    }
    return(f)
  }
  
  ##Prediction using laplace approx
  laplace_predict <- function(f, X, XX, y) {
    #k_new <- t(sapply(test$X, function(x_new) sapply(train$X, function(x) k(x_new, x))))
    f_new <- c()
    for (x_new in XX) {
      W <- diag(pnorm(y * as.numeric(f)) * (1-pnorm(y * as.numeric(f))))
      L <- chol(diag(n) + sqrt(W) %*% K %*% sqrt(W))
      k_new <- sapply(X, function(x) k(x_new, x))
      ###############################################
      jacob_diag <- (train$y+1)/2 - pnorm(train$y*f)
      f_new <- c(f_new, t(k_new)  %*% jacob_diag)
      ###############################################
    }
    return(f_new)
  }
  
  f <- get_mode(K, y)
  
  f_pred <- laplace_predict(f, X, XX, y)
  
  return(f_pred)
}

train_ind <- which(grey_df$Time < 9)
grey_df_train <- grey_df[train_ind,]
red_df_train <- grey_df[train_ind,]

#Test data, only the last time point.
test_ind <- which(grey_df$Time == 9)
grey_df_test <- grey_df[test_ind,]
red_df_test <- grey_df[test_ind,]


##Package data so it can be consumed by my custom function
X <- grey_df_train[,1:3]
XX <- grey_df_test[,1:3]
y <- grey_df_train$Presence
solx <- c()
for (i in 1:dim(X)[1]) {
  #Just to fix 0 indexing
  time_ind <- ifelse(X[i,1]==1, (X[i,1]-1)/delta_t+1 ,(X[i,1]-1)/delta_t)
  solx <- c(solx, G[X[i,3], X[i,2], time_ind])
}
solxx <- c()
for (i in 1:dim(XX)[1]) {
  #Just to fix 0 indexing
  time_ind <- ifelse(XX[i,1]==1, (XX[i,1]-1)/delta_t+1 ,(XX[i,1]-1)/delta_t)
  solxx <- c(solxx, G[XX[i,3], XX[i,2], time_ind])
}

f_pred <- nathans_gp_class(y, X, XX, solx, solxx)
gp_preds <- sapply(f_pred > 0, function(x) ifelse(x,1,-1))
gp_acc <- mean(gp_preds == test$y)

