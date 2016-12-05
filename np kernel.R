##Author -- Nathan Wycoff
##Description -- Nonparametric kernel function trained no 

#  -- I don't think we currently have a burn in
#  -- Do the nonparameteric kernel function model

#Requires spBayes
require(spBayes)
require(laGP)
require(kernlab)

##Problem Parameters
time_interval <- c(0,10)
dist_interval <- c(0,50)
D <- 10#Actual Diffusion Coefficient
alpha <- 0#Bonding Coefficient

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
R <- matrix(rep(0, x_size * t_size), nrow = x_size)#Initialize a matrix of zeros

#Enforce boundary and initial conditions
S[2:49,1] <- 0.1#Enforce Initial Conditions, point mass at 25
S[1,1] <- 0# Boundary Condition on one end
S[50,1] <- 0# Boundary Condition on the other
R[12,1] <- 200#Enforce Initial Conditions, point mass at 25
R[36,1] <- 200#Enforce Initial Conditions, point mass at 25

#Do iterations
for (j in 1:(t_size-1)) {
  for (i in 2:(x_size-1)) {
    #Update chemical S
    interaction <- - alpha * delta_t * R[i,j] * S[i,j]
    S[i,j+1] <- S[i,j] + D*(delta_t / delta_x^2) * (S[i+1,j] - 2 * S[i,j] + S[i-1,j]) + interaction
    
    #Update Chemical R
    interaction <- - alpha * delta_t * R[i,j] * S[i,j]
    R[i,j+1] <- R[i,j] + D*(delta_t / delta_x^2) * (R[i+1,j] - 2 * R[i,j] + R[i-1,j]) + interaction
  }
}

image(S)
image(R)

S_true <- S
S_used <- S


########
## Do our GP
########
#Params
sigma <- 0.1#Variance of noise corruption
iters <- 100

#Get some train and test data.
n_samples <- 100
n_train <- 50
results <- c()

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

#Normalize Data
sample_df$Estimated <- (sample_df$Estimated - mean(sample_df$Estimated)) / sd(sample_df$Estimated)
sample_df$Observed <- (sample_df$Observed - mean(sample_df$Observed)) / sd(sample_df$Observed)

#Train test split
train <- sample_df[1:n_train,]
test <- sample_df[(n_train+1):dim(sample_df)[1],]



#######
###Create hyper GP on dynamical sys
#######

##Get Covariance
D <- dist(matrix(c(train$x,test$x,test$t,test$t), ncol = 2))#Create a distance matrix
bins <- hist(D, plot = FALSE)$breaks#create bins

#Go over it all and get the covariances.
n <- length(c(train$x,test$x))
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
        y1 <- c(y1, c(train$Estimated, test$Estimated)[i])
        y2 <- c(y2, c(train$Estimated, test$Estimated)[j])
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
    k_new <- sapply(1:dim(X)[1], function(i) np_k(x_new, X[i,]))
    f <- c(f, t(k_new) %*% alpha)
  }
  
  return(f)
}


X <- cbind(train$x, train$t)
XX <- cbind(test$x, test$t)
y <- train$Observed
solx <- train$Estimated
solxx <- test$Estimated

pred <- nathans_gp_reg(y,X,XX,solx,solxx)
