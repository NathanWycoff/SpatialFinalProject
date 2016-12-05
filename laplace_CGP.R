## Algorithm to do a GP for classification using laplace approx
##
##Author -- Nathan Wycoff
##

##Generate some binary data
n <- 150
k1 <- 50
k2 <- 50
X <- c(rnorm(k1, -10, 7), rnorm(k2, 10,7), rnorm(n-k1-k2, 30,7))

hist(X[1:k1], col = rgb(0.8,0,0,0.5), xlim = range(X))
hist(X[(k1+1):(k1+k2)], col = rgb(0,0.8,0,0.5), add = TRUE)
hist(X[(k1+k2+1):n], col = rgb(0.8,0,0,0.5), add = TRUE)
y <- rep(c(-1,1,-1), c(k1,k2,n-k1-k2))

#Train-test split
t_ind <- sample(1:n, floor(n / 2))
train <- cbind(X,y)[t_ind,]
test <- cbind(X,y)[-t_ind,]

train <- data.frame(train)
colnames(train) <- c('X', 'y')
test <- data.frame(test)
colnames(test) <- c('X', 'y')

n <- floor(n / 2)


####
## Pre-built model
####
fit <- gausspr(I(as.factor(y)) ~ X, data = train, type = 'classification', kernel = 'rbfdot', kpar = list(sigma=1))
mean(predict(fit, test) == test$y)

##Look at logistic reg
fit_log <- glm(I(as.factor(y)) ~ X + 0, data = train, family = binomial)

preds <- sapply((predict(fit_log, newdata = test) > 0.5), function(x) ifelse(x,1,-1))

log_acc <- mean(preds == test$y)

##Now do a GP for classification with laplace approx
#Our covar func
k <- function(x1,x2) exp(-norm(x1 - x2, '2')^2)#Standard gaussian covar func

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


f_pred <- nathans_gp_class()
gp_preds <- sapply(f_pred > 0, function(x) ifelse(x,1,-1))
gp_acc <- mean(gp_preds == test$y)