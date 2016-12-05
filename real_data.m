%Load requirements
pkg load statistics
pkg load gpcf_sexp


%Generate Data
n = 50;
k = 25
X = vertcat(normrnd(-10,7,[k,1]),normrnd(10,7,[n-k,1]));
y = vertcat(repmat(0,k), repmat(1,n-k));

% Training test split
train_ind = randsample(n, floor(n/2));
X_train = X(train_ind);
y_train = y(train_ind);
test_ind = setdiff(1:n, train_ind);
X_test = X(test_ind);
y_test = y(test_ind);

% Create covariance functions
gpcf1 = gpcf_sexp('lengthScale', [0.9], 'magnSigma2', 10);

% Set the prior for the parameters of covariance functions 
pl = prior_unif();
pm = prior_sqrtunif();
gpcf1 = gpcf_sexp(gpcf1, 'lengthScale_prior', pl,'magnSigma2_prior', pm); %

% Create the GP structure (type is by default FULL)
gp = gp_set('cf', gpcf1, 'lik', lik_probit(), 'jitterSigma2', 1e-9);
%gp = gp_set('cf', gpcf1, 'lik', lik_logit(), 'jitterSigma2', 1e-9);

% ------- Laplace approximation --------
% Set the approximate inference method 
% (Laplace is default, so this could be skipped)
gp = gp_set(gp, 'latent_method', 'Laplace');

% Set the options for the scaled conjugate optimization
opt=optimset('TolFun',1e-3,'TolX',1e-3,'MaxIter',20,'Display','iter');
% Optimize with the scaled conjugate gradient method
gp=gp_optim(gp,x,y,'optimf',@fminscg,'opt',opt);

% Make predictions
[Ef_la, Varf_la, Ey_la, Vary_la, Py_la] = gp_pred(gp, x, y, xt, 'yt', ones(size(xt,1),1) );

