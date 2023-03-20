# chapter 8 in Blangiardo & Cameletti (2015)
library(INLA)
library(dplyr)

###---Bivariate model for spatially misaligned data---###

# To illustrate misalignment, we will simulate 15 and 20 locations for the response and covariate, respectively for a square domain between 0 and 1
n_y <- 15
n_x <- 20
set.seed(20)
loc_y <- cbind(runif(n_y), runif(n_y))
set.seed(1)
loc_x <- cbind(runif(n_x), runif(n_x))
plot(loc_y[,1:2])
points(loc_x[,1:2], col = "red", type = "p")

# joint bivariate model for both covariate and response, allows to account for all sources of uncertainty

# The same spatial process will be used in more than one linear predictor which we can used the `copy` feature to do

####################################################################################################
# joint model with gaussian distributions (both covaraite and response follow a normal distribution)
####################################################################################################

## Data simulation

# fix the parameters for both the random fields (the fields are the covariate mean (xi) and the additional spatial effect (u))
# These are the values of the spatial variance and the scale parameter k of the matern covariance

kappa_xi <- 5
sigma2_xi <- 0.5
kappa_u <- 7
sigma2_u <- 0.3

#In practice, to simulate values from a multivariate Normal distribution with covariance matrix deined by the Matérn covariance function, we use a standard algorithm based on the Cholesky factorization

simulate_GF <- function(coords, kappa, variance, lambda = 1) {
    #Compute the number of locations
    n <- nrow(coords)
    #Compute the distance matrix
    dist.m <- as.matrix(dist(coords))
    #Compute the Matern correlation matrix
    cor.m <- 2^(1-lambda)/gamma(lambda)*(kappa*dist.m)^lambda*besselK(x=dist.m*kappa, nu=lambda)
    diag(cor.m) <- 1
    #Compute the covariance matrix
    Sigma <- variance * cor.m
    #Simulate date using standard algorithm based on Cholesky fact.
    c(chol(Sigma) %*% rnorm(n=n,mean=0,sd=1))
}

# By using the simulate_GF function and the values chosen for the Matérn
# parameters, we simulate the values for the spatial process u(s) at the ny locations
# with:
set.seed(223)
u <- simulate_GF(coords = loc_y,
                 kappa = kappa_u,
                 variance = sigma2_u)
length(u)

# and the realizations of the process xi(s) at both set of spatial sites (first for the nx
# covariate locations) with:

set.seed(233)
xi <- simulate_GF(coords = rbind(loc_x, loc_y),
                  kappa=kappa_xi, variance=sigma2_xi)
length(xi)


# It is worth noting that we need the values xi also in the ny response locations in
# order to be able to simulate the observations using Eq. (8.3). However, when we it
# the model we will consider only the covariate values available at the nx locations.

# Now, we set the values for b0 (intercept), beta1 and the variances sigma2_e and sigma2_x:
b0 <- 10
beta1 <- -0.5
sigma2_e <- 0.16
sigma2_x <- 0.25

# Finally, the simulated values for the covariate and the response are obtained through
set.seed(4455)
x <- xi[1:n_x] + rnorm(n=n_x, mean = 0, sd = sqrt(sigma2_x))
set.seed(5544)
y <- b0 + beta1*xi[n_x + 1:n_y] + u + rnorm(n=n_y, mean = 0, sd=sqrt(sigma2_e))

##--Model fitting--##

# first step is to make mesh - this will be simple with only one outer extension
mesh <- inla.mesh.2d(loc=rbind(loc_x,loc_y), max.edge = 0.15,
                     cutoff = 0.03, offset = 0.1)
plot(mesh)

# Now define one SPDE model which will be used for both spatial effects xi(s) and u(s)
spde <- inla.spde2.matern(mesh=mesh, alpha=2)


# Due to the misalignment, in order to link the locations with the mesh we have to
# define two projection matrices:
A_x <- inla.spde.make.A(mesh = mesh, loc=loc_x)
dim(A_x)
A_y <- inla.spde.make.A(mesh = mesh, loc=loc_y)
dim(A_y)

# The model we are implementing is characterized by two likelihoods, one fo rthe covariate and noe for the response
# Here the response data can be rewriten the response data as a matrix where the # of columns is given by the # of likelihoods
# This data structure has to be taken into account using the inla.stack function - make a stack for each 

#stack for the covariate
stk_x <- inla.stack(data = list(y=cbind(x,NA)), # col 1 is the covariate values and col 2 is NA
                    effects = list(xi.field=1:spde$n.spde), # just an index of the spde 
                    A=list(A_x), # A matrix for covariate
                    tag="est.x"
                    )

# As we do not have replicates or any grouping, the index set for the spatial effect
# (xi.field) is just a sequence of integers from 1 to the number of nodes on the mesh.

# stack for the response is more complicated (it includes also the intercept and the spatial field u(s))
stk_y <- inla.stack(data = list(y=cbind(NA,y)), # col 1 is NAs and col 2 is the response variable values
                    effects=list(
                        list(u.field=1:spde$n.spde, x.field=1:spde$n.spde), #index for u(s) and xi(s) spde
                        list(intercept=rep(1,n_y)) #intercept (is there an intercept for each location?)
                    ),
                    A=list(A_y,1), #why does this have a "1" but not in stk_x
                    tag="est.y"
                    )

# observed response values yi in the second column of the data matrix. The
# effect argument is a list with two elements: the first includes two index sets for
# the two spatial random effects (with same projector matrix), while the second is the intercept vector

# now its possible to define the full stack object
stk <- inla.stack(stk_x,stk_y)


# The estimation procedure is carried out using the R-INLA "copy" feature which is adopted when a latent field is needed more than
# once in the model formulation.
# In this example, the spatial field xi(s) is used both in the speciication of the covariate and response linear predictor
# Need to define the "copy" feature in the formula

formula <- y ~ -1 + f(xi.field, model = spde) + # spatial field for covariate and -1 is added to remove R’s implicit intercept, which is replaced by the explicit +Intercept from when we created the stack
                    intercept + #now the response part which had an intercept and...
                    f(u.field, model=spde) + #the u(s)....
                    f(x.field, copy="xi.field", fixed = FALSE, #and the xi.field which we are copying from above 
                      hyper=list(theta=list(param=c(-1,10))))


# In this case, the random effect xi(s) is copied as it enters both in the linear predictor
# of xi and of yi. BUT the relationship between xi (covariate) and yi (response) involves B1 (coefficient of xi in the equation of yi) which is treated as hyperparameter in the estimation step. 
# we estimate the hyperparameter B1 (fixed=FALSE) and, in the prior speciication, we center the distribution on −1
# and assign a variance equal to 10 (this is suficient to include zero). The B1 hyperparameter
# is called theta in R-INLA and the default prior distribution is Gaussian

# We can also change the prior for sigma2_x and sigma2_e, choosing an informative logGamma with parameters 1 and 0.1 for both precisions
precprior <- list(theta=list(param=c(1, 0.1)))

# final step is the run of the inla function
output <- inla(formula, family = c("gaussian","gaussian"),#specifies a vector of two likelihoods
               data = inla.stack.data(stk),
               control.predictor = list(compute = TRUE,
                                        A=inla.stack.A(stk)),
               control.family = list(list(hyper=precprior), #same prior values for both the variances sigma2_x and sigma2_e
                                     list(hyper=precprior))
               )


##--Getting The Results--##

# the true values and the posterior summaries of the intercept (b0):
cbind(True = b0, output$summary.fixed[1, 1:5, drop=FALSE]) %>% round(4) #drop=FALSE inside[ ]is used to avoid the result to be coerced to a vector

# the true values and posterior summaries of the precisions (1/sigma2_x and 1/sigma2_e):
cbind(True = 1/c(Prec.x = sigma2_x, Prec.e = sigma2_e), output$summary.hyperpar[1:2,1:5])  %>% round(4)

# the true values and posterior summaries of the beta1 (i.e. coeifficent):

cbind(True=beta1, output$summary.hyperpar[5,1:5, drop=FALSE]) %>% round(4)

# For the precision, due to the asymmetry, it is better to look at the posterior marginal *distribution*. To plot these distributoins you can use the values from the  output$marginals.fixed or outputs$marginals.hyperpar (x is the values and y is the sd)


# The posterior summaries of te spatial parameters (for the process xi(s) and u(s)) can be extracted from the output using `inla.spde2.result`
xi_field <- inla.spde2.result(output, name="xi.field", spde)

u_field <- inla.spde2.result(output, name = "u.field", spde)

# In particular, the true values and posterior summaries for the spatial variances (sigma2_xi and sigma2_u) can be obtained:

cbind(True = c(sigma2.xi=sigma2_xi, sigma2.u=sigma2_u),
      mean=c(inla.emarginal(function(x) x, xi_field$marginals.var[[1]]),
             inla.emarginal(function(x) x, u_field$marginals.var[[1]])),
      rbind(inla.hpdmarginal(.95, xi_field$marginals.var[[1]]),
            inla.hpdmarginal(.95, u_field$marginals.var[[1]]))) %>% 
round(4)


# The true values for the range parameters for r_xi and r_u and the summary of the corresponding posterior:
cbind(True=c(range.xi=sqrt(8)/kappa_xi, range.u=sqrt(8)/kappa_u),
      mean=c(inla.emarginal(function(x) x, xi_field$marginals.range[[1]]),
             inla.emarginal(function(x) x, u_field$marginals.range[[1]])),
      rbind(inla.hpdmarginal(.95, xi_field$marginals.range[[1]]),
            inla.hpdmarginal(.95, u_field$marginals.range[[1]]))) %>% 
round(4)





#############################################
# joint model with non-gaussian distributions
#############################################

# Here we will consider both the covariate and the response to assume discrete values.
# The covariate is modeled as a Poisson distribution while the response is modeled as a Binomial 


##--Data simulation--##

# fix the parameters for both the random fields (the fields are the covariate mean (xi) and the additional spatial effect (u))
# These are the values of the spatial variance and the scale parameter k of the matern covariance

kappa_xi <- 5
sigma2_xi <- 0.5
kappa_u <- 7
sigma2_u <- 0.3

#In practice, to simulate values from a multivariate Normal distribution with covariance matrix deined by the Matérn covariance function, we use a standard algorithm based on the Cholesky factorization

simulate_GF <- function(coords, kappa, variance, lambda = 1) {
    #Compute the number of locations
    n <- nrow(coords)
    #Compute the distance matrix
    dist.m <- as.matrix(dist(coords))
    #Compute the Matern correlation matrix
    cor.m <- 2^(1-lambda)/gamma(lambda)*(kappa*dist.m)^lambda*besselK(x=dist.m*kappa, nu=lambda)
    diag(cor.m) <- 1
    #Compute the covariance matrix
    Sigma <- variance * cor.m
    #Simulate date using standard algorithm based on Cholesky fact.
    c(chol(Sigma) %*% rnorm(n=n,mean=0,sd=1))
}

# By using the simulate_GF function and the values chosen for the Matérn
# parameters, we simulate the values for the spatial process u(s) at the ny locations
# with:
set.seed(223)
u <- simulate_GF(coords = loc_y,
                 kappa = kappa_u,
                 variance = sigma2_u)
length(u)

# and the realizations of the process xi(s) at both set of spatial sites (first for the nx
# covariate locations) with:

set.seed(233)
xi <- simulate_GF(coords = rbind(loc_x, loc_y),
                  kappa=kappa_xi, variance=sigma2_xi)
length(xi)

# Now lets simulate the data for the covariate and the response
# To simulate the covaraite values (through the rpois function) the expected values Ei are required, which are generated using the Gamma distribution 
set.seed(134)
E <- rgamma(n=n_x, shape = 10, rate = 10) 
rho <- exp(xi[1:n_x])
set.seed(14)
x <- rpois(n=n_x, lambda=E*rho) #covariate values!

# To simulate the response data the intercept values has to be fixed (eual to 1)
b0 <- 1

# and the number of trials is randomly generated
set.seed(19)
ntrials <- 5 + rpois(n=n_y, lambda = 10)

# Finally to get the response values we use equation 8.4 and the antilogit transformation
eta_y <- b0 + #intercept
        beta1 + #coefficent
        xi[n_x + 1:n_y] + #this is just indexing the last 15 
        u
set.seed(553)
y <- rbinom(n=n_y, size=ntrials, prob=exp(eta_y)/(1 + exp(eta_y)))


##--Fitting the Model--##

# first step is to make mesh - this will be simple with only one outer extension
# The loc_x and loc_y are from the example above
mesh <- inla.mesh.2d(loc=rbind(loc_x,loc_y), max.edge = 0.15,
                     cutoff = 0.03, offset = 0.1)
plot(mesh)

# Now define one SPDE model which will be used for both spatial effects xi(s) and u(s)
spde <- inla.spde2.matern(mesh=mesh, alpha=2)


# Due to the misalignment, in order to link the locations with the mesh we have to
# define two projection matrices:
A_x <- inla.spde.make.A(mesh = mesh, loc=loc_x)
dim(A_x)
A_y <- inla.spde.make.A(mesh = mesh, loc=loc_y)
dim(A_y)

# Now need to create new data stack objects for each the covariate and the response
# First is for the covariate which will be modeled with a poisson distribution.
# NOTE in the `data` list I specify the type of link function, the expected # of cases (E) for the poisson likelihood
stk_x <- inla.stack(data = list(y=cbind(x, NA), E=E, link="log"),
                    effects = list(xi.field=1:spde$n.spde),
                    A=list(A_x),
                    tag="est.x")

stk_y <- inla.stack(data = list(y=cbind(NA,y), Ntrials = ntrials, link="logit"),
                    effects=list(
                        list(u.field=1:spde$n.spde, x.field=1:spde$n.spde), #index for u(s) and xi(s) spde
                        list(intercept=rep(1,length(y)))), #intercept (is there an intercept for each location?)
                    A=list(A_y,1), #why does this have a "1" but not in stk_x
                    tag="est.y"
                    )

# combine
stk <- inla.stack(stk_x, stk_y)


# The estimation procedure is carried out using the R-INLA "copy" feature which is adopted when a latent field is needed more than
# once in the model formulation.
# In this example, the spatial field xi(s) is used both in the speciication of the covariate and response linear predictor
# Need to define the "copy" feature in the formula

formula <- y ~ -1 + f(xi.field, model = spde) + # spatial field for covariate
                    intercept + #now the response part which had an intercept and...
                    f(u.field, model=spde) + #the u(s)....
                    f(x.field, copy="xi.field", fixed = FALSE, #and the xi.field which we are copying from above 
                      hyper=list(theta=list(param=c(-1,10))))

# To fit the model with the `inla` function we specify the poisson and the binomial families as well as the expected number of cases and of trials (which we can get from the stk)

output2 <- inla(formula, 
                data = inla.stack.data(stk),
                family = c("poisson","binomial"),
                E=inla.stack.data(stk)$E,
                Ntrials=inla.stack.data(stk)$Ntrials,
                control.predictor = list(compute=TRUE, A=inla.stack.A(stk)) # NOTE: Cam did it like this too
                )


##--Getting the results--## LOOK AT EXAMPLE ABOVE TO SEE HOW TO VIEW POSTERIOR SUMMARIES AND DISTRIBUTIONS FOR PARAMETERS








###########################################################################################
# Advanced Spatial Modeling with Stochastic Partial Differential Equations Using R and INLA
###########################################################################################

###---1.6  Advanced features---###

##--1.6.1 Several likeihoods--##

# This is a common approach to build joint models. By using a model with mode than one likelihood it is possible to build a 
# joint model with different types od outputs and the hyperparameteres in the likelihoods will be fitted separately

# Example using a toy dataset. This example will mimic a situation in which two sets of measurements are available but with different noise effects
data(SPDEtoy) 
SPDEtoy$y2 <- SPDEtoy$y + rnorm(nrow(SPDEtoy), sd = 2) #creating a new dataset

# Precision of the error terms between the two responses differ. Thus here we will fit a model with two gaussian likelihoods with different precisions to be estimated by each one of the sets of observations
# To fit a model with more than one likelihood the response variable must be a matrix with as many columns as the number of likelihoods. The number of rows in the total number of observations. The rows with data not associated with the likelihood are filled with NA values:A

n <- nrow(SPDEtoy) # number of locations

Y <- matrix(NA, ncol = 2, nrow = n *2)

# Add 'y' in the first column, rows 1 to 200
Y[1:n, 1] <- SPDEtoy$y

# Add 'y2' in the second column, rows 201 to 400
Y[n + 1:n, 2] <- SPDEtoy$y2

# The covaiates will be the same in both likelihoods! The 'family' agrument must be a vector with the names of the likelihoods used.
m0.2lik <- inla(Y ~ s1 + s2, family = c("gaussian","gaussian"),
                data = data.frame(Y = Y,
                                  s1 = rep(SPDEtoy$s1, 2), # why rep the responses twice I have no clue
                                  s2 = rep(SPDEtoy$s2, 2))
                )

summary(m0.2lik) # notice that  the summary of this model now shows estimates of two precisions. Also note how the precision of the second likelihood has a posterior mean smaller to the precision of the first likelihood. This is because the variance of the second set of observations is higher as they were generated by adding some noise to the original data.


##--1.6.2 Copy Model--##

# sometimes it is nessary to share an effect that is estimated from two or more parts of the dataset, so that all of them provide information about the effect when fitting the model. 
# This is the 'copy' effect as the new effect will be a copy of the original effect plus some tiny noise

# Here a new data.frame with the OG SPDEtoy data and the simulated one will be put together. This involves creates two new vectors of covariates by repeating the originial covariates. Then two indices will be created to ID to which group of observations a value belongs to:

y.vec <- c(SPDEtoy$y, SPDEtoy$y2) #responses combined
r <- rep(1:2, each = nrow(SPDEtoy)) #index to ID which group a value belongs to
s1.vec <- rep(SPDEtoy$s1, 2)
s2.vec <- rep(SPDEtoy$s2, 2)
i1 <- c(rep(1,n), rep(NA, n))
i2 <- c(rep(NA,n), rep(1,n))

d <- data.frame(y.vec, s1.vec, s2.vec, i1, i2)

# indices i1 and i2 have either 1 (i.e., the linear predictor includes the random effect) or NA (there is no random effect in the linear predictor). This ensures that the linear predictor includes the original random effect in the first 200 observations (using index i1), and the copied effect in the last 200 observations (using index i2)

tau.prior = list(prec = list(initial = 0.001, fixed = TRUE)) # This is to get similar results with previous examples (not shown). This is setting initial values for the precision of the random effects and fixing it in the prior definition. These values are passed to the prior defintion (in the call to the f() function below)

# The formula to fit the model is is defined using first an 'iid' mode with an index of i1 (first n values and NA for rest). The copied effect uses an index for the second set of observations and the values of covariates as weights.
f.copy <- y.vec ~ s1.vec + 
                    f(i1, s2.vec, model = "iid", hyper = tau.prior) +
                    f(i2, s2.vec, copy = "i1")

m0.copy <- inla(f.copy, data = d)
summary(m0.copy)


f.copy2 <- y.vec ~ s1.vec + f(i1, s2.vec, model = "iid") +
  f(i2, s2.vec, copy = "i1", fixed = FALSE)

m0.copy2 <- inla(f.copy2, data = d)
summary(m0.copy2)








###---3.1 Coregionalization Model---###

##--Data simulation--##

# Intercept on reparametrized model
alpha <- c(-5, 3, 10) 
# Random field marginal variances
m.var <- c(0.5, 0.4, 0.3) 
# GRF range parameters:
range <- c(4, 3, 2)
# Copy parameters: reparameterization of coregionalization 
# parameters
beta <- c(0.7, 0.5, -0.5) 
# Standard deviations of error terms
e.sd <- c(0.3, 0.2, 0.15)

# the number of observations of each response variable
n1 <- 99
n2 <- 100
n3 <- 101
 
# In this example there will be a different number of observations for each response variable and they will be observed at different locations
loc1 <- cbind(runif(n1)*10, runif(n1)*5)
loc2 <- cbind(runif(n2)*10, runif(n2)*5)
loc3 <- cbind(runif(n3)*10, runif(n3)*5)

# The `book.rMatern` function (practically same as simulate_GF as above) will be used to simulate independent random field realizations for each time.
# NOTE that we need, for example, the locations loc2 for the field u1 because this field is also used for the second variable.
book.rMatern <- function(n, coords, sigma=1, range, kappa = sqrt(8*nu)/range, variance = sigma^2, nu=1) {
    m <- as.matrix(dist(coords))
    m <- exp((1-nu)*log(2) + nu*log(kappa*m)-
             lgamma(nu))*besselK(m*kappa, nu)
    diag(m) <- 1
    return(drop(crossprod(chol(variance*m),
                          matrix(rnorm(nrow(coords)*n), ncol=n))))
}

set.seed(05101980)
z1 <- book.rMatern(1, rbind(loc1, loc2, loc3), range = range[1],
  sigma = sqrt(m.var[1]))
z2 <- book.rMatern(1, rbind(loc2, loc3), range = range[2],
  sigma = sqrt(m.var[2]))
z3 <- book.rMatern(1, loc3, range = range[3],
  sigma = sqrt(m.var[3]))

#Finally, we obtain samples from the observations:
set.seed(08011952)
y1 <- alpha[1] + z1[1:n1] + rnorm(n1, 0, e.sd[1])
y2 <- alpha[2] + beta[1] * z1[n1 + 1:n2] + z2[1:n2] + 
  rnorm(n2, 0, e.sd[2])
y3 <- alpha[3] + beta[2] * z1[n1 + n2 + 1:n3] + 
  beta[3] * z2[n2 + 1:n3] + z3 + rnorm(n3, 0, e.sd[3])


##--Model fitting--##

# This model only requires one mesh to fit all of the three spatial random field. Makes it easter to link the different effects across different outcomes at different spatial locations
mesh <- inla.mesh.2d(rbind(loc1,loc2,loc3),
                     max.edge = c(0.5, 1.25), offset = c(0.1,1.5), cutoff = 0.1)
plot(mesh)

# Next is defining the SPDE model -- they use the pcmatern but I probably won't (questions to ask!)
spde <- inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(0.5, 0.01), #P(range < 0.5) = 0.01
    prior.sigma = c(1,0.01) #P(sigma >1) = 0.01
)

# For each of the parameters (i.e., coefficients) of the copied effects, the prior is Gaussian with zero mean and precision 10.
hyper <- list(beta = list(prior = 'normal', param = c(0, 10))) # This is done in the above example as well but not sure if I need to (questions to ask!)


# make projection matrices for each set of locations
A1 <- inla.spde.make.A(mesh, loc1)
A2 <- inla.spde.make.A(mesh, loc2)
A3 <- inla.spde.make.A(mesh, loc3)

# Next lets make stacks for each and join them together
stack1 <- inla.stack(
    data = list(y = cbind(as.vector(y1), NA, NA)), #the response needs to be a vector
    A = list(A1),
    effects = list(
        list(intercept1 = 1, s1 = 1:spde$n.spde)
    )
)

stack2 <- inla.stack(
    data = list(y = cbind(NA, as.vector(y2), NA)),
    A = list(A2),
    effects = list(
        list(intercept2 = 1, s2 = 1:spde$n.spde, s12 = 1:spde$n.spde)
    )
)

stack3 <- inla.stack(
    data = list(y = cbind(NA, NA, as.vector(y3))),
    A = list(A3),
    effects = list(
        list(intercept3 = 1, s3 = 1:spde$n.spde, s13=1:spde$n.spde, s23=1:spde$n.spde)
    )
)

#combine
stack <- inla.stack(stack1, stack2, stack3)


# We use a PC prior for the error precision (see Section 1.6.5): -- don't worry about this either

hyper.eps <- list(hyper = list(prec = list(prior = 'pc.prec', 
                                           param = c(1, 0.01))))

# In this model there are 12 parameters in total: 2 hyperparameters for each pf the three spatial effects (i.e. range and sigma), 1 for each likelihood, and 3 copy parameters. To make the optimization process fast, the parameter values used in the simulation (plus some random noise) will be set as the initial values:

theta.ini <- c(log(1 / e.sd^2), 
  c(log(range), 
    log(sqrt(m.var)))[c(1, 4, 2, 5, 3, 6)], beta)
# We jitter the starting values to avoid artificially recovering 
# the true values
theta.ini = theta.ini + rnorm(length(theta.ini), 0, 0.1)

# The formula including all the terms in the model is defined as follows:

form <- y ~ 0 + intercept1 + intercept2 + intercept3 + 
  f(s1, model = spde) + f(s2, model = spde) + 
  f(s3, model = spde) + 
  f(s12, copy = "s1", fixed = FALSE, hyper = hyper) + 
  f(s13, copy = "s1", fixed = FALSE, hyper = hyper) + 
  f(s23, copy = "s2", fixed = FALSE, hyper = hyper) 


# Given that this model is complex and may take a long time to run, the empirical Bayes approach will be used, by setting int.strategy = 'eb' below, instead of integrating over the space of hyperparameters. The model is fitted with the following R code:

result <- inla(form, rep('gaussian', 3), 
  data = inla.stack.data(stack), 
  control.family = list(hyper.eps, hyper.eps, hyper.eps), 
  control.predictor = list(A = inla.stack.A(stack)),
  control.mode = list(theta = theta.ini, restart = TRUE),
  control.inla = list(int.strategy = 'eb'))


