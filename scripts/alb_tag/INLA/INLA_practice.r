# practicing INLA: A toy example

library(INLA)
data(SPDEtoy)

# this dataset has 3 columns
# first 2 are coordinates and the response is the third
str(SPDEtoy)

######################
# Making the INLA mesh
######################

# Here a domain is first defined to create the mesh
# (just named mesh5 because it was the 5th made in their example)
pl.dom <- cbind(c(0, 1, 1, 0.7, 0), c(0, 0, 0.7, 1, 1))
mesh5 <- inla.mesh.2d(loc.domain = pl.dom, 
                      max.edge = c(0.092,0.2)) 

###############################
# Parameterizing the SPDE model
###############################
spde5 <- inla.spde2.matern(mesh = mesh5,
                           alpha = 2 # smoothness parameter
                           )

##################
# Projector matrix
##################

# contains the basis function value for each basis
# with one basis function at each column
# This will be used to interpolate the random field
# that is being modeled at the mesh nodes

# Rows are obs and columns are vertices/nodes.
# The values in each row are how close the obs were to those vertices.
# There can be no more than 3 non-zero cells filled for each row and thus the matrix is sparse
# which makes is faster for computation

coords <- as.matrix(SPDEtoy[,1:2])
A5 <- inla.spde.make.A(mesh5, loc = coords)

# rows are obs and columns are vertices (row x vertices)
dim(A5) # 200 x 490

# For each obs the basis functiosn sums to 1
table(rowSums(A5))

################
# The Data Stack
################

# `inla.stack` is useful for organziing data, covariates, indices, and projector matrices.
# Helps control the way effects are projected in the linear predictor

# This function allows us to work with predictors that includes terms with different dimentions
# 3 main agruments are a vector list with the data, list of projector matrices, and list of effects, and a label (i.e. tag)

stk5 <- inla.stack(
    data = list(resp = SPDEtoy$y),
    A = list(A5, 1), # list of the projector matrix and a 1 (need to add the 1)
    effects = list(i = 1:spde5$n.spde, # This is another index. This can also be defined with inla.spde.make.index
                   beta0 = rep(1, nrow(SPDEtoy))),
    tag = 'est'
)

#or 

s_index<-inla.spde.make.index(name = "spatial.field",
                     n.spde = spde5$n.spde)

stk5 <- inla.stack(
    data = list(resp = SPDEtoy$y),
    A = list(A5, 1), # list of the projector matrix and a 1 (need to add the 1)
    effects = list(i = s_index, # This is another index. This can also be defined with inla.spde.make.index
                   beta0 = rep(1, nrow(SPDEtoy))),
    tag = 'est'
)
stk5 %>% str()
# The inla.stack() function automatically eliminates any column in a projector matrix that has a zero sum,
# and it generates a new and simplified matrix

###########################
# Model Fitting and Results
###########################

# To fit the model, the intercept in the formula must be removed and added as a covariate term in the list of
#  effects, so that all the covariate terms in the formula can be included in a projector matrix.
#  Then, the matrix of predictors is passed to the inla() function in its control.predictor argument, as follows:

res5 <- inla(resp ~ 0 + beta0 + f(i, model = spde5), #intercept (beta0) and the spatial random field (i)
             data = inla.stack.data(stk5),
             control.predictor = list(A = inla.stack.A(stk5)))

# inla function returns a object with several results:
# summaries, marginal posterior densities for each paramter,
# regression parameters, latent fields and all hyperparamters

# summary of the intercept (beta0):
res5$summary.fixed

# Similarly, the summary of the precision of the Gaussian likelihood
res5$summary.hyperpar[1,]

# *more on looking at the marginal distributions*

################################
# Projection of the random field
################################

# A common objective when dealing with spatial data collected at some locations is
# the prediction on a fine grid of the spatial model to get high resolution maps. 
# In this section we show how to do this prediction for the random field term only.

#Here we will predict at 3 locations: (0.1, 0.1), (0.5, 0.55), (0.7, 0.9)
pts3 <- rbind(c(0.1, 0.1), c(0.5, 0.55), c(0.7, 0.9))

A5pts3 <- inla.spde.make.A(mesh5, loc = pts3)
dim(A5pts3)

#to visualize only the columns with non-zero elements of this matrix:
jj3 <- which(colSums(A5pts3) > 0)
A5pts3[,jj3]

#this projector matrix can the be used to interpolating a function of the random field, for example the posterior mean
drop(A5pts3 %*% res5$summary.random$i$mean)

## prediction of the random field on a grid (mappping random field on a fine grid)
# we can use the inla.mesh.projector function to compute a projector matrix
# automatically for a grid of points over a squaire that contains the mesh
pgrid0 <- inla.mesh.projector(mesh5, xlim = 0:1, ylim = 0:1, 
                              dims = c(101,101) #101 x and y's
                              )

#posterior mean and sd for the random field in the grid
prd0.m <- inla.mesh.project(pgrid0,  res5$summary.random$i$mean)
prd0.s <- inla.mesh.project(pgrid0,  res5$summary.random$i$sd)
inla.mesh.project(pgrid0,  res5$summary.random$i$sd)
res5$summary.fixed$mean %>% head()

############
# Prediction
############

# Another quantity of interest when modeling spatially continuous processes is the 
# prediction of the expected value on a target location for which data have not been observed.
# aka computing the marginal distribution of the expected value at the target location

# have to into it in a stack like when model fitting but now
# all the fixed effects have to considered in the predictor and effects slots of the stack

stk5.pmu <- inla.stack(
  data = list(resp = NA), #set to NA to indicate that prediction should be carried out
  A = list(A5pts3, 1), 
  effects = list(i = 1:spde5$n.spde, beta0 = rep(1, 3)), 
  tag = 'prd5.mu')

# this stack is then joint to the data stack to fit the model again?
# I guess we can save time by considering the previous fitted model parameters here using `control.mode`
stk5.full <- inla.stack(stk5, stk5.pmu)

r5pmu <- inla(resp ~ 0 + beta0 + f(i, model = spde5), # intercept and random field 
  data = inla.stack.data(stk5.full), #stack for fitting the model and prediction
  control.mode = list(theta = res5$mode$theta, restart = FALSE), #just saying to use the same as before and not refit to save time
  control.predictor = list(A = inla.stack.A( 
    stk5.full), compute = TRUE))

# The fitted values for all the obs in the dataset are sumamrized in a data.frame
# To find the predicted values for the values with missing observations, the index to their rows in the data.frame must be found first
indd3r <- inla.stack.index(stk5.full, 'prd5.mu')$data
indd3r

# To get the summary of the posterior distributions of the marginal distribution at the target location's the
# index must be pased to the data.frame with the summary stats
r5pmu$summary.fitted.values[indd3r, c(1:3, 5)]

#Also, the posterior marginal distribution of the predicted values can be obtained:
marg3r <- r5pmu$marginals.fitted.values[indd3r]

# Finally, a 95% HPD interval for μ at the second target location can be computed with the function inla.hpdmarginal():
inla.hpdmarginal(0.95, marg3r[[2]])

######################
# Prediction on a grid
######################




###################################################################################################################
#                                  Doing the INLA turtorial through STSR book
###################################################################################################################
library("INLA")
library("dplyr")
library("tidyr")
library("ggplot2")
library("STRbook")
data("MOcarolinawren_long", package = "STRbook")

coords <- unique(MOcarolinawren_long[c("loc.ID", "lon", "lat")])

# non - convex hull method
boundary <- inla.nonconvex.hull(as.matrix(coords[, 2:3]))

# making the mesh
MOmesh <- inla.mesh.2d(boundary = boundary,
                       max.edge = c(0.8, 1.2), # max. edge length
                       cutoff = 0.1) # min. edge length

plot(MOmesh, asp = 1, main = "")
points(coords[c("lon", "lat")], col = "red", type = "p")

# making the SPDE model object
# NOTE: book version does use pcmatern to put priors on the range and sigma parameters but not going to do that here
spde <- inla.spde2.matern(mesh = MOmesh,
                            alpha = 2)

# For now lets try a spatial model ---- NOT spatial-temporal
n_years <- length(unique(MOcarolinawren_long$t))
n_spatial <- MOmesh$n
s_index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = n_spatial,
                                n.group = n_years)

# the s_index is a list with one vector (spatial.field) being a index from 1:n_spatial
# the second vector (spatial.field.group) being what group it belongs. So if it were spatial temporal 
# n.group agrument would be a vector of which year the observation belongs, and the third vector is
# just a vector of 1's repeating n_spatial * n.group times (if there is no n.group then it is just by n_spatial)

# making the A matrix
# again only doing a spatial model but if it were spatiotemporal there would also be a group and n.group agrument
coords.allyear <- MOcarolinawren_long[c("lon","lat")] %>% 
                    as.matrix()

PHI <- inla.spde.make.A(mesh = MOmesh,
                        loc = coords.allyear,
                        group = MOcarolinawren_long$t,
                        n.group = n_years)

dim(PHI) #1575 obs x 259 vertices ---- if it were spatiotemporal the 259 would be 5439 because there would be 259 vertices * number of years (21)

# making an INLA stack but both estimation (model fitting) and prediction
## First state: estimation
n_data <- nrow(MOcarolinawren_long)

stack_est <- inla.stack(
    data = list(cnt = MOcarolinawren_long$cnt),
    A = list(PHI, 1), #the 1 indicates the fixed effects are mapped one-to-one to the response
    effects = list(s_index,
                   list(Intercept = rep(1,n_data))),
    tag = "est"
)

## now lets construct a stack containing the matrices and vectors for defining the model at prediction locations
# here we will chose the vertices as the prediction locations
df_pred <- data.frame(lon = MOmesh$loc[,1],
                      lat = MOmesh$loc[,2]
                      )
n_pred <- nrow(df_pred)
PHI_pred <- Diagonal(n = n_pred) #this is an identify matrix

# predction stack is really similar to estimation but now data values are set to 'NA'
# to indicate that prediction should be carried out to these locations
## Second stack: Prediction
stack_pred <- inla.stack(
    data = list(cnt = NA),
    A = list(PHI_pred, 1),
    effects = list(s_index,
                   list(Intercept = rep(1, n_pred))),
    tag = "pred")


# The estimation and precition stack are combined
stack <- inla.stack(stack_est, stack_pred)

# Now all that remians before fitting the model is to define the formula 
# which is a combination of a standard R formula for the fixed effects and 
# an INLA formula for the spatial (or spatiotemporal) residual component. For the
# the spatial part we need to specify the name of the index as the first agrument "spatial.field" (s_index$spatial.field),
# the model (spde), the name of the grouping/time index (we do not need to do this because we are doing just a spatial model),
# and finally some model if you want it to be constructed across groups (e.g. ar1)

## Formula
formula <- cnt ~ -1 + Intercept +
                f(spatial.field,
                  model = spde)


# Finally we have everything to run the main function (inla) for fitting the model
# This needs the data from the stack (extracted through inla.stack.data) and the 
# family (in this case negative-binomial). Here inla is instructed to fit the model
# as well as compute the predictions at the required locations

# looks like inla.stack.data and inla.stack.A are both required in order to extract the
# data and needed for the control.predictor agrument 
output <- inla(formula,
               data = inla.stack.data(stack, spde = spde),
               family = "nbinomial",
               control.predictor = list(A = inla.stack.A(stack),
                                        compute = TRUE)
               )


# INlA provides approximate marginal posterior distributions for the spatial field (or for spatial fields if you had time)
# and the parameters (coefficients, 
#                     ρ = precision hyperparameter...how "precise" your measurements are in the sense of having larger or smaller error (dont think we have it here since we didnt include ar1), 
#                     τ = variance, 
#                     κ = scale i.e. spatial correlation length).

# From the posterior distributions over the precision parameter
# τ and scale parameter κ, we can readily obtain marginal posterior distributions over
# the more interpretable variance parameter σ2 and practical range parameter l.


output.field <- inla.spde2.result(inla = output,
                                  name = "spatial.field",
                                  spde = spde,
                                  do.transf = TRUE #I think this just means transform the data
                                  )

## plot p(beta0 | Z)
plot(output$marginals.fix$Intercept,
     type = 'l',
     xlab = expression(beta[0]),
     ylab = expression(p(beta[0]*"|"*Z)))

## plot p(rho | Z) #NOTE: NOT INCLUDED IN THIS MODEL SINCE ITS NOT AR1
plot(output$marginals.hyperpar$`GroupRho for spatial.field`,
     type = 'l',
     xlab = expression(rho),
     ylab = expression(p(rho*"|"*Z)))

## plot p(sigma^2 | Z)
plot(output.field$marginals.variance.nominal[[1]],
     type = 'l',
     xlab = expression(sigma^2),
     ylab = expression(p(sigma^2*"|"*Z)))

## plot p(range | Z)
plot(output.field$marginals.range.nominal[[1]],
     type = 'l',
     xlab = expression(l),
     ylab = expression(p(l*"|"*Z)))

## We provide the prediction (posterior mean) and prediction standard error (posterior standard deviation)
index_pred <- inla.stack.index(stack, "pred")$data
lp_mean <- output$summary.fitted.values$mean[index_pred] #mean
lp_sd <- output$summary.fitted.values$sd[index_pred] #sd

# Next, we need to create a spatial grid upon which we map the predictions and their associated prediction standard errors
grid_locs <- expand.grid(
    lon = seq(min(MOcarolinawren_long$lon) - 0.2,
              max(MOcarolinawren_long$lon) + 0.2,
              length.out = 80),
    lat = seq(min(MOcarolinawren_long$lat) - 0.2,
              max(MOcarolinawren_long$lat) + 0.2,
              length.out = 80))

# The function inla.mesh.projector provides all the information required, based on the created spatial grid, to carry out the mapping.
proj.grid <- inla.mesh.projector(MOmesh,
                                 xlim = c(min(MOcarolinawren_long$lon) - 0.2,
max(MOcarolinawren_long$lon) + 0.2),
                                 ylim = c(min(MOcarolinawren_long$lat) - 0.2,
max(MOcarolinawren_long$lat) + 0.2),
                                 dims = c(80, 80))


pred <- inla.mesh.project(proj.grid,  lp_mean)
sd <- inla.mesh.project(proj.grid,  lp_sd)


x<-pred %>% as.data.frame() %>% gather(loc,count) %>% cbind(grid_locs) %>% dplyr::select(-loc)

ggplot()+
geom_raster(x, mapping=aes(x = lon, y = lat,  fill = count))+
viridis::scale_fill_viridis()









###########################################################################################
# Advanced Spatial Modeling with Stochastic Partial Differential Equations Using R and INLA
###########################################################################################

###---Chapter 4 Point processess and preferential sampling---###

# point pattern records the occurence of events in a study region - these locations of the observed events depend on an underlying spatial process which is often modeled using an intensity function (lambda(s))
# lambda(s) - measures the average # of events per unit of space and it can be modeled to depend on covariates and other effects

# Modeling this point process is known as a log-gaussian Cox process (LGCP). 
# A Cox process is just a name for a Poisson process with vary intensity (think about dividing the study region into cells and counting the number of points in each one and modeling that)

# LGCP definition - The Cox process is a Poisson process with intensity (lambda(s)) that varies over space. Given some area (A - think about a grid cell), the probability of observing a certain number of points in that area follows a Poisson distribution with intensity (expected value)


###---Data Simulation---###

library(spatstat) # using rLGCP() function to sample from a log-cox point process
win <- owin(c(0,3),c(0,3)) # use a (0,3) x (0,3) simulation window
win


npix <- 300 # resolution of the grid
spatstat.options(npixel = npix)

# Modeling intensity as: log( lambda(s) ) = beta0 + S(s)
# beta0 is a fixed value and S(s) is a Gaussian spatial process with Matérn covariance and zero mean
# Parameter beta0 can be regarded as a global mean level for the log intensity; i.e. the log-intensity fluctuates about it according to the spatial process S(s)
# IF NO SPATIAL FIELD, the expected number of points is e^beta0 times the are of the window. Thus the expected number of points is:
beta0 <- 3
exp(beta0) * diff(range(win$xrange)) * diff(range(win$yrange)) #expected number


# We will use a Matern covariance function with nu = 1 (not sure what nu is!). The other parameters are the variance (sigma2x) and the scale (range)
# (look at the text to see why these values were chosen - TLDR will produce a smooth changes)
sigma2x <- 0.2 
range <- 1.2
nu <- 1

# Time to simulate points
library(RandomFields)

set.seed(1)
lg.s <- rLGCP('matern', beta0, var = sigma2x,
              scale = range / sqrt(8), nu = nu, win = win)

# The coordinates of the observed events of the point pattern can be obtained
xy <- cbind(lg.s$x, lg.s$y)[, 2:1]

n <- nrow(xy)

###---Inference---###

# The augmented data set is made of a binary response, with 1 for the observed points and 0 for some dummy observations. Both the observed and dummy observations will have associated ‘expected values’ or weights that will be included in the Poisson regression

# Building the mesh
# For point process we dont ususally use the locations points as mesh nodes. We will use the loc.docmain agrument to build the mesh.
loc.d <- 3 * cbind(c(0,1,1,0,0), c(0,0,1,1,0))
mesh <- inla.mesh.2d(loc.domain = loc.d, offset = c(0.3, 1),
                     max.edge = c(0.3, 0.7), cutoff = 0.05)
plot(mesh)

nv <- mesh$n

# spde model - here they set priors
spde <- inla.spde2.pcmatern(mesh = mesh,
  # PC-prior on range: P(practic.range < 0.05) = 0.01
  prior.range = c(0.05, 0.01),
  # PC-prior on sigma: P(sigma > 1) = 0.01
  prior.sigma = c(1, 0.01)) 

# create integration stack
# The SPDE approach for point pattern analysis defines the model at the nodes of the mesh. To fit the log-Cox point process model, these points are considered as integration points.
# The expected number of events is proportional to the area around the node (the areas of the polygons in the dual mesh). This means that at the nodes of the mesh with larger triangles, there are also larger expected values.

dmesh <- book.mesh.dual(mesh) # making a dual mesh

# converting the polygons into a `SpatialPolygons`
domain.polys <- Polygons(list(Polygon(loc.d)), '0')
domainSP <- SpatialPolygons(list(domain.polys))

# since the mesh is larger than the study area, we need to compute the interaction between each polygon in the dual mesh and the study area
library(rgeos)
w <- sapply(1:length(dmesh), function(i) {
  if (gIntersects(dmesh[i, ], domainSP))
    return(gArea(gIntersection(dmesh[i, ], domainSP)))
  else return(0)
})

# sum of the weights is the area of the study region
sum(w)

# checking how many integration points have zero weight
table(w > 0)

##--Data and projection matrices--##

# The vector of weights computed above (w) is what we need to used as the exposure (E) in the poisson likelihood in INLA
# We will make another vector of 1's (representing the observations) with a sequence of 0's (representing the mesh nodes)
y.pp <- rep(0:1, c(nv, n)) #nv was the number of nodes in mesh, n is the number of locations

# The exposure vecotr can be defined as:
e.pp <- c(w, rep(0, n)) #makes a vector that combines the weights with a sequence of 0's (the number of observations)

# Next we create the projection matrix. This takes two steps:
# 1) For the integration points this is just a diagonal because these locations are just the mesh vertices
imat <- Diagonal(nv, rep(1, nv))

# 2) For the observed points, another projection matrix is defined:
lmat <- inla.spde.make.A(mesh, xy)

# Then the entire projection matrix is:
A.pp <- rbind(imat, lmat)

# Now we can set up the data stack
stk.pp <- inla.stack(
  data = list(y = y.pp, e = e.pp),
  A = list(1, A.pp), # the 1 is for the intercept (b0). Since I put that first in A, the intercept needs to be first in the effects
  effects = list(
    list(b0 = rep(1, nv + n)),
    list(i = 1:spde$n.spde)
    ),
  tag = "pp"
)


##--Posterior Marginals--##

# Here we will fit the model to get the posterior marginals for all the parameters
pp.res <- inla(y ~ 0 + b0 + f(i, model = spde), # looks like putting 0 produces practically the same result as -1
               family = "poisson", 
               data = inla.stack.data(stk.pp),
               control.predictor = list(#compute = TRUE,
                                        A=inla.stack.A(stk.pp)),
               E = inla.stack.data(stk.pp)$e
               )
pp.res$summary.hyperpar



###---Chapter 4.2 Including a covariate in the log-Gaussian Cox process---###

# When adding a covariate to the linear predictor in a LGCP process we must know the values of the covariate at the location AND integration points

###---Simulation of the covaraite---###

# First need to define a covariate everywhere int he study area. 
# This will be computed at a grid using the `spatstat` package. 
# Locations simulated will cover both the study window and mesh points
library(spatstat)


# Use expanded range
x0 <- seq(min(mesh$loc[,1]), max(mesh$loc[,1]), length = npix) #using mesh and npix from above
y0 <- seq(min(mesh$loc[,2]), max(mesh$loc[,2]), length = npix)
gridcov <- outer(x0, y0, function(x,y) cos(x) - sin(y-2)) # covariate values at grid

# Now, the expected number of points is a function of the covariate
beta1 <- -0.5
sum(exp(beta0 + beta1 * gridcov) * diff(x0[1:2]) * diff(y0[1:2])) #beta0 from above

# simulate the point pattern
set.seed(1)
lg.s.c <- rLGCP('matern', im(beta0 + beta1 * gridcov, 
                             xcol = x0, yrow = y0), 
                var = sigma2x, scale = range / sqrt(8), 
                nu = 1, win = win)

# Both the spatial field and the point pattern are returned. The locations are:
xy.c <- cbind(lg.s.c$x, lg.s.c$y)[, 2:1]
n.c <- nrow(xy.c) # number of point pattern locations


###---Inference---###

# The covaraite values must be available at the point pattern locations and at the mesh nodes

# The covariate values will be collected by interpolating from the grid data. Need to make an `im` object and then we use the function `interp.im` 
covariate.im <- im(gridcov, x0, y0)
covariate.im
covariate <- interp.im(covariate.im,
                       x = c(mesh$loc[,1], xy.c[,1]), # the mesh and the point process lon positions
                       y = c(mesh$loc[, 2], xy.c[,2])) # the mesh and point process lat poistions

# The augmented data set is made of a binary response, with 1 for the observed points and 0 for the nodes.

y.pp.c <- rep(0:1, c(nv, n.c)) 
e.pp.c <- c(w, rep(0, n.c))

#  projection matrices for the integration points and the observation points
 
# diagonal matrix for integration point A matrix
imat <- Diagonal(nv, rep(1, nv))

# Projection matrix for the observed locations
lmat.c <- inla.spde.make.A(mesh, xy.c)
# The two projection matrics can be merged with rbind() for a total projection matrix

A.pp.c <- rbind(imat, lmat.c)


# Data stack is practically the same as above but it includes the covariates
stk.pp.c <- inla.stack(
  data = list(y = y.pp.c, e = e.pp.c),
  A = list(1 , A.pp.c),
  effects = list(
    list(b0 = 1,
         covariate = covariate),
    list(i = 1:nv)
  ),
  tag = 'pp.c'
)


# fitting the model
pp.c.res <- inla(y ~ 0 + b0 + covariate + f(i, model = spde),
                 family = 'poisson', 
                 data = inla.stack.data(stk.pp.c),
                 control.predictor = list(
                   A = inla.stack.A(stk.pp.c)),
                 E = inla.stack.data(stk.pp.c)$e)

# summary of the model hyperparameters
pp.c.res$summary.hyperpar # results differ from example














#########################
# BOOK DUAL MESH FUNCTION
#########################
book.mesh.dual <- function(mesh) {
    if (mesh$manifold=='R2') {
        ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
            colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
        library(parallel)
        pls <- mclapply(1:mesh$n, function(i) {
            p <- unique(Reduce('rbind', lapply(1:3, function(k) {
                j <- which(mesh$graph$tv[,k]==i)
                if (length(j)>0) 
                    return(rbind(ce[j, , drop=FALSE],
                                 cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                                       mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                                       mesh$loc[mesh$graph$tv[j, k], 2] +
                                       mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
                else return(ce[j, , drop=FALSE])
            })))
            j1 <- which(mesh$segm$bnd$idx[,1]==i)
            j2 <- which(mesh$segm$bnd$idx[,2]==i)
            if ((length(j1)>0) | (length(j2)>0)) {
                p <- unique(rbind(mesh$loc[i, 1:2], p,
                                  mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                                  mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                                  mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                                  mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
                yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
                xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
            }
            else {
                yy <- p[,2]-mesh$loc[i, 2]
                xx <- p[,1]-mesh$loc[i, 1]
            }
            Polygon(p[order(atan2(yy,xx)), ])
        })
        return(SpatialPolygons(lapply(1:mesh$n, function(i)
            Polygons(list(pls[[i]]), i))))
    }
    else stop("It only works for R2!")
}
