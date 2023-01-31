For references please see 
1) Illian, J.B., Martino, S., S�rbye, S.H., Gallego-Fern�ndez, J.B., Zunzunegui, M.,
Paz Esquivias, M. & Travis, J.M.J. (2013) Fitting complex ecological point process 
models with integrated nested Laplace approximation. Methods in Ecology and 
Evolution, 4(4),305-315.
2) Illian, J.B., S�rbye, S., Lindgren, F. & Rue, H.: workshop "spatial modelling with INLA"
3) Krainski, E.T., Lindgren, F., Simpson, D. & Rue, H. (2015) The R-INLA tutorial on 
SPDE models. 
4) Blangiardo, M., Cameletti, M., Baio, G. & Rue, H. (2013), Spatial and spatio-temporal models 
with R-INLA., Spatial and Spatio-temporal Epidemiology, 4, 33-49.

#APPENDIX 2

################################################################################################################
#Spatial single-species models (for data without excess zeros)
#and model selection
#example: Harbour porpoise density maps 
################################################################################################################

#Loading libraries
library(INLA)

#All the IGMRF-models are scaled automatically
inla.setOption(scale.model.default = TRUE)

#Loading the data file
data0 <- readRDS("data0.rds")

#Creating all possible models...
ColNam <- c("BT","CHL","NPP","PEA", "SP","DVV") 
n0 <- length(ColNam) 
id1 <- unlist(lapply(1:n0,function(x) combn(1:n0,x,simplify=F)),recursive=F) 

#...excluding the models with the highly correlated variables
ind <- rep(NA,length(id1))
for (i in 1:length(id1)) ind[i] <- any(grepl(1,id1[[i]])) & any(grepl(3,id1[[i]])) | 
  any(grepl(1,id1[[i]])) & any(grepl(4,id1[[i]])) | 
  any(grepl(3,id1[[i]])) & any(grepl(4,id1[[i]]))
id1 <- id1[-c(which(ind))]

#Besag-York-Mollie (BYM) models with default priors
#where the covariates' effects are modelled as second-order random walk processes (RW2)
all.models <- lapply(id1,function(x)paste(paste("density~1+",
  paste("f(inla.group(",ColNam[x],"),model='rw2')",collapse="+"),
  paste("+f(IDgrid,model='bym',graph='graph.dat')")))) 

#Run the models
DIC<-rep(NA,length(all.models))
for (i in 1:length(all.models)){
  formula <- as.formula(all.models[[i]])
  result <- inla(formula,family='gamma',
               data=data0, control.compute=list(dic=TRUE),
               control.predictor=list(compute=TRUE),
               verbose=TRUE
  )
  DIC[i] <- result$dic$dic
}

#The model with the lowest DIC values was considered as the best models. 
best.model <- all.models[[which.min(DIC)]]

################################################################################################################
#Spatial hurdle single-species models
#using the stochastic partial differential equation (SPDE) approach
#example: harbour seal usage density maps
################################################################################################################

#Loading the data file
data02 <- readRDS("data02.rds")
density <- data02$density

#define the two response variables
z <- (density>0) + 0
y <- ifelse(density>0, density, NA)

#Loading the boundary
bound <- readRDS("bound.rds")

#Creating a triangle mesh
(mesh <- inla.mesh.2d(loc=data02[,1:2], boundary=bound, max.edge=1.5, cut=0.1))$n

#Coordinates
coo <- as.matrix(data02[,1:2])

#SPDE model
spde <- inla.spde2.matern(mesh=mesh, alpha=2)

#Projector matrix for the SPDE model 
A <- inla.spde.make.A(mesh=mesh, loc=coo)

#The data stack (combining data)
#for the occurrence
stk.z <- inla.stack(tag='est.z', data=list(z=z, y=cbind(z, NA)), A=list(A, 1),
                    effects=list(list(i.z=1:spde$n.spde),
                                 list(z.b0=rep(1,length(z)),
                                      NPP=data02$NPP)))

#for the density
stk.y <- inla.stack(tag='est.y', data=list(r=y, y=cbind(NA, y)), A=list(A, 1),
                    effects=list(list(i.y=1:spde$n.spde),
                                 list(y.b0=rep(1, length(y)),
                                      BT=data02$BT,
                                      SP=data02$SP)))


#Fitting the model for each response separately
res.z <- inla(z ~ 0 + z.b0 + 
  f(inla.group(NPP), model='rw2',hyper = hyper.c) + 
  f(i.z, model=spde, hyper=hyper.spat), family='binomial',
  data=inla.stack.data(stk.z), control.compute=list(dic=TRUE),
  control.predictor=list(A=inla.stack.A(stk.z), compute=TRUE))

res.y <- inla(r ~ 0 + y.b0 + 
  f(inla.group(BT), model='rw2',hyper = hyper.c) +
  f(inla.group(SP), model='rw2',hyper = hyper.c) +
  f(i.y, model=spde, hyper=hyper.spat), family='gamma',
  data=inla.stack.data(stk.y), control.compute=list(dic=TRUE),
  control.predictor=list(A=inla.stack.A(stk.y), compute=TRUE))

#Defining a full data stack for the hurdle model
stk.zy <- inla.stack(stk.z, stk.y)

#Fitting the hurdle model
result.h <- inla(y ~ 0 + z.b0 + y.b0 +
  f(inla.group(NPP), model='rw2',hyper = hyper.c) + 
  f(inla.group(BT), model='rw2',hyper = hyper.c) +
  f(inla.group(SP), model='rw2',hyper = hyper.c) +
  f(i.z, model=spde, hyper=hyper.spat) + 
  f(i.y, copy='i.z', fixed=FALSE, hyper=hyper.spat),
  family=c('binomial', 'gamma'),
  data=inla.stack.data(stk.zy), control.compute=list(dic=TRUE),
  control.predictor=list(A=inla.stack.A(stk.zy), compute=TRUE))

################################################################################################################
#Spatiotemporal zero-inflated single-species models
#using the stochastic partial differential equation (SPDE) approach
#example: Black-legged kittiwake, observations
################################################################################################################

#Loading the data file
data03 <- readRDS("data03.rds")

#Creating a triangle mesh
mesh <- inla.mesh.2d(loc=data03[,1:2], boundary=bound, max.edge=1.5, cut=0.1)

#Coordinates
coo <- as.matrix(data03[,1:2])

#SPDE model
spde <- inla.spde2.matern(mesh=mesh, alpha=2)

#Projector matrix for the SPDE model 
A <- inla.spde.make.A(mesh=mesh, loc=coo)

#The data stack
stk <- inla.stack(tag='est', data=list(y=y), A=list(A, 1),
                    effects=list(list(i.y=1:spde$n.spde),
                                 list(b0=1, year = data03$year,
                                      effort = data03$effort,
                                      NPP=data03$NPP,
                                      PEA=data03$PEA,
                                      DVV=data03$DVV)))

#Fitting the zero-inflated Poisson model with default priors
#See "Prior choice" section about defining and comparing different priors
result.zip <- inla(y ~ 0 + b0 + f(year, model='rw1') + f(effort, model='rw2') +
                     f(inla.group(NPP), model='rw2') + 
                     f(inla.group(PEA), model='rw1') +
                     f(inla.group(DVV), model='rw2') +
                     f(i.y, model=spde), 
                   family="zeroinflatedpoisson1",E=E,
                   data=inla.stack.data(stk), control.compute=list(dic=TRUE),
                   control.predictor=list(A=inla.stack.A(stk), compute=TRUE))


################################################################################################################
#Joint spatial models
#using Besag-York-Mollie (BYM) specification
#example: Harbour porpoise density map (2005) and Herring abundance (ages 1 and 2&3) (across years)
#herring of different ages were regarded as separate species 
################################################################################################################

#Loading the data file
species1 <- readRDS("species1.rds")#Porpoise2005
species2 <- readRDS("species2.rds")#Herring (age 1)
species3 <- readRDS("species3.rds")#Herring (ages 2&3)

n1<-dim(species1)[1]
n2<-dim(species2)[1]
n3<-dim(species3)[1]

nothing1 = rep(NA, n1)
nothing2 = rep(NA, n2)
nothing3 = rep(NA, n3)

#Generating 3 vectors of response variables for the joint model 
y = as.vector(species1$density)
yNA = as.vector(c(y, nothing2, nothing3))

z = as.vector(species2$density)
zNA = as.vector(c(nothing1, z, nothing3))

w = as.vector(species3$density)
wNA = as.vector(c(nothing1, nothing2, w))

#Combine them in a matrix
outcome.matrix<-matrix(c(yNA,zNA, wNA), ncol=3)

#A factor vector with 3 levels indicating 3 outcome variables
mu = as.factor(c(rep(1,length(y)), rep(2,length(z)), rep(3,length(w))))

#Index vectors for the spatial effects 
i.spat1 = c(species1$IDgrid, nothing2, nothing3)
i.spat2 = c(nothing1, species2$IDgrid, nothing3)
i.spat3 = c(nothing1, nothing2, species3$IDgrid)

#The covariates 
BT1 <-c(species1$BT, nothing2, nothing3)

NPP2 <-c(nothing1, species2$NPP, nothing3)
SP2 <-c(nothing1, species2$SP, nothing3)
DVV2 <-c(nothing1, species2$DVV, nothing3)

CHL3 <-c(nothing1, nothing2, species3$CHL)
NPP3<-c(nothing1, nothing2, species3$NPP)
DVV3 <-c(nothing1, nothing2, species3$DVV)

#Data set
data=list(outcome.matrix=outcome.matrix, i.spat1=i.spat1, i.spat2=i.spat2, i.spat3=i.spat3,
          BT1=BT1, NPP2=NPP2, SP2=SP2, DVV2=DVV2, CHL3=CHL3, NPP3=NPP3, DVV3=DVV3
)

formula = outcome.matrix ~ mu -1 +
  f(inla.group(BT1), model="rw2", hyper=param.cc01) +
  f(inla.group(NPP2), model="rw2", hyper=param.cc01) +
  f(inla.group(SP2), model="rw2", hyper=param.cc01) +
  f(inla.group(DVV2), model="rw2", hyper=param.cc01) +
  f(inla.group(CHL3), model="rw2", hyper=param.cc01) +
  f(inla.group(NPP3), model="rw2", hyper=param.cc01) +
  f(inla.group(DVV3), model="rw2", hyper=param.cc01) +
  
  f(i.spat1, model="bym", graph = "graph2.dat", hyper= param.spat01) +
  f(i.spat2, copy="i.spat1", fixed=F) +
  f(i.spat3, copy="i.spat1", fixed=F)

result.j1 <- inla(formula,family=c("gamma","gamma","gamma"), data=data, 
               control.compute=list(dic=T), 
               control.predictor = list(compute = TRUE), verbose=TRUE,
               control.fixed = list(expand.factor.strategy = "inla")) 

################################################################################################################
#Joint hurdle spatial models
#using Besag-York-Mollie (BYM) specification
#example: Sandeels, density (include excess zeros) and Porpoises,2005 
################################################################################################################

#Loading the data file
species1 <- readRDS("species1.rds")#Porpoise (2005)
species2 <- readRDS("Sandeels.rds")#Sandeels, density
density <- species2$density

n1<-dim(species1)[1]
n2<-dim(species2)[1]

nothing1 = rep(NA, n1)
nothing2 = rep(NA, n2)

#define the two sandeels response variables
z <- as.vector((density>0) + 0)
y <- as.vector(ifelse(density>0, density, NA))

#Generating 3 vectors of response variables for the joint model 
w = as.vector(species1$density)
wNA = as.vector(c(w, nothing2, nothing2))

zNA = as.vector(c(nothing1, z, nothing2))
yNA = as.vector(c(nothing1, nothing2, y))

#Combine them in a matrix
outcome.matrix<-matrix(c(wNA,zNA, yNA), ncol=3)

#A factor vector with 3 levels indicating 3 outcome variables
mu = as.factor(c(rep(1,length(w)), rep(2,length(z)), rep(3,length(y))))

#Index vectors for the spatial effects 
i.spat1 = c(species1$IDgrid, nothing2, nothing2)#Porpoise
i.spat2 = c(nothing1, species2$IDgrid, nothing2)#Sandeels
i.spat3 = c(nothing1, nothing2, species2$IDgrid)#Sandeels

#The covariates 
BT1 <-c(species1$BT, nothing2, nothing2)#Porpoise
NPP2 <-c(nothing1, species2$NPP, nothing2)#Sandeels
NPP3<-c(nothing1, nothing2, species2$NPP)#Sandeels
DVV3 <-c(nothing1, nothing2, species2$DVV)#Sandeels

#Data set
data=list(outcome.matrix=outcome.matrix, i.spat1=i.spat1, i.spat2=i.spat2, i.spat3=i.spat3,
          BT1=BT1, NPP2=NPP2, NPP3=NPP3, DVV3=DVV3
)

formula = outcome.matrix ~ mu -1 +
  f(inla.group(BT1), model="rw2", hyper=param.cc02) +
  f(inla.group(NPP2), model="rw2", hyper=param.cc02) +
  f(inla.group(NPP3), model="rw2", hyper=param.cc02) +
  f(inla.group(DVV3), model="rw2", hyper=param.cc02) +  
  f(i.spat1, model="bym", graph = "graph3.dat", hyper= param.spat02) +
  f(i.spat2, copy="i.spat1", fixed=F) +
  f(i.spat3, copy="i.spat1", fixed=F)


result.j2 <- inla(formula,family=c('gamma','binomial', 'gamma'), data=data, 
               control.compute=list(dic=T), 
               control.predictor = list(compute = TRUE), verbose=TRUE,
               control.fixed = list(expand.factor.strategy = "inla")) 

################################################################################################################
#Joint modelling with misalignment 
#using the stochastic partial differential equation (SPDE) approach
#example: Black-legged kittiwake, observations
#
#The full R code is available here:
#Krainski, E.T., Lindgren, F., Simpson, D. & Rue, H. (2015) The R-INLA tutorial on SPDE models: Chapter 7
################################################################################################################

#Loading the data file
data0 <- readRDS("kittiwake_obs.rds")
data.c <- readRDS("biophysical_var.rds")

#Coordinates: bio/physical variables 
loc.c <- data.c[,1:2]
#Coordinates: species locations 
loc.y <- data0[,1:2]

#Projector matrices
Ac <- inla.spde.make.A(mesh, loc=loc.c)
Ay <- inla.spde.make.A(mesh, loc=loc.y)

#SPDE model
spde <- inla.spde2.matern(mesh=mesh)

#The data stack
stk.c <- inla.stack(data=list(y=cbind(data.c$NPP, NA)),
                    A=list(Ac), tag='dat.c',
                    effects=list(m=1:mesh$n))
stk.y <- inla.stack(data=list(y=cbind(NA, data0$Numbers)),
                    A=list(Ay, 1),
                    effects=list(list(c.y=1:mesh$n, x=1:mesh$n),
                                 list(a.y=rep(1,length(data0$Numbers)))))
stk <- inla.stack(stk.c, stk.y)

#The estimation of the regression coefficient in this approach is treated as 
#a hyperparameter, such as copy parameter of an latent field (Krainski et al. 2015).
formula <- y ~ 0 + a.y + 
  f(m, model=spde) + f(x, model=spde) + 
  f(c.y, copy='m', fixed=FALSE, hyper=hyper.y)

#Fitting the model 
result <- inla(form, data=inla.stack.data(stk), family=c('gaussian','poisson'),
            control.predictor=list(compute=TRUE, A=inla.stack.A(stk)))

#Posterior mean and posterior standard deviations:
mesh2locs <- rBind(Ac, Ay)
m.mprd <- drop(mesh2locs%*%res$summary.ran$m$mean)#mean
sd.mprd <- drop(mesh2locs%*%res$summary.ran$m$sd)#sd