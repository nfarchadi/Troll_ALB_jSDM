
#####
## Simple VAST model with covariates?
#####
library(VAST)
library(splines)
library(tidyverse)
example<- load_example( data_set="covariate_example" )
example$sampling_data %>% view()
example$covariate_data %>% view()

# Make settings (turning off bias.correct to save time for example)
settings <- make_settings(
    n_x = 100,
    Region = example$Region,
    purpose = "index2",
    FieldConfig = c("Omega1" = 0, "Omega2" = 0, "Epsilon1" = 0, "Epsilon2" = 0),
    RhoConfig = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0),
    ObsModel = c(1, 0),
    use_anisotropy = FALSE,
    bias.correct = FALSE,
    fine_scale = TRUE
)


# Define formula
X1_formula = ~ log(BOT_DEPTH)
X2_formula = X1_formula

# If all covariates as "static" (not changing among years), then set Year = NA to cause values to be duplicated internally for all values of Year
# If using a mix of static and dynamic covariates, then duplicate rows for static covariates for every value of Year
# Here, all covariates are static, so I'm using the first approach.
example$covariate_data[,'Year']<- NA

# Rescale covariates being used to have an SD >0.1 and <10 (for numerical stability)
example$covariate_data[,'BOT_DEPTH'] <- example$covariate_data[,'BOT_DEPTH'] / 100

# Run model
?fit_model()
fit_vast<- fit_model( "settings" = settings,
  Lat_i = example$sampling_data[,'Lat'],
  Lon_i = example$sampling_data[,'Lon'],
  t_i = example$sampling_data[,'Year'],
  b_i = example$sampling_data[,'Catch_KG'],
  a_i = rep(1, nrow(example$sampling_data)),
  X1_formula = X1_formula,
  X2_formula = X2_formula,
  covariate_data = example$covariate_data )


#####
## Can we fit this with sdmTMB
#####
# install.packages("sdmTMB", dependencies = TRUE)
library(sdmTMB)

data_sdmtmb <- data.frame("Year" = example$sampling_data[, "Year"], "Lat" = example$sampling_data[, "Lat"], "Lon" = example$sampling_data[, "Lon"], "Catch_KG" = example$sampling_data[, "Catch_KG"], "BOT_DEPTH" = example$covariate_data[, "BOT_DEPTH"])
data_sdmtmb$Presence = ifelse(example$sampling_data[,"Catch_KG"] > 0, 1, 0)


mesh_use<- fit_vast$spatial_list$MeshList$anisotropic_mesh 
mesh_sdmtmb <- sdmTMB::make_mesh(data_sdmtmb, c("Lon", "Lat"), mesh = mesh_use)

fit_tmb <- sdmTMB(
    Catch_KG ~ 0 + log(BOT_DEPTH),
    data = data_sdmtmb,
    mesh = mesh_sdmtmb,
    family = delta_lognormal(),
    time = "Year",
    time_varying = ~ 1,
    spatial = "off",
    spatiotemporal = "off"
)

#####
## Did that...work?
#####
summary(fit_tmb)
fit_vast$parameter_estimates$SD

# Coefficients for log(BOT_DEPTH) are spot on, the yearly intercepts are a bit off...not sure why..but also doesn't seem problematic.

#####
## sdmTMB presence absence?
#####
fit_tmb <- sdmTMB(
    Presence ~ 0 + log(BOT_DEPTH),
    data = data_sdmtmb,
    mesh = mesh_sdmtmb,
    family = binomial(),
    time = "Year",
    time_varying = ~ 1,
    spatial = "off",
    spatiotemporal = "off"
)

#####
## VAST?
#####
settings_pa<- make_settings(
    n_x = 100,
    Region = example$Region,
    purpose = "index2",
    FieldConfig = c("Omega1" = 0, "Omega2" = 0, "Epsilon1" = 0, "Epsilon2" = 0),
    RhoConfig = c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0),
    ObsModel = c(1, 0),
    use_anisotropy = FALSE,
    bias.correct = FALSE,
    fine_scale = TRUE
)

fit_vast<- fit_model( "settings" = settings,
  Lat_i = example$sampling_data[,'Lat'],
  Lon_i = example$sampling_data[,'Lon'],
  t_i = example$sampling_data[,'Year'],
  b_i = as_units(ifelse(example$sampling_data[, "Catch_KG"] > 0, 1, 0), unitless),
  a_i = rep(1, nrow(example$sampling_data)),
  X1_formula = X1_formula,
  covariate_data = example$covariate_data )

#####
## Did that...work?
#####
summary(fit_tmb)
fit_vast$parameter_estimates$SD

# Ignore the warning because its all the second linear predictor, which we don't care about. Otherwise, the coefficient estimates look pretty close to sdmTMB?
