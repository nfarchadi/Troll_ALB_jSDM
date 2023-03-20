calculate_n2 <- function(x) {
  if (all(is.na(x))){
    rep(NA, length.out=length(dep[1:23]))
  } else if(any(is.na(x[1:10]))){ ## if HYCOM is NA in top 25 m, not worth trying (and failing) to calc N2. (I changed 1:10 to 1:4)-NF
    rep(NA, length.out=length(dep[1:23]))
  } 
  else{
    oce::swN2(pressure = dep[1:23], sigmaTheta = oce::swSigmaTheta(x[24:46], x[1:23], pressure = dep[1:23]))#, referencePressure = median(depth[1:23], na.rm = TRUE))) ##24:26 is salinity layers, 1:23 is temp (or depth)
  }
}