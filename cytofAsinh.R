cytofAsinh <- function(value, cofactor = 5) { 
    value <- value-1 
    loID <- which(value < 0) 
    if(length(loID) > 0) 
        value[loID] <- rnorm(length(loID), mean = 0, sd = 0.01) 
    value <- value / cofactor 
    value <- asinh(value) # value <- log(value + sqrt(value^2 + 1)) 
    return(value) 
} 
