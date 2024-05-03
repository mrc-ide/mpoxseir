## clean workflow to run the model and plot results 

library(odin.dust)
library(tidyverse)

# compiles odin.dust model
seird <- odin.dust::odin_dust("inst/odin/simple_SEIR_odindust.R") 

# model inputs
# parameterising to DRC

# number of age classes (need this to be able to input anything)


# population size 
# https://www.populationpyramid.net/democratic-republic-of-the-congo/2019/
# this is entered into the model as S0

S0 <- c()
