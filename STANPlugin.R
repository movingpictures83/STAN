# replace "~/Desktop/Shared Code for Soon" by "folder-location"

# load libraries 
library(tidyverse)
library(rstan)
library(bayesplot)
library(loo)
library(bayestestR)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
n_chains <- 4

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix <<- prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }
}

run <- function() {}

output <- function(outputfile) {

df_stan <- readRDS(paste(pfix, parameters["RDS", 2], sep="/"))
m00 <- stan(file=str_c(paste(pfix, parameters["model", 2], sep="/")),
            data = df_stan , 
            chains=5, iter=11000, warmup=1000, thin=10,cores = 1)
saveRDS(m00, outputfile)
}
