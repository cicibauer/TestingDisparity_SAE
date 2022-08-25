#########################################################
# This R script applies the Spatial inseparable models 
# to obtained smoothed estimates of test positivity rate 
#########################################################

#-----------------------------load required packages ----------------------------------#
library(stringi); library(dplyr);library(tidyverse); library(tigris); library(tidycensus); library(tmap); 
library(readxl); library(INLA); library(spdep);library(ggcorrplot);library(rgdal);library(ggmap);library(ggplot2)
library("DClusterm"); library(yarrr)

#--------------------------------- Read in data --------------------------------------#
# read in processed data (based on the simulated dataset)
load("Data/processed_data_sp_model.rda")
### set the new test frequency to 1, if the positive_freq is NA
data_block_gap$positive_freq[data_block_gap$test_freq == 0] <- NA
data_block_gap$new_test_freq <- data_block_gap$test_freq 
data_block_gap$new_test_freq[data_block_gap$test_freq == 0] <- 1 

### sanity check
data_block_gap %>% dplyr::select(GEOID, week, test_freq, positive_freq, pop) %>% arrange(desc(test_freq))

### read in the adjacency matrix to fit the inseparable models
g <- inla.read.graph(filename="Data/adjacency_matrix.adj")
#data_block_gap$ID <- 1:nrow(data_block_gap)

#------------------------------------ Modeling ---------------------------------------#
######################################
### Models using rw2 latent effect ##
######################################
dat <- data_block_gap %>% mutate(ID.area = GEOID, ID.time=week, ID.time2=week, ID.area.time = 1:nrow(data_block_gap))
#--- Type I interaction ---#
# Formula
formula.intI <- positive_freq ~ f(ID.area, model="bym", graph= g) + f(ID.time, model="rw2") + 
  f(ID.time2,model="iid") + f(ID.area.time, model="iid")
# Run the model
mod.intI <- inla(formula.intI, data=dat, 
                 family="binomial", Ntrials= new_test_freq, control.family = list(link = "logit"),  
                 control.predictor=list(compute=TRUE, link=1), 
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))

#--- Type II interaction ---#
dat$ID.area.int <- dat$ID.area
dat$ID.time.int <- dat$ID.time
# Formula
formula.intII <- positive_freq ~ f(ID.area,model="bym",graph= g) +
  f(ID.time,model="rw2") + f(ID.time2,model="iid") +
  f(ID.area.int, model="iid", group=ID.time.int, control.group=list(model="rw2")) 
# Run the model
mod.intII <- inla(formula.intII, data=dat, 
                  family="binomial", Ntrials= new_test_freq, control.family = list(link = "logit"),  
                  control.predictor=list(compute=TRUE, link=1), 
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))

#--- Type III interaction ---#
# Formula
formula.intIII<- positive_freq ~  f(ID.area,model="bym",graph= g) +
  f(ID.time,model="rw2") + f(ID.time2,model="iid") +
  f(ID.time.int,model="iid", group=ID.area.int, control.group=list(model="besag", graph= g))
# Run the model
mod.intIII <- inla(formula.intIII, data=dat, 
                   family="binomial", Ntrials= new_test_freq, control.family = list(link = "logit"),  
                   control.predictor=list(compute=TRUE, link=1), 
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))

#--- Type IV interaction ---#
#Formula
formula.intIV <- positive_freq ~ f(ID.area,model="bym",graph= g) +
  f(ID.time,model="rw2") + f(ID.time2,model="iid") + 
  f(ID.area.int, model="besag", graph= g, group=ID.time.int,control.group=list(model="rw2"))
# Run the model
mod.intIV <- inla(formula.intIV, data=dat, 
                  family="binomial", Ntrials= new_test_freq, control.family = list(link = "logit"),  
                  control.predictor=list(compute=TRUE, link=1), 
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))

### Look at the summary of the models 
summary(mod.intI)
summary(mod.intII)
summary(mod.intIII)
summary(mod.intIV)

# calculate the observed test positivity rates
dat$rate_pos <- dat$positive_freq/dat$test_freq
# get the estimated test positivity rate using posterior mean
fitI <- mod.intI$summary.fitted.values[, 1]
fitII <- mod.intII$summary.fitted.values[, 1]
fitIII <- mod.intIII$summary.fitted.values[, 1]
fitIV <- mod.intIV$summary.fitted.values[, 1]

### Scatter plot compares the observed positivity rates and the estimated rates.
par(mfrow=c(2,2))
plot(dat$rate_pos, fitI, col=2, xlab="Observed Positivity Rate", ylab="Estimated Positivity Rate (Model 1 - rw2)", xlim=c(0, 1), ylim=c(0,1), cex.lab = 0.7)
abline(a=0, b=1, col="grey", lwd=2)

plot(dat$rate_pos, fitII, col=2, xlab="Observed Positivity Rate", ylab="Estimated Positivity Rate (Model 2 - rw2)", xlim=c(0, 1), ylim=c(0,1), cex.lab = 0.7)
abline(a=0, b=1, col="grey", lwd=2)

plot(dat$rate_pos, fitIII, col=2, xlab="Observed Positivity Rate", ylab="Estimated Positivity Rate (Model 3 - rw2)", xlim=c(0, 1), ylim=c(0,1), cex.lab = 0.7)
abline(a=0, b=1, col="grey", lwd=2)

plot(dat$rate_pos, fitIV, col=2, xlab="Observed Positivity Rate", ylab="Estimated Positivity Rate (Model 4 - rw2)", xlim=c(0, 1), ylim=c(0,1), cex.lab = 0.7)
abline(a=0, b=1, col="grey", lwd=2)

# save the model results (if needed)
save(list = c("mod.intIV","mod.intI","mod.intII","mod.intIII"), file="mods_rw2.Rdata")

# save the data file for testing disparity prediction
save(dat, file="Data/data_for_pred.rda")
