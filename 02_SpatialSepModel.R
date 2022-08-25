######################################################
#This R script is for smoothing positivity rate since
# its raw data contains lots of 0 and 1. 
# Spatial inseparable models
######################################################

#-----------------------------load required packages ----------------------------------#
library(stringi); library(dplyr);library(tidyverse); library(tigris); library(tidycensus); library(tmap); 
library(readxl); library(INLA); library(spdep);library(ggcorrplot);library(rgdal);library(ggmap);library(ggplot2)
library("DClusterm"); library(yarrr)

#--------------------------------- Read in data --------------------------------------#
# read in processed data
load("Data/processed_data_sp_model.rda")
#data_block_gap$GEOID <- as.character(data_block_gap$GEOID)
### set the new test frequency to 1, if the positive_freq is NA
data_block_gap$positive_freq[data_block_gap$test_freq == 0] <- NA
data_block_gap$new_test_freq <- data_block_gap$test_freq 
data_block_gap$new_test_freq[data_block_gap$test_freq == 0] <- 1 #to make prediction change 0 to 1

# sanity check
data_block_gap %>% dplyr::select(GEOID, week, test_freq, positive_freq, pop) %>% arrange(desc(test_freq))

# read in adjacency matrix
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

# calculate the observed posiitivity rates
dat$rate_pos <- dat$positive_freq/dat$test_freq
# get the estimated positivity rate
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

### Check CPO and dic, finish the code, so that we can decide on the best fitting model.
cpo.rw2 <- c(-sum(log(mod.intI$cpo$cpo), na.rm=T), -sum(log(mod.intII$cpo$cpo), na.rm=T), -sum(log(mod.intIII$cpo$cpo), na.rm=T), -sum(log(mod.intIV$cpo$cpo), na.rm=T))
pit.rw2 <- c(-sum(log(mod.intI$cpo$pit), na.rm=T), -sum(log(mod.intII$cpo$pit), na.rm=T), -sum(log(mod.intIII$cpo$pit), na.rm=T), -sum(log(mod.intIV$cpo$pit), na.rm=T))

dic.rw2 <- c(mod.intI$dic$dic, mod.intII$dic$dic, mod.intIII$dic$dic, mod.intIV$dic$dic)
waic.rw2 <-c(mod.intI$waic$waic, mod.intII$waic$waic, mod.intIII$waic$waic, mod.intIV$waic$waic)

#sigma^2 = 1/median
sigma2.rw2 <- cbind(1/mod.intI$summary.hyperpar$`0.5quant`,
                    1/mod.intII$summary.hyperpar$`0.5quant`,
                    1/mod.intIII$summary.hyperpar$`0.5quant`,
                    1/mod.intIV$summary.hyperpar$`0.5quant`)

# Combine together into a table
tab_rw2 <- rbind(cpo=cpo.rw2, waic=waic.rw2, dic=dic.rw2, sigma2 = sigma2.rw2)
tab_rw2 <- apply(tab_rw2, 2, function(x){sprintf("%.4f",x)})
rownames(tab_rw2) <- c("CPO","WAIC","DIC","Sigma^2: spatial_unstructured",
                       "Sigma^2: spatial_structured", "Sigma^2: time_structured",
                       "Sigma^2: time_unstructured", "Sigma^2: interaction")

tab_rw2

# inla.cpo(mod.intI,
#          force = FALSE,
#          verbose = TRUE,
#          recompute.mode = TRUE)

# save the model results (if needed)
save(list = c("mod.intIV","mod.intI","mod.intII","mod.intIII"), file="mods_rw2.Rdata")

# save the data file for prediction
save(dat, file="Data/data_for_pred.rda")

############################
#### plot of random effect
############################
###spatial
par(mfrow=c(2,2))
spatial1 <- mod.intI$summary.random$ID.area
spatial2 <- mod.intII$summary.random$ID.area
spatial3 <- mod.intIII$summary.random$ID.area
spatial4 <- mod.intIV$summary.random$ID.area

plot(spatial1$mean ~ spatial1$ID)
plot(spatial2$mean ~ spatial2$ID)
plot(spatial3$mean ~ spatial3$ID)
plot(spatial4$mean ~ spatial4$ID)

### structured time
time_str1 <- mod.intI$summary.random$ID.time
time_str2 <- mod.intII$summary.random$ID.time
time_str3 <- mod.intIII$summary.random$ID.time
time_str4 <- mod.intIV$summary.random$ID.time

plot(time_str1$mean ~ time_str1$ID)
plot(time_str2$mean ~ time_str2$ID)
plot(time_str3$mean ~ time_str3$ID)
plot(time_str4$mean ~ time_str4$ID)

### iid time
time_unstr1 <- mod.intI$summary.random$ID.time2
time_unstr2 <- mod.intII$summary.random$ID.time2
time_unstr3 <- mod.intIII$summary.random$ID.time2
time_unstr4 <- mod.intIV$summary.random$ID.time2

plot(time_unstr1$mean ~ time_unstr1$ID)
plot(time_unstr2$mean ~ time_unstr2$ID)
plot(time_unstr3$mean ~ time_unstr3$ID)
plot(time_unstr4$mean ~ time_unstr4$ID)

datout <- data.frame(dat3[1:11])
datout$estimated_pos_rate <- mod.intII$summary.fitted.values[, 1]
datout$lower_ci <- mod.intII$summary.fitted.values[, 3]
datout$upper_ci <- mod.intII$summary.fitted.values[, 5]

