######################################################
# This R script is for predicting positivity rate
######################################################


### data preparation
load("Data/data_for_pred.rda")
#dat$GEOID <- as.character(dat$GEOID)
dat$gap <- ceiling(dat$gap)

# adjacency matrix
g <- inla.read.graph("Data/adjacency_matrix.adj")


#### define function of calculating gap - using previous week testing intensity
calculate_gap <- function(dataset){
  dataset <- dataset %>% 
    group_by(week) %>%
    mutate(rank_test_intensity = as.numeric(factor(rank(-test_intensity, ties.method = "random"))),
           rank_est_pos_rate = as.numeric(factor(rank(-estimated_pos_rate))))
  
  ###calculate gap, !!!!!rank=1 means highest positive rate, highest test intensity!!!!!!
  dataset[,"rank_unmatch"] <- ifelse(dataset[, "rank_est_pos_rate"] < dataset[, "rank_test_intensity"], 1 ,0)
  dataset[,"rank_unmatch"][dataset[,"positive_rate"] == 0] <- 0 #rank_unmatch == 1 gap needs to be calculated
  
  #this will get the expected test intensity
  setDT(dataset)[, expect_test_intensity := test_intensity[match(rank_est_pos_rate, rank_test_intensity)], week]
  
  #calculate expected test frequency
  dataset[,"expect_test_freq"] <- dataset[,"expect_test_intensity"] * dataset[,"census_pop"] / 10000
  
  #calculate testing gap and round to integer
  dataset[, "gap"] <- dataset[, "expect_test_freq"] - dataset[, "test_freq"]
  dataset[, "gap"][dataset[,"rank_unmatch"] == 0] <- 0
  dataset[, "gap"] <- ceiling(dataset[, "gap"])
  
  return(dataset)
}


#### define function of calculating gap - using overall average testing intensity

calculate_gap2 <- function(dataset, pred_week){
  dataset <- as.data.frame(dataset)
  dat_total_test <- dataset[dataset[,"week"] != pred_week, ]
  total_test <- dat_total_test %>% group_by(block_grp) %>% summarise(across(c(test_freq), list(sum = sum))) #get total tests before pred week
  dataset <- merge(dataset, total_test, by = "block_grp") #merge dataset with total test num
  
  ###calculate overall average testing intensity
  dataset$ave_test <- dataset$test_freq_sum / (pred_week - 1) #average test for each cbg
  dataset$ave_test_intensity <- dataset$ave_test * 10000 / dataset$census_pop 
  
  dataset <- dataset %>% 
    group_by(week) %>%
    mutate(rank_ave_test_intensity = as.numeric(factor(rank(-ave_test_intensity, ties.method = "random"))),
           rank_est_pos_rate = as.numeric(factor(rank(-estimated_pos_rate))))
  
  ###calculate gap, !!!!!rank=1 means highest positive rate, highest test intensity!!!!!!
  dataset[,"rank_unmatch2"] <- ifelse(dataset[, "rank_est_pos_rate"] < dataset[, "rank_ave_test_intensity"], 1 ,0)
  dataset[,"rank_unmatch2"][dataset[,"positive_rate"] == 0] <- 0 #rank_unmatch == 1 gap needs to be calculated
  
  #this will get the expected test intensity
  setDT(dataset)[, expect_ave_test_intensity := ave_test_intensity[match(rank_est_pos_rate, rank_ave_test_intensity)], week]
  
  #calculate expected test frequency
  dataset[,"expect_ave_test_freq"] <- dataset[,"expect_ave_test_intensity"] * dataset[,"census_pop"] / 10000
  
  #calculate testing gap and round to integer
  dataset[, "gap2"] <- dataset[, "expect_ave_test_freq"] - dataset[, "ave_test"]
  dataset[, "gap2"][dataset[,"rank_unmatch2"] == 0] <- 0
  dataset[, "gap2"] <- ceiling(dataset[, "gap2"])
  
  return(dataset)
}


#### define function of calculating gap - using prev 4 weeks average testing intensity 
calculate_gap3 <- function(dataset, pred_week){
  dataset <- as.data.frame(dataset)
  dat_total_test <- dataset[dataset[,"week"] %in% (pred_week-4):(pred_week-1), ]
  total_test <- dat_total_test %>% group_by(block_grp) %>% summarise(across(c(test_freq), list(sum2 = sum))) #get total tests before pred week
  dataset <- merge(dataset, total_test, by = "block_grp") #merge dataset with total test num
  
  ###calculate overall average testing intensity
  dataset$ave_test.4w <- dataset$test_freq_sum2 / 4 #average test for each cbg
  dataset$ave_test_intensity.4w  <- dataset$ave_test.4w * 10000 / dataset$census_pop 
  
  dataset <- dataset %>% 
    group_by(week) %>%
    mutate(rank_ave_test_intensity.4w = as.numeric(factor(rank(-ave_test_intensity.4w, ties.method = "random"))),
           rank_est_pos_rate = as.numeric(factor(rank(-estimated_pos_rate))))
  
  ###calculate gap, !!!!!rank=1 means highest positive rate, highest test intensity!!!!!!
  dataset[,"rank_unmatch3"] <- ifelse(dataset[, "rank_est_pos_rate"] < dataset[, "rank_ave_test_intensity.4w"], 1 ,0)
  dataset[,"rank_unmatch3"][dataset[,"positive_rate"] == 0] <- 0 #rank_unmatch == 1 gap needs to be calculated
  
  #this will get the expected test intensity
  setDT(dataset)[, expect_ave_test_intensity.4w := ave_test_intensity[match(rank_est_pos_rate, rank_ave_test_intensity.4w)], week]
  
  #calculate expected test frequency
  dataset[,"expect_ave_test_freq.4w"] <- dataset[,"expect_ave_test_intensity.4w"] * dataset[,"census_pop"] / 10000
  
  #calculate testing gap and round to integer
  dataset[, "gap3"] <- dataset[, "expect_ave_test_freq.4w"] - dataset[, "ave_test.4w"]
  dataset[, "gap3"][dataset[,"rank_unmatch3"] == 0] <- 0
  dataset[, "gap3"] <- ceiling(dataset[, "gap3"])
  
  return(dataset)
}


#### define function to generate data set for prediction
dat_pred <- function(data, week_to_pred){
  temp_dat <- data
  ##generate data for week41
  #observed week 40
  week_prev <- temp_dat[temp_dat$week == week_to_pred-1,]
  pred <- week_prev %>%
    mutate(week = week_to_pred, new_test_freq = 1, 
           positive_freq = NA, positive_rate = NA, 
           ID.time = week_to_pred, ID.time2 = week_to_pred, 
           ID.time.int = week_to_pred)
  
  
  temp_dat <- rbind(subset(temp_dat, temp_dat$week %in% 1:(week_to_pred-1)), pred)
  
  #update variable ID.area.time 
  temp_dat$ID.area.time <- 1:nrow(temp_dat)
  
  return(temp_dat)
}

################################################
### Prediction using rw2 latend effect - week 14
################################################
dat1 <- dat_pred(data = dat, week_to_pred = 14)
#--- Type I interaction ---#
formula.intI <- positive_freq ~ f(ID.area, model="bym", graph=g) + f(ID.time, model="rw2") + 
  f(ID.time2,model="iid") + f(ID.area.time, model="iid")
pred.intI.rw2 <- inla(formula.intI, data=dat1, 
                      family="binomial", Ntrials= new_test_freq, control.family = list(link = "logit"),  
                      control.predictor=list(compute=TRUE, link=1), 
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))

#--- Type II interaction ---#
formula.intII <- positive_freq ~ f(ID.area,model="bym",graph=g) +
  f(ID.time,model="rw2") + f(ID.time2,model="iid") +
  f(ID.area.int, model="iid", group=ID.time.int, control.group=list(model="rw2")) 
pred.intII.rw2 <- inla(formula.intII, data=dat1, 
                       family="binomial", Ntrials= new_test_freq, control.family = list(link = "logit"),  
                       control.predictor=list(compute=TRUE, link=1), 
                       control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))

#--- Type III interaction ---#
formula.intIII<- positive_freq ~  f(ID.area,model="bym",graph=dat.mat) +
  f(ID.time,model="rw2") + f(ID.time2,model="iid") +
  f(ID.time.int,model="iid", group=ID.area.int, control.group=list(model="besag", graph=dat.mat))
pred.intIII.rw2 <- inla(formula.intIII, data=dat1, 
                        family="binomial", Ntrials= new_test_freq, control.family = list(link = "logit"),  
                        control.predictor=list(compute=TRUE, link=1), 
                        control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))

#--- Type IV interaction ---#
formula.intIV <- positive_freq ~ f(ID.area,model="bym",graph=dat.mat) +
  f(ID.time,model="rw2") + f(ID.time2,model="iid") + 
  f(ID.area.int, model="besag", graph=dat.mat, group=ID.time.int,control.group=list(model="rw2"))
pred.intIV.rw2 <- inla(formula.intIV, data=dat1, 
                       family="binomial", Ntrials= new_test_freq, control.family = list(link = "logit"),  
                       control.predictor=list(compute=TRUE, link=1), 
                       control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))


## Look at the summary of the models 
summary(pred.intI.rw2)
summary(pred.intII.rw2)
summary(pred.intIII.rw2)
summary(pred.intIV.rw2)

# Save the model results if needed
#save(list = c("pred.intI.rw2","pred.intII.rw2","pred.intIII.rw2","pred.intIV.rw2"), file="pred_w14_020722dat_rw2.Rdata")
#load("pred_w14_020722dat_rw2.Rdata")

# Get the estimated positivity rate and confidence interval
dat1$GEOID <- as.character(dat1$GEOID)
est_dat1 <- data.frame(dat1)
est_dat1$estimated_pos_rate <- pred.intII.rw2$summary.fitted.values[, 1]
est_dat1$lower_ci <- pred.intII.rw2$summary.fitted.values[, 3]
est_dat1$upper_ci <- pred.intII.rw2$summary.fitted.values[, 5]

##calculate testing gap
data_w14_int2  <- calculate_gap(est_dat1)
data_w14_int2  <- calculate_gap2(data_w14_int2)

# Plot to compare the estimated positivity rate and observed positivty rate for week 14.
plot(data_w14_int2$estimated_pos_rate[data_w14_int2$week == 14], 
     dat$positive_rate[dat$week == 14], col=2, 
     xlab="Observed Positivity Rate", 
     ylab="Estimated Positivity Rate (Model 2 - rw2)", 
     xlim=c(0,1),
     ylim= c(0,1),
     cex.lab = 0.7)
abline(a=0, b=1, col="grey", lwd=2)
