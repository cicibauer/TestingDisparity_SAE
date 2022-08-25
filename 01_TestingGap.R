################################################################
# This R script is to process data and calculate testing gaps
################################################################

#--------------------load required packages ---------------------#
library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)

#------------------------- Read in data ------------------------#
# Weekly counts of tests and cases and population size
data <- read_excel("Data/weekly_data.xlsx") # 52 weeks and 222 GEOIDs

#------------------------- calculate testing rank ------------------------#

# Calculate test intensity and positivity rate
#### rank test intensity per 10,000 population and positivity rate by week ########
#### rank = 1: highest test intensity ########
#### rank = 0: 0 test intensity ########
data_block_gap <- data %>%
  mutate(test_intensity = test_freq*10000/pop,
         positive_rate = positive_freq/test_freq) %>%
  group_by(week) %>%
  mutate(rank_test_intensity = as.numeric(factor(rank(-test_intensity))),
         rank_pos_rate = as.numeric(factor(rank(-positive_rate)))) %>%
  arrange(week)

# replace missing with 0
data_block_gap$rank_test_intensity[data_block_gap$test_intensity == 0] <- 0
data_block_gap$rank_pos_rate[is.na(data_block_gap$positive_rate)] <- 0

# calculate total testing intensity
data_block_gap$total_test <- ave(data_block_gap$test_freq, data_block_gap$GEOID, FUN=sum)
data_block_gap$total_intensity <- data_block_gap$total_test*10000/data_block_gap$pop

#----------------------------calculate testing gap------------------------#
### Calculate gap: rank=1 means highest positive rate and highest test intensity
### Gap exits if the rank of positivity rate is lower than rank of test intensity 
data_block_gap$rank_unmatch <- ifelse(data_block_gap$rank_pos_rate < data_block_gap$rank_test_intensity, 1 ,0)
data_block_gap$rank_unmatch[data_block_gap$positive_rate == 0] <- 0 #rank_unmatch == 1 gap needs to be calculated

# Expected test intensity by week
setDT(data_block_gap)[, expect_test_intensity := test_intensity[match(rank_pos_rate, rank_test_intensity)], week]
# calculate expected number of tests
data_block_gap$expect_test_freq <- data_block_gap$expect_test_intensity * data_block_gap$pop / 10000
# Gap is the difference between expected tests and real tests
data_block_gap$gap <- ifelse(data_block_gap$rank_unmatch == 1,
                             data_block_gap$expect_test_freq - data_block_gap$test_freq,
                             0)

#calculate cumulative gap by GEOID
data_block_gap$cum_gap <- ave(data_block_gap$gap, data_block_gap$GEOID, FUN=cumsum)
data_block_gap <- 
  data_block_gap %>% 
  group_by(GEOID) %>% 
  arrange(week) %>%
  mutate(cum_gap = cumsum(gap),mean_gap = max(cum_gap)/max(week)) %>%
  ungroup()

# save the processed data
save(data_block_gap, file="Data/processed_data_sp_model.rda")



