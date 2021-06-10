library(ipmr)
library(lme4)
library(dplyr)

data <- read.csv("ipmr_examples/data/lupine_all.csv") %>%
  select(location, year,
         log_area_t0, 
         flow_t0, 
         numrac_t0, 
         surv_t1, 
         log_area_t1,
         numab_t0)

seed <- read.csv("ipmr_examples/data/seedsperfruit.csv")
fruit <- read.csv("ipmr_examples/data/fruits_per_raceme.csv")

mean_fruit <- mean(fruit$NumFruits)
mean_seed  <- mean(seed$SEEDSPERFRUIT)

data <- mutate(data, seeds = #(
                 numrac_t0 
               #- numab_t0)
               * mean_fruit * mean_seed)

