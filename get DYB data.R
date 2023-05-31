library(tidyverse)
library(ddharmony)
library(censusAdjust)
library(DemoTools)
library(DDSQLtools)

load("C:/Users/ktang3/Desktop/Imperial/Life_Table_System/child mortality DHS BR.Rdata")
pyears_data <- readRDS("C:/Users/ktang3/Desktop/Imperial/SSA_mort/JE/data/pyears_data_smoothed.rds")
aggr.mat.cohort.0 <- pyears_data %>%
  rename(DHS = SurveyId) %>%
  group_by(country) %>% group_split %>%
  setNames(pyears_data$country %>% levels) 

joint.countries1 <- names(aggr.mat.br)

joint.countries2 <- names(aggr.mat.cohort.0)

joint.countries <- intersect(joint.countries1, joint.countries2) %>% setdiff(c("Rwanda", "Mali"))

locationid <- get_locations()

open.age <- 85

pop.m.list <- pop.f.list <- pop.m.5.list <- pop.f.5.list <- list()
deaths.m.list <- deaths.f.list <- deaths.m.5.list <- deaths.f.5.list <- list()

for(country in joint.countries){
  
  try(
    census_pop_counts <- DDharmonize_validate_PopCounts(locid = ifelse(country=="Cote d'Ivoire", 384,
                                                                       ifelse(country=="Tanzania",834,
                                                                              locationid$PK_LocID[which(locationid$Name == country)])),
                                                        times = 1950:2020,
                                                        process="census",
                                                        DataSourceShortName = "DYB")
  )
  
  try(
    census_deaths <- DDharmonize_validate_DeathCounts(locid = ifelse(country=="Cote d'Ivoire", 384,
                                                                     ifelse(country=="Tanzania",834,
                                                                            locationid$PK_LocID[which(locationid$Name == country)])),
                                                      times = 1950:2020,
                                                      process="census",
                                                      DataSourceShortName = "DYB",
                                                      retainKeys = TRUE)
  )
  
  
  if(exists("census_pop_counts") && !is.null(census_pop_counts)){
  ddharm_bf_census_m <- census_pop_counts %>%  
    filter(AgeSpan %in% c(-1, 1), AgeLabel != "Total", TimeLabel > 1960, SexID == 1, complete == TRUE) %>%
    select(TimeLabel, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID, DataSourceName) %>%
    distinct() %>%
    arrange(TimeLabel) %>%
    pivot_wider(names_from = TimeLabel, values_from = DataValue) %>%
    group_by(AgeStart, AgeLabel, AgeSpan, DataSourceName) %>%
    arrange(AgeStart, StatisticalConceptName) %>% ###De-facto first
    select(-c(StatisticalConceptName, SexID)) %>%
    summarise_all(function(y){first(na.omit(y))}) %>%
    ungroup() %>%
    select(-AgeSpan, -AgeStart) %>%
    mutate(DataSourceName = replace(DataSourceName, DataSourceName == "IPUMS International custom tabulations", "IPUMS"),
           DataSourceName = replace(DataSourceName, DataSourceName == "Demographic Yearbook", "DYB"),
           country = country,
           age = AgeLabel,
           .before = 3) %>%
    select(-AgeLabel)
  
  ddharm_bf_census_f <- census_pop_counts %>%  
    filter(AgeSpan %in% c(-1, 1), AgeLabel != "Total", TimeLabel > 1960, SexID == 2, complete == TRUE) %>%
    select(TimeLabel, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID, DataSourceName) %>%
    distinct() %>%
    arrange(TimeLabel) %>%
    pivot_wider(names_from = TimeLabel, values_from = DataValue) %>%
    group_by(AgeStart, AgeLabel, AgeSpan, DataSourceName) %>%
    arrange(AgeStart, StatisticalConceptName) %>% ###De-facto first
    select(-c(StatisticalConceptName, SexID)) %>%
    summarise_all(function(y){first(na.omit(y))}) %>%
    ungroup() %>%
    select(-AgeSpan, -AgeStart) %>%
    mutate(DataSourceName = replace(DataSourceName, DataSourceName == "IPUMS International custom tabulations", "IPUMS"),
           DataSourceName = replace(DataSourceName, DataSourceName == "Demographic Yearbook", "DYB"),
           country = country,
           age = AgeLabel,
           .before = 3) %>%
    select(-AgeLabel)
  
  ddharm_bf_census_m_5 <- census_pop_counts %>%  
    filter(AgeSpan %in% c(-1, 5), AgeLabel != "Total", TimeLabel > 1960, SexID == 1, five_year == TRUE) %>%
    select(TimeLabel, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID, DataSourceName) %>%
    distinct() %>%
    arrange(TimeLabel) %>%
    pivot_wider(names_from = TimeLabel, values_from = DataValue) %>%
    group_by(AgeStart, AgeLabel, AgeSpan, DataSourceName) %>%
    arrange(AgeStart, StatisticalConceptName) %>% ###De-facto first
    select(-c(StatisticalConceptName, SexID)) %>%
    summarise_all(function(y){first(na.omit(y))}) %>%
    ungroup() %>%
    select(-AgeSpan, -AgeStart) %>%
    mutate(DataSourceName = replace(DataSourceName, DataSourceName == "IPUMS International custom tabulations", "IPUMS"),
           DataSourceName = replace(DataSourceName, DataSourceName == "Demographic Yearbook", "DYB"),
           country = country,
           age = AgeLabel,
           .before = 3) %>%
    select(-AgeLabel)
  
  ddharm_bf_census_f_5 <- census_pop_counts %>%  
    filter(AgeSpan %in% c(-1, 5), AgeLabel != "Total", TimeLabel > 1960, SexID == 2, five_year == TRUE) %>%
    select(TimeLabel, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID, DataSourceName) %>%
    distinct() %>%
    arrange(TimeLabel) %>%
    pivot_wider(names_from = TimeLabel, values_from = DataValue) %>%
    group_by(AgeStart, AgeLabel, AgeSpan, DataSourceName) %>%
    arrange(AgeStart, StatisticalConceptName) %>% ###De-facto first
    select(-c(StatisticalConceptName, SexID)) %>%
    summarise_all(function(y){first(na.omit(y))}) %>%
    ungroup() %>%
    select(-AgeSpan, -AgeStart) %>%
    mutate(DataSourceName = replace(DataSourceName, DataSourceName == "IPUMS International custom tabulations", "IPUMS"),
           DataSourceName = replace(DataSourceName, DataSourceName == "Demographic Yearbook", "DYB"),
           country = country,
           age = AgeLabel,
           .before = 3) %>%
    select(-AgeLabel)
  } else {
    ddharm_bf_census_f <- ddharm_bf_census_m <- ddharm_bf_census_f_5 <- ddharm_bf_census_m_5 <- tibble(DataSourceName = character(),
                                                                                                       country = character(),
                                                                                                       age = numeric(),
                                                                                                       DataValue = numeric())
  }
  
  if(exists("census_deaths") && !is.null(census_deaths)){
  ddharm_bf_census_deaths_m <- census_deaths %>%  
    filter(AgeSpan %in% c(-1, 1), AgeLabel != "Total", TimeLabel > 1960, SexID == 1, complete == TRUE) %>%
    select(TimeLabel, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID, DataSourceName) %>%
    distinct() %>%
    arrange(TimeLabel) %>%
    pivot_wider(names_from = TimeLabel, values_from = DataValue) %>%
    group_by(AgeStart, AgeLabel, AgeSpan, DataSourceName) %>%
    arrange(AgeStart, StatisticalConceptName) %>% ###De-factor first
    select(-c(StatisticalConceptName, SexID)) %>%
    summarise_all(function(y){first(na.omit(y))}) %>%
    ungroup() %>%
    select(-AgeSpan, -AgeLabel) %>%
    mutate(aggr.age = ifelse(AgeStart >= open.age, open.age, AgeStart)) %>%
    group_by(aggr.age, DataSourceName) %>%
    summarise_at(vars(-AgeStart), function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
    rename(age = aggr.age) %>%
    mutate(DataSourceName = replace(DataSourceName, DataSourceName == "IPUMS International custom tabulations", "IPUMS"),
           DataSourceName = replace(DataSourceName, DataSourceName == "Demographic Yearbook", "DYB"),
           country = country,
           .before = 1) 
  
  ddharm_bf_census_deaths_f <- census_deaths %>%  
    filter(AgeSpan %in% c(-1, 1), AgeLabel != "Total", TimeLabel > 1960, SexID == 2, complete == TRUE) %>%
    select(TimeLabel, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID, DataSourceName) %>%
    distinct() %>%
    arrange(TimeLabel) %>%
    pivot_wider(names_from = TimeLabel, values_from = DataValue) %>%
    group_by(AgeStart, AgeLabel, AgeSpan, DataSourceName) %>%
    arrange(AgeStart, StatisticalConceptName) %>% ###De-factor first
    select(-c(StatisticalConceptName, SexID)) %>%
    summarise_all(function(y){first(na.omit(y))}) %>%
    ungroup() %>%
    select(-AgeSpan, -AgeLabel) %>%
    mutate(aggr.age = ifelse(AgeStart >= open.age, open.age, AgeStart)) %>%
    group_by(aggr.age, DataSourceName) %>%
    summarise_at(vars(-AgeStart), function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
    rename(age = aggr.age) %>%
    mutate(DataSourceName = replace(DataSourceName, DataSourceName == "IPUMS International custom tabulations", "IPUMS"),
           DataSourceName = replace(DataSourceName, DataSourceName == "Demographic Yearbook", "DYB"),
           country = country,
           .before = 1) 
  
  ddharm_bf_census_deaths_m_5 <- census_deaths %>%  
    filter(AgeSpan %in% c(-1, 5), AgeLabel != "Total", TimeLabel > 1960, SexID == 1, five_year == TRUE) %>%
    select(TimeLabel, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID, DataSourceName) %>%
    distinct() %>%
    arrange(TimeLabel) %>%
    pivot_wider(names_from = TimeLabel, values_from = DataValue) %>%
    group_by(AgeStart, AgeLabel, AgeSpan, DataSourceName) %>%
    arrange(AgeStart, StatisticalConceptName) %>% ###De-factor first
    select(-c(StatisticalConceptName, SexID)) %>%
    summarise_all(function(y){first(na.omit(y))}) %>%
    ungroup() %>%
    select(-AgeSpan, -AgeLabel) %>%
    mutate(aggr.age = ifelse(AgeStart >= open.age, open.age, AgeStart)) %>%
    group_by(aggr.age, DataSourceName) %>%
    summarise_at(vars(-AgeStart), function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
    rename(age = aggr.age) %>%
    mutate(DataSourceName = replace(DataSourceName, DataSourceName == "IPUMS International custom tabulations", "IPUMS"),
           DataSourceName = replace(DataSourceName, DataSourceName == "Demographic Yearbook", "DYB"),
           country = country,
           .before = 1) 
  
  ddharm_bf_census_deaths_f_5 <- census_deaths %>%  
    filter(AgeSpan %in% c(-1, 5), AgeLabel != "Total", TimeLabel > 1960, SexID == 2, five_year == TRUE) %>%
    select(TimeLabel, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID, DataSourceName) %>%
    distinct() %>%
    arrange(TimeLabel) %>%
    pivot_wider(names_from = TimeLabel, values_from = DataValue) %>%
    group_by(AgeStart, AgeLabel, AgeSpan, DataSourceName) %>%
    arrange(AgeStart, StatisticalConceptName) %>% ###De-factor first
    select(-c(StatisticalConceptName, SexID)) %>%
    summarise_all(function(y){first(na.omit(y))}) %>%
    ungroup() %>%
    select(-AgeSpan, -AgeLabel) %>%
    mutate(aggr.age = ifelse(AgeStart >= open.age, open.age, AgeStart)) %>%
    group_by(aggr.age, DataSourceName) %>%
    summarise_at(vars(-AgeStart), function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
    rename(age = aggr.age) %>%
    mutate(DataSourceName = replace(DataSourceName, DataSourceName == "IPUMS International custom tabulations", "IPUMS"),
           DataSourceName = replace(DataSourceName, DataSourceName == "Demographic Yearbook", "DYB"),
           country = country,
           .before = 1) 
  } else {
    ddharm_bf_census_deaths_f <- ddharm_bf_census_deaths_m <- ddharm_bf_census_deaths_f_5 <- ddharm_bf_census_deaths_m_5 <- tibble(DataSourceName = character(),
                                                                                                                                   country = character(),
                                                                                                                                   age = numeric(),
                                                                                                                                   DataValue = numeric())
    }
  
  pop.m.list[[country]] <- ddharm_bf_census_m
  pop.f.list[[country]] <- ddharm_bf_census_f
  pop.m.5.list[[country]] <- ddharm_bf_census_m_5
  pop.f.5.list[[country]] <- ddharm_bf_census_f_5
  
  deaths.m.list[[country]] <- ddharm_bf_census_deaths_m
  deaths.f.list[[country]] <- ddharm_bf_census_deaths_f
  deaths.m.5.list[[country]] <- ddharm_bf_census_deaths_m_5
  deaths.f.5.list[[country]] <- ddharm_bf_census_deaths_f_5
  
  census_pop_counts <- NULL
  census_deaths <- NULL
}

pop.m.list$Namibia %>% pivot_longer(cols = -(1:3)) %>%
  mutate(age = str_extract(age, "\\d*")) %>%
  drop_na(value) %>%
  pivot_wider(names_from = name, values_from = value) %>% view

save(pop.f.list, pop.m.list, pop.f.5.list, pop.m.5.list, deaths.f.list, deaths.f.5.list, deaths.m.list, deaths.m.5.list, file = "Census pop and deaths all countries 21-04-2023.RData")
