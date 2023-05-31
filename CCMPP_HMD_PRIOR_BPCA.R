params <- list(
  country = "Namibia", 
  age.knot.space = 2.5,
  year.knot.space = 2.5,
  supersmooth.pop = FALSE,
  download = FALSE,
  br.data = TRUE,
  census.deaths = FALSE,
  loghump = TRUE
)

# setwd("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_normalhump_trend_and_AR2")
setwd("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/CCMPP HMD BPCA")
load("~/more countries final avg sex Rwanda.RData")
# load("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/Census pop and deaths all countries.Rdata")
load("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/Census pop and deaths all countries 19-01-2023.RData")
load("C:/Users/ktang3/Desktop/Imperial/Life_Table_System/child mortality DHS BR.Rdata")

library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(tmbstan)
library(readr)
library(popReconstruct)
library(tidyverse)
library(TMB)
library(ggnewscale)
library(DDSQLtools)
library(ggforce)
library(stringr)
library(Matrix)
library(magic)
library(Hmisc)
library(gridExtra)
library(ddharmony)
library(censusAdjust)
library(DemoTools)
library(Rfast)

# require(devtools)
# require(httr)
# 
compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_bothsexes_thiele_normalhump_HMD_trend_and_err_BPCA.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_bothsexes_thiele_normalhump_HMD_trend_and_err_BPCA"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_bothsexes_thiele_normalhump_HMD_trend_and_err_BPCA_ARIMAhump.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_bothsexes_thiele_normalhump_HMD_trend_and_err_BPCA_ARIMAhump"))

# req <- GET("https://api.github.com/repos/sarahertog/ddharmony/git/trees/main?recursive=1")
# stop_for_status(req)
# filelist <- unlist(lapply(content(req)$tree, "[", "path"), use.names = F)
# 
# for (filename in filelist) {
#   one_function <- paste0("https://github.com/sarahertog/ddharmony/blob/main/", filename, "?raw=TRUE")
#   source_url(one_function)
#   rm(one_function)
# }
# rm(req, filelist, filename)

pyears_data <- readRDS("C:/Users/ktang3/Desktop/Imperial/SSA_mort/JE/data/pyears_data_smoothed.rds")
aggr.mat.cohort.0 <- pyears_data %>%
  rename(DHS = SurveyId) %>%
  group_by(country) %>% group_split %>%
  setNames(pyears_data$country %>% levels) 

projection_indices <- function(period_start,  period_end, interval, n_ages,
                               fx_idx, n_fx, n_sexes = 1) {
  
  stopifnot(n_sexes %in% 1:2)
  
  periods_out <- seq(period_start, period_end + interval, by = interval)
  periods <- periods_out[-length(periods_out)]
  ages <- 0:(n_ages - 1) * interval
  fertility_ages <- ages[fx_idx + 0:(n_fx - 1L)]
  sexes <- if(n_sexes == 1) "female" else c("female", "male")
  
  list(periods = periods,
       periods_out = periods_out,
       ages = ages,
       fertility_ages = fertility_ages,
       sexes = sexes,
       interval = interval,
       n_periods = length(periods),
       n_ages = n_ages,
       fx_idx = fx_idx,
       n_fx = n_fx,
       n_sexes = n_sexes)
}
projection_model_frames <- function(indices) {
  
  mf_population <- tidyr::crossing(period = indices$periods_out,
                                   sex = indices$sexes,
                                   age = indices$ages)
  
  mf_deaths <- tidyr::crossing(period = indices$periods,
                               sex = indices$sexes,
                               age = indices$ages)
  
  cohort_death_ages <- c(indices$ages, max(indices$ages) + indices$interval)
  mf_cohort_deaths <- tidyr::crossing(period = indices$periods,
                                      sex = indices$sexes,
                                      age = cohort_death_ages)
  
  mf_migrations <- tidyr::crossing(period = indices$periods,
                                   sex = indices$sexes,
                                   age = indices$ages)
  
  mf_births <- tidyr::crossing(period = indices$periods,
                               age = indices$fertility_ages)
  
  list(mf_population = mf_population,
       mf_deaths = mf_deaths,
       mf_cohort_deaths = mf_cohort_deaths,
       mf_migrations = mf_migrations,
       mf_births = mf_births)
}

proj_to_long <- function(proj, mf, value_col = "value") {
  
  val <- list()
  val$population <- mf$mf_population
  val$population[[value_col]] <- as.vector(proj$population)
  
  val$population <- mf$mf_cohort_deaths
  val$population[[value_col]] <- as.vector(proj$cohort_deathsg)
  
}

fit_tmb <- function(data,
                    par, 
                    outer_verbose = TRUE,
                    inner_verbose = FALSE,
                    random,
                    DLL = "leapfrog_TMBExports",
                    map = list(),
                    iter.max = 10000,
                    eval.max = 10000){
  
  obj <- TMB::MakeADFun(data = data,
                        par = par,
                        silent = !inner_verbose,
                        random = random,
                        DLL = DLL,
                        map = map)
  
  trace <- if(outer_verbose) 1 else 0
  f <- withCallingHandlers(
    stats::nlminb(obj$par, obj$fn, obj$gr,
                  control = list(trace = trace,
                                 rel.tol = 1e-8,
                                 iter.max = iter.max,
                                 eval.max = 10000)),
    warning = function(w) {
      if(grepl("NA/NaN function evaluation", w$message))
        invokeRestart("muffleWarning")
    }
  )
  
  if(f$convergence != 0)
    warning(paste("convergence error:", f$message))
  
  if(outer_verbose)
    message(paste("converged:", f$message))
  
  f$par.fixed <- f$par
  f$par.full <- obj$env$last.par
  
  objout <- TMB::MakeADFun(data = data,
                           par = par,
                           silent = !inner_verbose,
                           random = random,
                           DLL = DLL,
                           map = map)
  
  f$mode <- objout$report(f$par.full)
  
  val <- c(f, obj = list(objout))
  
  val
}

igme.5q0.m<-read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/5q0 male.csv")
igme.5q0.f<-read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/5q0 female.csv")
wpp.fx <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/WPP fx.csv")
#wpp.pop <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/WPP Pop estimates.csv")
wpp.pop.age.specific <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/WPP pop.csv")
#wpp.q4515 <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/WPP 45q15.csv")
#gbd.q4515 <- read.csv(file="C:/Users/ktang3/Desktop/Imperial/SSA_mort/GBD 45q15.csv",header=T)
wpp.qx <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/WPP age specific.csv")
#gbd.qx <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/GBD age specific.csv")

#igme.5q0.hivfree <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/Results (crisis-and-HIV-free)_u5mr.csv")
# 
# gbd.q4515$location_name<-str_replace(gbd.q4515$location_name,"Democratic Republic of the Congo","Congo Democratic Republic")
# gbd.q4515$location_name<-str_replace(gbd.q4515$location_name,"Côte d'Ivoire","Cote d'Ivoire")
# gbd.q4515$location_name<-str_replace(gbd.q4515$location_name,"United Republic of Tanzania","Tanzania")

igme.5q0.m$Country.Name <- str_replace(igme.5q0.m$Country.Name,"Democratic Republic of the Congo","Congo Democratic Republic")
igme.5q0.f$Country.Name <- str_replace(igme.5q0.f$Country.Name,"Democratic Republic of the Congo","Congo Democratic Republic")
igme.5q0.m$Country.Name<-str_replace(igme.5q0.m$Country.Name,"United Republic of Tanzania","Tanzania")
igme.5q0.f$Country.Name<-str_replace(igme.5q0.f$Country.Name,"United Republic of Tanzania","Tanzania")

wpp.qx$name <- str_replace(wpp.qx$name,"Democratic Republic of the Congo","Congo Democratic Republic")
wpp.fx$Name <- str_replace(wpp.fx$Name,"Democratic Republic of the Congo","Congo Democratic Republic")
wpp.qx$name<-str_replace(wpp.qx$name,"United Republic of Tanzania","Tanzania")
wpp.fx$Name<-str_replace(wpp.fx$Name,"United Republic of Tanzania","Tanzania")

# GBD.age <- gbd.qx %>%
#   dplyr::filter(sex != "both", 
#          measure_name == "Probability of death", 
#          name %in% joint.countries) %>%
#   mutate(
#     age = as.numeric(str_extract(age_group_name, "\\d+")) + 2)
# 
# GBD.mx <- GBD.age %>%
#   select(name, age, sex, year, val) %>%
#   mutate(sex = replace(sex, sex == "female", "f"),
#          sex = replace(sex, sex == "male", "m"),
#          GBD = val / (5 - 2.5 * val)) %>%
#   select(-val)

# WPP.mx <- wpp.qx %>%
#   filter(Sex != "Total") %>%
#   mutate(year = MidPeriod - 1,
#          sex = replace(Sex, Sex == "Male", "m"),
#          sex = replace(sex, Sex == "Female", "f"),
#          age = AgeGrpStart + 0.5 * (AgeGrpSpan - 1)) %>%
#   select(name, year, sex, age, mx) %>%
#   rename(WPP = mx)
# 
# WPP.mx.expand <- bind_rows(WPP.mx,
#                            WPP.mx %>% mutate(year = year - 2),
#                            WPP.mx %>% mutate(year = year - 1),
#                            WPP.mx %>% mutate(year = year + 1),
#                            WPP.mx %>% mutate(year = year + 2))  

open.age <- 85
n_ages <- open.age + 1
interval <- 1 #just in case

country <- params$country

bspline<-function (x,k,i,m=2) {
  if (m==-1) {basis<-as.numeric(x<k[i+1] & x>=k[i])} else {
    z0<-(x-k[i])/(k[i+m+1]-k[i])
    z1<-(k[i+m+2]-x)/(k[i+m+2]-k[i+1])
    basis<-z0*bspline(x,k,i,m-1)+z1*bspline(x,k,i+1,m-1) }
  basis
}

# library(MortCast)
# LQcoef.f <- LQcoef %>% filter(sex=="Female", !age%in%c("0","1-4")) %>% select(ax:vx)
# LQcoef.m <- LQcoef %>% filter(sex=="Male", !age%in%c("0","1-4")) %>% select(ax:vx)

##IGME 5q0
igme.5q0.df <- igme.5q0.f %>% filter(Country.Name == country) %>% reshape2::melt() %>% 
  mutate(year = as.numeric(gsub("X|\\.5","",variable)), child.mort = value/1000, Sex = "Female") %>%
  select(year,child.mort, Sex) %>%
  bind_rows(
    igme.5q0.m %>% filter(Country.Name == country) %>% reshape2::melt() %>% 
      mutate(year = as.numeric(gsub("X|\\.5","",variable)), child.mort = value/1000, Sex = "Male") %>%
      select(year,child.mort, Sex)
  )


#WPP 5q0
wpp.bf.qx <- wpp.qx %>% filter(name==country, Sex!="Total")

wpp.5q0 <- wpp.bf.qx %>% filter(AgeGrpStart %in% 0:1) %>% select(c(MidPeriod, Sex, px)) %>%
  group_by(MidPeriod, Sex) %>%
  summarise_at(vars(px),prod) %>%
  mutate(q50 = 1-px,
         year5 = MidPeriod - 3) %>%
  ungroup() %>% arrange(Sex, year5) %>%
  select(c(Sex, year5, q50))

wpp.5q0.interpolate <- bind_rows(
  approx(x = seq(1952, 2097, by = 5), y = wpp.5q0 %>% filter(Sex=="Female") %>% .$q50, xout=1952:2097)$y %>%
    as_tibble() %>%
    mutate(q50 = value,
           Sex = "Female",
           year = 1952:2097,
           .keep = "none"),
  
  approx(x = seq(1952, 2097, by = 5), y = wpp.5q0 %>% filter(Sex=="Male") %>% .$q50, xout=1952:2097)$y %>%
    as_tibble() %>%
    mutate(q50 = value,
           Sex = "Male",
           year = 1952:2097,
           .keep = "none")
)

if(params$download){
  locationid <- get_locations()
  #DDHarmonized smoothed
  try(
    census_pop_counts <- DDharmonize_validate_PopCounts(locid = ifelse(country=="Cote d'Ivoire", 384,
                                                                       ifelse(country=="Tanzania",834,
                                                                              locationid$PK_LocID[which(locationid$Name == country)])),
                                                        times = 1950:2020,
                                                        DataSourceShortName = "DYB") # time frame for censuses to extract from Demographic Yearbook
  )
  
  try(
    census_deaths <- DDharmonize_validate_DeathCounts(locid = ifelse(country=="Cote d'Ivoire", 384,
                                                                     ifelse(country=="Tanzania",834,
                                                                            locationid$PK_LocID[which(locationid$Name == country)])),
                                                      times = 1950:2020,
                                                      process="census",
                                                      DataSourceShortName = "DYB",
                                                      retainKeys = TRUE) # time frame for censuses to extract from Demographic Yearbook
  )
  
  ddharm_bf_census_m <- census_pop_counts %>%  
    filter(AgeSpan %in% c(-1, 1), AgeLabel != "Total", TimeLabel > 1960, SexID == 1, complete == TRUE) %>%
    select(TimeLabel, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID) %>%
    distinct() %>%
    arrange(TimeLabel) %>%
    pivot_wider(names_from = TimeLabel, values_from = DataValue) %>%
    group_by(AgeStart, AgeLabel, AgeSpan) %>%
    arrange(AgeStart, StatisticalConceptName) %>% ###De-facto first
    select(-c(StatisticalConceptName, SexID)) %>%
    summarise_all(function(y){first(na.omit(y))}) %>%
    select(-AgeSpan) %>%
    ungroup()
  
  ddharm_bf_census_f <- census_pop_counts %>%  
    filter(AgeSpan %in% c(-1, 1), AgeLabel != "Total", TimeLabel > 1960, SexID == 2, complete == TRUE) %>%
    select(TimeLabel, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID) %>%
    distinct() %>%
    arrange(TimeLabel) %>%
    pivot_wider(names_from = TimeLabel, values_from = DataValue) %>%
    group_by(AgeStart, AgeLabel, AgeSpan) %>%
    arrange(AgeStart, StatisticalConceptName) %>% ###De-factor first
    select(-c(StatisticalConceptName, SexID)) %>%
    summarise_all(function(y){first(na.omit(y))}) %>%
    select(-AgeSpan) %>%
    ungroup()
  
  ddharm_bf_census_m_5 <- census_pop_counts %>%  
    filter(AgeSpan %in% c(-1, 5), AgeLabel != "Total", TimeLabel > 1960, SexID == 1, five_year == TRUE) %>%
    select(TimeLabel, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID) %>%
    distinct() %>%
    arrange(TimeLabel) %>%
    pivot_wider(names_from = TimeLabel, values_from = DataValue) %>%
    group_by(AgeStart, AgeLabel, AgeSpan) %>%
    arrange(AgeStart, StatisticalConceptName) %>% ###De-factor first
    select(-c(StatisticalConceptName, SexID)) %>%
    summarise_all(function(y){first(na.omit(y))}) %>%
    select(-AgeSpan) %>%
    ungroup()
  
  ddharm_bf_census_f_5 <- census_pop_counts %>%  
    filter(AgeSpan %in% c(-1, 5), AgeLabel != "Total", TimeLabel > 1960, SexID == 2, five_year == TRUE) %>%
    select(TimeLabel, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID) %>%
    distinct() %>%
    arrange(TimeLabel) %>%
    pivot_wider(names_from = TimeLabel, values_from = DataValue) %>%
    group_by(AgeStart, AgeLabel, AgeSpan) %>%
    arrange(AgeStart, StatisticalConceptName) %>% ###De-factor first
    select(-c(StatisticalConceptName, SexID)) %>%
    summarise_all(function(y){first(na.omit(y))}) %>%
    select(-AgeSpan) %>%
    ungroup()
  
  if(exists("census_deaths") && !is.null(census_deaths)){
    ddharm_bf_census_deaths_m <- census_deaths %>%  
      filter(AgeSpan %in% c(-1, 1), AgeLabel != "Total", TimeLabel > 1960, SexID == 1, complete == TRUE) %>%
      select(TimeLabel, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID) %>%
      distinct() %>%
      arrange(TimeLabel) %>%
      pivot_wider(names_from = TimeLabel, values_from = DataValue) %>%
      group_by(AgeStart, AgeLabel, AgeSpan) %>%
      arrange(AgeStart, StatisticalConceptName) %>% ###De-factor first
      select(-c(StatisticalConceptName, SexID)) %>%
      summarise_all(function(y){first(na.omit(y))}) %>%
      ungroup() %>%
      select(-AgeSpan, -AgeLabel) %>%
      mutate(aggr.age = ifelse(AgeStart >= open.age, open.age, AgeStart)) %>%
      group_by(aggr.age) %>%
      summarise_at(vars(-AgeStart), function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
      rename(age = aggr.age)
    
    ddharm_bf_census_deaths_f <- census_deaths %>%  
      filter(AgeSpan %in% c(-1, 1), AgeLabel != "Total", TimeLabel > 1960, SexID == 2, complete == TRUE) %>%
      select(TimeLabel, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID) %>%
      distinct() %>%
      arrange(TimeLabel) %>%
      pivot_wider(names_from = TimeLabel, values_from = DataValue) %>%
      group_by(AgeStart, AgeLabel, AgeSpan) %>%
      arrange(AgeStart, StatisticalConceptName) %>% ###De-factor first
      select(-c(StatisticalConceptName, SexID)) %>%
      summarise_all(function(y){first(na.omit(y))}) %>%
      ungroup() %>%
      select(-AgeSpan, -AgeLabel) %>%
      mutate(aggr.age = ifelse(AgeStart >= open.age, open.age, AgeStart)) %>%
      group_by(aggr.age) %>%
      summarise_at(vars(-AgeStart), function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
      rename(age = aggr.age)
    
    ddharm_bf_census_deaths_m_5 <- census_deaths %>%  
      filter(AgeSpan %in% c(-1, 5), AgeLabel != "Total", TimeLabel > 1960, SexID == 1, five_year == TRUE) %>%
      select(TimeLabel, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID) %>%
      distinct() %>%
      arrange(TimeLabel) %>%
      pivot_wider(names_from = TimeLabel, values_from = DataValue) %>%
      group_by(AgeStart, AgeLabel, AgeSpan) %>%
      arrange(AgeStart, StatisticalConceptName) %>% ###De-factor first
      select(-c(StatisticalConceptName, SexID)) %>%
      summarise_all(function(y){first(na.omit(y))}) %>%
      ungroup() %>%
      select(-AgeSpan, -AgeLabel) %>%
      mutate(aggr.age = ifelse(AgeStart >= open.age, open.age, AgeStart)) %>%
      group_by(aggr.age) %>%
      summarise_at(vars(-AgeStart), function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
      rename(age = aggr.age)
    
    ddharm_bf_census_deaths_f_5 <- census_deaths %>%  
      filter(AgeSpan %in% c(-1, 5), AgeLabel != "Total", TimeLabel > 1960, SexID == 2, five_year == TRUE) %>%
      select(TimeLabel, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID) %>%
      distinct() %>%
      arrange(TimeLabel) %>%
      pivot_wider(names_from = TimeLabel, values_from = DataValue) %>%
      group_by(AgeStart, AgeLabel, AgeSpan) %>%
      arrange(AgeStart, StatisticalConceptName) %>% ###De-factor first
      select(-c(StatisticalConceptName, SexID)) %>%
      summarise_all(function(y){first(na.omit(y))}) %>%
      ungroup() %>%
      select(-AgeSpan, -AgeLabel) %>%
      mutate(aggr.age = ifelse(AgeStart >= open.age, open.age, AgeStart)) %>%
      group_by(aggr.age) %>%
      summarise_at(vars(-AgeStart), function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
      rename(age = aggr.age)
  } else {
    ddharm_bf_census_deaths_f <- ddharm_bf_census_deaths_m <- ddharm_bf_census_deaths_f_5 <- ddharm_bf_census_deaths_m_5 <- NULL
  }
  
} else {
  ddharm_bf_census_m <- if(!"DYB" %in% unique(pop.m.list[[country]]$DataSourceName)) tibble() else {pop.m.list[[country]] %>% 
      filter(DataSourceName == "DYB") %>%
      .[, colSums(is.na(.)) != nrow(.)] %>%
      select(-DataSourceName, -country) %>%
      pivot_longer(cols = -1) %>%
      mutate(age = str_extract(age, "\\d*")) %>%
      drop_na(value) %>%
      group_by(name) %>%
      mutate(age = as.numeric(age),
             drop = ifelse(max(age) < open.age, NA, 1)) %>%
      ungroup %>% drop_na(drop) %>% select(-drop) %>%
      pivot_wider(names_from = name, values_from = value) %>%
      rename(AgeStart = age) %>%
      mutate(AgeLabel = replace(AgeStart, AgeStart == max(AgeStart), paste0(max(AgeStart),"+")), .after = 1)
  }
  
  ddharm_bf_census_f <- if(!"DYB" %in% unique(pop.f.list[[country]]$DataSourceName)) tibble() else {pop.f.list[[country]] %>% 
      filter(DataSourceName == "DYB") %>%
      .[, colSums(is.na(.)) != nrow(.)] %>%
      select(-DataSourceName, -country) %>%
      pivot_longer(cols = -1) %>%
      mutate(age = str_extract(age, "\\d*")) %>%
      drop_na(value) %>%
      group_by(name) %>%
      mutate(age = as.numeric(age),
             drop = ifelse(max(age) < open.age, NA, 1)) %>%
      ungroup %>% drop_na(drop) %>% select(-drop) %>%
      pivot_wider(names_from = name, values_from = value) %>%
      rename(AgeStart = age) %>%
      mutate(AgeLabel = replace(AgeStart, AgeStart == max(AgeStart), paste0(max(AgeStart),"+")), .after = 1)
  }
  
  ddharm_bf_census_m_5 <- if(!"DYB" %in% unique(pop.m.5.list[[country]]$DataSourceName)) tibble() else {pop.m.5.list[[country]] %>% 
      filter(DataSourceName == "DYB") %>%
      .[, colSums(is.na(.)) != nrow(.)] %>%
      select(-DataSourceName, -country) %>%
      pivot_longer(cols = -1) %>%
      mutate(age = str_extract(age, "\\d*")) %>%
      drop_na(value) %>%
      group_by(name) %>%
      mutate(age = as.numeric(age),
             drop = ifelse(max(age) < open.age, NA, 1)) %>%
      ungroup %>% drop_na(drop) %>% select(-drop) %>%
      pivot_wider(names_from = name, values_from = value) %>%
      rename(AgeStart = age) %>%
      mutate(AgeLabel = paste0(AgeStart, "-", AgeStart + 4), .after = 1)
  }
  
  ddharm_bf_census_f_5 <- if(!"DYB" %in% unique(pop.f.5.list[[country]]$DataSourceName)) tibble() else {pop.f.5.list[[country]] %>% 
      filter(DataSourceName == "DYB") %>%
      .[, colSums(is.na(.)) != nrow(.)] %>%
      select(-DataSourceName, -country) %>%
      pivot_longer(cols = -1) %>%
      mutate(age = str_extract(age, "\\d*")) %>%
      drop_na(value) %>%
      group_by(name) %>%
      mutate(age = as.numeric(age),
             drop = ifelse(max(age) < open.age, NA, 1)) %>%
      ungroup %>% drop_na(drop) %>% select(-drop) %>%
      pivot_wider(names_from = name, values_from = value) %>%
      rename(AgeStart = age) %>%
      mutate(AgeLabel = paste0(AgeStart, "-", AgeStart + 4), .after = 1)
  }
  
  ddharm_bf_census_deaths_m <- if(is.null(deaths.m.list[[country]])) NULL else{deaths.m.list[[country]] %>%
      filter(DataSourceName == "DYB") %>%
      select(-DataSourceName, -country) %>%
      .[, colSums(is.na(.)) != nrow(.)]}
  
  ddharm_bf_census_deaths_f <- if(is.null(deaths.f.list[[country]])) NULL else{deaths.f.list[[country]] %>%
      filter(DataSourceName == "DYB") %>%
      select(-DataSourceName, -country) %>%
      .[, colSums(is.na(.)) != nrow(.)]}
  
  ddharm_bf_census_deaths_m_5 <- if(is.null(deaths.m.5.list[[country]])) NULL else{deaths.m.5.list[[country]] %>%
      filter(DataSourceName == "DYB") %>%
      select(-DataSourceName, -country) %>%
      .[, colSums(is.na(.)) != nrow(.)]}
  
  ddharm_bf_census_deaths_f_5 <- if(is.null(deaths.f.5.list[[country]])) NULL else{deaths.f.5.list[[country]] %>%
      filter(DataSourceName == "DYB") %>%
      select(-DataSourceName, -country) %>%
      .[, colSums(is.na(.)) != nrow(.)]}
}
#####################MIXING DE-FACTO AND DE-JURE HERE
if(nrow(ddharm_bf_census_m) != 0 & nrow(ddharm_bf_census_f) != 0){
  ddharm_smoothed <- tibble()
  
  for(i in 3:ncol(ddharm_bf_census_f)){
    dat.f <- ddharm_bf_census_f %>% select(1:2,i) %>% filter(!is.na(ddharm_bf_census_f[[i]])) %>% arrange(AgeStart)
    dat.m <- ddharm_bf_census_m %>% select(1:2,i) %>% filter(!is.na(ddharm_bf_census_m[[i]])) %>% arrange(AgeStart)
    
    stopifnot(length(unique(diff(dat.f$AgeStart))) == 1 | length(unique(diff(dat.m$AgeStart))) == 1)
    
    #census_workflow_adjust_smooth
    # 
    # smoothed.adult <- getSmoothedPop1(popM = dat.m[[3]], popF = dat.f[[3]], Age = dat.m[[1]], EduYrs = 2, subgroup = "adult")
    # smoothed.child <- getSmoothedPop1(popM = dat.m[[3]], popF = dat.f[[3]], Age = dat.m[[1]], EduYrs = 2, subgroup = "child")
    # 
    # adult.grad <- as.numeric(str_extract(smoothed.adult$best_smooth_method,"\\s\\d"))
    # child.grad <- as.numeric(str_extract(smoothed.child$best_smooth_method,"\\s\\d"))
    # 
    # cat("adult bestGrad5 =",adult.grad,"; child bestGrad5 =",child.grad, "\n")
    # 
    # if(adult.grad >= child.grad){
    #   ddharm_smoothed <- bind_rows(
    #     ddharm_smoothed,
    #     as_tibble(smoothed.adult$popM_smooth) %>% mutate(age = dat.m[[1]], AgeLabel = dat.m[[2]], 
    #                                                      sex="male", year = as.numeric(names(ddharm_bf_census_f)[i])),
    #     as_tibble(smoothed.adult$popF_smooth) %>% mutate(age = dat.f[[1]], AgeLabel = dat.f[[2]],
    #                                                      sex="female", year = as.numeric(names(ddharm_bf_census_m)[i]))
    #   )} else{
    #     ddharm_smoothed <- bind_rows(
    #       ddharm_smoothed,
    #       as_tibble(smoothed.child$popM_smooth) %>% mutate(age = dat.m[[1]], AgeLabel = dat.m[[2]], 
    #                                                        sex="male", year = as.numeric(names(ddharm_bf_census_f)[i])),
    #       as_tibble(smoothed.child$popF_smooth) %>% mutate(age = dat.f[[1]], AgeLabel = dat.f[[2]],
    #                                                        sex="female", year = as.numeric(names(ddharm_bf_census_m)[i]))
    #     )
    #   }
    
    smoothed <- census_workflow_adjust_smooth(popM = dat.m[[3]], popF = dat.f[[3]], Age = dat.m[[1]], EduYrs = 2)
    ddharm_smoothed <- bind_rows(
      ddharm_smoothed,
      as_tibble(smoothed$popM_smooth) %>% mutate(age = dat.m[[1]], AgeLabel = dat.m[[2]],
                                                 sex="male", year = as.numeric(names(ddharm_bf_census_m)[i])),
      as_tibble(smoothed$popF_smooth) %>% mutate(age = dat.f[[1]], AgeLabel = dat.f[[2]],
                                                 sex="female", year = as.numeric(names(ddharm_bf_census_f)[i]))
    )
    
  }
  
  ddharm_smoothed_mat <- ddharm_smoothed %>% pivot_wider(names_from = year, values_from = value) %>% arrange(sex, age)
  
  ddharm_census_f_smoothed <- ddharm_smoothed_mat %>% filter(sex=="female") %>% select(-sex)
  ddharm_census_m_smoothed <- ddharm_smoothed_mat %>% filter(sex=="male") %>% select(-sex)
  
  ddharm_bf_census_f_oag <- ddharm_census_f_smoothed %>% 
    #rename(age = AgeStart) %>%
    group_by(age) %>% 
    summarise_at(vars(-AgeLabel), function(y){first(na.omit(y))}) %>%
    mutate(aggr.age = ifelse(age >= open.age, open.age, age)) %>%
    group_by(aggr.age) %>%
    summarise_at(vars(-age), function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
    mutate(age = aggr.age) %>%
    group_by(age) %>%
    summarise_all(function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
    select(-2) 
  
  ddharm_bf_census_m_oag <- ddharm_census_m_smoothed %>% 
    #rename(age = AgeStart) %>%
    group_by(age) %>%
    summarise_at(vars(-AgeLabel), function(y){first(na.omit(y))}) %>%
    mutate(aggr.age = ifelse(age >= open.age, open.age, age)) %>%
    group_by(aggr.age) %>%
    summarise_at(vars(-age), function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
    mutate(age = aggr.age) %>%
    group_by(age) %>%
    summarise_all(function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
    select(-2) 
}

##census data 5 year age groups
if(nrow(ddharm_bf_census_m_5) != 0 & nrow(ddharm_bf_census_f_5) != 0){
  ddharm_smoothed_5 <- tibble()
  
  for(i in 3:ncol(ddharm_bf_census_f_5)){
    dat.f <- ddharm_bf_census_f_5 %>% select(1:2,i) %>% filter(!is.na(ddharm_bf_census_f_5[[i]])) %>% arrange(AgeStart)
    dat.m <- ddharm_bf_census_m_5 %>% select(1:2,i) %>% filter(!is.na(ddharm_bf_census_m_5[[i]])) %>% arrange(AgeStart)
    
    stopifnot(length(unique(diff(dat.f$AgeStart))) == 1 | length(unique(diff(dat.m$AgeStart))) == 1)
    
    # smoothed.adult <- getSmoothedPop5(popM = dat.m[[3]], popF = dat.f[[3]], Age = dat.m[[1]], EduYrs = 2, subgroup = "adult")
    # smoothed.child <- getSmoothedPop5(popM = dat.m[[3]], popF = dat.f[[3]], Age = dat.m[[1]], EduYrs = 2, subgroup = "child")
    # 
    # adult.grad <- as.numeric(str_extract(smoothed.adult$best_smooth_method,"\\s\\d"))
    # child.grad <- as.numeric(str_extract(smoothed.child$best_smooth_method,"\\s\\d"))
    # 
    # cat("adult bestGrad5 =",adult.grad,"; child bestGrad5 =",child.grad, "\n")
    # 
    # ddharm_smoothed_5 <- bind_rows(
    #   ddharm_smoothed_5,
    #   as_tibble(DemoTools::smooth_age_5(dat.m[[3]], dat.m[[1]], method="MAV", n=max(adult.grad, child.grad))) %>% mutate(age = dat.m[[1]], AgeLabel = dat.m[[2]], 
    #                                                                                                                      sex="male", year = as.numeric(names(ddharm_bf_census_f_5)[i])),
    #   as_tibble(DemoTools::smooth_age_5(dat.f[[3]], dat.f[[1]], method="MAV", n=max(adult.grad, child.grad))) %>% mutate(age = dat.f[[1]], AgeLabel = dat.f[[2]],
    #                                                                                                                      sex="female", year = as.numeric(names(ddharm_bf_census_m_5)[i]))
    # )
    
    smoothed <- census_workflow_adjust_smooth(popM = dat.m[[3]], popF = dat.f[[3]], Age = dat.m[[1]], EduYrs = 2)
    
    ddharm_smoothed_5 <- bind_rows(
      ddharm_smoothed_5,
      as_tibble(smoothed$popM_smooth) %>% mutate(age = 1:length(smoothed$popM_smoothed)-1, AgeLabel = 1:length(smoothed$popM_smoothed)-1,
                                                 sex="male", year = as.numeric(names(ddharm_bf_census_m_5)[i])),
      as_tibble(smoothed$popF_smooth) %>% mutate(age = 1:length(smoothed$popF_smoothed)-1, AgeLabel = 1:length(smoothed$popF_smoothed)-1,
                                                 sex="female", year = as.numeric(names(ddharm_bf_census_f_5)[i]))
    )
  }
  
  ddharm_smoothed_mat_5 <- ddharm_smoothed_5 %>% pivot_wider(names_from = year, values_from = value) %>% arrange(sex, age)
  
  ddharm_census_f_smoothed_5 <- ddharm_smoothed_mat_5 %>% filter(sex=="female") %>% select(-sex)
  ddharm_census_m_smoothed_5 <- ddharm_smoothed_mat_5 %>% filter(sex=="male") %>% select(-sex)
  
  ddharm_bf_census_f_oag_5 <- ddharm_census_f_smoothed_5 %>%
    group_by(age) %>%
    summarise_at(vars(-AgeLabel), function(y){first(na.omit(y))}) %>%
    mutate(aggr.age = ifelse(age >= open.age, open.age, age)) %>%
    group_by(aggr.age) %>%
    summarise_at(vars(-age), function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
    mutate(age = aggr.age) %>%
    group_by(age) %>%
    summarise_all(function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
    select(-2) 
  
  ddharm_bf_census_m_oag_5 <- ddharm_census_m_smoothed_5 %>%
    group_by(age) %>%
    summarise_at(vars(-AgeLabel), function(y){first(na.omit(y))}) %>%
    mutate(aggr.age = ifelse(age >= open.age, open.age, age)) %>%
    group_by(aggr.age) %>%
    summarise_at(vars(-age), function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
    mutate(age = aggr.age) %>%
    group_by(age) %>%
    summarise_all(function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
    select(-2) 
}

##WPP Pop Estimates
pop.singleyear <- wpp.pop.age.specific %>% filter(Location == country) %>% 
  # rename_with(.cols = X0:X100, ~str_extract(.x, "\\d+")) %>%
  # reshape2::melt(id.vars=c("Name", "Sex", "Reference")) %>% 
  select(AgeGrp, MidPeriod, PopMale, PopFemale) %>%
  pivot_longer(cols = 3:4) %>%
  mutate(year = MidPeriod,
         Sex = name,
         Sex = ifelse(Sex == 'PopMale', 'male', 'female'),
         value = value * 1000,
         age = AgeGrp) %>%
  select(Sex, year, age, value) %>%
  pivot_wider(values_from = value, names_from = year) %>%
  mutate(age = as.numeric(str_extract(age, "\\d+")))

pop.f.oag <- pop.singleyear  %>%
  filter(Sex=="female") %>%
  select(-Sex) %>%
  mutate(aggr.age = ifelse(age >= open.age, open.age, age)) %>%
  group_by(aggr.age) %>%
  summarise_at(vars(-age), function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
  mutate(age = aggr.age) %>%
  group_by(age) %>%
  summarise_all(function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
  select(-2) 

pop.m.oag <- pop.singleyear  %>%
  filter(Sex=="male") %>%
  select(-Sex) %>%
  mutate(aggr.age = ifelse(age >= open.age, open.age, age)) %>%
  group_by(aggr.age) %>%
  summarise_at(vars(-age), function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
  mutate(age = aggr.age) %>%
  group_by(age) %>%
  summarise_all(function(i){if(all(is.na(i))){NA} else {sum(i,na.rm=T)}}) %>%
  select(-2) 


bf.idx1<-projection_indices(period_start = 1960,
                            #period_end = max(floor(as.numeric(max(grep("\\d+", names(ddharm_bf_census_f), value=TRUE))) / 5)  * 5 + 5, 2015),
                            period_end = 2020,
                            interval = 1,
                            n_ages = n_ages,
                            fx_idx = 16L,
                            n_fx = 35L,
                            n_sexes = 2)

bf.idx5<-projection_indices(period_start = 1960,
                            #period_end = max(floor(as.numeric(max(grep("\\d+", names(ddharm_bf_census_f), value=TRUE))) / 5)  * 5 + 5, 2015),
                            period_end = 2020,
                            interval = 5,
                            n_ages = open.age %/% 5 + 1,
                            fx_idx = 4L,
                            n_fx = 4L,
                            n_sexes = 2)

#spline basis
knots.age <- seq(-3 * params$age.knot.space, open.age + 3 * params$age.knot.space, by = params$age.knot.space)
knots.time <- seq(-3 * params$year.knot.space + 1960, 2020 + 3 * params$year.knot.space, by = params$year.knot.space)
no.basis.age <- length(knots.age) - 4
no.basis.time <- length(knots.time) - 4

A.age<-c()
for(j in 1:no.basis.age) {
  A.age<-cbind(A.age, bspline(bf.idx1$ages, knots.age,j))
}

A.year<-c()
for(j in 1:no.basis.time) {
  A.year<-cbind(A.year, bspline(bf.idx1$periods, knots.time,j))
}

te.spline<-A.year %x% A.age

knots.fert <- seq(-3 * params$age.knot.space + min(bf.idx1$fertility_ages), 3 * params$age.knot.space + max(bf.idx1$fertility_ages), by = params$age.knot.space)
no.basis.fert <- length(knots.fert) - 4

A.age.fert<-c()
for(j in 1:no.basis.fert) {
  A.age.fert<-cbind(A.age.fert,bspline(bf.idx1$fertility_ages, knots.fert, j))
}

te.spline.fert <- A.year %x% A.age.fert 

##WPP fx
fx <- wpp.fx %>% filter(Name == country) %>% reshape2::melt() %>% 
  mutate(year = as.numeric(str_extract(Reference,"\\d{4}")), 
         fx = value/1000,
         age = as.numeric(str_extract(variable,"\\d{2}"))) %>%
  select(year,age, fx)

fx_init.singleyear <-  fx %>% filter(year %in% bf.idx1$periods) %>%
  pivot_wider(names_from=year, values_from=fx) %>% 
  select(-age) %>% as.matrix() %>%
  apply(1, function(i){approx(x=seq(1960, 2015, by = 5)+2, y=i, xout=1960:2020, rule=2)$y}) %>%
  apply(1, function(i){approx(x=seq(15,45,by=5)+2 , y=i, xout=15:49, rule=2)$y}) %>%
  `colnames<-`(1960:2020)

log_fx_mean <- as.vector(log(fx_init.singleyear))

##DHS data cohort smoothed
bf5.smooth <- aggr.mat.cohort.0[[country]] %>%
  filter(period %in% bf.idx1$periods) %>%
  group_by(mm1, tips, DHS, agegr, period) %>%
  summarise_at(vars(pyears, event, pyears2, adjusted), sum) %>%
  ungroup()

bf5.smooth$period <- factor(bf5.smooth$period,levels=bf.idx1$periods)
#bf5.smooth$tips <- factor(bf5.smooth$tips,levels=0:14)

dhs.start.age <- 15
dhs.end.age <- 59

bf5.f.no0.smooth <- bf5.smooth %>% filter(mm1=="female", agegr >= dhs.start.age, agegr <= dhs.end.age) %>% arrange(period,tips,agegr)
bf5.m.no0.smooth <- bf5.smooth %>% filter(mm1=="male", agegr >= dhs.start.age, agegr <= dhs.end.age) %>% arrange(period,tips,agegr)

#BR Data
if(params$br.data){
  bf5.smooth.br <- aggr.mat.br[[country]] %>%
    filter(period %in% bf.idx1$periods) %>%
    group_by(b4, tips, DHS, agegr, period) %>%
    summarise_at(vars(pyears, event), sum) %>%
    ungroup() %>%
    mutate(tips = as.numeric(levels(tips))[tips])
} else {
  bf5.smooth.br <- aggr.mat.br[[country]] %>% slice(0)
}

bf5.smooth.br$period <- factor(bf5.smooth.br$period,levels=bf.idx1$periods)
#bf5.smooth.br$tips <- factor(bf5.smooth.br$tips,levels=0:14)

br.start.age <- 0
br.end.age <- 20

bf5.f.no0.br <- bf5.smooth.br %>% filter(b4 == "female", agegr >= br.start.age, agegr <= br.end.age) %>% arrange(period, tips, agegr)
bf5.m.no0.br <- bf5.smooth.br %>% filter(b4 == "male", agegr >= br.start.age, agegr <= br.end.age) %>% arrange(period, tips, agegr)

# if(is.null(aggr.mat.br[[joint.countries[[i]]]]) || nrow(aggr.mat.br[[joint.countries[i]]]) == 0){
#   bf5.smooth.br <- data.mat.br.m  <- data.mat.br.f <- aggr.mat.br$Angola %>% slice(0)
# } else {
#   bf5.smooth.br <- aggr.mat.br[[joint.countries[i]]] %>%
#     #  filter(period < 2012) %>%
#     filter(period %in% year.start:year.end) %>%
#     group_by(b4, tips, DHS, agegr, period) %>%
#     summarise_at(vars(pyears, event), sum) %>%
#     ungroup() %>%
#     mutate(tips = as.numeric(levels(tips))[tips])
#   
#   data.mat.br.m <- bf5.smooth.br %>% filter(b4 == "male", agegr >= br.start.age, agegr <= br.end.age) %>% arrange(period, tips, agegr)
#   data.mat.br.f <- bf5.smooth.br %>% filter(b4 == "female", agegr >= br.start.age, agegr <= br.end.age) %>% arrange(period, tips, agegr)
# }

igme.h.mean.f <- igme.5q0.df %>% filter(Sex=="Female", year %in% bf.idx1$periods,  !year %in% unique(bf5.smooth.br$period)) %>% right_join(as_tibble(bf.idx1$periods)%>%mutate(year=value)) %>% .$child.mort %>% log()
h.mean.f <- wpp.5q0.interpolate %>% filter(Sex=="Female", year %in% bf.idx1$periods, !year %in% unique(bf5.smooth.br$period)) %>% .$q50 %>% log() ##Using WPP 5q0 estimates

igme.h.mean.m <- igme.5q0.df %>% filter(Sex=="Male", year %in% bf.idx1$periods, !year %in% unique(bf5.smooth.br$period)) %>% right_join(as_tibble(bf.idx1$periods)%>%mutate(year=value)) %>% .$child.mort %>% log()
h.mean.m <- wpp.5q0.interpolate %>% filter(Sex=="Male", year %in% bf.idx1$periods, !year %in% unique(bf5.smooth.br$period)) %>% .$q50 %>% log() ##Using WPP 5q0 estimates

#smooth monotone
ddharm_bf_census_f_oag_mono <- ddharm_bf_census_m_oag_mono <-
  ddharm_bf_census_f_oag_5_mono <- ddharm_bf_census_m_oag_5_mono <-
  ddharm_bf_census_f_oag[1]

ddharm_bf_census_f_oag_mono <- ddharm_bf_census_f_oag %>% 
  mutate(age.5 = age %/% 5) %>%
  group_by(age.5) %>%
  summarise_at(vars(-age), sum) %>%
  .[,-1] %>%
  as.matrix %>%
  apply(2, function(i){c(DemoTools::agesmth1(Value = i[!is.na(i)], Age = 0:(length(i[!is.na(i)]) - 1) * 5, method = "loess", OAG= TRUE), i[is.na(i)])}) %>%
  apply(2, function(i){c(DemoTools::graduate(Value = i[!is.na(i)], Age = 0:(length(i[!is.na(i)]) - 1) * 5, method = "mono", OAG= TRUE), rep(NA, length(which(is.na(i))) * 5))}) %>%
  as_tibble %>%
  mutate(age = 0:open.age, .before = 1)

ddharm_bf_census_m_oag_mono <- ddharm_bf_census_m_oag %>% 
  mutate(age.5 = age %/% 5) %>%
  group_by(age.5) %>%
  summarise_at(vars(-age), sum) %>%
  .[,-1] %>%
  as.matrix %>%
  apply(2, function(i){c(DemoTools::agesmth1(Value = i[!is.na(i)], Age = 0:(length(i[!is.na(i)]) - 1) * 5, method = "loess", OAG= TRUE), i[is.na(i)])}) %>%
  apply(2, function(i){c(DemoTools::graduate(Value = i[!is.na(i)], Age = 0:(length(i[!is.na(i)]) - 1) * 5, method = "mono", OAG= TRUE), rep(NA, length(which(is.na(i))) * 5))}) %>%
  as_tibble %>%
  mutate(age = 0:open.age, .before = 1)

ddharm_bf_census_f_oag_5_mono <- ddharm_bf_census_f_oag_5 %>% 
  mutate(age.5 = age %/% 5) %>%
  group_by(age.5) %>%
  summarise_at(vars(-age), sum) %>%
  .[,-1] %>%
  as.matrix %>%
  apply(2, function(i){c(DemoTools::agesmth1(Value = i[!is.na(i)], Age = 0:(length(i[!is.na(i)]) - 1) * 5, method = "loess", OAG= TRUE), i[is.na(i)])}) %>%
  apply(2, function(i){c(DemoTools::graduate(Value = i[!is.na(i)], Age = 0:(length(i[!is.na(i)]) - 1) * 5, method = "mono", OAG= TRUE), rep(NA, length(which(is.na(i))) * 5))}) %>%
  as_tibble %>%
  mutate(age = 0:open.age, .before = 1)

ddharm_bf_census_m_oag_5_mono <- ddharm_bf_census_m_oag_5 %>% 
  mutate(age.5 = age %/% 5) %>%
  group_by(age.5) %>%
  summarise_at(vars(-age), sum) %>%
  .[,-1] %>%
  as.matrix %>%
  apply(2, function(i){c(DemoTools::agesmth1(Value = i[!is.na(i)], Age = 0:(length(i[!is.na(i)]) - 1) * 5, method = "loess", OAG= TRUE), i[is.na(i)])}) %>%
  apply(2, function(i){c(DemoTools::graduate(Value = i[!is.na(i)], Age = 0:(length(i[!is.na(i)]) - 1) * 5, method = "mono", OAG= TRUE), rep(NA, length(which(is.na(i))) * 5))}) %>%
  as_tibble %>%
  mutate(age = 0:open.age, .before = 1)

if(params$supersmooth.pop){
  ddharm_bf_census_f_oag <- ddharm_bf_census_f_oag_mono; ddharm_bf_census_m_oag <- ddharm_bf_census_m_oag_mono
  ddharm_bf_census_f_oag_5 <- ddharm_bf_census_f_oag_5_mono; ddharm_bf_census_m_oag_5 <- ddharm_bf_census_m_oag_5_mono
}

data.f <- if(nrow(ddharm_bf_census_f)!=0) {as.matrix(log(ddharm_bf_census_f_oag[,-1]))} else {as.matrix(tibble(.rows=open.age+1))}; data.m <- if(nrow(ddharm_bf_census_m)!=0) {as.matrix(log(ddharm_bf_census_m_oag[,-1]))} else {as.matrix(tibble(.rows=open.age+1))}
data.f.5 <- if(nrow(ddharm_bf_census_f)!=0) {as.matrix(log(ddharm_bf_census_f_oag_5[,-1] %>% select(!matches(colnames(data.f)))))} else if(exists("ddharm_bf_census_f_oag_5") && nrow(ddharm_bf_census_f_oag_5) != 0) {
  as.matrix(log(ddharm_bf_census_f_oag_5[,-1]))} else {as.matrix(tibble(.rows=open.age+1))}
data.m.5 <- if(nrow(ddharm_bf_census_m)!=0){as.matrix(log(ddharm_bf_census_m_oag_5[,-1] %>% select(!matches(colnames(data.f)))))} else if(exists("ddharm_bf_census_f_oag_5") && nrow(ddharm_bf_census_f_oag_5) != 0) {
  as.matrix(log(ddharm_bf_census_m_oag_5[,-1]))} else {as.matrix(tibble(.rows=open.age+1))}

# if(country == "Zimbabwe"){
#   data.f <- as.matrix(log(ddharm_bf_census_f_oag[,-(1:2)])); data.m <- as.matrix(log(ddharm_bf_census_m_oag[,-(1:2)]))
#   data.f.5 <- as.matrix(log(ddharm_bf_census_f_oag_5[,-(1:2)] %>% select(!matches(colnames(data.f))))); data.m.5 <- as.matrix(log(ddharm_bf_census_m_oag_5[,-(1:2)] %>% select(!matches(colnames(data.f)))))
#   #data.f.55 <- as.matrix(log(ddharm_bf_census_f_oag_5[,-(1:2)])); data.m <- as.matrix(log(ddharm_bf_census_m_oag_5[,-(1:2)]))
# }

# if(!length(unique(bf5.smooth$DHS)) == 1 || min(diff(as.numeric(str_extract(unique(bf5.smooth$DHS), "\\d+")))) < 14) {
#   skip<-c("Ethiopia","Central African Republic","Comoros","Sao Tome and Principe","Botswana","Cape Verde","Equatorial Guinea","Eritrea","Nigeria (Ondo State)","Ghana","Mauritania","Sudan", "Guinea-Bissau")
#   joint.countries<-names(aggr.mat.cohort.0)[!names(aggr.mat.cohort.0)%in%skip]
# 
#   par.all <- split(unname(tmb.full.joint.common$env$last.par.best), names(tmb.full.joint.common$env$last.par.best))
#   
#   tp.common <- par.all$tips_params_common
#   tp.common[6]<-tp.common[6] + par.all$tp_common_5
#   tp.common[11]<-tp.common[11] + par.all$tp_common_10
#   
#   tp <- split(par.all$tips_params, rep(joint.countries, each = length(par.all$tips_params_common)))[params$country]
# 
#   tp.init <- tp[[1]] + tp.common; map <- list(log_lambda_tp = factor(NA),
#                                               tp_params = factor(rep(NA,15)),
#                                               tp_slope = factor(NA),
#                                               tp_params_5 = factor(NA),
#                                               tp_params_10 = factor(NA))
# } else {tp.init <- rep(0,15); map <- list()}

if(!length(unique(bf5.smooth$DHS)) == 1 || min(diff(as.numeric(str_extract(unique(bf5.smooth$DHS), "\\d+")))) < 14) {
  fit.life.table <- readRDS("C:/Users/ktang3/Desktop/Imperial/Life_Table_System/JE/fit_2022-07-08.rds")
  
  br.skip1 <- names(aggr.mat.br)[which(lapply(aggr.mat.br, length) == 0)]
  br.joint.countries1 <- names(aggr.mat.br)[!names(aggr.mat.br) %in% br.skip1]
  br.joint.countries2 <- names(aggr.mat.cohort.0)
  br.joint.countries <- intersect(br.joint.countries1, br.joint.countries2) %>% setdiff(c("Rwanda", "Mali"))
  
  tp.params <- fit.life.table$mode$tp_params[, match(params$country, br.joint.countries)] + fit.life.table$mode$tp_params_common
  tp.params.age0 <- fit.life.table$mode$tp_params_br_age0[, match(params$country, br.joint.countries)] + fit.life.table$mode$tp_params_br_common_age0
  tp.params.age1 <- fit.life.table$mode$tp_params_br_age1[, match(params$country, br.joint.countries)] + fit.life.table$mode$tp_params_br_common_age1
  tp.params.age5 <- fit.life.table$mode$tp_params_br_age5[, match(params$country, br.joint.countries)] + fit.life.table$mode$tp_params_br_common_age5
  tp.params.age10 <- fit.life.table$mode$tp_params_br_age10[, match(params$country, br.joint.countries)] + fit.life.table$mode$tp_params_br_common_age10
  
  tp.list <- append(
    list(tp_params = tp.params,
         tp_params_br_age0 = tp.params.age0,
         tp_params_br_age1 = tp.params.age1,
         tp_params_br_age5 = tp.params.age5,
         tp_params_br_age10 = tp.params.age10),
    
    split(unname(fit.life.table$par.full), names(fit.life.table$par.full)) %>%
      .[names(.) %>% str_detect("tp_slope|_5|_10")] %>%
      `names<-`(names(.) %>% str_replace("_common", ""))) %>%
    append( split(unname(fit.life.table$par.full), names(fit.life.table$par.full)) %>%
              .[names(.) %>% str_detect("lambda_tp")] %>%
              .[!names(.) %>% str_detect("common")])
  
} else {
  tp.list <- list(
    log_lambda_tp = 5,
    tp_params = rep(0, 15),
    tp_slope = 0,
    tp_params_5 = 0,
    tp_params_10 = 0,
    
    log_lambda_tp_br_age0 = 5,
    tp_params_br_age0 = rep(0, 15),
    tp_slope_br_age0 = 0,
    tp_params_5_br_age0 = 0,
    tp_params_10_br_age0 = 0,
    
    log_lambda_tp_br_age1 = 5,
    tp_params_br_age1 = rep(0, 15),
    tp_slope_br_age1 = 0,
    tp_params_5_br_age1 = 0,
    tp_params_10_br_age1 = 0,
    
    log_lambda_tp_br_age5 = 5,
    tp_params_br_age5 = rep(0, 15),
    tp_slope_br_age5 = 0,
    tp_params_5_br_age5 = 0,
    tp_params_10_br_age5 = 0,
    
    log_lambda_tp_br_age10 = 5,
    tp_params_br_age10 = rep(0, 15),
    tp_slope_br_age10 = 0,
    tp_params_5_br_age10 = 0,
    tp_params_10_br_age10 = 0,
  )}

basepop.f <- ifelse(pop.f.oag$`1960`==0, 1, pop.f.oag$'1960')
basepop.m <- ifelse(pop.m.oag$`1960`==0, 1, pop.m.oag$'1960')

# tp.init <- rep(0, 15)

#cor matrix from HMD/LTs
load('../../Life_Table_System/smoothed pca log new bpca cov inflated err young ages pure prior partial cov.RData')

if(params$loghump){
  par.df <- readRDS("../../Life_Table_System/fit single LT HP PCA new bpca cov inflated err young ages pure prior partial cov loghump.rds")
  rep.par <- readRDS("../../Life_Table_System/fit LT par cov PCA new bpca cov inflated err young ages pure prior partial loghump constrained.rds")
} else {
  par.df <- readRDS("../../Life_Table_System/fit single LT HP PCA new bpca cov inflated err young ages pure prior partial cov.rds")
  rep.par <- readRDS("../../Life_Table_System/fit LT par cov PCA new bpca cov inflated err young ages pure prior partial cov.rds")
}

V_mf_mean <- rep.par$V_mf_mean
V_mf_slope <- rep.par$V_mf_slope
V_mf <- rep.par$V_mf

log.spline.prior.m <- par.df %>% filter(sex == "Male") %>% select('phi':'B') %>% apply(2, mean)
log.spline.prior.f <- par.df %>% filter(sex == "Female") %>% select('phi':'B') %>% apply(2, mean)
log.spline.prior.m[c(1:3, 7, 8)] <- rep.par$log_spline_means_m_common
log.spline.prior.f[c(1:3, 7, 8)] <- rep.par$log_spline_means_f_common

init_lambda_f <- exp(log.spline.prior.f[4]); init_lambda_m <- exp(log.spline.prior.m[4]); 
init_delta_f <- exp(log.spline.prior.f[5]); init_delta_m <- exp(log.spline.prior.m[5]); 
init_epsilon_f <- exp(log.spline.prior.f[6]); init_epsilon_m <- exp(log.spline.prior.m[6]); 

full.penal.gx <- as(0.5 * diag(no.basis.time) %x% crossprod(diff(diag(no.basis.age),differences=1)) +
                      0.5 * crossprod(diff(diag(no.basis.time),differences=1)) %x% diag(no.basis.age) +
                      #0.5 * crossprod(diff(diag(no.basis.time))%x%diff(diag(no.basis.age))) +
                      1e-3 * diag(no.basis.time * no.basis.age), "sparseMatrix")

AR1.PREC <- function(n, rho) {
  P <- diag(n)
  for(i in 2:n){
    P[i, i-1] <- -rho
  }
  P <- P / sqrt(1-rho^2)
  P[1,1] <- 1
  
  P
}

AR2.PREC <- function(n, phi){
  #sigma2 <- 1
  #gamma0 <- (1-phi[2]) / (1+phi[2]) * sigma2 / ((1-phi[2])^2 - phi[1]^2)
  
  #gamma0 = 1, all variance = 1
  sigma2 <- (1+phi[2]) / (1-phi[2]) * ((1-phi[2])^2 - phi[1]^2)
  gamma.init <- matrix(c(1, phi[1] / (1-phi[2]), phi[1]/(1-phi[2]),1),2,2)
  
  VAR <- adiag(gamma.init, diag(rep(sigma2, n-2)))
  T <- diag(n) - phi[1] * rbind(matrix(0, 2, n), cbind(rep(0, n-2), diag(n-2), rep(0, n-2))) - phi[2] * rbind(matrix(0, 2, n), cbind(diag(n-2), matrix(0, n-2, 2)))
  
  t(T) %*% solve(VAR) %*% T
}

full.penal.gx.AR2.tensor <- as(AR2.PREC(no.basis.time, c(2*0.5, -0.5)) %x% AR2.PREC(no.basis.age, c(2*0.5, -0.5)), "sparseMatrix")
full.penal.fx.AR2.tensor <- as(AR2.PREC(no.basis.time, c(2*0.5, -0.5)) %x% AR2.PREC(no.basis.fert, c(2*0.5, -0.5)), "sparseMatrix")

gumbel.theta.AR2.marginal.fx <- -log(0.01) *  sqrt(mean(diag(te.spline.fert %*% solve(full.penal.fx.AR2.tensor) %*% t(te.spline.fert)))) * 1.96 / log(1.1)
gumbel.theta.AR2.marginal.gx <- -log(0.01) *  sqrt(mean(diag(te.spline %*% solve(full.penal.gx.AR2.tensor) %*% t(te.spline)))) * 1.96 / 0.08

ARIMA.invD <- solve(rbind(c(1, rep(0, no.basis.time-1)), diff(diag(no.basis.time))))
ARIMA.vec <- diff(diag(2)) %*% A.year[1:2,2:4] %*% ARIMA.invD[2:4,2:4] #omitting initial variance as coefficient = 0

gumbel.theta.ARIMA.marginal.lambda <- 
  gumbel.theta.ARIMA.marginal.delta <- 
  gumbel.theta.ARIMA.marginal.epsilon <- -log(0.01) *  
  sqrt(ARIMA.vec %*% solve(AR2.PREC(no.basis.time-1, c(0.5, 0)))[1:3, 1:3] %*% t(ARIMA.vec)) * 
  1.96 / log(1.05)

# RW2.invD <- round(solve(rbind(cbind(diag(2), matrix(0, 2, no.basis.time - 2)), diff(diag(no.basis.time), differences = 2))), digits = 10)
# RW2.vec <- round(diff(diag(bf.idx1$n_periods), differences = 2) %*% A.year %*% RW2.invD, digits = 10)
# 
# gumbel.theta.RW2 <- -log(0.01) *  
#   sqrt(mean(diag(tcrossprod(RW2.vec)))) / RW2.sd  ## = log-ratio of yearly (percentage change + 1), i.e. log( (y_2 / y_1)  / (y_1 / y_0))
# #1.96 / log(1.005) ## = log-ratio of yearly (percentage change + 1), i.e. log( (y_2 / y_1)  / (y_1 / y_0))
# 
# gumbel.theta.RW2.standardised <- -log(0.01) *  
#   sqrt(mean(diag(tcrossprod(RW2.vec)))) / RW2.sd.standardised  ## = log-ratio of yearly (percentage change + 1), i.e. log( (y_2 / y_1)  / (y_1 / y_0))
# 
# gumbel.theta.RW2.manual <- -log(0.01) *  
#   sqrt(mean(diag(tcrossprod(RW2.vec)))) *
#   1.96 / log(1.1) 

d.rho <- function(rho, n) {
  sqrt( (1-n)*log(1+3*rho) + (3-2*n)*log(1-rho) + n*log(1+rho) )
}

upper.gamma <- 40
hiv.cut <- 12

get.ar1.rho.alpha.func <- function(theta, L, alpha) {
  ( (exp(-theta*sqrt(1-L)) - exp(-theta)) / (1 - exp(-theta)) - alpha)^2
}
get.ar2.rho.alpha.func <- function(theta, L, alpha) {
  ( (exp(-theta * sqrt((1-L)^2 / (1+L))) - exp(-theta)) / (1 - exp(-theta) ) - alpha)^2
}

D2 <- diff(diag(no.basis.time),differences=2)
inv.spline.err.mat <- as(t(D2) %*% solve(tcrossprod(D2)), "sparseMatrix")
spline.trend <- (1:no.basis.time - 1) * params$year.knot.space

if(params$census.deaths){
  census.deaths.f <- if(is.null(ddharm_bf_census_deaths_f) || nrow(ddharm_bf_census_deaths_f) == 0) {as.matrix(tibble(.rows=open.age+1))} else {as.matrix(ddharm_bf_census_deaths_f[,-1])}
  census.deaths.m <- if(is.null(ddharm_bf_census_deaths_m) || nrow(ddharm_bf_census_deaths_m) == 0) {as.matrix(tibble(.rows=open.age+1))} else {as.matrix(ddharm_bf_census_deaths_m[,-1])} 
  census.deaths.f.5 <- if(is.null(ddharm_bf_census_deaths_f_5) || nrow(ddharm_bf_census_deaths_f_5) == 0) {as.matrix(tibble(.rows=open.age / 5 + 1))} else if(!is.null(colnames(census.deaths.f))) {as.matrix(full_join(tibble(age = seq(0, open.age, by = 5)), ddharm_bf_census_deaths_f_5)[,-1] %>% select(!matches(colnames(census.deaths.f))))} else {as.matrix(full_join(tibble(age = seq(0, open.age, by = 5)), ddharm_bf_census_deaths_f_5)[,-1])}
  census.deaths.m.5 <- if(is.null(ddharm_bf_census_deaths_f_5) || nrow(ddharm_bf_census_deaths_m_5) == 0) {as.matrix(tibble(.rows=open.age / 5 + 1))} else if(!is.null(colnames(census.deaths.m))) {as.matrix(full_join(tibble(age = seq(0, open.age, by = 5)), ddharm_bf_census_deaths_m_5)[,-1] %>% select(!matches(colnames(census.deaths.f))))} else {as.matrix(full_join(tibble(age = seq(0, open.age, by = 5)), ddharm_bf_census_deaths_m_5)[,-1])}
} else {
  census.deaths.f <- as.matrix(tibble(.rows=open.age+1))
  census.deaths.m <- as.matrix(tibble(.rows=open.age+1))
  census.deaths.f.5 <- as.matrix(tibble(.rows=open.age / 5 + 1))
  census.deaths.m.5 <- as.matrix(tibble(.rows=open.age / 5 + 1))
}


A.year.1900 <- c()
for(j in 1:51) {
  A.year.1900 <- cbind(A.year.1900, bspline(1900:2020, seq(-3 * 2.5 + 1900, 2020 + 3 * 2.5, by = 2.5),j))
}

C <- rbind(1, 1:nrow(A.year.1900)) %*% A.year.1900
Z <- qr.Q(qr(t(C)), complete = T)[, -(1:2)]
P <- diff(diag(ncol(A.year.1900)), differences = 2)
inv.VAR <- solve(P %*% Z) 
A.err.const <- A.year.1900 %*% Z %*% inv.VAR

theta_rho <- optim(3, get.ar2.rho.alpha.func, method="L-BFGS-B", L = 0.1, alpha = 0.01)$par
thiele_age <- c(0.01, 1:(nrow(smoothed.pca) - 1))

#create dat and par list####
# data.f <- data.f.5 <- data.m <- data.m.5 <- as.matrix(tibble(.rows=open.age+1))
# census.deaths.m.5 <- census.deaths.f.5 <-  as.matrix(tibble(.rows=open.age %/% 5+1))
# bf5.f.no0.smooth <- bf5.f.no0.smooth[1,]; bf5.m.no0.smooth <- bf5.m.no0.smooth[1,]
# skip<-c("Ethiopia","Central African Republic","Comoros","Sao Tome and Principe","Botswana","Cape Verde","Equatorial Guinea","Eritrea","Nigeria (Ondo State)","Ghana","Mauritania","Sudan")
# joint.countries<-names(aggr.mat.cohort.0)[!names(aggr.mat.cohort.0)%in%skip]
# 
# all.list<-function(x,age.start,age.end,year.start,year.end){
#   no.basis = 15
#   knots<-seq(0,1,length=no.basis-2)
#   dk<-knots[2]-knots[1]
#   knots<-c(knots[1]-dk*(3:1),knots,knots[no.basis-2]+dk*(1:3))
#   age<-seq(0,1,length=age.end-age.start+1)
#   year<-seq(0,1,length=year.end-year.start+1)
# 
#   A.age<-c()
#   for(j in 1:no.basis) {
#     A.age<-cbind(A.age,bspline(age,knots,j))
#   }
# 
#   A.year<-c()
#   for(j in 1:no.basis) {
#     A.year<-cbind(A.year,bspline(year,knots,j))
#   }
# 
#   te.spline<-A.year%x%A.age
# 
#   ind<-function(x){
#     aggr.mat.reduced<-x
#     data.mat.m<-aggr.mat.reduced[aggr.mat.reduced$mm1=="male" & aggr.mat.reduced$agegr %in% age.start:age.end & aggr.mat.reduced$period %in% year.start:year.end,c(7,5,3,2,4)]
#     data.mat.f<-aggr.mat.reduced[aggr.mat.reduced$mm1=="female" & aggr.mat.reduced$agegr %in% age.start:age.end & aggr.mat.reduced$period %in% year.start:year.end,c(7,5,3,2,4)]
# 
#     age.start.col<-min(which(A.age[min(data.mat.m[,3])-age.start+1,]!=0))
#     age.end.col<-max(which(A.age[max(data.mat.m[,3])-age.start+1,]!=0))
#     year.start.col<-min(which(A.year[min(data.mat.m[,4])-year.start+1,]!=0))
#     year.end.col<-max(which(A.year[max(data.mat.m[,4])-year.start+1,]!=0))
#     age.start<-min(data.mat.m[,3])
#     age.end<-max(data.mat.m[,3])
#     year.start<-min(data.mat.m[,4])
#     year.end<-max(data.mat.m[,4])
#     tp<-max(data.mat.m[,5],data.mat.f[,5])
# 
#     ind.place<-c(age.start.col,age.end.col,year.start.col,year.end.col,age.start,age.end,year.start,year.end,tp)
#     names(ind.place)<-c("age.start.col", "age.end.col","year.start.col","year.end.col","age.start","age.end","year.start","year.end","tp")
#     ind.place
#   }
# 
# 
#   everything.func<-function(i){
#     age.seq<-as.numeric(ind.list[[i]][["age.start"]]:ind.list[[i]][["age.end"]])
#     year.seq<-as.numeric(ind.list[[i]][["year.start"]]:ind.list[[i]][["year.end"]])
# 
#     age.m<-A.age[age.seq-age.start+1,]%*%joint.countries.age.m[[i]]
#     age.f<-A.age[age.seq-age.start+1,]%*%joint.countries.age.f[[i]]
#     time.m<-A.year[year.seq-year.start+1,]%*%joint.countries.time.m[[i]]
#     time.f<-A.year[year.seq-year.start+1,]%*%joint.countries.time.f[[i]]
#     agetime.m<-matrix(te.spline[age.seq-age.start+1+rep((age.end-age.start+1)*(year.seq-year.start),each=length(age.seq)),]%*%joint.countries.2d.m[[i]],length(age.seq),length(year.seq))
#     agetime.f<-matrix(te.spline[age.seq-age.start+1+rep((age.end-age.start+1)*(year.seq-year.start),each=length(age.seq)),]%*%joint.countries.2d.f[[i]],length(age.seq),length(year.seq))
# 
#     mort.m<-joint.countries.avg.m[[i]]+t(apply(apply(agetime.m,2,function(x){x+age.m}),1,function(x){x+time.m}))
#     mort.f<-joint.countries.avg.f[[i]]+t(apply(apply(agetime.f,2,function(x){x+age.f}),1,function(x){x+time.f}))
#     rownames(age.m)<-rownames(age.f)<-rownames(agetime.m)<-rownames(agetime.f)<-rownames(mort.m)<-rownames(mort.f)<-age.seq
#     rownames(time.m)<-rownames(time.f)<-colnames(agetime.m)<-colnames(agetime.f)<-colnames(mort.m)<-colnames(mort.f)<-year.seq
# 
# 
#     if(i=="Rwanda"){
#       mort.m[,4:6]<-mort.m[,4:6]+matrix(rep.par.list$rwanda_geno_m,56,3,byrow=T)
#       mort.f[,4:6]<-mort.f[,4:6]+matrix(rep.par.list$rwanda_geno_f,56,3,byrow=T)
#     }
# 
#     tips<-joint.countries.tp.c[[i]]+tp.common
# 
#     list(tp=tips,avg.m=joint.countries.avg.m[[i]],avg.f=joint.countries.avg.f[[i]],age.m=age.m,age.f=age.f,time.m=time.m,time.f=time.f,agetime.m=agetime.m,agetime.f=agetime.f,mort.m=mort.m,mort.f=mort.f)
#   }
# 
#   ind.list<-lapply(aggr.mat.cohort.0[joint.countries],ind)
#   rep.par.list<-split(unname(x$env$last.par.best),names(x$env$last.par.best))
#   joint.countries.avg.m<-split(rep.par.list$avg_m,joint.countries)
#   joint.countries.avg.f<-split(rep.par.list$avg_f,joint.countries)
#   joint.countries.age.m<-split(rep.par.list$spline_params_m_age,rep(joint.countries,each=no.basis))
#   joint.countries.age.f<-split(rep.par.list$spline_params_f_age,rep(joint.countries,each=no.basis))
#   joint.countries.time.m<-split(rep.par.list$spline_params_m_time,rep(joint.countries,each=no.basis))
#   joint.countries.time.f<-split(rep.par.list$spline_params_f_time,rep(joint.countries,each=no.basis))
#   joint.countries.2d.m<-split(rep.par.list$spline_params_m_2d,rep(joint.countries,each=no.basis*no.basis))
#   joint.countries.2d.f<-split(rep.par.list$spline_params_f_2d,rep(joint.countries,each=no.basis*no.basis))
# 
#   joint.countries.tp.c<-split(rep.par.list$tips_params,rep(joint.countries,each=length(rep.par.list$tips_params_common)))
#   tp.common<-rep.par.list$tips_params_common
#   tp.common[6]<-tp.common[6]+rep.par.list$tp_common_5
#   tp.common[11]<-tp.common[11]+rep.par.list$tp_common_10
# 
# 
#   setNames(lapply(joint.countries,everything.func),joint.countries)
# }
# 
# more.countries.avg<-all.list(tmb.full.joint.common,age.start=10,age.end=65,year.start=1990,year.end=2017)
# 
# for(i in 1:length(more.countries.avg)){
#   more.countries.avg[[i]]$avg.m<-more.countries.avg[[i]]$avg.m+unname(tmb.full.joint.common$env$last.par.best["avg_common_m"])
#   more.countries.avg[[i]]$avg.f<-more.countries.avg[[i]]$avg.f+unname(tmb.full.joint.common$env$last.par.best["avg_common_f"])
#   more.countries.avg[[i]]$mort.m<-more.countries.avg[[i]]$mort.m+unname(tmb.full.joint.common$env$last.par.best["avg_common_m"])
#   more.countries.avg[[i]]$mort.f<-more.countries.avg[[i]]$mort.f+unname(tmb.full.joint.common$env$last.par.best["avg_common_f"])
# }
# 
# tp.init <- more.countries.avg[[params$country]]$tp; map <- list(log_lambda_tp = factor(NA),
#                                                                 tp_params = factor(rep(NA,15)),
#                                                                 tp_slope = factor(NA),
#                                                                 tp_params_5 = factor(NA),
#                                                                 tp_params_10 = factor(NA))

# 
# data.loghump.vec.RW <- list(log_basepop_mean_f = log(basepop.f), log_basepop_mean_m = log(basepop.m),
#                             log_fx_mean = log_fx_mean,
#                             srb = rep(1.05, bf.idx1$n_periods),
#                             interval = bf.idx1$interval,
#                             n_periods = bf.idx1$n_periods,
#                             fx_idx = bf.idx1$fx_idx,
#                             n_fx = bf.idx1$n_fx,
#                             census_log_pop_f = data.f, census_log_pop_m = data.m,
#                             census_year_idx = match(bf.idx1$interval * floor(as.numeric(colnames(data.f)) / bf.idx1$interval), bf.idx1$periods_out),
#                             #census_year_grow_idx = as.numeric(colnames(data.f)) - bf.idx1$interval * floor(as.numeric(colnames(data.f)) / bf.idx1$interval),
#                             
#                             open_idx = bf.idx1$n_ages,
#                             #oag = apply(data.f, 2, function(i){length(na.omit(i))}),
#                             pop_start = rep(6, ncol(data.f)), 
#                             pop_end = ifelse(apply(data.f, 2, function(i){length(na.omit(i))})-1 > 75 + 1, 75 + 1, apply(data.f, 2, function(i){length(na.omit(i))})-1),
#                             #pop_end = apply(data.f, 2, function(i){length(na.omit(i))})-1,
#                             
#                             census_log_pop_f_5 = data.f.5, census_log_pop_m_5 = data.m.5,
#                             census_year_idx_5 = match(bf.idx1$interval * floor(as.numeric(colnames(data.f.5)) / bf.idx1$interval), bf.idx1$periods_out),
#                             #census_year_group_idx_5 = bf.idx1$ages %/% 5 + 1,
#                             #n_agegrp_5 = bf.idx5$n_ages,
#                             #pop_start_5 = rep(2, ncol(data.f.5)),
#                             #pop_end_5 =  ifelse(apply(data.f.5, 2, function(i){length(na.omit(i))})-1 > 15, 15, apply(data.f.5, 2, function(i){length(na.omit(i))})-1),
#                             
#                             census_year_group_idx_5 = bf.idx1$ages + 1,
#                             n_agegrp_5 = bf.idx1$n_ages,
#                             pop_start_5 = rep(6, ncol(data.f.5)),
#                             pop_end_5 =  ifelse(apply(data.f.5, 2, function(i){length(na.omit(i))})-1 > 75 + 1, 75 + 1, apply(data.f.5, 2, function(i){length(na.omit(i))})-1),
#                             
#                             census_deaths_f = census.deaths.f, census_deaths_m = census.deaths.m, 
#                             census_deaths_year_idx = match(bf.idx1$interval * floor(as.numeric(colnames(census.deaths.f)) / bf.idx1$interval), bf.idx1$periods_out),
#                             deaths_start = rep(6, ncol(census.deaths.f)),
#                             deaths_end = ifelse(apply(census.deaths.f, 2, function(i){length(na.omit(i))})-1 > 75 + 1, 75 + 1, apply(census.deaths.f, 2, function(i){length(na.omit(i))})-1),
#                             
#                             census_deaths_f_5 = ifelse(is.na(census.deaths.f.5), 1, census.deaths.f.5), census_deaths_m_5 = ifelse(is.na(census.deaths.m.5), 1, census.deaths.m.5), 
#                             census_deaths_year_idx_5 = match(bf.idx1$interval * floor(as.numeric(colnames(census.deaths.f.5)) / bf.idx1$interval), bf.idx1$periods_out),
#                             census_deaths_group_idx_5 = bf.idx1$ages %/% 5 + 1,
#                             deaths_start_5 = rep(2, ncol(census.deaths.f.5)),
#                             deaths_end_5 = ifelse(apply(census.deaths.f.5, 2, function(i){length(na.omit(i))})-1 > 15, 15, apply(census.deaths.f.5, 2, function(i){length(na.omit(i))})-1),
#                             n_agegrp_5_deaths = bf.idx5$n_ages,
#                             
#                             df = bf5.f.no0.smooth$adjusted, dm = bf5.m.no0.smooth$adjusted,
#                             Ef = bf5.f.no0.smooth$pyears2, Em = bf5.m.no0.smooth$pyears2,
#                             df_age = as.integer(bf5.f.no0.smooth$agegr + 1), dm_age = as.integer(bf5.m.no0.smooth$agegr + 1),
#                             df_time = match(bf5.f.no0.smooth$period, levels(bf5.f.no0.smooth$period)), dm_time = match(bf5.m.no0.smooth$period, levels(bf5.m.no0.smooth$period)),
#                             df_tp = as.integer(bf5.f.no0.smooth$tips), dm_tp = as.integer(bf5.m.no0.smooth$tips),
#                             
#                             df_br = bf5.f.no0.br$event, dm = bf5.m.no0.br$event,
#                             Ef_br = bf5.f.no0.br$pyears, Em = bf5.m.no0.br$pyears,
#                             df_age_br = as.integer(bf5.f.no0.br$agegr + 1), dm_age = as.integer(bf5.m.no0.br$agegr + 1),
#                             df_time_br = match(bf5.f.no0.br$period, levels(bf5.f.no0.br$period)), dm_time = match(bf5.m.no0.br$period, levels(bf5.m.no0.br$period)),
#                            
#                             thiele_age = thiele_age,
#                             D_time = as(A.year, "sparseMatrix"),
#                             
#                             tp_mat_br_m_age0 <- as(model.matrix(~ factor(bf5.m.no0.br$tips, levels = 0:14) - 1) * as.numeric(bf5.m.no0.br$agegr == 0), "sparseMatrix"),
#                             tp_mat_br_m_age1 <- as(model.matrix(~ factor(bf5.m.no0.br$tips, levels = 0:14) - 1) * as.numeric(bf5.m.no0.br$agegr %in% 1:4), "sparseMatrix"),
#                             tp_mat_br_m_age5 <- as(model.matrix(~ factor(bf5.m.no0.br$tips, levels = 0:14) - 1) * as.numeric(bf5.m.no0.br$agegr %in% 5:9), "sparseMatrix"),
#                             tp_mat_br_m_age10 <- as(model.matrix(~ factor(bf5.m.no0.br$tips, levels = 0:14) - 1) * as.numeric(bf5.m.no0.br$agegr >= 10), "sparseMatrix"),
#                             
#                             tp_mat_br_f_age0 <- as(model.matrix(~ factor(bf5.f.no0.br$tips, levels = 0:14) - 1) * as.numeric(bf5.f.no0.br$agegr == 0), "sparseMatrix"),
#                             tp_mat_br_f_age1 <- as(model.matrix(~ factor(bf5.f.no0.br$tips, levels = 0:14) - 1) * as.numeric(bf5.f.no0.br$agegr %in% 1:4), "sparseMatrix"),
#                             tp_mat_br_f_age5 <- as(model.matrix(~ factor(bf5.f.no0.br$tips, levels = 0:14) - 1) * as.numeric(bf5.f.no0.br$agegr %in% 5:9), "sparseMatrix"),
#                             tp_mat_br_f_age10 <- as(model.matrix(~ factor(bf5.f.no0.br$tips, levels = 0:14) - 1) * as.numeric(bf5.f.no0.br$agegr >= 10), "sparseMatrix"),
#                             
#                             penal_tp = as(crossprod(diff(diag(15))), "sparseMatrix"), #only penalising the latter 14 coefficients 
#                             penal_tp_0 = as(diag(c(1, rep(0, 14))), "sparseMatrix"),
#                             null_penal_tp = as(exp(15)*tcrossprod(c(0,1,1,1,rep(0, 11))),"sparseMatrix"),
#                             tp_mean = -2:12,
#                             
#                             D_agetime = as(te.spline, "sparseMatrix"),
#                             D_agetime_fert = as(te.spline.fert, "sparseMatrix"),
#                             
#                             upper_marginal_sd_fx = log(1.1)/1.96,
#                             upper_marginal_sd_gx = 0.08/1.96,
#                             
#                             D_firstrow = A.age[1,1:3],
#                             D_firstrow_ARIMA = ARIMA.vec,
#                             
#                             theta_rho_phi = -log(0.01)/d.rho(0.9, no.basis.time),
#                             theta_rho_lambda = -log(0.01 * 2)/sqrt(-log(1 - 0.9^2)),
#                             theta_rho_fx_age = -log(0.01)/d.rho(0.9, no.basis.time),
#                             theta_rho_fx_time = -log(0.01)/d.rho(0.9, no.basis.time),
#                             theta_rho_gx_age = -log(0.01)/d.rho(0.9, no.basis.time),
#                             theta_rho_gx_time = -log(0.01)/d.rho(0.9, no.basis.time),
#                             
#                             theta_tp = gumbel.theta.tp,
#                             theta_logpop = 3,
#                             
#                             hiv_var_cut = hiv.cut,
#                             
#                             pad=3L,
#                             
#                             upper_marginal_sd2_lambda = log(1.1)/1.96,
#                             upper_marginal_sd2_delta = log(1.1)/1.96,
#                             upper_marginal_sd2_epsilon = log(1.1)/1.96,
#                             
#                             theta_marginal_phi = -log(0.01) *  sqrt(mean(diag(t(D2) %*% solve(tcrossprod(D2)) %*% t(t(D2) %*% solve(tcrossprod(D2)))))) * 1.96 / log(1.2),
#                             spline_trend = spline.trend,
#                             inv_spline_err_mat = inv.spline.err.mat,
#                             
#                             theta_marginal_phi_child_posterior = gumbel.theta.RW2.manual,
#                             
#                             theta_marginal_phi_child_posterior_phi = gumbel.theta.RW2[[1]],
#                             theta_marginal_phi_child_posterior_psi = gumbel.theta.RW2[[2]],
#                             theta_marginal_phi_child_posterior_A = gumbel.theta.RW2[[3]],
#                             theta_marginal_phi_child_posterior_B = gumbel.theta.RW2[[4]],
#                             
#                             theta_marginal_phi_child_posterior_phi_standardised = gumbel.theta.RW2.standardised[[1]],
#                             theta_marginal_phi_child_posterior_psi_standardised = gumbel.theta.RW2.standardised[[2]],
#                             theta_marginal_phi_child_posterior_A_standardised = gumbel.theta.RW2.standardised[[3]],
#                             theta_marginal_phi_child_posterior_B_standardised = gumbel.theta.RW2.standardised[[4]],
#                             
#                             # theta_marginal_phi_child_posterior_phi = gumbel.theta.RW2.manual,
#                             # theta_marginal_phi_child_posterior_psi = gumbel.theta.RW2.manual,
#                             # theta_marginal_phi_child_posterior_A = gumbel.theta.RW2.manual,
#                             # theta_marginal_phi_child_posterior_B = gumbel.theta.RW2.manual,
#                             
#                             RW2_PREC = as(crossprod(diff(diag(no.basis.time), differences = 2)), "sparseMatrix"),
#                             RW2_null = as(1e-5*diag(no.basis.time), "sparseMatrix"),
#                             
#                             theta_rho_phi_base1 = optim(3, get.ar2.rho.alpha.func, method="L-BFGS-B", L = 0.1, alpha = 1-(1 - 0.1)/sqrt(1+0.1)-0.0001)$par,
#                             theta_rho_lambda_base1 = optim(3, get.ar1.rho.alpha.func, method="L-BFGS-B", L = 0.1, alpha = 1-sqrt(1 - 0.1)-0.0001)$par,
#                             theta_rho_fx_age_base1 = optim(3, get.ar2.rho.alpha.func, method="L-BFGS-B", L = 0.1, alpha = 1-(1 - 0.1)/sqrt(1+0.1)-0.0001)$par,
#                             theta_rho_fx_time_base1 = optim(3, get.ar2.rho.alpha.func, method="L-BFGS-B", L = 0.1, alpha = 1-(1 - 0.1)/sqrt(1+0.1)-0.0001)$par,
#                             theta_rho_gx_age_base1 = optim(3, get.ar2.rho.alpha.func, method="L-BFGS-B", L = 0.1, alpha = 0.1)$par,
#                             theta_rho_gx_time_base1 = optim(3, get.ar2.rho.alpha.func, method="L-BFGS-B", L = 0.1, alpha = 0.1)$par,
#                             
#                             cor_mat_m = as(cor.mat.m, "sparseMatrix"),
#                             cor_mat_f = as(cor.mat.f, "sparseMatrix"),
#                             
#                             cor_mat_m_standardised = as(cor.mat.m.standardised, "sparseMatrix"),
#                             cor_mat_f_standardised = as(cor.mat.f.standardised, "sparseMatrix"),
#                             
#                             theta_5q0 = -log(0.01) * 1.96/ log(2),
#                             log_dat_v5q0_m = h.mean.m, log_dat_v5q0_f = h.mean.f,
#                             year_ind_5q0 = match(wpp.5q0.interpolate %>% filter(Sex=="Female", year %in% bf.idx1$periods) %>% .$year, bf.idx1$periods) - 1,
#                             
#                             scale_sd_m = unlist(sd.m),
#                             scale_sd_f = unlist(sd.f)
# )

dat.list <- list(log_basepop_mean_f = log(basepop.f), log_basepop_mean_m = log(basepop.m),
                 log_fx_mean = log_fx_mean,
                 srb = rep(1.05, bf.idx1$n_periods),
                 interval = bf.idx1$interval,
                 n_periods = bf.idx1$n_periods,
                 fx_idx = bf.idx1$fx_idx,
                 n_fx = bf.idx1$n_fx,
                 census_log_pop_f = data.f, census_log_pop_m = data.m,
                 census_year_idx = match(bf.idx1$interval * floor(as.numeric(colnames(data.f)) / bf.idx1$interval), bf.idx1$periods_out),
                 #census_year_grow_idx = as.numeric(colnames(data.f)) - bf.idx1$interval * floor(as.numeric(colnames(data.f)) / bf.idx1$interval),
                 
                 open_idx = bf.idx1$n_ages,
                 #oag = apply(data.f, 2, function(i){length(na.omit(i))}),
                 pop_start = rep(6, ncol(data.f)), 
                 pop_end = ifelse(apply(data.f, 2, function(i){length(na.omit(i))})-1 > 75 + 1, 75 + 1, apply(data.f, 2, function(i){length(na.omit(i))})-1),
                 #pop_end = apply(data.f, 2, function(i){length(na.omit(i))})-1,
                 
                 census_log_pop_f_5 = data.f.5, census_log_pop_m_5 = data.m.5,
                 census_year_idx_5 = match(bf.idx1$interval * floor(as.numeric(colnames(data.f.5)) / bf.idx1$interval), bf.idx1$periods_out),
                 #census_year_group_idx_5 = bf.idx1$ages %/% 5 + 1,
                 #n_agegrp_5 = bf.idx5$n_ages,
                 #pop_start_5 = rep(2, ncol(data.f.5)),
                 #pop_end_5 =  ifelse(apply(data.f.5, 2, function(i){length(na.omit(i))})-1 > 15, 15, apply(data.f.5, 2, function(i){length(na.omit(i))})-1),
                 
                 census_year_group_idx_5 = bf.idx1$ages + 1,
                 n_agegrp_5 = bf.idx1$n_ages,
                 pop_start_5 = rep(6, ncol(data.f.5)),
                 pop_end_5 =  ifelse(apply(data.f.5, 2, function(i){length(na.omit(i))})-1 > 75 + 1, 75 + 1, apply(data.f.5, 2, function(i){length(na.omit(i))})-1),
                 
                 census_deaths_f = census.deaths.f, census_deaths_m = census.deaths.m, 
                 census_deaths_year_idx = match(bf.idx1$interval * floor(as.numeric(colnames(census.deaths.f)) / bf.idx1$interval), bf.idx1$periods_out),
                 deaths_start = rep(6, ncol(census.deaths.f)),
                 deaths_end = ifelse(apply(census.deaths.f, 2, function(i){length(na.omit(i))})-1 > 75 + 1, 75 + 1, apply(census.deaths.f, 2, function(i){length(na.omit(i))})-1),
                 
                 census_deaths_f_5 = ifelse(is.na(census.deaths.f.5), 1, census.deaths.f.5), census_deaths_m_5 = ifelse(is.na(census.deaths.m.5), 1, census.deaths.m.5), 
                 census_deaths_year_idx_5 = match(bf.idx1$interval * floor(as.numeric(colnames(census.deaths.f.5)) / bf.idx1$interval), bf.idx1$periods_out),
                 census_deaths_group_idx_5 = bf.idx1$ages %/% 5 + 1,
                 deaths_start_5 = rep(2, ncol(census.deaths.f.5)),
                 deaths_end_5 = ifelse(apply(census.deaths.f.5, 2, function(i){length(na.omit(i))})-1 > 15, 15, apply(census.deaths.f.5, 2, function(i){length(na.omit(i))})-1),
                 n_agegrp_5_deaths = bf.idx5$n_ages,
                 
                 df = bf5.f.no0.smooth$adjusted, dm = bf5.m.no0.smooth$adjusted,
                 Ef = bf5.f.no0.smooth$pyears2, Em = bf5.m.no0.smooth$pyears2,
                 df_age = as.integer(bf5.f.no0.smooth$agegr + 1), dm_age = as.integer(bf5.m.no0.smooth$agegr + 1),
                 df_time = match(bf5.f.no0.smooth$period, levels(bf5.f.no0.smooth$period)), dm_time = match(bf5.m.no0.smooth$period, levels(bf5.m.no0.smooth$period)),
                 df_tp = as.integer(bf5.f.no0.smooth$tips), dm_tp = as.integer(bf5.m.no0.smooth$tips),
                 
                 df_br = bf5.f.no0.br$event, dm_br = bf5.m.no0.br$event,
                 Ef_br = bf5.f.no0.br$pyears, Em_br = bf5.m.no0.br$pyears,
                 df_age_br = as.integer(bf5.f.no0.br$agegr + 1), dm_age_br = as.integer(bf5.m.no0.br$agegr + 1),
                 df_time_br = match(bf5.f.no0.br$period, levels(bf5.f.no0.br$period)), dm_time_br = match(bf5.m.no0.br$period, levels(bf5.m.no0.br$period)),
                 
                 thiele_age = thiele_age,
                 D_time = as(A.year, "sparseMatrix"),
                 D_err =  as(A.err.const[bf.idx1$periods - 1900 + 1,], "dgCMatrix"),
                 slope_ind = bf.idx1$periods - mean(bf.idx1$periods),
                 
                 tp_mat_br_m_age0 = as(model.matrix(~ factor(bf5.m.no0.br$tips, levels = 0:14) - 1) * as.numeric(bf5.m.no0.br$agegr == 0), "sparseMatrix"),
                 tp_mat_br_m_age1 = as(model.matrix(~ factor(bf5.m.no0.br$tips, levels = 0:14) - 1) * as.numeric(bf5.m.no0.br$agegr %in% 1:4), "sparseMatrix"),
                 tp_mat_br_m_age5 = as(model.matrix(~ factor(bf5.m.no0.br$tips, levels = 0:14) - 1) * as.numeric(bf5.m.no0.br$agegr %in% 5:9), "sparseMatrix"),
                 tp_mat_br_m_age10 = as(model.matrix(~ factor(bf5.m.no0.br$tips, levels = 0:14) - 1) * as.numeric(bf5.m.no0.br$agegr >= 10), "sparseMatrix"),
                 
                 tp_mat_br_f_age0 = as(model.matrix(~ factor(bf5.f.no0.br$tips, levels = 0:14) - 1) * as.numeric(bf5.f.no0.br$agegr == 0), "sparseMatrix"),
                 tp_mat_br_f_age1 = as(model.matrix(~ factor(bf5.f.no0.br$tips, levels = 0:14) - 1) * as.numeric(bf5.f.no0.br$agegr %in% 1:4), "sparseMatrix"),
                 tp_mat_br_f_age5 = as(model.matrix(~ factor(bf5.f.no0.br$tips, levels = 0:14) - 1) * as.numeric(bf5.f.no0.br$agegr %in% 5:9), "sparseMatrix"),
                 tp_mat_br_f_age10 = as(model.matrix(~ factor(bf5.f.no0.br$tips, levels = 0:14) - 1) * as.numeric(bf5.f.no0.br$agegr >= 10), "sparseMatrix"),
                 
                 penal_tp = as(crossprod(diff(diag(15))), "dgCMatrix"), #only penalising the latter 14 coefficients 
                 penal_tp_0 = as(diag(c(1, rep(0, 14))), "dgCMatrix"),
                 null_penal_tp = as(exp(15)*tcrossprod(c(0,1,1,1,rep(0, 11))),"dgCMatrix"),
                 tp_mean = -2:12,
                 
                 null_penal_tp_br = as(exp(15) * tcrossprod(c(rep(0, 6), 1, 1, 1, rep(0, 6))),"dgCMatrix"),
                 tp_mean_br = -7:7,
                
                 log_spline_prior_means_m = log.spline.prior.m,
                 log_spline_prior_means_f = log.spline.prior.f,

                 # RW2_PREC = as(crossprod(diff(diag(no.basis.time), differences = 2)), "dgCMatrix"),
                 # RW2_null = as(1e-6 * diag(no.basis.time), "dgCMatrix"),

                 theta_rho_hump = theta_rho,
                 
                 V_mf_mean = V_mf_mean + (1960 - mean(bf.idx1$periods))^2 * V_mf_slope,
                 V_mf_slope = V_mf_slope,
                 V_mf = V_mf,
                 
                 D_agetime = as(te.spline, "sparseMatrix"),
                 D_agetime_fert = as(te.spline.fert, "sparseMatrix"),
                 D_firstrow = A.age[1,1:3],
                 D_firstrow_ARIMA = ARIMA.vec,
                 
                 theta_rho_fx_age = -log(0.01)/d.rho(0.9, no.basis.time),
                 theta_rho_fx_time = -log(0.01)/d.rho(0.9, no.basis.time),
                 theta_rho_gx_age = -log(0.01)/d.rho(0.9, no.basis.time),
                 theta_rho_gx_time = -log(0.01)/d.rho(0.9, no.basis.time),
                 
                 upper_marginal_sd_fx = log(1.2)/1.96,
                 upper_marginal_sd_gx = 0.2/1.96,
                 
                 theta_tp = -log(0.01) * 1.96/ 0.5,
                 theta_logpop = -log(0.01) * 1.96/ log(1.8),
                 
                 theta_rho_fx_age_base1 = optim(3, get.ar2.rho.alpha.func, method="L-BFGS-B", L = 0.1, alpha = 1-(1 - 0.1)/sqrt(1+0.1)-0.001)$par,
                 theta_rho_fx_time_base1 = optim(3, get.ar2.rho.alpha.func, method="L-BFGS-B", L = 0.1, alpha = 1-(1 - 0.1)/sqrt(1+0.1)-0.001)$par,
                 theta_rho_gx_age_base1 = optim(3, get.ar2.rho.alpha.func, method="L-BFGS-B", L = 0.1, alpha = 1-(1 - 0.1)/sqrt(1+0.1)-0.001)$par,
                 theta_rho_gx_time_base1 = optim(3, get.ar2.rho.alpha.func, method="L-BFGS-B", L = 0.1, alpha = 1-(1 - 0.1)/sqrt(1+0.1)-0.001)$par,
                 # theta_rho_gx_age_base1 = optim(3, get.ar2.rho.alpha.func, method="L-BFGS-B", L = 0.1, alpha = 0.01)$par,
                 # theta_rho_gx_time_base1 = optim(3, get.ar2.rho.alpha.func, method="L-BFGS-B", L = 0.1, alpha = 0.01)$par,
                 
                 theta_5q0 = -log(0.1) * 1.96 / log(1.4),
                 log_dat_v5q0_m = h.mean.m, log_dat_v5q0_f = h.mean.f,
                 year_ind_5q0 = match(wpp.5q0.interpolate %>% filter(Sex=="Female", year %in% bf.idx1$periods, !year %in% unique(bf5.smooth.br$period)) %>% .$year, bf.idx1$periods) - 1,
                 
                 loghump = ifelse(params$loghump, 1L, 0L),
                 # hump_arima_initial = 0L,
                 
                 lambda_prior_sd = 0.3,
                 delta_prior_sd = 0.1,
                 epsilon_prior_sd = 0.1,
                 
                 # log_V_factor_prior_sd = c(10, 5, 1)
                 log_V_factor_prior_sd = c(1, 1, 1),
                 
                 calc_5q0 = as.integer(!params$br.data),
                 # calc_5q0 = 1L,
                 
                 PCA = smoothed.pca)

par.list <- list(log_basepop_f = log(basepop.f), log_basepop_m = log(basepop.m),
                 log_fx_spline_params = rep(0, no.basis.fert * no.basis.time),
                 gx_f_spline_params = rep(0, no.basis.time * no.basis.age), gx_m_spline_params = rep(0, no.basis.time * no.basis.age),
                 
                 #log_tau2_logpop = c(log(1.96^2/log(1.2)^2), log(1.96^2/log(1.5)^2), log(1.96^2/log(1.2)^2), log(1.96^2/log(1.5)^2)),
                 log_tau2_logpop = c(log(1.96^2/log(2)^2), log(1.96^2/log(1.5)^2)),
                 
                 log_dispersion = 1,
                 log_dispersion_br = 1,
                 
                 # log_lambda_tp = 5,
                 # tp_params = tp.init,
                 # tp_slope = 0,
                 # tp_params_5 = 0,
                 # tp_params_10 = 0,
                 # 
                 # log_lambda_tp_br_age0 = 5,
                 # tp_params_br_age0 = tp.init,
                 # tp_slope_br_age0 = 0,
                 # tp_params_5_br_age0 = 0,
                 # tp_params_10_br_age0 = 0,
                 # 
                 # log_lambda_tp_br_age1 = 5,
                 # tp_params_br_age1 = tp.init,
                 # tp_slope_br_age1 = 0,
                 # tp_params_5_br_age1 = 0,
                 # tp_params_10_br_age1 = 0,
                 # 
                 # log_lambda_tp_br_age5 = 5,
                 # tp_params_br_age5 = tp.init,
                 # tp_slope_br_age5 = 0,
                 # tp_params_5_br_age5 = 0,
                 # tp_params_10_br_age5 = 0,
                 # 
                 # log_lambda_tp_br_age10 = 5,
                 # tp_params_br_age10 = tp.init,
                 # tp_slope_br_age10 = 0,
                 # tp_params_5_br_age10 = 0,
                 # tp_params_10_br_age10 = 0,
                 
                 log_spline_means_m = c(-4, log.spline.prior.m[c(2,3,7,8)]),
                 log_spline_means_f = c(-4, log.spline.prior.f[c(2,3,7,8)]),
                 
                 log_spline_slopes_m = rep(0, 5),
                 log_spline_slopes_f = rep(0, 5),
                 
                 log_spline_err_m = matrix(0, ncol(A.err.const), 5),
                 log_spline_err_f = matrix(0, ncol(A.err.const), 5),
                 
                 log_spline_hump_m = matrix(log.spline.prior.m[4:6], 3, ncol(A.year)) %>% t,
                 log_spline_hump_f = matrix(log.spline.prior.f[4:6], 3, ncol(A.year)) %>% t,
                 
                 log_prec_all_RW2_hump = c(3, 4, 4),
                 logit_rho_lambda = 1,
                 logit_rho_delta = 1,
                 logit_rho_epsilon = 1,
                 
                 logit_rho_fx_age = 0,
                 logit_rho_fx_time = 0,
                 logit_rho_gx_age = 0,
                 logit_rho_gx_time = 0,
                 
                 log_marginal_lambda_gx = log((gumbel.theta.AR2.marginal.gx/-log(0.01))^2) + 1,
                 log_marginal_lambda_fx = log((gumbel.theta.AR2.marginal.fx/-log(0.01))^2) + 1,
                 
                 log_dispersion_census_deaths = c(0, 0),
                 log_dispersion_census_deaths_5 = c(0, 0),
                 
                 log_prec_5q0 = log(1.96^2/log(1.5)^2),
                 
                 log_V_factor_mean = log(1),
                 log_V_factor_slope = log(1),
                 log_V_factor = log(1)) %>%
  append(tp.list)
  #               
# TMB.fun <- MakeADFun(data = dat.list,
#                      parameters = par.list,
#                      
#                      random = c("log_basepop_f","log_basepop_m",
#                                 "log_fx_spline_params",
#                                 "gx_f_spline_params","gx_m_spline_params",
#                                 
#                                 "tp_params", "tp_slope",
#                                 "tp_params_5", "tp_params_10",
#                                 
#                                 "tp_params_br_age0", "tp_slope_br_age0",
#                                 "tp_params_5_br_age0", "tp_params_10_br_age0",
#                                 
#                                 "tp_params_br_age1", "tp_slope_br_age1",
#                                 "tp_params_5_br_age1", "tp_params_10_br_age1",
#                                 
#                                 "tp_params_br_age5", "tp_slope_br_age5",
#                                 "tp_params_5_br_age5", "tp_params_10_br_age5",
#                                 
#                                 "tp_params_br_age10", "tp_slope_br_age10",
#                                 "tp_params_5_br_age10", "tp_params_10_br_age10",
#                                 
#                                 "log_spline_means_m", "log_spline_means_f",
#                                 "log_spline_slopes_m", "log_spline_slopes_f",
#                                 "log_spline_err_m", "log_spline_err_f",
#                                 "log_spline_hump_m", "log_spline_hump_f"),
#                      
#                      DLL = "ccmpp_bothsexes_thiele_normalhump_HMD_trend_and_err",
#                      
#                      map = list(
#                        log_V_factor_mean = factor(NA),
#                        log_V_factor_slope = factor(NA),
#                        log_V_factor = factor(NA)
#                        # log_basepop_f = factor(matrix(NA, dim(basepop.f)[1], dim(basepop.f)[2])),
#                        # log_basepop_m = factor(matrix(NA, dim(basepop.f)[1], dim(basepop.f)[2]))
#                      ),
#                      
#                      silent = F)
# 
# f <- stats::nlminb(TMB.fun$par, TMB.fun$fn, TMB.fun$gr,
#                    control = list(trace = 1,
#                                   iter.max = 10000,
#                                   eval.max = 10000))
# 
# f.report <- TMB.fun$report(TMB.fun$env$last.par.best)
# 
# save(TMB.fun, f, f.report, file = 'fit Namibia new.RData')


fit <- fit_tmb(data = dat.list,
               
               par = par.list,
               
               inner_verbose=TRUE,
               random = c("log_basepop_f","log_basepop_m",
                          "log_fx_spline_params",
                          "gx_f_spline_params","gx_m_spline_params",
                          
                          "tp_params", "tp_slope",
                          "tp_params_5", "tp_params_10",
                          
                          "tp_params_br_age0", "tp_slope_br_age0",
                          "tp_params_5_br_age0", "tp_params_10_br_age0",
                          
                          "tp_params_br_age1", "tp_slope_br_age1",
                          "tp_params_5_br_age1", "tp_params_10_br_age1",
                          
                          "tp_params_br_age5", "tp_slope_br_age5",
                          "tp_params_5_br_age5", "tp_params_10_br_age5",
                          
                          "tp_params_br_age10", "tp_slope_br_age10",
                          "tp_params_5_br_age10", "tp_params_10_br_age10",
                          
                          "log_spline_means_m", "log_spline_means_f",
                          "log_spline_slopes_m", "log_spline_slopes_f",
                          "log_spline_err_m", "log_spline_err_f",
                          "log_spline_hump_m", "log_spline_hump_f"),
               
               DLL = "ccmpp_bothsexes_thiele_normalhump_HMD_trend_and_err_BPCA_ARIMAhump",
               
               map = list(
                 # log_V_factor_mean = factor(NA),
                 log_V_factor_slope = factor(NA)
                 # log_V_factor = factor(NA)
                 # log_basepop_f = factor(matrix(NA, dim(basepop.f)[1], dim(basepop.f)[2])),
                 # log_basepop_m = factor(matrix(NA, dim(basepop.f)[1], dim(basepop.f)[2]))
               )
)

save(fit, file = paste(params$country, ' BPCA.RData'))

#variance estimates####
rmvnorm_sparseprec <- function(n, mean = rep(0, nrow(prec)), prec = diag(lenth(mean))) {
  z = matrix(rnorm(n * length(mean)), ncol = n)
  L_inv = Matrix::Cholesky(prec)
  v <- mean + Matrix::solve(as(L_inv, "pMatrix"), Matrix::solve(Matrix::t(as(L_inv, "Matrix")), z))
  as.matrix(Matrix::t(v))
}  

var.sample <- function(fit, nsample){
  r <- fit$obj$env$random
  par_f <- fit$par.full[-r]
  
  par_r <- fit$par.full[r]
  hess_r <- fit$obj$env$spHess(fit$par.full, random = TRUE)
  smp_r <- rmvnorm_sparseprec(nsample, par_r, hess_r)
  
  smp <- matrix(0, nsample, length(fit$par.full))
  smp[ , r] <- smp_r
  smp[ ,-r] <- matrix(par_f, nsample, length(par_f), byrow = TRUE)
  colnames(smp) <- rep("NA", length(fit$par.full))
  colnames(smp)[r] <- names(par_r)
  colnames(smp)[-r] <- names(par_f)
  smp
}

fit.par.sim <- var.sample(fit, 1000)
fit.var.sim <-  apply(fit.par.sim, 1, fit$obj$report)

fit.var.sim.CI.list <- list(lower = list(),
                            mean = list(),
                            upper = list())

temp <- c()

for(i in 1:length(fit.var.sim[[1]])){
  for(j in 1:length(fit.var.sim)){
    temp <- cbind(temp, c(fit.var.sim[[j]][[i]]))
  }
  
  fit.var.sim.CI.list[[1]][[ names(fit.var.sim[[1]])[i] ]] <- array(apply(temp, 1, quantile, probs = c(0.025)), 
                                                                    `if`(is.null(dim(fit.var.sim[[1]][[i]])), length(fit.var.sim[[1]][[i]]), dim(fit.var.sim[[1]][[i]])), 
                                                                    dimnames = dimnames(fit.var.sim[[1]][[i]]))
  
  fit.var.sim.CI.list[[2]][[ names(fit.var.sim[[1]])[i] ]] <- array(apply(temp, 1, quantile, probs = c(0.5)), 
                                                                    `if`(is.null(dim(fit.var.sim[[1]][[i]])), length(fit.var.sim[[1]][[i]]), dim(fit.var.sim[[1]][[i]])), 
                                                                    dimnames = dimnames(fit.var.sim[[1]][[i]]))
  
  fit.var.sim.CI.list[[3]][[ names(fit.var.sim[[1]])[i] ]] <- array(apply(temp, 1, quantile, probs = c(0.975)), 
                                                                    `if`(is.null(dim(fit.var.sim[[1]][[i]])), length(fit.var.sim[[1]][[i]]), dim(fit.var.sim[[1]][[i]])), 
                                                                    dimnames = dimnames(fit.var.sim[[1]][[i]]))
  
  temp <- c()
}

fit.par.sim.CI <- apply(fit.par.sim, 2, quantile, probs = c(0.025, 0.5, 0.975))

fit.par.sim.CI.list <- list(lower = split(unname(fit.par.sim.CI[1,]), names(fit.par.sim.CI[1,])),
                            mean = split(unname(fit.par.sim.CI[2,]), names(fit.par.sim.CI[2,])),
                            upper = split(unname(fit.par.sim.CI[3,]), names(fit.par.sim.CI[3,])))

for(i in 1:length(fit.par.sim.CI.list)){
  for(j in 1:length(fit.par.sim.CI.list[[1]])){
    fit.par.sim.CI.list[[i]][[j]] <- array(fit.par.sim.CI.list[[i]][[j]], dim = dim(as.matrix(fit$obj$env$parList()[[names(fit.par.sim.CI.list[[i]][j])]])))
  }
}

fit$fit.par.sim.CI.list <- fit.par.sim.CI.list
fit$fit.var.sim.CI.list <- fit.var.sim.CI.list
fit$fit.par.sim <- fit.par.sim
fit$fit.var.sim <- fit.var.sim

saveRDS(fit, file = paste0(params$country," BPCA no 5q0.rds"))

#plot par estimates####
fit <- readRDS(paste0(params$country, ' BPCA 5q0.rds'))

par.df <- bind_rows(
  fit$mode$log_par_f %>%
    `colnames<-`(c("phi", "d", "psi", "lambda", "delta", "epsilon", "A", "B")) %>%
    as_tibble() %>%
    mutate(sex = "female",
           year = bf.idx1$periods),
  
  fit$mode$log_par_m %>%
    `colnames<-`(c("phi", "d", "psi", "lambda", "delta", "epsilon", "A", "B")) %>%
    as_tibble() %>%
    mutate(sex = "male",
           year = bf.idx1$periods)
  ) %>%
  mutate(across(phi:epsilon, exp)) %>%
  pivot_longer(phi:B) %>%
  full_join(
    bind_rows(
      fit$fit.var.sim.CI.list$upper$log_par_f %>%
        `colnames<-`(c("phi", "d", "psi", "lambda", "delta", "epsilon", "A", "B")) %>%
        as_tibble() %>%
        mutate(sex = "female",
               year = bf.idx1$periods),
      
      fit$fit.var.sim.CI.list$upper$log_par_m %>%
        `colnames<-`(c("phi", "d", "psi", "lambda", "delta", "epsilon", "A", "B")) %>%
        as_tibble() %>%
        mutate(sex = "male",
               year = bf.idx1$periods)
    ) %>%
      mutate(across(phi:epsilon, exp)) %>%
      pivot_longer(phi:B) %>%
      rename(upper = value) %>%
      full_join(
        bind_rows(
          fit$fit.var.sim.CI.list$lower$log_par_f %>%
            `colnames<-`(c("phi", "d", "psi", "lambda", "delta", "epsilon", "A", "B")) %>%
            as_tibble() %>%
            mutate(sex = "female",
                   year = bf.idx1$periods),
          
          fit$fit.var.sim.CI.list$lower$log_par_m %>%
            `colnames<-`(c("phi", "d", "psi", "lambda", "delta", "epsilon", "A", "B")) %>%
            as_tibble() %>%
            mutate(sex = "male",
                   year = bf.idx1$periods)
        ) %>%
          mutate(across(phi:epsilon, exp)) %>%
          pivot_longer(phi:B) %>%
          rename(lower = value)
      )) %>%
  mutate(name = fct_relevel(name, "phi", "d", "psi", "lambda", "delta", "epsilon", "A", "B"))

par.df %>%
  ggplot() +
  geom_line(aes(x = year, y = value, col = sex)) +
  geom_ribbon(aes(x = year, ymax = upper, ymin = lower, fill = sex), alpha = 0.2) +
  facet_wrap(~name, scales = "free_y")

#plot tips####
cbind(fit$mode$tp_params,
      fit$mode$tp_params_br_age0,
      fit$mode$tp_params_br_age1,
      fit$mode$tp_params_br_age5,
      fit$mode$tp_params_br_age10) %>%
  `colnames<-`(c("Adult", "Age 0", "Age 1-4", "Age 5-9", "Age 10+")) %>%
  `rownames<-`(0:14) %>%
  reshape2::melt() %>%
  ggplot() +
  geom_line(aes(x = Var1, y = value)) +
  coord_flip() +
  facet_grid(~Var2)

# if(length(unique(bf5.smooth$DHS)) == 1 || diff(as.numeric(str_extract(unique(bf5.smooth$DHS), "\\d+"))) > 14) {
#   load("~/more countries final avg sex Rwanda.RData")
#   skip<-c("Ethiopia","Central African Republic","Comoros","Sao Tome and Principe","Botswana","Cape Verde","Equatorial Guinea","Eritrea","Nigeria (Ondo State)","Ghana","Mauritania","Sudan")
#   joint.countries<-names(aggr.mat.cohort.0)[!names(aggr.mat.cohort.0)%in%skip]
#   
#   all.list<-function(x,age.start,age.end,year.start,year.end){
#     no.basis = 15
#     knots<-seq(0,1,length=no.basis-2)
#     dk<-knots[2]-knots[1]	
#     knots<-c(knots[1]-dk*(3:1),knots,knots[no.basis-2]+dk*(1:3))
#     age<-seq(0,1,length=age.end-age.start+1)
#     year<-seq(0,1,length=year.end-year.start+1)
#     
#     A.age<-c()
#     for(j in 1:no.basis) {
#       A.age<-cbind(A.age,bspline(age,knots,j))
#     }
#     
#     A.year<-c()
#     for(j in 1:no.basis) {
#       A.year<-cbind(A.year,bspline(year,knots,j))
#     }
#     
#     te.spline<-A.year%x%A.age
#     
#     ind<-function(x){
#       aggr.mat.reduced<-x
#       data.mat.m<-aggr.mat.reduced[aggr.mat.reduced$mm1=="male" & aggr.mat.reduced$agegr %in% age.start:age.end & aggr.mat.reduced$period %in% year.start:year.end,c(7,5,3,2,4)]
#       data.mat.f<-aggr.mat.reduced[aggr.mat.reduced$mm1=="female" & aggr.mat.reduced$agegr %in% age.start:age.end & aggr.mat.reduced$period %in% year.start:year.end,c(7,5,3,2,4)]
#       
#       age.start.col<-min(which(A.age[min(data.mat.m[,3])-age.start+1,]!=0))
#       age.end.col<-max(which(A.age[max(data.mat.m[,3])-age.start+1,]!=0))
#       year.start.col<-min(which(A.year[min(data.mat.m[,4])-year.start+1,]!=0))
#       year.end.col<-max(which(A.year[max(data.mat.m[,4])-year.start+1,]!=0))
#       age.start<-min(data.mat.m[,3])
#       age.end<-max(data.mat.m[,3])
#       year.start<-min(data.mat.m[,4])
#       year.end<-max(data.mat.m[,4])
#       tp<-max(data.mat.m[,5],data.mat.f[,5])
#       
#       ind.place<-c(age.start.col,age.end.col,year.start.col,year.end.col,age.start,age.end,year.start,year.end,tp)
#       names(ind.place)<-c("age.start.col", "age.end.col","year.start.col","year.end.col","age.start","age.end","year.start","year.end","tp")
#       ind.place
#     }  
#     
#     
#     everything.func<-function(i){
#       age.seq<-as.numeric(ind.list[[i]][["age.start"]]:ind.list[[i]][["age.end"]])
#       year.seq<-as.numeric(ind.list[[i]][["year.start"]]:ind.list[[i]][["year.end"]])
#       
#       age.m<-A.age[age.seq-age.start+1,]%*%joint.countries.age.m[[i]]
#       age.f<-A.age[age.seq-age.start+1,]%*%joint.countries.age.f[[i]]
#       time.m<-A.year[year.seq-year.start+1,]%*%joint.countries.time.m[[i]]
#       time.f<-A.year[year.seq-year.start+1,]%*%joint.countries.time.f[[i]]
#       agetime.m<-matrix(te.spline[age.seq-age.start+1+rep((age.end-age.start+1)*(year.seq-year.start),each=length(age.seq)),]%*%joint.countries.2d.m[[i]],length(age.seq),length(year.seq))
#       agetime.f<-matrix(te.spline[age.seq-age.start+1+rep((age.end-age.start+1)*(year.seq-year.start),each=length(age.seq)),]%*%joint.countries.2d.f[[i]],length(age.seq),length(year.seq))
#       
#       mort.m<-joint.countries.avg.m[[i]]+t(apply(apply(agetime.m,2,function(x){x+age.m}),1,function(x){x+time.m}))
#       mort.f<-joint.countries.avg.f[[i]]+t(apply(apply(agetime.f,2,function(x){x+age.f}),1,function(x){x+time.f}))
#       rownames(age.m)<-rownames(age.f)<-rownames(agetime.m)<-rownames(agetime.f)<-rownames(mort.m)<-rownames(mort.f)<-age.seq
#       rownames(time.m)<-rownames(time.f)<-colnames(agetime.m)<-colnames(agetime.f)<-colnames(mort.m)<-colnames(mort.f)<-year.seq
#       
#       
#       if(i=="Rwanda"){
#         mort.m[,4:6]<-mort.m[,4:6]+matrix(rep.par.list$rwanda_geno_m,56,3,byrow=T)
#         mort.f[,4:6]<-mort.f[,4:6]+matrix(rep.par.list$rwanda_geno_f,56,3,byrow=T)
#       }
#       
#       tips<-joint.countries.tp.c[[i]]+tp.common
#       
#       list(tp=tips,avg.m=joint.countries.avg.m[[i]],avg.f=joint.countries.avg.f[[i]],age.m=age.m,age.f=age.f,time.m=time.m,time.f=time.f,agetime.m=agetime.m,agetime.f=agetime.f,mort.m=mort.m,mort.f=mort.f)
#     }
#     
#     ind.list<-lapply(aggr.mat.cohort.0[joint.countries],ind)
#     rep.par.list<-split(unname(x$env$last.par.best),names(x$env$last.par.best))
#     joint.countries.avg.m<-split(rep.par.list$avg_m,joint.countries)
#     joint.countries.avg.f<-split(rep.par.list$avg_f,joint.countries)
#     joint.countries.age.m<-split(rep.par.list$spline_params_m_age,rep(joint.countries,each=no.basis))
#     joint.countries.age.f<-split(rep.par.list$spline_params_f_age,rep(joint.countries,each=no.basis))
#     joint.countries.time.m<-split(rep.par.list$spline_params_m_time,rep(joint.countries,each=no.basis))
#     joint.countries.time.f<-split(rep.par.list$spline_params_f_time,rep(joint.countries,each=no.basis))
#     joint.countries.2d.m<-split(rep.par.list$spline_params_m_2d,rep(joint.countries,each=no.basis*no.basis))
#     joint.countries.2d.f<-split(rep.par.list$spline_params_f_2d,rep(joint.countries,each=no.basis*no.basis))
#     
#     joint.countries.tp.c<-split(rep.par.list$tips_params,rep(joint.countries,each=length(rep.par.list$tips_params_common)))
#     tp.common<-rep.par.list$tips_params_common
#     tp.common[6]<-tp.common[6]+rep.par.list$tp_common_5
#     tp.common[11]<-tp.common[11]+rep.par.list$tp_common_10
#     
#     
#     setNames(lapply(joint.countries,everything.func),joint.countries)
#   }
#   
#   more.countries.avg<-all.list(tmb.full.joint.common,age.start=10,age.end=65,year.start=1990,year.end=2017)
#   
#   tp.all<- exp(more.countries.avg[[params$country]]$tp) %>% as_tibble() %>%
#     mutate(tips = 0:14)
# } else {
#   tp.all <- apply(thiele.par.sim, 1, function(i){
#     tp <- i %>% split(names(.)) %>% .$tp_params
#     tp[6] <- tp[6] + i %>% split(names(.)) %>% .$tp_params_5
#     tp[11] <- tp[11] + i %>% split(names(.)) %>% .$tp_params_10
#     exp(tp)
#   }) %>%
#     `rownames<-`(0:14) %>%
#     reshape2::melt() %>% as_tibble() %>% 
#     select(tips = Var1, sim = Var2, value = value)}
# 
# if(length(unique(bf5.smooth$DHS)) == 1 || diff(as.numeric(str_extract(unique(bf5.smooth$DHS), "\\d+"))) > 14) {
#   tp.all %>%
#     ggplot + geom_line(aes(x = tips, y = value)) + coord_flip()
# } else{
#   tp.all %>%
#     ggplot() + stat_boxplot(geom="errorbar", aes(y = tips, x = value, group = tips), position=position_dodge(0), lwd = 0.8) +
#     geom_boxplot(aes(y = tips, x = value, group = tips), lwd = 0.8) +
#     ggtitle(paste(country, "Estimated TiPS")) +
#     theme_bw() +
#     theme(text = element_text(size=25), legend.key.size = unit(1.5,"cm"), strip.text = element_text(size=30),
#           plot.title = element_text(hjust = 0.5, face = "bold", size = 35))
# }

#plot mort schedules vs DHS spline####
DHS.plot <- aggr.mat.cohort.0[[params$country]] %>%
  filter(agegr %in% 15:59) %>%
  group_by(country, mm1, tips, DHS, agegr, period) %>%
  summarise_at(vars(pyears, event, pyears2, adjusted), sum) %>%
  ungroup() %>%
  mutate(rate = adjusted / pyears2) %>%
  select(country, mm1, tips, DHS, agegr, period, rate, adjusted, pyears2) %>%
  rename(sex = mm1,
         age = agegr,
         year = period) %>%
  full_join(aggr.mat.br[[params$country]] %>%
              mutate(country = params$country) %>%
              filter(agegr %in% 0:20) %>%
              group_by(country, b4, tips, DHS, agegr, period) %>%
              summarise_at(vars(pyears, event), sum) %>%
              ungroup() %>%
              mutate(rate.br = event / pyears,
                     tips = as.numeric(levels(tips))[tips]) %>%
              select(country, b4, tips, DHS, agegr, period, rate.br, event, pyears) %>%
              rename(sex = b4,
                     age = agegr,
                     year = period)) %>%
  full_join(
    bind_rows(
      fit$mode$mx_mat_f %>%
        `colnames<-`(bf.idx1$periods) %>%
        `rownames<-`(1:nrow(fit$mode$mx_mat_f) - 1) %>%
        reshape2::melt() %>%
        mutate(sex = "female") %>%
        rename(age = Var1,
               year = Var2,
               est = value),
      
      fit$mode$mx_mat_m %>%
        `colnames<-`(bf.idx1$periods) %>%
        `rownames<-`(1:nrow(fit$mode$mx_mat_f) - 1) %>%
        reshape2::melt() %>%
        mutate(sex = "male") %>%
        rename(age = Var1,
               year = Var2,
               est = value)
      )
    ) %>%
  full_join(
    bind_rows(
      fit$fit.var.sim.CI.list$upper$mx_mat_f %>%
        `colnames<-`(bf.idx1$periods) %>%
        `rownames<-`(1:nrow(fit$mode$mx_mat_f) - 1) %>%
        reshape2::melt() %>%
        mutate(sex = "female") %>%
        rename(age = Var1,
               year = Var2,
               upper = value),
      
      fit$fit.var.sim.CI.list$upper$mx_mat_m %>%
        `colnames<-`(bf.idx1$periods) %>%
        `rownames<-`(1:nrow(fit$mode$mx_mat_f) - 1) %>%
        reshape2::melt() %>%
        mutate(sex = "male") %>%
        rename(age = Var1,
               year = Var2,
               upper = value)
    )
  ) %>%
  full_join(
    bind_rows(
      fit$fit.var.sim.CI.list$lower$mx_mat_f %>%
        `colnames<-`(bf.idx1$periods) %>%
        `rownames<-`(1:nrow(fit$mode$mx_mat_f) - 1) %>%
        reshape2::melt() %>%
        mutate(sex = "female") %>%
        rename(age = Var1,
               year = Var2,
               lower = value),
      
      fit$fit.var.sim.CI.list$lower$mx_mat_m %>%
        `colnames<-`(bf.idx1$periods) %>%
        `rownames<-`(1:nrow(fit$mode$mx_mat_f) - 1) %>%
        reshape2::melt() %>%
        mutate(sex = "male") %>%
        rename(age = Var1,
               year = Var2,
               lower = value)
    )
  )

DHS.plot %>%
  filter(year %in% seq(min(year), max(year), length = 6)) %>%
  ggplot() +
  geom_line(aes(x = age, y = est, col = sex)) +
  scale_y_continuous(trans = "log") +
  facet_wrap(~year)

par.mat <- par.df %>%
  select(-upper, -lower) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  arrange(sex, year) %>%
  select(phi:B)

component.df.all <- full_join(
  par.df %>%
    select(-upper, -lower) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    arrange(sex, year) %>%
    select(1:2) %>%
    bind_cols(
      apply(par.mat, 1, function(i){i[1] ^ ((0:95 + i[2])^i[3])}) %>%
        t %>%
        `colnames<-`(0:95)) %>%
    pivot_longer(cols = -(1:2), names_to = "age") %>%
    rename(Child = value),
  
  par.df %>%
    select(-upper, -lower) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    arrange(sex, year) %>%
    select(1:2) %>%
    bind_cols(
      apply(par.mat, 1, function(i){i[4] * exp(-i[5] * (log(c(0.01, 1:95)) - log(i[6]))^2 )}) %>%
        t %>%
        `colnames<-`(0:95)) %>%
    pivot_longer(cols = -(1:2), names_to = "age") %>%
    rename(Hump = value)) %>%
  full_join(
    par.df %>%
      select(-upper, -lower) %>%
      pivot_wider(names_from = name, values_from = value) %>%
      arrange(sex, year) %>%
      select(1:2) %>%
      bind_cols(
        apply(par.mat, 1, function(i){exp(i[7] * smoothed.pca[1:96,1] + i[8] * smoothed.pca[1:96,2])}) %>%
          t %>%
          `colnames<-`(0:95)) %>%
      pivot_longer(cols = -(1:2), names_to = "age") %>%
      rename(Senescence = value)) %>%
  mutate(All = Child + Hump + Senescence,
         age = as.numeric(age))

component.df.all %>%
  filter(year %in% seq(min(.$year), max(.$year), length = 5)) %>%
  pivot_longer(cols = Child:All, names_to = "component") %>%
  mutate(value = value * 1000) %>%
  filter(value > 1e-8,
         age <= 95) %>%
  {ggplot(.) + 
      geom_line(aes(x = age, y = value, col = sex), lwd = 1.2) +
      # scale_color_manual(values = c("female" = f.col, "male" = m.col), name = "") +
      theme_bw() + 
      # ggplot.theme +
      theme(strip.text.y = element_text(size = 25),
            legend.text = element_text(size = 25)) +
      ggtitle(.$name) + ylab(bquote(m[x]~"*"~1000)) + xlab("Age") +
      scale_y_continuous(trans = "log", labels = function(i){as.character(round(i, digits = 3))}, expand = c(0.25, 0.25)) +
      scale_x_continuous(breaks = seq(0, 90, by = 15)) + 
      facet_grid(component ~ year , scales = "free_y", switch = "y")}



component.df.all %>%
  filter(sex == "male") %>%
  {plot_ly(y = 0:95, x = .$year %>% unique %>% sort, z = matrix(.$Hump, nrow = 96), scene = "scene") %>% 
      add_surface(colorbar = list(title = "", x = 0.5, y = 0.9, orientation = "h"),
                  showscale = FALSE) %>%
      layout(
        scene = list(yaxis = list(title = " Age ", tickfont = list(size = 15), titlefont = list(size = 25)),
                     xaxis = list(title = " Year ", tickfont = list(size = 15), titlefont = list(size = 25)),
                     zaxis = list(title = "", range = c(0, range.hump[2])),
                     camera = list(eye = list(x = 1.5, y = -2.3, z = 1.7))))}
#plot mort schedules vs census deaths####
load(paste0("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_RW_Gumbel_1_and_5_common_AR2_phi_hiv_var_separate_normal_ARIMA_hump_base1_common_logpop/RW2/tight_prior/", country, " base1 common logpop RW2 child posterior.RData"))

thiele.no.census.mx.df <- bind_rows(
  thiele.f.loghump.oag.RW.ori$mode$mx_mat_f %>%
    `colnames<-`(bf.idx1$periods) %>%  
    `rownames<-`(thiele_age ) %>%
    reshape2::melt() %>%
    mutate(sex = "female"),
  
  thiele.f.loghump.oag.RW.ori$mode$mx_mat_m%>%
    `colnames<-`(bf.idx1$periods) %>%  
    `rownames<-`(thiele_age ) %>%
    reshape2::melt() %>%
    mutate(sex = "male")
) %>%
  rename(age = Var1,
         year = Var2)

census.mx.df <- inner_join(
  bind_rows(ddharm_bf_census_deaths_f %>% mutate(sex = "female", .after = 1), ddharm_bf_census_deaths_m %>% mutate(sex = "male", .after = 1)) %>%
    pivot_longer(cols = -(1:2), values_to = "death", names_to = "year"),
  
  bind_rows(ddharm_bf_census_f_oag %>% mutate(sex = "female", .after = 1), ddharm_bf_census_m_oag %>% mutate(sex = "male", .after = 1)) %>%
    pivot_longer(cols = -(1:2), values_to = "pop", names_to = "year")
) %>%
  mutate(mx = death / (pop + 0.5 * death),
         year = as.numeric(year),
         year = year - 0.5,
         age = age - 0.5)

thiele.mx.df <- full_join(
  bind_rows(
    full_join(
      thiele.f %>% filter(period %in% unique(census.mx.df$year - 0.5)) %>% mutate(period = period + 0.5) %>% rename(value1 = value),
      thiele.f %>% filter(period %in% unique(census.mx.df$year + 0.5)) %>% mutate(period = period - 0.5) %>% rename(value2 = value)
    ),
    
    full_join(
      thiele.m %>% filter(period %in% unique(census.mx.df$year - 0.5)) %>% mutate(period = period + 0.5) %>% rename(value1 = value),
      thiele.m %>% filter(period %in% unique(census.mx.df$year + 0.5)) %>% mutate(period = period - 0.5) %>% rename(value2 = value)
    )
  ) %>%
    mutate(value = 0.5 * value1 + 0.5 * value2),
  
  
  bind_rows(
    full_join(
      mx.var.df.f %>% filter(period %in% unique(census.mx.df$year - 0.5)) %>% mutate(period = period + 0.5) %>% rename(`2.5%1` = `2.5%`, `97.5%1` = `97.5%`),
      mx.var.df.f %>% filter(period %in% unique(census.mx.df$year + 0.5)) %>% mutate(period = period - 0.5) %>% rename(`2.5%2` = `2.5%`, `97.5%2` = `97.5%`)
    ),
    
    full_join(
      mx.var.df.m %>% filter(period %in% unique(census.mx.df$year - 0.5)) %>% mutate(period = period + 0.5) %>% rename(`2.5%1` = `2.5%`, `97.5%1` = `97.5%`),
      mx.var.df.m %>% filter(period %in% unique(census.mx.df$year + 0.5)) %>% mutate(period = period - 0.5) %>% rename(`2.5%2` = `2.5%`, `97.5%2` = `97.5%`)
    )
  ) %>%
    mutate(`2.5%` = 0.5 * `2.5%1` + 0.5 * `2.5%2`,
           `97.5%` = 0.5 * `97.5%1` + 0.5 * `97.5%2`)
) %>%
  rename(year = period) %>%
  select(age, year, model, sex, hump, value, `2.5%`, `97.5%`) %>%
  bind_rows(  
    full_join(
      GBD.mx %>% filter(name == country, year %in% unique(census.mx.df$year - 0.5)) %>% mutate(year = year + 0.5) %>% rename(value1 = GBD),
      GBD.mx %>% filter(name == country, year %in% unique(census.mx.df$year + 0.5)) %>% mutate(year = year - 0.5) %>% rename(value2 = GBD)
    ) %>%
      mutate(value = 0.5 * value1 + 0.5 * value2,
             model = "GBD",
             sex = replace(sex, sex == "m", "male"),
             sex = replace(sex, sex == "f", "female"),
             hump = "log-Normal hump") %>%
      select(-name, -value1, -value2),
    
    full_join(
      WPP.mx.expand %>% filter(name == country, year %in% unique(census.mx.df$year - 0.5)) %>% mutate(year = year + 0.5) %>% rename(value1 = WPP),
      WPP.mx.expand %>% filter(name == country, year %in% unique(census.mx.df$year + 0.5)) %>% mutate(year = year - 0.5) %>% rename(value2 = WPP)
    ) %>%
      mutate(value = 0.5 * value1 + 0.5 * value2,
             model = "WPP",
             sex = replace(sex, sex == "m", "male"),
             sex = replace(sex, sex == "f", "female"),
             hump = "log-Normal hump") %>%
      select(-name, -value1, -value2),
    
    full_join(
      thiele.no.census.mx.df %>% filter(year %in% unique(census.mx.df$year - 0.5)) %>% mutate(year = year + 0.5) %>% rename(value1 = value),
      thiele.no.census.mx.df %>% filter(year %in% unique(census.mx.df$year + 0.5)) %>% mutate(year = year - 0.5) %>% rename(value2 = value)
    ) %>%
      mutate(value = 0.5 * value1 + 0.5 * value2,
             model = "CCMPP no census deaths",
             sex = replace(sex, sex == "m", "male"),
             sex = replace(sex, sex == "f", "female"),
             hump = "log-Normal hump") %>%
      select(-value1, -value2)
  )

ggplot(thiele.mx.df %>% filter(model == "Thiele RW", hump=="log-Normal hump", sex == "female")) + geom_line(aes(x = age, y = value, col = model), lwd=1.2, linetype = 1) +
  geom_line(data = filter(thiele.mx.df, model != "Thiele RW", sex == "female"), aes(x = age, y = value, col = model), lwd=1.2, linetype = 2) +
  geom_point(data = filter(census.mx.df, mx!=0, sex == "female"), aes(x = age, y = mx), size=2, stroke=1) + 
  geom_ribbon(aes(x = age, ymin = `2.5%`, ymax = `97.5%`), alpha=0.2, fill = 4) +
  scale_shape_manual(values = 0:14) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  ggtitle(paste(country, "Females (Normal hump)")) + ylab(bquote(""*m[x])) + xlab("Age") +
  theme_bw() +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
        legend.position="bottom",
        legend.box = "vertical",
        axis.text = element_text(size=10)) +
  facet_wrap(~ year) + 
  guides(shape = guide_legend(nrow=1, label.position = "top", override.aes = list(size=5)))

ggplot(thiele.mx.df %>% filter(model == "Thiele RW", hump=="log-Normal hump", sex == "male")) + geom_line(aes(x = age, y = value, col = model), lwd=1.2, linetype = 1) +
  geom_line(data = filter(thiele.mx.df, model != "Thiele RW", sex == "male"), aes(x = age, y = value, col = model), lwd=1.2, linetype = 2) +
  geom_point(data = filter(census.mx.df, mx!=0, sex == "male"), aes(x = age, y = mx), size=2, stroke=1) + 
  geom_ribbon(aes(x = age, ymin = `2.5%`, ymax = `97.5%`), alpha=0.2, fill = 4) +
  scale_shape_manual(values = 0:14) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  ggtitle(paste(country, "Males (Normal hump)")) + ylab(bquote(""*m[x])) + xlab("Age") +
  theme_bw() +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
        legend.position="bottom",
        legend.box = "vertical",
        axis.text = element_text(size=10)) +
  facet_wrap(~ year) + 
  guides(shape = guide_legend(nrow=1, label.position = "top", override.aes = list(size=5)))

#5 year census
census.mx.5.df <- inner_join(
  bind_rows(ddharm_bf_census_deaths_f_5 %>% mutate(sex = "female", .after = 1), ddharm_bf_census_deaths_m_5 %>% mutate(sex = "male", .after = 1)) %>%
    pivot_longer(cols = -(1:2), values_to = "death", names_to = "year") %>%
    mutate(age = age + 2),
  
  bind_rows(ddharm_bf_census_f_oag_5 %>% mutate(sex = "female", .after = 1,
                                                age = age%/% 5) %>%
              group_by(sex, age) %>%
              summarise_all(sum) %>%
              mutate(age = age * 5 + 2),
            
            ddharm_bf_census_m_oag_5 %>% mutate(sex = "male", .after = 1,
                                                age = age%/% 5) %>%
              group_by(sex, age) %>%
              summarise_all(sum) %>%
              mutate(age = age * 5 + 2)
  ) %>%
    pivot_longer(cols = -(1:2), values_to = "pop", names_to = "year")
) %>%
  mutate(mx = death / (pop + 0.5 * death),
         year = as.numeric(year),
         year = year - 0.5,
         age = age - 0.5)


thiele.mx.5.df <- full_join(
  bind_rows(
    full_join(
      thiele.f %>% filter(period %in% unique(census.mx.5.df$year - 0.5)) %>% mutate(period = period + 0.5) %>% rename(value1 = value),
      thiele.f %>% filter(period %in% unique(census.mx.5.df$year + 0.5)) %>% mutate(period = period - 0.5) %>% rename(value2 = value)
    ),
    
    full_join(
      thiele.m %>% filter(period %in% unique(census.mx.5.df$year - 0.5)) %>% mutate(period = period + 0.5) %>% rename(value1 = value),
      thiele.m %>% filter(period %in% unique(census.mx.5.df$year + 0.5)) %>% mutate(period = period - 0.5) %>% rename(value2 = value)
    )
  ) %>%
    mutate(value = 0.5 * value1 + 0.5 * value2),
  
  
  bind_rows(
    full_join(
      mx.var.df.f %>% filter(period %in% unique(census.mx.5.df$year - 0.5)) %>% mutate(period = period + 0.5) %>% rename(`2.5%1` = `2.5%`, `97.5%1` = `97.5%`),
      mx.var.df.f %>% filter(period %in% unique(census.mx.5.df$year + 0.5)) %>% mutate(period = period - 0.5) %>% rename(`2.5%2` = `2.5%`, `97.5%2` = `97.5%`)
    ),
    
    full_join(
      mx.var.df.m %>% filter(period %in% unique(census.mx.5.df$year - 0.5)) %>% mutate(period = period + 0.5) %>% rename(`2.5%1` = `2.5%`, `97.5%1` = `97.5%`),
      mx.var.df.m %>% filter(period %in% unique(census.mx.5.df$year + 0.5)) %>% mutate(period = period - 0.5) %>% rename(`2.5%2` = `2.5%`, `97.5%2` = `97.5%`)
    )
  ) %>%
    mutate(`2.5%` = 0.5 * `2.5%1` + 0.5 * `2.5%2`,
           `97.5%` = 0.5 * `97.5%1` + 0.5 * `97.5%2`)
) %>%
  rename(year = period) %>%
  select(age, year, model, sex, hump, value, `2.5%`, `97.5%`) %>%
  bind_rows(  
    full_join(
      GBD.mx %>% filter(name == country, year %in% unique(census.mx.5.df$year - 0.5)) %>% mutate(year = year + 0.5) %>% rename(value1 = GBD),
      GBD.mx %>% filter(name == country, year %in% unique(census.mx.5.df$year + 0.5)) %>% mutate(year = year - 0.5) %>% rename(value2 = GBD)
    ) %>%
      mutate(value = 0.5 * value1 + 0.5 * value2,
             model = "GBD",
             sex = replace(sex, sex == "m", "male"),
             sex = replace(sex, sex == "f", "female"),
             hump = "log-Normal hump") %>%
      select(-name, -value1, -value2),
    
    full_join(
      WPP.mx.expand %>% filter(name == country, year %in% unique(census.mx.5.df$year - 0.5)) %>% mutate(year = year + 0.5) %>% rename(value1 = WPP),
      WPP.mx.expand %>% filter(name == country, year %in% unique(census.mx.5.df$year + 0.5)) %>% mutate(year = year - 0.5) %>% rename(value2 = WPP)
    ) %>%
      mutate(value = 0.5 * value1 + 0.5 * value2,
             model = "WPP",
             sex = replace(sex, sex == "m", "male"),
             sex = replace(sex, sex == "f", "female"),
             hump = "log-Normal hump") %>%
      select(-name, -value1, -value2),
    
    full_join(
      thiele.no.census.mx.df %>% filter(year %in% unique(census.mx.5.df$year - 0.5)) %>% mutate(year = year + 0.5) %>% rename(value1 = value),
      thiele.no.census.mx.df %>% filter(year %in% unique(census.mx.5.df$year + 0.5)) %>% mutate(year = year - 0.5) %>% rename(value2 = value)
    ) %>%
      mutate(value = 0.5 * value1 + 0.5 * value2,
             model = "CCMPP no census deaths",
             sex = replace(sex, sex == "m", "male"),
             sex = replace(sex, sex == "f", "female"),
             hump = "log-Normal hump") %>%
      select(-value1, -value2)
  )

ggplot(thiele.mx.5.df %>% filter(model == "Thiele RW", hump=="log-Normal hump", sex == "female")) + geom_line(aes(x = age, y = value, col = model), lwd=1.2, linetype = 1) +
  geom_line(data = filter(thiele.mx.5.df, model != "Thiele RW", sex == "female"), aes(x = age, y = value, col = model), lwd=1.2, linetype = 2) +
  geom_point(data = filter(census.mx.5.df, mx!=0, sex == "female"), aes(x = age, y = mx), size=2, stroke=1) + 
  geom_ribbon(aes(x = age, ymin = `2.5%`, ymax = `97.5%`), alpha=0.2, fill = 4) +
  scale_shape_manual(values = 0:14) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  ggtitle(paste(country, "Females (Normal hump)")) + ylab(bquote(""*m[x])) + xlab("Age") +
  theme_bw() +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
        legend.position="bottom",
        legend.box = "vertical",
        axis.text = element_text(size=10)) +
  facet_wrap(~ year) + 
  guides(shape = guide_legend(nrow=1, label.position = "top", override.aes = list(size=5)))

ggplot(thiele.mx.5.df %>% filter(model == "Thiele RW", hump=="log-Normal hump", sex == "male")) + geom_line(aes(x = age, y = value, col = model), lwd=1.2, linetype = 1) +
  geom_line(data = filter(thiele.mx.5.df, model != "Thiele RW", sex == "male"), aes(x = age, y = value, col = model), lwd=1.2, linetype = 2) +
  geom_point(data = filter(census.mx.5.df, mx!=0, sex == "male"), aes(x = age, y = mx), size=2, stroke=1) + 
  geom_ribbon(aes(x = age, ymin = `2.5%`, ymax = `97.5%`), alpha=0.2, fill = 4) +
  scale_shape_manual(values = 0:14) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  ggtitle(paste(country, "Males (Normal hump)")) + ylab(bquote(""*m[x])) + xlab("Age") +
  theme_bw() +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
        legend.position="bottom",
        legend.box = "vertical",
        axis.text = element_text(size=10)) +
  facet_wrap(~ year) + 
  guides(shape = guide_legend(nrow=1, label.position = "top", override.aes = list(size=5)))

#plot thiele decomposed####
mx.mat.func <- function(x, hump="normal"){
  age <- 0.5:110.5
  if(hump=="normal"){
    child.mort <- sapply(x$mode$psi_f, function(y){exp(-y*(age-2))}) %*% diag(x$mode$phi_f)
    hump.mort <- sapply(x$mode$epsilon_f, function(y){(y-age)^2}) %*% diag(-x$mode$delta_f) %>% exp() %*% diag(x$mode$lambda_f)
    old.mort <- sapply(x$mode$B_f, function(y){exp(y*(age-92))}) %*% diag(x$mode$A_f)
    child.mort.m <- sapply(x$mode$psi_m, function(y){exp(-y*(age-2))}) %*% diag(x$mode$phi_m)
    hump.mort.m <- sapply(x$mode$epsilon_m, function(y){(y-age)^2}) %*% diag(-x$mode$delta_m) %>% exp() %*% diag(x$mode$lambda_m)
    old.mort.m <- sapply(x$mode$B_m, function(y){exp(y*(age-92))}) %*% diag(x$mode$A_m)
  } else if(hump=="log") {
    child.mort <- sapply(x$mode$psi_f, function(y){exp(-y*(age-2))}) %*% diag(x$mode$phi_f)
    hump.mort <- sapply(x$mode$epsilon_f, function(y){(log(y)-log(age))^2}) %*% diag(-x$mode$delta_f) %>% exp() %*% diag(x$mode$lambda_f)
    old.mort <- sapply(x$mode$B_f, function(y){exp(y*(age-92))}) %*% diag(x$mode$A_f)
    child.mort.m <- sapply(x$mode$psi_m, function(y){exp(-y*(age-2))}) %*% diag(x$mode$phi_m)
    hump.mort.m <- sapply(x$mode$epsilon_m, function(y){(log(y)-log(age))^2}) %*% diag(-x$mode$delta_m) %>% exp() %*% diag(x$mode$lambda_m)
    old.mort.m <- sapply(x$mode$B_m, function(y){exp(y*(age-92))}) %*% diag(x$mode$A_m)
  }
  
  lapply(list(Child = child.mort, Hump = hump.mort, Senescent = old.mort, Total = child.mort + hump.mort + old.mort), function(i){
    i %>%
      `colnames<-`(bf.idx1$periods) %>%  
      `rownames<-`(0:110) %>%
      reshape2::melt()
  })%>% 
    map2(names(.), ~ add_column(.x, stage = rep(.y, nrow(.x)))) %>% 
    bind_rows() %>%
    mutate(sex="female") %>%
    bind_rows(
      lapply(list(Child = child.mort.m, Hump = hump.mort.m, Senescent = old.mort.m, Total = child.mort.m + hump.mort.m + old.mort.m), function(i){
        i %>%
          `colnames<-`(bf.idx1$periods) %>%  
          `rownames<-`(0:110) %>%
          reshape2::melt()
      })%>% 
        map2(names(.), ~ add_column(.x, stage = rep(.y, nrow(.x)))) %>% 
        bind_rows() %>%
        mutate(sex="male")
    )
}

thiele.decomp.loghump <- lapply(loghump.models.list, mx.mat.func, hump="normal") %>%
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>%
  bind_rows(
    lapply(
      list(Child = exp(hmd.mu.f[-1])[1] * exp(-exp(hmd.mu.f[-1])[2] * (0:110 - 2)),
           Hump = init_lambda_f * exp(- init_delta_f * ((0:110 - init_epsilon_f)^2)),
           Senescent = exp(hmd.mu.f[-1])[6] * exp(exp(hmd.mu.f[-1])[7] * (0.:110-92)),
           Total = exp(hmd.mu.f[-1])[1] * exp(-exp(hmd.mu.f[-1])[2] * (0:110 - 2)) + 
             init_lambda_f * exp(- init_delta_f * ((0:110 - init_epsilon_f)^2)) +
             exp(hmd.mu.f[-1])[6] * exp(exp(hmd.mu.f[-1])[7] * (0:110 - 92))
      ), function(i) {
        i %>%
          reshape2::melt() %>%
          mutate(Var1 = 0:110)
      }) %>% 
      map2(names(.), ~ add_column(.x, stage = rep(.y, nrow(.x)))) %>% 
      bind_rows() %>%
      mutate(sex="female", model = "Initial Values") %>%
      crossing(Var2 = bf.idx1$periods),
    
    lapply(
      list(Child = exp(hmd.mu.m[-1])[1] * exp(-exp(hmd.mu.m[-1])[2] * (0:110 - 2)),
           Hump = init_lambda_m * exp(- init_delta_m * ((0:110 - init_epsilon_m)^2)),
           Senescent = exp(hmd.mu.m[-1])[6] * exp(exp(hmd.mu.m[-1])[7] * (0.:110-92)),
           Total = exp(hmd.mu.m[-1])[1] * exp(-exp(hmd.mu.m[-1])[2] * (0:110 - 2)) + 
             init_lambda_m * exp(- init_delta_m * ((0:110 - init_epsilon_m)^2)) +
             exp(hmd.mu.m[-1])[6] * exp(exp(hmd.mu.m[-1])[7] * (0:110 - 92))
      ), function(i) {
        i %>%
          reshape2::melt() %>%
          mutate(Var1 = 0:110)
      }) %>% 
      map2(names(.), ~ add_column(.x, stage = rep(.y, nrow(.x)))) %>% 
      bind_rows() %>%
      mutate(sex="male", model = "Initial Values") %>%
      crossing(Var2 = bf.idx1$periods)
  ) %>%
  setNames(c("age", "period", "value", "stage", "sex", "model")) %>%
  mutate(period = as.factor(period),
         age = as.numeric(age),
         stage = fct_relevel(stage, c("Total", "Child", "Hump", "Senescent")),
         model = fct_relevel(model, c("Initial Values", "Thiele RW"))
  )

ggplot(filter(thiele.decomp.loghump, sex=="female")) + geom_line(aes(x = age, y = value, col = period), lwd = 1.2) +
  scale_color_manual(values = colorRampPalette(colors = c("blue", "red"))(bf.idx1$n_periods)) + 
  scale_y_continuous(trans="log", labels = function(x){as.character(round(x, 5))}) +
  theme_bw() +
  theme(text = element_text(size=20),
        strip.text = element_text(size=20)) + 
  ggtitle(paste(country, "Females")) +
  facet_grid(stage~model, scales="free_y", switch = "y")

ggplot(filter(thiele.decomp.loghump, sex=="male")) + geom_line(aes(x = age, y = value, col = period), lwd = 1.2) +
  scale_color_manual(values = colorRampPalette(colors = c("blue", "red"))(bf.idx1$n_periods)) + 
  theme_bw() +
  scale_y_continuous(trans="log", labels = function(x){as.character(round(x, 5))}) +
  theme(text = element_text(size=20),
        strip.text = element_text(size=20)) + 
  ggtitle(paste(country, "Males")) +
  facet_grid(stage~model, scales="free_y", switch = "y")

#plot 45q15 estimates####
q4515.func <- function(x){
  as_tibble(apply(x$mode$mx_mat_f[16:59,],2,function(x){1-prod((1-0.5*x)/(1+0.5*x))})) %>%
    #as_tibble(apply(x$mode$mx_mat_f[4:12,],2,function(x){1-exp(sum(-5*x))})) %>%
    mutate(sex="female", year = bf.idx1$periods) %>%
    bind_rows(
      as_tibble(apply(x$mode$mx_mat_m[16:59,],2,function(x){1-prod((1-0.5*x)/(1+0.5*x))})) %>%
        #as_tibble(apply(x$mode$mx_mat_m[4:12,],2,function(x){1-exp(sum(-5*x))})) %>%
        mutate(sex="male", year = bf.idx1$periods)
    )
}
wpp.bf.q4515 <- filter(wpp.q4515, Name==country)

gbd.bf.q4515 <- filter(gbd.q4515, location_name == country, year_id %in% bf.idx1$periods, sex_name!="both") %>%
  select(c(4,7,12))

q4515.df.loghump <-  lapply(loghump.models.list, q4515.func) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>%
  bind_rows(
    reshape2::melt(wpp.bf.q4515[,-c(1,3,4)]) %>% 
      mutate(model = "WPP Estimates",
             year = as.numeric(gsub("X|\\.\\d{4}","",variable)),
             value = value/1000,
             sex = Sex) %>%
      select(-c(variable, Sex)),
    
    gbd.bf.q4515 %>% 
      mutate(model = "GBD Estimates",
             year = year_id-2,
             value = val,
             sex = sex_name) %>%
      select(model:sex),
    
    tibble(value = 1-prod((1-0.5 * exp(hmd.mu.f[-1])[1] * exp(-exp(hmd.mu.f[-1])[2] * (15:59-2)) +
                             init_lambda_f * exp(- init_delta_f * ((15:59 - init_epsilon_f)^2)) +
                             exp(hmd.mu.f[-1])[6] * exp(exp(hmd.mu.f[-1])[7] * (15:59 -92))) / 
                            (1+0.5 * exp(hmd.mu.f[-1])[1] * exp(-exp(hmd.mu.f[-1])[2] * (15:59-2)) + 
                               init_lambda_f * exp(- init_delta_f * ((15:59 - init_epsilon_f)^2)) +
                               exp(hmd.mu.f[-1])[6] * exp(exp(hmd.mu.f[-1])[7] * (15:59 -92)))), 
           sex="female", model = "Initial Values", year = bf.idx1$periods),
    
    tibble(value = 1-prod((1-0.5 * exp(hmd.mu.m[-1])[1] * exp(-exp(hmd.mu.m[-1])[2] * (15:59-2)) +
                             init_lambda_m * exp(- init_delta_m * ((15:59 - init_epsilon_m)^2)) +
                             exp(hmd.mu.m[-1])[6] * exp(exp(hmd.mu.m[-1])[7] * (15:59 -92))) / 
                            (1+0.5 * exp(hmd.mu.m[-1])[1] * exp(-exp(hmd.mu.m[-1])[2] * (15:59-2)) + 
                               init_lambda_m * exp(- init_delta_m * ((15:59 - init_epsilon_m)^2)) +
                               exp(hmd.mu.m[-1])[6] * exp(exp(hmd.mu.m[-1])[7] * (15:59 -92)))), 
           sex="male", model = "Initial Values", year = bf.idx1$periods)
  ) %>%
  mutate(model = fct_relevel(model, c("GBD Estimates","WPP Estimates", "Thiele RW")),
         hump = "log-Normal hump")

q4515.var.df <- bind_rows(
  as_tibble(t(apply(sapply(thiele.var.sim, function(i){apply(i$mx_mat_f[16:59,],2,function(x){1-prod((1-0.5*x)/(1+0.5*x))})}), 1, quantile, c(0.025, 0.975))))%>% 
    mutate(year = bf.idx1$periods, sex = "female"),
  as_tibble(t(apply(sapply(thiele.var.sim, function(i){apply(i$mx_mat_m[16:59,],2,function(x){1-prod((1-0.5*x)/(1+0.5*x))})}), 1, quantile, c(0.025, 0.975))))%>% 
    mutate(year = bf.idx1$periods, sex = "male")
)

a <- q4515.df.loghump %>% 
  ggplot() + geom_line(data = filter(q4515.df.loghump, model!="Initial Values"), aes(x = year, y = value, col = sex, linetype = model), lwd = 1.2) + ylab(bquote(""[45]*q[15])) +
  geom_ribbon(data = q4515.var.df, aes(x = year, ymin = `2.5%`, ymax = `97.5%`, fill = sex), alpha=0.2) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_linetype_manual(values=c("dashed", "dotted", "solid", "dotdash")) +
  #geom_point(data = filter(q4515.df.loghump, model=="Initial Values"), aes(x = year, y = value, col = sex), size=3) + 
  ggtitle(bquote(.(country)~"Estimated"[45]*q[15])) +
  theme_bw() +
  theme(text = element_text(size=25), legend.key.size = unit(1.5,"cm"), strip.text = element_text(size=30),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
        legend.text = element_text(size=25))

#plot 5q0 estimates####
q50.func <- function(x){
  as_tibble(apply(x$mode$mx_mat_f[1:5,],2,function(x){1-prod((1-0.5*x)/(1+0.5*x))})) %>%
    mutate(sex="female", year = bf.idx1$periods) %>%
    bind_rows(
      as_tibble(apply(x$mode$mx_mat_m[1:5,],2,function(x){1-prod((1-0.5*x)/(1+0.5*x))})) %>%
        mutate(sex="male", year = bf.idx1$periods)
    )
}

# q50.func <- function(x){
#   as_tibble(apply(x$mode$mx_mat_f[1:5,],2,function(x){1 - exp(-sum(x))})) %>%
#     mutate(sex="female", year = bf.idx1$periods) %>%
#     bind_rows(
#       as_tibble(apply(x$mode$mx_mat_m[1:5,],2,function(x){1 - exp(-sum(x))})) %>%
#         mutate(sex="male", year = bf.idx1$periods)
#     )
# }

q50.df.loghump <- lapply(loghump.models.list, q50.func) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>%
  bind_rows(
    as_tibble(exp(igme.h.mean.f)) %>% mutate(model="IGME Estimates", year=bf.idx1$periods, variable = NULL, sex = "female"),
    as_tibble(exp(igme.h.mean.m)) %>% mutate(model="IGME Estimates", year=bf.idx1$periods, variable = NULL, sex = "male"),
    
    as_tibble(exp(h.mean.f)) %>% mutate(model="WPP Estimates", year=bf.idx1$periods, variable = NULL, sex = "female"),
    as_tibble(exp(h.mean.m)) %>% mutate(model="WPP Estimates", year=bf.idx1$periods, variable = NULL, sex = "male"),
    
    tibble(value = 1-prod((1-0.5 * exp(hmd.mu.f[-1])[1] * exp(-exp(hmd.mu.f[-1])[2] * (0:4-2)) +
                             init_lambda_f * exp(- init_delta_f * ((log(0:4) - log(init_epsilon_f))^2)) +
                             exp(hmd.mu.f[-1])[6] * exp(exp(hmd.mu.f[-1])[7] * (0:4 -92))) / 
                            (1+0.5 * exp(hmd.mu.f[-1])[1] * exp(-exp(hmd.mu.f[-1])[2] * (0:4-2)) + 
                               init_lambda_f * exp(- init_delta_f * ((log(0:4) - log(init_epsilon_f))^2)) +
                               exp(hmd.mu.f[-1])[6] * exp(exp(hmd.mu.f[-1])[7] * (0:4 -92)))), 
           sex="female", model = "Initial Values", year = bf.idx1$periods),
    
    tibble(value = 1-prod((1-0.5 * exp(hmd.mu.m[-1])[1] * exp(-exp(hmd.mu.m[-1])[2] * (0:4-2)) +
                             init_lambda_m * exp(- init_delta_m * ((log(0:4) - log(init_epsilon_m))^2)) +
                             exp(hmd.mu.m[-1])[6] * exp(exp(hmd.mu.m[-1])[7] * (0:4 -92))) / 
                            (1+0.5 * exp(hmd.mu.m[-1])[1] * exp(-exp(hmd.mu.m[-1])[2] * (0:4-2)) + 
                               init_lambda_m * exp(- init_delta_m * ((log(0:4) - log(init_epsilon_m))^2)) +
                               exp(hmd.mu.m[-1])[6] * exp(exp(hmd.mu.m[-1])[7] * (0:4 -92)))), 
           sex="male", model = "Initial Values", year = bf.idx1$periods)
  ) %>%
  mutate(model = fct_relevel(model, "IGME Estimates", "WPP Estimates", "Thiele RW"),
         hump = "log-Normal hump")

q50.var.df <- bind_rows(
  as_tibble(t(apply(sapply(thiele.var.sim, function(i){apply(i$mx_mat_f[1:5,],2,function(x){1-prod((1-0.5*x)/(1+0.5*x))})}), 1, quantile, c(0.025, 0.975))))%>% 
    mutate(year = bf.idx1$periods, sex = "female"),
  as_tibble(t(apply(sapply(thiele.var.sim, function(i){apply(i$mx_mat_m[1:5,],2,function(x){1-prod((1-0.5*x)/(1+0.5*x))})}), 1, quantile, c(0.025, 0.975))))%>% 
    mutate(year = bf.idx1$periods, sex = "male")
)

# q50.var.df <- bind_rows(
#   as_tibble(t(apply(sapply(thiele.var.sim, function(i){apply(i$mx_mat_f[1:5,],2,function(x){1 - exp(-sum(x))})}), 1, quantile, c(0.025, 0.975))))%>% 
#     mutate(year = bf.idx1$periods, sex = "female"),
#   as_tibble(t(apply(sapply(thiele.var.sim, function(i){apply(i$mx_mat_m[1:5,],2,function(x){1 - exp(-sum(x))})}), 1, quantile, c(0.025, 0.975))))%>% 
#     mutate(year = bf.idx1$periods, sex = "male")
# )

b <- q50.df.loghump %>% mutate(hump = fct_inorder(hump)) %>%
  ggplot() + geom_line(data=filter(q50.df.loghump, model!="Initial Values"), aes(x = year, y = value, col = sex, linetype = model), lwd = 1.2) + ylab(bquote(""[5]*q[0])) +
  #geom_point(data = filter(q50.df.loghump, model=="Initial Values"), aes(x = year, y = value, col = sex), size=3) + 
  scale_linetype_manual(values=c("dashed", "dotted", "solid", "dotdash")) +
  geom_ribbon(data = q50.var.df, aes(x = year, ymin = `2.5%`, ymax = `97.5%`, fill = sex), alpha=0.2) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  ggtitle(bquote(.(country)~"Estimated"[5]*q[0])) +
  theme_bw() +
  theme(text = element_text(size=25), legend.key.size = unit(1.5,"cm"), strip.text = element_text(size=30),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
        legend.text = element_text(size=25))

grid.arrange(a,b)

# ggsave(filename="BF 45q15.png",device="png",width=15.4686,height=1 * 15.4686,unit="cm",scale=3.15,
#               plot={grid.arrange(a,b)
#               })


#plot population counts#####
mf5 <- projection_model_frames(bf.idx1)
get.pop<- function(x){
  pop.mat <- matrix(x$mode$population_f, bf.idx1$n_ages, bf.idx1$n_periods+1)
  pop.mat.m <- matrix(x$mode$population_m, bf.idx1$n_ages, bf.idx1$n_periods+1)
  rownames(pop.mat) <-  rownames(pop.mat.m) <- c(bf.idx1$ages)
  colnames(pop.mat) <- colnames(pop.mat.m) <- bf.idx1$periods_out
  bind_rows(
    reshape2::melt(pop.mat) %>% mutate(sex="female"),
    reshape2::melt(pop.mat.m) %>% mutate(sex="male")
  ) %>%
    select(age = Var1, year = Var2, value = value, sex)
}

pop.df <- lapply(loghump.models.list, get.pop) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>%
  bind_rows(
    bind_rows(ddharm_bf_census_f_oag %>% mutate(age = c(ddharm_bf_census_f_oag$age), sex = "female"),
              ddharm_bf_census_m_oag %>% mutate(age = c(ddharm_bf_census_m_oag$age), sex = "male")
    ) %>%
      pivot_longer(!age & !sex) %>% mutate(year=as.numeric(name), name=NULL, model="UNPD Census smoothed"),
    
    bind_rows(mutate(pop.f.oag[,-(1:3)], age = pop.f.oag$age, sex = "female"),
              mutate(pop.m.oag[,-(1:3)], age = pop.m.oag$age, sex = "male")
    )%>%
      pivot_longer(!age & !sex) %>% mutate(year=as.numeric(name), name=NULL, model="WPP Estimates")
  ) %>%
  mutate(cohort = year - age,
         model = fct_relevel(model, c("UNPD Census smoothed", "WPP Estimates", "Thiele RW")),
         hump = "log-Normal hump"
  )

pop.var.df <- bind_rows(
  as_tibble(t(apply(sapply(thiele.var.sim, function(i){i$population_f}), 1, quantile, c(0.025, 0.5, 0.975))))%>%
    mutate(sex = "female",
           age = rep(bf.idx1$ages, bf.idx1$n_periods + 1),
           year = rep(bf.idx1$periods_out, each = bf.idx1$n_ages),
           period5 = sprintf("%d-%d", year, year + 4)),
  
  as_tibble(t(apply(sapply(thiele.var.sim, function(i){i$population_m}), 1, quantile, c(0.025, 0.5, 0.975))))%>%
    mutate(sex = "male",
           age = rep(bf.idx1$ages, bf.idx1$n_periods + 1),
           year = rep(bf.idx1$periods_out, each = bf.idx1$n_ages),
           period5 = sprintf("%d-%d", year, year + 4))
)

# get.census.pop<- function(x){
#   pop.mat <- matrix(x$mode$census_proj_mat_f, bf.idx1$n_ages, ncol(data.f))
#   pop.mat.m <- matrix(x$mode$census_proj_mat_m, bf.idx1$n_ages, ncol(data.f))
#   rownames(pop.mat) <-  rownames(pop.mat.m) <- c(bf.idx1$ages)
#   colnames(pop.mat) <- colnames(pop.mat.m) <- colnames(data.f)
#   bind_rows(
#     reshape2::melt(pop.mat) %>% mutate(sex="female"),
#     reshape2::melt(pop.mat.m) %>% mutate(sex="male")
#   ) %>%
#     select(age = Var1, year = Var2, value = value, sex)
# }

#censpop.df <- lapply(loghump.models.list, get.census.pop) %>% 
# map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
# bind_rows() %>%
# bind_rows(
#   bind_rows(ddharm_bf_census_f_oag %>% mutate(age = c(ddharm_bf_census_f_oag$age), sex = "female"),
#             ddharm_bf_census_m_oag %>% mutate(age = c(ddharm_bf_census_m_oag$age), sex = "male")
#   ) %>%
#     pivot_longer(!age & !sex) %>% mutate(year=as.numeric(name), name=NULL, model="UNPD Census smoothed"),
#   
#   bind_rows(mutate(pop.f.oag[,-(1:3)], age = pop.f.oag$age, sex = "female"),
#             mutate(pop.m.oag[,-(1:3)], age = pop.m.oag$age, sex = "male")
#   )%>%
#     pivot_longer(!age & !sex) %>% mutate(year=as.numeric(name), name=NULL, model="WPP Estimates")
# ) %>%
# mutate(cohort = year - age,
#        model = fct_relevel(model, c("UNPD Census smoothed", "WPP Estimates", "Thiele RW")),
#        hump = "log-Normal hump"
# )
# 
# censpop.var.df <- bind_rows(
#   as_tibble(t(apply(sapply(thiele.var.sim, function(i){i$census_proj_mat_f}), 1, quantile, c(0.025, 0.975))))%>%
#     mutate(sex = "female",
#            age = rep(bf.idx1$ages, ncol(data.f)),
#            year = rep(as.numeric(colnames(data.f)), each = bf.idx1$n_ages),
#            period5 = sprintf("%d-%d", year, year + 4)),
#   
#   as_tibble(t(apply(sapply(thiele.var.sim, function(i){i$census_proj_mat_m}), 1, quantile, c(0.025, 0.975))))%>%
#     mutate(sex = "male",
#            age = rep(bf.idx1$ages, ncol(data.m)),
#            year = rep(as.numeric(colnames(data.m)), each = bf.idx1$n_ages),
#            period5 = sprintf("%d-%d", year, year + 4))
# )
# 
# 
# get.census.pop.5<- function(x){
#   pop.mat <- matrix(x$mode$census_proj_mat_f, bf.idx1$n_ages, ncol(data.f.5))
#   pop.mat.m <- matrix(x$mode$census_proj_mat_m, bf.idx1$n_ages, ncol(data.f.5))
#   rownames(pop.mat) <-  rownames(pop.mat.m) <- c(bf.idx1$ages)
#   colnames(pop.mat) <- colnames(pop.mat.m) <- colnames(data.f)
#   bind_rows(
#     reshape2::melt(pop.mat) %>% mutate(sex="female"),
#     reshape2::melt(pop.mat.m) %>% mutate(sex="male")
#   ) %>%
#     select(age = Var1, year = Var2, value = value, sex) %>%
#     mutate(age5 = age %/% 5) %>%
#     group_by(age5, year, sex) %>%
#     summarise_at(vars(value), sum)
# }
# 
# censpop.df <- lapply(loghump.models.list, get.census.pop) %>% 
#   map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
#   bind_rows() %>%
#   bind_rows(
#     bind_rows(ddharm_bf_census_f_oag %>% mutate(age = c(ddharm_bf_census_f_oag$age), sex = "female"),
#               ddharm_bf_census_m_oag %>% mutate(age = c(ddharm_bf_census_m_oag$age), sex = "male")
#     ) %>%
#       pivot_longer(!age & !sex) %>% mutate(year=as.numeric(name), name=NULL, model="UNPD Census smoothed"),
#     
#     bind_rows(mutate(pop.f.oag[,-(1:3)], age = pop.f.oag$age, sex = "female"),
#               mutate(pop.m.oag[,-(1:3)], age = pop.m.oag$age, sex = "male")
#     )%>%
#       pivot_longer(!age & !sex) %>% mutate(year=as.numeric(name), name=NULL, model="WPP Estimates")
#   ) %>%
#   mutate(cohort = year - age,
#          model = fct_relevel(model, c("UNPD Census smoothed", "WPP Estimates", "Thiele RW")),
#          hump = "log-Normal hump"
#   )
# 
# censpop.var.df <- bind_rows(
#   as_tibble(t(apply(sapply(thiele.var.sim, function(i){i$census_proj_mat_f}), 1, quantile, c(0.025, 0.975))))%>%
#     mutate(sex = "female",
#            age = rep(bf.idx1$ages, ncol(data.f)),
#            year = rep(as.numeric(colnames(data.f)), each = bf.idx1$n_ages),
#            period5 = sprintf("%d-%d", year, year + 4)),
#   
#   as_tibble(t(apply(sapply(thiele.var.sim, function(i){i$census_proj_mat_m}), 1, quantile, c(0.025, 0.975))))%>%
#     mutate(sex = "male",
#            age = rep(bf.idx1$ages, ncol(data.m)),
#            year = rep(as.numeric(colnames(data.m)), each = bf.idx1$n_ages),
#            period5 = sprintf("%d-%d", year, year + 4))
# )
# 
# 
# allpop.df <- full_join(pop.df, censpop.df) %>% mutate(model = fct_relevel(model, c("UNPD Census smoothed", "WPP Estimates", "Thiele RW")))
# allpop.var.df <- full_join(pop.var.df, censpop.var.df)

pop.df %>% filter(year %in% c(1960, colnames(data.f)), !model %in% c("UNPD Census smoothed", "WPP Estimates")) %>%
  ggplot() + geom_line(aes(x = age, y = value, col = sex, linetype=model), lwd = 1) +
  geom_point(data = pop.df %>% filter(year == 1960, model %in% c("WPP Estimates")),
             aes(x = age, y = value, pch = model, col = sex), size=1.5) +
  geom_point(data = pop.df %>% filter(year %in% c(1960, colnames(data.f)), model %in% c("UNPD Census smoothed")),
             aes(x = age, y = value, pch = model, col = sex), size=1.5) +
  scale_shape_manual(values = c(16,0,3)) +
  geom_ribbon(data=pop.var.df %>% filter(year %in% c(1960, colnames(data.f))), 
              aes(x = age, ymin = `2.5%`, ymax = `97.5%`, fill = sex), alpha = 0.2) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_linetype(guide = "none") +
  theme_bw() + 
  theme(text = element_text(size=25)) + ylab("Population counts") + 
  scale_y_continuous(label = function(x) format(x, scientific = TRUE)) +
  ggtitle(country) +
  facet_grid(year~sex, scale="free_y", switch = "y")
#facet_grid_paginate(year~sex, scale="free_y", switch = "y", nrow = 3, ncol = 2, page=1)

pop.df %>% filter(!model %in% c("UNPD Census smoothed", "WPP Estimates"), age %in% seq(0, 80, by = 5)) %>%
  ggplot() + geom_line(aes(x = year, y = value, col = sex, linetype = model), lwd = 1.2) +
  geom_point(data = pop.df %>% filter(year %in% c(1960, colnames(data.f)), model %in% c("UNPD Census smoothed", "WPP Estimates"), age %in% seq(0, 80, by = 5)),
             aes(x = year, y = value, pch = model, col = sex), size=3) +
  scale_shape_manual(values = c(0, 16, 3)) +
  geom_ribbon(data= pop.var.df %>% filter(age %in% seq(0, 80, by = 5)), 
              aes(x = year, ymin = `2.5%`, ymax = `97.5%`, fill = sex), alpha = 0.2) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_linetype(guide = "none") +
  theme_bw() + 
  theme(text = element_text(size=25)) +
  scale_y_continuous(label = function(x) format(x, scientific = TRUE)) +
  ggtitle(country) +
  facet_grid_paginate(age~sex, scale="free", switch="y", nrow = 8, ncol = 2, page=1)

pop.df %>% filter(!model %in% c("UNPD Census smoothed", "WPP Estimates"), age %in% seq(0, 80, by = 5)) %>%
  ggplot() + geom_line(aes(x = year, y = value, col = sex, linetype = model), lwd = 1.2) +
  geom_point(data = pop.df %>% filter(year %in% c(1960, colnames(data.f)), model %in% c("UNPD Census smoothed", "WPP Estimates"), age %in% seq(0, 80, by = 5)),
             aes(x = year, y = value, pch = model, col = sex), size=3) +
  scale_shape_manual(values = c(0, 16, 3)) +
  geom_ribbon(data= pop.var.df %>% filter(age %in% seq(0, 80, by = 5)), 
              aes(x = year, ymin = `2.5%`, ymax = `97.5%`, fill = sex), alpha = 0.2) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  theme(text = element_text(size=25)) +
  scale_linetype(guide = "none") +
  theme_bw() + 
  theme(text = element_text(size=25)) +
  scale_y_continuous(label = function(x) format(x, scientific = TRUE)) +
  ggtitle(country) +
  facet_grid_paginate(age~sex, scale="free", switch="y", nrow = 8, ncol = 2, page=2)

#plot migration####
get.mig<- function(x){
  pop.mat <- matrix(x$mode$migrations_f, bf.idx1$n_ages, bf.idx1$n_periods)
  pop.mat.m <- matrix(x$mode$migrations_m, bf.idx1$n_ages, bf.idx1$n_periods)
  rownames(pop.mat) <- rownames(pop.mat.m) <- c(bf.idx1$ages)
  colnames(pop.mat) <- colnames(pop.mat.m) <- bf.idx1$periods
  bind_rows(reshape2::melt(pop.mat) %>% mutate(sex = "female"),
            reshape2::melt(pop.mat.m) %>% mutate(sex = "male")
  ) %>% select(age = Var1, year = Var2, mig = value, sex)
}

mig.df <- lapply(loghump.models.list, get.mig) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>%
  mutate(hump = "log-Normal hump")

mig.var.df <- bind_rows(
  as_tibble(t(apply(sapply(thiele.var.sim, function(i){i$migrations_f/i$population_f[1:(bf.idx1$n_ages*bf.idx1$n_periods)]}), 1, quantile, c(0.025, 0.975)))) %>%
    mutate(sex = "female",
           age = rep(bf.idx1$ages, bf.idx1$n_periods),
           year = rep(bf.idx1$periods, each = bf.idx1$n_ages),
           period5 = sprintf("%d-%d", year, year + 4)),
  
  as_tibble(t(apply(sapply(thiele.var.sim, function(i){i$migrations_m/i$population_m[1:(bf.idx1$n_ages*bf.idx1$n_periods)]}), 1, quantile, c(0.025, 0.975)))) %>%
    mutate(sex = "male",
           age = rep(bf.idx1$ages, bf.idx1$n_periods),
           year = rep(bf.idx1$periods, each = bf.idx1$n_ages),
           period5 = sprintf("%d-%d", year, year + 4))
) %>%
  mutate(cohort = year - age)

inner_join(mig.df, pop.df) %>% mutate(gx = mig/value, model = fct_relevel(model, "Thiele RW")) %>%
  filter(year %in% c(1960, 1975, 1985, 1995, 2005, 2015)) %>%
  ggplot() + geom_line(aes(x = age, y = gx, col = sex, linetype=model), lwd = 1.2) +
  geom_ribbon(data = mig.var.df %>% filter(year %in% c(1960, 1975, 1985, 1995, 2005, 2015)),
              aes(x = age, ymin = `2.5%`, ymax = `97.5%`, fill = sex), alpha=0.2) +
  ggtitle(paste(country,"Estimated migration Proportions")) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_linetype(guide = "none") +
  theme_bw() + 
  theme(text = element_text(size=25), legend.key.size = unit(1.5,"cm"), strip.text = element_text(size=30),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35)) +
  facet_wrap_paginate(~year, nrow=3, ncol=2, page=1)

inner_join(mig.df, pop.df) %>% mutate(gx = mig/value, model = fct_relevel(model, "Thiele RW")) %>% filter(age %in% seq(0,70, by = 5)) %>%
  ggplot() + geom_line(aes(x = year, y = gx, col = sex, linetype=model), lwd = 1.2) +
  geom_ribbon(data = mig.var.df %>% filter(age %in% seq(0,70, by = 5)),
              aes(x = year, ymin = `2.5%`, ymax = `97.5%`, fill = sex), alpha=0.2) +
  ggtitle(paste(country,"Estimated migration Proportions")) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_linetype(guide = "none") +
  theme_bw() + 
  theme(text = element_text(size=20), legend.key.size = unit(1.5,"cm"), strip.text = element_text(size=30),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35)) +
  facet_wrap_paginate(~age, nrow=3, ncol=3, page=1)

inner_join(mig.df, pop.df) %>% mutate(gx = mig/value, model = fct_relevel(model, "Thiele RW")) %>% filter(age %in% seq(0,70, by = 5)) %>%
  ggplot() + geom_line(aes(x = year, y = gx, col = sex, linetype=model), lwd = 1.2) +
  geom_ribbon(data = mig.var.df %>% filter(age %in% seq(0,70, by = 5)),
              aes(x = year, ymin = `2.5%`, ymax = `97.5%`, fill = sex), alpha=0.2) +
  ggtitle(paste(country,"Estimated migration Proportions")) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_linetype(guide = "none") +
  theme_bw() + 
  theme(text = element_text(size=25), legend.key.size = unit(1.5,"cm"), strip.text = element_text(size=30),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35)) +
  facet_wrap_paginate(~age, nrow=3, ncol=3, page=2)

cohort.label <- paste("Cohort born in", unique(pop.df$cohort))
names(cohort.label) <- unique(pop.df$cohort)

inner_join(mig.df, pop.df) %>% mutate(gx = mig/value, model = fct_relevel(model, "Thiele RW")) %>% filter(cohort %in% seq(1945, 1970, by = 5)) %>%
  ggplot() + geom_line(aes(x = age, y = gx, col = sex, linetype = model), lwd = 1.2) +
  geom_ribbon(data = mig.var.df %>% filter(cohort %in% seq(1945, 1970, by = 5)),
              aes(x = age, ymin = `2.5%`, ymax = `97.5%`, fill = sex), alpha=0.2) +
  ggtitle(paste(country,"Estimated migration Proportions")) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_linetype(guide = "none") +
  theme_bw() + 
  theme(text = element_text(size=25), legend.key.size = unit(1.5,"cm"), strip.text = element_text(size=30),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35)) +
  facet_wrap_paginate(~cohort, nrow=3, ncol=2, page = 1, labeller = labeller(cohort = cohort.label))

#total migration####
total.mig.scale <- format_format(big.mark = " ", decimal.mark = ",", scientific = FALSE)

format(3.0, total.mig.scale)

mig.df %>% group_by(sex, year, model, hump) %>%
  summarise_at(vars(mig), sum) %>%
  ggplot() + geom_line(aes(x = year, y = mig, col = sex), lwd = 2) + 
  scale_y_continuous(labels = total.mig.scale, n.breaks = 10) +
  theme_bw() + 
  theme(text = element_text(size = 25),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ylab("Value") +
  xlab("Year") +
  ggtitle(paste(country, "Total Net Migration"))

#plot tfr####
get.fert<- function(x){
  pop.mat <- matrix(x$mode$fx, bf.idx1$n_fx, bf.idx1$n_periods)
  rownames(pop.mat) <- c(bf.idx1$fertility_ages)
  colnames(pop.mat) <- bf.idx1$periods
  reshape2::melt(pop.mat) %>% select(age = Var1, year = Var2, value = value)
}

fx.df <- lapply(loghump.models.list, get.fert) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>%
  bind_rows(
    reshape2::melt(fx_init.singleyear %>% `rownames<-`(bf.idx1$fertility_ages)) %>%
      select(age = Var1, year = Var2, value = value) %>%
      mutate(model = "Initial Values")
  ) %>%
  mutate(hump = "log-Normal hump",
         model = fct_relevel(model, c("Initial Values", "Thiele RW"))) 

fx.var.df <- as_tibble(t(apply(sapply(thiele.var.sim, function(i){i$fx}), 1, quantile, c(0.025, 0.975)))) %>%
  mutate(age = rep(bf.idx1$fertility_ages, bf.idx1$n_periods),
         year = rep(bf.idx1$periods, each = bf.idx1$n_fx),
         period5 = sprintf("%d-%d", year, year + 4))


tfr.df <- fx.df %>%
  group_by(year, model, hump) %>%
  summarise_at(vars(value), sum)

tfr.var.df <- as_tibble(t(apply(sapply(thiele.var.sim, function(i){apply(matrix(i$fx, bf.idx1$n_fx, bf.idx1$n_periods), 2, function(j){sum(j)})}), 1, quantile, c(0.025, 0.975)))) %>%
  mutate(year = bf.idx1$periods)

tfr.df %>% filter(model != "Initial Values") %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) +
  geom_point(data = filter(tfr.df, model=="Initial Values"), aes(x = year, y = value), size = 3) +
  geom_ribbon(data = tfr.var.df, aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha=0.2, fill = 2) +
  ggtitle(paste(country, "Estimated Total Fertility Rates")) +
  scale_color_discrete(guide = "none") +
  theme_bw() + 
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
        legend.key.width = unit(2,"cm")) 

aa <- tfr.df %>% filter(model != "Initial Values") %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) +
  geom_point(data = filter(tfr.df, model=="Initial Values"), aes(x = year, y = value), size = 3) +
  geom_ribbon(data = tfr.var.df, aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha=0.2, fill = 2) +
  ggtitle(paste(country, "Estimated Total Fertility Rates")) +
  theme_bw() +  
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
        legend.key.width = unit(2,"cm")) +
  scale_color_discrete(guide = "none")


#plot mean age of birth####
mfr.df <- fx.df %>%
  group_by(year, model, hump) %>%
  summarise_at(vars(value), function(i){sum(i *15:49)/ sum(i)})

mfr.var.df <- as_tibble(t(apply(sapply(thiele.var.sim, function(i){apply(matrix(i$fx, bf.idx1$n_fx, bf.idx1$n_periods), 2, function(j){sum(j*15:49)/sum(j)})}), 1, quantile, c(0.025, 0.975)))) %>%
  mutate(year = bf.idx1$periods)


mfr.df %>% filter(model != "Initial Values") %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) +
  geom_point(data = filter(mfr.df, model=="Initial Values"), aes(x = year, y = value), size = 3) +
  geom_ribbon(data = mfr.var.df, aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha=0.2, fill = 2) +
  ggtitle(paste(country, "Estimated Mean Fertility Age")) +
  scale_color_discrete(guide = "none") +
  theme_bw() +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
        legend.key.width = unit(2,"cm")) 

#fertility####
fx.df %>% filter(model != "Initial Values", year %in% c(1960, 1975, 1985, 1995, 2005, 2015)) %>%
  ggplot() + geom_line(aes(x = age, y = value, col = model), lwd = 1.2) +
  geom_point(data = filter(fx.df, model=="Initial Values", year %in% c(1960, 1975, 1985, 1995, 2005, 2015)), aes(x = age, y = value), size = 3) +
  geom_ribbon(data = filter(fx.var.df, year %in% c(1960, 1975, 1985, 1995, 2005, 2015)), aes(x = age, ymin = `2.5%`, ymax = `97.5%`), alpha=0.2, fill = 2) +
  ggtitle(paste(country, "Estimated Fertility Rates")) +
  scale_color_discrete(guide = "none") +
  theme_bw() + 
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
        legend.key.width = unit(2,"cm")) +
  facet_wrap_paginate(~year, nrow = 3, ncol = 2, page=1)

fx.df %>% filter(model != "Initial Values", age %in% c(15,25,35,45)) %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) +
  geom_point(data = filter(fx.df, model=="Initial Values", age %in% c(15,25,35,45)), aes(x = year, y = value), size = 3) +
  geom_ribbon(data = filter(fx.var.df, age %in% c(15,25,35,45)), aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha=0.2, fill = 2) +
  scale_color_discrete(guide = "none") +
  theme_bw() + 
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
        legend.key.width = unit(2,"cm")) +
  facet_wrap_paginate(~age, scale="free_y", nrow = 2, ncol = 2, page=1)


bb <- fx.df %>% filter(model != "Initial Values", age %in% c(15,25,35,45)) %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) +
  geom_point(data = filter(fx.df, model=="Initial Values", age %in% c(15,25,35,45)), aes(x = year, y = value), size = 3) +
  geom_ribbon(data = filter(fx.var.df, age %in% c(15,25,35,45)), aes(x = year, ymin = `2.5%`, ymax = `97.5%`), alpha=0.2, fill = 2) +
  scale_color_discrete(guide = "none") +
  theme_bw() + 
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
        legend.key.width = unit(2,"cm")) +
  facet_wrap_paginate(~age, scale="free_y", nrow = 2, ncol = 2, page=1)

# ggsave(filename="BF fert.png",device="png",width=15.4686,height=1 * 15.4686,unit="cm",scale=3.15,
#               plot={grid.arrange(aa, bb)})

#raw pop plot####
ddharm_bf_census_f %>% pivot_longer(cols = 3:ncol(ddharm_bf_census_f), names_to = "year") %>%
  ggplot() + geom_point(aes(x = AgeStart, y = value), size = 2) +
  geom_point(data = subset(ddharm_bf_census_f %>% pivot_longer(cols = 3:ncol(ddharm_bf_census_f), names_to = "year"), AgeStart%%5==0), 
             aes(x = AgeStart, y = value), col="red", size = 2) + 
  scale_x_continuous(breaks = seq(0, 100, by = 10), name = "Age") +
  ggtitle(paste(params$country,"Raw Population Counts (Females)")) + 
  theme_bw() +
  theme(strip.text = element_text(size=20),
        plot.title = element_text(size=30, hjust=0.5),
        text = element_text(size = 30)) +
  facet_grid(~ year)

ddharm_bf_census_f_oag %>% pivot_longer(cols = 2:ncol(ddharm_bf_census_f_oag), names_to = "year") %>%
  ggplot() + geom_point(aes(x = age, y = value), size = 2) +
  # geom_point(data = subset(ddharm_bf_census_f_oag %>% pivot_longer(cols = 2:ncol(ddharm_bf_census_f_oag), names_to = "year"), age%%5==0), 
  #            aes(x = age, y = value), col="red") + 
  scale_x_continuous(breaks = seq(0, 100, by = 10), name = "Age") +
  ggtitle(paste(params$country,"Smoothed Population Counts (Females)")) + 
  theme_bw() +
  theme(strip.text = element_text(size=20),
        plot.title = element_text(size=30, hjust=0.5),
        text = element_text(size = 30)) +
  facet_grid(~ year)


#ggsave plots####

ggsave(filename="BF mort schedule female.png",device="png",width=15.4686,height=1.2 * 15.4686,unit="cm",scale=3.15,
       plot={ggplot(thiele.f %>% filter(age %in% 0:70, hump=="log-Normal hump")) + geom_line(aes(x = age, y = value, col=model), lwd=1.2, linetype = 2) +
           geom_line(data = filter(thiele.f, model == "Spline average"), aes(x = age, y = value, col = model), lwd=1.2, linetype=1) +
           geom_point(data = filter(DHS.plot.f, value!=0), aes(x = age, y = value, shape = tips), alpha = 0) +
           geom_point(data = filter(DHS.plot.f, value!=0, tips %in% 1:3), aes(x = age, y = value, shape = tips), size=2, stroke=1) +
           geom_point(data = filter(DHS.plot.f, value!=0, !tips %in% 1:3), aes(x = age, y = value, shape = tips), size=1, alpha = 0.3) +
           geom_ribbon(data = mx.var.df.f, aes(x = age, ymin = `2.5%`, ymax = `97.5%`), alpha=0.2, fill = 4) +
           scale_shape_manual(values = 0:14) +
           scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
           ggtitle(paste(country, "Females")) + ylab(bquote(""[5]*m[x])) + xlab("Age") +
           theme(text = element_text(size=25),
                 plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
                 axis.text  = element_text(size=10),
                 legend.position="bottom",
                 legend.box = "vertical") +
           facet_wrap(~period) +
           guides(shape = guide_legend(nrow=1, label.position = "top", override.aes = list(size=5)))
       })

ggsave(filename="BF mort schedule male.png",device="png",width=15.4686,height=1.2 * 15.4686,unit="cm",scale=3.15,
       plot={
         ggplot(thiele.m %>% filter(age %in% 0:70, hump == "log-Normal hump")) + geom_line(aes(x = age, y = value, col=model), lwd=1.2, linetype = 2) +
           geom_line(data = filter(thiele.m, model == "Spline average"), aes(x = age, y = value, col = model), lwd=1.2, linetype=1)+
           geom_point(data = filter(DHS.plot.m, value!=0), aes(x = age, y = value, shape = tips), alpha = 0) +
           geom_point(data = filter(DHS.plot.m, value!=0, tips %in% 1:3), aes(x = age, y = value, shape = tips), size=2, stroke=1) +
           geom_point(data = filter(DHS.plot.m, value!=0, !tips %in% 1:3), aes(x = age, y = value, shape = tips), size=1, alpha = 0.3) +
           geom_ribbon(data = mx.var.df.m, aes(x = age, ymin = `2.5%`, ymax = `97.5%`), alpha=0.2, fill = 4) +
           scale_shape_manual(values = 0:14) +
           scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
           facet_wrap(~period) +
           ggtitle(paste(country, "Males")) + ylab(bquote(""[5]*m[x])) + xlab("Age") +
           theme(text = element_text(size=25),
                 plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
                 legend.position="bottom",
                 legend.box = "vertical",
                 axis.text = element_text(size=10)) +
           facet_wrap(~period) +
           guides(shape = guide_legend(nrow=1, label.position = "top", override.aes = list(size=5)))
         
       })


ggsave(filename="BF mort schedule time.png",device="png",width=15.4686,height=1.2 * 15.4686,unit="cm",scale=3.15,
       plot={
         ggplot(bind_rows(thiele.f, thiele.m) %>% filter(age %in% seq(0,70, by = 5), hump=="log-Normal hump")) + geom_line(aes(x = period, y = value, col=sex, linetype=model), lwd=1.2) +
           scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
           geom_ribbon(data = bind_rows(mx.var.df.m, mx.var.df.f) %>% filter(age %in% seq(0,70, by = 5)), aes(x = period, ymin = `2.5%`, ymax = `97.5%`, fill = sex), alpha=0.2) +
           ggtitle(paste(country, "Estimate Mortality Schedules")) + ylab(bquote(""[5]*m[x])) + xlab("5-year age groups") +
           scale_color_manual(values = c("red", "blue")) +
           scale_fill_manual(values = c("red", "blue")) +
           theme(text = element_text(size=25),
                 plot.title = element_text(hjust = 0.5, face = "bold", size = 35)) +
           facet_wrap(~age)
       })


tp.all %>% mutate(country = "Namibia")

a <- tp.all %>%
  ggplot() + stat_boxplot(geom="errorbar", aes(y = tips, x = value, group = tips), position=position_dodge(0), lwd = 0.8) +
  geom_boxplot(aes(y = tips, x = value, group = tips), lwd = 0.8) +
  ggtitle(paste(country, "Estimated TiPS")) +
  theme(text = element_text(size=25), legend.key.size = unit(1.5,"cm"), strip.text = element_text(size=30),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35))

b <- tp.all.bf %>%
  ggplot() + stat_boxplot(geom="errorbar", aes(y = tips, x = value, group = tips), position=position_dodge(0), lwd = 0.8) +
  geom_boxplot(aes(y = tips, x = value, group = tips), lwd = 0.8) +
  ggtitle("Burkina Faso Estimated TiPS") +
  theme(text = element_text(size=25), legend.key.size = unit(1.5,"cm"), strip.text = element_text(size=30),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35))

grid.arrange(a,b, nrow = 1)

thiele.par.loghump %>% mutate(country = "Namibia") %>%
  bind_rows(thiele.par.loghump.bf %>% mutate(country = "Burkina Faso")) -> thiele.par.loghump.all

par.var.df %>% mutate(country = "Namibia") %>%
  bind_rows(par.var.df.bf %>% mutate(country = "Burkina Faso")) -> par.var.df.all

thiele.par.loghump.all <- thiele.par.loghump.all %>% mutate(country = fct_relevel(country, "Namibia"))
par.var.df.all <- par.var.df.all %>% mutate(country = fct_relevel(country, "Namibia"))

ggplot(thiele.par.loghump.all %>% filter(model == "Thiele RW", name %in% c("lambda", "delta", "epsilon"))) + 
  geom_line(aes(x = period5, y = value, linetype = model, col = sex), lwd = 1.2) +
  geom_ribbon(data = par.var.df.all %>% filter(name %in% c("lambda", "delta", "epsilon")), aes(x = period5, ymin = `2.5%`, ymax = `97.5%`,fill=sex),alpha=0.2) + 
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
        legend.key.width = unit(2,"cm"),
        legend.position = "bottom",
        legend.direction="vertical",
        legend.text.align = 0,
        strip.text =  element_text(size=30, margin = margin(l = 10, t = 10, b = 10, r = 10))) +
  xlab("Year") +
  ggtitle(bquote(.(country)~"Estimated Parameters")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash"), guide = "none") + 
  facet_grid(name~country, scales="free_y", labeller = labeller(name=label_parsed), switch="y")

ggplot(thiele.par.loghump.all %>% filter(model == "Thiele RW")) + 
  geom_line(aes(x = period5, y = value, linetype = model, col = sex), lwd = 1.2) +
  geom_ribbon(data = par.var.df.all, aes(x = period5, ymin = `2.5%`, ymax = `97.5%`,fill=sex),alpha=0.2) + 
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35),
        legend.key.width = unit(2,"cm"),
        legend.position = "bottom",
        legend.direction="vertical",
        legend.text.align = 0,
        strip.text =  element_text(size=30, margin = margin(l = 10, t = 10, b = 10, r = 10))) +
  xlab("Year") +
  ggtitle(bquote(.(country)~"Estimated Parameters")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash"), guide = "none") + 
  facet_grid_paginate(name~country, scales="free_y", labeller = labeller(name=label_parsed), switch="y", nrow = 4, ncol = 2, page = 2)


