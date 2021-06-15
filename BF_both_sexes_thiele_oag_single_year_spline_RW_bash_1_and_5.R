args <- commandArgs(T)
print(args) #list the command line arguments. 

myvar <- as.numeric(args[1])

.libPaths(c("C:/Users/ktang3/Documents/R/win-library/4.0", .libPaths()))
#setwd("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_RW_Gumbel_1_and_5")
setwd("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_RW_Gumbel_1_and_5_common_AR2_all")
load("C:/Users/ktang3/Documents/cohort smooth 1900-2017.RData")
skip<-c("Ethiopia","Central African Republic","Comoros","Sao Tome and Principe","Botswana","Cape Verde","Equatorial Guinea","Eritrea","Nigeria (Ondo State)","Ghana","Mauritania","Sudan")
joint.countries<-names(aggr.mat.cohort.0)[!names(aggr.mat.cohort.0)%in%skip]

params <- list(
  country = joint.countries[myvar], 
  age.knot.space = 2.5,
  year.knot.space = 2.5
)

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

require(devtools)
require(httr)

req <- GET("https://api.github.com/repos/sarahertog/ddharmony/git/trees/main?recursive=1")
stop_for_status(req)
filelist <- unlist(lapply(content(req)$tree, "[", "path"), use.names = F)

for (filename in filelist) {
  one_function <- paste0("https://github.com/sarahertog/ddharmony/blob/main/", filename, "?raw=TRUE")
  source_url(one_function)
  rm(one_function)
}
rm(req, filelist, filename)

#dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_bothsexes_thiele_loghump_oag_RW_originalscale_spline_RW_aggr_gumbel_common"))
#dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_bothsexes_thiele_loghump_oag_RW_originalscale_spline_RW_aggr_gumbel_common_AR2"))
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_bothsexes_thiele_loghump_oag_RW_originalscale_spline_RW_aggr_gumbel_common_AR2_all"))

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

make_tmb_obj <- function(data,
                         par,
                         model = "ccmpp_tmb",
                         inner_verbose = FALSE,
                         calc_outputs = TRUE,
                         random,
                         DLL = "leapfrog_TMBExports",
                         map = list()) {
  
  data$model <- model
  data$calc_outputs <- as.integer(calc_outputs)
  
  obj <- TMB::MakeADFun(data = data,
                        parameters = par,
                        DLL = DLL,
                        silent = !inner_verbose,
                        random = random,
                        map = map)
  class(obj) <- "tmb_obj"
  
  obj
}

fit_tmb <- function(tmb_input,
                    outer_verbose = TRUE,
                    inner_verbose = FALSE,
                    max_iter = 250,
                    random,
                    DLL = "leapfrog_TMBExports",
                    map = list(),
                    stepmin = 1,
                    stepmax = 1
) {
  
  ## stopifnot(inherits(tmb_input, "tmb_input"))
  
  if(is.null(tmb_input$model))
    tmb_input$model <- "ccmpp_tmb"
  
  obj <- make_tmb_obj(data = tmb_input$data,
                      par = tmb_input$par_init,
                      model = tmb_input$model,
                      inner_verbose = inner_verbose,
                      calc_outputs = 0L,
                      random = random,
                      DLL = DLL,
                      map = map)
  
  trace <- if(outer_verbose) 1 else 0
  f <- withCallingHandlers(
    stats::nlminb(obj$par, obj$fn, obj$gr,
                  control = list(trace = trace,
                                 iter.max = max_iter,
                                 step.min = stepmin,
                                 step.max = stepmax)),
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
  
  objout <- make_tmb_obj(tmb_input$data,
                         tmb_input$par_init,
                         model = tmb_input$model,
                         inner_verbose = inner_verbose,
                         calc_outputs = 1L,
                         DLL = DLL,
                         map = map,
                         random = random)
  f$mode <- objout$report(f$par.full)
  
  val <- c(f, obj = list(objout))
  class(val) <- "leapfrog_fit"
  
  val
}

igme.5q0.m<-read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/5q0 male.csv")
igme.5q0.f<-read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/5q0 female.csv")
wpp.fx <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/WPP fx.csv")
wpp.pop <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/WPP Pop estimates.csv")
wpp.pop.age.specific <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/WPP pop age specific.csv")
wpp.q4515 <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/WPP 45q15.csv")
gbd.q4515 <- read.csv(file="C:/Users/ktang3/Desktop/Imperial/SSA_mort/GBD 45q15.csv",header=T)
wpp.qx <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/WPP age specific.csv")

gbd.q4515$location_name<-str_replace(gbd.q4515$location_name,"Democratic Republic of the Congo","Congo Democratic Republic")
gbd.q4515$location_name<-str_replace(gbd.q4515$location_name,"Côte d'Ivoire","Cote d'Ivoire")
gbd.q4515$location_name<-str_replace(gbd.q4515$location_name,"United Republic of Tanzania","Tanzania")

igme.5q0.m$Country.Name <- str_replace(igme.5q0.m$Country.Name,"Democratic Republic of the Congo","Congo Democratic Republic")
igme.5q0.f$Country.Name <- str_replace(igme.5q0.f$Country.Name,"Democratic Republic of the Congo","Congo Democratic Republic")
igme.5q0.m$Country.Name<-str_replace(igme.5q0.m$Country.Name,"United Republic of Tanzania","Tanzania")
igme.5q0.f$Country.Name<-str_replace(igme.5q0.f$Country.Name,"United Republic of Tanzania","Tanzania")

wpp.qx$name <- str_replace(wpp.qx$name,"Democratic Republic of the Congo","Congo Democratic Republic")
wpp.qx$name <- str_replace(wpp.qx$name,"Democratic Republic of the Congo","Congo Democratic Republic")
wpp.qx$name<-str_replace(wpp.qx$name,"United Republic of Tanzania","Tanzania")
wpp.qx$name<-str_replace(wpp.qx$name,"United Republic of Tanzania","Tanzania")

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

library(MortCast)
LQcoef.f <- LQcoef %>% filter(sex=="Female", !age%in%c("0","1-4")) %>% select(ax:vx)
LQcoef.m <- LQcoef %>% filter(sex=="Male", !age%in%c("0","1-4")) %>% select(ax:vx)

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

#DDHarmonized smoothed
census_pop_counts <- DDharmonize_validate_PopCounts(locid = ifelse(country=="Cote d'Ivoire", 384, 
                                                                   ifelse(country=="Tanzania",834, country)),       
                                                    times = 1950:2020,
                                                    DataSourceShortName = "DYB") # time frame for censuses to extract from Demographic Yearbook

#####################MIXING DE-FACTO AND DE-JURE HERE
ddharm_bf_census_m <- census_pop_counts %>%  
  filter(AgeSpan %in% c(-1, 1), AgeLabel != "Total", ReferencePeriod > 1960, SexID == 1, complete == TRUE) %>%
  select(ReferencePeriod, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID) %>%
  distinct() %>%
  pivot_wider(names_from = ReferencePeriod, values_from = DataValue) %>%
  group_by(AgeStart, AgeLabel, AgeSpan) %>%
  arrange(AgeStart, StatisticalConceptName) %>% ###De-factor first
  select(-c(StatisticalConceptName, SexID)) %>%
  summarise_all(function(y){first(na.omit(y))}) %>%
  select(-AgeSpan) %>%
  ungroup()

ddharm_bf_census_f <- census_pop_counts %>%  
  filter(AgeSpan %in% c(-1, 1), AgeLabel != "Total", ReferencePeriod > 1960, SexID == 2, complete == TRUE) %>%
  select(ReferencePeriod, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID) %>%
  distinct() %>%
  pivot_wider(names_from = ReferencePeriod, values_from = DataValue) %>%
  group_by(AgeStart, AgeLabel, AgeSpan) %>%
  arrange(AgeStart, StatisticalConceptName) %>% ###De-factor first
  select(-c(StatisticalConceptName, SexID)) %>%
  summarise_all(function(y){first(na.omit(y))}) %>%
  select(-AgeSpan) %>%
  ungroup()

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

##census data 5 year age groups
ddharm_bf_census_m_5 <- census_pop_counts %>%  
  filter(AgeSpan %in% c(-1, 5), AgeLabel != "Total", ReferencePeriod > 1960, SexID == 1, five_year == TRUE) %>%
  select(ReferencePeriod, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID) %>%
  distinct() %>%
  pivot_wider(names_from = ReferencePeriod, values_from = DataValue) %>%
  group_by(AgeStart, AgeLabel, AgeSpan) %>%
  arrange(AgeStart, StatisticalConceptName) %>% ###De-factor first
  select(-c(StatisticalConceptName, SexID)) %>%
  summarise_all(function(y){first(na.omit(y))}) %>%
  select(-AgeSpan) %>%
  ungroup()

ddharm_bf_census_f_5 <- census_pop_counts %>%  
  filter(AgeSpan %in% c(-1, 5), AgeLabel != "Total", ReferencePeriod > 1960, SexID == 2, five_year == TRUE) %>%
  select(ReferencePeriod, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID) %>%
  distinct() %>%
  pivot_wider(names_from = ReferencePeriod, values_from = DataValue) %>%
  group_by(AgeStart, AgeLabel, AgeSpan) %>%
  arrange(AgeStart, StatisticalConceptName) %>% ###De-factor first
  select(-c(StatisticalConceptName, SexID)) %>%
  summarise_all(function(y){first(na.omit(y))}) %>%
  select(-AgeSpan) %>%
  ungroup()

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


##WPP Pop Estimates
pop.singleyear <- wpp.pop.age.specific %>% filter(Name == country) %>% 
  setNames(c(names(wpp.pop)[1:3], 0:100)) %>%
  reshape2::melt(id.vars=c("Name", "Sex", "Reference")) %>% 
  mutate(year = Reference, 
         value = as.numeric(str_replace_all(value, "\\s", "")) * 1000,
         age = variable) %>%
  select(Sex, year, age, value) %>%
  pivot_wider(values_from = value, names_from = year) %>%
  mutate(age = as.numeric(levels(age))[age])

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
  A.age<-cbind(A.age,bspline(bf.idx1$ages,knots.age,j))
}

A.year<-c()
for(j in 1:no.basis.time) {
  A.year<-cbind(A.year,bspline(bf.idx1$periods,knots.time,j))
}

te.spline<-A.year%x%A.age

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

thiele.loghump.prior <- function(h, sex) {
  cat(h,"\n")
  log.m0 <- h - log(1 - 0.5 * exp(h)) + log(0.2)
  if(sex == "male"){
    LQ.mx <- c(log.m0, as.matrix(LQcoef.m[,1:3]) %*% c(1, h, h^2))
  } else if (sex == "female") {
    LQ.mx <- c(log.m0, as.matrix(LQcoef.f[,1:3]) %*% c(1, h, h^2))
  }
  
  #get priors for phi and psi using LQ mx at 0-4 to 10-14
  child.coef <- lm((LQ.mx[2:3]-log.m0) ~ c(5, 10)-1)$coef
  psi <- -child.coef
  
  #get priors for lambda, delta and epsilon using 10-14 to 40-44
  hump.coef <- lm(LQ.mx[3:9] ~  log(seq(12, 42, by = 5)) + I(log(seq(12, 42, by = 5))^2))$coef
  delta <- -hump.coef[3]
  epsilon <- exp(-hump.coef[2] / (2 * hump.coef[3]))
  lambda <- exp(hump.coef[1]+delta*(log(epsilon)^2))
  
  lambda <- ifelse(lambda<1e-5 | lambda>0.015, 8e-3, lambda)
  delta <- ifelse(delta<1e-10, 1e-3, delta)
  epsilon <- ifelse(epsilon>40 | epsilon<10, 25, epsilon)
  
  #get priors for A and B using 65-69 to 90-94
  old.coef <- lm(LQ.mx[14:19] ~ seq(67, 92, by = 5))$coef
  A <- exp(old.coef[1])
  B <- old.coef[2]
  
  thiele.min <- function(par, dat){
    all.age <- seq(2, 92, by = 5)
    par <- exp(par)
    est <- par[1] * exp(-par[2] * (all.age-2)) + par[3] * exp(-par[4]*((log(all.age) - log(par[5]))^2)) + par[6] * exp(par[7] * (all.age-92))
    ess <- mean(sum((log(est[-1]) - dat[-1])^2) + 1e3 * (log(est[1]) - dat[1])^2)
    return(ess)
  }
  nlm <- nlminb(start = log(c(exp(log.m0), psi, 0.02, 2, 22, exp(LQ.mx[19])-2e-3*exp(-0.5*(log(92)-log(22))^2), 0.1)), thiele.min, dat = LQ.mx[1:19], control = list(eval.max = 8000, iter.max = 8000, step.min = 1e-10, step.max = 1e-4))
  stopifnot(nlm$convergence ==0)
  cat(nlm$message,"\n")
  return(setNames(exp(nlm$par), c("phi", "psi", "lambda", "delta", "epsilon", "A", "B")))
}

igme.h.mean.f <- igme.5q0.df %>% filter(Sex=="Female", year %in% bf.idx1$periods) %>% right_join(as_tibble(bf.idx1$periods)%>%mutate(year=value)) %>% .$child.mort %>% log()
h.mean.f <- wpp.5q0.interpolate %>% filter(Sex=="Female", year %in% bf.idx1$periods) %>% .$q50 %>% log() ##Using WPP 5q0 estimates
#thiele.prior.f <- sapply(h.mean.f, thiele.prior, sex="female")
thiele.loghump.prior.f <- sapply(h.mean.f, thiele.loghump.prior, sex="female")

igme.h.mean.m <- igme.5q0.df %>% filter(Sex=="Male", year %in% bf.idx1$periods) %>% right_join(as_tibble(bf.idx1$periods)%>%mutate(year=value)) %>% .$child.mort %>% log()
h.mean.m <- wpp.5q0.interpolate %>% filter(Sex=="Male", year %in% bf.idx1$periods) %>% .$q50 %>% log() ##Using WPP 5q0 estimates
#thiele.prior.m <- sapply(h.mean.m, thiele.prior, sex="male")
thiele.loghump.prior.m <- sapply(h.mean.m, thiele.loghump.prior, sex="male")

thiele_age <- 0.5:97.5

##DHS data cohort smoothed
bf5.smooth <- aggr.mat.cohort.0[[country]] %>%
  filter(period %in% bf.idx1$periods) %>%
  group_by(mm1, tips, DHS, agegr, period) %>%
  summarise_at(vars(pyears, event, pyears2, adjusted), sum) %>%
  ungroup()

bf5.smooth$period <- factor(bf5.smooth$period,levels=bf.idx1$periods)
bf5.smooth$tips <- factor(bf5.smooth$tips,levels=0:14)

dhs.start.age <- 15
dhs.end.age <- 59

bf5.f.no0.smooth <- bf5.smooth %>% filter(mm1=="female", agegr >= dhs.start.age, agegr <= dhs.end.age) %>% arrange(period,tips,agegr)
bf5.m.no0.smooth <- bf5.smooth %>% filter(mm1=="male", agegr >= dhs.start.age, agegr <= dhs.end.age) %>% arrange(period,tips,agegr)

basepop.f <- ifelse(pop.f.oag$`1960`==0, 1, pop.f.oag$'1960')
basepop.m <- ifelse(pop.m.oag$`1960`==0, 1, pop.m.oag$'1960')

data.f <- if(nrow(ddharm_bf_census_f)!=0) {as.matrix(log(ddharm_bf_census_f_oag[,-1]))} else {as.matrix(tibble(.rows=open.age+1))}; data.m <- if(nrow(ddharm_bf_census_m)!=0) {as.matrix(log(ddharm_bf_census_m_oag[,-1]))} else {as.matrix(tibble(.rows=open.age+1))}
data.f.5 <- if(!is.null(data.f)) {as.matrix(log(ddharm_bf_census_f_oag_5[,-1] %>% select(!matches(colnames(data.f)))))} else {as.matrix(log(ddharm_bf_census_f_oag_5[,-1]))}; data.m.5 <- if(!is.null(data.m)){as.matrix(log(ddharm_bf_census_m_oag_5[,-1] %>% select(!matches(colnames(data.f)))))} else {as.matrix(log(ddharm_bf_census_m_oag_5[,-1]))}

if(country == "Zimbabwe"){
  data.f <- as.matrix(log(ddharm_bf_census_f_oag[,-(1:2)])); data.m <- as.matrix(log(ddharm_bf_census_m_oag[,-(1:2)]))
  data.f.5 <- as.matrix(log(ddharm_bf_census_f_oag_5[,-(1:2)] %>% select(!matches(colnames(data.f))))); data.m.5 <- as.matrix(log(ddharm_bf_census_m_oag_5[,-(1:2)] %>% select(!matches(colnames(data.f)))))
  #data.f.55 <- as.matrix(log(ddharm_bf_census_f_oag_5[,-(1:2)])); data.m <- as.matrix(log(ddharm_bf_census_m_oag_5[,-(1:2)]))
}

init_lambda_f <- thiele.loghump.prior.f[3,1]; init_lambda_m <- thiele.loghump.prior.m[3,1]; 
init_delta_f <- thiele.loghump.prior.f[4,1]; init_delta_m <- thiele.loghump.prior.m[4,1]; 
init_epsilon_f <- thiele.loghump.prior.f[5,1]; init_epsilon_m <- thiele.loghump.prior.m[5,1]; 

if(length(unique(bf5.smooth$DHS)) == 1 || diff(as.numeric(str_extract(unique(bf5.smooth$DHS), "\\d+"))) > 14) {
  load("C:/Users/ktang3/Documents/more countries final avg sex Rwanda.RData")
  skip<-c("Ethiopia","Central African Republic","Comoros","Sao Tome and Principe","Botswana","Cape Verde","Equatorial Guinea","Eritrea","Nigeria (Ondo State)","Ghana","Mauritania","Sudan")
  joint.countries<-names(aggr.mat.cohort.0)[!names(aggr.mat.cohort.0)%in%skip]
  
  all.list<-function(x,age.start,age.end,year.start,year.end){
    no.basis = 15
    knots<-seq(0,1,length=no.basis-2)
    dk<-knots[2]-knots[1]	
    knots<-c(knots[1]-dk*(3:1),knots,knots[no.basis-2]+dk*(1:3))
    age<-seq(0,1,length=age.end-age.start+1)
    year<-seq(0,1,length=year.end-year.start+1)
    
    A.age<-c()
    for(j in 1:no.basis) {
      A.age<-cbind(A.age,bspline(age,knots,j))
    }
    
    A.year<-c()
    for(j in 1:no.basis) {
      A.year<-cbind(A.year,bspline(year,knots,j))
    }
    
    te.spline<-A.year%x%A.age
    
    ind<-function(x){
      aggr.mat.reduced<-x
      data.mat.m<-aggr.mat.reduced[aggr.mat.reduced$mm1=="male" & aggr.mat.reduced$agegr %in% age.start:age.end & aggr.mat.reduced$period %in% year.start:year.end,c(7,5,3,2,4)]
      data.mat.f<-aggr.mat.reduced[aggr.mat.reduced$mm1=="female" & aggr.mat.reduced$agegr %in% age.start:age.end & aggr.mat.reduced$period %in% year.start:year.end,c(7,5,3,2,4)]
      
      age.start.col<-min(which(A.age[min(data.mat.m[,3])-age.start+1,]!=0))
      age.end.col<-max(which(A.age[max(data.mat.m[,3])-age.start+1,]!=0))
      year.start.col<-min(which(A.year[min(data.mat.m[,4])-year.start+1,]!=0))
      year.end.col<-max(which(A.year[max(data.mat.m[,4])-year.start+1,]!=0))
      age.start<-min(data.mat.m[,3])
      age.end<-max(data.mat.m[,3])
      year.start<-min(data.mat.m[,4])
      year.end<-max(data.mat.m[,4])
      tp<-max(data.mat.m[,5],data.mat.f[,5])
      
      ind.place<-c(age.start.col,age.end.col,year.start.col,year.end.col,age.start,age.end,year.start,year.end,tp)
      names(ind.place)<-c("age.start.col", "age.end.col","year.start.col","year.end.col","age.start","age.end","year.start","year.end","tp")
      ind.place
    }  
    
    
    everything.func<-function(i){
      age.seq<-as.numeric(ind.list[[i]][["age.start"]]:ind.list[[i]][["age.end"]])
      year.seq<-as.numeric(ind.list[[i]][["year.start"]]:ind.list[[i]][["year.end"]])
      
      age.m<-A.age[age.seq-age.start+1,]%*%joint.countries.age.m[[i]]
      age.f<-A.age[age.seq-age.start+1,]%*%joint.countries.age.f[[i]]
      time.m<-A.year[year.seq-year.start+1,]%*%joint.countries.time.m[[i]]
      time.f<-A.year[year.seq-year.start+1,]%*%joint.countries.time.f[[i]]
      agetime.m<-matrix(te.spline[age.seq-age.start+1+rep((age.end-age.start+1)*(year.seq-year.start),each=length(age.seq)),]%*%joint.countries.2d.m[[i]],length(age.seq),length(year.seq))
      agetime.f<-matrix(te.spline[age.seq-age.start+1+rep((age.end-age.start+1)*(year.seq-year.start),each=length(age.seq)),]%*%joint.countries.2d.f[[i]],length(age.seq),length(year.seq))
      
      mort.m<-joint.countries.avg.m[[i]]+t(apply(apply(agetime.m,2,function(x){x+age.m}),1,function(x){x+time.m}))
      mort.f<-joint.countries.avg.f[[i]]+t(apply(apply(agetime.f,2,function(x){x+age.f}),1,function(x){x+time.f}))
      rownames(age.m)<-rownames(age.f)<-rownames(agetime.m)<-rownames(agetime.f)<-rownames(mort.m)<-rownames(mort.f)<-age.seq
      rownames(time.m)<-rownames(time.f)<-colnames(agetime.m)<-colnames(agetime.f)<-colnames(mort.m)<-colnames(mort.f)<-year.seq
      
      
      if(i=="Rwanda"){
        mort.m[,4:6]<-mort.m[,4:6]+matrix(rep.par.list$rwanda_geno_m,56,3,byrow=T)
        mort.f[,4:6]<-mort.f[,4:6]+matrix(rep.par.list$rwanda_geno_f,56,3,byrow=T)
      }
      
      tips<-joint.countries.tp.c[[i]]+tp.common
      
      list(tp=tips,avg.m=joint.countries.avg.m[[i]],avg.f=joint.countries.avg.f[[i]],age.m=age.m,age.f=age.f,time.m=time.m,time.f=time.f,agetime.m=agetime.m,agetime.f=agetime.f,mort.m=mort.m,mort.f=mort.f)
    }
    
    ind.list<-lapply(aggr.mat.cohort.0[joint.countries],ind)
    rep.par.list<-split(unname(x$env$last.par.best),names(x$env$last.par.best))
    joint.countries.avg.m<-split(rep.par.list$avg_m,joint.countries)
    joint.countries.avg.f<-split(rep.par.list$avg_f,joint.countries)
    joint.countries.age.m<-split(rep.par.list$spline_params_m_age,rep(joint.countries,each=no.basis))
    joint.countries.age.f<-split(rep.par.list$spline_params_f_age,rep(joint.countries,each=no.basis))
    joint.countries.time.m<-split(rep.par.list$spline_params_m_time,rep(joint.countries,each=no.basis))
    joint.countries.time.f<-split(rep.par.list$spline_params_f_time,rep(joint.countries,each=no.basis))
    joint.countries.2d.m<-split(rep.par.list$spline_params_m_2d,rep(joint.countries,each=no.basis*no.basis))
    joint.countries.2d.f<-split(rep.par.list$spline_params_f_2d,rep(joint.countries,each=no.basis*no.basis))
    
    joint.countries.tp.c<-split(rep.par.list$tips_params,rep(joint.countries,each=length(rep.par.list$tips_params_common)))
    tp.common<-rep.par.list$tips_params_common
    tp.common[6]<-tp.common[6]+rep.par.list$tp_common_5
    tp.common[11]<-tp.common[11]+rep.par.list$tp_common_10
    
    
    setNames(lapply(joint.countries,everything.func),joint.countries)
  }
  
  more.countries.avg<-all.list(tmb.full.joint.common,age.start=10,age.end=65,year.start=1990,year.end=2017)
  
  for(i in 1:length(more.countries.avg)){
    more.countries.avg[[i]]$avg.m<-more.countries.avg[[i]]$avg.m+unname(tmb.full.joint.common$env$last.par.best["avg_common_m"])
    more.countries.avg[[i]]$avg.f<-more.countries.avg[[i]]$avg.f+unname(tmb.full.joint.common$env$last.par.best["avg_common_f"])
    more.countries.avg[[i]]$mort.m<-more.countries.avg[[i]]$mort.m+unname(tmb.full.joint.common$env$last.par.best["avg_common_m"])
    more.countries.avg[[i]]$mort.f<-more.countries.avg[[i]]$mort.f+unname(tmb.full.joint.common$env$last.par.best["avg_common_f"])
  }
  
  tp.init <- more.countries.avg[[params$country]]$tp; map <- list(log_lambda_tp = factor(NA),
                                                                  tp_params = factor(rep(NA,15)),
                                                                  tp_slope = factor(NA),
                                                                  tp_params_5 = factor(NA),
                                                                  tp_params_10 = factor(NA))
} else {tp.init <- rep(0,15); map <- list()}


full.penal.gx <- as(0.5 * diag(no.basis.time) %x% crossprod(diff(diag(no.basis.age),differences=1)) +
                      0.5 * crossprod(diff(diag(no.basis.time),differences=1)) %x% diag(no.basis.age) +
                      #0.5 * crossprod(diff(diag(no.basis.time))%x%diff(diag(no.basis.age))) +
                      1e-3 * diag(no.basis.time * no.basis.age), "sparseMatrix")

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
 
# full.penal.gx.AR2 <- as(0.5 * diag(no.basis.time) %x% AR2.PREC(no.basis.age, c(2*0.9, -0.9)) +
#                           0.5 * AR2.PREC(no.basis.time, c(2*0.9, -0.9)) %x% diag(no.basis.age), "sparseMatrix")

full.penal.gx.AR2.tensor <- as(AR2.PREC(no.basis.age, c(2*0.9, -0.9)) %x% AR2.PREC(no.basis.time, c(2*0.9, -0.9)), "sparseMatrix")

full.penal.fx <- as(0.5 * diag(no.basis.time) %x% crossprod(diff(diag(no.basis.fert))) +
                      0.5 * crossprod(diff(diag(no.basis.time))) %x% diag(no.basis.fert) +
                      1e-3 * diag(no.basis.time * no.basis.fert), "sparseMatrix")

# full.penal.fx.sum2zero <- as(0.5 * diag(no.basis.time) %x% crossprod(diff(diag(no.basis.fert))) +
#                                0.5 * crossprod(diff(diag(no.basis.time))) %x% diag(no.basis.fert) +
#                                1e5 * tcrossprod(rep(1,no.basis.time * no.basis.fert)), "sparseMatrix")

full.penal.fx.AR2.tensor <- as(AR2.PREC(no.basis.time, c(2*0.9, -0.9)) %x% AR2.PREC(no.basis.fert, c(2*0.9, -0.9)), "sparseMatrix")

#full.penal.time <- as(crossprod(diff(diag(no.basis.time), differences = 2)) + 1e-3 * diag(no.basis.time),"sparseMatrix")
full.penal.time <- as(crossprod(diff(diag(no.basis.time), differences = 1)) + 1e-3 * diag(no.basis.time),"sparseMatrix")


gumbel.theta.fx <- -log(0.01) * sqrt(mean(diag(te.spline.fert %*% solve(full.penal.fx) %*% t(te.spline.fert)))) * 1.96 / log(1.1)
#gumbel.theta.fx.sum2zero <- -log(0.01) * sqrt(mean(diag(te.spline.fert %*% solve(full.penal.fx.sum2zero) %*% t(te.spline.fert)))) * 1.96 / log(1.1)

gumbel.theta.gx <- -log(0.01) *  sqrt(mean(diag(te.spline %*% solve(full.penal.gx) %*% t(te.spline)))) * 1.96 / 0.08

gumbel.theta.phi <- -log(0.01) * sqrt(mean(A.year %*% solve(full.penal.time) %*% t(A.year))) * 1.96 / log(1.1)
gumbel.theta.psi <- -log(0.01) * sqrt(mean(A.year %*% solve(full.penal.time) %*% t(A.year))) * 1.96 / log(1.1)
gumbel.theta.A <- -log(0.01) * sqrt(mean(A.year %*% solve(full.penal.time) %*% t(A.year))) * 1.96 / log(1.1)
gumbel.theta.B <- -log(0.01) * sqrt(mean(A.year %*% solve(full.penal.time) %*% t(A.year))) * 1.96 / log(1.1)

gumbel.theta.lambda <- -log(0.01) * sqrt(mean(A.year %*% solve(full.penal.time) %*% t(A.year))) * 1.96 / log(3)
gumbel.theta.delta <- -log(0.01) * sqrt(mean(A.year %*% solve(full.penal.time) %*% t(A.year))) * 1.96 / log(3)
gumbel.theta.epsilon <- -log(0.01) * sqrt(mean(A.year %*% solve(full.penal.time) %*% t(A.year))) * 1.96 / log(3)

gumbel.theta.tp <- -log(0.01) * 1.96/1

gumbel.theta.AR2.marginal.fx <- sqrt(mean(diag(te.spline.fert %*% solve(full.penal.fx.AR2.tensor) %*% t(te.spline.fert)))) * 1.96 / log(1.1)
gumbel.theta.AR2.marginal.gx <- sqrt(mean(diag(te.spline %*% solve(full.penal.gx.AR2.tensor) %*% t(te.spline)))) * 1.96 / 0.08 
#gumbel.theta.AR2.marginal.gx <- sqrt(mean(diag(te.spline %*% solve(full.penal.gx.AR2) %*% t(te.spline)))) * 1.96 / 0.08 

gumbel.theta.AR2.marginal.phi <- 
  gumbel.theta.AR2.marginal.psi <-
  gumbel.theta.AR2.marginal.A <-
  gumbel.theta.AR2.marginal.B <- -log(0.01) *  sqrt(mean(diag(A.year %*% solve(AR2.PREC(no.basis.time, c(2*0.9, -0.9))) %*% t(A.year))))  * 1.96 / log(1.1)
# gumbel.theta.AR2.marginal.psi <- -log(0.01) * sqrt(mean(A.year %*% solve(full.penal.time) %*% t(A.year))) * 1.96 / log(1.1)
# gumbel.theta.AR2.marginal.A <- -log(0.01) * sqrt(mean(A.year %*% solve(full.penal.time) %*% t(A.year))) * 1.96 / log(1.1)
# gumbel.theta.AR2.marginal.B <- -log(0.01) * sqrt(mean(A.year %*% solve(full.penal.time) %*% t(A.year))) * 1.96 / log(1.1)

gumbel.theta.AR2.marginal.lambda <- 
  gumbel.theta.AR2.marginal.delta <- 
  gumbel.theta.AR2.marginal.epsilon <- -log(0.01) *  sqrt(mean(diag(A.year %*% solve(AR2.PREC(no.basis.time, c(2*0.9, -0.9))) %*% t(A.year)))) * 1.96 / log(3)

data.loghump.vec.RW <- list(log_basepop_mean_f = log(basepop.f), log_basepop_mean_m = log(basepop.m),
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
                            
                            df = bf5.f.no0.smooth$adjusted, dm = bf5.m.no0.smooth$adjusted,
                            Ef = bf5.f.no0.smooth$pyears2, Em = bf5.m.no0.smooth$pyears2,
                            df_age = bf5.f.no0.smooth$agegr + 1, dm_age = bf5.m.no0.smooth$agegr + 1,
                            df_time = match(bf5.f.no0.smooth$period, levels(bf5.f.no0.smooth$period)), dm_time = match(bf5.m.no0.smooth$period, levels(bf5.m.no0.smooth$period)),
                            df_tp = c(bf5.f.no0.smooth$tips)-1, dm_tp = c(bf5.m.no0.smooth$tips)-1,
                            
                            log_phi_mean_f = log(thiele.loghump.prior.f[1,]), log_phi_mean_m = log(thiele.loghump.prior.m[1,]),
                            log_psi_mean_f = log(thiele.loghump.prior.f[2,]), log_psi_mean_m = log(thiele.loghump.prior.m[2,]),
                            log_lambda_mean_f = log(init_lambda_f), log_lambda_mean_m = log(init_lambda_m),
                            log_delta_mean_f = log(init_delta_f), log_delta_mean_m = log(init_delta_m),
                            log_epsilon_mean_f = log(init_epsilon_f), log_epsilon_mean_m = log(init_epsilon_m),
                            log_A_mean_f = log(thiele.loghump.prior.f[6,]), log_A_mean_m = log(thiele.loghump.prior.m[6,]),
                            log_B_mean_f = log(thiele.loghump.prior.f[7,]), log_B_mean_m = log(thiele.loghump.prior.m[7,]),
                            
                            thiele_age = thiele_age,
                            
                            full_penal_time = full.penal.time,
                            
                            full_penal_gx = full.penal.gx.AR2.tensor,
                            
                            full_penal_fx = full.penal.fx,
                            #full_penal_fx = full.penal.fx.sum2zero,
                            
                            penal_tp = as(crossprod(diff(diag(15))), "sparseMatrix"), #only penalising the latter 14 coefficients 
                            penal_tp_0 = as(diag(c(1, rep(0, 14))), "sparseMatrix"),
                            null_penal_tp = as(exp(15)*tcrossprod(c(0,1,1,1,rep(0, 11))),"sparseMatrix"),
                            tp_mean = -2:12,
                            
                            D_time = as(A.year, "sparseMatrix"),
                            D_agetime = as(te.spline, "sparseMatrix"),
                            D_agetime_fert = as(te.spline.fert, "sparseMatrix"),
                            
                            theta_fx = gumbel.theta.fx,
                            #theta_fx = gumbel.theta.fx.sum2zero,
                            theta_gx = gumbel.theta.gx,
                            theta_phi = gumbel.theta.phi,
                            theta_psi = gumbel.theta.psi,
                            theta_lambda = gumbel.theta.lambda,
                            theta_delta = gumbel.theta.delta,
                            theta_epsilon = gumbel.theta.epsilon,
                            theta_A = gumbel.theta.A,
                            theta_B = gumbel.theta.B,
                            theta_tp = gumbel.theta.tp,
                            
                            theta_marginal_gx = gumbel.theta.AR2.marginal.gx,
                            theta_marginal_fx = gumbel.theta.AR2.marginal.fx,
                            theta_marginal_phi = gumbel.theta.AR2.marginal.phi,
                            theta_marginal_psi = gumbel.theta.AR2.marginal.psi,
                            theta_marginal_lambda = gumbel.theta.AR2.marginal.lambda,
                            theta_marginal_delta = gumbel.theta.AR2.marginal.delta,
                            theta_marginal_epsilon = gumbel.theta.AR2.marginal.epsilon,
                            theta_marginal_A = gumbel.theta.AR2.marginal.A,
                            theta_marginal_B = gumbel.theta.AR2.marginal.B
)

par.vec <- list(log_tau2_logpop_f = c(2,3), log_tau2_logpop_m = c(2,3),
                
                log_basepop_f = log(basepop.f), log_basepop_m = log(basepop.m),
                log_fx_spline_params = rep(0, no.basis.fert * no.basis.time),
                log_fx_intercept = 0,
                gx_f_spline_params = rep(0, no.basis.time * no.basis.age), gx_m_spline_params = rep(0, no.basis.time * no.basis.age),
                
                log_lambda_tp = 5,
                tp_params = tp.init,
                
                tp_slope = 0,
                tp_params_5 = 0,
                tp_params_10 = 0,
                
                log_phi_f_spline_params = rep(0, no.basis.time), log_phi_m_spline_params = rep(0, no.basis.time),
                log_psi_f_spline_params = rep(0, no.basis.time), log_psi_m_spline_params = rep(0, no.basis.time),
                log_lambda_f_spline_params = rep(0, no.basis.time), log_lambda_m_spline_params = rep(0, no.basis.time),
                log_delta_f_spline_params = rep(0, no.basis.time), log_delta_m_spline_params = rep(0, no.basis.time),
                log_epsilon_f_spline_params = rep(0, no.basis.time), log_epsilon_m_spline_params = rep(0, no.basis.time),
                log_A_f_spline_params = rep(0, no.basis.time), log_A_m_spline_params = rep(0, no.basis.time),
                log_B_f_spline_params = rep(0, no.basis.time), log_B_m_spline_params = rep(0, no.basis.time),
                
                log_lambda_fx = log((gumbel.theta.fx/-log(0.01))^2) + 0.1, 
                #log_lambda_fx = log((gumbel.theta.fx.sum2zero/-log(0.01))^2) + 0.1, 
                log_lambda_gx = log((gumbel.theta.gx/-log(0.01))^2) + 0.1,
 
                log_lambda_phi = log((gumbel.theta.phi/-log(0.01))^2) + 0.1,
                log_lambda_psi = log((gumbel.theta.psi/-log(0.01))^2) + 0.1,
                log_lambda_A = log((gumbel.theta.A/-log(0.01))^2) + 0.1,
                log_lambda_B = log((gumbel.theta.B/-log(0.01))^2) + 0.1,
                
                log_lambda_lambda = log((gumbel.theta.lambda/-log(0.01))^2) + 0.1,
                log_lambda_delta = log((gumbel.theta.delta/-log(0.01))^2) + 0.1,
                log_lambda_epsilon = log((gumbel.theta.epsilon/-log(0.01))^2) + 0.1,
                
                log_tau2_logpop = c(log(1.96^2/log(1.5)^2), log(1.96^2/log(2)^2), log(1.96^2/log(1.5)^2), log(1.96^2/log(2)^2)),
                log_dispersion = c(1.3, 1.3),
                
                log_marginal_lambda_gx = log((gumbel.theta.AR2.marginal.gx/-log(0.01))^2) + 0.1,
                log_marginal_lambda_fx = log((gumbel.theta.AR2.marginal.fx/-log(0.01))^2) + 0.1,
                log_marginal_lambda_phi = log((gumbel.theta.AR2.marginal.phi/-log(0.01))^2) + 0.1,
                log_marginal_lambda_psi = log((gumbel.theta.AR2.marginal.psi/-log(0.01))^2) + 0.1,
                log_marginal_lambda_A = log((gumbel.theta.AR2.marginal.A/-log(0.01))^2) + 0.1,
                log_marginal_lambda_B = log((gumbel.theta.AR2.marginal.B/-log(0.01))^2) + 0.1,
                log_marginal_lambda_lambda = log((gumbel.theta.AR2.marginal.lambda/-log(0.01))^2) + 0.1,
                log_marginal_lambda_delta = log((gumbel.theta.AR2.marginal.delta/-log(0.01))^2) + 0.1,
                log_marginal_lambda_epsilon = log((gumbel.theta.AR2.marginal.epsilon/-log(0.01))^2) + 0.1
                )

input.thiele.loghump.oag.vec.RW <- list(data = data.loghump.vec.RW, par_init = par.vec, model = "ccmpp_vr_tmb")

system.time(thiele.f.loghump.oag.RW.ori <- fit_tmb(input.thiele.loghump.oag.vec.RW, inner_verbose=TRUE,
                                                   random = c("log_basepop_f","log_basepop_m",
                                                              "log_fx_spline_params",
                                                              "gx_f_spline_params","gx_m_spline_params",
                                                              "tp_params",
                                                              "log_phi_f_spline_params", "log_phi_m_spline_params",
                                                              "log_psi_f_spline_params", "log_psi_m_spline_params",
                                                              "log_lambda_f_spline_params", "log_lambda_m_spline_params",
                                                              "log_delta_f_spline_params", "log_delta_m_spline_params",
                                                              "log_epsilon_f_spline_params", "log_epsilon_m_spline_params",
                                                              "log_A_f_spline_params", "log_A_m_spline_params",
                                                              "log_B_f_spline_params", "log_B_m_spline_params"),
                                                   DLL="ccmpp_bothsexes_thiele_loghump_oag_RW_originalscale_spline_RW_aggr_gumbel_common_AR2_all",
                                                   map = map,
                                                   stepmin = 1e-10, stepmax = 1)
)

save(thiele.f.loghump.oag.RW.ori,file=paste0(params$country," tau Gumbel AR2.RData"))
