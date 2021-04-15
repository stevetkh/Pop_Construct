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

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_f_thiele.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_f_thiele"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_f_thiele_noDHS.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_f_thiele_noDHS"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_AR_f_no0.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_AR_f_no0"))

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


igme.5q0.f<-read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/5q0 male.csv")
igme.5q0.m<-read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/5q0 female.csv")
wpp.fx <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/WPP fx.csv")
wpp.pop <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/WPP Pop estimates.csv")
wpp.q4515 <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/WPP 45q15.csv")

country <- "Uganda"

library(MortCast)
load("~/cohort smooth 1900-2017.RData")

open.age <- 75
n_ages = open.age / 5 + 1

##IGME priors
igme.5q0.df <- igme.5q0.f %>% filter(Country.Name == country) %>% reshape2::melt() %>% 
  mutate(year = as.numeric(gsub("X|\\.5","",variable)), child.mort = value/1000, Sex = "Female") %>%
  select(year,child.mort, Sex) %>%
  bind_rows(
    igme.5q0.m %>% filter(Country.Name == country) %>% reshape2::melt() %>% 
      mutate(year = as.numeric(gsub("X|\\.5","",variable)), child.mort = value/1000, Sex = "Male") %>%
      select(year,child.mort, Sex)
  )

igme.5q0.5 <- igme.5q0.df %>% mutate(year5 = 5 * floor(year/5)) %>% group_by(Sex, year5) %>%
  summarise_at(vars(child.mort),mean) %>% ungroup()

#DDHarmonized smoothed
census_pop_counts <- DDharmonize_validate_PopCounts(locid = country,       
                                                    times = 1950:2020) # time frame for censuses to extract

#####################MIXING DE-FACTO AND DE-JURE HERE
ddharm_bf_census_m <- census_pop_counts %>% select(ReferencePeriod, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID) %>%
  filter(AgeSpan %in% c(-1, 5), AgeLabel != "Total") %>%
  distinct() %>%
  pivot_wider(names_from = ReferencePeriod, values_from = DataValue) %>%
  filter(SexID == 1) %>%
  group_by(AgeStart, AgeLabel, AgeSpan) %>%
  arrange(AgeStart, StatisticalConceptName) %>% ###De-factor first
  select(-c(StatisticalConceptName, SexID)) %>%
  summarise_all(function(y){first(na.omit(y))}) %>%
  select(-AgeSpan) %>%
  ungroup()

ddharm_bf_census_f <- census_pop_counts %>% select(ReferencePeriod, StatisticalConceptName, AgeStart, AgeLabel, AgeSpan, DataValue, SexID) %>%
  filter(AgeSpan %in% c(-1, 5), AgeLabel != "Total") %>%
  distinct() %>%
  pivot_wider(names_from = ReferencePeriod, values_from = DataValue) %>%
  filter(SexID == 2) %>%
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
  
  smoothed.adult <- getSmoothedPop5(popM = dat.m[[3]], popF = dat.f[[3]], Age = dat.m[[1]], EduYrs = 2, subgroup = "adult")
  smoothed.child <- getSmoothedPop5(popM = dat.m[[3]], popF = dat.f[[3]], Age = dat.m[[1]], EduYrs = 2, subgroup = "child")
  
  adult.grad <- as.numeric(str_extract(smoothed.adult$best_smooth_method,"\\s\\d"))
  child.grad <- as.numeric(str_extract(smoothed.child$best_smooth_method,"\\s\\d"))
  
  cat("adult bestGrad5 =",adult.grad,"; child bestGrad5 =",child.grad, "\n")
  
  ddharm_smoothed <- bind_rows(
    ddharm_smoothed,
    as_tibble(DemoTools::smooth_age_5(dat.m[[3]], dat.m[[1]], method="MAV", n=adult.grad)) %>% mutate(age = dat.m[[1]], AgeLabel = dat.m[[2]], 
                                                                                                      sex="male", year = as.numeric(names(ddharm_bf_census_f)[i])),
    as_tibble(DemoTools::smooth_age_5(dat.f[[3]], dat.f[[1]], method="MAV", n=adult.grad)) %>% mutate(age = dat.f[[1]], AgeLabel = dat.f[[2]],
                                                                                                      sex="female", year = as.numeric(names(ddharm_bf_census_m)[i]))
  )
}

ddharm_smoothed_mat <- ddharm_smoothed %>% pivot_wider(names_from = year, values_from = value) %>% arrange(sex, age)

ddharm_census_f_smoothed <- ddharm_smoothed_mat %>% filter(sex=="female") %>% select(-sex)
ddharm_census_m_smoothed <- ddharm_smoothed_mat %>% filter(sex=="male") %>% select(-sex)

ddharm_bf_census_f_smoothed_aggr <- ddharm_census_f_smoothed %>%
  arrange(AgeLabel) %>%
  mutate(age.aggr = ifelse(AgeLabel %in% c("0","0-4","1-4"),"0-4", 
                           ifelse(AgeLabel %in% c(sprintf("%i-%i", seq(open.age, 90, by=5), seq(open.age+4, 94, by=5)), "80+", "85+", "90+", "95+", "98+"), paste0(open.age,"+"), AgeLabel))) %>%
  group_by(age.aggr) %>%
  summarise_at(vars(-AgeLabel),sum,na.rm=T) %>%
  mutate(age = as.numeric(regmatches(age.aggr, regexpr("\\d+",age.aggr))), .before = 1) %>%
  arrange(age) %>% select(-age.aggr)

ddharm_bf_census_m_smoothed_aggr <- ddharm_census_m_smoothed %>%
  arrange(AgeLabel) %>%
  mutate(age.aggr = ifelse(AgeLabel %in% c("0","0-4","1-4"),"0-4", 
                           ifelse(AgeLabel %in% c(sprintf("%i-%i", seq(open.age, 90, by=5), seq(open.age+4, 94, by=5)), "80+", "85+", "90+", "95+", "98+"), paste0(open.age,"+"), AgeLabel))) %>%
  group_by(age.aggr) %>%
  summarise_at(vars(-AgeLabel),sum,na.rm=T) %>%
  mutate(age = as.numeric(regmatches(age.aggr, regexpr("\\d+",age.aggr))), .before = 1) %>%
  arrange(age) %>% select(-age.aggr)

ddharm_bf_census_m_raw <- ddharm_bf_census_m %>%
  mutate(age.aggr = ifelse(AgeLabel %in% c("0","0-4","1-4"),"0-4", 
                           ifelse(AgeLabel %in% c(sprintf("%i-%i", seq(open.age, 90, by=5), seq(open.age+4, 94, by=5)), "80+", "85+", "90+", "95+", "98+"), paste0(open.age,"+"), AgeLabel))) %>%
  group_by(age.aggr) %>%
  summarise_at(vars(-(AgeStart:AgeLabel)),sum,na.rm=T) %>%
  mutate(age = as.numeric(regmatches(age.aggr, regexpr("\\d+",age.aggr))), .before = 1) %>%
  arrange(age) %>% select(-age.aggr)

ddharm_bf_census_f_raw <- ddharm_bf_census_f %>%
  mutate(age.aggr = ifelse(AgeLabel %in% c("0","0-4","1-4"),"0-4", 
                           ifelse(AgeLabel %in% c(sprintf("%i-%i", seq(open.age, 90, by=5), seq(open.age+4, 94, by=5)), "80+", "85+", "90+", "95+", "98+"), paste0(open.age,"+"), AgeLabel))) %>%
  group_by(age.aggr) %>%
  summarise_at(vars(-(AgeStart:AgeLabel)),sum,na.rm=T) %>%
  mutate(age = as.numeric(regmatches(age.aggr, regexpr("\\d+",age.aggr))), .before = 1) %>%
  arrange(age) %>% select(-age.aggr)

##WPP Pop Estimates
pop <- wpp.pop %>% filter(Name == country) %>% 
  setNames(c(names(wpp.pop)[1:3], seq(0, 100, by = 5))) %>%
  reshape2::melt(id.vars=c("Name", "Sex", "Reference")) %>% 
  mutate(year = as.numeric(str_extract(Reference,"\\d{4}")), 
         value = as.numeric(value) * 1000,
         age = variable) %>%
  select(Sex, year, age, value) %>%
  pivot_wider(values_from = value, names_from = year)

pop.f <- filter(pop, Sex=="female") %>% select(-Sex)
pop.m <- filter(pop, Sex=="male") %>% select(-Sex)

pop.f.aggr <- pop.f[1:n_ages,]
pop.m.aggr <- pop.m[1:n_ages,]

pop.f.aggr[n_ages, -1] <- t(pop.f %>% slice(-(1:(n_ages-1))) %>% select(-1) %>% colSums())
pop.m.aggr[n_ages, -1] <- t(pop.m %>% slice(-(1:(n_ages-1))) %>% select(-1) %>% colSums())

##Pop data from Mark's package
#data(burkina_faso_females)

bf.idx5<-projection_indices(period_start = 1960,
                            period_end = 2015,
                            interval = 5,
                            n_ages = n_ages,
                            fx_idx = 4L,
                            n_fx = 7L,
                            n_sexes = 1)

##DHS data cohort smoothed
bf5.smooth <- aggr.mat.cohort.0[[country]] %>%
  mutate(age5 = 5 * floor(agegr / 5),
         period5 = 5 * floor(period / 5)) %>%
  group_by(mm1, tips, DHS, age5, period5) %>%
  summarise_at(vars(pyears, event, pyears2, adjusted), sum) %>%
  ungroup()

bf5.smooth$period5 <- factor(bf5.smooth$period5,levels=bf.idx5$periods)
bf5.smooth$tips <- factor(bf5.smooth$tips,levels=0:14)
bf5.f.smooth <- bf5.smooth %>% filter(mm1=="female") %>% arrange(period5,tips,age5)
bf5.f.no0.smooth <- filter(bf5.f.smooth, age5  >= 15)
bf5.f.0.smooth <- filter(bf5.f.smooth, age5 == 0)

tips.DX.smooth <- as(model.matrix(event~as.factor(tips)-1,data=bf5.f.no0.smooth),"sparseMatrix")
tips.DX.0.smooth <- as(model.matrix(event~as.factor(tips)-1,data=bf5.f.0.smooth),"sparseMatrix")

LQcoef.f <- LQcoef %>% filter(sex=="Female", !age%in%c("0","1-4")) %>% select(ax:vx)
LQcoef.f$age5 <- seq(5,110,by=5)

LQ.baseline.DX.smooth <- model.matrix(adjusted~factor(age5, levels = seq(5, 110, by = 5)) - 1, data=bf5.f.no0.smooth) %*% as.matrix(LQcoef.f[,1:4])
LQ.baseline.DX.ax.smooth <- LQ.baseline.DX.smooth[,1]
LQ.baseline.DX.bx.smooth <- data.frame(period5 = as.factor(bf.idx5$periods)) %>% 
  bind_rows(bind_cols(n = 1 : nrow(bf5.f.no0.smooth), ax=LQ.baseline.DX.smooth[,2], period5=bf5.f.no0.smooth$period5)) %>% 
  pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
  slice(-1) %>% select(-n) %>%
  as.matrix() %>% as("sparseMatrix")
LQ.baseline.DX.cx.smooth <- data.frame(period5 = as.factor(bf.idx5$periods)) %>% 
  bind_rows(bind_cols(n = 1 : nrow(bf5.f.no0.smooth), ax=LQ.baseline.DX.smooth[,3], period5=bf5.f.no0.smooth$period5)) %>% 
  pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
  slice(-1) %>% select(-n) %>%
  as.matrix() %>% as("sparseMatrix")
LQ.baseline.DX.vx.smooth <- data.frame(period5 = as.factor(bf.idx5$periods)) %>% 
  bind_rows(bind_cols(n = 1 : nrow(bf5.f.no0.smooth), ax=LQ.baseline.DX.smooth[,4], period5=bf5.f.no0.smooth$period5)) %>% 
  pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
  slice(-1) %>% select(-n) %>%
  as.matrix() %>% as("sparseMatrix")

basepop_init <- as.numeric(pop.f.aggr$'1960') %>% ifelse(.==0, 0.5, .)

#fx_init <- burkina.faso.females$fertility.rates[4:10, ]
#fx_init <- reshape2::melt(fx_init) %>% 
#  slice(c(rep(0:(ncol(fx_init)-1) * nrow(fx_init), each=nrow(fx_init) * 5) + 1:nrow(fx_init),
#          rep((ncol(fx_init)-1) * nrow(fx_init) + 1:nrow(fx_init), max(bf.idx5$periods - 2004)) )) %>%
#  mutate(year = rep(1960:max(bf.idx5$periods), each = nrow(fx_init))) %>%
#  select(-2) %>%
#  pivot_wider(names_from = year, values_from = value) %>%
#  select(num_range("",bf.idx5$periods)) %>%
#  as.matrix()

##WPP fx
fx <- wpp.fx %>% filter(Name == country) %>% reshape2::melt() %>% 
  mutate(year = as.numeric(str_extract(Reference,"\\d{4}")), 
         fx = value/1000,
         age = as.numeric(str_extract(variable,"\\d{2}"))) %>%
  select(year,age, fx)

fx_init <-  fx %>% filter(year %in% bf.idx5$periods) %>%
  pivot_wider(names_from=year, values_from=fx) %>% 
  select(-age) %>% as.matrix()

log_basepop_mean <- as.vector(log(basepop_init))
log_fx_mean <- as.vector(log(fx_init))

thiele.prior <- function(h, sex) {
  cat(h,"\n")
  log.m0 <- h - log(1 - 0.5 * exp(h)) + log(0.2)
  if(sex == "male"){
    LQ.mx <- c(log.m0, as.matrix(LQcoef.m[,1:3]) %*% c(1, h, h^2))
  } else if (sex == "female") {
    LQ.mx <- c(log.m0, as.matrix(LQcoef.f[,1:3]) %*% c(1, h, h^2))
  }
  
  #get priors for phi and psi using LQ mx at 0-4 to 10-14
  child.coef <- lm(LQ.mx[1:3] ~ c(2, 7, 12), weights = c(1e15, 1, 1))$coef
  phi <- exp(child.coef[1])
  psi <- -child.coef[2]
  
  #get priors for lambda, delta and epsilon using 10-14 to 40-44
  hump.coef <- lm(LQ.mx[3:9] ~  seq(12, 42, by = 5) + I(seq(12, 42, by = 5)^2))$coef
  delta <- -hump.coef[3]
  epsilon <- -hump.coef[2] / (2 * hump.coef[3])
  lambda <- exp(hump.coef[1] - hump.coef[2]^2 / (4 * hump.coef[3]))
  
  #get priors for A and B using 65-69 to 90-94
  old.coef <- lm(LQ.mx[14:19] ~ seq(67, 92, by = 5))$coef
  A <- exp(old.coef[1])
  B <- old.coef[2]
  
  thiele.min <- function(par, dat){
    all.age <- seq(2, 92, by = 5)
    par <- exp(par)
    est <- par[1] * exp(-par[2] * all.age) + par[3] * exp(-par[4]*(all.age - par[5])^2) + par[6] * exp(par[7] * all.age)
    ess <- sum((log(est[-1]) - dat[-1])^2) + 1e4 * (log(est[1]) - dat[1])^2
    return(ess)
  }
  nlm <- nlminb(start = log(c(phi, psi, lambda, delta, epsilon, A, B)), thiele.min, dat = LQ.mx[1:19], control = list(eval.max = 1000, iter.max = 1000))
  stopifnot(nlm$convergence ==0)
  return(setNames(exp(nlm$par), c("phi", "psi", "lambda", "delta", "epsilon", "A", "B")))
}

igme.h.mean.f <- igme.5q0.5 %>% filter(Sex=="Female", year5 %in% levels(bf5.f.no0.smooth$period5)) %>% .$child.mort %>% log()
thiele.prior.f <- sapply(igme.h.mean.f, thiele.prior, sex="female")

thiele_age <- seq(2, 112, by = 5)

data.f <- as.matrix(log(ddharm_bf_census_f_smoothed_aggr[,-1]))
#data.f <- as.matrix(log(ddharm_bf_census_f_raw[,-1]))
#data.f <- as.matrix(log(bf_census_f[,-1]))
#data.f <- as.matrix(log(bf.pop.aggr.f %>% select('1975', '1985', '1995', '2005', '2015')))

data.vec <- list(log_basepop_mean_f = log_basepop_mean,
                 log_fx_mean = log_fx_mean,
                 srb = rep(1.05, bf.idx5$n_periods),
                 interval = bf.idx5$interval,
                 n_periods = bf.idx5$n_periods,
                 fx_idx = 4L,
                 n_fx = 7L,
                 census_log_pop_f = data.f,
                 census_year_idx = match(bf.idx5$interval * floor(as.numeric(colnames(data.f)) / bf.idx5$interval), bf.idx5$periods),
                 census_year_grow_idx = as.numeric(colnames(data.f)) - bf.idx5$interval * floor(as.numeric(colnames(data.f)) / bf.idx5$interval),

                 open_idx = bf.idx5$n_ages,
                 pop_start = 3, pop_end = open.age / 5,
                
                 df = bf5.f.no0.smooth$adjusted,
                 Ef = bf5.f.no0.smooth$pyears2,
                 df_age = match(bf5.f.no0.smooth$age5+2, thiele_age),
                 df_time = match(bf5.f.no0.smooth$period5, levels(bf5.f.no0.smooth$period5)),
                 df_tp = c(bf5.f.no0.smooth$tips)-1,
                 
                 log_phi_mean = log(thiele.prior.f[1,]),
                 log_psi_mean = log(thiele.prior.f[2,]),
                 log_lambda_mean = log(thiele.prior.f[3,]),
                 log_delta_mean = log(thiele.prior.f[4,]),
                 log_epsilon_mean = log(thiele.prior.f[5,]),
                 log_A_mean = log(thiele.prior.f[6,]),
                 log_B_mean = log(thiele.prior.f[7,]),
                
                 #log_phi_mean = rep(mean(log(thiele.prior.f[1,])), bf.idx5$n_periods),
                 #log_psi_mean =  rep(mean(log(thiele.prior.f[2,])), bf.idx5$n_periods),
                 #log_lambda_mean =  rep(mean(log(thiele.prior.f[3,])), bf.idx5$n_periods),
                 #log_delta_mean =  rep(mean(log(thiele.prior.f[4,])), bf.idx5$n_periods),
                 #log_epsilon_mean =  rep(mean(log(thiele.prior.f[5,])), bf.idx5$n_periods),
                 #log_A_mean =  rep(mean(log(thiele.prior.f[6,])), bf.idx5$n_periods),
                 #log_B_mean =  rep(mean(log(thiele.prior.f[7,])), bf.idx5$n_periods),
                
                 thiele_age = thiele_age,
                 
                 penal_tp = as(crossprod(diff(diag(15))),"sparseMatrix"),
                 null_penal_tp = as(exp(15)*tcrossprod(c(0,1,1,1,rep(0,11))),"sparseMatrix"),
                 penal_tp_0 = as(tcrossprod(c(1,rep(0,14))),"sparseMatrix"),
                 
                 LQ_baseline_mx_DX_f = LQ.baseline.DX.ax.smooth,
                 h_DX_f = LQ.baseline.DX.bx.smooth,
                 h2_DX_f = LQ.baseline.DX.cx.smooth,
                 k_DX_f = LQ.baseline.DX.vx.smooth,
                 tp_DX_f = tips.DX.smooth,
                 LQ_baseline_f = as.matrix(LQcoef.f[,1:4]),
                 h_constant_f = igme.h.mean.f
)

par.vec <- list(log_tau2_logpop_f = c(1,0),
                #log_tau2_logpop_f_base=0,
                log_tau2_fx = 2,
                log_tau2_gx_f = 4,
                log_basepop_f = log_basepop_mean,
                log_fx = log_fx_mean,
                gx_f = rep(0, bf.idx5$n_ages * bf.idx5$n_periods),
                logit_rho_g_x = 4,
                logit_rho_g_t = 2,
                
                log_lambda_tp = 0,
                log_lambda_tp_0_inflated_sd = 0.3,
                tp_params = rep(0,15),
                
                log_dispersion = 1,
                
                log_phi_innov = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                log_psi_innov = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                log_lambda_innov = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                log_delta_innov = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                log_epsilon_innov = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                log_A_innov = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                log_B_innov = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                
                log_marginal_prec_phi = 0,
                log_marginal_prec_psi = 0,
                log_marginal_prec_lambda = 0,
                log_marginal_prec_delta = 0,
                log_marginal_prec_epsilon = 0,
                log_marginal_prec_A = 0,
                log_marginal_prec_B = 0,
                
                logit_rho_phi = 1,
                logit_rho_psi = 1,
                logit_rho_lambda = 1,
                logit_rho_delta = 1,
                logit_rho_epsilon = 1,
                logit_rho_A = 1,
                logit_rho_B = 1,
                
                log_marginal_prec_h = 0,
                log_marginal_prec_k = 0,
                h_params_f = igme.h.mean.f, 
                k_params_f = rep(0,bf.idx5$n_periods),
                logit_rho_h = 0,
                logit_rho_k = 0
)


input.LQ.both.vec <- list(data = data.vec, par_init = par.vec, model = "ccmpp_vr_tmb")

system.time(LQ.f.no0 <- fit_tmb(input.LQ.both.vec,inner_verbose=TRUE, random = c("log_basepop_f",
                                                                                 "log_fx",
                                                                                 "gx_f",
                                                                                 "tp_params",
                                                                                 "h_params_f",
                                                                                 "k_params_f"
                                                                                 ),
                                DLL="ccmpp_LQuad_bothsexes_AR_f_no0"
                                )
            ) 

system.time(thiele.f.no0 <- fit_tmb(input.LQ.both.vec,inner_verbose=TRUE, random = c("log_basepop_f",
                                                                                 "log_fx",
                                                                                 "gx_f",
                                                                                 "tp_params",
                                                                                 "log_phi_innov",
                                                                                 "log_psi_innov",
                                                                                 "log_lambda_innov",
                                                                                 "log_delta_innov",
                                                                                 "log_epsilon_innov",
                                                                                 "log_A_innov",
                                                                                 "log_B_innov"
                                                                                 ),
                                DLL="ccmpp_f_thiele"
                                )
            ) 

system.time(thiele.f.no0.MVN <- fit_tmb(input.LQ.both.vec,inner_verbose=TRUE, random = c("log_basepop_f",
                                                                                     "log_fx",
                                                                                     "gx_f",
                                                                                     "tp_params",
                                                                                     "log_phi_innov",
                                                                                     "log_psi_innov",
                                                                                     "log_lambda_innov",
                                                                                     "log_delta_innov",
                                                                                     "log_epsilon_innov",
                                                                                     "log_A_innov",
                                                                                     "log_B_innov"
                                                                                     ),
                                        DLL="ccmpp_f_thiele",
                                        map = list(logit_rho_phi = factor(NA),
                                                   logit_rho_psi = factor(NA),
                                                   logit_rho_lambda = factor(NA),
                                                   logit_rho_delta = factor(NA),
                                                   logit_rho_epsilon = factor(NA),
                                                   logit_rho_A = factor(NA),
                                                   logit_rho_B = factor(NA)
                                                   )
                                        )
            ) 


system.time(thiele.f.noDHS <- fit_tmb(input.LQ.both.vec,inner_verbose=TRUE, random = c("log_basepop_f",
                                                                                     "log_fx",
                                                                                     "gx_f",
                                                                                     "log_phi_innov",
                                                                                     "log_psi_innov",
                                                                                     "log_lambda_innov",
                                                                                     "log_delta_innov",
                                                                                     "log_epsilon_innov",
                                                                                     "log_A_innov",
                                                                                     "log_B_innov"
                                                                                     ),
                                    DLL="ccmpp_f_thiele_noDHS"
                                    )
            ) 

input.LQ.both.vec$par_init$logit_rho_g_t <- 0
system.time(thiele.f.fixgt <- fit_tmb(input.LQ.both.vec,inner_verbose=TRUE, random = c("log_basepop_f",
                                                                                       "log_fx",
                                                                                       "gx_f",
                                                                                       "tp_params",
                                                                                       "log_phi_innov",
                                                                                       "log_psi_innov",
                                                                                       "log_lambda_innov",
                                                                                       "log_delta_innov",
                                                                                       "log_epsilon_innov",
                                                                                       "log_A_innov",
                                                                                       "log_B_innov"
                                                                                       ),
                                      DLL="ccmpp_f_thiele",
                                      map=list(logit_rho_g_t = factor(NA))
                                      )
            ) 

models.list <- list("LQ" = LQ.f.no0, "Thiele" = thiele.f.no0, "Thiele MVN" = thiele.f.no0.MVN, "Thiele no DHS" = thiele.f.noDHS,
                     "Thiele fixgt" = thiele.f.fixgt)

#q4515####
q4515.func <- function(x){
  as_tibble(apply(x$mode$mx_mat_f[4:12,],2,function(x){1-prod((1-2.5*x)/(1+2.5*x))})) %>%
    mutate(sex="female", year = bf.idx5$periods)
  #%>%
  #  bind_rows(
  #    as_tibble(apply(x$mode$mx_mat_m[4:12,],2,function(x){1-prod((1-2.5*x)/(1+2.5*x))})) %>%
  #      mutate(sex="male", year = bf.idx5$periods)
  #  )
}
#wpp.bf.q4515<-read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/WPP_BF_45q15.csv")

wpp.bf.q4515 <- filter(wpp.q4515, Name==country)

q4515.df <- lapply(models.list, q4515.func) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>%
  bind_rows(
    as_tibble(reshape2::melt(wpp.bf.q4515[2,-c(1:4)]/1000)) %>% mutate(model="WPP Estimates", year=bf.idx5$periods, variable = NULL, sex = "female"),
    as_tibble(reshape2::melt(wpp.bf.q4515[1,-c(1:4)]/1000)) %>% mutate(model="WPP Estimates", year=bf.idx5$periods, variable = NULL, sex = "male")
  ) %>%
  mutate(model = fct_relevel(model, "WPP Estimates"))

q4515.df %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model, linetype = sex), lwd = 1.2) + ylab(bquote(""[45]*q[15])) +
  theme(text = element_text(size=25))

#child mx####
q50.func <- function(x){
  as_tibble(5 * x$mode$mx_mat_f[1,] / (1 + 2.5 * x$mode$mx_mat_f[1,])) %>%
    mutate(sex="female", year = bf.idx5$periods)
  #%>%
  #  bind_rows(
  #    as_tibble(apply(x$mode$mx_mat_m[4:12,],2,function(x){1-prod((1-2.5*x)/(1+2.5*x))})) %>%
  #      mutate(sex="male", year = bf.idx5$periods)
  #  )
}

q50.df <- lapply(models.list, q50.func) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>%
  bind_rows(
    as_tibble(exp(igme.h.mean.f)) %>% mutate(model="IGME Estimates", year=bf.idx5$periods, variable = NULL, sex = "female"),
  ) %>%
  mutate(model = fct_relevel(model, "IGME Estimates"))

q50.df %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model, linetype = sex), lwd = 1.2) + ylab(bquote(""[5]*q[0])) +
  theme(text = element_text(size=25))

##Comparison to DHS-spline####
load("~/more countries final avg sex Rwanda.RData")
skip<-c("Ethiopia","Central African Republic","Comoros","Sao Tome and Principe","Botswana","Cape Verde","Equatorial Guinea","Eritrea","Nigeria (Ondo State)","Ghana","Mauritania","Sudan")
joint.countries<-names(aggr.mat.cohort.0)[!names(aggr.mat.cohort.0)%in%skip]

no.basis = 15

bspline<-function (x,k,i,m=2) {
  if (m==-1) {basis<-as.numeric(x<k[i+1] & x>=k[i])} else {
    z0<-(x-k[i])/(k[i+m+1]-k[i])
    z1<-(k[i+m+2]-x)/(k[i+m+2]-k[i+1])
    basis<-z0*bspline(x,k,i,m-1)+z1*bspline(x,k,i+1,m-1) }
  basis
}

all.list<-function(x,age.start,age.end,year.start,year.end){
  no.basis=no.basis
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

thiele <- lapply(models.list, function(y){y$mode$mx_mat_f %>%
    `colnames<-`(bf.idx5$periods) %>%  
    `rownames<-`(seq(0, nrow(y$mode$mx_mat_f)*5-5, by = 5)) %>%
    reshape2::melt()}) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>%
  setNames(c("age5", "period5", "value", "model")) %>%
  mutate(period5 = sprintf("%d-%d",period5, period5+4)) %>%
  bind_rows(
    reshape2::melt(more.countries.avg[[country]]$mort.f) %>% 
      mutate(age5 = 5 * floor(Var1 / 5),
             period5 = 5 * floor(Var2 / 5)) %>% 
      group_by(age5, period5) %>%
      summarise_at(vars(value), mean) %>%
      mutate(value = exp(value),
             period5 = sprintf("%d-%d",period5, period5+4),
             model = "Spline average")
  ) %>%
  mutate(model = str_wrap(model, 20))

DHS.plot <- bf5.f.no0.smooth %>% mutate(value = adjusted / pyears2, period5 = as.numeric(levels(period5)[period5])) %>%
  select(tips, age5, period5, value) %>%
  mutate(period5 = sprintf("%d-%d", period5, period5 + 4))

ggplot(thiele %>% filter(age5 %in% 0:70)) + geom_line(aes(x = age5, y = value, col=model), lwd=1.2, linetype = 2) +
  geom_line(data = filter(thiele, model == "Spline average"), aes(x = age5, y = value, col = model), lwd=1.2, linetype=1)+
  geom_point(data = filter(DHS.plot, value!=0), aes(x = age5, y = value, shape = tips), alpha = 0) + 
  geom_point(data = filter(DHS.plot, value!=0, tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=3, stroke=2) +
  geom_point(data = filter(DHS.plot, value!=0, !tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=2, alpha = 0.3) +
  scale_shape_manual(values = 0:14) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  facet_wrap(~period5) + 
  ggtitle(country) + ylab(bquote(""[5]*m[x])) + xlab("5-year age groups") +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35))

par.deviation <- bind_rows(
  as_tibble(log(thiele.f.no0$mode$phi) - data.vec$log_phi_mean) %>% 
    mutate(par = "phi", period5 = bf.idx5$periods),
  
  as_tibble(log(thiele.f.no0$mode$psi) - data.vec$log_psi_mean) %>% 
    mutate(par = "psi", period5 = bf.idx5$periods),
  
  as_tibble(log(thiele.f.no0$mode$lambda) - data.vec$log_lambda_mean) %>% 
    mutate(par = "lambda", period5 = bf.idx5$periods),
  
  as_tibble(log(thiele.f.no0$mode$delta) - data.vec$log_delta_mean) %>% 
    mutate(par = "delta", period5 = bf.idx5$periods),
  
  as_tibble(log(thiele.f.no0$mode$epsilon) - data.vec$log_epsilon_mean) %>% 
    mutate(par = "epsilon", period5 = bf.idx5$periods),
  
  as_tibble(log(thiele.f.no0$mode$A) - data.vec$log_A_mean) %>% 
    mutate(par = "A", period5 = bf.idx5$periods),
  
  as_tibble(log(thiele.f.no0$mode$B) - data.vec$log_B_mean) %>% 
    mutate(par = "B", period5 = bf.idx5$periods)
) %>%
  mutate(par = fct_inorder(par))


ggplot(par.deviation) + geom_line(aes(x = period5, y = value, col = par), lwd = 1.2) +
  theme(text = element_text(size=20)) +
  ggtitle("Deviation from prior means")


#populaton counts####
#ADD AR.f$mode$census_proj
#ddharm_bf_census_f_smoothed_aggr$'1995' <- ddharm_bf_census_f_smoothed_aggr$'1985'^(1-10/11) *  ddharm_bf_census_f_smoothed_aggr$'1996'^(10/11)
#ddharm_bf_census_f_smoothed_aggr$'2005' <- ddharm_bf_census_f_smoothed_aggr$'1996'^(1-9/10) *  ddharm_bf_census_f_smoothed_aggr$'2006'^(9/10)

#ddharm_bf_census_f_raw$'1995' <- ddharm_bf_census_f_raw$'1985'^(1-10/11) *  ddharm_bf_census_f_raw$'1996'^(10/11)
#ddharm_bf_census_f_raw$'2005' <- ddharm_bf_census_f_raw$'1996'^(1-9/10) *  ddharm_bf_census_f_raw$'2006'^(9/10)


mf5 <- projection_model_frames(bf.idx5)
get.pop<- function(x){
  pop.mat <- matrix(x$mode$population_f, bf.idx5$n_ages, bf.idx5$n_periods+1)
  #pop.mat <- rbind(pop.mat, apply(pop.mat,2,sum))
  #rownames(pop.mat) <- c(bf.idx5$ages, "All")
  rownames(pop.mat) <- c(bf.idx5$ages)
  colnames(pop.mat) <- bf.idx5$periods_out
  reshape2::melt(pop.mat) %>% select(age = Var1, year = Var2, value = value)
}

pop.df <- lapply(models.list, get.pop) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>%
  bind_rows(
    ddharm_bf_census_f_raw %>%
      mutate(age = c(ddharm_bf_census_f_raw$age), model="UNPD Census") %>%
      pivot_longer(!age & !model) %>% mutate(year=as.numeric(name), name=NULL),
    
    ddharm_bf_census_f_smoothed_aggr %>%
      mutate(age = c(ddharm_bf_census_f_raw$age), model="UNPD Census smoothed") %>%
      pivot_longer(!age & !model) %>% mutate(year=as.numeric(name), name=NULL),
    
    mutate(pop.f.aggr[,-(1:3)], age = ddharm_bf_census_f_raw$age, model="WPP Estimates") %>%
      pivot_longer(!age & !model) %>% mutate(year=as.numeric(name), name=NULL)
  ) %>%
  mutate(cohort = year - age,
         cohort = 5 * floor(cohort / 5),
         model = fct_relevel(model, c("UNPD Census", "UNPD Census smoothed", "WPP Estimates"))
         )


pop.df %>% filter(year%%5==0) %>%
  ggplot() + geom_line(aes(x = age, y = value, col = model), lwd = 1.2) +
  theme(text = element_text(size=25)) +
  facet_wrap_paginate(~cohort, scale="free_y", nrow = 2, ncol = 3, page=3)


mf5 <- projection_model_frames(bf.idx5)
get.census.pop<- function(x){
  pop.mat <- matrix(x$mode$census_proj_mat_f, bf.idx5$n_ages, ncol(data.f))
  #pop.mat <- rbind(pop.mat, apply(pop.mat,2,sum))
  #rownames(pop.mat) <- c(bf.idx5$ages, "All")
  rownames(pop.mat) <- c(bf.idx5$ages)
  colnames(pop.mat) <- colnames(data.f)
  reshape2::melt(pop.mat) %>% select(age = Var1, year = Var2, value = value)
}

censpop.df <- lapply(models.list, get.census.pop) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>%
  bind_rows(
    ddharm_bf_census_f_raw %>%
      mutate(age = c(ddharm_bf_census_f_raw$age), model="UNPD Census") %>%
      pivot_longer(!age & !model) %>% mutate(year=as.numeric(name), name=NULL),
    
    ddharm_bf_census_f_smoothed_aggr %>%
      mutate(age = c(ddharm_bf_census_f_smoothed_aggr$age), model="UNPD Census smoothed") %>%
      pivot_longer(!age & !model) %>% mutate(year=as.numeric(name), name=NULL),
    
    mutate(pop.f.aggr[,-(1:3)], age = ddharm_bf_census_f_raw$age, model="WPP Estimates") %>%
      pivot_longer(!age & !model) %>% mutate(year=as.numeric(name), name=NULL)
  ) %>%
  mutate(cohort = year - age,
         model = fct_relevel(model, c("UNPD Census", "UNPD Census smoothed", "WPP Estimates"))
  )

allpop.df <- full_join(pop.df, censpop.df) %>% mutate(model = fct_relevel(model, c("UNPD Census", "UNPD Census smoothed", "WPP Estimates")))


allpop.df %>% filter(year %in% c(1960,colnames(data.f))) %>%
  ggplot() + geom_line(aes(x = age, y = value, col = model), lwd = 1.2) +
  theme(text = element_text(size=25)) + ylab("Population counts") + 
  facet_wrap_paginate(~year, scale="free_y", nrow = 2, ncol = 3, page=1)

allpop.df %>% filter(model != "WPP Estimates") %>%
  ggplot() + geom_line(aes(x = age, y = value, col = model), lwd = 1.2) +
  theme(text = element_text(size=25)) +
  facet_wrap_paginate(~cohort, scale="free_y", nrow = 2, ncol = 2, page=2)

censpop.df %>% filter(model != "WPP Estimates") %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) +
  theme(text = element_text(size=25)) + ylab("Population counts") +
  facet_wrap_paginate(~age, scale="free", nrow = 4, ncol = 3, page=1)

allpop.df %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) +
  theme(text = element_text(size=25)) +
  facet_wrap_paginate(~age, scale="free", nrow = 4, ncol = 3, page=1)

#migration####
get.mig<- function(x){
  pop.mat <- matrix(x$mode$migrations_f, bf.idx5$n_ages, bf.idx5$n_periods)
  #pop.mat <- rbind(pop.mat, apply(pop.mat,2,sum))
  #rownames(pop.mat) <- c(bf.idx5$ages, "All")
  rownames(pop.mat) <- c(bf.idx5$ages)
  colnames(pop.mat) <- bf.idx5$periods
  reshape2::melt(pop.mat) %>% select(age = Var1, year = Var2, mig = value)
}

mig.df <- lapply(models.list, get.mig) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows()

inner_join(mig.df, pop.df) %>% mutate(gx = mig/value) %>%
  filter(year %in% c(1960, 1975, 1985, 1995, 2005, 2015)) %>%
  ggplot() + geom_line(aes(x = age, y = gx, col = model), lwd = 1.2) +
  theme(text = element_text(size=25)) +
  facet_wrap_paginate(~ year, nrow=2, ncol=3, page=1, scales = "free_y")

inner_join(mig.df, pop.df) %>% mutate(gx = mig/value) %>%
  ggplot() + geom_line(aes(x = year, y = gx, col = model), lwd = 1.2) +
  theme(text = element_text(size=25)) +
  facet_wrap_paginate(~ age, scales="free_y", nrow=2, ncol=3, page=2)

inner_join(mig.df, pop.df) %>% mutate(gx = mig/value) %>%
  ggplot() + geom_line(aes(x = age, y = gx, col = model), lwd = 1.2) +
  theme(text = element_text(size=25)) +
  facet_wrap_paginate(~ cohort, nrow=2, ncol=3, page=3)



#fertility####
get.fert<- function(x){
  pop.mat <- matrix(x$mode$fx, bf.idx5$n_fx, bf.idx5$n_periods)
  rownames(pop.mat) <- c(bf.idx5$fertility_ages)
  colnames(pop.mat) <- bf.idx5$periods
  reshape2::melt(pop.mat) %>% select(age = Var1, year = Var2, value = value)
}

fx.df <- lapply(models.list, get.fert) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>%
  bind_rows(
    reshape2::melt(fx_init %>% `rownames<-`(bf.idx5$fertility_ages)) %>%
      select(age = Var1, year = Var2, value = value) %>%
      mutate(model = "Initial Values")
  ) %>%
  mutate(model = fct_relevel(model, "Initial Values")) 

fx.df %>%
  ggplot() + geom_line(aes(x = age, y = value, col = model), lwd = 1.2) +
  theme(text = element_text(size=25)) +
  facet_wrap_paginate(~year, scale="free_y", nrow = 2, ncol = 4, page=1)

fx.df %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) +
  theme(text = element_text(size=25)) +
  facet_wrap_paginate(~age, scale="free_y", nrow = 2, ncol = 4, page=1)
