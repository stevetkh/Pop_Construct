#devtools::load_all("C:/Users/ktang3/Documents/GitHub/leapfrog")
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

load(file="C:/Users/ktang3/Desktop/Imperial/Pop_Construct/swe_fert.rda")
load(file="C:/Users/ktang3/Desktop/Imperial/Pop_Construct/swe_lt.rda")
mlt <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/MLT_UN2011_130_1y_complete.csv")

#compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/trial.cpp")
#dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/trial"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_singlesex.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_singlesex"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_AR.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_AR"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_AR_vec.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_AR_vec"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_weighted_hmean.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_weighted_hmean"))

projection_indices <- function(period_start, interval, n_periods, n_ages,
                               fx_idx, n_fx, n_sexes = 1) {
  
  stopifnot(n_sexes %in% 1:2)
  
  periods_out <- period_start + 0:n_periods * interval
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
       n_periods = n_periods,
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

#compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_vr_tmb.cpp")
#dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_vr_tmb"))

#' Convert projection output to long format
#'
#' @param proj CCMPP projection output from [`ccmppR()`].
#' @param mf model frames list from [`projection_model_frames()`].
#' @param value_col column name for projection output.
#'
proj_to_long <- function(proj, mf, value_col = "value") {
  
  val <- list()
  val$population <- mf$mf_population
  val$population[[value_col]] <- as.vector(proj$population)
  
  val$population <- mf$mf_cohort_deaths
  val$population[[value_col]] <- as.vector(proj$cohort_deathsg)
  
}

#' ## Baseline 5x5 popualtion projection
#'
#' Construct a population projection from 1900 to 2010 for age 0 to 90+ by
#' 5 year age group with no migration.
#'
#' The first step is to calculate life table sx values. Since
#' $s_x = _nL_x / _nL_{x-n}$ and for the open age group is approximated as
#' $s_x = T_x / T_{x-n}$, truncate the life table to one age older (age 95)
#' than the oldest age in the population projection (90+).
#'
#' There is a bit of imprecision here by using `exposure` as population size
#' instead of getting the start of year population data.

###Sweden #####
##+ swe_lt_sx
swe_lt5 <- swe_lt %>%
  arrange(iso3, sex, year, age) %>%
  mutate(population = exposure * as.integer(year %% 5 == 0),
         year = 5 * floor(year / 5),
         age = pmin(5 * floor(age / 5), 95)) %>%
  group_by(iso3, sex, year, age) %>%
  summarise_at(vars(population, exposure, deaths, Lx, dx ,lx), sum) %>%
  group_by(iso3, sex, year) %>%
  mutate(Tx = rev(cumsum(rev(Lx))),
         sx = if_else(age < max(age),
                      Lx / lag(Lx, default = 25 * 1e5),
                      Tx / lag(Tx)),
         mx = dx/Lx,
         mx5 = dx/population,
         qx5 = mx5/(1+0.5*mx5)) %>%  
  ungroup()

## mlt5
mlt5 <- as_tibble(mlt) %>%
  mutate(age5 =  pmin(5 * floor(age / 5), 95)) %>%
#  subset(!is.na(sx1)) %>%
  group_by(Type_MLT, Family, Type, Sex, E0, age5) %>%
  summarise_at(vars(Lx1,dx1), sum) %>%
  group_by(Type_MLT, Family, Type, Sex, E0) %>%
  mutate(Tx = rev(cumsum(rev(Lx1))),
         sx = if_else(age5 < max(age5),
                      Lx1 / lag(Lx1, default = 5 * 1e5),
                      Tx / lag(Tx)),
         mx = dx1/Lx1,
        mx5 = 5*dx1/Lx1) %>%
  ungroup()


#' For inputs, use the life table 'exposure' for 1900 as the base population.
#' Assume that fertility occurs for women age 15 to 49, and the sex ratio at
#' birth is 1.05

##+ ccmpp_inputs

idx5 <- projection_indices(period_start = 1900,
                           interval = 5,
                           n_periods = 22,
                           n_ages = 19,
                           fx_idx = 4L,
                           n_fx = 7,
                           n_sexes = 1)

mf5 <- projection_model_frames(idx5)


swe_lt5_sx <- swe_lt5 %>%
  filter(sex == "female",
         year %in% idx5$periods) %>%
  arrange(iso3, sex, year, age)

swe_lt5_pop <- swe_lt5 %>%
  filter(sex %in% idx5$sexes,
         year %in% idx5$periods_out) %>%
  mutate(age = pmin(age, max(idx5$ages))) %>%
  group_by(iso3, sex, year, age) %>%
  summarise_at(vars(population, exposure, deaths), sum) %>%
  ungroup()

swe_fert5 <- swe_fert %>%
  mutate(year = 5 * floor(year / 5),
         age = 5 * floor(age / 5)) %>%
  group_by(iso3, year, age) %>%
  summarise_at(vars(births, exposure), sum) %>%
  mutate(asfr = births / exposure) %>%
  filter(year %in% idx5$periods,
         age %in% idx5$fertility_ages) %>%
  arrange(iso3, year, age)


basepop <- filter(swe_lt5_pop, year == min(idx5$periods))$population
sx <- matrix(swe_lt5_sx$sx, length(idx5$ages)+1)
fx <- matrix(swe_fert5$asfr, length(idx5$fertility_ages))
gx <- matrix(0, length(idx5$ages), length(idx5$periods))
srb <- rep(1.05, idx5$n_periods)

baseline_sx<-subset(mlt5, Type=="CD West" & match(tolower(Sex) , tolower(idx5$sexes)) & E0==45)

#' # PopReconstruct
#'
basepop_init <- basepop
sx_init <- sx
fx_init <- pmax(fx, min(fx[fx > 0]))
gx_init <- gx

log_basepop_mean <- as.vector(log(basepop_init))
logit_sx_mean <- as.vector(qlogis(sx_init))
log_fx_mean <- as.vector(log(fx_init))
gx_mean <- as.vector(gx_init)

census_log_pop <- matrix(log(swe_lt5_pop$population),
                         nrow = length(idx5$ages))
census_year_idx <- seq.int(1, idx5$n_periods + 1, by = 5 / idx5$interval)

deaths_obs <- swe_lt5_pop %>%
  filter(year != max(year)) %>%
  .$deaths %>%
  matrix(nrow = idx5$n_ages)

births_obs <- matrix(swe_fert5$births, nrow = idx5$n_fx)

data <- list(log_basepop_mean = log_basepop_mean,
             logit_sx_mean = logit_sx_mean,
             log_fx_mean = log_fx_mean,
             gx_mean = gx_mean,
             srb = srb,
             interval = idx5$interval,
             n_periods = idx5$n_periods,
             fx_idx = idx5$fx_idx,
             n_fx = idx5$n_fx,
             census_log_pop = census_log_pop[ , census_year_idx, drop = FALSE],
             census_year_idx = census_year_idx,
             sx_params_mean = rep(1,idx5$n_periods),
             deaths_obs = deaths_obs,
             logit_baseline_sx = as.vector(qlogis(baseline_sx$sx)),
             #age 20-40
             hiv_mort_start = 6, hiv_mort_end = 10,
             log_hiv_mort_mean = rep(0, idx5$n_periods)
             )

par <- list(log_tau2_logpop = 0,
            log_tau2_sx = 0,
            log_tau2_fx = 0,
            log_tau2_gx = 0,
            log_basepop = log_basepop_mean,
            logit_sx = logit_sx_mean,
            log_fx = log_fx_mean,
            gx = gx_mean,
            sx_params = rep(1,idx5$n_periods),
            log_dispersion=2,
            log_tau2_hiv_mort = 0,
            log_hiv_mort = rep(0,idx5$n_periods))

make_tmb_obj <- function(data,
                         par,
                         model = "ccmpp_tmb",
                         inner_verbose = FALSE,
                         calc_outputs = TRUE,
                         random,
                         DLL = "leapfrog_TMBExports") {
  
  data$model <- model
  data$calc_outputs <- as.integer(calc_outputs)
  
  obj <- TMB::MakeADFun(data = data,
                        parameters = par,
                        DLL = DLL,
                        silent = !inner_verbose,
                        random = random)
  class(obj) <- "tmb_obj"
  
  obj
}

fit_tmb <- function(tmb_input,
                    outer_verbose = TRUE,
                    inner_verbose = FALSE,
                    max_iter = 250,
                    random
) {
  
  ## stopifnot(inherits(tmb_input, "tmb_input"))
  
  if(is.null(tmb_input$model))
    tmb_input$model <- "ccmpp_tmb"
  
  obj <- make_tmb_obj(data = tmb_input$data,
                      par = tmb_input$par_init,
                      model = tmb_input$model,
                      inner_verbose = inner_verbose,
                      calc_outputs = 0L,
                      random = random)
  
  trace <- if(outer_verbose) 1 else 0
  f <- withCallingHandlers(
    stats::nlminb(obj$par, obj$fn, obj$gr,
                  control = list(trace = trace,
                                 iter.max = max_iter)),
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
                         random = random)
  f$mode <- objout$report(f$par.full)
  
  val <- c(f, obj = list(objout))
  class(val) <- "leapfrog_fit"
  
  val
}


obj <- make_tmb_obj(data, par, inner_verbose = TRUE , model = "ccmpp_vr_tmb", 
                    random =  c("log_basepop", "sx_params", "log_hiv_mort", "log_fx", "gx"))

input <- list(data = data, par_init = par, model = "ccmpp_vr_tmb")
input2 <- list(data = data, par_init = par, model = "ccmpp_tmb")

init_sim <- obj$report()

#obj$fn()
#obj$gr()

fit <- fit_tmb(input,inner_verbose=TRUE, random = c("log_basepop", "sx_params", "log_hiv_mort", "log_fx", "gx"))

fit2 <- fit_tmb(input2,inner_verbose=TRUE)

#fit <- sample_tmb(fit)


mf5$mf_population %>%
  mutate(population = fit$mode$population,
         source = "projection (no migration)") %>%
  bind_rows(
    swe_lt5_pop %>%
      select(age, period = year, population) %>%
      mutate(source = "HMD exposure denominator")
  ) %>%
  mutate(source = fct_inorder(source)) %>%
  count(period, source, wt = population) %>%
  ggplot(aes(period, n, linetype = source)) +
  geom_line() +
  scale_y_continuous("Population (millions)",
                     labels = scales::label_number(scale = 1e-6)) +
  expand_limits(y = 0) +
  ggtitle("Sweden population 1900-2010")


if( !file.exists("swe5yr_fitstan.rds") ) {
  fitstan <- fit_tmbstan(input, chains = 4, iter = 2000,
                         use_inits = TRUE, refresh = 100)
  saveRDS(fitstan, "swe5yr_fitstan.rds")
} else {
  fitstan <- readRDS("swe5yr_fitstan.rds")
}

fitstan <- sample_tmbstan(fitstan, TRUE)


n_samples <- ncol(fit$sample$population)

colnames(fit$sample$population) <- seq_len(n_samples)
colnames(fit$sample$migrations) <- seq_len(n_samples)
colnames(fit$sample$fx) <- seq_len(n_samples)

init_pop <- ccmppR(basepop_init, sx_init, fx_init, gx_init,
                   srb = srb, interval = 5, fx_idx = 4)

df <- crossing(year = c(years, max(years)+1),
               sex = "female",
               age = ages) %>%
  mutate(init_pop = as.vector(init_pop_mat))

census_pop <- crossing(sex = "female",
                       age_group = c(sprintf("%02d-%02d", 0:15*5, 0:15*5+4), "80+")) %>%
  bind_cols(as_tibble(burkina.faso.females$census.pop.counts)) %>%
  gather(year, census_pop, `1975`:`2005`) %>%
  type_convert(cols(year = col_double()))

df <- df %>%
  left_join(census_pop) %>%
  bind_cols(as_tibble(fit$sample$population)) %>%
  gather(sample, value, `1`:last_col())

agepop <- df %>%
  group_by(year, sex, age_group) %>%
  summarise(init_pop = mean(init_pop),
            census_pop = mean(census_pop),
            mean = mean(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

totalpop <- df %>%
  group_by(year, sample) %>%
  summarise(init_pop = sum(init_pop),
            census_pop = sum(census_pop),
            value = sum(value)) %>%
  group_by(year) %>%
  summarise(init_pop = mean(init_pop),
            census_pop = mean(census_pop),
            mean = mean(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))


##+ Migrations

migrations <- crossing(year = c(years, years),
                       age_lower = ages) %>%
  mutate(init_migrations = init_sim$migrations) %>%
  bind_cols(as_tibble(fit$sample$migrations)) %>%
  gather(sample, value, `1`:last_col()) 

total_migrations <- migrations %>%
  group_by(year, sample) %>%
  summarise(init_migrations = sum(init_migrations),
            value = sum(value)) %>%
  group_by(year) %>%
  summarise(init_migrations = mean(init_migrations),
            mean = mean(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

total_migrations <- migrations %>%
  group_by(year, sample) %>%
  summarise(init_migrations = sum(init_migrations),
            value = sum(value)) %>%
  group_by(year) %>%
  summarise(init_migrations = mean(init_migrations),
            mean = mean(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))


ggplot(total_migrations, aes(year, mean, ymin = lower, ymax = upper)) +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  geom_line(aes(y = init_migrations), linetype = "dashed") +
  scale_y_continuous("Net migration (000s)",
                     labels = scales::number_format(scale=1e-3)) +
  labs(x = NULL) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  ggtitle("SWE: total net migrations (thousands)")


migrations %>%
  mutate(age_group = cut(age_lower, c(0, 15, 50, Inf),
                         c("0-14", "15-49", "50+"), right = FALSE)) %>%
  group_by(year, age_group, sample) %>%
  summarise(init_migrations = sum(init_migrations),
            value = sum(value)) %>%
  summarise(init_migrations = mean(init_migrations),
            mean = mean(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975)) %>%
  ggplot(aes(year, mean, ymin = lower, ymax = upper,
             color = age_group, fill = age_group)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line() +
  geom_line(aes(y = init_migrations), linetype = "dashed") +
  scale_y_continuous("Net migration (000s)",
                     labels = scales::number_format(scale=1e-3)) +
  labs(x = NULL) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  ggtitle("SWE: net migrations by age (thousands)")

###Burkina Faso####
library(MortCast)
load("~/cohort smooth 1900-2017.RData")

bf.pop <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/Burkina Faso Pop.csv")
colnames(bf.pop) <- c("Name", "Sex", "Year", sprintf("%d-%d", seq(0, 95, by=5), seq(4, 99, by=5)), "100+")

bf.pop.f <- bf.pop %>% filter(Sex=="Female") %>% select(-c("Name", "Sex")) %>% apply(2,function(x)as.numeric(gsub(" ","",x))) %>% as_tibble() %>%
            pivot_longer(!Year,names_to="Age",values_to="counts") %>% pivot_wider(names_from=Year,values_from="counts")
bf.pop.f[,-1] <- bf.pop.f[,-1]*1000
  
bf.pop.m <- bf.pop %>% filter(Sex=="Male") %>% select(-c("Name", "Sex")) %>% apply(2,function(x)as.numeric(gsub(" ","",x))) %>% as_tibble() %>%
  pivot_longer(!Year,names_to="Age",values_to="counts") %>% pivot_wider(names_from=Year,values_from="counts")
bf.pop.m[,-1] <- bf.pop.m[,-1]*1000

data(burkina_faso_females)

n_ages = nrow(burkina.faso.females$census.pop.counts)

projection_indices <- function(period_start, interval, n_periods, n_ages,
                               fx_idx, n_fx, n_sexes = 1) {
  
  stopifnot(n_sexes %in% 1:2)
  
  periods_out <- period_start + 0:n_periods * interval
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
       n_periods = n_periods,
       n_ages = n_ages,
       fx_idx = fx_idx,
       n_fx = n_fx,
       n_sexes = n_sexes)
}

bf.idx5<-projection_indices(period_start = 1960,
                            interval = 1,
                            age_interval = 5,
                            n_ages = n_ages,
                            n_periods = 2015-1960+1,
                            fx_idx = 4L,
                            n_fx = 7L,
                            n_sexes = 2)

bf.pop.aggr.f <- bf.pop.f[1:bf.idx5$n_ages,]
bf.pop.aggr.f[bf.idx5$n_ages,-1] <- t(bf.pop.f %>% slice(-(1:(bf.idx5$n_ages-1))) %>% select(-1) %>% colSums())
bf.pop.aggr.f[bf.idx5$n_ages,1] <- "80+"
bf.pop.aggr.f <- bf.pop.aggr.f %>% select(num_range("",bf.idx5$periods))

bf.pop.aggr.m <- bf.pop.m[1:bf.idx5$n_ages,]
bf.pop.aggr.m[bf.idx5$n_ages,-1] <- t(bf.pop.m %>% slice(-(1:(bf.idx5$n_ages-1))) %>% select(-1) %>% colSums())
bf.pop.aggr.m[bf.idx5$n_ages,1] <- "80+"
bf.pop.aggr.m <- bf.pop.aggr.m %>% select(num_range("",bf.idx5$periods))

## model life tables
{mlt <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/MLT_UN2011_130_1y_complete.csv")
mlt5 <- as_tibble(mlt) %>%
  mutate(age5 =  pmin(5 * floor(age / 5), max(bf.idx5$ages))) %>%
  #  subset(!is.na(sx1)) %>%
  group_by(Type_MLT, Family, Type, Sex, E0, age5) %>%
  summarise_at(vars(Lx1,dx1), sum) %>%
  group_by(Type_MLT, Family, Type, Sex, E0) %>%
  mutate(Tx = rev(cumsum(rev(Lx1))),
         sx = if_else(age5 < max(age5),
                      Lx1 / lag(Lx1, default = 5 * 1e5),
                      Tx / lag(Tx)),
         mx = dx1/Lx1,
         mx5 = 5*dx1/Lx1) %>%
  ungroup()

log_baseline_mx<-log(subset(mlt5, Type=="CD West" & match(tolower(Sex) , "female") & E0==45)$mx)
}

##Log Quad
LQcoef.f<-rbind(data.frame(ax=0, bx=1, cx=0, vx=0) , LQcoef %>% filter(sex=="Female", !age%in%c("0","1-4")) %>% select(ax:vx))
LQcoef.f$age5 <- seq(0,110,by=5)

LQcoef.m<-rbind(data.frame(ax=0, bx=1, cx=0, vx=0) , LQcoef %>% filter(sex=="Male", !age%in%c("0","1-4")) %>% select(ax, bx, cx, vx))
LQcoef.m$age5 <- seq(0,110,by=5)

library(rdhs)
library(demogsurv)
library(data.table)
library(survival)
countries <- dhs_countries()
surveys <- dhs_surveys(countryIds = "BF", surveyYearStart=1900, surveyType = "DHS",surveyCharacteristicIds = 1)
ird <- dhs_datasets(fileType = "IR", fileFormat = "flat")[SurveyId %in% surveys$SurveyId]
ird$path <- unlist(get_datasets(ird$FileName))
ir <- list()
for(survid in ird$SurveyId){
  print(survid)
  dat <- readRDS(ird[SurveyId == survid]$path)
  dat <- dat[grep("caseid|^v0|^v1|^b|^mm", names(dat))]
  ir[[survid]] <- dat
}
ir <- lapply(ir, haven::as_factor)
ir <- Map(data.frame,
          SurveyId = surveys$SurveyId,
          CountryName = surveys$CountryName,
          SurveyYear = surveys$SurveyYear,
          ir)

zw<-lapply(ir,reshape_sib_data)
for(i in 1:length(zw)){
  zw[[i]]$death<-factor(zw[[i]]$mm2,c("dead","alive"))=="dead"
  zw[[i]]$tstop<-ifelse(zw[[i]]$death,zw[[i]]$mm8,zw[[i]]$v008)
  zw[[i]]$weight<-zw[[i]]$v005/1e6
}

aggr<-list()
for(i in 1:length(zw)){
  aggr[[names(zw)[i]]] <- demog_pyears(~mm1, zw[[i]],
                                       period = bf.idx5$periods,
                                       agegr = seq(0, 70, 5),
                                       tips = 0:15,
                                       dob = "mm4",
                                       intv = "v008",
                                       event = "death",
                                       tstart = "mm4",
                                       tstop = "tstop",
                                       weights = "weight")$data
}
aggr.mat.form.list <- function(x){
  aggr.mat<-x  
  aggr.mat<-aggr.mat[!aggr.mat$mm1 %in% c("missing",9),]
  aggr.mat.reduced<-aggr.mat[aggr.mat$pyears!=0,]
  aggr.mat.reduced
}

bf5.list <- lapply(aggr,aggr.mat.form.list)
bf5 <- c()
for(i in 1:length(bf5.list)) {bf5 <- rbind(bf5,cbind(bf5.list[[i]],DHS=names(bf5.list)[i]))}
bf5$period5 <- factor(gsub("-\\d{4,}","",bf5$period),levels=bf.idx5$periods)
bf5$age5 <- factor(gsub("-\\d+","",bf5$agegr),levels=bf.idx5$ages)
bf5$tips <- factor(bf5$tips,levels=0:14)
bf5.f <- bf5 %>% filter(mm1=="female") %>% arrange(period5,tips,age5)
bf5.m <- bf5 %>% filter(mm1=="male") %>% arrange(period5,tips,age5)

baseline.DX <- model.matrix(event~as.factor(age5)-1,data=bf5.f) %*% log_baseline_mx
year.DX <- as(model.matrix(event~as.factor(period5)-1,data=bf5.f),"sparseMatrix")
tips.DX <- as(model.matrix(event~as.factor(tips)-1,data=bf5.f),"sparseMatrix")
tips.DX.m <- as(model.matrix(event~as.factor(tips)-1,data=bf5.m),"sparseMatrix")

LQ.baseline.DX <- model.matrix(event~factor(age5, levels = seq(0, 110, by = 5)) - 1, data=bf5.f) %*% as.matrix(LQcoef.f[,1:4])
LQ.baseline.DX.ax <- LQ.baseline.DX[,1]
LQ.baseline.DX.bx <- data.frame(period5 = as.factor(bf.idx5$periods)) %>% 
                      bind_rows(bind_cols(n = 1 : nrow(bf5.f), ax=LQ.baseline.DX[,2], period5=bf5.f$period5)) %>% 
                      pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
                      slice(-1) %>% select(-n) %>%
                      as.matrix() %>% as("sparseMatrix")
LQ.baseline.DX.cx <- data.frame(period5 = as.factor(bf.idx5$periods)) %>% 
                      bind_rows(bind_cols(n = 1 : nrow(bf5.f), ax=LQ.baseline.DX[,3], period5=bf5.f$period5)) %>% 
                      pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
                      slice(-1) %>% select(-n) %>%
                      as.matrix() %>% as("sparseMatrix")
LQ.baseline.DX.vx <- data.frame(period5 = as.factor(bf.idx5$periods)) %>% 
                      bind_rows(bind_cols(n = 1 : nrow(bf5.f), ax=LQ.baseline.DX[,4], period5=bf5.f$period5)) %>% 
                      pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
                      slice(-1) %>% select(-n) %>%
                      as.matrix() %>% as("sparseMatrix")


LQ.baseline.DX.m <- model.matrix(event~factor(age5, levels = seq(0, 110, by = 5)) - 1, data=bf5.m) %*% as.matrix(LQcoef.m[,1:4])
LQ.baseline.DX.ax.m <- LQ.baseline.DX.m[,1]
LQ.baseline.DX.bx.m <- data.frame(period5 = as.factor(bf.idx5$periods)) %>% 
  bind_rows(bind_cols(n = 1 : nrow(bf5.m), ax=LQ.baseline.DX.m[,2], period5=bf5.m$period5)) %>% 
  pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
  slice(-1) %>% select(-n) %>%
  as.matrix() %>% as("sparseMatrix")
LQ.baseline.DX.cx.m <- data.frame(period5 = as.factor(bf.idx5$periods)) %>% 
  bind_rows(bind_cols(n = 1 : nrow(bf5.m), ax=LQ.baseline.DX.m[,3], period5=bf5.m$period5)) %>% 
  pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
  slice(-1) %>% select(-n) %>%
  as.matrix() %>% as("sparseMatrix")
LQ.baseline.DX.vx.m <- data.frame(period5 = as.factor(bf.idx5$periods)) %>% 
  bind_rows(bind_cols(n = 1 : nrow(bf5.m), ax=LQ.baseline.DX.m[,4], period5=bf5.m$period5)) %>% 
  pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
  slice(-1) %>% select(-n) %>%
  as.matrix() %>% as("sparseMatrix")

basepop_init <- as.numeric(bf.pop.aggr.f$'1960')
basepop_init.m <- as.numeric(bf.pop.aggr.m$'1960')

sx_init <- burkina.faso.females$survival.proportions
sx_init <- cbind(sx_init, `2005` = sx_init[ , "2000"], `2010` = sx_init[ , "2000"], `2015` = sx_init[ , "2000"])

fx_init <- burkina.faso.females$fertility.rates[4:10, ]
fx_init <- reshape2::melt(fx_init) %>% 
  slice(c(rep(0:(ncol(fx_init)-1) * nrow(fx_init), each=nrow(fx_init) * 5) + 1:nrow(fx_init),
          rep((ncol(fx_init)-1) * nrow(fx_init) + 1:nrow(fx_init), max(bf.idx5$periods - 2004)) )) %>%
  mutate(year = rep(1960:max(bf.idx5$periods), each = nrow(fx_init))) %>%
  select(-2) %>%
  pivot_wider(names_from = year, values_from = value) %>%
  as.matrix()

gx_init <- burkina.faso.females$migration.proportions
gx_init <- cbind(gx_init, `2005` = gx_init[ , "2000"], `2010` = gx_init[ , "2000"], `2015` = gx_init[ , "2000"])

log_basepop_mean <- as.vector(log(basepop_init))
logit_sx_mean <- as.vector(qlogis(sx_init))
log_fx_mean <- as.vector(log(fx_init[,-1]))
gx_mean <- as.vector(gx_init)

##prior means for LQ h from IGME 5q0
igme.5q0<-read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/IGME BF 5q0.csv")

igme.5q0.df <- reshape2::melt(igme.5q0) %>% mutate(year = as.numeric(gsub("X|\\.5","",variable)), child.mort = value/1000) %>%
  select(Sex,year,child.mort)


igme.5q0.5<-reshape2::melt(igme.5q0) %>% mutate(year = as.numeric(gsub("X","",variable)), child.mort = value/1000) %>%
  select(Sex,year,child.mort) %>%
  mutate(year5 = 5 * floor(year/5)) %>% group_by(Sex, year5) %>%
  summarise_at(vars(child.mort),mean)


data <- list(log_basepop_mean = log_basepop_mean,
             log_fx_mean = log_fx_mean,
             gx_mean = gx_mean,
             srb = rep(1.05, bf.idx5$n_periods),
             interval = 5,
             n_periods = bf.idx5$n_periods,
             fx_idx = 4L,
             n_fx = 7L,
             census_log_pop = log(burkina.faso.females$census.pop.counts),
             census_year_idx = c(4L, 6L, 8L, 10L),
            
             baseline_mx_DX = baseline.DX,
             mx_params_DX = year.DX,
             tp_DX = tips.DX,
             penal_mx = as(crossprod(diff(diag(bf.idx5$n_periods))),"sparseMatrix"),
             penal_mx_0 = as(tcrossprod(c(1,rep(0,bf.idx5$n_periods-1))),"sparseMatrix"),
             penal_tp = as(crossprod(diff(diag(15))),"sparseMatrix"),
             null_penal_tp = as(exp(15)*tcrossprod(c(0,1,1,1,rep(0,11))),"sparseMatrix"),
             penal_tp_0 = as(tcrossprod(c(1,rep(0,14))),"sparseMatrix"),
             dx = bf5.f$event,
             Ex = bf5.f$pyears,
            
             mx_params_mean = rep(0,bf.idx5$n_periods),
             log_baseline_mx = log_baseline_mx,
            
             open_idx = bf.idx5$n_ages, #open age starting at 80
             LQ_baseline_mx_DX = LQ.baseline.DX.ax,
             h_DX = LQ.baseline.DX.bx,
             h2_DX = LQ.baseline.DX.cx,
             k_DX = LQ.baseline.DX.vx,
             LQ_baseline = as.matrix(LQcoef.f[,1:4])
)

par <- list(log_tau2_logpop = 0,
            log_tau2_sx = 0,
            log_tau2_fx = 0,
            log_tau2_gx = 0,
            log_basepop = log_basepop_mean,
            logit_sx = logit_sx_mean,
            log_fx = log_fx_mean,
            gx = gx_mean,
            log_lambda_tp = 0,
            log_lambda_tp_0_inflated_sd = 0,
            log_lambda_mx = 0,
            log_lambda_mx_0_inflated_sd = 0,
            log_dispersion=2,
            
            mx_params = rep(0,bf.idx5$n_periods),
            tp_params = rep(0,15),
            
            log_lambda_h = 0,
            log_lambda_h_0_inflated_sd = 0,
            log_lambda_k = 0,
            log_lambda_k_0_inflated_sd = 0, 
            h_params = rep(log(0.2),bf.idx5$n_periods),
            k_params = rep(log(0.6),bf.idx5$n_periods),
            
            log_tau2_hiv_mort = 0,
            log_hiv_mort = rep(0,bf.idx5$n_periods))

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
                    map = list()
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
                                 iter.max = max_iter)),
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


obj <- make_tmb_obj(data, par, inner_verbose = TRUE , model = "ccmpp_vr_tmb", 
                    random =  c("log_basepop", "mx_params", "tp_params", "log_fx", "gx"),
                    DLL="trial2")

input <- list(data = data, par_init = par, model = "ccmpp_vr_tmb")
fit <- fit_tmb(input,inner_verbose=TRUE, random = c("log_basepop", "mx_params", "tp_params", "log_fx", "gx"), DLL="trial2")
est.mx <- matrix(fit$mode$mx, bf.idx5$n_ages, bf.idx5$n_periods)
colnames(est.mx) <- bf.idx5$periods

input.LQ <- list(data = data, par_init = par, model = "ccmpp_vr_tmb")
fit.LQ <- fit_tmb(input.LQ,inner_verbose=TRUE, random = c("log_basepop", "h_params", "k_params", "tp_params", "log_fx", "gx"), DLL="ccmpp_LQuad_singlesex")


data.bothsexes <- list(log_basepop_mean_f = log_basepop_mean, log_basepop_mean_m = as.vector(log(basepop_init.m)),
             log_fx_mean = log_fx_mean,
             gx_mean_f = gx_mean, gx_mean_m = gx_mean, 
             srb = rep(1.05, bf.idx5$n_periods),
             interval = 5,
             n_periods = bf.idx5$n_periods,
             fx_idx = 4L,
             n_fx = 7L,
             census_log_pop_f = as.matrix(log(bf.pop.aggr.f[,-1])), census_log_pop_m = as.matrix(log(bf.pop.aggr.m[,-1])),
             census_year_idx = 2:bf.idx5$n_periods,
             
             open_idx = bf.idx5$n_ages, #open age starting at 80
             LQ_baseline_mx_DX_f = LQ.baseline.DX.ax, LQ_baseline_mx_DX_m = LQ.baseline.DX.ax.m,
             h_DX_f = LQ.baseline.DX.bx, h_DX_m = LQ.baseline.DX.bx.m,
             h2_DX_f = LQ.baseline.DX.cx, h2_DX_m = LQ.baseline.DX.cx.m,
             k_DX_f = LQ.baseline.DX.vx, k_DX_m = LQ.baseline.DX.vx.m,
             tp_DX_f = tips.DX, tp_DX_m = tips.DX.m,
             
             penal_mx = as(crossprod(diff(diag(bf.idx5$n_periods))),"sparseMatrix"),
             penal_mx_0 = as(tcrossprod(c(1,rep(0,bf.idx5$n_periods-1))),"sparseMatrix"),
             penal_tp = as(crossprod(diff(diag(15))),"sparseMatrix"),
             null_penal_tp = as(exp(15)*tcrossprod(c(0,1,1,1,rep(0,11))),"sparseMatrix"),
             penal_tp_0 = as(tcrossprod(c(1,rep(0,14))),"sparseMatrix"),
             
             LQ_baseline_f = as.matrix(LQcoef.f[,1:4]), LQ_baseline_m = as.matrix(LQcoef.m[,1:4]),
             df = bf5.f$event, dm = bf5.m$event,
             Ef = bf5.f$pyears, Em = bf5.m$pyears,
             
             rho_mat_ones = as(cbind(rep(0,bf.idx5$n_periods - 1), diag(bf.idx5$n_periods - 1)), "sparseMatrix"),
             rho_mat = as(cbind(diag(bf.idx5$n_periods - 1), rep(0,bf.idx5$n_periods - 1)), "sparseMatrix")
)



par.bothsexes <- list(log_tau2_logpop_f = 0, log_tau2_logpop_m = 0,
            log_tau2_fx = 0,
            log_tau2_gx_f = 0, log_tau2_gx_m = 0,
            log_lambda_h_f = 0, log_lambda_h_m = 0,
            log_lambda_h_0_inflated_sd_f = 0, log_lambda_h_0_inflated_sd_m = 0,
            log_lambda_k_f = 0, log_lambda_k_m = 0,
            log_lambda_k_0_inflated_sd_f = 0, log_lambda_k_0_inflated_sd_m = 0, 
            log_lambda_tp = 0,
            log_lambda_tp_0_inflated_sd = 0,
            
            log_dispersion_f = 2, log_dispersion_m = 2,
            
            log_basepop_f = log_basepop_mean, log_basepop_m = as.vector(log(basepop_init.m)),
            log_fx = log_fx_mean,
            gx_f = gx_mean, gx_m = gx_mean,
            
            h_params_f = rep(log(0.2),bf.idx5$n_periods), h_params_m = rep(log(0.2),bf.idx5$n_periods),
            k_params_f = rep(0,bf.idx5$n_periods), k_params_m = rep(0,bf.idx5$n_periods),
            
            tp_params = rep(0,15),
            
            logit_rho_h_f = 0.5, logit_rho_h_m = 0.5,
            logit_rho_k_f = 0.5, logit_rho_k_m = 0.5,
            h_constant_f = log(0.1), h_constant_m = log(0.1),
            k_constant_f = 4, k_constant_m = 4
)


input.LQ.both <- list(data = data.bothsexes, par_init = par.bothsexes, model = "ccmpp_vr_tmb")
fit.LQ.both <- fit_tmb(input.LQ.both,inner_verbose=TRUE, random = c("log_basepop_f", "log_basepop_m",
                                                          "log_fx",
                                                          "gx_f", "gx_m",
                                                          "h_params_f", "h_params_m",
                                                          "k_params_f", "k_params_m",
                                                          "tp_params"), DLL="ccmpp_LQuad_bothsexes")

fit.LQ.both.AR <- fit_tmb(input.LQ.both,inner_verbose=TRUE, map = list(k_constant_f = factor(NA), k_constant_m = factor(NA)), random = c("log_basepop_f", "log_basepop_m",
                                                                    "log_fx",
                                                                    "gx_f", "gx_m",
                                                                    "h_params_f", "h_params_m",
                                                                    "k_params_f", "k_params_m",
                                                                    "h_constant_f", "h_constant_m",
                                                                    "tp_params"), DLL="ccmpp_LQuad_bothsexes_AR")


data.vec <- list(log_basepop_mean_f = log_basepop_mean, log_basepop_mean_m = as.vector(log(basepop_init.m)),
                       log_fx_mean = log_fx_mean,
                       srb = rep(1.05, bf.idx5$n_periods),
                       interval = bf.idx5$interval,
                       n_periods = bf.idx5$n_periods,
                       fx_idx = 4L,
                       n_fx = 7L,
                       census_log_pop_f = as.matrix(log(bf.pop.aggr.f[,-1])), census_log_pop_m = as.matrix(log(bf.pop.aggr.m[,-1])),
                       census_year_idx = match(colnames(bf.pop.aggr.f)[-1], bf.idx5$periods),
                       
                       open_idx = bf.idx5$n_ages, #open age starting at 80
                       LQ_baseline_mx_DX_f = LQ.baseline.DX.ax, LQ_baseline_mx_DX_m = LQ.baseline.DX.ax.m,
                       h_DX_f = LQ.baseline.DX.bx, h_DX_m = LQ.baseline.DX.bx.m,
                       h2_DX_f = LQ.baseline.DX.cx, h2_DX_m = LQ.baseline.DX.cx.m,
                       k_DX_f = LQ.baseline.DX.vx, k_DX_m = LQ.baseline.DX.vx.m,
                       tp_DX_f = tips.DX, tp_DX_m = tips.DX.m,
                       
                       penal_tp = as(crossprod(diff(diag(15))),"sparseMatrix"),
                       null_penal_tp = as(exp(15)*tcrossprod(c(0,1,1,1,rep(0,11))),"sparseMatrix"),
                       penal_tp_0 = as(tcrossprod(c(1,rep(0,14))),"sparseMatrix"),
                       
                       LQ_baseline_f = as.matrix(LQcoef.f[,1:4]), LQ_baseline_m = as.matrix(LQcoef.m[,1:4]),
                       df = bf5.f$event, dm = bf5.m$event,
                       Ef = bf5.f$pyears, Em = bf5.m$pyears,
                 
                       h_mean_f = igme.5q0.df %>% filter(Sex=="Female", year %in% bf.idx5$periods) %>% .$child.mort %>% log(),
                       h_mean_m = igme.5q0.df %>% filter(Sex=="Male", year %in% bf.idx5$periods) %>% .$child.mort %>% log()
                 )



par.vec <- list(log_tau2_logpop_f = 5, log_tau2_logpop_m = 5,
                      log_tau2_fx = 5,
                      log_tau2_gx_f = 1.8, log_tau2_gx_m = 1.8,
                      log_marginal_var_h = rep(log(0.5^2), 2),
                      log_marginal_var_k = rep(log(0.5^2), 2),
                      log_lambda_tp = 1.7,
                      log_lambda_tp_0_inflated_sd = 0.3,
                      
                      log_dispersion = rep(1.3, 2),
                      
                      log_basepop_f = log_basepop_mean, log_basepop_m = as.vector(log(basepop_init.m)),
                      log_fx = log_fx_mean,
                      gx_f = rep(0, bf.idx5$n_ages * bf.idx5$n_periods), gx_m = rep(0, bf.idx5$n_ages * bf.idx5$n_periods),
                
                      h_params_f = igme.5q0.df %>% filter(Sex=="Female", year %in% bf.idx5$periods) %>% .$child.mort %>% log(), 
                      h_params_m = igme.5q0.df %>% filter(Sex=="Male", year %in% bf.idx5$periods) %>% .$child.mort %>% log(),
                      k_params_f = rep(0,bf.idx5$n_periods), k_params_m = rep(0,bf.idx5$n_periods),
                      
                      tp_params = rep(0,15),
                      
                      logit_rho_h = rep(2, 2),
                      logit_rho_k = rep(0, 2),
                      h_constant_f = igme.5q0.df %>% filter(Sex=="Female", year %in% bf.idx5$periods) %>% .$child.mort %>% log(),
                      h_constant_m = igme.5q0.df %>% filter(Sex=="Male", year %in% bf.idx5$periods) %>% .$child.mort %>% log(),
                      k_constant_f = 0, k_constant_m = 0,
                
                      logit_rho_g_x = rep(0, 2),
                      logit_rho_g_t = rep(4, 2)
                )


input.LQ.both.vec <- list(data = data.vec, par_init = par.vec, model = "ccmpp_vr_tmb")
fit.LQ.both.vec <- fit_tmb(input.LQ.both.vec,inner_verbose=TRUE, random = c("log_basepop_f", "log_basepop_m",
                                                                    "log_fx",
                                                                    "gx_f", "gx_m",
                                                                    "h_params_f", "h_params_m",
                                                                    "k_params_f", "k_params_m",
                                                                    "h_constant_f", "h_constant_m",
                                                                    "k_constant_f", "k_constant_m",
                                                                    "tp_params"), DLL="ccmpp_LQuad_bothsexes_AR_vec",
                           map = list(h_constant_f = factor(rep(1,bf.idx5$n_periods)), h_constant_m = factor(rep(1,bf.idx5$n_periods))))

fit.LQ.both.vec.fixk <- fit_tmb(input.LQ.both.vec,inner_verbose=TRUE, random = c("log_basepop_f", "log_basepop_m",
                                                                            "log_fx",
                                                                            "gx_f", "gx_m",
                                                                            "h_params_f", "h_params_m",
                                                                            "k_params_f", "k_params_m",
                                                                            "h_constant_f", "h_constant_m",
                                                                            "k_constant_f", "k_constant_m",
                                                                            "tp_params"), DLL="ccmpp_LQuad_bothsexes_AR_vec",
                           map = list(h_constant_f = factor(rep(1,bf.idx5$n_periods)), h_constant_m = factor(rep(1,bf.idx5$n_periods)),
                                      k_constant_f = factor(NA), k_constant_m = factor(NA)))


fit.LQ.both.vec.igme.prior <- fit_tmb(input.LQ.both.vec,inner_verbose=TRUE, random = c("log_basepop_f", "log_basepop_m",
                                                                                 "log_fx",
                                                                                 "gx_f", "gx_m",
                                                                                 "h_params_f", "h_params_m",
                                                                                 "k_params_f", "k_params_m",
                                                                                 "k_constant_f", "k_constant_m",
                                                                                 "tp_params"), DLL="ccmpp_LQuad_bothsexes_AR_vec",
                                map = list(h_constant_f = factor(rep(NA,bf.idx5$n_periods)), h_constant_m = factor(rep(NA,bf.idx5$n_periods))))

fit.LQ.both.weighted.hmean <- fit_tmb(input.LQ.both.vec,inner_verbose=TRUE, random = c("log_basepop_f", "log_basepop_m",
                                                                                 "log_fx",
                                                                                 "gx_f", "gx_m",
                                                                                 "h_params_f", "h_params_m",
                                                                                 "k_params_f", "k_params_m",
                                                                                 "k_constant_f", "k_constant_m",
                                                                                 "tp_params"), DLL="ccmpp_LQuad_bothsexes_weighted_hmean")

###plots####
mf5 <- projection_model_frames(bf.idx5)
haha <-fit.LQ.both.weighted.hmean

#population
mf5$mf_population %>% arrange(sex,period,age) %>%
  mutate(population = c(haha$mode$population_f,haha$mode$population_m),
         source = "Fitted") %>%
  bind_rows(
    bind_rows(
      cbind.data.frame(age = bf.idx5$ages, bf.pop.aggr.f) %>% pivot_longer(!age) %>%
        select(age = age, period = name, population = value) %>%
        mutate(source = "Data", period = as.numeric(period), sex="female"),
      
      cbind.data.frame(age = bf.idx5$ages, bf.pop.aggr.m) %>% pivot_longer(!age) %>%
        select(age = age, period = name, population = value) %>%
        mutate(source = "Data", period = as.numeric(period), sex="male")
    )
  ) %>%
  mutate(source = fct_inorder(source)) %>%
  count(period, source, sex, wt = population) %>%
  ggplot(aes(period, n, linetype = source, color=sex)) +
  geom_line() +
  scale_y_continuous("Population (millions)",
                     labels = scales::label_number(scale = 1e-6)) +
  expand_limits(y = 0) +
  ggtitle("BF population 1960-2020")

#mx
persp(x = seq(0,110,by=5), y = bf.idx5$periods, z = log(haha$mode$mx_mat_f), theta = -10, phi = 25)
persp(x = seq(0,110,by=5), y = bf.idx5$periods, z = log(haha$mode$mx_mat_m), theta = -10, phi = 25)

bf.ccmpp.mx <- crossing(age=seq(0,110,by=5),year=bf.idx5$periods, sex=c("male","female")) %>% arrange(sex,year) %>% 
  mutate(mx = c(haha$mode$mx_mat_f, haha$mode$mx_mat_m))

bf.ccmpp.mx.rw <- crossing(age=seq(0,110,by=5),year=bf.idx5$periods, sex=c("male","female")) %>% arrange(sex,year) %>% 
  mutate(mx = c(fit.LQ.both$mode$mx_mat_f, fit.LQ.both$mode$mx_mat_m))

plot.df <- bf.ccmpp.mx %>% mutate(model = "AR") %>% 
  bind_rows(mutate(bf.ccmpp.mx.rw, model = "RW"))

ggsave(filename="compare RW males year.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,
       plot={
         ggplot(plot.df) +
    geom_line(data = subset(plot.df, model=="AR" & sex=="male"), aes(x = age, y = mx, group = year, color = as.factor(year), linetype = model), lwd=1.2) +
    geom_line(data = subset(plot.df, model=="RW" & sex=="male"), aes(x = age, y = mx, group = year, color = as.factor(year), linetype = model), lwd=1.2) +
    scale_y_continuous(trans="log") + ylab(bquote(""[5]*m[x])) + ggtitle(bquote("AR  vs  RW  male "~""[5]*m[x]~"(on log scale)")) +
    theme(plot.title = element_text(hjust = 0.5, size = 30), title = element_text(size = 20), axis.text = element_text(size = 15))
       })

ggsave(filename="compare RW females year.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,
       plot={
         ggplot(plot.df) +
    geom_line(data = subset(plot.df, model=="AR" & sex=="female"), aes(x = age, y = mx, group = year, color = as.factor(year), linetype = model), lwd=1.2) +
    geom_line(data = subset(plot.df, model=="RW" & sex=="female"), aes(x = age, y = mx, group = year, color = as.factor(year), linetype = model), lwd=1.2) +
    scale_y_continuous(trans="log") + ylab(bquote(""[5]*m[x])) + ggtitle(bquote("AR  vs  RW  female "~""[5]*m[x]~"(on log scale)")) +
    theme(plot.title = element_text(hjust = 0.5, size = 30), title = element_text(size = 20), axis.text = element_text(size = 15))
       })
#45q15
q4515.f<-apply(haha$mode$mx_mat_f[4:12,],2,function(x){1-prod((1-2.5*x)/(1+2.5*x))})
q4515.m<-apply(haha$mode$mx_mat_m[4:12,],2,function(x){1-prod((1-2.5*x)/(1+2.5*x))})
ggsave(filename="CCMPP q4515.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,
       plot={
         ggplot() + geom_line(aes(x=bf.idx5$periods, y=q4515.f,color="female")) + geom_line(aes(x=bf.idx5$periods, y=q4515.m,color="male")) +
           xlab("Year") + ylab(bquote(""[45]*q[15])) + ggtitle(bquote("CCMPP-LQ Model "[45]*q[15]~"(assuming UDD)")) + 
           scale_color_manual(values = c("red", "blue"), name = "sex", guide = guide_legend(reverse=TRUE)) + 
           theme(plot.title = element_text(hjust = 0.5, size = 30), title = element_text(size = 20), axis.text = element_text(size = 15))
       })


##compared to only RW on mx and WPP estimates
wpp.bf.q4515<-read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/WPP_BF_45q15.csv")
q4515.f.rw<-apply(fit.LQ.both$mode$mx_mat_f[4:12,],2,function(x){1-prod((1-2.5*x)/(1+2.5*x))})
q4515.m.rw<-apply(fit.LQ.both$mode$mx_mat_m[4:12,],2,function(x){1-prod((1-2.5*x)/(1+2.5*x))})

ggsave(filename="CCMPP q4515 compare WPP and RW.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,
       plot={
         ggplot() + geom_line(aes(x=bf.idx5$periods, y=q4515.f,color="female",linetype="AR on mx + AR(1):AR(1) on gx"),lwd=1.2) + geom_line(aes(x=bf.idx5$periods, y=q4515.m,color="male",linetype="AR on mx + AR(1):AR(1) on gx"),lwd=1.2) +
  geom_line(aes(x=bf.idx5$periods, y=q4515.f.rw, color="female", linetype="RW on mx (previous fit)"),lwd=1.2) + geom_line(aes(x=bf.idx5$periods, y=q4515.m.rw,color="male", linetype="RW on mx (previous fit)"),lwd=1.2) + 
  geom_line(aes(x=bf.idx5$periods, y=as.numeric(wpp.bf.q4515[1,6:17])/1000, color="female", linetype="WPP estimates"),lwd=1.2) + geom_line(aes(x=bf.idx5$periods, y=as.numeric(wpp.bf.q4515[2,6:17])/1000, color="male", linetype="WPP estimates"),lwd=1.2) + 
  xlab("Year") + ylab(bquote(""[45]*q[15])) + ggtitle(bquote("CCMPP-LQ Model "[45]*q[15]~"(assuming UDD)")) + 
  scale_linetype_manual(values=c("solid","dotted","dashed")) +
  scale_color_manual(values = c("red", "blue"), name = "sex", guide = guide_legend(reverse=TRUE)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 30), title = element_text(size = 20), axis.text = element_text(size = 15))
       })

#empirical sx
emp.sx.m <- emp.sx.f <- matrix(0, nrow(bf.pop.aggr.m)-1, ncol(bf.pop.aggr.m)-1)

for(j in 1:(ncol(bf.pop.aggr.m)-1)) {
  for (i in 1:(nrow(bf.pop.aggr.m)-1)) {
    emp.sx.m[i,j] <- as.numeric(bf.pop.aggr.m[i+1,j+1] / bf.pop.aggr.m[i,j])
    emp.sx.f[i,j] <- as.numeric(bf.pop.aggr.f[i+1,j+1] / bf.pop.aggr.f[i,j])
  }
}

rownames(emp.sx.m) <- rownames(emp.sx.f) <- seq(0,75,by=5)
colnames(emp.sx.m) <- colnames(emp.sx.f) <- bf.idx5$periods[-bf.idx5$n_periods]

s4515.df <- reshape2::melt(t(1-apply(emp.sx.m[4:12,],2,prod))) %>% dplyr::select(year = Var2, s4515 = value) %>% mutate(source="empirical", sex="male") %>%
  bind_rows(reshape2::melt(t(1-apply(emp.sx.f[4:12,],2,prod))) %>% dplyr::select(year = Var2, s4515 = value) %>% mutate(source="empirical", sex="female")) %>%
  bind_rows(
    crossing(year = bf.idx5$periods, sex=c("female","male")) %>% arrange(sex,year) %>%
      mutate(s4515 = c(1-apply(haha$mode$sx_mat_f[5:13,],2,prod),1-apply(haha$mode$sx_mat_m[5:13,],2,prod)), source = "LQ"))

ggsave(filename="CCMPP 1-s4515.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,
       plot={
         ggplot(s4515.df) + geom_line(aes(x = year, y = s4515, color = sex, linetype = source), lwd = 1.2) +
           scale_color_manual(values = c("red", "blue")) + scale_linetype(labels = c("Empirical (deaths + migration)", "LogQuad Estimates")) +
           ylab(bquote(1 - ""[45]*s[15])) + ggtitle(bquote(1-""[45]*s[15]~"Empirical vs LQ Modelled")) + 
           theme(plot.title = element_text(hjust = 0.5, size = 30), title = element_text(size = 20), axis.text = element_text(size = 15))
       })

#migration 
gx_m <- matrix(haha$mode$gx_m, bf.idx5$n_ages, bf.idx5$n_periods)
gx_f <- matrix(haha$mode$gx_f, bf.idx5$n_ages, bf.idx5$n_periods)
rownames(gx_m) <- rownames(gx_f) <- bf.idx5$ages
colnames(gx_m) <- colnames(gx_f) <- bf.idx5$periods

reshape2::melt(gx_m) %>% ggplot +
  geom_line(aes(x = Var2, y = value, color = as.factor(Var1), group = as.factor(Var1)))

reshape2::melt(gx_f) %>% ggplot +
  geom_line(aes(x = Var2, y = value, color = as.factor(Var1), group = as.factor(Var1)))

g4515.df <- mf5$mf_migrations %>% arrange(sex,period,age) %>%
  mutate(migration = c(haha$mode$migrations_f,haha$mode$migrations_m)) %>%
  left_join(mf5$mf_population %>% arrange(sex,period,age) %>%
              mutate(population = c(haha$mode$population_f,haha$mode$population_m)), by=c("period","sex","age")) %>%
  filter(age %in% 15:55) %>% group_by(sex,period) %>%
  summarise_at(vars(migration, population), sum) %>% mutate(g4515 = migration / population)

gx.df <-reshape2::melt(gx_f) %>% mutate(sex="female") %>%
  bind_rows(reshape2::melt(gx_m) %>% mutate(sex="male")) %>% 
  setNames(c("age","year","gx","sex"))

ggsave(filename="gx 15 to 60.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,
       plot = {
         ggplot(filter(gx.df, age %in% 15:55)) + 
           geom_line(data = filter(gx.df, age %in% 15:55, sex=="female"), aes(x = year, y = gx, color = age, group = age), lwd = 1.2) +
           scale_color_gradientn(colors=c("gray40","red","magenta")) + new_scale_color() +
           geom_line(data = filter(gx.df, age %in% 15:55, sex=="male"), aes(x = year, y = gx, color = age, group = age), lwd = 1.2, linetype = 2) +
           scale_color_gradientn(colors=c("gray40","blue","cyan")) +
           ggtitle(bquote("Estimated"~""[5]*g[x])) + geom_hline(yintercept = 0, lwd=1.2) + 
           theme(plot.title = element_text(hjust = 0.5, size = 30), title = element_text(size = 20), axis.text = element_text(size = 15))
         
       })

#tips
tp<-split(unname(fit.LQ.both$par.full),names(fit.LQ.both$par.full))$tp_params
plot(x=exp(tp), y=0:14, type="l")

###stan####
#' # PopReconstruct VR
#'


deaths_obs <- swe_lt5_pop %>%
  filter(year != max(year)) %>%
  .$deaths %>%
  matrix(nrow = idx5$n_ages)

births_obs <- matrix(swe_fert5$births, nrow = idx5$n_fx)


data_vr <- list(log_basepop_mean = log_basepop_mean,
                logit_sx_mean = logit_sx_mean,
                log_fx_mean = log_fx_mean,
                gx_mean = gx_mean,
                srb = srb,
                interval = idx5$interval,
                n_periods = idx5$n_periods,
                fx_idx = idx5$fx_idx,
                n_fx = idx5$n_fx,
                census_log_pop = census_log_pop[ , census_year_idx, drop = FALSE],
                census_year_idx = census_year_idx,
                deaths_obs = deaths_obs,
                births_obs = births_obs)

par <- list(log_tau2_logpop = 0,
            log_tau2_sx = 0,
            log_tau2_fx = 0,
            log_tau2_gx = 0,
            log_basepop = log_basepop_mean,
            logit_sx = logit_sx_mean,
            log_fx = log_fx_mean,
            gx = gx_mean)

obj_vr <- make_tmb_obj(data_vr,
                       par,
                       model = "ccmpp_vr_tmb",
                       inner_verbose = TRUE)
init_sim <- obj_vr$report()

obj_vr$fn()
obj_vr$gr()

input_vr <- list(data = data_vr, par_init = par, model = "ccmpp_vr_tmb")

fit_vr <- fit_tmb(input_vr, inner_verbose = TRUE)
fit_vr[1:6]

fit_vr <- sample_tmb(fit_vr)


if(!file.exists("swe5yr_fitstan_vr.rds")) {
  fitstan_vr <- fit_tmbstan(input_vr, chains = 4, iter = 1000,
                            use_inits = TRUE, refresh = 100)
  saveRDS(fitstan_vr, "swe5yr_fitstan_vr.rds")
} else {
  fitstan_vr <- readRDS("swe5yr_fitstan_vr.rds")
}


fitstan_vr <- sample_tmbstan(fitstan_vr, TRUE)

fitstan_vr$sample$migrations %>%
  `colnames<-`(seq_len(ncol(.))) %>%
  as_tibble() %>%
  bind_cols(mf5$mf_migrations, .) %>%
  gather(sample, value, `1`:last_col()) %>%
  group_by(period, sample) %>%
  summarise(value = sum(value)) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975)) %>%
  print(n = Inf)

fitstan$sample$migrations %>%
  `colnames<-`(seq_len(ncol(.))) %>%
  as_tibble() %>%
  bind_cols(mf5$mf_migrations, .) %>%
  gather(sample, value, `1`:last_col()) %>%
  group_by(period, sample) %>%
  summarise(value = sum(value)) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975)) %>%
  print(n = Inf)


fitstan_vr$sample$population %>%
  `colnames<-`(seq_len(ncol(.))) %>%
  as_tibble() %>%
  bind_cols(mf5$mf_population, .) %>%
  gather(sample, value, `1`:last_col()) %>%
  group_by(period, sample) %>%
  summarise(value = sum(value)) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

fitstan$sample$population %>%
  `colnames<-`(seq_len(ncol(.))) %>%
  as_tibble() %>%
  bind_cols(mf5$mf_population, .) %>%
  gather(sample, value, `1`:last_col()) %>%
  group_by(period, sample) %>%
  summarise(value = sum(value)) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

fitstan_vr$sample$period_deaths %>%
  `colnames<-`(seq_len(ncol(.))) %>%
  as_tibble() %>%
  bind_cols(mf5$mf_deaths, .) %>%
  gather(sample, value, `1`:last_col()) %>%
  group_by(period, sample) %>%
  summarise(value = sum(value)) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

fitstan$sample$period_deaths %>%
  `colnames<-`(seq_len(ncol(.))) %>%
  as_tibble() %>%
  bind_cols(mf5$mf_deaths, .) %>%
  gather(sample, value, `1`:last_col()) %>%
  group_by(period, sample) %>%
  summarise(value = sum(value)) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))


fitstan_vr$sample$births %>%
  `colnames<-`(seq_len(ncol(.))) %>%
  as_tibble() %>%
  bind_cols(mf5$mf_births, .) %>%
  gather(sample, value, `1`:last_col()) %>%
  group_by(period, sample) %>%
  summarise(value = sum(value)) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

fitstan$sample$births %>%
  `colnames<-`(seq_len(ncol(.))) %>%
  as_tibble() %>%
  bind_cols(mf5$mf_births, .) %>%
  gather(sample, value, `1`:last_col()) %>%
  group_by(period, sample) %>%
  summarise(value = sum(value)) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))


fitstan$sample$period_deaths %>%
  `colnames<-`(seq_len(ncol(.))) %>%
  as_tibble() %>%
  bind_cols(mf5$mf_period_deaths, .) %>%
  gather(sample, value, `1`:last_col()) %>%
  group_by(period, sample) %>%
  summarise(value = sum(value)) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

fitstan$sample$period_deaths %>%
  `colnames<-`(seq_len(ncol(.))) %>%
  as_tibble()