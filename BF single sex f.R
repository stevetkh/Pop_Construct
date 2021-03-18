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
library(DDSQLtools)
library(ggforce)

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_AR_f_no0.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_AR_f_no0"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_MVN_f_no0.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_MVN_f_no0"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_weighted_hmean_f_no0.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_weighted_hmean_f_no0"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_ARIMAk_f_no0.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_ARIMAk_f_no0"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_VAR_f_no0.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/ccmpp_LQuad_bothsexes_VAR_f_no0"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/trial.cpp")

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

###Burkina Faso####
library(MortCast)
load("~/cohort smooth 1900-2017.RData")

###UNPD Census Data
bf_census <- get_recorddata(dataProcessTypeIds = "Census",
                            indicatorTypeIds = "Population by age and sex",
                            locIds = "Burkina Faso",
                            locAreaTypeIds = "Whole Area",
                            subGroupIds = "Total or All groups",
                            isComplete = "Abridged")
open.age <- 85

bf_census_f <- bf_census %>% filter(SexName=="Female",!AgeLabel %in% c("Total", "Unknown")) %>%
  select(AgeLabel, TimeMid, TimeLabel, StatisticalConceptName, DataValue) %>%
  group_by(AgeLabel, TimeMid) %>%
  filter(StatisticalConceptName == "De-jure") %>%
  arrange(TimeMid, AgeLabel) %>%
  summarise_at(vars(DataValue),first) %>% 
  arrange(TimeMid, AgeLabel) %>% ungroup() %>%
  pivot_wider(names_from = TimeMid, values_from = DataValue) %>%
  arrange(AgeLabel) %>% select(c(1,2,4,5,8)) %>%
  mutate(age.aggr = ifelse(AgeLabel %in% c("0","0-4","1-4"),"0-4", 
                           ifelse(AgeLabel %in% c(sprintf("%i-%i", seq(open.age, 90, by=5), seq(open.age+4, 94, by=5)), "80+", "85+", "95+", "98+"), paste0(open.age,"+"), AgeLabel))) %>%
  group_by(age.aggr) %>%
  summarise_at(vars(-AgeLabel),sum,na.rm=T) %>%
  mutate(age = as.numeric(regmatches(age.aggr, regexpr("\\d+",age.aggr))), .before = 1) %>%
  arrange(age) %>% select(-age.aggr) %>%
  setNames(c("age",1975,1985,1996,2006))

##WPP Pop Estimates
bf.pop <- read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/Burkina Faso Pop.csv")
colnames(bf.pop) <- c("Name", "Sex", "Year", sprintf("%d-%d", seq(0, 95, by=5), seq(4, 99, by=5)), "100+")

bf.pop.f <- bf.pop %>% filter(Sex=="Female") %>% select(-c("Name", "Sex")) %>% apply(2,function(x)as.numeric(gsub(" ","",x))) %>% as_tibble() %>%
  pivot_longer(!Year,names_to="Age",values_to="counts") %>% pivot_wider(names_from=Year,values_from="counts")
bf.pop.f[,-1] <- bf.pop.f[,-1]*1000

##Pop data from Mark's package
data(burkina_faso_females)

##CHANGE
#n_ages = nrow(burkina.faso.females$census.pop.counts)
n_ages = open.age / 5 + 1

bf.idx5<-projection_indices(period_start = 1960,
                            period_end = 2015,
                            interval = 5,
                            n_ages = n_ages,
                            fx_idx = 4L,
                            n_fx = 7L,
                            n_sexes = 1)

bf.pop.aggr.f <- bf.pop.f[1:bf.idx5$n_ages,]
bf.pop.aggr.f[bf.idx5$n_ages,-1] <- t(bf.pop.f %>% slice(-(1:(bf.idx5$n_ages-1))) %>% select(-1) %>% colSums())
#bf.pop.aggr.f[bf.idx5$n_ages,1] <- "80+"
#bf.pop.aggr.f <- bf.pop.aggr.f %>% select(sprintf("%d",bf.idx5$periods))

##Log Quad
LQcoef.f <- LQcoef %>% filter(sex=="Female", !age%in%c("0","1-4")) %>% select(ax:vx)
LQcoef.f$age5 <- seq(5,110,by=5)

##DHS data
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
bf5$age5 <- as.numeric(gsub("-\\d+","",bf5$agegr))
bf5$tips <- factor(bf5$tips,levels=0:14)
bf5.f <- bf5 %>% filter(mm1=="female") %>% arrange(period5,tips,age5)
bf5.f.no0 <- filter(bf5.f, age5  >= 15 & age5 <=45)
bf5.f.0 <- filter(bf5.f, age5 == 0)

tips.DX <- as(model.matrix(event~as.factor(tips)-1,data=bf5.f.no0),"sparseMatrix")
tips.DX.0 <- as(model.matrix(event~as.factor(tips)-1,data=bf5.f.0),"sparseMatrix")

LQ.baseline.DX <- model.matrix(event~factor(age5, levels = seq(5, 110, by = 5)) - 1, data=bf5.f.no0) %*% as.matrix(LQcoef.f[,1:4])
LQ.baseline.DX.ax <- LQ.baseline.DX[,1]
LQ.baseline.DX.bx <- data.frame(period5 = as.factor(bf.idx5$periods)) %>% 
  bind_rows(bind_cols(n = 1 : nrow(bf5.f.no0), ax=LQ.baseline.DX[,2], period5=bf5.f.no0$period5)) %>% 
  pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
  slice(-1) %>% select(-n) %>%
  as.matrix() %>% as("sparseMatrix")
LQ.baseline.DX.cx <- data.frame(period5 = as.factor(bf.idx5$periods)) %>% 
  bind_rows(bind_cols(n = 1 : nrow(bf5.f.no0), ax=LQ.baseline.DX[,3], period5=bf5.f.no0$period5)) %>% 
  pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
  slice(-1) %>% select(-n) %>%
  as.matrix() %>% as("sparseMatrix")
LQ.baseline.DX.vx <- data.frame(period5 = as.factor(bf.idx5$periods)) %>% 
  bind_rows(bind_cols(n = 1 : nrow(bf5.f.no0), ax=LQ.baseline.DX[,4], period5=bf5.f.no0$period5)) %>% 
  pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
  slice(-1) %>% select(-n) %>%
  as.matrix() %>% as("sparseMatrix")
LQ.baseline.0 <- as(model.matrix(event~as.factor(period5)-1,data=bf5.f.0),"sparseMatrix")
  

##DHS data cohort smoothed
bf5.smooth <- aggr.mat.cohort.0$`Burkina Faso` %>%
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
LQ.baseline.0.smooth <- as(model.matrix(adjusted~as.factor(period5)-1,data=bf5.f.0.smooth),"sparseMatrix")


basepop_init <- as.numeric(bf.pop.aggr.f$'1960') %>% ifelse(.==0, 0.5, .)

fx_init <- burkina.faso.females$fertility.rates[4:10, ]
fx_init <- reshape2::melt(fx_init) %>% 
  slice(c(rep(0:(ncol(fx_init)-1) * nrow(fx_init), each=nrow(fx_init) * 5) + 1:nrow(fx_init),
          rep((ncol(fx_init)-1) * nrow(fx_init) + 1:nrow(fx_init), max(bf.idx5$periods - 2004)) )) %>%
  mutate(year = rep(1960:max(bf.idx5$periods), each = nrow(fx_init))) %>%
  select(-2) %>%
  pivot_wider(names_from = year, values_from = value) %>%
  select(num_range("",bf.idx5$periods)) %>%
  as.matrix()

log_basepop_mean <- as.vector(log(basepop_init))
log_fx_mean <- as.vector(log(fx_init))

##prior means for LQ h from IGME 5q0
igme.5q0<-read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/IGME BF 5q0.csv")
igme.5q0.df <- reshape2::melt(igme.5q0) %>% mutate(year = as.numeric(gsub("X|\\.5","",variable)), child.mort = value/1000) %>%
  select(Sex,year,child.mort)

igme.5q0.5 <- reshape2::melt(igme.5q0) %>% mutate(year = as.numeric(gsub("X","",variable)), child.mort = value/1000) %>%
  select(Sex,year,child.mort) %>%
  mutate(year5 = 5 * floor(year/5)) %>% group_by(Sex, year5) %>%
  summarise_at(vars(child.mort),mean) %>% ungroup()


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


data.f <- as.matrix(log(bf_census_f[,-1]))
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
                 #census_year_grow_idx = rep(0,4),
                 
                 open_idx = bf.idx5$n_ages, #open age starting at 70
                 LQ_baseline_mx_DX_f = LQ.baseline.DX.ax,
                 h_DX_f = LQ.baseline.DX.bx,
                 h2_DX_f = LQ.baseline.DX.cx,
                 k_DX_f = LQ.baseline.DX.vx,
                 tp_DX_f = tips.DX,

                 penal_tp = as(crossprod(diff(diag(15))),"sparseMatrix"),
                 null_penal_tp = as(exp(15)*tcrossprod(c(0,1,1,1,rep(0,11))),"sparseMatrix"),
                 penal_tp_0 = as(tcrossprod(c(1,rep(0,14))),"sparseMatrix"),
                 
                 LQ_baseline_f = as.matrix(LQcoef.f[,1:4]),
                 df = bf5.f.no0$event,
                 Ef = bf5.f.no0$pyears,
                 
                 #LQ_baseline_0 = LQ.baseline.0,
                 #tp_DX_f_0 = tips.DX.0,
                 #df0 = bf5.f.0$event,
                 #Ef0 = bf5.f.0$pyears,
                 
                 h_mean_f = igme.5q0.5 %>% filter(Sex=="Female", year5 %in% bf.idx5$periods) %>% .$child.mort %>% log(),
                 
                 pop_start = 1, pop_end = open.age / 5 + 1
)

data.vec.smooth <- list(log_basepop_mean_f = log_basepop_mean,
                        log_fx_mean = log_fx_mean,
                        srb = rep(1.05, bf.idx5$n_periods),
                        interval = bf.idx5$interval,
                        n_periods = bf.idx5$n_periods,
                        fx_idx = 4L,
                        n_fx = 7L,
                        census_log_pop_f = data.f,
                        census_year_idx = match(bf.idx5$interval * floor(as.numeric(colnames(data.f)) / bf.idx5$interval), bf.idx5$periods),
                        census_year_grow_idx = as.numeric(colnames(data.f)) - bf.idx5$interval * floor(as.numeric(colnames(data.f)) / bf.idx5$interval),
                        #census_year_grow_idx = rep(0,4),
                        
                        open_idx = bf.idx5$n_ages,
                        LQ_baseline_mx_DX_f = LQ.baseline.DX.ax.smooth,
                        h_DX_f = LQ.baseline.DX.bx.smooth,
                        h2_DX_f = LQ.baseline.DX.cx.smooth,
                        k_DX_f = LQ.baseline.DX.vx.smooth,
                        tp_DX_f = tips.DX.smooth,
                        LQ_baseline_0 = LQ.baseline.0.smooth,
                        tp_DX_f_0 = tips.DX.0.smooth,
                        
                        penal_tp = as(crossprod(diff(diag(15))),"sparseMatrix"),
                        null_penal_tp = as(exp(15)*tcrossprod(c(0,1,1,1,rep(0,11))),"sparseMatrix"),
                        penal_tp_0 = as(tcrossprod(c(1,rep(0,14))),"sparseMatrix"),
                        
                        LQ_baseline_f = as.matrix(LQcoef.f[,1:4]),
                        df = bf5.f.no0.smooth$adjusted,
                        Ef = bf5.f.no0.smooth$pyears2,
                        df0 = bf5.f.0.smooth$adjusted,
                        Ef0 = bf5.f.0.smooth$pyears2,
                        
                        h_mean_f = igme.5q0.5 %>% filter(Sex=="Female", year5 %in% bf.idx5$periods) %>% .$child.mort %>% log(),
                        
                        pop_start = 1, pop_end = open.age / 5 + 1
)

par.vec <- list(log_tau2_logpop_f = 2,
                log_tau2_fx = 1,
                log_tau2_gx_f = 1,
                log_marginal_var_h = 2,
                log_marginal_var_k = 2,
                log_marginal_prec_h = 2,
                log_marginal_prec_k = 2,
                log_lambda_tp = 1.6,
                log_lambda_tp_0_inflated_sd = 0.3,
                
                log_dispersion = 1,
                
                log_basepop_f = log_basepop_mean,
                log_fx = log_fx_mean,
                gx_f = rep(0.5, bf.idx5$n_ages * bf.idx5$n_periods), gx_m = rep(0, bf.idx5$n_ages * bf.idx5$n_periods),
                
                h_params_f = igme.5q0.5 %>% filter(Sex=="Female", year5 %in% bf.idx5$periods) %>% .$child.mort %>% log(), 
                k_params_f = rep(0,bf.idx5$n_periods),
                
                tp_params = rep(0,15),
                
                logit_rho_h = 0,
                logit_rho_k = 0,
                h_constant_f = igme.5q0.5 %>% filter(Sex=="Female", year5 %in% bf.idx5$periods) %>% .$child.mort %>% log(),
                k_constant_f = 0,
                
                logit_rho_g_x = 0,
                logit_rho_g_t = 0,
                
                logit_eigens = c(1,1),
                angles = c(0,pi/2),
                hk_err_corr = c(0)
)


input.LQ.both.vec <- list(data = data.vec, par_init = par.vec, model = "ccmpp_vr_tmb")
input.LQ.both.vec.smooth <- list(data = data.vec.smooth, par_init = par.vec, model = "ccmpp_vr_tmb")


system.time(AR.f.no0 <- fit_tmb(input.LQ.both.vec,inner_verbose=TRUE, random = c("log_basepop_f",
                                                                                 "log_fx",
                                                                                 "gx_f",
                                                                                 "h_params_f",
                                                                                 "k_params_f",
                                                                                 "h_constant_f",
                                                                                 "k_constant_f",
                                                                                 "tp_params"), 
                                DLL="ccmpp_LQuad_bothsexes_AR_f_no0",
                                map = list(h_constant_f = factor(rep(NA,bf.idx5$n_periods)),
                                           k_constant_f = factor(NA)
                                ))
) 



system.time(hmean.f.no0 <- fit_tmb(input.LQ.both.vec,inner_verbose=TRUE, random = c("log_basepop_f",
                                                                                    "log_fx",
                                                                                    "gx_f",
                                                                                    "h_params_f",
                                                                                    "k_params_f",
                                                                                    "h_constant_f",
                                                                                    "k_constant_f",
                                                                                    "tp_params"), 
                                   DLL="ccmpp_LQuad_bothsexes_weighted_hmean_f_no0",
                                   map = list(h_constant_f = factor(rep(NA,bf.idx5$n_periods)),
                                              k_constant_f = factor(NA)
                                   ))
) 

system.time(MVN.f.no0 <- fit_tmb(input.LQ.both.vec,inner_verbose=TRUE, random = c("log_basepop_f",
                                                                                 "log_fx",
                                                                                 "gx_f",
                                                                                 "h_params_f",
                                                                                 "k_params_f",
                                                                                 "h_constant_f",
                                                                                 "k_constant_f",
                                                                                 "tp_params"), 
                                DLL="ccmpp_LQuad_bothsexes_MVN_f_no0",
                                map = list(h_constant_f = factor(rep(NA,bf.idx5$n_periods)),
                                           k_constant_f = factor(NA)
                                ))
) 


system.time(ARIMAk.f.no0 <- fit_tmb(input.LQ.both.vec,inner_verbose=TRUE, random = c("log_basepop_f",
                                                                                    "log_fx",
                                                                                    "gx_f",
                                                                                    "h_params_f",
                                                                                    "k_params_f",
                                                                                    "h_constant_f",
                                                                                    "tp_params"), 
                                   DLL="ccmpp_LQuad_bothsexes_ARIMAk_f_no0",
                                   map = list(h_constant_f = factor(rep(NA,bf.idx5$n_periods))
                                   ))
) 

system.time(VAR.f.no0 <- fit_tmb(input.LQ.both.vec,inner_verbose=TRUE, random = c("log_basepop_f",
                                                                                 "log_fx",
                                                                                 "gx_f",
                                                                                 "h_params_f",
                                                                                 "k_params_f",
                                                                                 "h_constant_f",
                                                                                 "k_constant_f",
                                                                                 "tp_params"), 
                                DLL="ccmpp_LQuad_bothsexes_VAR_f_no0",
                                map = list(h_constant_f = factor(rep(NA,bf.idx5$n_periods)),
                                           k_constant_f = factor(NA)
                                ))
) 


#45q15####
q4515.func <- function(x){
  as_tibble(apply(x$mode$mx_mat_f[4:12,],2,function(x){1-prod((1-2.5*x)/(1+2.5*x))})) %>%
    mutate(sex="female", year = bf.idx5$periods)
}

wpp.bf.q4515<-read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/WPP_BF_45q15.csv")

q4515.df <- lapply(models.list, q4515.func) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>%
  bind_rows(
    as_tibble(reshape2::melt(wpp.bf.q4515[1,-c(1:5)]/1000)) %>% mutate(model="WPP Estimates", year=bf.idx5$periods, variable = NULL)
  ) %>%
  mutate(model = fct_relevel(model, "WPP Estimates"))


q4515.df %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) + ylab(bquote(""[45]*q[15])) +
  theme(text = element_text(size=25))


#h params####
get.h <- function(x) {
  x$par.full %>% split(names(.)) %>% .$h_params_f %>% as_tibble() %>% mutate(sex="female", year = bf.idx5$periods)
}

models.list <- list("AR-h AR-k" = AR.f.no0, "hmean AR-k" = hmean.f.no0, "MVN-h MVN-k" = MVN.f.no0, "AR-h ARIMA-k" = ARIMAk.f.no0)

h.df <- lapply(models.list, get.h) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>%
  bind_rows(
    filter(igme.5q0.5, year5 %in% bf.idx5$periods, Sex=="Female") %>%
      mutate(model="IGME Estimates", value=log(child.mort), year=year5, year5=NULL, child.mort=NULL, Sex=NULL)
  ) %>%
  mutate(model = fct_relevel(model, "IGME Estimates"))

ggplot(h.df) + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) + ylab(bquote("log("[5]*q[0]*") or LogQuad h")) + 
  theme(text = element_text(size=25))

#k params####
get.k <- function(x) {
  x$par.full %>% split(names(.)) %>% .$k_params_f %>% as_tibble() %>% mutate(sex="female", year = bf.idx5$periods)
}

k.df <- lapply(models.list, get.k) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows()

ggplot(k.df) + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) + ylab("LogQuad k") +
  theme(text = element_text(size=25))

#population####
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
    bf_census_f %>%
      mutate(age = c(bf_census_f$age), model="UNPD Census") %>%
      pivot_longer(!age & !model) %>% mutate(year=as.numeric(name), name=NULL)
    ) %>%
  bind_rows(
      mutate(bf.pop.aggr.f[,-1], age = bf_census_f$age, model="WPP Estimates") %>%
      pivot_longer(!age & !model) %>% mutate(year=as.numeric(name), name=NULL)
  ) %>% 
  mutate(cohort = year - age,
         cohort = 5 * floor(cohort / 5),
         model = fct_relevel(model, "UNPD Census"))


pop.df %>% filter(year %in% c(1960, 1975, 1985, 1995, 2005, 2015)) %>%
  ggplot() + geom_line(aes(x = age, y = value, col = model), lwd = 1.2) +
  theme(text = element_text(size=25)) +
  facet_wrap_paginate(~year, scale="free_y", nrow = 2, ncol = 3, page=1)

pop.df %>% filter(model != "WPP Estimates") %>%
  ggplot() + geom_line(aes(x = age, y = value, col = model), lwd = 1.2) +
  theme(text = element_text(size=25)) +
  facet_wrap_paginate(~cohort, scale="free_y", nrow = 2, ncol = 3, page=3)

pop.df %>% filter(model != "WPP Estimates") %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) +
  theme(text = element_text(size=25)) +
  facet_wrap_paginate(~age, scale="free", nrow = 2, ncol = 3, page=1)


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
  facet_wrap_paginate(~ year, nrow=2, ncol=3, page=1)

inner_join(mig.df, pop.df) %>% mutate(gx = mig/value) %>%
  ggplot() + geom_line(aes(x = age, y = gx, col = model), lwd = 1.2) +
  theme(text = element_text(size=25)) +
  facet_wrap_paginate(~ cohort, nrow=2, ncol=3, page=3)

inner_join(mig.df, pop.df) %>% mutate(gx = mig/value) %>%
  ggplot() + geom_line(aes(x = year, y = gx, col = model), lwd = 1.2) +
  theme(text = element_text(size=25)) +
  facet_wrap_paginate(~ age, scales="free_y", nrow=2, ncol=3, page=3)


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
  facet_wrap_paginate(~year, scale="free_y", nrow = 2, ncol = 3, page=2)

fx.df %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) +
  theme(text = element_text(size=25)) +
  facet_wrap_paginate(~age, scale="free_y", nrow = 2, ncol = 3, page=1)

#births####
get.births<- function(x){
  pop.mat <- matrix(x$mode$births, bf.idx5$n_fx, bf.idx5$n_periods)
  pop.mat <- rbind(pop.mat, apply(pop.mat,2,sum))
  rownames(pop.mat) <- c(bf.idx5$fertility_ages, "All")
  colnames(pop.mat) <- bf.idx5$periods
  reshape2::melt(pop.mat) %>% select(age = Var1, year = Var2, value = value)
}

births.df <- bind_rows(
    mutate(get.births(AR.f), model="ARh ARk"),
    mutate(get.births(ARIMA.f), model="ARh ARIMAk"),
    mutate(get.births(ARh.MVNk.f), model="ARh MVNk")
  )

births.df %>% mutate(age = factor(age, levels = c(bf.idx5$fertility_ages, "All"))) %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) + 
  facet_wrap_paginate( ~ age, scales = "free", nrow = 2, ncol = 3, page = 1)

#deaths####
get.period.deaths<- function(x){
  pop.mat <- matrix(x$mode$period_deaths_f, bf.idx5$n_ages, bf.idx5$n_periods)
  pop.mat <- rbind(pop.mat, apply(pop.mat,2,sum))
  rownames(pop.mat) <- c(bf.idx5$ages, "All")
  colnames(pop.mat) <- bf.idx5$periods
  reshape2::melt(pop.mat) %>% select(age = Var1, year = Var2, value = value)
}
get.cohort.deaths<- function(x){
  pop.mat <- matrix(x$mode$cohort_deaths_f, bf.idx5$n_ages+1, bf.idx5$n_periods)
  pop.mat <- rbind(pop.mat, apply(pop.mat,2,sum))
  rownames(pop.mat) <- c("0-",bf.idx5$ages, "All")
  colnames(pop.mat) <- bf.idx5$periods
  reshape2::melt(pop.mat) %>% select(age = Var1, year = Var2, value = value)
}

period.deaths.df <- bind_rows(
  mutate(get.period.deaths(AR.f), model="ARh ARk"),
  mutate(get.period.deaths(ARIMA.f), model="ARh ARIMAk"),
  mutate(get.period.deaths(ARh.MVNk.f), model="ARh MVNk")
)

cohort.deaths.df <- bind_rows(
  mutate(get.cohort.deaths(AR.f), model="ARh ARk"),
  mutate(get.cohort.deaths(ARIMA.f), model="ARh ARIMAk"),
  mutate(get.cohort.deaths(ARh.MVNk.f), model="ARh MVNk")
)


period.deaths.df %>% mutate(age = factor(age, levels = c(bf.idx5$ages, "All"))) %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) + 
  facet_wrap_paginate( ~ age, scales = "free", nrow = 2, ncol = 3, page = 1)

cohort.deaths.df %>% mutate(age = factor(age, levels = c("0-", bf.idx5$ages, "All"))) %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) + 
  facet_wrap_paginate( ~ age, scales = "free", nrow = 2, ncol = 3, page = 1)

#sx####
get.sx <- function(x){
  pop.mat <- matrix(x$mode$sx_mat_f, bf.idx5$n_ages+1, bf.idx5$n_periods)
  rownames(pop.mat) <- c("0-",bf.idx5$ages)
  colnames(pop.mat) <- bf.idx5$periods
  reshape2::melt(pop.mat) %>% select(age = Var1, year = Var2, value = value)
}


#empirical sx
emp.sx.f <- matrix(0, nrow(data.vec$census_log_pop_f)-1, ncol(data.vec$census_log_pop_f)-1)

for(j in 1:(ncol(data.vec$census_log_pop_f)-1)) {
  for (i in 1:(nrow(data.vec$census_log_pop_f)-2)) {
    emp.sx.f[i,j] <- sqrt(as.numeric(exp(data.vec$census_log_pop_f)[i+2,j+1] / exp(data.vec$census_log_pop_f)[i,j]))
  }
}

emp.sx.f.wpp <- matrix(0, nrow(bf.pop.aggr.f)-1, ncol(bf.pop.aggr.f)-1-1)

for(j in 1:(ncol(bf.pop.aggr.f)-1-1)) {
  for (i in 1:(nrow(bf.pop.aggr.f)-1)) {
    emp.sx.f.wpp[i,j] <- as.numeric(bf.pop.aggr.f[,-1][i+1,j+1] / bf.pop.aggr.f[,-1][i,j])
  }
}

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

matrix(AR.f$mode$sx_mat_f, bf.idx5$n_ages+1, bf.idx5$n_periods)
matrix(AR.f$mode$population_f, bf.idx5$n_ages, bf.idx5$n_periods+1)

AR.f$mode$infants_f*(1-matrix(AR.f$mode$sx_mat_f, bf.idx5$n_ages+1, bf.idx5$n_periods)[1,])
matrix(AR.f$mode$cohort_deaths_f, bf.idx5$n_ages+1, bf.idx5$n_periods)[1,]

#mx####
get.mx<- function(x){
  pop.mat <- matrix(x$mode$mx_mat_f, nrow(LQcoef.f)+1, bf.idx5$n_periods)
  rownames(pop.mat) <- seq(0, 110, by = 5)
  colnames(pop.mat) <- bf.idx5$periods
  reshape2::melt(pop.mat) %>% select(age = Var1, year = Var2, value = value)
}

mx.df <- lapply(models.list, get.mx) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows()


bf5.f.plot <- bf5.f.no0 %>% select(year = period5, age = age5, tips = tips, event = event, pyears = pyears) %>%
  mutate(year = as.numeric(as.character(year)), raw.mort = event/pyears)


mx.df %>%
  ggplot() + geom_line(aes(x = age, y = value, col = model), lwd = 1.2) +
  geom_point(data = filter(bf5.f.plot, raw.mort!=0), aes(x = age, y = raw.mort, shape = tips), alpha = 0) + 
  geom_point(data = filter(bf5.f.plot, raw.mort!=0, tips %in% 1:3), aes(x = age, y = raw.mort, shape = tips), size=3, stroke=2) +
  geom_point(data = filter(bf5.f.plot, raw.mort!=0, !tips %in% 1:3), aes(x = age, y = raw.mort, shape = tips), size=2, alpha = 0.15) +
  scale_shape_manual(values = 0:14) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  theme(text = element_text(size=25)) +
  facet_wrap_paginate(~year, nrow = 2, ncol = 3, page=2)



mx.df %>% filter(age %in% 15:55) %>%
  ggplot() + geom_line(aes(x = year, y = value, col = model), lwd = 1.2) +
  geom_point(data = filter(bf5.f.plot, raw.mort!=0), aes(x = year, y = raw.mort, shape = tips), alpha = 0) + 
  geom_point(data = filter(bf5.f.plot, raw.mort!=0, tips %in% 1:3), aes(x = year, y = raw.mort, shape = tips), size=3, stroke=2) +
  geom_point(data = filter(bf5.f.plot, raw.mort!=0, !tips %in% 1:3), aes(x = year, y = raw.mort, shape = tips), size=2, alpha = 0.15) +
  scale_shape_manual(values = 0:14) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  theme(text = element_text(size=25)) +
  facet_wrap_paginate(~age, scales = "free_y", nrow = 2, ncol = 3, page=1)

