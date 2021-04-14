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
library(MortCast)
library(shinystan)

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_innov.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_innov"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_LQ_DEF.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_LQ_DEF"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_LQ_DEF_innov.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_LQ_DEF_innov"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_LQ_DEF_innov_TN.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_LQ_DEF_innov_TN"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_LQ_DEF_innov_TN_onek.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_LQ_DEF_innov_TN_onek"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_LQ_DEF_innov_onek_varyDEF.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_LQ_DEF_innov_onek_varyDEF"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_LQ_DEF_innov_onek_varyDEF_noprior.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_LQ_DEF_innov_onek_varyDEF_noprior"))

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_LQ_DEF_innov_varyk_oneDEF.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_LQ_DEF_innov_varyk_oneDEF"))

load("~/cohort smooth 1900-2017.RData")
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


##Log Quad
LQcoef.f <- LQcoef %>% filter(sex=="Female", !age%in%c("0","1-4")) %>% select(ax:vx)
LQcoef.f$age5 <- seq(5,110,by=5)

LQcoef.m <- LQcoef %>% filter(sex=="Male", !age%in%c("0","1-4")) %>% select(ax:vx)
LQcoef.m$age5 <- seq(5,110,by=5)

igme.5q0.f<-read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/5q0 male.csv")
igme.5q0.m<-read.csv("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/5q0 female.csv")

country <- "Zimbabwe"

##prior means for LQ h from IGME 5q0
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

##DHS data cohort smoothed
bf5.smooth <- aggr.mat.cohort.0[[country]] %>%
  mutate(age5 = 5 * floor(agegr / 5),
         period5 = 5 * floor(period / 5)) %>%
  group_by(mm1, tips, DHS, age5, period5) %>%
  summarise_at(vars(pyears, event, pyears2, adjusted), sum) %>%
  ungroup()

bf5.smooth$period5 <- factor(bf5.smooth$period5)
bf5.smooth$tips <- factor(bf5.smooth$tips,levels=0:14)
bf5.f.smooth <- bf5.smooth %>% filter(mm1=="female") %>% arrange(period5,tips,age5)
bf5.f.no0.smooth <- filter(bf5.f.smooth, age5  >= 15)

bf5.m.smooth <- bf5.smooth %>% filter(mm1=="male") %>% arrange(period5,tips,age5)
bf5.m.no0.smooth <- filter(bf5.m.smooth, age5 >= 15)

tips.DX.smooth <- as(model.matrix(event~as.factor(tips)-1,data=bf5.f.no0.smooth),"sparseMatrix")
tips.DX.smooth.m <- as(model.matrix(event~as.factor(tips)-1,data=bf5.m.no0.smooth),"sparseMatrix")
#tips.DX.0.smooth <- as(model.matrix(event~as.factor(tips)-1,data=bf5.f.0.smooth),"sparseMatrix")

LQ.baseline.DX.smooth <- model.matrix(adjusted~factor(age5, levels = seq(5, 110, by = 5)) - 1, data=bf5.f.no0.smooth) %*% as.matrix(LQcoef.f[,1:4])
LQ.baseline.DX.ax.smooth <- LQ.baseline.DX.smooth[,1]
LQ.baseline.DX.bx.smooth <- data.frame(period5 = levels(bf5.smooth$period5)) %>% 
  bind_rows(bind_cols(n = 1 : nrow(bf5.f.no0.smooth), ax=LQ.baseline.DX.smooth[,2], period5=bf5.f.no0.smooth$period5)) %>% 
  pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
  slice(-1) %>% select(-n) %>%
  as.matrix() %>% as("sparseMatrix")
LQ.baseline.DX.cx.smooth <- data.frame(period5 = levels(bf5.smooth$period5)) %>% 
  bind_rows(bind_cols(n = 1 : nrow(bf5.f.no0.smooth), ax=LQ.baseline.DX.smooth[,3], period5=bf5.f.no0.smooth$period5)) %>% 
  pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
  slice(-1) %>% select(-n) %>%
  as.matrix() %>% as("sparseMatrix")
LQ.baseline.DX.vx.smooth <- data.frame(period5 = levels(bf5.smooth$period5)) %>% 
  bind_rows(bind_cols(n = 1 : nrow(bf5.f.no0.smooth), ax=LQ.baseline.DX.smooth[,4], period5=bf5.f.no0.smooth$period5)) %>% 
  pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
  slice(-1) %>% select(-n) %>%
  as.matrix() %>% as("sparseMatrix")

LQ.baseline.DX.smooth.m <- model.matrix(adjusted~factor(age5, levels = seq(5, 110, by = 5)) - 1, data=bf5.m.no0.smooth) %*% as.matrix(LQcoef.m[,1:4])
LQ.baseline.DX.ax.smooth.m <- LQ.baseline.DX.smooth.m[,1]
LQ.baseline.DX.bx.smooth.m <- data.frame(period5 = levels(bf5.smooth$period5)) %>% 
  bind_rows(bind_cols(n = 1 : nrow(bf5.m.no0.smooth), ax=LQ.baseline.DX.smooth.m[,2], period5=bf5.m.no0.smooth$period5)) %>% 
  pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
  slice(-1) %>% select(-n) %>%
  as.matrix() %>% as("sparseMatrix")
LQ.baseline.DX.cx.smooth.m <- data.frame(period5 = levels(bf5.smooth$period5)) %>% 
  bind_rows(bind_cols(n = 1 : nrow(bf5.m.no0.smooth), ax=LQ.baseline.DX.smooth.m[,3], period5=bf5.m.no0.smooth$period5)) %>% 
  pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
  slice(-1) %>% select(-n) %>%
  as.matrix() %>% as("sparseMatrix")
LQ.baseline.DX.vx.smooth.m <- data.frame(period5 = levels(bf5.smooth$period5)) %>% 
  bind_rows(bind_cols(n = 1 : nrow(bf5.m.no0.smooth), ax=LQ.baseline.DX.smooth.m[,4], period5=bf5.m.no0.smooth$period5)) %>% 
  pivot_wider(names_from=period5, values_from=ax, values_fill=0) %>% 
  slice(-1) %>% select(-n) %>%
  as.matrix() %>% as("sparseMatrix")

data.vec.f <- list(LQ_baseline_mx_DX_f = LQ.baseline.DX.ax.smooth,
                   h_DX_f = LQ.baseline.DX.bx.smooth,
                   h2_DX_f = LQ.baseline.DX.cx.smooth,
                   k_DX_f = LQ.baseline.DX.vx.smooth,
                   tp_DX_f = tips.DX.smooth,
                   
                   penal_tp = as(crossprod(diff(diag(15))),"sparseMatrix"),
                   null_penal_tp = as(exp(15)*tcrossprod(c(0,1,1,1,rep(0,11))),"sparseMatrix"),
                   penal_tp_0 = as(tcrossprod(c(1,rep(0,14))),"sparseMatrix"),
                   
                   LQ_baseline_f = as.matrix(LQcoef.f[,1:4]),
                   df = bf5.f.no0.smooth$adjusted,
                   Ef = bf5.f.no0.smooth$pyears2,
                   
                   h_mean_f = igme.5q0.5 %>% filter(Sex=="Female", year5 %in% levels(bf5.f.no0.smooth$period5)) %>% .$child.mort %>% log(),
                   h_constant_f = igme.5q0.5 %>% filter(Sex=="Female", year5 %in% levels(bf5.f.no0.smooth$period5)) %>% .$child.mort %>% log(),
                   
                   DEF_log_age = log(bf5.f.no0.smooth$age5 + 2),
                   DEF_time = as(model.matrix(adjusted~period5-1,data = bf5.f.no0.smooth), "sparseMatrix"),
                   
                   sigma_D_mean = 100,
                   sigma_E_mean = 100,
                   sigma_F_mean = 100
)

par.vec.f <- list(log_marginal_prec_h = 3,
                  log_marginal_prec_k = 3,
                  log_lambda_tp = 0,
                  log_lambda_tp_0_inflated_sd = 0.3,
                  
                  log_dispersion = 1,
                  
                  h_params_f = igme.5q0.5 %>% filter(Sex=="Female", year5 %in% levels(bf5.m.no0.smooth$period5)) %>% .$child.mort %>% log(), 
                  k_params_f = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                  
                  tp_params = rep(0,15),
                  
                  logit_rho_h = 0,
                  logit_rho_k = 0,
                  #h_constant_f = igme.5q0.5 %>% filter(Sex=="Female", year5 %in% levels(bf5.m.no0.smooth$period5)) %>% .$child.mort %>% log(),
                  #k_constant_f = 0,
                  
                  log_D = rep(log(0.005), length(levels(bf5.f.no0.smooth$period5))),
                  log_E = rep(log(5), length(levels(bf5.f.no0.smooth$period5))),
                  log_F = rep(log(35), length(levels(bf5.f.no0.smooth$period5))),
                  
                  log_D_mean = -6,
                  log_E_mean = 1.5,
                  log_F_mean = 3.5,
                  
                  log_marginal_prec_D = 3,
                  log_marginal_prec_E = 3,
                  log_marginal_prec_F = 3,
                  
                  logit_rho_D = 0,
                  logit_rho_E = 0,
                  logit_rho_F = 0,
                  
                  h_params_f_innov = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                  k_params_f_innov = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                  
                  log_D_innov = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                  log_E_innov = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                  log_F_innov = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                  k_param_f = 0
)


data.vec.m <- list(LQ_baseline_mx_DX_f = LQ.baseline.DX.ax.smooth.m,
                   h_DX_f = LQ.baseline.DX.bx.smooth.m,
                   h2_DX_f = LQ.baseline.DX.cx.smooth.m,
                   k_DX_f = LQ.baseline.DX.vx.smooth.m,
                   tp_DX_f = tips.DX.smooth.m,
                   
                   penal_tp = as(crossprod(diff(diag(15))),"sparseMatrix"),
                   null_penal_tp = as(exp(15)*tcrossprod(c(0,1,1,1,rep(0,11))),"sparseMatrix"),
                   penal_tp_0 = as(tcrossprod(c(1,rep(0,14))),"sparseMatrix"),
                   
                   LQ_baseline_f = as.matrix(LQcoef.m[,1:4]),
                   df = bf5.m.no0.smooth$adjusted,
                   Ef = bf5.m.no0.smooth$pyears2,
                   
                   h_mean_f = igme.5q0.5 %>% filter(Sex=="Male", year5 %in% levels(bf5.m.no0.smooth$period5)) %>% .$child.mort %>% log(),
                   
                   DEF_log_age = log(bf5.m.no0.smooth$age5 + 2),
                   DEF_time = as(model.matrix(adjusted ~ period5-1,data = bf5.m.no0.smooth), "sparseMatrix")
)


par.vec.m <- list(log_marginal_prec_h = 1,
                  log_marginal_prec_k = 1,
                  log_lambda_tp = 1,
                  log_lambda_tp_0_inflated_sd = 0.3,
                  
                  log_dispersion = 1,
                  
                  h_params_f = igme.5q0.5 %>% filter(Sex=="Male", year5 %in% levels(bf5.m.no0.smooth$period5)) %>% .$child.mort %>% log(), 
                  k_params_f = rep(0,length(levels(bf5.m.no0.smooth$period5))),
                  
                  tp_params = rep(0,15),
                  
                  logit_rho_h = 0,
                  logit_rho_k = 0,
                  h_constant_f = igme.5q0.5 %>% filter(Sex=="Male", year5 %in% levels(bf5.m.no0.smooth$period5)) %>% .$child.mort %>% log(),
                  k_constant_f = 0,
                  
                  log_D = rep(log(0.005), length(levels(bf5.m.no0.smooth$period5))),
                  log_E = rep(log(5), length(levels(bf5.m.no0.smooth$period5))),
                  log_F = rep(log(35), length(levels(bf5.m.no0.smooth$period5))),
                  
                  log_D_mean = log(0.005),
                  log_E_mean = log(5),
                  log_F_mean = log(35),
                  
                  log_marginal_prec_D = 1,
                  log_marginal_prec_E = 1,
                  log_marginal_prec_F = 1,
                  
                  logit_rho_D = 1,
                  logit_rho_E = 1,
                  logit_rho_F = 1,
                  
                  h_params_f_innov = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                  k_params_f_innov = rep(0, length(levels(bf5.f.no0.smooth$period5)))
)

input.LQ.both.vec.f <- list(data = data.vec.f, par_init = par.vec.f, model = "ccmpp_vr_tmb")
input.LQ.both.vec.m <- list(data = data.vec.m, par_init = par.vec.m, model = "ccmpp_vr_tmb")

just.DHS.m <- fit_tmb(input.LQ.both.vec.m,inner_verbose=TRUE, random = c("h_params_f",
                                                                         "k_params_f",
                                                                         "h_constant_f",
                                                                         "k_constant_f",
                                                                         "tp_params"), 
                      DLL="just_DHS_f",
                      map = list(h_constant_f = factor(rep(NA,length(levels(bf5.m.no0.smooth$period5)))),
                                 k_constant_f = factor(NA)
                      ))



just.DHS.fixh.m <- fit_tmb(input.LQ.both.vec.m,inner_verbose=TRUE, random = c("h_params_f",
                                                                              "k_params_f",
                                                                              "h_constant_f",
                                                                              "k_constant_f",
                                                                              "tp_params"), 
                           DLL="just_DHS_f",
                           map = list(h_constant_f = factor(rep(NA,length(levels(bf5.m.no0.smooth$period5)))),
                                      k_constant_f = factor(NA),
                                      h_params_f = factor(rep(NA,length(levels(bf5.m.no0.smooth$period5)))),
                                      logit_rho_h = factor(NA),
                                      log_marginal_prec_h = factor(NA)
                           ))

just.DHS.m.innov <- fit_tmb(input.LQ.both.vec.m,inner_verbose=TRUE, random = c("h_params_f_innov",
                                                                         "k_params_f_innov",
                                                                         "h_constant_f",
                                                                         "tp_params"), 
                      DLL="just_DHS_f_innov",
                      map = list(h_constant_f = factor(rep(NA,length(levels(bf5.m.no0.smooth$period5))))
                      ))

just.DHS.m.DEF <- fit_tmb(input.LQ.both.vec.m,inner_verbose=TRUE, random = c("h_params_f",
                                                                         "k_params_f",
                                                                         "h_constant_f",
                                                                         "k_constant_f",
                                                                         "tp_params",
                                                                         "log_D",
                                                                         "log_E",
                                                                         "log_F",
                                                                         "log_D_mean",
                                                                         "log_E_mean",
                                                                         "log_F_mean"), 
                      DLL="just_DHS_f_LQ_DEF",
                      map = list(h_constant_f = factor(rep(NA,length(levels(bf5.m.no0.smooth$period5)))),
                                 k_constant_f = factor(NA),
                                 log_D = factor(rep(1, length(levels(bf5.m.no0.smooth$period5)))),
                                 log_E = factor(rep(1, length(levels(bf5.m.no0.smooth$period5)))),
                                 log_F = factor(rep(1, length(levels(bf5.m.no0.smooth$period5)))),
                                 logit_rho_D = factor(NA),
                                 logit_rho_E = factor(NA),
                                 logit_rho_F = factor(NA),
                                 log_D_mean = factor(NA),
                                 log_E_mean = factor(NA),
                                 log_F_mean = factor(NA),
                                 log_marginal_prec_D = factor(NA),
                                 log_marginal_prec_E = factor(NA),
                                 log_marginal_prec_F = factor(NA)
                                 ))


just.DHS.m.DEF <- fit_tmb(input.LQ.both.vec.m,inner_verbose=TRUE, random = c("h_params_f",
                                                                             "k_params_f",
                                                                             "h_constant_f",
                                                                             "k_constant_f",
                                                                             "tp_params",
                                                                             "log_D",
                                                                             "log_E",
                                                                             "log_F",
                                                                             "log_D_mean",
                                                                             "log_E_mean",
                                                                             "log_F_mean"), 
                          DLL="just_DHS_f_LQ_DEF",
                          map = list(h_constant_f = factor(rep(NA,length(levels(bf5.m.no0.smooth$period5)))),
                                     k_constant_f = factor(NA)
                                     
                                     #log_D = factor(rep(1, length(levels(bf5.m.no0.smooth$period5)))),
                                     #logit_rho_D = factor(NA),
                                     #log_D_mean = factor(NA),
                                     #log_marginal_prec_D = factor(1),
                                     
                                     #log_E = factor(rep(1, length(levels(bf5.m.no0.smooth$period5)))),
                                     #logit_rho_E = factor(NA),
                                     #log_E_mean = factor(NA),
                                     #log_marginal_prec_E = factor(1),
                                     
                                     #log_F = factor(rep(1, length(levels(bf5.m.no0.smooth$period5)))),
                                     #logit_rho_F = factor(NA),
                                     #log_F_mean = factor(NA),
                                     #log_marginal_prec_F = factor(1)
                          ))



just.DHS.f <- fit_tmb(input.LQ.both.vec.f,inner_verbose=TRUE, random = c("h_params_f",
                                                                         "k_params_f",
                                                                         "h_constant_f",
                                                                         "k_constant_f",
                                                                         "tp_params"), 
                      DLL="just_DHS_f",
                      map = list(h_constant_f = factor(rep(NA,length(levels(bf5.m.no0.smooth$period5)))),
                                 k_constant_f = factor(NA)
                      ))



just.DHS.fixh.f <- fit_tmb(input.LQ.both.vec.f,inner_verbose=TRUE, random = c("h_params_f",
                                                                              "k_params_f",
                                                                              "h_constant_f",
                                                                              "k_constant_f",
                                                                              "tp_params"), 
                           DLL="just_DHS_f",
                           map = list(h_constant_f = factor(rep(NA,length(levels(bf5.m.no0.smooth$period5)))),
                                      k_constant_f = factor(NA),
                                      h_params_f = factor(rep(NA,length(levels(bf5.m.no0.smooth$period5)))),
                                      logit_rho_h = factor(NA),
                                      log_marginal_prec_h = factor(NA)
                           ))


just.DHS.f.innov <- fit_tmb(input.LQ.both.vec.f,inner_verbose=TRUE, random = c("h_params_f_innov",
                                                                               "k_params_f_innov",
                                                                               "h_constant_f",
                                                                               "tp_params"), 
                            DLL="just_DHS_f_innov",
                            map = list(h_constant_f = factor(rep(NA,length(levels(bf5.m.no0.smooth$period5))))
                            ))



just.DHS.f.DEF <- fit_tmb(input.LQ.both.vec.f,inner_verbose=TRUE, random = c("h_params_f",
                                                                         "k_params_f",
                                                                         "h_constant_f",
                                                                         "k_constant_f",
                                                                         "tp_params",
                                                                         "log_D",
                                                                         "log_E",
                                                                         "log_F",
                                                                         "log_D_mean",
                                                                         "log_E_mean",
                                                                         "log_F_mean"), 
                      DLL="just_DHS_f_LQ_DEF",
                      map = list(h_constant_f = factor(rep(NA,length(levels(bf5.f.no0.smooth$period5)))),
                                 k_constant_f = factor(NA),
                                 log_D = factor(rep(1, length(levels(bf5.m.no0.smooth$period5)))),
                                 log_E = factor(rep(1, length(levels(bf5.m.no0.smooth$period5)))),
                                 log_F = factor(rep(1, length(levels(bf5.m.no0.smooth$period5)))),
                                 logit_rho_D = factor(NA),
                                 logit_rho_E = factor(NA),
                                 logit_rho_F = factor(NA),
                                 log_D_mean = factor(NA),
                                 log_E_mean = factor(NA),
                                 log_F_mean = factor(NA),
                                 log_marginal_prec_D = factor(NA),
                                 log_marginal_prec_E = factor(NA),
                                 log_marginal_prec_F = factor(NA)
                      ))



just.DHS.f.DEF <- fit_tmb(input.LQ.both.vec.f,inner_verbose=TRUE, random = c("h_params_f",
                                                                             "k_params_f",
                                                                             "h_constant_f",
                                                                             "k_constant_f",
                                                                             "tp_params",
                                                                             "log_D",
                                                                             "log_E",
                                                                             "log_F",
                                                                             "log_D_mean",
                                                                             "log_E_mean",
                                                                             "log_F_mean"), 
                          DLL="just_DHS_f_LQ_DEF",
                          map = list(h_constant_f = factor(rep(NA,length(levels(bf5.f.no0.smooth$period5)))),
                                     k_constant_f = factor(NA),
                                     
                                     log_D = factor(rep(1, length(levels(bf5.m.no0.smooth$period5)))),
                                     logit_rho_D = factor(NA),
                                     log_D_mean = factor(NA),
                                     log_marginal_prec_D = factor(1),
                                     
                                     log_E = factor(rep(1, length(levels(bf5.m.no0.smooth$period5)))),
                                     logit_rho_E = factor(NA),
                                     log_E_mean = factor(NA),
                                     log_marginal_prec_E = factor(1),
                                     
                                     log_F = factor(rep(1, length(levels(bf5.m.no0.smooth$period5)))),
                                     logit_rho_F = factor(NA),
                                     log_F_mean = factor(NA),
                                     log_marginal_prec_F = factor(1)
                          ))



just.DHS.f.1k.DEF <- fit_tmb(input.LQ.both.vec.f,inner_verbose=TRUE, random = c("h_params_f_innov",
                                                                             "k_param_f",
                                                                             "tp_params",
                                                                             "log_D_innov",
                                                                             "log_E_innov",
                                                                             "log_F_innov",
                                                                             "log_D_mean",
                                                                             "log_E_mean",
                                                                             "log_F_mean"
                                                                             ), 
                          DLL="just_DHS_f_LQ_DEF_innov_onek_varyDEF",
                          map = list(log_D_mean = factor(NA),
                                     log_E_mean = factor(NA),
                                     log_F_mean = factor(NA)
                          )
                          )

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

##males
h <- just.DHS.m$par.full %>% split(names(.)) %>% .$h_params_f
k <- just.DHS.m$par.full %>% split(names(.)) %>% .$k_params_f
fixh.k <- just.DHS.fixh.m$par.full %>% split(names(.)) %>% .$k_params_f
h.DEF <-just.DHS.m.DEF$par.full %>% split(names(.)) %>% .$h_params_f
k.DEF <-just.DHS.m.DEF$par.full %>% split(names(.)) %>% .$k_params_f
  
LQ <- as.matrix(LQcoef.m)[,1:4] %*% rbind(rep(1, length(unique(bf5.m.no0.smooth$period5))),
                                          h,
                                          h^2,
                                          k) %>% 
  `colnames<-`(unique(bf5.m.no0.smooth$period5)) %>%
  `rownames<-`(seq(5, 110, by = 5)) %>%
  reshape2::melt() %>%
  mutate(value = exp(value),
         model="LQ") %>%
  setNames(c("age5", "period5", "value", "model")) %>%
  bind_rows(
    as.matrix(LQcoef.m)[,1:4] %*% rbind(rep(1, length(unique(bf5.m.no0.smooth$period5))),
                                        data.vec.m$h_mean_f,
                                        data.vec.m$h_mean_f^2,
                                        fixh.k) %>% 
      `colnames<-`(unique(bf5.m.no0.smooth$period5)) %>%
      `rownames<-`(seq(5, 110, by = 5)) %>%
      reshape2::melt() %>%
      mutate(value = exp(value),
             model="LQ fixh") %>%
      setNames(c("age5", "period5", "value", "model"))
  )  %>%
  bind_rows(
    (exp(as.matrix(LQcoef.m)[,1:4] %*% rbind(rep(1, length(unique(bf5.m.no0.smooth$period5))),
                                             h.DEF,
                                             h.DEF^2,
                                             k.DEF)) +
       exp(just.DHS.m.DEF$par.full%>%split(names(.))%>%.$log_D -exp(just.DHS.m.DEF$par.full%>%split(names(.))%>%.$log_E) * (log(seq(5, 110, by= 5) + 2) - just.DHS.m.DEF$par.full%>%split(names(.))%>%.$log_F) ^ 2)) %>% 
      log() %>% 
      `colnames<-`(unique(bf5.m.no0.smooth$period5)) %>%
      `rownames<-`(seq(5, 110, by = 5)) %>%
      reshape2::melt() %>%
      mutate(value = exp(value),
             model="LQ DEF") %>%
      setNames(c("age5", "period5", "value", "model"))
  ) %>%
  mutate(period5 = sprintf("%d-%d",period5, period5+4)) %>%
  bind_rows(
    reshape2::melt(more.countries.avg[[country]]$mort.m) %>% 
      mutate(age5 = 5 * floor(Var1 / 5),
             period5 = 5 * floor(Var2 / 5)) %>% 
      group_by(age5, period5) %>%
      summarise_at(vars(value), mean) %>%
      mutate(value = exp(value),
             period5 = sprintf("%d-%d",period5, period5+4),
             model = "Spline average")
  ) %>%
  mutate(model = fct_relevel(model, "LQ", "LQ fixh", "LQ DEF", "Spline average"))


DHS.plot <- bf5.m.no0.smooth %>% mutate(value = adjusted / pyears2, period5 = as.numeric(levels(period5)[period5])) %>%
  select(tips, age5, period5, value) %>%
  mutate(period5 = sprintf("%d-%d", period5, period5 + 4))

ggplot(LQ %>% filter(age5 %in% 0:70)) + geom_line(aes(x = age5, y = value, linetype = model), lwd=1.2, col="blue") +
  geom_line(data = LQ %>% filter(model == "LQ DEF", age5 %in% 0:70), aes(x = age5, y = value, linetype = model), col = "orange", lwd = 1.2) +
  geom_line(data = LQ %>% filter(model == "Spline average"), aes(x = age5, y = value, linetype = model), col = "green3", lwd = 1.2) +
  geom_point(data = filter(DHS.plot, value!=0), aes(x = age5, y = value, shape = tips), alpha = 0) + 
  geom_point(data = filter(DHS.plot, value!=0, tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=3, stroke=2) +
  geom_point(data = filter(DHS.plot, value!=0, !tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=2, alpha = 0.3) +
  scale_shape_manual(values = 0:14) +
  scale_linetype_manual(values = c("solid", "dashed", "solid", "twodash"), guide = guide_legend(override.aes = list(col = c("blue", "blue", "orange", "green3")))) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  facet_wrap(~period5) + 
  ggtitle(country) + ylab(bquote(""[5]*m[x])) + xlab("5-year age groups") +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35))

hk.df <- as_tibble(h) %>% mutate(period5 = unique(bf5.m.no0.smooth$period5), par ="h", model = "LQ") %>%
  bind_rows(as_tibble(igme.5q0.5 %>% filter(Sex=="Male", year5 %in% levels(bf5.m.no0.smooth$period5)) %>% .$child.mort %>% log()) %>%
              mutate(period5 = unique(bf5.m.no0.smooth$period5), par ="h", model = "IGME"),
            as_tibble(h.DEF) %>%
              mutate(period5 = unique(bf5.m.no0.smooth$period5), par ="h", model = "LQ DEF")) %>%
  bind_rows(
    as_tibble(k) %>% mutate(period5 = unique(bf5.m.no0.smooth$period5), par ="k", model = "LQ") %>%
      bind_rows(as_tibble(fixh.k) %>%
                  mutate(period5 = unique(bf5.m.no0.smooth$period5), par ="k", model = "LQ fixh"),
                as_tibble(k.DEF) %>%
                  mutate(period5 = unique(bf5.m.no0.smooth$period5), par ="k", model = "LQ DEF"))
    ) %>%
  mutate(period5 = as.numeric(levels(period5)[period5]),
         model = fct_relevel(model, "IGME","LQ fixh", "LQ", "LQ DEF"))
  
ggplot(hk.df) +geom_line(aes(x = period5, y = value, col = model), lwd = 1.2) +
  scale_color_manual(values = c("red", "red", "blue", "green3"))+
  facet_wrap(~par, scales = "free_y") +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35))


just.DHS.m.DEF$par.full%>%split(names(.))%>%.$log_D - exp(just.DHS.m.DEF$par.full%>%split(names(.))%>%.$log_E) * (log(seq(5, 110, by= 5) + 2) - just.DHS.m.DEF$par.full%>%split(names(.))%>%.$log_F) ^ 2 %>%
  exp() %>%
  mutate()
  


 
##females
h <- just.DHS.f$par.full %>% split(names(.)) %>% .$h_params_f
k <- just.DHS.f$par.full %>% split(names(.)) %>% .$k_params_f
fixh.k <- just.DHS.fixh.f$par.full %>% split(names(.)) %>% .$k_params_f
h.DEF <-just.DHS.f.DEF$par.full %>% split(names(.)) %>% .$h_params_f
k.DEF <-just.DHS.f.DEF$par.full %>% split(names(.)) %>% .$k_params_f

sapply(just.DHS.f.1k.DEF$mode$log_F, function(x){(log(seq(5, 110, by = 5) + 2) - x) ^ 2})
  
LQ <- as.matrix(LQcoef.f)[,1:4] %*% rbind(rep(1, length(unique(bf5.f.no0.smooth$period5))),
                                          h,
                                          h^2,
                                          k) %>% 
  `colnames<-`(unique(bf5.f.no0.smooth$period5)) %>%
  `rownames<-`(seq(5, 110, by = 5)) %>%
  reshape2::melt() %>%
  mutate(value = exp(value),
         model="LQ") %>%
  setNames(c("age5", "period5", "value", "model")) %>%
  bind_rows(
    as.matrix(LQcoef.f)[,1:4] %*% rbind(rep(1, length(unique(bf5.f.no0.smooth$period5))),
                                        data.vec.f$h_mean_f,
                                        data.vec.f$h_mean_f^2,
                                        fixh.k) %>% 
      `colnames<-`(unique(bf5.f.no0.smooth$period5)) %>%
      `rownames<-`(seq(5, 110, by = 5)) %>%
      reshape2::melt() %>%
      mutate(value = exp(value),
             model="LQ fixh") %>%
      setNames(c("age5", "period5", "value", "model"))
  )  %>%
  bind_rows(
    (exp(as.matrix(LQcoef.f)[,1:4] %*% rbind(rep(1, length(unique(bf5.f.no0.smooth$period5))),
                                             h.DEF,
                                             h.DEF^2,
                                             k.DEF)) +
       exp(just.DHS.f.DEF$par.full%>%split(names(.))%>%.$log_D -exp(just.DHS.f.DEF$par.full%>%split(names(.))%>%.$log_E) * (log(seq(5, 110, by= 5) + 2) - just.DHS.f.DEF$par.full%>%split(names(.))%>%.$log_F) ^ 2)) %>% 
      log() %>% 
      `colnames<-`(unique(bf5.f.no0.smooth$period5)) %>%
      `rownames<-`(seq(5, 110, by = 5)) %>%
      reshape2::melt() %>%
      mutate(value = exp(value),
             model="LQ DEF") %>%
      setNames(c("age5", "period5", "value", "model"))
  ) %>%
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
  mutate(model = fct_relevel(model, "LQ", "LQ fixh", "LQ DEF", "Spline average"))


DHS.plot <- bf5.f.no0.smooth %>% mutate(value = adjusted / pyears2, period5 = as.numeric(levels(period5)[period5])) %>%
  select(tips, age5, period5, value) %>%
  mutate(period5 = sprintf("%d-%d", period5, period5 + 4))

ggplot(LQ %>% filter(age5 %in% 0:70)) + geom_line(aes(x = age5, y = value, linetype = model), lwd=1.2, col="red") +
  geom_line(data = LQ %>% filter(model == "LQ DEF", age5 %in% 0:70), aes(x = age5, y = value, linetype = model), col = "orange", lwd = 1.2) +
  geom_line(data = LQ %>% filter(model == "Spline average"), aes(x = age5, y = value, linetype = model), col = "green3", lwd = 1.2) +
  geom_point(data = filter(DHS.plot, value!=0), aes(x = age5, y = value, shape = tips), alpha = 0) + 
  geom_point(data = filter(DHS.plot, value!=0, tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=3, stroke=2) +
  geom_point(data = filter(DHS.plot, value!=0, !tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=2, alpha = 0.3) +
  scale_shape_manual(values = 0:14) +
  scale_linetype_manual(values = c("solid", "dashed", "solid", "twodash"), guide = guide_legend(override.aes = list(col = c("red", "red", "orange", "green3")))) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  facet_wrap(~period5) + 
  ggtitle(country) + ylab(bquote(""[5]*m[x])) + xlab("5-year age groups") +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35))

hk.df <- as_tibble(h) %>% mutate(period5 = unique(bf5.f.no0.smooth$period5), par ="h", model = "LQ") %>%
  bind_rows(as_tibble(igme.5q0.5 %>% filter(Sex=="Female", year5 %in% levels(bf5.f.no0.smooth$period5)) %>% .$child.mort %>% log()) %>%
              mutate(period5 = unique(bf5.f.no0.smooth$period5), par ="h", model = "IGME"),
            as_tibble(h.DEF) %>%
              mutate(period5 = unique(bf5.f.no0.smooth$period5), par ="h", model = "LQ DEF")) %>%
  bind_rows(
    as_tibble(k) %>% mutate(period5 = unique(bf5.f.no0.smooth$period5), par ="k", model = "LQ") %>%
      bind_rows(as_tibble(fixh.k) %>%
                  mutate(period5 = unique(bf5.f.no0.smooth$period5), par ="k", model = "LQ fixh"),
                as_tibble(k.DEF) %>%
                  mutate(period5 = unique(bf5.f.no0.smooth$period5), par ="k", model = "LQ DEF"))
  ) %>%
  mutate(period5 = as.numeric(levels(period5)[period5]),
         model = fct_relevel(model, "IGME","LQ fixh", "LQ", "LQ DEF"))

ggplot(hk.df) +geom_line(aes(x = period5, y = value, col = model), lwd = 1.2) +
  scale_color_manual(values = c("red", "red", "blue", "green3"))+
  facet_wrap(~par, scales = "free_y") +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35))


##No DHS####

#males
h <- just.DHS.m$par.full %>% split(names(.)) %>% .$h_params_f
k <- just.DHS.m$par.full %>% split(names(.)) %>% .$k_params_f
fixh.k <- just.DHS.fixh.m$par.full %>% split(names(.)) %>% .$k_params_f

LQ <- as.matrix(LQcoef.m)[,1:4] %*% rbind(rep(1, length(unique(bf5.m.no0.smooth$period5))),
                                          h,
                                          h^2,
                                          k) %>% 
  `colnames<-`(unique(bf5.m.no0.smooth$period5)) %>%
  `rownames<-`(seq(5, 110, by = 5)) %>%
  reshape2::melt() %>%
  mutate(value = exp(value),
         model="LQ") %>%
  setNames(c("age5", "period5", "value", "model")) %>%
  bind_rows(
    as.matrix(LQcoef.m)[,1:4] %*% rbind(rep(1, length(unique(bf5.m.no0.smooth$period5))),
                                        data.vec.m$h_mean_f,
                                        data.vec.m$h_mean_f^2,
                                        fixh.k) %>% 
      `colnames<-`(unique(bf5.m.no0.smooth$period5)) %>%
      `rownames<-`(seq(5, 110, by = 5)) %>%
      reshape2::melt() %>%
      mutate(value = exp(value),
             model="LQ fixh") %>%
      setNames(c("age5", "period5", "value", "model"))
  ) %>%
  mutate(period5 = sprintf("%d-%d",period5, period5+4))


DHS.plot <- bf5.m.no0.smooth %>% mutate(value = adjusted / pyears2, period5 = as.numeric(levels(period5)[period5])) %>%
  select(tips, age5, period5, value) %>%
  mutate(period5 = sprintf("%d-%d", period5, period5 + 4))

ggplot(LQ %>% filter(age5 %in% 10:70)) + geom_line(aes(x = age5, y = value, linetype = model), lwd=1.2, col="blue") +
  geom_point(data = filter(DHS.plot, value!=0), aes(x = age5, y = value, shape = tips), alpha = 0) + 
  geom_point(data = filter(DHS.plot, value!=0, tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=3, stroke=2) +
  geom_point(data = filter(DHS.plot, value!=0, !tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=2, alpha = 0.3) +
  scale_shape_manual(values = 0:14) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  facet_wrap(~period5) + 
  ggtitle(country) + ylab(bquote(""[5]*m[x])) + xlab("5-year age groups") +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35))


#females
h <- just.DHS.f$par.full %>% split(names(.)) %>% .$h_params_f
k <- just.DHS.f$par.full %>% split(names(.)) %>% .$k_params_f
fixh.k <- just.DHS.fixh.f$par.full %>% split(names(.)) %>% .$k_params_f

LQ <- as.matrix(LQcoef.f)[,1:4] %*% rbind(rep(1, length(unique(bf5.f.no0.smooth$period5))),
                                          h,
                                          h^2,
                                          k) %>% 
  `colnames<-`(unique(bf5.f.no0.smooth$period5)) %>%
  `rownames<-`(seq(5, 110, by = 5)) %>%
  reshape2::melt() %>%
  mutate(value = exp(value),
         model="LQ") %>%
  setNames(c("age5", "period5", "value", "model")) %>%
  bind_rows(
    as.matrix(LQcoef.f)[,1:4] %*% rbind(rep(1, length(unique(bf5.f.no0.smooth$period5))),
                                        data.vec.f$h_mean_f,
                                        data.vec.f$h_mean_f^2,
                                        fixh.k) %>% 
      `colnames<-`(unique(bf5.f.no0.smooth$period5)) %>%
      `rownames<-`(seq(5, 110, by = 5)) %>%
      reshape2::melt() %>%
      mutate(value = exp(value),
             model="LQ fixh") %>%
      setNames(c("age5", "period5", "value", "model"))
  ) %>%
  mutate(period5 = sprintf("%d-%d",period5, period5+4))

DHS.plot <- bf5.f.no0.smooth %>% mutate(value = adjusted / pyears2, period5 = as.numeric(levels(period5)[period5])) %>%
  select(tips, age5, period5, value) %>%
  mutate(period5 = sprintf("%d-%d", period5, period5 + 4))

ggplot(LQ %>% filter(age5 %in% 10:70)) + geom_line(aes(x = age5, y = value, linetype = model), lwd=1.2, col="red") +
  geom_point(data = filter(DHS.plot, value!=0), aes(x = age5, y = value, shape = tips), alpha = 0) + 
  geom_point(data = filter(DHS.plot, value!=0, tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=3, stroke=2) +
  geom_point(data = filter(DHS.plot, value!=0, !tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=2, alpha = 0.3) +
  scale_shape_manual(values = 0:14) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  facet_wrap(~period5) + 
  ggtitle(country) + ylab(bquote(""[5]*m[x])) + xlab("5-year age groups") +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35))


#stan####

data.vec.f <- list(LQ_baseline_mx_DX_f = LQ.baseline.DX.ax.smooth,
                   h_DX_f = LQ.baseline.DX.bx.smooth,
                   h2_DX_f = LQ.baseline.DX.cx.smooth,
                   k_DX_f = LQ.baseline.DX.vx.smooth,
                   tp_DX_f = tips.DX.smooth,
                   
                   penal_tp = as(crossprod(diff(diag(15))),"sparseMatrix"),
                   null_penal_tp = as(exp(15)*tcrossprod(c(0,1,1,1,rep(0,11))),"sparseMatrix"),
                   penal_tp_0 = as(tcrossprod(c(1,rep(0,14))),"sparseMatrix"),
                   
                   LQ_baseline_f = as.matrix(LQcoef.f[,1:4]),
                   df = bf5.f.no0.smooth$adjusted,
                   Ef = bf5.f.no0.smooth$pyears2,
                   
                   h_mean_f = igme.5q0.5 %>% filter(Sex=="Female", year5 %in% levels(bf5.f.no0.smooth$period5)) %>% .$child.mort %>% log(),
                   
                   DEF_log_age = log(bf5.f.no0.smooth$age5 + 2),
                   DEF_time = as(model.matrix(adjusted~period5-1,data = bf5.f.no0.smooth), "sparseMatrix")
)

par.vec.f <- list(log_marginal_prec_h = 1,
                  log_marginal_prec_k = 1,
                  log_lambda_tp = 0,
                  log_lambda_tp_0_inflated_sd = 0.3,
                  
                  log_dispersion = 1,
                  
                  h_params_f = igme.5q0.5 %>% filter(Sex=="Female", year5 %in% levels(bf5.m.no0.smooth$period5)) %>% .$child.mort %>% log(), 
                  k_params_f = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                  
                  tp_params = rep(0,15),
                  
                  logit_rho_h = 0,
                  logit_rho_k = 0,
                  h_constant_f = igme.5q0.5 %>% filter(Sex=="Female", year5 %in% levels(bf5.m.no0.smooth$period5)) %>% .$child.mort %>% log(),
                  k_constant_f = 0,
                  
                  log_D = -5,
                  log_E = 3,
                  log_F = 4,
                  
                  h_params_f_innov = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                  k_params_f_innov = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                  
                  rho_h = 0,
                  rho_k = 0,
                  
                  k_param_f = 0
)


m.stan <- make_tmb_obj(data.vec.m, par.vec.m, inner_verbose = TRUE, calc_outputs = TRUE, 
                       random = c("h_params_f",
                                  "k_params_f",
                                  "h_constant_f",
                                  "k_constant_f",
                                  "tp_params",
                                  "log_D",
                                  "log_E",
                                  "log_F",
                                  "log_D_mean",
                                  "log_E_mean",
                                  "log_F_mean"), 
                       DLL="just_DHS_f_LQ_DEF",
                       map = list(log_D = factor(rep(1, length(levels(bf5.m.no0.smooth$period5)))),           
                                  logit_rho_D = factor(NA),
                                  log_D_mean = factor(NA),
                                  log_marginal_prec_D = factor(1),
                                  
                                  log_E = factor(rep(1, length(levels(bf5.m.no0.smooth$period5)))),
                                  logit_rho_E = factor(NA),
                                  log_E_mean = factor(NA),
                                  log_marginal_prec_E = factor(1),
                        
                                  log_F = factor(rep(1, length(levels(bf5.m.no0.smooth$period5)))),
                                  logit_rho_F = factor(NA),
                                  log_F_mean = factor(NA),
                                  log_marginal_prec_F = factor(1),
                                  
                                  h_constant_f = factor(rep(NA,length(levels(bf5.f.no0.smooth$period5)))),
                                  k_constant_f = factor(NA)
                                  )
                       )


f.stan <- make_tmb_obj(data.vec.f, par.vec.f, inner_verbose = TRUE, calc_outputs = TRUE, 
                       random = c("h_params_f",
                                  "k_params_f",
                                  "h_constant_f",
                                  "k_constant_f",
                                  "tp_params"
                                  ), 
                       DLL="just_DHS_f",
                       map = list(h_constant_f = factor(rep(NA,length(levels(bf5.f.no0.smooth$period5)))),
                                  k_constant_f = factor(NA)
                       )
)

f.stan.innov <- make_tmb_obj(data.vec.f, par.vec.f, inner_verbose = TRUE, calc_outputs = TRUE, 
                       random = c("h_params_f_innov",
                                  "k_params_f_innov",
                                  "h_constant_f",
                                  "tp_params"
                       ), 
                       DLL="just_DHS_f_innov",
                       map = list(h_constant_f = factor(rep(NA,length(levels(bf5.f.no0.smooth$period5))))
                       )
)


f.stan.DEF.innov <- make_tmb_obj(data.vec.f, par.vec.f, inner_verbose = TRUE, calc_outputs = TRUE, 
                             random = c("h_params_f_innov",
                                        "k_params_f_innov",
                                        "h_constant_f",
                                        "tp_params",
                                        "log_D",
                                        "log_E",
                                        "log_F"
                             ), 
                             DLL="just_DHS_f_LQ_DEF_innov",
                             map = list(h_constant_f = factor(rep(NA,length(levels(bf5.f.no0.smooth$period5))))
                             )
)

f.stan.DEF.innov.TN <- make_tmb_obj(data.vec.f, par.vec.f, inner_verbose = TRUE, calc_outputs = TRUE, 
                                 random = c("h_params_f_innov",
                                            "k_params_f_innov",
                                            "h_constant_f",
                                            "tp_params",
                                            "log_D",
                                            "log_E",
                                            "log_F"
                                 ), 
                                 DLL="just_DHS_f_LQ_DEF_innov_TN",
                                 map = list(h_constant_f = factor(rep(NA,length(levels(bf5.f.no0.smooth$period5))))
                                 )
)

f.stan.DEF.innov.TN.onek <- make_tmb_obj(data.vec.f, par.vec.f, inner_verbose = TRUE, calc_outputs = TRUE, 
                                         random = c("h_params_f_innov",
                                                    "k_param_f",
                                                    "h_constant_f",
                                                    "tp_params",
                                                    "log_D",
                                                    "log_E",
                                                    "log_F"
                                         ), 
                                         DLL="just_DHS_f_LQ_DEF_innov_TN_onek",
                                         map = list(h_constant_f = factor(rep(NA,length(levels(bf5.f.no0.smooth$period5))))
                                         )
)

f.stan.DEF.innov.1k.varyDEF <- make_tmb_obj(data.vec.f, par.vec.f, inner_verbose = TRUE, calc_outputs = TRUE, 
                                         random = c("h_params_f_innov",
                                                    "k_param_f",
                                                    "tp_params",
                                                    "log_D_innov",
                                                    "log_E_innov",
                                                    "log_F_innov",
                                                    "log_D_mean",
                                                    "log_E_mean",
                                                    "log_F_mean"
                                         ), 
                                         DLL="just_DHS_f_LQ_DEF_innov_onek_varyDEF"
                                         )



fit.stan.f <- tmbstan(f.stan, chains = 1, iter = 2000, control = list(adapt_delta = 0.9, max_treedepth = 12))
fit.stan.f.innov <- tmbstan(f.stan.innov, chains = 1, iter = 2000, control = list(adapt_delta = 0.9, max_treedepth = 12))

fit.stan.f.innov.lp <- tmbstan(f.stan.innov, chains = 1, iter = 2000, laplace=TRUE)
fit.stan.f.DEF.innov <- tmbstan(f.stan.DEF.innov, chains = 1, iter = 2000)
fit.stan.f.DEF.innov.TN <- tmbstan(f.stan.DEF.innov.TN, chains = 1, iter = 2000, lower = c(rep(-Inf, 5), -1, -1), upper = c(rep(Inf, 5), 1, 1), control = list(adapt_delta = 0.95, max_treedepth = 15))
fit.stan.f.DEF.innov.TN.onek <- tmbstan(f.stan.DEF.innov.TN.onek, chains = 1, iter = 2000, lower = c(rep(-Inf, 4), -1), upper = c(rep(Inf, 4), 1), control = list(adapt_delta = 0.95, max_treedepth = 15))


fit.stan.f.1k.varyDEF <- tmbstan(f.stan.DEF.innov.1k.varyDEF, chains = 1, iter = 2000, control = list(adapt_delta = 0.9, max_treedepth = 12))
extract(fit.stan.f.innov)$logit_rho_h


launch_shinystan(fit.stan.f.DEF.innov.TN)

fit.stan.m <- tmbstan(m.stanlaplace = TRUE)