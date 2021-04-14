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

compile("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_thiele.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/just_DHS_f_thiele"))

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

country <- "Burkina Faso"

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

thiele_age <- seq(2, 92, by = 5)

data.vec.f <- list(df = bf5.f.no0.smooth$adjusted,
                   Ef = bf5.f.no0.smooth$pyears2,
                   df_age = match(bf5.f.no0.smooth$age5+2, thiele_age),
                   df_time = match(bf5.f.no0.smooth$period5, levels(bf5.f.no0.smooth$period5)),
                   df_tp = c(bf5.f.no0.smooth$tips)-1,
                     
                   log_phi_mean = log(thiele.prior.f[1,]),
                   log_psi_mean = log(thiele.prior.f[2,]),
                   log_lambda_mean = log(thiele.prior.f[3,]),
                   log_delta_mean = log(thiele.prior.f[4,]),
                   log_epsilon_mean = log(thiele.prior.f[5,]),
                   #log_epsilon_mean = rep(log(15), length(levels(bf5.f.no0.smooth$period5))),
                   log_A_mean = log(thiele.prior.f[6,]),
                   log_B_mean = log(thiele.prior.f[7,]),
                   
                   thiele_age = thiele_age,
                   
                   penal_tp = as(crossprod(diff(diag(15))),"sparseMatrix"),
                   null_penal_tp = as(exp(15)*tcrossprod(c(0,1,1,1,rep(0,11))),"sparseMatrix"),
                   penal_tp_0 = as(tcrossprod(c(1,rep(0,14))),"sparseMatrix")
                   )

par.vec.f <- list(log_lambda_tp = 0,
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
                  
                  log_marginal_prec_phi = 1,
                  log_marginal_prec_psi = 1,
                  log_marginal_prec_lambda = 1,
                  log_marginal_prec_delta = 1,
                  log_marginal_prec_epsilon = 1,
                  log_marginal_prec_A = 1,
                  log_marginal_prec_B = 1,
                 
                  logit_rho_phi = 0,
                  logit_rho_psi = 0,
                  logit_rho_lambda = 0,
                  logit_rho_delta = 0,
                  logit_rho_epsilon = 0,
                  logit_rho_A = 0,
                  logit_rho_B = 0
                  )

input.LQ.both.vec.f <- list(data = data.vec.f, par_init = par.vec.f, model = "ccmpp_vr_tmb")

just.DHS.f.1k.DEF.fixmean <- fit_tmb(input.LQ.both.vec.f,inner_verbose=TRUE, random = c("tp_params",
                                                                                        "log_phi_innov",
                                                                                        "log_psi_innov",
                                                                                        "log_lambda_innov",
                                                                                        "log_delta_innov",
                                                                                        "log_epsilon_innov",
                                                                                        "log_A_innov",
                                                                                        "log_B_innov"
                                                                                        ), 
                                     DLL="just_DHS_f_thiele"
                                     )


just.DHS.f.1k.DEF.fixmean.log35 <- fit_tmb(input.LQ.both.vec.f,inner_verbose=TRUE, random = c("tp_params",
                                                                                        "log_phi_innov",
                                                                                        "log_psi_innov",
                                                                                        "log_lambda_innov",
                                                                                        "log_delta_innov",
                                                                                        "log_epsilon_innov",
                                                                                        "log_A_innov",
                                                                                        "log_B_innov"
                                                                                        ), 
                                           DLL="just_DHS_f_thiele"
                                           )

just.DHS.f.1k.DEF.fixmean.log15 <- fit_tmb(input.LQ.both.vec.f,inner_verbose=TRUE, random = c("tp_params",
                                                                                              "log_phi_innov",
                                                                                              "log_psi_innov",
                                                                                              "log_lambda_innov",
                                                                                              "log_delta_innov",
                                                                                              "log_epsilon_innov",
                                                                                              "log_A_innov",
                                                                                              "log_B_innov"
                                                                                              ), 
                                           DLL="just_DHS_f_thiele"
                                           )


mx.mat.func <- function(x){
  age <- seq(2, 112, by = 5)
  
  child.mort <- sapply(x$mode$psi, function(y){exp(-y*age)}) %*% diag(x$mode$phi)
  
  hump.mort <- sapply(x$mode$epsilon, function(y){(y-age)^2}) %*% diag(-x$mode$delta) %>% exp() %*% diag(x$mode$lambda)
  
  old.mort <- sapply(x$mode$B, function(y){exp(y*age)}) %*% diag(x$mode$A)
  
  list(child.mort = child.mort, hump.mort = hump.mort, old.mort = old.mort, mx = child.mort + hump.mort + old.mort)
}

models.list <- list("Epsilon ~ Regressed" = just.DHS.f.1k.DEF.fixmean)

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

thiele <- lapply(models.list, function(y){mx.mat.func(y)$mx %>%
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%  
    `rownames<-`(seq(0, 110, by = 5)) %>%
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
  as_tibble(log(just.DHS.f.1k.DEF.fixmean$mode$phi) - data.vec.f$log_phi_mean) %>% 
    mutate(par = "phi", period5 = sort(unique(bf5.f.no0.smooth$period5))),
  
  as_tibble(log(just.DHS.f.1k.DEF.fixmean$mode$psi) - data.vec.f$log_psi_mean) %>% 
    mutate(par = "psi", period5 = sort(unique(bf5.f.no0.smooth$period5))),
  
  as_tibble(log(just.DHS.f.1k.DEF.fixmean$mode$lambda) - data.vec.f$log_lambda_mean) %>% 
    mutate(par = "lambda", period5 = sort(unique(bf5.f.no0.smooth$period5))),
  
  as_tibble(log(just.DHS.f.1k.DEF.fixmean$mode$delta) - data.vec.f$log_delta_mean) %>% 
    mutate(par = "delta", period5 = sort(unique(bf5.f.no0.smooth$period5))),
  
  as_tibble(log(just.DHS.f.1k.DEF.fixmean$mode$epsilon) - data.vec.f$log_epsilon_mean) %>% 
    mutate(par = "epsilon", period5 = sort(unique(bf5.f.no0.smooth$period5))),
  
  as_tibble(log(just.DHS.f.1k.DEF.fixmean$mode$A) - data.vec.f$log_A_mean) %>% 
    mutate(par = "A", period5 = sort(unique(bf5.f.no0.smooth$period5))),
  
  as_tibble(log(just.DHS.f.1k.DEF.fixmean$mode$B) - data.vec.f$log_B_mean) %>% 
    mutate(par = "B", period5 = sort(unique(bf5.f.no0.smooth$period5)))
  ) %>%
  mutate(period5 = as.numeric(levels(period5)[period5]),
         par = fct_inorder(par))

ggplot(par.deviation) + geom_line(aes(x = period5, y = value, col = par), lwd = 1.2) +
  theme(text = element_text(size=20)) +
  ggtitle("Deviation from prior means")

lapply(models.list, function(y){mx.mat.func(y)$child.mort %>%
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%  
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt()}) %>% 
  map2(names(.), ~ add_column(.x, prior = rep(.y, nrow(.x)), model = "Child")) %>% 
  bind_rows()

thiele.decomp <- lapply(models.list, function(y){mx.mat.func(y)$child.mort %>%
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%  
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt()}) %>% 
  map2(names(.), ~ add_column(.x, prior = rep(.y, nrow(.x)), model = "Child")) %>% 
  bind_rows() %>%
  bind_rows(
    lapply(models.list, function(y){mx.mat.func(y)$hump.mort %>%
        `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%  
        `rownames<-`(seq(0, 110, by = 5)) %>%
        reshape2::melt()}) %>% 
      map2(names(.), ~ add_column(.x, prior = rep(.y, nrow(.x)), model = "Hump")) %>% 
      bind_rows()
  ) %>%
  bind_rows(
    lapply(models.list, function(y){mx.mat.func(y)$old.mort %>%
        `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%  
        `rownames<-`(seq(0, 110, by = 5)) %>%
        reshape2::melt()}) %>% 
      map2(names(.), ~ add_column(.x, prior = rep(.y, nrow(.x)), model = "Senescent")) %>% 
      bind_rows()
  ) %>%
  bind_rows(
    lapply(models.list, function(y){mx.mat.func(y)$mx %>%
        `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%  
        `rownames<-`(seq(0, 110, by = 5)) %>%
        reshape2::melt()}) %>% 
      map2(names(.), ~ add_column(.x, prior = rep(.y, nrow(.x)), model = "Total")) %>% 
      bind_rows()
  ) %>%
  setNames(c("age5", "period5", "value", "prior", "model")) %>%
  mutate(period5 = as.factor(period5),
         age5 = as.numeric(age5),
         model = fct_relevel(model, c("Total", "Child", "Hump", "Senescent")),
         prior = fct_inorder(prior))

ggplot(thiele.decomp) + geom_line(aes(x = age5, y = value, col = period5, linetype = prior), lwd = 1.2) +
  scale_color_manual(values = colorRampPalette(colors = c("blue", "red"))(length(unique(bf5.f.no0.smooth$period5)))) + 
  scale_y_continuous(trans="log", labels = function(x){as.character(round(x,3))}) +
  theme(text = element_text(size=20),
        strip.text = element_text(size=20)) + 
  ggtitle("Burkina Faso females") +
  facet_grid(model~prior, scales="free_y", switch = "y")

#males ####
igme.h.mean.m <- igme.5q0.5 %>% filter(Sex=="Male", year5 %in% levels(bf5.m.no0.smooth$period5)) %>% .$child.mort %>% log()
thiele.prior.m <- sapply(igme.h.mean.m, thiele.prior, sex="male")

thiele_age <- seq(2, 92, by = 5)

data.vec.m <- list(df = bf5.m.no0.smooth$adjusted,
                   Ef = bf5.m.no0.smooth$pyears2,
                   df_age = match(bf5.m.no0.smooth$age5+2, thiele_age),
                   df_time = match(bf5.m.no0.smooth$period5, levels(bf5.m.no0.smooth$period5)),
                   df_tp = c(bf5.m.no0.smooth$tips)-1,
                   
                   log_phi_mean = log(thiele.prior.m[1,]),
                   log_psi_mean = log(thiele.prior.m[2,]),
                   log_lambda_mean = log(thiele.prior.m[3,]),
                   log_delta_mean = log(thiele.prior.m[4,]),
                   log_epsilon_mean = log(thiele.prior.m[5,]),
                   #log_epsilon_mean = rep(log(15), length(levels(bf5.f.no0.smooth$period5))),
                   log_A_mean = log(thiele.prior.m[6,]),
                   log_B_mean = log(thiele.prior.m[7,]),
                   
                   thiele_age = thiele_age,
                   
                   penal_tp = as(crossprod(diff(diag(15))),"sparseMatrix"),
                   null_penal_tp = as(exp(15)*tcrossprod(c(0,1,1,1,rep(0,11))),"sparseMatrix"),
                   penal_tp_0 = as(tcrossprod(c(1,rep(0,14))),"sparseMatrix")
)

par.vec.m <- list(log_lambda_tp = 0,
                  log_lambda_tp_0_inflated_sd = 0.3,
                  tp_params = rep(0,15),
                  
                  log_dispersion = 1,
                  
                  log_phi_innov = rep(0, length(levels(bf5.m.no0.smooth$period5))),
                  log_psi_innov = rep(0, length(levels(bf5.m.no0.smooth$period5))),
                  log_lambda_innov = rep(0, length(levels(bf5.m.no0.smooth$period5))),
                  log_delta_innov = rep(0, length(levels(bf5.m.no0.smooth$period5))),
                  log_epsilon_innov = rep(0, length(levels(bf5.m.no0.smooth$period5))),
                  log_A_innov = rep(0, length(levels(bf5.m.no0.smooth$period5))),
                  log_B_innov = rep(0, length(levels(bf5.m.no0.smooth$period5))),
                  
                  log_marginal_prec_phi = 3,
                  log_marginal_prec_psi = 3,
                  log_marginal_prec_lambda = 3,
                  log_marginal_prec_delta = 3,
                  log_marginal_prec_epsilon = 3,
                  log_marginal_prec_A = 3,
                  log_marginal_prec_B = 3,
                  
                  logit_rho_phi = 0,
                  logit_rho_psi = 0,
                  logit_rho_lambda = 0,
                  logit_rho_delta = 0,
                  logit_rho_epsilon = 0,
                  logit_rho_A = 0,
                  logit_rho_B = 0
)

input.LQ.both.vec.m <- list(data = data.vec.m, par_init = par.vec.m, model = "ccmpp_vr_tmb")

just.DHS.m.1k.DEF.fixmean <- fit_tmb(input.LQ.both.vec.m,inner_verbose=TRUE, random = c("tp_params",
                                                                                        "log_phi_innov",
                                                                                        "log_psi_innov",
                                                                                        "log_lambda_innov",
                                                                                        "log_delta_innov",
                                                                                        "log_epsilon_innov",
                                                                                        "log_A_innov",
                                                                                        "log_B_innov"
                                                                                        ), 
                                     DLL="just_DHS_f_thiele"
                                     )


just.DHS.m.1k.DEF.fixmean.log35 <- fit_tmb(input.LQ.both.vec.m,inner_verbose=TRUE, random = c("tp_params",
                                                                                              "log_phi_innov",
                                                                                              "log_psi_innov",
                                                                                              "log_lambda_innov",
                                                                                              "log_delta_innov",
                                                                                              "log_epsilon_innov",
                                                                                              "log_A_innov",
                                                                                              "log_B_innov"
                                                                                              ), 
                                           DLL="just_DHS_f_thiele"
                                           )

just.DHS.m.1k.DEF.fixmean.log15 <- fit_tmb(input.LQ.both.vec.m,inner_verbose=TRUE, random = c("tp_params",
                                                                                              "log_phi_innov",
                                                                                              "log_psi_innov",
                                                                                              "log_lambda_innov",
                                                                                              "log_delta_innov",
                                                                                              "log_epsilon_innov",
                                                                                              "log_A_innov",
                                                                                              "log_B_innov"
                                                                                              ), 
                                           DLL="just_DHS_f_thiele"
                                           )


models.list <- list("Epsilon ~ Regressed" = just.DHS.m.1k.DEF.fixmean)

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

thiele <- lapply(models.list, function(y){mx.mat.func(y)$mx %>%
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%  
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt()}) %>% 
  map2(names(.), ~ add_column(.x, model = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>%
  setNames(c("age5", "period5", "value", "model")) %>%
  mutate(period5 = sprintf("%d-%d",period5, period5+4)) %>%
  bind_rows(
    reshape2::melt(more.countries.avg[[country]]$mort.m) %>% 
      filter(Var1 %in% seq(2, 72, by = 5), Var2 %in% seq(1982, 2017, by = 5)) %>%
      mutate(age5 = Var1-2,
             period5 = Var2-2) %>%
      #mutate(age5 = 5 * floor(Var1 / 5),
      #       period5 = 5 * floor(Var2 / 5)) %>% 
      #group_by(age5, period5) %>%
      #summarise_at(vars(value), mean) %>%
      mutate(value = exp(value),
             period5 = sprintf("%d-%d",period5, period5+4),
             model = "Spline average")
  ) %>%
  mutate(model = str_wrap(model, 20))

DHS.plot <- bf5.m.no0.smooth %>% mutate(value = adjusted / pyears2, period5 = as.numeric(levels(period5)[period5])) %>%
  select(tips, age5, period5, value) %>%
  mutate(period5 = sprintf("%d-%d", period5, period5 + 4))

ggplot(thiele %>% filter(age5 %in% 10:70)) + geom_line(aes(x = age5, y = value, col=model), lwd=1.2, linetype = 2) +
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
  as_tibble(log(just.DHS.m.1k.DEF.fixmean$mode$phi) - data.vec.m$log_phi_mean) %>% 
    mutate(par = "phi", period5 = sort(unique(bf5.m.no0.smooth$period5))),
  
  as_tibble(log(just.DHS.m.1k.DEF.fixmean$mode$psi) - data.vec.m$log_psi_mean) %>% 
    mutate(par = "psi", period5 = sort(unique(bf5.m.no0.smooth$period5))),
  
  as_tibble(log(just.DHS.m.1k.DEF.fixmean$mode$lambda) - data.vec.m$log_lambda_mean) %>% 
    mutate(par = "lambda", period5 = sort(unique(bf5.m.no0.smooth$period5))),
  
  as_tibble(log(just.DHS.m.1k.DEF.fixmean$mode$delta) - data.vec.m$log_delta_mean) %>% 
    mutate(par = "delta", period5 = sort(unique(bf5.m.no0.smooth$period5))),
  
  as_tibble(log(just.DHS.m.1k.DEF.fixmean$mode$epsilon) - data.vec.m$log_epsilon_mean) %>% 
    mutate(par = "epsilon", period5 = sort(unique(bf5.m.no0.smooth$period5))),
  
  as_tibble(log(just.DHS.m.1k.DEF.fixmean$mode$A) - data.vec.m$log_A_mean) %>% 
    mutate(par = "A", period5 = sort(unique(bf5.m.no0.smooth$period5))),
  
  as_tibble(log(just.DHS.m.1k.DEF.fixmean$mode$B) - data.vec.m$log_B_mean) %>% 
    mutate(par = "B", period5 = sort(unique(bf5.m.no0.smooth$period5)))
) %>%
  mutate(period5 = as.numeric(levels(period5)[period5]),
         par = fct_inorder(par))

ggplot(par.deviation) + geom_line(aes(x = period5, y = value, col = par), lwd = 1.2) +
  theme(text = element_text(size=20)) +
  ggtitle("Deviation from prior means")


thiele.decomp <-
  thiele.decomp <- lapply(models.list, function(y){mx.mat.func(y)$child.mort %>%
      `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%  
      `rownames<-`(seq(0, 110, by = 5)) %>%
      reshape2::melt()}) %>% 
  map2(names(.), ~ add_column(.x, prior = rep(.y, nrow(.x)), model = "Child")) %>% 
  bind_rows() %>%
  bind_rows(
    lapply(models.list, function(y){mx.mat.func(y)$hump.mort %>%
        `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%  
        `rownames<-`(seq(0, 110, by = 5)) %>%
        reshape2::melt()}) %>% 
      map2(names(.), ~ add_column(.x, prior = rep(.y, nrow(.x)), model = "Hump")) %>% 
      bind_rows()
  ) %>%
  bind_rows(
    lapply(models.list, function(y){mx.mat.func(y)$old.mort %>%
        `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%  
        `rownames<-`(seq(0, 110, by = 5)) %>%
        reshape2::melt()}) %>% 
      map2(names(.), ~ add_column(.x, prior = rep(.y, nrow(.x)), model = "Senescent")) %>% 
      bind_rows()
  ) %>%
  bind_rows(
    lapply(models.list, function(y){mx.mat.func(y)$mx %>%
        `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%  
        `rownames<-`(seq(0, 110, by = 5)) %>%
        reshape2::melt()}) %>% 
      map2(names(.), ~ add_column(.x, prior = rep(.y, nrow(.x)), model = "Total")) %>% 
      bind_rows()
  ) %>%
  setNames(c("age5", "period5", "value", "prior", "model")) %>%
  mutate(period5 = as.factor(period5),
         age5 = as.numeric(age5),
         model = fct_relevel(model, c("Total", "Child", "Hump", "Senescent")),
         prior = fct_inorder(prior))

ggplot(thiele.decomp) + geom_line(aes(x = age5, y = value, col = period5, linetype = prior), lwd = 1.2) +
  scale_color_manual(values = colorRampPalette(colors = c("blue", "red"))(length(unique(bf5.m.no0.smooth$period5)))) + 
  scale_y_continuous(trans="log", labels = function(x){as.character(round(x, 5))}) +
  theme(text = element_text(size=20),
        strip.text = element_text(size=20)) + 
  ggtitle("Burkina Faso males") +
  facet_grid(model~prior, scales="free_y", switch = "y")
