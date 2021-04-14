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
                   
                   sigma_D_mean = 5,
                   sigma_E_mean = 3,
                   sigma_F_mean = 0.25
                   )

par.vec.f <- list(log_marginal_prec_h = 3,
                  log_lambda_tp = 0,
                  log_lambda_tp_0_inflated_sd = 0.3,
                  
                  log_dispersion = 1,
                  
                  h_params_f = igme.5q0.5 %>% filter(Sex=="Female", year5 %in% levels(bf5.m.no0.smooth$period5)) %>% .$child.mort %>% log(), 
                  k_params_f = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                  
                  tp_params = rep(0,15),
                  
                  logit_rho_h = 0,
                  #h_constant_f = igme.5q0.5 %>% filter(Sex=="Female", year5 %in% levels(bf5.m.no0.smooth$period5)) %>% .$child.mort %>% log(),
                  #k_constant_f = 0,
                  
                  log_D = rep(log(0.0005), length(levels(bf5.f.no0.smooth$period5))),
                  log_E = rep(log(25), length(levels(bf5.f.no0.smooth$period5))),
                  log_F = rep(log(35), length(levels(bf5.f.no0.smooth$period5))),
                  
                  #log_D_mean = log(0.0005),
                  #log_E_mean = log(4.2),
                  #log_F_mean = log(35),
                   
                  log_D_mean = log(0.0005),
                  log_E_mean = log(8),
                  log_F_mean = log(35),
                  
                  log_marginal_prec_D = -1,
                  log_marginal_prec_E = -1,
                  log_marginal_prec_F = -1,
                  
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

input.LQ.both.vec.f <- list(data = data.vec.f, par_init = par.vec.f, model = "ccmpp_vr_tmb")

just.DHS.f.1k.DEF.fixmean <- fit_tmb(input.LQ.both.vec.f,inner_verbose=TRUE, random = c("h_params_f_innov",
                                                                                "k_param_f",
                                                                                "tp_params",
                                                                                "log_D_innov",
                                                                                "log_E_innov",
                                                                                "log_F_innov"
                                                                                ), 
                                     DLL="just_DHS_f_LQ_DEF_innov_onek_varyDEF",
                                     map = list(log_D_mean = factor(NA),
                                                log_E_mean = factor(NA),
                                                log_F_mean = factor(NA)
                                                )
                                     )


just.DHS.f.1k.DEF.noprior <- fit_tmb(input.LQ.both.vec.f,inner_verbose=TRUE, random = c("h_params_f_innov",
                                                                                        "k_param_f",
                                                                                        "tp_params",
                                                                                        "log_D_innov",
                                                                                        "log_E_innov",
                                                                                        "log_F_innov"
                                                                                        ), 
                                     DLL="just_DHS_f_LQ_DEF_innov_onek_varyDEF_noprior"
                                     )


just.DHS.f.1k.DEF.vaguemean <- fit_tmb(input.LQ.both.vec.f,inner_verbose=TRUE, random = c("h_params_f_innov",
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


just.DHS.f.1k.DEF.vaguemean.fe <- fit_tmb(input.LQ.both.vec.f,inner_verbose=TRUE, random = c("h_params_f_innov",
                                                                                          "k_param_f",
                                                                                          "tp_params",
                                                                                          "log_D_innov",
                                                                                          "log_E_innov",
                                                                                          "log_F_innov"
                                                                                          ),
                                          DLL="just_DHS_f_LQ_DEF_innov_onek_varyDEF"
                                          )



just.DHS.f.1k.DEF.tightmean <- fit_tmb(input.LQ.both.vec.f,inner_verbose=TRUE, random = c("h_params_f_innov",
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


just.DHS.f.1k.DEF.tightmean.fe <- fit_tmb(input.LQ.both.vec.f,inner_verbose=TRUE, random = c("h_params_f_innov",
                                                                                             "k_param_f",
                                                                                             "tp_params",
                                                                                             "log_D_innov",
                                                                                             "log_E_innov",
                                                                                             "log_F_innov"
                                                                                             ), 
                                          DLL="just_DHS_f_LQ_DEF_innov_onek_varyDEF"
                                          )



mx.mat.func <- function(x, sex){
  logage <- log(seq(2,112, by=5))
  DEF.hump <- sapply(x$mode$log_F, function(y){-(logage-y)^2}) %*%
    diag(exp(x$mode$log_E)) %>%
    exp() %*% diag(exp(x$mode$log_D))
  
  if(sex=="female"){LQcoef <- as.matrix(LQcoef.f[,1:4])} else if(sex=="male"){LQcoef <- as.matrix(LQcoef.m[,1:4])}
  
  LQ <- LQcoef %*% rbind(rep(1, length(x$mode$h_params_f)),
                         x$mode$h_params_f,
                         x$mode$h_params_f^2,
                         x$mode$k_params_f
  ) %>% exp()
  LQ <- rbind(exp(x$mode$h_params_f - log(1 - 0.5*exp(x$mode$h_params_f)) + log(0.2)), LQ)
  
  list(LQ=LQ, DEF = DEF.hump, mx = LQ+DEF.hump)
  
}

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

LQDEF <- bind_rows(
  mx.mat.func(just.DHS.f.1k.DEF.noprior, sex="female")$mx %>% 
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="no prior fixed effects DEF mean"),
  
  mx.mat.func(just.DHS.f.1k.DEF.fixmean, sex="female")$mx %>% 
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="fixed DEF mean"),
  
  mx.mat.func(just.DHS.f.1k.DEF.vaguemean, sex="female")$mx %>% 
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="vague prior random effects DEF mean"),
  
  mx.mat.func(just.DHS.f.1k.DEF.vaguemean.fe, sex="female")$mx %>% 
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="vague prior fixed effects DEF mean"),
  
  mx.mat.func(just.DHS.f.1k.DEF.tightmean, sex="female")$mx %>% 
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="tight prior random effects DEF mean")
  
  #mx.mat.func(just.DHS.f.1k.DEF.tightmean.fe, sex="female")$mx %>% 
  #  `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%
  #  `rownames<-`(seq(0, 110, by = 5)) %>%
  #  reshape2::melt() %>%
  #  mutate( model="tight prior fixed effects DEF mean")
) %>%
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

ggplot(LQDEF %>% filter(age5 %in% 0:70)) + geom_line(aes(x = age5, y = value, col=model), lwd=1.2, linetype = 2) +
  geom_line(data = filter(LQDEF, model == "Spline average"), aes(x = age5, y = value, col = model), lwd=1.2, linetype=1)+
  geom_point(data = filter(DHS.plot, value!=0), aes(x = age5, y = value, shape = tips), alpha = 0) + 
  geom_point(data = filter(DHS.plot, value!=0, tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=3, stroke=2) +
  geom_point(data = filter(DHS.plot, value!=0, !tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=2, alpha = 0.3) +
  scale_shape_manual(values = 0:14) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  facet_wrap(~period5) + 
  ggtitle(country) + ylab(bquote(""[5]*m[x])) + xlab("5-year age groups") +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35))

LQ <- bind_rows(
  mx.mat.func(just.DHS.f.1k.DEF.noprior, sex="female")$LQ %>% 
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="no prior fixed effects DEF mean"),
  
  mx.mat.func(just.DHS.f.1k.DEF.fixmean, sex="female")$LQ %>% 
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="fixed DEF mean"),
  
  mx.mat.func(just.DHS.f.1k.DEF.vaguemean, sex="female")$LQ %>% 
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="vague prior random effects DEF mean"),
  
  mx.mat.func(just.DHS.f.1k.DEF.vaguemean.fe, sex="female")$LQ %>% 
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="vague prior fixed effects DEF mean"),
  
  mx.mat.func(just.DHS.f.1k.DEF.tightmean, sex="female")$LQ %>% 
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="tight prior random effects DEF mean")
) %>%
  setNames(c("age5", "period5", "value", "model")) %>%
  mutate(period5 = sprintf("%d-%d",period5, period5+4),
         model = str_wrap(model, 20))

DEF <- bind_rows(
  mx.mat.func(just.DHS.f.1k.DEF.noprior, sex="female")$DEF %>% 
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="no prior fixed effects DEF mean"),
  
  mx.mat.func(just.DHS.f.1k.DEF.fixmean, sex="female")$DEF %>% 
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="fixed DEF mean"),
  
  mx.mat.func(just.DHS.f.1k.DEF.vaguemean, sex="female")$DEF %>% 
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="vague prior random effects DEF mean"),
  
  mx.mat.func(just.DHS.f.1k.DEF.vaguemean.fe, sex="female")$DEF %>% 
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="vague prior fixed effects DEF mean"),
  
  mx.mat.func(just.DHS.f.1k.DEF.tightmean, sex="female")$DEF %>% 
    `colnames<-`(sort(unique(bf5.f.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="tight prior random effects DEF mean")
) %>%
  setNames(c("age5", "period5", "value", "model")) %>%
  mutate(period5 = sprintf("%d-%d",period5, period5+4),
         model = str_wrap(model, 20))

ggplot(DEF %>% filter(age5 %in% 10:70)) + geom_line(aes(x = age5, y = value, col=model), lwd=1.2, linetype = 2) +
  geom_point(data = filter(DHS.plot, value!=0), aes(x = age5, y = value, shape = tips), alpha = 0) + 
  geom_point(data = filter(DHS.plot, value!=0, tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=3, stroke=2) +
  geom_point(data = filter(DHS.plot, value!=0, !tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=2, alpha = 0.3) +
  scale_shape_manual(values = 0:14) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  facet_wrap(~period5) + 
  ggtitle(country) + ylab(bquote(""[5]*m[x])) + xlab("5-year age groups") +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35))

ggplot(LQ %>% filter(age5 %in% 0:70)) + geom_line(aes(x = age5, y = value, col=model), lwd=1.2, linetype = 2) +
  geom_point(data = filter(DHS.plot, value!=0), aes(x = age5, y = value, shape = tips), alpha = 0) + 
  geom_point(data = filter(DHS.plot, value!=0, tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=3, stroke=2) +
  geom_point(data = filter(DHS.plot, value!=0, !tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=2, alpha = 0.3) +
  scale_shape_manual(values = 0:14) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  facet_wrap(~period5) + 
  ggtitle(country) + ylab(bquote(""[5]*m[x])) + xlab("5-year age groups") +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35))



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

                   h_constant_f = igme.5q0.5 %>% filter(Sex=="Female", year5 %in% levels(bf5.f.no0.smooth$period5)) %>% .$child.mort %>% log(),
                   
                   DEF_log_age = log(bf5.f.no0.smooth$age5 + 2)
)

par.vec.f <- list(log_marginal_prec_h = 2,
                  log_marginal_prec_k = 3,
                  log_lambda_tp = 0,
                  log_lambda_tp_0_inflated_sd = 0.3,
                  
                  log_dispersion = 1,
                  
                  tp_params = rep(0,15),
                  
                  logit_rho_h = 0,
                  logit_rho_k = 0,
                  
                  log_D = log(0.05),
                  log_E = log(4.2),
                  log_F = log(35),
                  
                  h_params_f_innov = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                  k_params_f_innov = rep(0, length(levels(bf5.f.no0.smooth$period5)))
)

input.LQ.both.vec.f <- list(data = data.vec.f, par_init = par.vec.f, model = "ccmpp_vr_tmb")

just.DHS.f.k.1DEF <- fit_tmb(input.LQ.both.vec.f,inner_verbose=TRUE, random = c("h_params_f_innov",
                                                                                "k_params_f_innov",
                                                                                "tp_params"
                                                                                ),
                             DLL="just_DHS_f_LQ_DEF_innov_varyk_oneDEF"
                             #map = list(k_params_f_innov = factor(rep(NA, length(levels(bf5.f.no0.smooth$period5)))))
                             )

#males####
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
                   
                   h_constant_f = igme.5q0.5 %>% filter(Sex=="Male", year5 %in% levels(bf5.m.no0.smooth$period5)) %>% .$child.mort %>% log(),
                   
                   DEF_log_age = log(bf5.m.no0.smooth$age5 + 2),
                   DEF_time = as(model.matrix(adjusted~period5-1,data = bf5.m.no0.smooth), "sparseMatrix"),
                   
                   sigma_D_mean = 5,
                   sigma_E_mean = 3,
                   sigma_F_mean = 0.1
)

par.vec.m <- list(log_marginal_prec_h = 3,
                  log_lambda_tp = 0,
                  log_lambda_tp_0_inflated_sd = 0.3,
                  
                  log_dispersion = 1,
                  
                  h_params_f = igme.5q0.5 %>% filter(Sex=="Male", year5 %in% levels(bf5.m.no0.smooth$period5)) %>% .$child.mort %>% log(), 
                  k_params_f = rep(0, length(levels(bf5.f.no0.smooth$period5))),
                  
                  tp_params = rep(0,15),
                  
                  logit_rho_h = 0,
                  
                  log_D = rep(log(0.0005), length(levels(bf5.m.no0.smooth$period5))),
                  log_E = rep(log(25), length(levels(bf5.m.no0.smooth$period5))),
                  log_F = rep(log(35), length(levels(bf5.m.no0.smooth$period5))),
                  
                  #log_D_mean = log(0.0005),
                  #log_E_mean = log(4.2),
                  #log_F_mean = log(35),
                  
                  log_D_mean = log(0.0005),
                  log_E_mean = log(8),
                  log_F_mean = log(35),
                  
                  log_marginal_prec_D = 1,
                  log_marginal_prec_E = 1,
                  log_marginal_prec_F = 1,
                  
                  logit_rho_D = 0,
                  logit_rho_E = 0,
                  logit_rho_F = 0,
                  
                  h_params_f_innov = rep(0, length(levels(bf5.m.no0.smooth$period5))),
                  k_params_f_innov = rep(0, length(levels(bf5.m.no0.smooth$period5))),
                  
                  log_D_innov = rep(0, length(levels(bf5.m.no0.smooth$period5))),
                  log_E_innov = rep(0, length(levels(bf5.m.no0.smooth$period5))),
                  log_F_innov = rep(0, length(levels(bf5.m.no0.smooth$period5))),
                  k_param_f = 0
)

input.LQ.both.vec.m <- list(data = data.vec.m, par_init = par.vec.m, model = "ccmpp_vr_tmb")


just.DHS.m.1k.DEF.fixmean <- fit_tmb(input.LQ.both.vec.m,inner_verbose=TRUE, random = c("h_params_f_innov",
                                                                                        "k_param_f",
                                                                                        "tp_params",
                                                                                        "log_D_innov",
                                                                                        "log_E_innov",
                                                                                        "log_F_innov"
                                                                                        ), 
                                     DLL="just_DHS_f_LQ_DEF_innov_onek_varyDEF",
                                     map = list(log_D_mean = factor(NA),
                                                log_E_mean = factor(NA),
                                                log_F_mean = factor(NA)
                                                )
)


just.DHS.m.1k.DEF.noprior <- fit_tmb(input.LQ.both.vec.m,inner_verbose=TRUE, random = c("h_params_f_innov",
                                                                                        "k_param_f",
                                                                                        "tp_params",
                                                                                        "log_D_innov",
                                                                                        "log_E_innov",
                                                                                        "log_F_innov"
                                                                                        ), 
                                     DLL="just_DHS_f_LQ_DEF_innov_onek_varyDEF_noprior"
                                     )


just.DHS.m.1k.DEF.vaguemean <- fit_tmb(input.LQ.both.vec.m, inner_verbose=TRUE, random = c("h_params_f_innov",
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

just.DHS.m.1k.DEF.vaguemean.fe <- fit_tmb(input.LQ.both.vec.m,inner_verbose=TRUE, random = c("h_params_f_innov",
                                                                                             "k_param_f",
                                                                                             "tp_params",
                                                                                             "log_D_innov",
                                                                                             "log_E_innov",
                                                                                             "log_F_innov"
                                                                                             ),
                                          DLL="just_DHS_f_LQ_DEF_innov_onek_varyDEF"
                                          )



just.DHS.m.1k.DEF.tightmean <- fit_tmb(input.LQ.both.vec.m,inner_verbose=TRUE, random = c("h_params_f_innov",
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


just.DHS.m.1k.DEF.tightmean.fe <- fit_tmb(input.LQ.both.vec.m,inner_verbose=TRUE, random = c("h_params_f_innov",
                                                                                             "k_param_f",
                                                                                             "tp_params",
                                                                                             "log_D_innov",
                                                                                             "log_E_innov",
                                                                                             "log_F_innov"
                                                                                             ), 
                                          DLL="just_DHS_f_LQ_DEF_innov_onek_varyDEF"
                                          )


mx.mat.func <- function(x, sex){
  logage <- log(seq(2,112, by=5))
  DEF.hump <- sapply(x$mode$log_F, function(y){-(logage-y)^2}) %*%
                diag(exp(x$mode$log_E)) %>%
                exp() %*% diag(exp(x$mode$log_D))
  
  if(sex=="female"){LQcoef <- as.matrix(LQcoef.f[,1:4])} else if(sex=="male"){LQcoef <- as.matrix(LQcoef.m[,1:4])}
  
  LQ <- LQcoef %*% rbind(rep(1, length(x$mode$h_params_f)),
                         x$mode$h_params_f,
                         x$mode$h_params_f^2,
                         x$mode$k_params_f
                         ) %>% exp()
  LQ <- rbind(exp(x$mode$h_params_f - log(1 - 0.5*exp(x$mode$h_params_f)) + log(0.2)), LQ)
  
  list(LQ=LQ, DEF = DEF.hump, mx = LQ+DEF.hump)
  
}

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

LQDEF <- bind_rows(
  mx.mat.func(just.DHS.m.1k.DEF.noprior, sex="male")$mx %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="no prior fixed effects DEF mean"),
  
  mx.mat.func(just.DHS.m.1k.DEF.fixmean, sex="male")$mx %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="fixed DEF mean"),
  
  mx.mat.func(just.DHS.m.1k.DEF.vaguemean, sex="male")$mx %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="vague prior random effects DEF mean"),
  
  mx.mat.func(just.DHS.m.1k.DEF.vaguemean.fe, sex="male")$mx %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="vague prior fixed effects DEF mean"),
  
  mx.mat.func(just.DHS.m.1k.DEF.tightmean, sex="male")$mx %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="tight prior random effects DEF mean"),

  mx.mat.func(just.DHS.m.1k.DEF.tightmean.fe, sex="male")$mx %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="tight prior fixed effects DEF mean")
  ) %>%
  setNames(c("age5", "period5", "value", "model")) %>%
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
  mutate(model = str_wrap(model, 20))

DHS.plot <- bf5.m.no0.smooth %>% mutate(value = adjusted / pyears2, period5 = as.numeric(levels(period5)[period5])) %>%
  select(tips, age5, period5, value) %>%
  mutate(period5 = sprintf("%d-%d", period5, period5 + 4))

ggplot(LQDEF %>% filter(age5 %in% 0:70)) + geom_line(aes(x = age5, y = value, col=model), lwd=1.2, linetype = 2) +
  geom_line(data = filter(LQDEF, model == "Spline average"), aes(x = age5, y = value, col = model), lwd=1.2, linetype=1)+
  geom_point(data = filter(DHS.plot, value!=0), aes(x = age5, y = value, shape = tips), alpha = 0) + 
  geom_point(data = filter(DHS.plot, value!=0, tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=3, stroke=2) +
  geom_point(data = filter(DHS.plot, value!=0, !tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=2, alpha = 0.3) +
  scale_shape_manual(values = 0:14) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  facet_wrap(~period5) + 
  ggtitle(country) + ylab(bquote(""[5]*m[x])) + xlab("5-year age groups") +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35))

LQ <- bind_rows(
  mx.mat.func(just.DHS.m.1k.DEF.noprior, sex="male")$LQ %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="no prior fixed effects DEF mean"),
  
  mx.mat.func(just.DHS.m.1k.DEF.fixmean, sex="male")$LQ %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="fixed DEF mean"),
  
  mx.mat.func(just.DHS.m.1k.DEF.vaguemean, sex="male")$LQ %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="vague prior random effects DEF mean"),
  
  mx.mat.func(just.DHS.m.1k.DEF.vaguemean.fe, sex="male")$LQ %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="vague prior fixed effects DEF mean"),
  
  mx.mat.func(just.DHS.m.1k.DEF.tightmean, sex="male")$LQ %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="tight prior random effects DEF mean"),
  
  mx.mat.func(just.DHS.m.1k.DEF.tightmean.fe, sex="male")$LQ %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="tight prior fixed effects DEF mean")
) %>%
  setNames(c("age5", "period5", "value", "model")) %>%
  mutate(period5 = sprintf("%d-%d",period5, period5+4),
         model = str_wrap(model, 20))

DEF <- bind_rows(
  mx.mat.func(just.DHS.m.1k.DEF.noprior, sex="male")$DEF %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="no prior fixed effects DEF mean"),
  
  mx.mat.func(just.DHS.m.1k.DEF.fixmean, sex="male")$DEF %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="fixed DEF mean"),
  
  mx.mat.func(just.DHS.m.1k.DEF.vaguemean, sex="male")$DEF %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="vague prior random effects DEF mean"),
  
  mx.mat.func(just.DHS.m.1k.DEF.vaguemean.fe, sex="male")$DEF %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="vague prior fixed effects DEF mean"),
  
  mx.mat.func(just.DHS.m.1k.DEF.tightmean, sex="male")$DEF %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="tight prior random effects DEF mean"),
  
  mx.mat.func(just.DHS.m.1k.DEF.tightmean.fe, sex="male")$DEF %>% 
    `colnames<-`(sort(unique(bf5.m.no0.smooth$period5))) %>%
    `rownames<-`(seq(0, 110, by = 5)) %>%
    reshape2::melt() %>%
    mutate( model="tight prior fixed effects DEF mean")
) %>%
  setNames(c("age5", "period5", "value", "model")) %>%
  mutate(period5 = sprintf("%d-%d",period5, period5+4),
         model = str_wrap(model, 20))


ggplot(DEF %>% filter(age5 %in% 10:70, model != "vague prior random\neffects DEF mean")) + geom_line(aes(x = age5, y = value, col=model), lwd=1.2, linetype = 2) +
  geom_point(data = filter(DHS.plot, value!=0), aes(x = age5, y = value, shape = tips), alpha = 0) + 
  geom_point(data = filter(DHS.plot, value!=0, tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=3, stroke=2) +
  geom_point(data = filter(DHS.plot, value!=0, !tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=2, alpha = 0.3) +
  scale_shape_manual(values = 0:14) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  facet_wrap(~period5) + 
  ggtitle(paste(country, "(vague prior+random effects not shown as scale distorted)")) + ylab(bquote(""[5]*m[x])) + xlab("5-year age groups") +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35))


ggplot(LQ %>% filter(age5 %in% 0:70)) + geom_line(aes(x = age5, y = value, col=model), lwd=1.2, linetype = 2) +
  geom_point(data = filter(DHS.plot, value!=0), aes(x = age5, y = value, shape = tips), alpha = 0) + 
  geom_point(data = filter(DHS.plot, value!=0, tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=3, stroke=2) +
  geom_point(data = filter(DHS.plot, value!=0, !tips %in% 1:3), aes(x = age5, y = value, shape = tips), size=2, alpha = 0.3) +
  scale_shape_manual(values = 0:14) +
  scale_y_continuous(trans = "log", labels = function(x){as.character(round(x,3))}) +
  facet_wrap(~period5) + 
  ggtitle(country) + ylab(bquote(""[5]*m[x])) + xlab("5-year age groups") +
  theme(text = element_text(size=25),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 35))

