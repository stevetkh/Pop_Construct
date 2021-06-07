load("~/cohort smooth 1900-2017.RData")
skip<-c("Ethiopia","Central African Republic","Comoros","Sao Tome and Principe","Botswana","Cape Verde","Equatorial Guinea","Eritrea","Nigeria (Ondo State)","Ghana","Mauritania","Sudan")
joint.countries<-names(aggr.mat.cohort.0)[!names(aggr.mat.cohort.0)%in%skip]

for(i in joint.countries){
  tryCatch({
    rmarkdown::render(
      "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_oag.rmd",
      params=list(country=i),
      output_file= paste0(gsub("\\s|'","_",i),".pdf"),
      output_dir = "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/Rmd_RW/"
      )
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
}


for(i in joint.countries[-11]){
  tryCatch({
    rmarkdown::render(
      "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_oag_RW_only.rmd",
      params=list(country=i),
      output_file= paste0(gsub("\\s|'","_",i),".pdf"),
      output_dir = "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/Rmd_RW_only/"
    )
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
}


setwd("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele spline RData")
for(i in joint.countries[-(1:11)]){
  tryCatch({
    rmarkdown::render(
      "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_oag_RW_only.rmd",
      params=list(country=i,
                  lambda_f = 0.003,
                  lambda_m = 0.003,
                  delta_f = 0.4,
                  delta_m = 0.4,
                  epsilon_f = 25,
                  epsilon_m = 30),
      output_file= paste0(gsub("\\s|'","_",i),".pdf"),
      output_dir = "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/Rmd_RW_only_new/"
    )
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
}


for(i in joint.countries[-(11)]){
  tryCatch({
    rmarkdown::render(
      "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_oag_RW_only_5-74.rmd",
      params=list(country=i,
                  lambda_f = 0.003,
                  lambda_m = 0.003,
                  delta_f = 0.4,
                  delta_m = 0.4,
                  epsilon_f = 25,
                  epsilon_m = 30),
      output_file= paste0(gsub("\\s|'","_",i),".pdf"),
      output_dir = "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/Rmd_RW_only_5-74/"
    )
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
}


for(i in joint.countries[-(11)]){
  tryCatch({
    rmarkdown::render(
      "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_oag_RW_only_5-74_LQinit.rmd",
      params=list(country=i,
                  log_phi_hyperprec = log(2/0.01),
                  log_psi_hyperprec = log(2/0.01),
                  log_lambda_hyperprec = log(2/0.1),
                  log_delta_hyperprec = log(2/0.1),
                  log_epsilon_hyperprec = log(2/0.1),
                  log_A_hyperprec = log(2/0.001),
                  log_B_hyperprec = log(2/0.001)),
      output_file= paste0(gsub("\\s|'","_",i),".pdf"),
      output_dir = "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/Rmd_RW_only_5-74_LQinit/"
    )
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
}

for(i in joint.countries[-(11)]){
  tryCatch({
    rmarkdown::render(
      "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_oag_RW_only_5-74_LQinit_Ad.rmd",
      params=list(country=i,
                  log_phi_hyperprec = log(2/0.01),
                  log_psi_hyperprec = log(2/0.01),
                  log_lambda_hyperprec = log(2/0.1),
                  log_delta_hyperprec = log(2/0.1),
                  log_epsilon_hyperprec = log(2/0.1),
                  log_A_hyperprec = log(2/0.001),
                  log_B_hyperprec = log(2/0.001)),
      output_file= paste0(gsub("\\s|'","_",i),".pdf"),
      output_dir = "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/Rmd_RW_only_5-74_LQinit_Ad/"
    )
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
}

for(i in joint.countries[-(11)]){
  tryCatch({
    rmarkdown::render(
      "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_oag_RW_only_5-74_LQinit_spline_novar.rmd",
      params=list(country=i),
      output_file= paste0(gsub("\\s|'","_",i),"_no_var.pdf"),
      output_dir = "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_spline_RW_reports/"
    )
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
}

for(i in joint.countries[-(11)]){
  tryCatch({
    rmarkdown::render(
      "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_oag_RW_only_5-74_LQinit_spline_1_and_5.rmd",
      params=list(country=i),
      output_file= paste0(gsub("\\s|'","_",i),"_Gumble_wrong_tips_and_priors.pdf"),
      output_dir = "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_RW_Gumbel_1_and_5_wrong_tips_and_priors/"
    )
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
}

for(i in joint.countries[-(1:13)]){
  tryCatch({
    rmarkdown::render(
      "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_oag_RW_only_5-74_LQinit_spline_1_and_5_common_sp.rmd",
      params=list(country=i),
      output_file= paste0(gsub("\\s|'","_",i),"_tau_Gumble.pdf"),
      output_dir = "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_RW_Gumbel_1_and_5/"
    )
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
}


for(i in joint.countries[-11]){
  tryCatch({
    rmarkdown::render(
      "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_oag_RW_only_5-74_LQinit_spline_1_and_5_common_sp.rmd",
      params=list(country=i),
      output_file= paste0(gsub("\\s|'","_",i),"_tau_Gumble.pdf"),
      output_dir = "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_RW_Gumbel_1_and_5_common_AR2/"
    )
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
}

