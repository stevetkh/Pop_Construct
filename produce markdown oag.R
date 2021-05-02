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


for(i in joint.countries){
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
