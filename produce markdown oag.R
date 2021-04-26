load("~/more countries final avg sex Rwanda.RData")
skip<-c("Ethiopia","Central African Republic","Comoros","Sao Tome and Principe","Botswana","Cape Verde","Equatorial Guinea","Eritrea","Nigeria (Ondo State)","Ghana","Mauritania","Sudan")
joint.countries<-names(aggr.mat.cohort.0)[!names(aggr.mat.cohort.0)%in%skip]


for(i in joint.countries){
  try(
    rmarkdown::render(
      "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/thiele_oag.rmd",
      params=list(country=i),
      output_file= paste0(gsub(" \\s | ' ","_",i),".pdf"),
      output_dir = "C:/Users/ktang3/Desktop/Imperial/Pop_Construct/Rmd/"
      )
  )
}

for(i in joint.countries){
  print(i)
}