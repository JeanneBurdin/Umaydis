#########
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
#########
library(stringr)
library(stringi)

#Lecture des fichers SAM
Umaydis.K9_15k <- read.delim("Umaydis.K9_40k-55k.sam", header=FALSE)

#Enlever toutes les cases NA 
Umaydis.K9_15k <- na.omit(Umaydis.K9_15k)

##clip
#SACE
Umaydis.K9_15k$clip <- str_detect(Umaydis.K9_15k$V6, "S")
clip <- vector()
left <- vector()
right <- vector()
for(i in 1:nrow(Umaydis.K9_15k)){
  df <- data.frame(number = stri_remove_empty(str_split(Umaydis.K9_15k[i, "V6"], "[A-Z]+")[[1]]),
                   letter = stri_remove_empty(str_split(Umaydis.K9_15k[i, "V6"], "[0-9]+")[[1]]))
  if(df[1,2] == "S"){
    left <- c(left, df[1,1])
  }else{
    left <- c(left, 0)
  }
  if(df[nrow(df),2] == "S"){
    right <- c(right, df[nrow(df),1])
  }else{
    right <- c(right, 0)
  }
}
Umaydis.K9_15k$left <- left
Umaydis.K9_15k$right <- right

#Export du data frame dans un tableau excel
write_csv(Umaydis.K9_15k,"Umaydis.K9_40k-55kk.csv")
