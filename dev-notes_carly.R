use_version()
usethis::use_gpl_license(version = 3, include_future = TRUE)

##Make new functions and run them
use_r("my_summaries")
#now edit my_summaries.R to include the functions you want
load_all()
cvec <- c("Boston", NA, "Brookline", "Brighton", NA, "Boston")
nvec <- c(2022, 2021, NA, 2021, 2021, NA, 2022)
char_summary(cvec)
numeric_summary(nvec, na.rm=TRUE)

##check if your package runs correctly
devtools::check()

##publish the package on Github
use_git() #adds Git support
use_github()

#################################################################
##make tiny sample of the jackson fischer data for dev purposes #
#################################################################

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("imcdatasets")
library(imcdatasets)

##100-image subset of breast dataset
breast.sce <- JacksonFischer_2020_BreastCancer("sce", full_dataset = F)
#colData(breast.sce)
data <- data.frame(breast.sce$patient_id, breast.sce$image_number, breast.sce$cell_id, breast.sce$cell_x, breast.sce$cell_y, breast.sce$cell_metacluster,
                   breast.sce$patient_age, breast.sce$tumor_grade)
names(data) <- c('patient_id','image_number','cell_id','cell_x','cell_y','cell_metacluster',
                 'patient_age','tumor_grade')
data_example <- data[data$patient_id %in% c(1:42),]
##remove patients with unknown clinical subtype
#data <- data[data$tumor_clinical_type %in% c('HR-HER2+','HR+HER2-','HR+HER2+','TripleNeg'),]
save(data_example, file='data_example.RData')





