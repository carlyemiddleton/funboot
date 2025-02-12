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

######################################################################
##make tiny sample of the jackson fischer data for tutorial purposes #
######################################################################

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
data_example <- data[data$patient_id %in% c(1:42),] #trim down to just the first 10 images
##remove patients with unknown clinical subtype
#data <- data[data$tumor_clinical_type %in% c('HR-HER2+','HR+HER2-','HR+HER2+','TripleNeg'),]
#save(data_example, file='data_example.RData')

usethis::use_data(data_example)
usethis::use_r("data") #create a data.R file in the R subdirectory
#write the data.R file description
devtools::document() #do this to export your functions

###################################
##make preprocess_data() function #
###################################

data("data_example")

sumfun.data <- preprocess_data(data=data_example, from.cell=7, to.cell=3, qc.cellcount.cutoff=20, P=50, perm.yn=T,
                       R=200, inc=1, image.dims=c(0,1000,0,1000), summary.function='L')

pffrmodel <- fit_model(formula=outcome ~ patient_age + tumor_grade + L.obs,
                      data=sumfun.data, spatial.covars = c('L.obs'))



detach("package:phantem", unload = TRUE)







