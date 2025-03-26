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
names(data) <- c('patient_id','image_number','cell_id','cell_x','cell_y','cell_type',
                 'patient_age','tumor_grade')
breastcancer_data_subset <- data[data$patient_id %in% c(1:42),] #trim down to just the first 10 images
##remove patients with unknown clinical subtype
#data <- data[data$tumor_clinical_type %in% c('HR-HER2+','HR+HER2-','HR+HER2+','TripleNeg'),]
#save(data_example, file='data_example.RData')

usethis::use_data(breastcancer_data_subset)
#usethis::use_r("data") #create a data.R file in the R subdirectory
#write the data.R file description
devtools::document() #do this to export your functions
devtools::build_vignettes()

######################################################################
##make tiny sample of the melanoma data for tutorial purposes #
######################################################################

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("imcdatasets")
library(imcdatasets)

##100-image subset of breast dataset
melanoma.sce <- HochSchulz_2022_Melanoma("sce", full_dataset = F)
#View(colData(melanoma.sce))
data <- data.frame(melanoma.sce$patient_id, melanoma.sce$image_number, melanoma.sce$cell_id, melanoma.sce$cell_x, melanoma.sce$cell_y, melanoma.sce$cell_type,
                   melanoma.sce$patient_age, melanoma.sce$patient_gender)
names(data) <- c('patient_id','image_number','cell_id','cell_x','cell_y','cell_type',
                 'patient_age','patient_gender')
melanoma_data <- data[!is.na(data$patient_id),] #trim down to just patients with known IDs
melanoma_data_subset <- melanoma_data[melanoma_data$patient_id %in%
                                        c(64, 36, 84, 51, 71),] #trim to just patients 64, 36, 84, /51, 71

usethis::use_data(melanoma_data_subset)
usethis::use_r("data") #create a data.R file in the R subdirectory
#write the data.R file description
devtools::document() #do this to export your functions
usethis::use_vignette("locacola-vignette")
devtools::load_all()

###################################
##make preprocess_data() function #
###################################
library(spatstat.data)
library(spatstat.univar)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(refund)
library(dplyr)

detach("package:funboot", unload = TRUE)
library(devtools)
devtools::install_github('carlyemiddleton/funboot')
library(funboot)

data("melanoma_data_subset")
data("breastcancer_data_subset")

##preprocess data
sumfun.data <- preprocess_data(data=melanoma_data_subset, from.cell="Macrophage", to.cell="CD8- T cell",
                      qc.cellcount.cutoff=20, P=50, perm.yn=F,
                       R=200, inc=1, image.dims=c(0,1200,0,1200), summary.function='g',seed=456)

#sumfun.data$tumor_grade <- ifelse(sumfun.data$tumor_grade=='1', 1, 0)
sumfun.data$patient_id <- factor(sumfun.data$patient_id)
sumfun.data$patient_gender <- ifelse(sumfun.data$patient_gender=='f', 1, 0)
sumfun.data$patient_age <- ifelse(sumfun.data$patient_age=='<45', 0, 1)

##fit model
#formula <- outcome ~ tumor_grade + patient_age + g.obs
formula <- outcome ~  patient_gender + patient_age + s(patient_id, bs='re')

pffrmodel <- fit_model(formula=formula,spatial.covars = 'g.obs',
                      data=sumfun.data)
summary(pffrmodel)

##Calculate CBs
CBs <- wildBS_CB(formula=formula,
                   data=sumfun.data, spatial.covars = NULL, re=c('patient_id'),
                   B=3,alpha=.05,seed=456)
#Plot CBs
plot.wildBS_CB(CBs)

##Adjust pointwise CIs for multiplicity
p.adjusted <- pointwise_test(data=sumfun.data, formula = formula, r.star = c(50,100), re=c('patient_id'),
                             patient.id=c('patient_id'),alpha=.05,image.id=c('image_number'),
                             spatial.covars = c('g.obs'),method = c('bonferroni'))
p.adjusted

##Run the F test
#formula.full <- outcome ~ tumor_grade #+ patient_age + L.obs #+ s(patient_id, bs='re')
#formula.red <- outcome ~ 1#tumor_grade + patient_age  #+ s(patient_id, bs='re')
formula.full <- outcome ~ patient_gender + patient_age + g.obs + s(patient_id, bs='re')
formula.red <- outcome ~ patient_gender + s(patient_id, bs='re')

F.p.value <- Ftest(formula.full, formula.red,image.id='image_number',patient.id='patient_id',
                  data=sumfun.data, spatial.covars = c('g.obs'),
                  B=50,alpha=.05,re=c('patient_id'),seed=456)
F.p.value


