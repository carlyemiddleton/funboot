#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("imcdatasets")
library(imcdatasets)

##100-image subset of breast cancer dataset
breast.sce <- imcdatasets::JacksonFischer_2020_BreastCancer("sce", full_dataset = F)
View(colData(breast.sce))
range(colData(breast.sce)$cell_x) #note the window dimensions
range(colData(breast.sce)$cell_y)

colData(breast.sce)
breastcancer_data_subset_100images <- data.frame(breast.sce$patient_id, breast.sce$image_number, breast.sce$cell_id, breast.sce$cell_x, breast.sce$cell_y, breast.sce$cell_metacluster,
                   breast.sce$patient_age, breast.sce$tumor_grade)
names(breastcancer_data_subset_100images) <- c('patient_id','image_number','cell_id','cell_x','cell_y','cell_type',
                 'patient_age','tumor_grade')
usethis::use_data(breastcancer_data_subset_100images)

#################################################
##make tiny sample of the data for the vignette #
#################################################
library(usethis)

breastcancer_data_subset_10images <- breastcancer_data_subset_100images[breastcancer_data_subset_100images$patient_id %in% c(1:42),] #trim down to just the first 10 images
#save(breastcancer_data_subset, file='breastcancer_data_subset.RData')
usethis::use_data(breastcancer_data_subset_10images)












