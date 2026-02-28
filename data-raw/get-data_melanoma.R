#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("imcdatasets")
library(imcdatasets)

melanoma.sce <- HochSchulz_2022_Melanoma(
  data_type = c("sce"),
  panel = "rna",
  full_dataset = FALSE,
  version = "latest",
  metadata = FALSE,
  on_disk = FALSE,
  h5FilesPath = NULL,
  force = FALSE
)
View(colData(melanoma.sce))
#50 38-channel images

melanoma_data <- data.frame(melanoma.sce$patient_id, melanoma.sce$image_number, melanoma.sce$cell_id, melanoma.sce$cell_x, melanoma.sce$cell_y, melanoma.sce$cell_type,
                   melanoma.sce$patient_age)
names(melanoma_data) <- c('patient_id','image_number','cell_id','cell_x','cell_y','cell_type',
                 'patient_age')

##Remove patients with missing patient_id (there are 2 of them)
melanoma_data <- melanoma_data[!is.na(melanoma_data$patient_id),]

usethis::use_data(melanoma_data, overwrite = TRUE)
