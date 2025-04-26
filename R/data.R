#' Jackson-Fischer Breast Cancer Data Subset
#'
#' Eight variables taken from the 100-image subset of the 2020 [Jackson-Fischer breast cancer dataset](https://bodenmillergroup.github.io/imcdatasets/), and re-named to follow the input format for *funboot.*  Each row of the data corresponds to one cell.
#'
#' @format A data frame with 285851 rows and 8 variables:
#' \describe{
#'   \item{patient_id}{Patient ID}
#'   \item{image_number}{Image Number}
#'   \item{cell_id}{Cell ID}
#'   \item{cell_x}{X coordinate of the cell on the image}
#'   \item{cell_y}{Y coordinate of the cell on the image}
#'   \item{cell_type}{Cell metacluster number}
#'   \item{patient_age}{Patient age in years}
#'   \item{tumor_grade}{Grade of tumor (either 1, 2 or 3)}
#' }
"breastcancer_data_subset_100images"
