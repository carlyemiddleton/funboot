#' Hoch-Schulz Melanoma Data Subset
#'
#' 7 variables taken from the 2022 \href{https://bodenmillergroup.github.io/imcdatasets/}{Hoch-Schulz melanoma data subset}, and re-named to follow the input format for *funboot.*  Each row of the data corresponds to one cell.
#'
#' @format A data frame with 307931 rows and 7 variables:
#' \describe{
#'   \item{patient_id}{Patient ID}
#'   \item{image_number}{Image Number}
#'   \item{cell_id}{Cell ID}
#'   \item{cell_x}{X coordinate of the cell on the image}
#'   \item{cell_y}{Y coordinate of the cell on the image}
#'   \item{cell_type}{Cell type}
#'   \item{patient_age}{Patient age in years}
#' }
"melanoma_data"
