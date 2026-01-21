#' An S4 class to represent Breeding Value results
#'
#' @slot result A data.frame containing the calculated BLUP or BLUE values.
#' @slot raw_data A data.frame containing the cleaned raw input data.
#' @slot trait A character string indicating the trait name.
#' @slot method A character string indicating the method used ("blup" or "blue").
#' @slot heritability A numeric value representing broad-sense heritability (H2).
#' @slot model_info A list containing model convergence info and variance components.
#' @exportClass BreedingValue
#' @importFrom methods setClass
setClass("BreedingValue",
         slots = c(
           result = "data.frame",
           raw_data = "data.frame",
           trait = "character",
           method = "character",
           heritability = "numeric",
           model_info = "list"
         )
)

#' An S4 class to represent Descriptive Statistics
#'
#' @slot statistics A data.frame containing the summary metrics (Mean, SD, CV, etc.).
#' @slot raw_data A numeric vector containing the input data (cleaned).
#' @slot variable_name A character string indicating the name of the analyzed variable.
#' @exportClass DescriptiveStat
setClass("DescriptiveStat",
         slots = c(
           statistics = "data.frame",
           raw_data = "numeric",
           variable_name = "character"
         )
)
