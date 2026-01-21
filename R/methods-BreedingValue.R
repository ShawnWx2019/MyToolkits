#' Show method for BreedingValue object
#'
#' @param object A BreedingValue object.
#' @export
#' @importFrom methods show
setMethod("show", "BreedingValue", function(object) {
  cat("\n=== Breeding Value Analysis (", toupper(object@method), ") ===\n", sep = "")
  cat("Trait       : ", object@trait, "\n", sep = "")

  # Show Heritability if available
  if (!is.na(object@heritability)) {
    cat("Heritability: ", sprintf("%.2f", object@heritability), "\n", sep = "")
  }

  cat("Dimensions  : ", nrow(object@result), " lines across ",
      object@model_info$n_year, " Years, ",
      object@model_info$n_loc, " Locations.\n", sep = "")

  # Warning for Singular fit in BLUP
  if (object@model_info$is_singular && object@method == "blup") {
    cat("\n[WARNING]: Model is Singular (Low Heritability)!\n")
    cat("           Genetic variance is estimated to be near zero.\n")
    cat("           BLUP values have shrunk to the mean.\n")
  } else {
    cat("Status      : ", object@model_info$message, "\n", sep = "")
  }

  cat("\nPreview of Results:\n")
  print(head(object@result, 5))

  if(nrow(object@result) > 5) {
    cat("... (", nrow(object@result)-5, " more lines) \n", sep="")
  }
  cat("\nTip: Use 'get_complete_table(object)' to get full data with raw values.\n")
})


