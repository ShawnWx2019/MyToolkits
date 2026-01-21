#' Calculate BLUP or BLUE Values for Multi-Environment Trials
#'
#' @description
#' Calculates BLUP or BLUE values and returns an S4 object containing results,
#' metadata, and estimated broad-sense heritability.
#'
#' @param x Input data frame.
#' @param Tag Trait name to analyze.
#' @param method Method to use: "blup" (Random effects) or "blue" (Fixed effects).
#' @param year.name Column name for year.
#' @param loc.name Column name for location.
#' @param trait.name Column name for trait.
#' @param line.name Column name for line/genotype.
#' @param value.name Column name for phenotypic value.
#' @param use_interaction Logical. Include GxE interaction terms in the model?
#'
#' @return An object of class \code{\linkS4class{BreedingValue}}.
#' @export
#'
#' @importFrom dplyr %>% rename filter select mutate all_of distinct
#' @importFrom tibble column_to_rownames
#' @importFrom lme4 lmer fixef ranef isSingular VarCorr
#' @importFrom emmeans lsmeans
#' @importFrom stats setNames as.formula na.omit rnorm
#' @importFrom methods new
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create synthetic data for demonstration
#' set.seed(123)
#' n_lines <- 20
#' n_loc <- 3
#' n_year <- 2
#' rep_per_loc <- 2
#'
#' df_test <- expand.grid(
#'   line = paste0("L", 1:n_lines),
#'   location = paste0("Loc", 1:n_loc),
#'   year = paste0("202", 1:n_year),
#'   rep = 1:rep_per_loc
#' )
#'
#' # Simulate Phenotypic Value: G + E + GxE + error
#' # Genetic effect (true value)
#' g_eff <- setNames(rnorm(n_lines, mean = 0, sd = 5), paste0("L", 1:n_lines))
#'
#' # Construct value
#' df_test$value <- 100 +
#'   g_eff[df_test$line] +                 # G
#'   rnorm(nrow(df_test), 0, 10) +         # E + GxE combined noise
#'   rnorm(nrow(df_test), 0, 2)            # Residual
#'
#' df_test$trait_name <- "Yield"
#'
#' # 2. Run BLUP Analysis
#' # This will calculate BLUPs and estimate Heritability
#' res_blup <- get_unbiase_by_year_loc(
#'   x = df_test,
#'   Tag = "Yield",
#'   method = "blup",
#'   trait.name = "trait_name"
#' )
#'
#' # Show the object (Auto-print summary)
#' res_blup
#'
#' # 3. Run BLUE Analysis
#' res_blue <- get_unbiase_by_year_loc(
#'   x = df_test,
#'   Tag = "Yield",
#'   method = "blue",
#'   trait.name = "trait_name"
#' )
#'
#' res_blue
#'
#' # 4. Get the complete table (Results + Raw Data)
#' full_table <- get_complete_table(res_blup)
#' head(full_table)
#'
get_unbiase_by_year_loc <- function(
    x,
    Tag = "FBL",
    method = 'blup',
    year.name = 'year',
    loc.name = "location",
    trait.name = "trait",
    line.name = 'line',
    value.name = "value",
    use_interaction = TRUE
) {

  # --- 1. Data Preprocessing ---
  req_cols <- c(year.name, loc.name, trait.name, line.name, value.name)
  if (!all(req_cols %in% names(x))) {
    stop("Error: Input data frame is missing specified columns.")
  }

  # Clean and format data
  df_clean <- x %>%
    filter(.data[[trait.name]] == Tag) %>%
    select(
      year = all_of(year.name),
      location = all_of(loc.name),
      line = all_of(line.name),
      value = all_of(value.name)
    ) %>%
    mutate(
      year = as.factor(year),
      location = as.factor(location),
      line = as.factor(line),
      value = as.numeric(value)
    )

  if (any(is.na(df_clean$value))) {
    warning("NA values detected and removed from analysis.")
    df_clean <- na.omit(df_clean)
  }

  n_year <- length(unique(df_clean$year))
  n_loc  <- length(unique(df_clean$location))
  n_total_obs <- nrow(df_clean)
  n_lines_check <- length(unique(df_clean$line))

  method <- tolower(method)

  # --- 2. Model Construction ---
  # Environmental Effects
  env_effs <- c()
  if (n_year > 1) env_effs <- c(env_effs, "(1|year)")
  if (n_loc > 1)  env_effs <- c(env_effs, "(1|location)")

  # Interaction Effects
  gxe_effs <- c()
  if (use_interaction && n_year > 1) gxe_effs <- c(gxe_effs, "(1|line:year)")
  if (use_interaction && n_loc > 1)  gxe_effs <- c(gxe_effs, "(1|line:location)")

  # Initialize variables
  is_singular <- FALSE
  msg <- "Converged"
  h2_val <- NA_real_
  final_df <- data.frame()
  f_str <- ""

  # --- 3. Calculation Logic ---
  if (method == 'blup') {
    # BLUP: Genotype is Random
    gen_eff <- "(1|line)"
    rhs <- paste(c(gen_eff, env_effs, gxe_effs), collapse = " + ")
    f_str <- paste("value ~", rhs)

    # Run Model
    model <- lmer(as.formula(f_str), data = df_clean)

    # Check Singularity
    is_singular <- lme4::isSingular(model)
    if (is_singular) {
      msg <- "Singular fit (Genetic variance approx. 0)"
      h2_val <- 0
    } else {
      # Calculate Heritability
      h2_val <- calculate_h2(model, n_year, n_loc, n_total_obs, n_lines_check)
    }

    # Extract BLUPs
    mu <- fixef(model)["(Intercept)"]
    blup_vals <- ranef(model)$line + mu
    final_df <- blup_vals %>% setNames(Tag)

  } else if (method == 'blue') {
    # BLUE: Genotype is Fixed
    gen_eff <- "line"
    # Note: Environment remains random to account for env variance
    rhs <- paste(c(gen_eff, env_effs), collapse = " + ")
    f_str <- paste("value ~", rhs)

    # Run Model
    model <- lmer(as.formula(f_str), data = df_clean)
    is_singular <- lme4::isSingular(model)

    # Calculate LSMeans
    blue_out <- lsmeans(model, "line") %>% as.data.frame()
    final_df <- blue_out %>%
      select(line, lsmean) %>%
      column_to_rownames("line") %>%
      setNames(Tag)

  } else {
    stop("Invalid method. Please use 'blup' or 'blue'.")
  }

  # --- 4. Return S4 Object ---
  new("BreedingValue",
      result = final_df,
      raw_data = df_clean,
      trait = Tag,
      method = method,
      heritability = h2_val,
      model_info = list(
        is_singular = is_singular,
        message = msg,
        n_year = n_year,
        n_loc = n_loc,
        formula = f_str
      )
  )
}
