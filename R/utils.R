# Internal helper function to calculate H2
calculate_h2 <- function(model, n_year, n_loc, n_total_obs, n_line) {

  # 1. Extract Variance Components
  vc <- lme4::VarCorr(model)
  vc_df <- as.data.frame(vc)

  # 2. Get Genotypic Variance (Vg)
  vg <- vc_df[vc_df$grp == "line", "vcov"]
  if (length(vg) == 0) vg <- 0

  # 3. Get Residual Variance (Ve)
  ve <- attr(vc, "sc")^2

  # 4. Get Interaction Variance (V_gxe) - Sum of all line interactions
  # Finds rows where 'grp' contains "line" but is not just "line" (e.g., "line:year")
  vgxe_indices <- grep("line:", vc_df$grp)
  vgxe <- sum(vc_df[vgxe_indices, "vcov"])

  # 5. Calculate Effective Replicates (Approximate for unbalanced data)
  # Total observations / (Num Lines * Num Environments)
  n_env <- max(1, n_year * n_loc) # Avoid division by zero

  # If locations are nested in years or just separate, n_env logic might vary.
  # Here we treat distinct year-loc combinations as total environments if possible,
  # or just use the interaction structure.
  # Simple approach:
  n_rep_eff <- n_total_obs / (n_line * n_env)

  # 6. Calculate H2
  # Formula: Vg / (Vg + Vgxe/n_env + Ve/(n_env * n_rep))
  # Which simplifies to: Vg / (Vg + Vgxe/n_env + Ve/(Total_Obs/n_line))

  denom <- vg + (vgxe / n_env) + (ve / (n_env * n_rep_eff))

  if (denom == 0) return(0)
  h2 <- vg / denom
  return(h2)
}
