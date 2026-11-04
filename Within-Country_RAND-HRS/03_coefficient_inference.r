# =========================================================
# 03_coefficient_inference.R
# Within-HRS transfer learning across race/ethnicity groups:
# coefficient inference
#
# For the best source combination per target (from 02):
#   - Point estimate: glmtrans (ridge, alpha = 0) on full data
#   - 95% CI and SE: 1000 bootstrap refits of the target data
#   - Coefficients reported on the log-odds and odds-ratio scales
#   - Wealth reported per $100,000 (rescaled after fitting)
#
# Input:  output/01_detection_results.RData
# Output: output/03_coefficient_inference.xlsx      (log-odds scale)
#         output/03_coefficient_inference_OR.xlsx   (odds-ratio scale)
#         output/03_coefficient_inference_OR.tex
#
# Bootstrap row indices are pre-drawn with a fixed seed and do not
# depend on parallel settings. Random cross-validation folds inside
# glmtrans use parallel RNG streams (L'Ecuyer-CMRG), so results may
# vary slightly with the number of cores and the software environment.
# =========================================================

library(dplyr)
library(tidyr)
library(glmtrans)
library(openxlsx)
library(parallel)
library(kableExtra)

RNGkind("L'Ecuyer-CMRG")

out_dir <- "output"

load(file.path(out_dir, "01_detection_results.RData"))

n_cores <- 16
n_boot  <- 1000


# ---------------------------------------------------------
# 1. Best source combination per target (from 02, Best_Per_Target)
# ---------------------------------------------------------

best_pairs <- list(
  Black    = c("Hispanic", "White"),
  Hispanic = c("White", "Black"),
  Other    = c("White", "Black"),
  White    = c("Hispanic", "Black")
)


# ---------------------------------------------------------
# 2. Variable labels
# ---------------------------------------------------------

label_map <- c(
  intercept                      = "Intercept",
  age                            = "Age",
  female                         = "Female",
  educationGED                   = "Education: GED",
  educationHS_Graduate           = "Education: HS Graduate",
  educationSome_College          = "Education: Some College",
  educationCollege_Plus          = "Education: College+",
  wealth                         = "Wealth (per $100k)",
  quit_smokingFormer             = "Smoking: Former",
  quit_smokingNever              = "Smoking: Never",
  bmi                            = "BMI",
  diabetes_treat_typeOral_Only   = "Diabetes Tx: Oral only",
  diabetes_treat_typeInsulin_Any = "Diabetes Tx: Any insulin",
  medicare                       = "Medicare",
  medicaid                       = "Medicaid",
  private_ins                    = "Private insurance"
)


# ---------------------------------------------------------
# 3. Bootstrap inference for one target-source pair
# ---------------------------------------------------------

bootstrap_inference <- function(target_name, source_names,
                                n_boot = 1000, seed = 99) {
  
  target_data <- group_glmtrans_clean[[target_name]]
  source_data <- group_glmtrans_clean[source_names]
  n_target    <- nrow(target_data$x)
  
  # Point estimate
  set.seed(seed)
  capture.output(
    fit_point <- glmtrans(
      target = target_data,
      source = source_data,
      family = "binomial",
      alpha  = 0,
      transfer.source.id = "all"
    )
  )
  coef_point <- fit_point$beta
  coef_names <- names(coef_point)
  n_coef     <- length(coef_point)
  
  # Pre-drawn bootstrap indices
  set.seed(seed)
  boot_idx_list <- lapply(seq_len(n_boot), function(b)
    sample(seq_len(n_target), size = n_target, replace = TRUE))
  
  # Bootstrap refits; failed fits are recorded as NA
  beta_list <- mclapply(boot_idx_list, function(idx) {
    boot_target <- list(x = target_data$x[idx, , drop = FALSE],
                        y = target_data$y[idx])
    tryCatch(
      suppressWarnings({
        capture.output(
          f <- glmtrans(target = boot_target, source = source_data,
                        family = "binomial", alpha = 0,
                        transfer.source.id = "all")
        )
        f$beta
      }),
      error = function(e) rep(NA_real_, n_coef)
    )
  }, mc.cores = n_cores)
  
  beta_mat <- do.call(rbind, beta_list)
  colnames(beta_mat) <- coef_names
  
  n_used  <- sum(!is.na(beta_mat[, 1]))
  boot_se <- apply(beta_mat, 2, sd,       na.rm = TRUE)
  ci_low  <- apply(beta_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
  ci_high <- apply(beta_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  out <- data.frame(
    Target        = target_name,
    Source        = paste(source_names, collapse = "+"),
    Variable      = coef_names,
    Label         = ifelse(coef_names %in% names(label_map),
                           label_map[coef_names], coef_names),
    Coefficient   = as.numeric(coef_point),
    Boot_SE       = as.numeric(boot_se),
    CI_Low        = as.numeric(ci_low),
    CI_High       = as.numeric(ci_high),
    N_Boot_Used   = n_used,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  # Wealth: per $1 -> per $100,000
  w <- out$Variable == "wealth"
  out[w, c("Coefficient", "Boot_SE", "CI_Low", "CI_High")] <-
    out[w, c("Coefficient", "Boot_SE", "CI_Low", "CI_High")] * 1e5
  
  out$CI_Excludes_0 <- ifelse(out$CI_Low > 0 | out$CI_High < 0, "Yes", "No")
  
  out
}


# ---------------------------------------------------------
# 4. Run all targets
# ---------------------------------------------------------

all_inference <- bind_rows(lapply(names(best_pairs), function(tg) {
  bootstrap_inference(tg, best_pairs[[tg]], n_boot = n_boot)
}))


# ---------------------------------------------------------
# 5. Log-odds scale: side-by-side table and export
# ---------------------------------------------------------

group_n  <- sapply(group_glmtrans_clean, function(d) nrow(d$x))
group_ev <- sapply(group_glmtrans_clean, function(d) sum(d$y == 1))
hdr <- sapply(names(best_pairs), function(tg)
  sprintf("%s (N=%d, events=%d)", tg, group_n[[tg]], group_ev[[tg]]))
names(hdr) <- names(best_pairs)

fmt_cell <- function(b, lo, hi, excl)
  sprintf("%.3f [%.3f, %.3f]%s", b, lo, hi, ifelse(excl == "Yes", "*", ""))

side_by_side <- all_inference %>%
  mutate(Cell = fmt_cell(Coefficient, CI_Low, CI_High, CI_Excludes_0),
         Target_hdr = hdr[Target]) %>%
  select(Label, Target_hdr, Cell) %>%
  pivot_wider(names_from = Target_hdr, values_from = Cell)

var_order <- all_inference %>% filter(Target == "Black") %>% pull(Label)
side_by_side <- side_by_side[match(var_order, side_by_side$Label), ]

wb <- createWorkbook()

addWorksheet(wb, "Side_by_Side")
writeData(wb, "Side_by_Side", side_by_side)
writeData(wb, "Side_by_Side",
          "* 95% bootstrap CI excludes 0 (stability flag, not a p-value). Wealth rescaled to per $100,000.",
          startRow = nrow(side_by_side) + 3, startCol = 1)

addWorksheet(wb, "All_Inference")
writeData(wb, "All_Inference", all_inference)

for (tg in names(best_pairs)) {
  addWorksheet(wb, tg)
  writeData(wb, tg, all_inference %>% filter(Target == tg))
}

saveWorkbook(wb, file.path(out_dir, "03_coefficient_inference.xlsx"),
             overwrite = TRUE)


# ---------------------------------------------------------
# 6. Odds-ratio scale
# ---------------------------------------------------------

or_table <- all_inference %>%
  filter(Variable != "intercept") %>%
  mutate(
    Significant = CI_Excludes_0,
    OR_fmt = sprintf("%.3f", exp(Coefficient)),
    CI_fmt = sprintf("%.3f–%.3f", exp(CI_Low), exp(CI_High))
  )

or_export <- or_table %>%
  select(Target, Source, Variable, Label,
         OR = OR_fmt, `95% CI` = CI_fmt,
         N_Boot_Used, Significant)

wb_or <- createWorkbook()
addWorksheet(wb_or, "All_Inference_OR")
writeData(wb_or, "All_Inference_OR", or_export)
saveWorkbook(wb_or, file.path(out_dir, "03_coefficient_inference_OR.xlsx"),
             overwrite = TRUE)


# ---------------------------------------------------------
# 7. LaTeX table (grouped by target)
# ---------------------------------------------------------

targets <- unique(or_table$Target)

group_labels <- character(0)
group_sizes  <- integer(0)
tex_variable <- character(0)
tex_or_ci    <- character(0)

for (tg in targets) {
  sub <- or_table %>% filter(Target == tg)
  src <- unique(sub$Source)
  
  group_labels <- c(group_labels, sprintf("%s $\\leftarrow$ %s", tg, src))
  group_sizes  <- c(group_sizes, nrow(sub))
  
  tex_variable <- c(tex_variable, gsub("$", "\\$", sub$Label, fixed = TRUE))
  tex_or_ci    <- c(tex_or_ci,
                    sprintf("%s (%s)%s", sub$OR_fmt, sub$CI_fmt,
                            ifelse(sub$Significant == "Yes", "*", "")))
}

tex_df <- data.frame(
  Variable      = tex_variable,
  `OR (95% CI)` = tex_or_ci,
  check.names   = FALSE
)

kbl_out <- kbl(
  tex_df,
  format = "latex", booktabs = TRUE, escape = FALSE, align = "lr",
  caption = "Bootstrap coefficient stability analysis: odds ratios and 95\\% CIs across best target-source pairs (HRS race/ethnicity transfer learning).",
  label = "tab:hrs_coef_or"
) %>%
  kable_styling(latex_options = "hold_position", font_size = 9) %>%
  footnote(symbol = "95\\% bootstrap CI excludes 1 (stability flag, not a p-value).",
           escape = FALSE)

start <- 1
for (i in seq_along(targets)) {
  end <- start + group_sizes[i] - 1
  kbl_out <- pack_rows(kbl_out, group_labels[i], start, end, bold = TRUE)
  start <- end + 1
}

writeLines(as.character(kbl_out),
           file.path(out_dir, "03_coefficient_inference_OR.tex"))