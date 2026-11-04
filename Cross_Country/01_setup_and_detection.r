# =========================================================
# 01_setup_and_detection.R
# Cross-national transfer learning: data setup and source detection
#
# Input:  data/*_final.xlsx (harmonized cohort files)
# Output: output/01_detection_results.RData   (used by 02)
#         output/01_detection_summary_country.xlsx
# =========================================================

library(readxl)
library(dplyr)
library(glmtrans)
library(openxlsx)

data_dir <- "data"
out_dir  <- "output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ---------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------

file_map <- c(
  HRS    = "HRS_final.xlsx",
  CHARLS = "CHARLS_final.xlsx",
  KLoSA  = "KLoSA_final.xlsx",
  LASI   = "LASI_final.xlsx",
  MHAS   = "MAHS_final.xlsx",
  HAALSI = "SA_final.xlsx",
  SHARE  = "SHARE_final.xlsx",
  ELSA   = "UK_final.xlsx"
)

for (f in file_map) {
  if (!file.exists(file.path(data_dir, f))) {
    stop("Data file not found: ", file.path(data_dir, f))
  }
}

country_list <- lapply(file_map, function(f) read_excel(file.path(data_dir, f)))
names(country_list) <- names(file_map)


# ---------------------------------------------------------
# 2. Outcome: burden_binary = 1 if pain + vision >= 1
# ---------------------------------------------------------

for (country in names(country_list)) {
  country_list[[country]]$burden_score <-
    country_list[[country]]$pain + country_list[[country]]$vision
  country_list[[country]]$burden_binary <-
    ifelse(country_list[[country]]$burden_score == 0, 0, 1)
}


# ---------------------------------------------------------
# 3. Predictors
# ---------------------------------------------------------

predictors <- c(
  "age", "female", "ever_smoke", "diabetes_treat",
  "hypertension_med", "medical_insurance", "social_eco_ppp"
)

model_vars <- c("burden_binary", predictors)

for (country in names(country_list)) {
  missing_vars <- setdiff(model_vars, names(country_list[[country]]))
  if (length(missing_vars) > 0) {
    stop("Missing variables in ", country, ": ",
         paste(missing_vars, collapse = ", "))
  }
}


# ---------------------------------------------------------
# 4. Convert to glmtrans format; drop missing outcome
# ---------------------------------------------------------

country_glmtrans_clean <- lapply(country_list, function(df) {
  idx <- !is.na(df$burden_binary)
  list(
    x = as.matrix(df[idx, predictors]),
    y = df$burden_binary[idx]
  )
})


# ---------------------------------------------------------
# 5. Sample size summary
# ---------------------------------------------------------

sample_summary <- data.frame(
  Country       = names(country_glmtrans_clean),
  N             = sapply(country_glmtrans_clean, function(x) nrow(x$x)),
  Outcome_1_N   = sapply(country_glmtrans_clean, function(x) sum(x$y == 1)),
  Outcome_1_Pct = round(sapply(country_glmtrans_clean, function(x) mean(x$y == 1) * 100), 2),
  Outcome_0_N   = sapply(country_glmtrans_clean, function(x) sum(x$y == 0)),
  Outcome_0_Pct = round(sapply(country_glmtrans_clean, function(x) mean(x$y == 0) * 100), 2)
)


# ---------------------------------------------------------
# 6. Burden score dispersion (descriptive)
# ---------------------------------------------------------

burden_summary <- data.frame(
  Country  = names(country_list),
  Mean     = sapply(country_list, function(df) mean(df$burden_score, na.rm = TRUE)),
  SD       = sapply(country_list, function(df) sd(df$burden_score, na.rm = TRUE)),
  Variance = sapply(country_list, function(df) var(df$burden_score, na.rm = TRUE))
)
burden_summary$Dispersion_Ratio <- burden_summary$Variance / burden_summary$Mean


# ---------------------------------------------------------
# 7. Source detection
# Descriptive only; does not restrict source combinations.
# ---------------------------------------------------------

run_source_detection <- function(target_name) {
  set.seed(42)
  source_names <- setdiff(names(country_glmtrans_clean), target_name)
  
  res <- source_detection(
    target = country_glmtrans_clean[[target_name]],
    source = country_glmtrans_clean[source_names],
    family = "binomial",
    detection.info = FALSE
  )
  
  loss_diff <- res$source.loss - res$threshold
  
  data.frame(
    Target       = target_name,
    Source       = source_names,
    Source_Loss  = as.numeric(res$source.loss),
    Threshold    = as.numeric(res$threshold),
    Loss_Diff    = as.numeric(loss_diff),
    Transferable = as.numeric(loss_diff) < 0,
    stringsAsFactors = FALSE
  )
}

transferability_matrix <- bind_rows(
  lapply(names(country_glmtrans_clean), run_source_detection)
)

transferable_sources_table <- transferability_matrix %>%
  filter(Transferable) %>%
  arrange(Target, Loss_Diff)

source_popularity <- transferable_sources_table %>%
  group_by(Source) %>%
  summarise(
    Accepted_Count = n(),
    Accepted_By    = paste(sort(Target), collapse = ", "),
    .groups        = "drop"
  ) %>%
  arrange(desc(Accepted_Count), Source)


# ---------------------------------------------------------
# 8. All non-empty source combinations
# 7 candidate sources per target: 2^7 - 1 = 127; 1016 in total.
# ---------------------------------------------------------

generate_all_source_combinations <- function(target_name) {
  sources <- setdiff(names(country_glmtrans_clean), target_name)
  combos <- unlist(
    lapply(seq_along(sources), function(k) combn(sources, k, simplify = FALSE)),
    recursive = FALSE
  )
  data.frame(
    Target             = target_name,
    Source_Combination = sapply(combos, paste, collapse = "+"),
    N_Source           = sapply(combos, length),
    stringsAsFactors   = FALSE
  )
}

all_source_combination_table <- bind_rows(
  lapply(names(country_glmtrans_clean), generate_all_source_combinations)
)


# ---------------------------------------------------------
# 9. Save
# ---------------------------------------------------------

write.xlsx(
  list(
    Sample_Summary          = sample_summary,
    Burden_Distribution     = burden_summary,
    Transferability_Matrix  = transferability_matrix,
    Transferable_Sources    = transferable_sources_table,
    Source_Popularity       = source_popularity,
    All_Source_Combinations = all_source_combination_table
  ),
  file = file.path(out_dir, "01_detection_summary_country.xlsx"),
  overwrite = TRUE
)

save(
  country_glmtrans_clean, predictors,
  sample_summary, burden_summary,
  transferability_matrix, transferable_sources_table, source_popularity,
  all_source_combination_table,
  file = file.path(out_dir, "01_detection_results.RData")
)