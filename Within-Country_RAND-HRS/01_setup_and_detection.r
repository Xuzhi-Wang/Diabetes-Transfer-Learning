# =========================================================
# 01_setup_and_detection.R
# Within-HRS transfer learning across race/ethnicity groups:
# data setup and source detection
#
# Input:  data/hrs_2018_cross_sectional.rds (cleaned HRS 2018)
# Output: output/01_detection_results.RData   (used by 02)
#         output/01_detection_summary_hrs_race.xlsx
# =========================================================

library(dplyr)
library(glmtrans)
library(openxlsx)

data_path <- "data/hrs_2018_cross_sectional.rds"
out_dir   <- "output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hrs <- readRDS(data_path)


# ---------------------------------------------------------
# 1. Outcome: 1 = Poor/Fair eyesight; 0 = Good/Very Good/Excellent
# ---------------------------------------------------------

hrs_model <- hrs %>%
  mutate(
    eyesight_binary = case_when(
      eyesight %in% c(1, 2)    ~ 1,
      eyesight %in% c(3, 4, 5) ~ 0,
      TRUE                     ~ NA_real_
    )
  )


# ---------------------------------------------------------
# 2. Predictors (race_eth defines the groups, not a predictor)
# ---------------------------------------------------------

predictors <- c(
  "age", "female", "education", "wealth", "quit_smoking", "bmi",
  "diabetes_treat_type", "medicare", "medicaid", "private_ins"
)

model_vars <- c("eyesight_binary", "race_eth", predictors)

missing_vars <- setdiff(model_vars, names(hrs_model))
if (length(missing_vars) > 0) {
  stop("Missing variables in dataset: ", paste(missing_vars, collapse = ", "))
}


# ---------------------------------------------------------
# 3. Race/ethnicity groups
# ---------------------------------------------------------

race_labels <- c("1" = "Hispanic", "2" = "White", "3" = "Black", "4" = "Other")
group_order <- c("Hispanic", "White", "Black", "Other")

hrs_race <- hrs_model %>%
  select(all_of(model_vars)) %>%
  filter(race_eth %in% c(1, 2, 3, 4)) %>%
  mutate(group = race_labels[as.character(race_eth)])


# ---------------------------------------------------------
# 4. Complete cases
# ---------------------------------------------------------

hrs_race <- hrs_race %>%
  mutate(complete_model_case = complete.cases(select(., all_of(model_vars))))

complete_case_summary <- hrs_race %>%
  group_by(group) %>%
  summarise(
    Original_N      = n(),
    Complete_Case_N = sum(complete_model_case),
    Removed_N       = Original_N - Complete_Case_N,
    Removed_Pct     = round(100 * Removed_N / Original_N, 2),
    .groups = "drop"
  ) %>%
  arrange(match(group, group_order))

hrs_model_cc <- hrs_race %>%
  filter(complete_model_case) %>%
  select(-complete_model_case)


# ---------------------------------------------------------
# 5. Categorical variables
# diabetes_treat_type: insulin only and oral + insulin are
# combined into Insulin_Any.
# ---------------------------------------------------------

hrs_model_cc <- hrs_model_cc %>%
  mutate(
    education = factor(
      education,
      levels = 1:5,
      labels = c("Less_HS", "GED", "HS_Graduate", "Some_College", "College_Plus")
    ),
    quit_smoking = factor(
      quit_smoking,
      levels = 0:2,
      labels = c("Current", "Former", "Never")
    ),
    diabetes_treat_type = factor(
      case_when(
        diabetes_treat_type %in% c(2, 3) ~ "Insulin_Any",
        diabetes_treat_type == 1         ~ "Oral_Only",
        diabetes_treat_type == 0         ~ "No_Treatment",
        TRUE                             ~ NA_character_
      ),
      levels = c("No_Treatment", "Oral_Only", "Insulin_Any")
    )
  )


# ---------------------------------------------------------
# 6. Design matrix
# Built once on the full sample so all groups share the same columns.
# ---------------------------------------------------------

model_formula <- as.formula(
  paste("eyesight_binary ~", paste(predictors, collapse = " + "))
)

X_all <- model.matrix(model_formula, data = hrs_model_cc)[, -1, drop = FALSE]
y_all <- hrs_model_cc$eyesight_binary
grp   <- hrs_model_cc$group


# ---------------------------------------------------------
# 7. Split into glmtrans format by group
# group_row_labels is used by 02 to build the pooled model.
# ---------------------------------------------------------

group_glmtrans_clean <- lapply(group_order, function(g) {
  idx <- grp == g
  list(x = X_all[idx, , drop = FALSE], y = y_all[idx])
})
names(group_glmtrans_clean) <- group_order

group_row_labels <- lapply(group_order, function(g) rep(g, sum(grp == g)))
names(group_row_labels) <- group_order


# ---------------------------------------------------------
# 8. Checks
# ---------------------------------------------------------

x_colnames <- lapply(group_glmtrans_clean, function(d) colnames(d$x))
if (!all(sapply(x_colnames, identical, x_colnames[[1]]))) {
  stop("Column mismatch across groups.")
}

y_values <- unique(unlist(lapply(group_glmtrans_clean, function(d) d$y)))
if (!all(y_values %in% c(0, 1))) {
  stop("Outcome is not binary 0/1.")
}


# ---------------------------------------------------------
# 9. Sample size summary
# ---------------------------------------------------------

sample_summary <- data.frame(
  Group         = names(group_glmtrans_clean),
  N             = sapply(group_glmtrans_clean, function(d) nrow(d$x)),
  Outcome_1_N   = sapply(group_glmtrans_clean, function(d) sum(d$y == 1)),
  Outcome_1_Pct = round(sapply(group_glmtrans_clean, function(d) mean(d$y == 1) * 100), 2),
  Outcome_0_N   = sapply(group_glmtrans_clean, function(d) sum(d$y == 0)),
  Outcome_0_Pct = round(sapply(group_glmtrans_clean, function(d) mean(d$y == 0) * 100), 2)
)


# ---------------------------------------------------------
# 10. Source detection
# Descriptive only; does not restrict source combinations.
# ---------------------------------------------------------

run_source_detection <- function(target_name) {
  set.seed(42)
  source_names <- setdiff(names(group_glmtrans_clean), target_name)
  
  res <- source_detection(
    target = group_glmtrans_clean[[target_name]],
    source = group_glmtrans_clean[source_names],
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
  lapply(names(group_glmtrans_clean), run_source_detection)
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
# 11. All non-empty source combinations
# 3 candidate sources per target: 2^3 - 1 = 7; 28 in total.
# ---------------------------------------------------------

generate_all_source_combinations <- function(target_name) {
  sources <- setdiff(names(group_glmtrans_clean), target_name)
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
  lapply(names(group_glmtrans_clean), generate_all_source_combinations)
)


# ---------------------------------------------------------
# 12. Save
# ---------------------------------------------------------

write.xlsx(
  list(
    Complete_Case_Summary   = complete_case_summary,
    Sample_Summary          = sample_summary,
    Transferability_Matrix  = transferability_matrix,
    Transferable_Sources    = transferable_sources_table,
    Source_Popularity       = source_popularity,
    All_Source_Combinations = all_source_combination_table
  ),
  file = file.path(out_dir, "01_detection_summary_hrs_race.xlsx"),
  overwrite = TRUE
)

save(
  group_glmtrans_clean, group_row_labels,
  predictors, model_formula,
  complete_case_summary, sample_summary,
  transferability_matrix, transferable_sources_table, source_popularity,
  all_source_combination_table,
  file = file.path(out_dir, "01_detection_results.RData")
)