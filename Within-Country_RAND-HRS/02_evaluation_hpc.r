# =========================================================
# 02_evaluation_hpc.R
# Within-HRS transfer learning across race/ethnicity groups:
# paired evaluation
#
# For each target and source combination, compares:
#   (1) Baseline : ridge logistic (cv.glmnet, alpha = 0), target only
#   (2) Pooled   : ridge logistic, target train + sources,
#                  with group indicators as covariates
#                  (following Lee et al. 2023)
#   (3) Transfer : glmtrans, same source combination
# 70/30 train/test split stratified by outcome; 300 repeats.
#
# Input:  output/01_detection_results.RData
# Output: output/02_evaluation_results.RData   (used by 03)
#         output/02_evaluation_results.xlsx
# =========================================================

library(dplyr)
library(glmtrans)
library(glmnet)
library(pROC)
library(openxlsx)
library(parallel)

out_dir <- "output"

load(file.path(out_dir, "01_detection_results.RData"))


# ---------------------------------------------------------
# 1. Settings
# ---------------------------------------------------------

n_cores   <- 16
n_repeat  <- 300
test_prop <- 0.3

evaluation_seeds <- 1:n_repeat


# ---------------------------------------------------------
# 2. Single-source pairs
# ---------------------------------------------------------

single_source_table <- all_source_combination_table %>%
  filter(N_Source == 1) %>%
  rename(Source = Source_Combination) %>%
  arrange(Target, Source)


# ---------------------------------------------------------
# 3. Helper functions
# ---------------------------------------------------------

safe_auc <- function(y_true, prob) {
  if (length(unique(y_true)) < 2) return(NA_real_)
  as.numeric(pROC::auc(y_true, prob, quiet = TRUE))
}

calc_brier <- function(y_true, prob) mean((prob - y_true)^2)

calc_metrics <- function(y_true, prob, threshold = 0.5) {
  pred <- ifelse(prob >= threshold, 1, 0)
  data.frame(
    ACC   = mean(pred == y_true),
    AUC   = safe_auc(y_true, prob),
    Brier = calc_brier(y_true, prob)
  )
}

# Stratified test-set indices
stratified_test_id <- function(y, prop = test_prop) {
  idx1 <- which(y == 1)
  idx0 <- which(y == 0)
  test1 <- sample(idx1, size = round(prop * length(idx1)))
  test0 <- sample(idx0, size = round(prop * length(idx0)))
  sort(c(test1, test0))
}

# Ridge logistic (lambda.min); used by baseline and pooled
fit_ridge_prob <- function(x_train, y_train, x_test) {
  cvfit <- cv.glmnet(
    x      = x_train,
    y      = y_train,
    family = "binomial",
    alpha  = 0
  )
  as.numeric(
    predict(cvfit, newx = x_test, s = "lambda.min", type = "response")
  )
}

# Pooled training design: target train + sources + group indicators
build_pooled_train <- function(target_name, source_names,
                               target_train_x, target_train_y) {
  
  x_list <- c(list(target_train_x),
              lapply(source_names, function(s) group_glmtrans_clean[[s]]$x))
  y_list <- c(list(target_train_y),
              lapply(source_names, function(s) group_glmtrans_clean[[s]]$y))
  grp_list <- c(
    list(rep(target_name, nrow(target_train_x))),
    lapply(source_names, function(s) group_row_labels[[s]])
  )
  
  X_core <- do.call(rbind, x_list)
  y_pool <- unlist(y_list)
  grp    <- factor(unlist(grp_list))
  
  G <- model.matrix(~ grp)[, -1, drop = FALSE]
  list(x = cbind(X_core, G), y = y_pool, grp_levels = levels(grp))
}

# Pooled test design for target rows
build_pooled_test <- function(te_x, target_name, grp_levels, train_colnames) {
  
  dummy_cols <- paste0("grp", grp_levels[-1])
  G_test <- matrix(0, nrow = nrow(te_x), ncol = length(dummy_cols),
                   dimnames = list(NULL, dummy_cols))
  
  if (target_name %in% grp_levels[-1]) {
    G_test[, paste0("grp", target_name)] <- 1
  }
  
  te_pooled <- cbind(te_x, G_test)
  te_pooled[, train_colnames, drop = FALSE]
}


# ---------------------------------------------------------
# 4. Paired evaluation for one target-source combination
# ---------------------------------------------------------

run_paired_evaluation <- function(target_name,
                                  source_names,
                                  seeds = evaluation_seeds) {
  
  target_data <- group_glmtrans_clean[[target_name]]
  target_n    <- nrow(target_data$x)
  
  results <- lapply(seeds, function(seed) {
    
    set.seed(seed)
    
    test_id  <- stratified_test_id(target_data$y, prop = test_prop)
    train_id <- setdiff(seq_len(target_n), test_id)
    
    tr_x <- target_data$x[train_id, , drop = FALSE]
    tr_y <- target_data$y[train_id]
    te_x <- target_data$x[test_id, , drop = FALSE]
    te_y <- target_data$y[test_id]
    
    # Baseline
    prob_base <- fit_ridge_prob(tr_x, tr_y, te_x)
    m_base    <- calc_metrics(te_y, prob_base)
    
    # Pooled
    pooled    <- build_pooled_train(target_name, source_names, tr_x, tr_y)
    te_pool   <- build_pooled_test(te_x, target_name,
                                   pooled$grp_levels, colnames(pooled$x))
    prob_pool <- fit_ridge_prob(pooled$x, pooled$y, te_pool)
    m_pool    <- calc_metrics(te_y, prob_pool)
    
    # Transfer
    capture.output(
      fit_tr <- glmtrans(
        target = list(x = tr_x, y = tr_y),
        source = group_glmtrans_clean[source_names],
        family = "binomial",
        transfer.source.id = "all"
      )
    )
    prob_tr <- as.numeric(predict(fit_tr, newx = te_x, type = "response"))
    m_tr    <- calc_metrics(te_y, prob_tr)
    
    data.frame(
      Target             = target_name,
      Source_Combination = paste(source_names, collapse = "+"),
      N_Source           = length(source_names),
      Repeat             = seed,
      Train_N            = length(train_id),
      Test_N             = length(test_id),
      Train_Event_N      = sum(tr_y == 1),
      Test_Event_N       = sum(te_y == 1),
      Test_Event_Pct     = mean(te_y == 1) * 100,
      
      Baseline_ACC = m_base$ACC, Baseline_AUC = m_base$AUC, Baseline_Brier = m_base$Brier,
      Pooled_ACC   = m_pool$ACC, Pooled_AUC   = m_pool$AUC, Pooled_Brier   = m_pool$Brier,
      Transfer_ACC = m_tr$ACC,   Transfer_AUC = m_tr$AUC,   Transfer_Brier = m_tr$Brier,
      
      # Positive gain = transfer better
      Gain_vs_Baseline_ACC   = m_tr$ACC     - m_base$ACC,
      Gain_vs_Baseline_AUC   = m_tr$AUC     - m_base$AUC,
      Gain_vs_Baseline_Brier = m_base$Brier - m_tr$Brier,
      
      Gain_vs_Pooled_ACC   = m_tr$ACC     - m_pool$ACC,
      Gain_vs_Pooled_AUC   = m_tr$AUC     - m_pool$AUC,
      Gain_vs_Pooled_Brier = m_pool$Brier - m_tr$Brier,
      
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(results)
}


# ---------------------------------------------------------
# 5. Run all source combinations
# ---------------------------------------------------------

raw_results <- bind_rows(
  mclapply(
    seq_len(nrow(all_source_combination_table)),
    function(i) {
      target_name  <- all_source_combination_table$Target[i]
      source_names <- unlist(strsplit(
        all_source_combination_table$Source_Combination[i], "\\+"))
      run_paired_evaluation(target_name, source_names, evaluation_seeds)
    },
    mc.cores = n_cores
  )
)

raw_results_single <- raw_results %>%
  filter(N_Source == 1) %>%
  rename(Source = Source_Combination)


# ---------------------------------------------------------
# 6. Summary statistics
# ---------------------------------------------------------

summarise_results <- function(df, group_vars) {
  
  ci_mean <- function(v) mean(v, na.rm = TRUE)
  ci_lo   <- function(v) {
    m <- mean(v, na.rm = TRUE); s <- sd(v, na.rm = TRUE); n <- sum(!is.na(v))
    m - 1.96 * s / sqrt(n)
  }
  ci_hi   <- function(v) {
    m <- mean(v, na.rm = TRUE); s <- sd(v, na.rm = TRUE); n <- sum(!is.na(v))
    m + 1.96 * s / sqrt(n)
  }
  pval    <- function(v) tryCatch(t.test(v)$p.value, error = function(e) NA_real_)
  
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      N_Repeat            = n(),
      Mean_Test_Event_Pct = mean(Test_Event_Pct),
      
      Baseline_Mean_AUC = mean(Baseline_AUC, na.rm = TRUE),
      Pooled_Mean_AUC   = mean(Pooled_AUC,   na.rm = TRUE),
      Transfer_Mean_AUC = mean(Transfer_AUC, na.rm = TRUE),
      
      Baseline_Mean_Brier = mean(Baseline_Brier, na.rm = TRUE),
      Pooled_Mean_Brier   = mean(Pooled_Brier,   na.rm = TRUE),
      Transfer_Mean_Brier = mean(Transfer_Brier, na.rm = TRUE),
      
      GvB_AUC_Mean   = ci_mean(Gain_vs_Baseline_AUC),
      GvB_AUC_Lo     = ci_lo(Gain_vs_Baseline_AUC),
      GvB_AUC_Hi     = ci_hi(Gain_vs_Baseline_AUC),
      GvB_AUC_P      = pval(Gain_vs_Baseline_AUC),
      GvB_Brier_Mean = ci_mean(Gain_vs_Baseline_Brier),
      
      GvP_AUC_Mean   = ci_mean(Gain_vs_Pooled_AUC),
      GvP_AUC_Lo     = ci_lo(Gain_vs_Pooled_AUC),
      GvP_AUC_Hi     = ci_hi(Gain_vs_Pooled_AUC),
      GvP_AUC_P      = pval(Gain_vs_Pooled_AUC),
      GvP_Brier_Mean = ci_mean(Gain_vs_Pooled_Brier),
      
      .groups = "drop"
    )
}

summary_all <- summarise_results(
  raw_results, c("Target", "Source_Combination", "N_Source")
) %>% arrange(Target, desc(GvB_AUC_Mean))

summary_single <- summary_all %>% filter(N_Source == 1)


# ---------------------------------------------------------
# 7. Source detection labels (single-source pairs)
# ---------------------------------------------------------

detection_label <- transferability_matrix %>%
  select(Target, Source, Source_Loss, Threshold, Loss_Diff, Transferable)

summary_single_with_detection <- summary_single %>%
  rename(Source = Source_Combination) %>%
  left_join(detection_label, by = c("Target", "Source")) %>%
  arrange(Target, desc(GvB_AUC_Mean))

detection_vs_performance <- summary_single_with_detection %>%
  group_by(Transferable) %>%
  summarise(
    N_Pairs      = n(),
    Mean_GvB_AUC = mean(GvB_AUC_Mean, na.rm = TRUE),
    Mean_GvP_AUC = mean(GvP_AUC_Mean, na.rm = TRUE),
    .groups = "drop"
  )


# ---------------------------------------------------------
# 8. Best combination per target; effect of number of sources
# ---------------------------------------------------------

best_combo <- summary_all %>%
  group_by(Target) %>%
  slice_max(GvB_AUC_Mean, n = 1, with_ties = FALSE) %>%
  ungroup()

source_number_effect <- summary_all %>%
  group_by(Target, N_Source) %>%
  summarise(
    Best_GvB_AUC = max(GvB_AUC_Mean, na.rm = TRUE),
    Best_GvP_AUC = max(GvP_AUC_Mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Target, N_Source)


# ---------------------------------------------------------
# 9. Save
# ---------------------------------------------------------

save(
  raw_results, raw_results_single,
  summary_all, summary_single, summary_single_with_detection,
  best_combo, source_number_effect, detection_vs_performance,
  single_source_table, all_source_combination_table,
  file = file.path(out_dir, "02_evaluation_results.RData")
)

wb <- createWorkbook()
addWorksheet(wb, "Raw_All");               writeData(wb, "Raw_All", raw_results)
addWorksheet(wb, "Summary_All");           writeData(wb, "Summary_All", summary_all)
addWorksheet(wb, "Single_With_Detection"); writeData(wb, "Single_With_Detection", summary_single_with_detection)
addWorksheet(wb, "Best_Per_Target");       writeData(wb, "Best_Per_Target", best_combo)
addWorksheet(wb, "Detection_vs_Perf");     writeData(wb, "Detection_vs_Perf", detection_vs_performance)
addWorksheet(wb, "Source_Number_Effect");  writeData(wb, "Source_Number_Effect", source_number_effect)
addWorksheet(wb, "All_Combinations");      writeData(wb, "All_Combinations", all_source_combination_table)
saveWorkbook(wb, file.path(out_dir, "02_evaluation_results.xlsx"), overwrite = TRUE)