# =========================================================
# Figure3.r
# Figure 3: relative AUC (vs. target-only baseline) of the pooled
# and transfer models, for the best source combination per target.
#
# Input:  Cross_Country/output/02_evaluation_results.xlsx
#         Within-Country_RAND-HRS/output/02_evaluation_results.xlsx
# Output: Figures/Figure3.pdf, Figures/Figure3.png
#
# Run from the repository root.
# =========================================================

library(openxlsx)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(patchwork)

gw_file  <- "Cross_Country/output/02_evaluation_results.xlsx"
hrs_file <- "Within-Country_RAND-HRS/output/02_evaluation_results.xlsx"
fig_dir  <- "Figures"
dir.create(fig_dir, showWarnings = FALSE)


# ---------------------------------------------------------
# 1. Read evaluation results
# ---------------------------------------------------------

gw_best  <- read.xlsx(gw_file,  sheet = "Best_Per_Target")
gw_raw   <- read.xlsx(gw_file,  sheet = "Raw_All")
hrs_best <- read.xlsx(hrs_file, sheet = "Best_Per_Target")
hrs_raw  <- read.xlsx(hrs_file, sheet = "Raw_All")

if (colnames(hrs_raw)[1] != "Target") {
  colnames(hrs_raw)[1] <- "Target"
}


# ---------------------------------------------------------
# 2. Paired AUC ratios across repeats
# ---------------------------------------------------------

compute_ratios <- function(best_df, raw_df) {
  map_dfr(seq_len(nrow(best_df)), function(i) {
    tgt   <- best_df$Target[i]
    combo <- best_df$Source_Combination[i]

    sub <- raw_df %>%
      filter(Target == tgt, Source_Combination == combo)

    pooled_ratio   <- sub$Pooled_AUC   / sub$Baseline_AUC
    transfer_ratio <- sub$Transfer_AUC / sub$Baseline_AUC

    tibble(
      Target = tgt, Source_Combination = combo, N_Repeat = nrow(sub),
      Pooled_est   = mean(pooled_ratio),
      Pooled_lo    = quantile(pooled_ratio, 0.025),
      Pooled_hi    = quantile(pooled_ratio, 0.975),
      Transfer_est = mean(transfer_ratio),
      Transfer_lo  = quantile(transfer_ratio, 0.025),
      Transfer_hi  = quantile(transfer_ratio, 0.975)
    )
  })
}

to_long <- function(ratios_df) {
  bind_rows(
    ratios_df %>% transmute(Target, Method = "Pooled",
                            estimate = Pooled_est, lb_95 = Pooled_lo, ub_95 = Pooled_hi),
    ratios_df %>% transmute(Target, Method = "Transfer",
                            estimate = Transfer_est, lb_95 = Transfer_lo, ub_95 = Transfer_hi)
  )
}

gw_long  <- to_long(compute_ratios(gw_best,  gw_raw))
hrs_long <- to_long(compute_ratios(hrs_best, hrs_raw))

# Display label only
gw_long <- gw_long %>%
  mutate(Target = recode(as.character(Target), "KLOSA" = "KLoSA"))

order_by_transfer <- function(long_df) {
  ord <- long_df %>% filter(Method == "Transfer") %>% arrange(desc(estimate)) %>% pull(Target)
  long_df$Target <- factor(long_df$Target, levels = rev(ord))
  long_df
}

gw_long  <- order_by_transfer(gw_long)
hrs_long <- order_by_transfer(hrs_long)


# ---------------------------------------------------------
# 3. Shared x-axis range
# ---------------------------------------------------------

all_vals <- c(gw_long$lb_95, gw_long$ub_95, hrs_long$lb_95, hrs_long$ub_95)
x_lo <- min(all_vals) * 0.93
x_hi <- max(all_vals) * 1.07
x_breaks <- c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5)
x_breaks <- x_breaks[x_breaks >= x_lo & x_breaks <= x_hi]


# ---------------------------------------------------------
# 4. Forest plot
# ---------------------------------------------------------

col_values <- c("Pooled" = "#6A2C96", "Transfer" = "#2C5D96")

make_forest <- function(df, panel_tag) {
  ggplot(df, aes(x = estimate, y = Target, group = Method, colour = Method)) +
    geom_vline(xintercept = 1, color = "gray50", linewidth = 1, linetype = "dashed") +
    geom_errorbarh(aes(xmin = lb_95, xmax = ub_95), linewidth = 1,
                   position = position_dodge(width = 0.5), height = 0.3) +
    geom_point(aes(shape = Method), fill = "white", size = 1.8, stroke = 1,
               position = position_dodge(width = 0.5)) +
    scale_shape_manual(values = c("Pooled" = 22, "Transfer" = 21)) +
    scale_x_log10(limits = c(x_lo, x_hi), breaks = x_breaks, labels = sprintf("%.1f", x_breaks)) +
    scale_color_manual(values = col_values) +
    labs(x = "Relative AUC (Baseline = 1)", y = "", title = panel_tag) +
    theme_bw() +
    theme(text = element_text(size = 13, face = "bold"),
          plot.title = element_text(size = 10, face = "bold", hjust = 0),
          legend.position = "bottom",
          legend.title = element_blank())
}

p_gw  <- make_forest(gw_long,  "A. International Aging Cohorts")
p_hrs <- make_forest(hrs_long, "B. RAND HRS (Health and Retirement Study)")


# ---------------------------------------------------------
# 5. Combine and save
# ---------------------------------------------------------

figure3 <- (p_gw | p_hrs) +
  plot_layout(ncol = 2, guides = "collect", widths = c(1.15, 1)) &
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "Figure3.pdf"), figure3, width = 9, height = 5, units = "in")
ggsave(file.path(fig_dir, "Figure3.png"), figure3, width = 9, height = 5, units = "in", dpi = 300)
