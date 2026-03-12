##library packages
library(readxl)
library(dplyr)
library(sandwich)
library(lmtest)

## Read data with SUC, ACE consumption and other explanatory factors
test <- read_excel("test.xlsx")
names(test)

## log-linear regression model+ Z-score + HC3 (heteroskedasticity-consistent standard errors)
fit_log_linear <- function(data, outcome, covariates, keep_fit = FALSE) {
  
  ## Select outcome and covariates, then remove missing values
  sub <- data %>%
    dplyr::select(all_of(c(outcome, covariates))) %>%
    na.omit()
  
  ## Check whether log-transformation is applicable
  if (any(sub[[outcome]] <= 0, na.rm = TRUE)) {
    stop("Outcome contains non-positive values; log() not applicable")
  }
  
  ## Log-transform outcome
  y <- log(sub[[outcome]])
  
  ## Standardize covariates using Z-scores
  X_scaled <- scale(sub[, covariates, drop = FALSE])
  X <- as.data.frame(X_scaled)
  colnames(X) <- covariates
  
  ## Store scaling parameters for back-transformation if needed
  x_center <- attr(X_scaled, "scaled:center")
  x_scale  <- attr(X_scaled, "scaled:scale")
  
  ## Build model dataset
  df_model <- cbind(logY = y, X)
  
  ## Fit OLS model
  fit <- lm(logY ~ ., data = df_model)
  
  ## Compute HC3 robust covariance matrix and coefficient tests
  vcov_hc3 <- vcovHC(fit, type = "HC3")
  ct <- coeftest(fit, vcov. = vcov_hc3, df = Inf)
  
  ## Remove intercept
  ct <- ct[rownames(ct) != "(Intercept)", , drop = FALSE]
  
  beta <- ct[, "Estimate"]
  se   <- ct[, "Std. Error"]
  pval <- ct[, "Pr(>|z|)"]
  
  ## Calculate 95% confidence intervals on log scale
  zcrit <- 1.96
  lower <- beta - zcrit * se
  upper <- beta + zcrit * se
  
  ## Back-transform to change ratio scale
  CR      <- exp(beta)
  CR_low  <- exp(lower)
  CR_high <- exp(upper)
  pct     <- (CR - 1) * 100
  
  out_tbl <- tibble::tibble(
    Variable          = names(beta),
    Beta              = as.numeric(beta),
    SE_HC3            = as.numeric(se),
    CR_per_1SD        = as.numeric(CR),
    CR_Lower_95CI     = as.numeric(CR_low),
    CR_Upper_95CI     = as.numeric(CR_high),
    PctChange_per_1SD = as.numeric(pct),
    p_value           = as.numeric(pval),
    N_used            = nrow(sub)
  )
  
  ## Return only the result table if keep_fit = FALSE
  if (!keep_fit) return(out_tbl)
  
  ## Return both result table and fitted objects if keep_fit = TRUE
  list(
    table      = out_tbl,
    fit        = fit,
    vcov_hc3   = vcov_hc3,
    covariates = covariates,
    outcome    = outcome,
    sub_raw    = sub,
    X_z        = X,
    x_center   = x_center,
    x_scale    = x_scale
  )
}

##SUC model
cov_suc <- c("aged.60","high",
             "college.above",
             "D.income.percapita","pop.density",
             "AQI_GDR")


res_suc <- fit_log_linear(test, outcome = "SUC", covariates = cov_suc)

print(res_suc)
##ACE model
cov_ace <- c("aged.60","high",
             "college.above",
             "D.income.percapita","pop.density",
             "AQI_GDR")

res_ace <- fit_log_linear(test, outcome = "ACE", covariates = cov_ace)

print(res_ace)

##forest plot
library(tidyverse)
library(forcats)
library(cowplot)
library(ggplot2)

# =========================
# 0) merge SUC + ACE results table
# =========================
mk_dat <- function(res, sweet) {
  res %>%
    transmute(
      factor    = Variable,
      CR        = CR_per_1SD,
      CI_lower  = CR_Lower_95CI,
      CI_upper  = CR_Upper_95CI,
      p_value   = p_value,
      sweetener = sweet
    ) %>%
    mutate(across(c(CR, CI_lower, CI_upper, p_value), as.numeric))
}

dat_suc <- mk_dat(res_suc, "SUC")
dat_ace <- mk_dat(res_ace, "ACE")
ord <- dat_suc$factor

dat <- bind_rows(dat_suc, dat_ace) %>%
  mutate(
    factor = factor(factor, levels = ord)
  )

# =========================
# 1)labels and headline
# =========================
fmt2 <- function(x) sprintf("%.2f", x)
fmt_p <- function(p) case_when(
  is.na(p) ~ NA_character_,
  p < .001 ~ "<0.001",
  p <  .01 ~ sprintf("%.3f", p),
  TRUE     ~ sprintf("%.2f", p)
)

plot_df <- dat %>%
  mutate(
    CR_lab = paste0(fmt2(CR), " (", fmt2(CI_lower), "-", fmt2(CI_upper), ")"),
    p_lab  = fmt_p(p_value),
    p_star = case_when(
      is.na(p_value) ~ "",
      p_value < 0.01 ~ " **",
      p_value < 0.05 ~ " *",
      TRUE ~ ""
    ),
    p_lab_star = paste0(p_lab, p_star)
  ) %>%
  bind_rows(tibble(
    factor     = factor("Z", levels = c(ord, "Z")),
    CR         = NA_real_,
    CI_lower   = NA_real_,
    CI_upper   = NA_real_,
    p_value    = NA_real_,
    sweetener  = NA_character_,
    CR_lab     = "Consumption Ratio (95% CI)",
    p_lab      = "p-value",
    p_star     = "",
    p_lab_star = "p-value"
  )) %>%
  mutate(
    model = factor(as.character(factor), levels = c(ord, "Z")) |> fct_rev()
  )

# log coordinate

plot_df <- plot_df %>%
  mutate(
    logx    = log(CR),
    log_low = log(CI_lower),
    log_hi  = log(CI_upper)
  )

# =========================
# 2) x-axis range and scale
# =========================
ticks_vals <- c(0.5, 0.75, 1, 1.5, 2)
ticks_log  <- log(ticks_vals)

xr  <- range(plot_df$log_low, plot_df$log_hi, na.rm = TRUE)
pad <- diff(xr) * 0.10
xlim_mid <- c(min(xr[1], min(ticks_log)) - pad, max(xr[2], max(ticks_log)) + pad)

cols <- c("SUC" = "#bd0026", "ACE" = "#2a6fdf")
dodge_w <- 0.60

# =========================
# 3) left：Model
# =========================
left_df <- plot_df %>%
  filter(factor != "Z") %>%
  group_by(factor) %>% slice(1) %>% ungroup() %>%
  bind_rows(plot_df %>% filter(factor == "Z"))

mix_left <- ggplot(left_df, aes(y = model)) +
  geom_text(data = filter(left_df, factor != "Z"),
            aes(x = 0, label = as.character(factor)), hjust = 0) +
  geom_text(data = filter(left_df, factor == "Z"),
            aes(x = 0, label = "Model"), hjust = 0, fontface = "bold") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, .02))) +
  theme_void(base_size = 11) +
  theme(plot.margin = margin(r = 2)) +
  coord_cartesian(clip = "off")

# =========================
# 4) PR text column
# =========================
mix_txt <- ggplot(plot_df, aes(y = model)) +
  geom_text(
    data = filter(plot_df, factor != "Z"),
    aes(x = 0, label = CR_lab, color = sweetener),
    hjust = 0,
    position = position_dodge(width = dodge_w),
    show.legend = FALSE
  ) +
  geom_text(
    data = filter(plot_df, factor == "Z"),
    aes(x = 0, label = "Consumption Ratio (95% CI)"),
    hjust = 0,
    fontface = "bold"
  ) +
  scale_color_manual(values = cols) +
  scale_x_continuous(limits = c(0, 1.8), expand = expansion(mult = c(0, .02))) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(l = 2, r = 2)
  ) +
  coord_cartesian(clip = "off")

# =========================
# 5) middle forest plot
# =========================
mix_mid <- ggplot(filter(plot_df, factor != "Z"), aes(y = model, color = sweetener)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  geom_errorbarh(
    aes(xmin = log_low, xmax = log_hi),
    height = 0,
    position = position_dodge(width = dodge_w),
    linewidth = 1
  ) +
  geom_point(
    aes(x = logx),
    size = 2,
    position = position_dodge(width = dodge_w)
  ) +
  scale_color_manual(values = cols) +
  scale_x_continuous(breaks = ticks_log, labels = ticks_vals, limits = xlim_mid) +
  labs(x = "Consumption ratio per 1 SD (log scale)", y = NULL, color = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.line.x  = element_line(color = "black"),
    legend.position = "top",
    plot.margin = margin(l = 2, r = 2)
  )

# =========================
# 6) p-value colunmn
# =========================
mix_pval <- ggplot(plot_df, aes(y = model)) +
  geom_text(
    data = filter(plot_df, factor != "Z"),
    aes(x = 0, label = p_lab_star, color = sweetener),
    hjust = 0, position = position_dodge(width = dodge_w), show.legend = FALSE
  ) +
  geom_text(
    data = filter(plot_df, factor == "Z"),
    aes(x = 0, label = "p-value"), hjust = 0, fontface = "bold"
  ) +
  scale_color_manual(values = cols) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, .02))) +
  theme_minimal(base_size = 11) +
  theme(axis.title = element_blank(),
        axis.text  = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(l = 2)) +
  coord_cartesian(clip = "off")

# =========================
# 7) Merge all subplots into a complete forest plot
# =========================
mix_forest <- cowplot::plot_grid(
  mix_left, mix_txt, mix_mid, mix_pval,
  ncol = 4,
  rel_widths = c(1.2, 1.3, 3.5, 0.8),
  align = "h",
  axis  = "tb"
)

mix_forest

##partial effect curve
library(dplyr)
library(tibble)
library(lmtest)
library(sandwich)


##partial effect curve function
library(ggplot2)

partial_effect_curve <- function(fit_obj, focal,
                                 grid_n = 100,
                                 grid_q = c(0, 1),
                                 x_axis = c("raw", "z"),
                                 y_scale = c("CR", "Y")) {
  
  x_axis  <- match.arg(x_axis)
  y_scale <- match.arg(y_scale)
  
  if (!(focal %in% fit_obj$covariates)) stop("focal is not in covariates")
  
  fit  <- fit_obj$fit
  V    <- fit_obj$vcov_hc3
  Xz   <- fit_obj$X_z
  sub  <- fit_obj$sub_raw
  
  z_lo <- unname(quantile(Xz[[focal]], probs = grid_q[1], na.rm = TRUE))
  z_hi <- unname(quantile(Xz[[focal]], probs = grid_q[2], na.rm = TRUE))
  z_grid <- seq(z_lo, z_hi, length.out = grid_n)
  
  newX <- as.data.frame(matrix(0, nrow = grid_n, ncol = length(fit_obj$covariates)))# Keep all other variables fixed at 0 (=mean), only let focal change.
  colnames(newX) <- fit_obj$covariates
  newX[[focal]] <- z_grid
  
  # Linear prediction: eta = Xb
  Xmat <- model.matrix(~ ., data = newX)  # add intercept
  eta  <- as.numeric(Xmat %*% coef(fit))
  se_eta <- sqrt(diag(Xmat %*% V %*% t(Xmat)))
  eta_lo <- eta - 1.96 * se_eta
  eta_hi <- eta + 1.96 * se_eta
  eta0 <- as.numeric(coef(fit)[["(Intercept)"]])
  
  if (y_scale == "CR") {
    y    <- exp(eta    - eta0)
    y_lo <- exp(eta_lo - eta0)
    y_hi <- exp(eta_hi - eta0)
    ylab <- "Consumption ratio (relative to mean; CR)"
  } else {
    y    <- exp(eta)
    y_lo <- exp(eta_lo)
    y_hi <- exp(eta_hi)
    ylab <- "Predicted consumption (geometric mean)"
  }
  
  # x-axis uses the original scale (raw data)
  x_raw <- z_grid * unname(fit_obj$x_scale[[focal]]) + unname(fit_obj$x_center[[focal]])
  
  df <- tibble(
    focal   = focal,
    z       = z_grid,
    x_raw   = x_raw,
    y       = y,
    y_lo    = y_lo,
    y_hi    = y_hi
  )
  
  p <- ggplot(df, aes(x = if (x_axis == "raw") x_raw else z, y = y)) +
    geom_ribbon(aes(ymin = y_lo, ymax = y_hi), alpha = 0.20) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(
      x = if (x_axis == "raw") paste0(focal, " (raw scale)") else paste0(focal, " (z-score)"),
      y = ylab,
      title = paste0("Partial (marginal) effect of ", focal, " (others held at mean)")
    ) +
    theme_classic()
  
  list(data = df, plot = p)
}

covs <- c("aged.60","high","D.income.percapita","AQI_GDR","college.above","pop.density")

fit_suc <- fit_log_linear(test, outcome = "SUC", covariates = covs, keep_fit = TRUE)
fit_ace <- fit_log_linear(test, outcome = "ACE", covariates = covs, keep_fit = TRUE)

# take AQI_GDR as example
pe_suc <- partial_effect_curve(fit_suc, focal = "AQI_GDR", x_axis = "raw", y_scale = "CR")
pe_ace <- partial_effect_curve(fit_ace, focal = "AQI_GDR", x_axis = "raw", y_scale = "CR")

# Merge data
df_plot <- bind_rows(
  pe_suc$data %>% mutate(sweetener = "SUC"),
  pe_ace$data %>% mutate(sweetener = "ACE")
)

# rug data
rug_df <- bind_rows(
  fit_suc$sub_raw %>% transmute(x_raw = .data[["AQI_GDR"]], sweetener = "SUC"),
  fit_ace$sub_raw %>% transmute(x_raw = .data[["AQI_GDR"]], sweetener = "ACE")
) %>%
  filter(!is.na(x_raw))

#partial effect curve plot
ggplot(df_plot, aes(x = x_raw, y = y, colour = sweetener, fill = sweetener)) +
  geom_ribbon(aes(ymin = y_lo, ymax = y_hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_rug(
    data = rug_df,
    aes(x = x_raw, colour = sweetener),
    sides = "b", alpha = 0.4,linewidth=1,
    inherit.aes = FALSE
  )+
  scale_colour_manual(values = c(SUC = "#bd0026", ACE = "#2a6fdf")) +
  scale_fill_manual(values   = c(SUC = "#f4b3bf", ACE = "#9EC9F3")) +
  labs(
    x = "Good AQI days (%)",
    y = "Consumption ratio"
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),legend.position = c(0.8,0.8))