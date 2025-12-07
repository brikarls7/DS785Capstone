#  Setup & Load libraries 
library(tidyverse)
library(vip)
library(xgboost)
library(cluster)
library(factoextra)
library(corrplot)
library(skimr)
library(rsample)
library(yardstick)
library(pROC)
library(PRROC)
library(kknn)
library(ranger)
library(broom)
library(officer)
library(cowplot)
library(gridExtra)
library(janitor)
library(ggplot2)
library(dplyr)
library(tibble)

set.seed(123)
dir.create("plots", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)

# Load & clean data 

#If need to reload data:
# pcos <- read_csv("PCOS_data.csv")

# Coerce numeric-like lab columns to numeric (fix beta-HCG issue) 
beta_cols <- names(pcos)[stringr::str_detect(names(pcos), "beta-HCG\\(mIU/mL\\)")]
pcos <- pcos %>% mutate(across(all_of(beta_cols), 
                               ~ readr::parse_number(as.character(.x))))

pcos <- pcos %>% mutate(`AMH(ng/mL)` = as.numeric(`AMH(ng/mL)`))
# LH/FSH ratio creation due to external research
# Source: https://jofem.org/index.php/jofem/article/view/716/284284501
pcos <- pcos %>% mutate(`LH/FSH Ratio` = if_else(is.na(`FSH(mIU/mL)`) | `FSH(mIU/mL)` == 0, NA_real_, `LH(mIU/mL)` / `FSH(mIU/mL)`))

# Train/Test Split (70/30 Stratified)
set.seed(123)
split <- rsample::initial_split(pcos, prop = 0.7, strata = `PCOS (Y/N)`)
train <- rsample::training(split)
test  <- rsample::testing(split)

# Imputation
train_imp <- train %>% 
  mutate(across(where(is.numeric), ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)))

test_imp  <- test %>% 
  mutate(across(where(is.numeric), ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)))

# Baseline model: Logisitc Regression 
baseline_rate <- mean(train_imp$`PCOS (Y/N)` == "Yes")
readr::write_csv(tibble(baseline_positive_rate = baseline_rate), "outputs/baseline_rate.csv")

logit_model <- glm(`PCOS (Y/N)` ~ ., data = train_imp, family = binomial())
saveRDS(logit_model, "outputs/logit_model.rds")

logit_test_probs <- as.numeric(predict(logit_model, newdata = as.data.frame(test_imp), type = "response"))
logit_test_class <- factor(ifelse(logit_test_probs >= 0.5, "Yes", "No"), levels = c("No","Yes"))

#  KNN 
# tune k via 5-fold CV on AUC 
set.seed(123)
train_knn <- train_imp %>%
  mutate(`PCOS (Y/N)` = factor(`PCOS (Y/N)`)) %>%
  select(`PCOS (Y/N)`, where(is.numeric))

test_knn <- test_imp %>%
  mutate(`PCOS (Y/N)` = factor(`PCOS (Y/N)`)) %>%
  select(`PCOS (Y/N)`, where(is.numeric))

# 5-fold CV with stratification
folds  <- vfold_cv(train_knn, v = 5, strata = `PCOS (Y/N)`)

# k values to try
k_grid <- seq(3, 51, by = 2)

pos_level <- levels(train_knn$`PCOS (Y/N)`)[2]

impute_median <- function(df) {
  df %>%
    mutate(across(
      .cols = where(is.numeric),
      .fns  = ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)
    ))
}

cv_auc_for_k <- function(k) {
  aucs <- purrr::map_dbl(folds$splits, function(s) {
    tr <- analysis(s)   %>% impute_median()
    va <- assessment(s) %>% impute_median()
    
    # Make sure both folds have both classes
    tr$`PCOS (Y/N)` <- factor(tr$`PCOS (Y/N)`, levels = levels(train_knn$`PCOS (Y/N)`))
    va$`PCOS (Y/N)` <- factor(va$`PCOS (Y/N)`, levels = levels(train_knn$`PCOS (Y/N)`))
    
    if (nlevels(droplevels(tr$`PCOS (Y/N)`)) < 2)
      return(NA_real_)
    if (nlevels(droplevels(va$`PCOS (Y/N)`)) < 2)
      return(NA_real_)
    
    # KNN model – kknn does scaling for us
    fit <- kknn::kknn(
      formula = `PCOS (Y/N)` ~ .,
      train   = tr,
      test    = va,
      k       = k,
      distance = 2,
      kernel   = "optimal",
      scale    = TRUE
    )
    probs_pos <- fit$prob[, pos_level]
    as.numeric(pROC::auc(va$`PCOS (Y/N)`, probs_pos))
  })
  
  mean(aucs, na.rm = TRUE)
}

knn_cv <- tibble(k= k_grid,
                 mean_cv_auc = purrr::map_dbl(k_grid, cv_auc_for_k))
readr::write_csv(knn_cv, "outputs/knn_tuning.csv")
k_best <- knn_cv$k[which.max(knn_cv$mean_cv_auc)]

# Final KNN using k best
train_knn$`PCOS (Y/N)` <- factor(train_knn$`PCOS (Y/N)`,
                                 levels = c("No", "Yes"))

knn_fit <- kknn::kknn(
  formula = `PCOS (Y/N)` ~ .,
  train   = train_knn,
  test    = test_knn,
  k       = k_best,
  distance = 2,
  kernel   = "optimal",
  scale    = TRUE   
)

# Predicted probabilities 
knn_test_probs <- knn_fit$prob[, "Yes"]

knn_test_class <- factor(
  ifelse(knn_test_probs >= 0.5, "Yes", "No"),
  levels = c("No", "Yes")
)


#  XGBoost
y_train <- train_imp$`PCOS (Y/N)`
y_test <- test_imp$`PCOS (Y/N)`

train_x <- subset(train_imp, select = -`PCOS (Y/N)`) 
x_train_mat <- model.matrix(~ . - 1, data = train_x)
dtrain <- xgb.DMatrix(data = x_train_mat, label = as.numeric(y_train) - 1)

test_x <- subset(test_imp, select = -`PCOS (Y/N)`) 
x_test_mat <- model.matrix(~ . - 1, data = test_x)
dtest  <- xgb.DMatrix(data = x_test_mat,  label = as.numeric(y_test) - 1)

xgb_grid <- expand.grid(
  eta = c(0.03, 0.07, 0.1),
  max_depth = c(3,4,5),
  nrounds = c(200, 400, 600)
)
xgb_cv_results <- purrr::pmap_df(xgb_grid, function(eta, max_depth, nrounds) {
  cv <- xgb.cv(
    params = list(
      objective = "binary:logistic",
      eval_metric = "auc",
      eta = eta,
      max_depth = max_depth,
      subsample = 0.9,
      colsample_bytree = 0.9
    ),
    data = dtrain,
    nrounds = nrounds,
    nfold = 5,
    verbose = 0,
    stratified = TRUE,
    early_stopping_rounds = 25
  )
  
  tibble(
    eta = eta,
    max_depth = max_depth,
    nrounds = cv$best_iteration,
    cv_auc = max(cv$evaluation_log$test_auc_mean)
  )
})

readr::write_csv(xgb_cv_results, "outputs/xgb_tuning.csv")
xgb_best <- xgb_cv_results %>% arrange(desc(cv_auc))

xgb_best <- xgb_best[1,]
xgb_model <- xgboost::xgboost(
  params = list(
    objective = "binary:logistic", eval_metric = "auc",
    eta = xgb_best$eta, max_depth = xgb_best$max_depth,
    subsample = 0.9, colsample_bytree = 0.9
  ),
  data = dtrain, nrounds = xgb_best$nrounds, verbose = 0
)
saveRDS(xgb_model, "outputs/xgb_model.rds")
xgb_test_probs <- predict(xgb_model, dtest)
xgb_test_class <- factor(ifelse(xgb_test_probs >= 0.5, "Yes", "No"), levels = c("No","Yes"))

# Metrics & Plots
calc_metrics <- function(truth, probs, cls) {
  pr <- PRROC::pr.curve(scores.class0 = probs[truth == "Yes"],
                        scores.class1 = probs[truth == "No"], curve = FALSE)
  tibble(
    ROC_AUC = as.numeric(pROC::auc(truth, probs)),
    PR_AUC  = unname(pr$auc.integral),
    Accuracy = yardstick::accuracy_vec(truth, cls),
    Sensitivity = yardstick::sens_vec(truth, cls),
    Specificity = yardstick::spec_vec(truth, cls)
  )
}

metrics_tbl <- bind_rows(
  calc_metrics(y_test, logit_test_probs, logit_test_class) %>% mutate(Model = "Logistic"),
  calc_metrics(y_test, knn_test_probs,   knn_test_class)   %>% mutate(Model = paste0("KNN (k=", k_best, ")")),
  calc_metrics(y_test, xgb_test_probs,   xgb_test_class)   %>% mutate(Model = "XGBoost (tuned)")
) %>% select(Model, everything()) %>% arrange(desc(ROC_AUC))
readr::write_csv(metrics_tbl, "outputs/test_metrics_summary.csv")

# ROC curves
# Ensure factor levels: negative first, positive second
y_test <- factor(y_test, levels = c("No","Yes"))

roc_pcos <- bind_rows(
  yardstick::roc_curve(
    data  = tibble(y_test = y_test, .pred_Yes = logit_test_probs),
    truth = y_test, .pred_Yes, event_level = "second"
  ) %>% mutate(Model = "Logistic"),
  
  yardstick::roc_curve(
    data  = tibble(y_test = y_test, .pred_Yes = knn_test_probs),
    truth = y_test, .pred_Yes, event_level = "second"
  ) %>% mutate(Model = paste0("KNN (k=", k_best, ")")),
  
  yardstick::roc_curve(
    data  = tibble(y_test = y_test, .pred_Yes = xgb_test_probs),
    truth = y_test, .pred_Yes, event_level = "second"
  ) %>% mutate(Model = "XGBoost")
)

g_roc <- ggplot(roc_pcos, aes(x = 1 - specificity, y = sensitivity, color = Model)) +
  geom_path() + geom_abline(lty = 3) + theme_minimal() +
  labs(title = "ROC Curves (Test)", x = "1 - Specificity", y = "Sensitivity")
ggsave("plots/ROC_curves.png", g_roc, width = 7, height = 5, dpi = 150)


# Confusion matrices
cm_plot <- function(truth, cls, title) {
  pcosc <- tibble(truth, pred = cls)
  cm <- yardstick::conf_mat(pcosc, truth, pred)
  autoplot(cm, type = "heatmap") + ggtitle(title)
}
cm_grid <- cowplot::plot_grid(
  cm_plot(y_test, logit_test_class, "Logistic"),
  cm_plot(y_test, knn_test_class,   paste0("KNN (k=", k_best, ")")),
  # cm_plot(y_test, rf_test_class,    "Random Forest"),
  cm_plot(y_test, xgb_test_class,   "XGBoost"),
  ncol = 2
)
ggsave("plots/Confusion_matrices.png", cm_grid, width = 9, height = 6.5, dpi = 150)

# VIPs
log_coefs <- broom::tidy(logit_model) %>%
  dplyr::filter(term != "(Intercept)") %>%
  mutate(importance = abs(estimate)) %>%
  arrange(desc(importance)) %>% slice_head(n = 15)
g_log_vip <- ggplot(log_coefs, aes(x = reorder(term, importance), y = importance)) +
  geom_col() + coord_flip() + theme_minimal() +
  labs(title = "Variable Importance — Logistic (|coef|)", x = "", y = "|Coefficient|")
ggsave("plots/Logistic_VIP.png", g_log_vip, width = 7, height = 6, dpi = 150)

png("plots/XGB_VIP.png", width = 900, height = 700, res = 120)
vip::vip(xgb_model, num_features = 10)
dev.off()

# KNN decision regions on two clinical features
f1 <- if ("AMH" %in% colnames(pcos)) "AMH" else if ("BMI" %in% colnames(pcos)) "BMI" else NULL
f2 <- if ("LH/FSH Ratio" %in% colnames(pcos)) "LH/FSH Ratio" else if ("Age (yrs)" %in% colnames(pcos)) "Age (yrs)" else NULL

if (!is.null(f1) && !is.null(f2)) {
  viz <- pcos %>%
    select(`PCOS (Y/N)`, all_of(c(f1, f2))) %>%
    tidyr::drop_na()
  
  viz_scaled <- viz %>%
    mutate(across(all_of(c(f1, f2)), scale))
  
  gx <- seq(min(viz_scaled[[f1]]), max(viz_scaled[[f1]]), length.out = 250)
  gy <- seq(min(viz_scaled[[f2]]), max(viz_scaled[[f2]]), length.out = 250)

  grid <- expand.grid(setNames(list(gx, gy), c(f1, f2)))
  
  fit2 <- kknn::kknn(`PCOS (Y/N)` ~ ., train = viz_scaled, test = grid, k = k_best)
  
  prob_mat <- as.data.frame(fit2$prob)
  if ("Yes" %in% colnames(prob_mat)) {
    grid$pred_prob <- prob_mat[["Yes"]]
  } else if ("No" %in% colnames(prob_mat)) {
    grid$pred_prob <- 1 - prob_mat[["No"]]
  } else {
    stop("Could not find Yes/No probability columns in kknn output.")
  }
  grid$pred_class <- factor(ifelse(grid$pred_prob >= 0.5, "Yes", "No"), levels = c("No","Yes"))
  
  g_knn <- ggplot() +
    # decision surface by probability
    geom_raster(
      data = grid,
      aes(x = .data[[f1]], y = .data[[f2]], fill = pred_prob),
      alpha = 0.35
    ) +
    # training points (true labels)
    geom_point(
      data = viz_scaled,
      aes(x = .data[[f1]], y = .data[[f2]], shape = `PCOS (Y/N)`, color = `PCOS (Y/N)`),
      size = 1.6, alpha = 0.9
    ) +
    scale_fill_continuous(name = "P(Yes)") +
    theme_minimal() +
    labs(
      title = paste0("KNN Decision Regions (", f1, " vs ", f2, ")"),
      x = f1, y = f2
    )
  
  dir.create("plots", showWarnings = FALSE, recursive = TRUE)
  ggsave("plots/KNN_decision_regions.png", g_knn, width = 7, height = 5, dpi = 150)
}
