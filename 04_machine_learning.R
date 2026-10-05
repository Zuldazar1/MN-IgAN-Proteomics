library(data.table)
library(lightgbm)
library(pROC)
ml_proteins <- unique(cox_results$Protein)
ml_proteins <- unique(ml_proteins)
ml_ids <- intersect(
  setdiff(names(proteomics_scaled), "RowName"),
  cox_clinical[!is.na(event), ID]
)
X <- t(as.matrix(proteomics_scaled[match(ml_proteins, RowName),..ml_ids]))
colnames(X) <- ml_proteins
rownames(X) <- ml_ids
storage.mode(X) <- "numeric"

y <- as.integer(cox_clinical[match(ml_ids, ID), event])
clin_ml <- rbindlist(
  list(mn_clinical, igan_clinical),
  use.names = TRUE,
  fill = TRUE
)

clin_use <- clin_ml[
  match(ml_ids, ID),
  .(ID, Age, Gender, eGFR, UTP24h)
]
scale_pos_weight <- sum(y == 0) / sum(y == 1)

lgb_params <- list(
  objective = "binary",
  metric = "auc",
  learning_rate = 0.01,
  num_leaves = 10L,
  max_depth = 15L,
  feature_fraction = 1,
  min_gain_to_split = 0.01,
  scale_pos_weight = scale_pos_weight,
  verbosity = -1,
  seed = 1L
)
nrounds <- 500L

nfold <- 5L

idx_pos <- sample(which(y == 1))
idx_neg <- sample(which(y == 0))
folds <- vector("list", nfold)

for (k in seq_len(nfold)) {
  folds[[k]] <- c(
    idx_pos[seq(k, length(idx_pos), by = nfold)],
    idx_neg[seq(k, length(idx_neg), by = nfold)]
  )
}
##90% Gain####
feat_names <- colnames(X)
gain_sum <- setNames(rep(0, length(feat_names)), feat_names)
for (k in seq_len(nfold)) {
te <- folds[[k]]
  tr <- setdiff(seq_len(nrow(X)), te)
  
  dtrain <- lightgbm::lgb.Dataset(data = X[tr, , drop = FALSE], label = y[tr])
  model <- lightgbm::lgb.train(params = lgb_params, data = dtrain, nrounds = nrounds, verbose = -1)
  imp <- lightgbm::lgb.importance(model, percentage = FALSE)
  
  gain_k <- setNames(rep(0, length(feat_names)), feat_names)
  gain_k[imp$Feature] <- imp$Gain
  gain_k <- gain_k / sum(gain_k)
  
  gain_sum <- gain_sum + gain_k
}
imp_df <- data.table(Protein = feat_names, TotalGain_cv = as.numeric(gain_sum / nfold))
setorder(imp_df, -TotalGain_cv)
imp_df[, CumGain := cumsum(TotalGain_cv)]
imp_df[, CumProp := CumGain / sum(TotalGain_cv)]
k90 <- which(imp_df$CumProp >= 0.90)[1]
top90_df <- imp_df[1:k90]
top90_feats <- top90_df$Protein
ordered_feats <- top90_feats
##SFS####
X_sfs <- X[, ordered_feats, drop = FALSE]
delta_auc <- 0
patience <- 3L
feats_now <- character(0)
sfs_res <- list()
pred_list <- vector("list", length(ordered_feats))
best_auc <- -Inf
best_k <- NA_integer_
no_improve <- 0L
for (i in seq_along(ordered_feats)) {
  feats_now <- c(feats_now, ordered_feats[i])
  pred_now <- rep(NA_real_, length(y))
  auc_fold <- rep(NA_real_, nfold)
  for (k in seq_len(nfold)) {
te <- folds[[k]]
tr <- setdiff(seq_len(nrow(X_sfs)), te)
dtrain <- lightgbm::lgb.Dataset(data = X_sfs[tr, feats_now, drop = FALSE], label = y[tr])
model <- lightgbm::lgb.train(params = lgb_params, data = dtrain, nrounds = nrounds, verbose = -1)
pred_te <- as.numeric(predict(model, X_sfs[te, feats_now, drop = FALSE]))
pred_now[te] <- pred_te
roc_k <- pROC::roc(y[te], pred_te, levels = c(0, 1), direction = "<", quiet = TRUE)
auc_fold[k] <- as.numeric(pROC::auc(roc_k))
  }
roc_all <- pROC::roc(y, pred_now, levels = c(0, 1), direction = "<", quiet = TRUE)
auc_all <- as.numeric(pROC::auc(roc_all))
sfs_res[[i]] <- data.table(
    k = i,
    added = ordered_feats[i],
    n_feats = length(feats_now),
    AUC_all = auc_all,
    AUC_mean = mean(auc_fold),
    AUC_std = sd(auc_fold)
  )
pred_list[[i]] <- pred_now
if (auc_all > best_auc + delta_auc) {
    best_auc <- auc_all
    best_k <- i
    no_improve <- 0L
  } else {
    no_improve <- no_improve + 1L
  }
  
  if (no_improve >= patience) break
}
sfs_results <- rbindlist(sfs_res, fill = TRUE)
best_features <- ordered_feats[1:best_k]
best_pred <- pred_list[[best_k]]
final_panel <- ordered_feats[seq_len(best_k)]
##
##OOF####
protein_df <- data.table(ID = ml_ids, event = y)
protein_df <- cbind(protein_df, as.data.table(X))
df_model <- cbind(protein_df,
  clin_use[, .(Age, Gender, eGFR, UTP24h)]
)
model_defs <- list(
  Protein = final_panel,
  Protein_Clinical = c(final_panel, "eGFR", "Age", "Gender", "UTP24h"),
  Clinical = c("eGFR", "Age", "Gender", "UTP24h")
)
y_model <- y
folds_model <- folds
lgb_params_model <- lgb_params

pred_mat <- data.table(ID = df_model$ID, event = y_model)

for (nm in names(model_defs)) {
  feats <- model_defs[[nm]]
  X_now <- as.matrix(df_model[, ..feats])
  storage.mode(X_now) <- "numeric"
  pred_now <- rep(NA_real_, length(y_model))
  
  for (k in seq_len(nfold)) {
    te <- folds_model[[k]]
    tr <- setdiff(seq_len(nrow(X_now)), te)
    
    dtrain <- lightgbm::lgb.Dataset(data = X_now[tr, , drop = FALSE], label = y_model[tr])
    model <- lightgbm::lgb.train(params = lgb_params_model, data = dtrain, nrounds = nrounds, verbose = -1)
    pred_now[te] <- as.numeric(predict(model, X_now[te, , drop = FALSE]))
  }
  
  pred_mat[[nm]] <- pred_now
}


inner_folds <- vector("list", nfold)

for (k in seq_len(nfold)) {
  te <- folds_model[[k]]
  tr <- setdiff(seq_along(y_model), te)
  y_train <- y_model[tr]
  
  idx_pos <- sample(which(y_train == 1))
  idx_neg <- sample(which(y_train == 0))
  folds_now <- vector("list", nfold)
  
  for (j in seq_len(nfold)) {
    folds_now[[j]] <- c(
      idx_pos[seq(j, length(idx_pos), by = nfold)],
      idx_neg[seq(j, length(idx_neg), by = nfold)]
    )
  }
  
  inner_folds[[k]] <- folds_now
}


for (nm in names(model_defs)) {
  feats <- model_defs[[nm]]
  X_cal <- as.matrix(df_model[, ..feats])
  storage.mode(X_cal) <- "numeric"
  pred_cal <- rep(NA_real_, length(y_model))
  
  for (k in seq_len(nfold)) {
    te <- folds_model[[k]]
    tr <- setdiff(seq_len(nrow(X_cal)), te)
    y_train <- y_model[tr]
    inner_oof <- rep(NA_real_, length(tr))
    
    for (j in seq_len(nfold)) {
      va_local <- inner_folds[[k]][[j]]
      tr_local <- setdiff(seq_along(tr), va_local)
      
      tr_inner <- tr[tr_local]
      va_inner <- tr[va_local]
      
      dtrain <- lightgbm::lgb.Dataset(data = X_cal[tr_inner, , drop = FALSE], label = y_model[tr_inner])
      model <- lightgbm::lgb.train(params = lgb_params_model, data = dtrain, nrounds = nrounds, verbose = -1)
      inner_oof[va_local] <- as.numeric(predict(model, X_cal[va_inner, , drop = FALSE]))
    }
    
    inner_logit <- qlogis(pmin(pmax(inner_oof, 1e-6), 1 - 1e-6))
    sigmoid_model <- glm(event ~ score, data = data.frame(event = y_train, score = inner_logit), family = binomial())
    
    outer_raw <- pred_mat[[nm]][te]
    outer_logit <- qlogis(pmin(pmax(outer_raw, 1e-6), 1 - 1e-6))
pred_cal[te] <- as.numeric(
predict(sigmoid_model, newdata = data.frame(score = outer_logit), type = "response")
    )
  }
pred_mat[[paste0(nm, "_Calibrated")]] <- pred_cal
}


pred_calibrated <- pred_mat[, .(
  ID,
  event,
  Protein_Calibrated,
  Protein_Clinical_Calibrated,
  Clinical_Calibrated
)]


auc_calibrated <- rbindlist(lapply(names(model_defs), function(nm) {
  pred_name <- paste0(nm, "_Calibrated")
  roc_now <- pROC::roc(y_model, pred_mat[[pred_name]], levels = c(0, 1), direction = "<", quiet = TRUE)
  data.table(Model = nm, AUC = as.numeric(pROC::auc(roc_now)))
}))


##Calibration results####
B_cal <- 500L
calibrated_columns <- c(
  Protein = "Protein_Calibrated",
  Protein_Clinical = "Protein_Clinical_Calibrated",
  Clinical = "Clinical_Calibrated"
)
calibration_results <- list()
for (nm in names(calibrated_columns)) {
p <- pred_mat[[calibrated_columns[[nm]]]]
p_use <- pmin(pmax(p, 1e-6), 1 - 1e-6)
lp <- qlogis(p_use)
fit_intercept <- glm(y_model ~ 1, offset = lp, family = binomial())
fit_slope <- glm(y_model ~ lp, family = binomial())
cal_intercept <- as.numeric(coef(fit_intercept)[1])
cal_slope <- as.numeric(coef(fit_slope)[2])
brier <- mean((y_model - p)^2)
boot_intercept <- rep(NA_real_, B_cal)
boot_slope <- rep(NA_real_, B_cal)
boot_brier <- rep(NA_real_, B_cal)
  for (b in seq_len(B_cal)) {
    idx <- sample(seq_along(y_model), replace = TRUE)
    
    y_boot <- y_model[idx]
    p_boot <- p[idx]
    lp_boot <- qlogis(pmin(pmax(p_boot, 1e-6), 1 - 1e-6))
    
    fit_intercept_boot <- glm(y_boot ~ 1, offset = lp_boot, family = binomial())
    fit_slope_boot <- glm(y_boot ~ lp_boot, family = binomial())
    
    boot_intercept[b] <- as.numeric(coef(fit_intercept_boot)[1])
    boot_slope[b] <- as.numeric(coef(fit_slope_boot)[2])
    boot_brier[b] <- mean((y_boot - p_boot)^2)
  }
  
  calibration_results[[nm]] <- data.table(
    Model = nm,
    Calibration_intercept = cal_intercept,
    Intercept_lower95 = quantile(boot_intercept, 0.025, na.rm = TRUE),
    Intercept_upper95 = quantile(boot_intercept, 0.975, na.rm = TRUE),
    Calibration_slope = cal_slope,
    Slope_lower95 = quantile(boot_slope, 0.025, na.rm = TRUE),
    Slope_upper95 = quantile(boot_slope, 0.975, na.rm = TRUE),
    Brier_score = brier,
    Brier_lower95 = quantile(boot_brier, 0.025, na.rm = TRUE),
    Brier_upper95 = quantile(boot_brier, 0.975, na.rm = TRUE)
  )
}

calibration_results <- rbindlist(calibration_results)
calibration_results[, Calibration_intercept_95CI := sprintf(
  "%.3f [%.3f–%.3f]",
  Calibration_intercept, Intercept_lower95, Intercept_upper95
)]
calibration_results[, Calibration_slope_95CI := sprintf(
  "%.3f [%.3f–%.3f]",
  Calibration_slope, Slope_lower95, Slope_upper95
)]
calibration_results[, Brier_score_95CI := sprintf(
  "%.3f [%.3f–%.3f]",
  Brier_score, Brier_lower95, Brier_upper95
)]
calibration_results <- calibration_results[, .(
  Model,
  Calibration_intercept_95CI,
  Calibration_slope_95CI,
  Brier_score_95CI
)]

##ROC analysis;Bootstrap;DeLong####
library(data.table)
library(pROC)
y_cal <- as.integer(pred_mat$event)
roc_protein <- pROC::roc(y_cal, pred_mat$Protein_Calibrated, levels = c(0, 1), direction = "<", quiet = TRUE)
roc_protein_clinical <- pROC::roc(y_cal, pred_mat$Protein_Clinical_Calibrated, levels = c(0, 1), direction = "<", quiet = TRUE)
roc_clinical <- pROC::roc(y_cal, pred_mat$Clinical_Calibrated, levels = c(0, 1), direction = "<", quiet = TRUE)
auc_protein <- as.numeric(pROC::auc(roc_protein))
auc_protein_clinical <- as.numeric(pROC::auc(roc_protein_clinical))
auc_clinical <- as.numeric(pROC::auc(roc_clinical))
B <- 500L
idx_event <- which(y_cal == 1)
idx_nonevent <- which(y_cal == 0)
boot_auc <- data.table(
  Protein = rep(NA_real_, B),
  Protein_Clinical = rep(NA_real_, B),
  Clinical = rep(NA_real_, B)
)
for (b in seq_len(B)) {
idx <- c(
sample(idx_event, length(idx_event), replace = TRUE),
sample(idx_nonevent, length(idx_nonevent), replace = TRUE)
  )
  
y_boot <- y_cal[idx]
boot_auc$Protein[b] <- as.numeric(pROC::auc(
pROC::roc(y_boot, pred_mat$Protein_Calibrated[idx], levels = c(0, 1), direction = "<", quiet = TRUE)
  ))
boot_auc$Protein_Clinical[b] <- as.numeric(pROC::auc(
pROC::roc(y_boot, pred_mat$Protein_Clinical_Calibrated[idx], levels = c(0, 1), direction = "<", quiet = TRUE)
  ))
boot_auc$Clinical[b] <- as.numeric(pROC::auc(
pROC::roc(y_boot, pred_mat$Clinical_Calibrated[idx], levels = c(0, 1), direction = "<", quiet = TRUE)
  ))
}

auc_results <- data.table(
Model = c("Protein", "Protein_Clinical", "Clinical"),
AUC = c(auc_protein, auc_protein_clinical, auc_clinical),
  Lower95 = c(
    quantile(boot_auc$Protein, 0.025),
    quantile(boot_auc$Protein_Clinical, 0.025),
    quantile(boot_auc$Clinical, 0.025)
  ),
  Upper95 = c(
    quantile(boot_auc$Protein, 0.975),
    quantile(boot_auc$Protein_Clinical, 0.975),
    quantile(boot_auc$Clinical, 0.975)
  )
)
auc_results[, AUC_95CI := sprintf("%.3f [%.3f–%.3f]", AUC, Lower95, Upper95)]
boot_auc[, Delta_PC_vs_C := Protein_Clinical - Clinical]
boot_auc[, Delta_P_vs_C := Protein - Clinical]
boot_auc[, Delta_PC_vs_P := Protein_Clinical - Protein]
delong_pc_c <- pROC::roc.test(roc_protein_clinical, roc_clinical, method = "delong", paired = TRUE)
delong_p_c <- pROC::roc.test(roc_protein, roc_clinical, method = "delong", paired = TRUE)
delong_pc_p <- pROC::roc.test(roc_protein_clinical, roc_protein, method = "delong", paired = TRUE)
auc_comparison <- data.table(
  Comparison = c(
    "Protein + Clinical vs Clinical",
    "Protein vs Clinical",
    "Protein + Clinical vs Protein"
  ),
  Delta_AUC = c(
    auc_protein_clinical - auc_clinical,
    auc_protein - auc_clinical,
    auc_protein_clinical - auc_protein
  ),
  Lower95 = c(
    quantile(boot_auc$Delta_PC_vs_C, 0.025),
    quantile(boot_auc$Delta_P_vs_C, 0.025),
    quantile(boot_auc$Delta_PC_vs_P, 0.025)
  ),
  Upper95 = c(
    quantile(boot_auc$Delta_PC_vs_C, 0.975),
    quantile(boot_auc$Delta_P_vs_C, 0.975),
    quantile(boot_auc$Delta_PC_vs_P, 0.975)
  ),
  P_Delong = c(
    delong_pc_c$p.value,
    delong_p_c$p.value,
    delong_pc_p$p.value
  )
)
auc_comparison[, Delta_AUC_95CI := sprintf("%.3f [%.3f–%.3f]", Delta_AUC, Lower95, Upper95)]
##DCA####
y_dca <- as.integer(pred_mat$event)
N_dca <- length(y_dca)
event_rate <- mean(y_dca)
thresholds <- seq(0, 0.99, by = 0.01)
dca_columns <- c(
  "Protein" = "Protein_Calibrated",
  "Protein + Clinical" = "Protein_Clinical_Calibrated",
  "Clinical" = "Clinical_Calibrated"
)

dca_list <- list()

for (model_name in names(dca_columns)) {
  p <- pred_mat[[dca_columns[[model_name]]]]
  net_benefit <- numeric(length(thresholds))
  
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]
    positive <- p >= threshold
    TP <- sum(positive & y_dca == 1)
    FP <- sum(positive & y_dca == 0)
    net_benefit[i] <- TP / N_dca - FP / N_dca * threshold / (1 - threshold)
  }
  
  dca_list[[model_name]] <- data.table(
    Strategy = model_name,
    Threshold = thresholds,
    Net_benefit = net_benefit
  )
}

dca_models <- rbindlist(dca_list)

dca_treat_all <- data.table(
  Strategy = "Treat all",
  Threshold = thresholds,
  Net_benefit = event_rate - (1 - event_rate) * thresholds / (1 - thresholds)
)

dca_treat_none <- data.table(
  Strategy = "Treat none",
  Threshold = thresholds,
  Net_benefit = 0
)

dca_results <- rbindlist(list(dca_models, dca_treat_all, dca_treat_none))

dca_results[, Strategy := factor(
  Strategy,
  levels = c("Protein", "Protein + Clinical", "Clinical", "Treat all", "Treat none")
)]