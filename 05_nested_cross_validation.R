library(data.table)
library(survival)
library(lightgbm)
library(pROC)

outer_nfold <- 10
inner_nfold <- 5
gain_cutoff <- 0.90
delta_auc <- 0
patience <- 3
nrounds <- 500

lgb_params <- list(
  objective = "binary",
  metric = "auc",
  learning_rate = 0.01,
  num_leaves = 10,
  max_depth = 15,
  feature_fraction = 1,
  min_gain_to_split = 0.01,
  verbosity = -1,
  seed = 1
)

cv_seed <- seed
set.seed(cv_seed)

a <- as.matrix(a_MNIGA[, -1, with = FALSE])
storage.mode(a) <- "numeric"
a <- t(scale(t(a)))
pro_data <- cbind(data.table(RowName = a_MNIGA$RowName), as.data.table(a))

b <- rbindlist(list(as.data.table(amn), as.data.table(aiga)), use.names = TRUE, fill = TRUE)

c <- rbindlist(list(
  MNann[, .(ID, Age, Gender, eGFR, UTP24 = `24hUTP`)],
  IGAann[, .(ID, Age, Gender, eGFR, UTP24 = `24hUTP`)]
), use.names = TRUE, fill = TRUE)

c <- unique(c, by = "ID")
c[, Gender01 := fifelse(tolower(trimws(as.character(Gender))) == "male", 1L,
                        fifelse(tolower(trimws(as.character(Gender))) == "female", 0L, NA_integer_))]

pro <- intersect(unique(overlap_sig$Protein), pro_data$RowName)
ids <- intersect(setdiff(names(pro_data), "RowName"), b$ID)

x <- t(as.matrix(pro_data[match(pro, RowName), ..ids]))
colnames(x) <- pro
rownames(x) <- ids
storage.mode(x) <- "numeric"

surv <- b[match(ids, ID)]
clin <- c[match(ids, ID)]

keep <- complete.cases(
  surv$time_years, surv$event, surv$eGFR, surv$Age, surv$Gender,
  clin$Age, clin$Gender01, clin$eGFR, clin$UTP24, x
)

x <- x[keep, , drop = FALSE]
ids <- ids[keep]
surv <- surv[keep]
clin <- clin[keep]
y <- as.integer(surv$event)

pos <- sample(which(y == 1))
neg <- sample(which(y == 0))

folds <- lapply(seq_len(outer_nfold), function(k) {
  c(pos[seq(k, length(pos), by = outer_nfold)],
    neg[seq(k, length(neg), by = outer_nfold)])
})

pred_all <- list()
sum_all <- list()
cox_all <- list()
gain_all <- list()
panel_all <- list()
sfs_all <- list()

for (i in seq_len(outer_nfold)) {
  
  te <- folds[[i]]
  tr <- setdiff(seq_len(nrow(x)), te)
  
  xtr <- x[tr, , drop = FALSE]
  xte <- x[te, , drop = FALSE]
  ytr <- y[tr]
  yte <- y[te]
  
  surv_tr <- surv[tr]
  clin_tr <- clin[tr]
  clin_te <- clin[te]
  
  # Cox
  cox1 <- rbindlist(lapply(pro, function(p) {
    
    dat <- data.table(
      time_years = surv_tr$time_years,
      event = ytr,
      Protein = xtr[, p],
      eGFR = surv_tr$eGFR,
      Age = surv_tr$Age,
      Gender = surv_tr$Gender
    )
    
    fit <- coxph(Surv(time_years, event) ~ Protein + eGFR + Age + Gender, data = dat)
    s <- summary(fit)
    
    data.table(
      Fold = i,
      Protein = p,
      HR = s$conf.int["Protein", "exp(coef)"],
      P = s$coefficients["Protein", "Pr(>|z|)"]
    )
  }))
  
  cox1[, Selected := P < 0.05]
  cox_all[[i]] <- cox1
  sel <- cox1[Selected == TRUE, Protein]
  
  # Inner folds
  pos1 <- sample(which(ytr == 1))
  neg1 <- sample(which(ytr == 0))
  
  inner <- lapply(seq_len(inner_nfold), function(k) {
    c(pos1[seq(k, length(pos1), by = inner_nfold)],
      neg1[seq(k, length(neg1), by = inner_nfold)])
  })
  
  # Gain
  gsum <- setNames(rep(0, length(sel)), sel)
  
  for (j in seq_len(inner_nfold)) {
    
    te1 <- inner[[j]]
    tr1 <- setdiff(seq_len(nrow(xtr)), te1)
    
    par <- lgb_params
    par$scale_pos_weight <- sum(ytr[tr1] == 0) / sum(ytr[tr1] == 1)
    
    m <- lgb.train(
      par,
      lgb.Dataset(xtr[tr1, sel, drop = FALSE], label = ytr[tr1]),
      nrounds = nrounds
    )
    
    imp <- lgb.importance(m, percentage = FALSE)
    g <- setNames(rep(0, length(sel)), sel)
    g[imp$Feature] <- imp$Gain
    
    if (sum(g) > 0) g <- g / sum(g)
    gsum <- gsum + g
  }
  
  gain1 <- data.table(Fold = i, Protein = sel, Gain = as.numeric(gsum / inner_nfold))
  setorder(gain1, -Gain)
  gain1[, Cum := cumsum(Gain) / sum(Gain)]
  
  k90 <- which(gain1$Cum >= gain_cutoff)[1]
  gain1[, Top90 := seq_len(.N) <= k90]
  gain_all[[i]] <- gain1
  
  ord <- gain1[Top90 == TRUE, Protein]
  
  # SFS
  cur <- character()
  tmp <- list()
  best <- -Inf
  best_k <- 1
  stop_n <- 0
  
  for (j in seq_along(ord)) {
    
    cur <- c(cur, ord[j])
    pred <- rep(NA_real_, length(ytr))
    auc1 <- rep(NA_real_, inner_nfold)
    
    for (k in seq_len(inner_nfold)) {
      
      te1 <- inner[[k]]
      tr1 <- setdiff(seq_len(nrow(xtr)), te1)
      
      par <- lgb_params
      par$scale_pos_weight <- sum(ytr[tr1] == 0) / sum(ytr[tr1] == 1)
      
      m <- lgb.train(
        par,
        lgb.Dataset(xtr[tr1, cur, drop = FALSE], label = ytr[tr1]),
        nrounds = nrounds
      )
      
      p1 <- predict(m, xtr[te1, cur, drop = FALSE])
      pred[te1] <- p1
      auc1[k] <- as.numeric(auc(roc(ytr[te1], p1, levels = c(0, 1), direction = "<", quiet = TRUE)))
    }
    
    auc2 <- as.numeric(auc(roc(ytr, pred, levels = c(0, 1), direction = "<", quiet = TRUE)))
    
    tmp[[j]] <- data.table(
      Fold = i,
      Step = j,
      Protein = ord[j],
      OOF_AUC = auc2,
      Mean_AUC = mean(auc1),
      SD_AUC = sd(auc1)
    )
    
    if (auc2 > best + delta_auc) {
      best <- auc2
      best_k <- j
      stop_n <- 0
    } else {
      stop_n <- stop_n + 1
    }
    
    if (stop_n >= patience) break
  }
  
  sfs_all[[i]] <- rbindlist(tmp)
  
  panel <- ord[seq_len(best_k)]
  panel_all[[i]] <- data.table(Fold = i, Rank = seq_along(panel), Protein = panel)
  
  # Final models
  cv <- c("eGFR", "Age", "Gender01", "UTP24")
  
  xtr_p <- xtr[, panel, drop = FALSE]
  xte_p <- xte[, panel, drop = FALSE]
  
  xtr_c <- as.matrix(clin_tr[, ..cv])
  xte_c <- as.matrix(clin_te[, ..cv])
  storage.mode(xtr_c) <- "numeric"
  storage.mode(xte_c) <- "numeric"
  
  xtr_pc <- cbind(xtr_p, xtr_c)
  xte_pc <- cbind(xte_p, xte_c)
  
  par <- lgb_params
  par$scale_pos_weight <- sum(ytr == 0) / sum(ytr == 1)
  
  m1 <- lgb.train(par, lgb.Dataset(xtr_p, label = ytr), nrounds = nrounds)
  m2 <- lgb.train(par, lgb.Dataset(xtr_pc, label = ytr), nrounds = nrounds)
  m3 <- lgb.train(par, lgb.Dataset(xtr_c, label = ytr), nrounds = nrounds)
  
  p1_raw <- predict(m1, xte_p)
  p2_raw <- predict(m2, xte_pc)
  p3_raw <- predict(m3, xte_c)
  
  # Sigmoid calibration
  z1 <- rep(NA_real_, length(ytr))
  z2 <- rep(NA_real_, length(ytr))
  z3 <- rep(NA_real_, length(ytr))
  
  for (k in seq_len(inner_nfold)) {
    
    te1 <- inner[[k]]
    tr1 <- setdiff(seq_len(length(ytr)), te1)
    
    par <- lgb_params
    par$scale_pos_weight <- sum(ytr[tr1] == 0) / sum(ytr[tr1] == 1)
    
    a1 <- lgb.train(par, lgb.Dataset(xtr_p[tr1, , drop = FALSE], label = ytr[tr1]), nrounds = nrounds)
    a2 <- lgb.train(par, lgb.Dataset(xtr_pc[tr1, , drop = FALSE], label = ytr[tr1]), nrounds = nrounds)
    a3 <- lgb.train(par, lgb.Dataset(xtr_c[tr1, , drop = FALSE], label = ytr[tr1]), nrounds = nrounds)
    
    z1[te1] <- predict(a1, xtr_p[te1, , drop = FALSE])
    z2[te1] <- predict(a2, xtr_pc[te1, , drop = FALSE])
    z3[te1] <- predict(a3, xtr_c[te1, , drop = FALSE])
  }
  
  s1 <- qlogis(pmin(pmax(z1, 1e-6), 1 - 1e-6))
  s2 <- qlogis(pmin(pmax(z2, 1e-6), 1 - 1e-6))
  s3 <- qlogis(pmin(pmax(z3, 1e-6), 1 - 1e-6))
  
  f1 <- glm(ytr ~ s1, family = binomial())
  f2 <- glm(ytr ~ s2, family = binomial())
  f3 <- glm(ytr ~ s3, family = binomial())
  
  q1 <- qlogis(pmin(pmax(p1_raw, 1e-6), 1 - 1e-6))
  q2 <- qlogis(pmin(pmax(p2_raw, 1e-6), 1 - 1e-6))
  q3 <- qlogis(pmin(pmax(p3_raw, 1e-6), 1 - 1e-6))
  
  p1 <- plogis(coef(f1)[1] + coef(f1)[2] * q1)
  p2 <- plogis(coef(f2)[1] + coef(f2)[2] * q2)
  p3 <- plogis(coef(f3)[1] + coef(f3)[2] * q3)
  
  pred_all[[i]] <- data.table(
    ID = ids[te],
    event = yte,
    Fold = i,
    Pro = p1,
    Pro_Clin = p2,
    Clin = p3
  )
  
  sum_all[[i]] <- data.table(
    Fold = i,
    Cox_N = length(sel),
    Panel_N = length(panel),
    Panel = paste(panel, collapse = ";"),
    Inner_AUC = best,
    AUC_Pro = as.numeric(auc(roc(yte, p1, levels = c(0, 1), direction = "<", quiet = TRUE))),
    AUC_Pro_Clin = as.numeric(auc(roc(yte, p2, levels = c(0, 1), direction = "<", quiet = TRUE))),
    AUC_Clin = as.numeric(auc(roc(yte, p3, levels = c(0, 1), direction = "<", quiet = TRUE)))
  )
}

# Results
pred <- rbindlist(pred_all)
fold_res <- rbindlist(sum_all)
cox_res <- rbindlist(cox_all)
gain_res <- rbindlist(gain_all)
panels <- rbindlist(panel_all)
sfs <- rbindlist(sfs_all)

# Pooled nested OOF AUC
r1 <- roc(pred$event, pred$Pro, levels = c(0, 1), direction = "<", quiet = TRUE)
r2 <- roc(pred$event, pred$Pro_Clin, levels = c(0, 1), direction = "<", quiet = TRUE)
r3 <- roc(pred$event, pred$Clin, levels = c(0, 1), direction = "<", quiet = TRUE)

auc_oof <- data.table(
  Model = c("Pro", "Pro_Clin", "Clin"),
  AUC = c(
    as.numeric(auc(r1)),
    as.numeric(auc(r2)),
    as.numeric(auc(r3))
  )
)

# Paired DeLong
delong <- data.table(
  Comparison = c("Pro vs Clin", "Pro_Clin vs Clin", "Pro_Clin vs Pro"),
  Delta_AUC = c(
    as.numeric(auc(r1) - auc(r3)),
    as.numeric(auc(r2) - auc(r3)),
    as.numeric(auc(r2) - auc(r1))
  ),
  P = c(
    roc.test(r1, r3, method = "delong", paired = TRUE)$p.value,
    roc.test(r2, r3, method = "delong", paired = TRUE)$p.value,
    roc.test(r2, r1, method = "delong", paired = TRUE)$p.value
  )
)

# Calibration
mods <- c(Pro = "Pro", Pro_Clin = "Pro_Clin", Clin = "Clin")

cal <- rbindlist(lapply(names(mods), function(nm) {
  
  p <- pred[[mods[nm]]]
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  
  f1 <- glm(pred$event ~ 1, offset = lp, family = binomial())
  f2 <- glm(pred$event ~ lp, family = binomial())
  
  data.table(
    Model = nm,
    Intercept = coef(f1)[1],
    Slope = coef(f2)[2],
    Brier = mean((pred$event - p)^2)
  )
}))

# DCA
th <- seq(0.01, 0.50, by = 0.01)

dca <- rbindlist(lapply(names(mods), function(nm) {
  
  p <- pred[[mods[nm]]]
  
  rbindlist(lapply(th, function(t) {
    
    pos <- p >= t
    tp <- sum(pos & pred$event == 1)
    fp <- sum(pos & pred$event == 0)
    
    data.table(
      Model = nm,
      Threshold = t,
      Net_benefit = tp / nrow(pred) - fp / nrow(pred) * t / (1 - t)
    )
  }))
}))

er <- mean(pred$event)

dca <- rbindlist(list(
  dca,
  data.table(
    Model = "Treat all",
    Threshold = th,
    Net_benefit = er - (1 - er) * th / (1 - th)
  ),
  data.table(
    Model = "Treat none",
    Threshold = th,
    Net_benefit = 0
  )
))

# Stability
cox_stability <- cox_res[Selected == TRUE, .(N = uniqueN(Fold)), by = Protein]
cox_stability[, Frequency := N / outer_nfold]
setorder(cox_stability, -N)
panel_stability <- panels[, .(N = uniqueN(Fold)), by = Protein]
panel_stability[, Frequency := N / outer_nfold]
setorder(panel_stability, -N)