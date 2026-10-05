library(data.table)
library(MASS)

MN <- fread(
  file.path(proteomics_dir, "MN.txt")
)

IGA <- fread(
  file.path(proteomics_dir, "IGA.txt")
)

mn_clinical <- fread(
  file.path(clinical_dir, "MNann.csv")
)

igan_clinical <- fread(
  file.path(clinical_dir, "IGAann.csv")
)


# CKD stage
mn_clinical[, CKD_stage := fcase(
  eGFR >= 90, "G1",
  eGFR >= 60, "G2",
  eGFR >= 45, "G3a",
  eGFR >= 30, "G3b",
  eGFR >= 15, "G4",
  eGFR < 15,  "G5"
)]

igan_clinical[, CKD_stage := fcase(
  eGFR >= 90, "G1",
  eGFR >= 60, "G2",
  eGFR >= 45, "G3a",
  eGFR >= 30, "G3b",
  eGFR >= 15, "G4",
  eGFR < 15,  "G5"
)]

CKD_levels <- c("G1", "G2", "G3a", "G3b", "G4", "G5")

mn_clinical[, CKD_stage := ordered(CKD_stage, levels = CKD_levels)]
igan_clinical[, CKD_stage := ordered(CKD_stage, levels = CKD_levels)]

mn_clinical[, `:=`(
  Gender = factor(Gender),
  Age = as.numeric(Age),
  Alb = as.numeric(Alb)
)]

igan_clinical[, `:=`(
  Gender = factor(Gender),
  Age = as.numeric(Age),
  Alb = as.numeric(Alb)
)]


# Ordinal logistic regression
run_ordinal_panel <- function(expr_dt, clinical_dt) {
  
  sample_ids <- intersect(
    setdiff(names(expr_dt), "RowName"),
    clinical_dt$ID
  )
  
  meta <- clinical_dt[
    ID %in% sample_ids,
    .(sample = ID, CKD_stage, Age, Gender, Alb)
  ]
  
  res <- vector("list", nrow(expr_dt))
  
  for (i in seq_len(nrow(expr_dt))) {
    
    dat <- data.table(
      sample = sample_ids,
      protein = as.numeric(unlist(expr_dt[i, ..sample_ids]))
    )
    
    dat <- merge(dat, meta, by = "sample")
    dat[, protein := as.numeric(scale(protein))]
    
    fit <- MASS::polr(
      CKD_stage ~ protein + Age + Gender + Alb,
      data = dat,
      Hess = TRUE,
      method = "logistic"
    )
    
    co <- coef(summary(fit))
    
    beta <- co["protein", "Value"]
    se <- co["protein", "Std. Error"]
    p <- 2 * pnorm(abs(beta / se), lower.tail = FALSE)
    
    res[[i]] <- data.table(
      Protein = expr_dt$RowName[i],
      N = nrow(dat),
      Beta = beta,
      SE = se,
      OR = exp(beta),
      LCL = exp(beta - 1.96 * se),
      UCL = exp(beta + 1.96 * se),
      P = p
    )
  }
  
  res <- rbindlist(res)
  res[, FDR := p.adjust(P, method = "BH")]
  setorder(res, FDR)
  
  res
}


# MN
res_mn <- run_ordinal_panel(
  MN,
  mn_clinical
)

# IgAN
res_igan <- run_ordinal_panel(
  IGA,
  igan_clinical
)

# Significant proteins
mn_sig <- res_mn[FDR < 0.05]
igan_sig <- res_igan[FDR < 0.05]

# Shared significant proteins
shared_sig <- merge(
  mn_sig,
  igan_sig,
  by = "Protein",
  suffixes = c("_MN", "_IgAN")
)