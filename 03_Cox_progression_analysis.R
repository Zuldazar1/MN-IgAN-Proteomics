library(data.table)
library(survival)

# Standardize each protein across the combined MN and IgAN cohort
proteomics <- fread("MN_IGA.txt")
protein_matrix <- as.matrix(proteomics[, -1, with = FALSE])
storage.mode(protein_matrix) <- "numeric"
protein_matrix <- t(scale(t(protein_matrix)))

proteomics_scaled <- cbind(
  data.table(RowName = proteomics$RowName),
  as.data.table(protein_matrix)
)

# Cox regression adjusted for baseline eGFR, age and sex
cox_results <- vector("list", length(candidate_proteins))

for (i in seq_along(candidate_proteins)) {
  protein <- candidate_proteins[i]
  expr <- as.numeric(unlist(proteomics_scaled[RowName == protein, ..sample_ids],
                            use.names = FALSE))
  
  dat <- copy(cox_clinical[match(sample_ids, ID)])
  dat[, Protein := expr]
  dat <- dat[complete.cases(time_years, event, Protein, eGFR, Age, Gender)]
  
  fit <- coxph(Surv(time_years, event) ~ Protein + eGFR + Age + Gender, data = dat)
  s <- summary(fit)
  co <- s$coefficients["Protein", ]
  ci <- s$conf.int["Protein", ]
  
  cox_results[[i]] <- data.table(
    Protein = protein, HR = ci["exp(coef)"],
    Lower95 = ci["lower .95"], Upper95 = ci["upper .95"],
    P = co["Pr(>|z|)"]
  )
}

cox_results <- rbindlist(cox_results)
cox_results[, FDR := p.adjust(P, method = "BH")]
cox_significant <- cox_results[P < 0.05]