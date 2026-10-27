########################################################################
## Ki-67, 1p36 loss and 9p21/CDKN2A loss as prognostic factors in skull
## base chordoma: systematic review and meta-analysis
##
## Reproduces every pooled estimate, subgroup analysis, sensitivity
## analysis and small-study-effect test reported in the article and its
## supplementary material.
##
## Input : data/estimates.csv, data/cohort_attributes.csv
## Output: results/syntheses.csv, results/subgroups.csv,
##         results/sensitivity.csv, results/small_study_effects.csv
##
## Usage : source("R/meta_analysis.R"); res <- run_meta_analysis()
##         or, from the command line: Rscript R/meta_analysis.R
## Requires: metafor (>= 4.0)
########################################################################

suppressPackageStartupMessages(library(metafor))

## Syntheses defined for the review. "pool" marks syntheses eligible for a
## pooled estimate; the others are displayed without a summary estimate.
SYNTHESES <- read.table(header = TRUE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE, quote = "", text = "
id|biomarker|outcome|analysis|pool|scale
KI67-BIN-PFS-ADJ|Ki-67, binary|Recurrence/progression|Adjusted HRs (primary analysis)|TRUE|HR
KI67-BIN-OS-ADJ|Ki-67, binary|Overall survival|Adjusted HRs|TRUE|HR
KI67-BIN-PFS-UNADJ|Ki-67, binary|Recurrence/progression|Unadjusted and derived HRs|TRUE|HR
KI67-BIN-OS-UNADJ|Ki-67, binary|Overall survival|Unadjusted and derived HRs|TRUE|HR
KI67-CONT-PFS-ADJ|Ki-67, continuous (per 1 percentage point)|Recurrence/progression|Adjusted HRs|TRUE|HR per 1 percentage point
KI67-CONT-OS-ADJ|Ki-67, continuous (per 1 percentage point)|Overall survival|Adjusted HRs|TRUE|HR per 1 percentage point
KI67-CONT-PFS-UNADJ|Ki-67, continuous (per 1 percentage point)|Recurrence/progression|Unadjusted HRs|TRUE|HR per 1 percentage point
KI67-CONT-OS-UNADJ|Ki-67, continuous (per 1 percentage point)|Overall survival|Unadjusted HRs|TRUE|HR per 1 percentage point
KI67-CONT-DSS-DISPLAY|Ki-67, continuous (per 1 percentage point)|Disease-specific survival|Display|FALSE|HR per 1 percentage point
KI67-BIN-DSS-DISPLAY|Ki-67, binary|Disease-specific survival|Display|FALSE|HR
KI67-HSCORE-PFS-ADJ|Ki-67, H-score (per unit, 0-300)|Recurrence/progression|Display|FALSE|HR per unit
KI67-ORD-PFS-ADJ|Ki-67, ordinal categories|Recurrence/progression|Display|FALSE|HR
ZENONOS-PFS-AFTER-RT-DISPLAY|Ki-67, 1p36 and 9p21 (Zenonos 2018)|Progression after radiotherapy|Display|FALSE|HR
1P36-PFS-DISPLAY|1p36 loss|Recurrence/progression|Display|FALSE|HR
1P36-OS-DISPLAY|1p36 loss|Overall survival|Display|FALSE|HR
9P21-PFS-DISPLAY|9p21/CDKN2A loss|Recurrence/progression|Display|FALSE|HR
9P21-OS-DSS-DISPLAY|9p21/CDKN2A loss|Overall or disease-specific survival|Display|FALSE|HR
")

MIN_K <- 3         # minimum number of independent cohorts for a pooled estimate
MIN_K_EGGER <- 10  # minimum number of cohorts for funnel plot and Egger's test

## ------------------------------------------------------------ helpers
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
hr_ci <- function(b, lb, ub) paste0(f2(exp(b)), " (", f2(exp(lb)), "–", f2(exp(ub)), ")")

## Random-effects model: REML estimate of tau^2 and Hartung-Knapp-Sidik-Jonkman
## confidence interval with the ad hoc variance correction (metafor: test = "adhoc").
fit_re <- function(dd, method = "REML") {
  if (nrow(dd) < 2) return(NULL)
  tryCatch(rma(yi = dd$yi, sei = dd$sei, method = method, test = "adhoc",
               control = list(stepadj = 0.5, maxiter = 1000)), error = function(z) NULL)
}

read_inputs <- function(data_dir = "data") {
  e <- read.csv(file.path(data_dir, "estimates.csv"), stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  at <- read.csv(file.path(data_dir, "cohort_attributes.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
  e$yi <- log(e$hr); e$sei <- e$se_log_hr
  e$ci_lower[!is.na(e$ci_lower) & e$ci_lower <= 0] <- NA
  e$ci_upper[!is.na(e$ci_upper) & e$ci_upper <= 0] <- NA
  e <- merge(e, at[, setdiff(names(at), "study")], by = "study_id", all.x = TRUE, sort = FALSE)
  e$direct <- e$source_code %in% c("reported", "ci_from_p")
  e
}

## Rows of a synthesis used in the main analysis or display
synthesis_rows <- function(e, sid) {
  dd <- e[e$synthesis_id == sid & e$role %in% c("main", "display") & !is.na(e$yi) & !is.na(e$sei), ]
  dd[order(dd$year, dd$study), ]
}

sens_row <- function(sid, label, dd) {
  k <- nrow(dd)
  if (k < MIN_K) return(data.frame(synthesis_id = sid, analysis = label, k = k, hr_95ci = paste0("Not pooled (k = ", k, ")"),
                                   i2 = NA, studies = paste(dd$study, collapse = "; ")))
  m <- fit_re(dd)
  data.frame(synthesis_id = sid, analysis = label, k = k, hr_95ci = hr_ci(m$b, m$ci.lb, m$ci.ub),
             i2 = round(m$I2), studies = paste(dd$study, collapse = "; "))
}

## ------------------------------------------------------------ main function
run_meta_analysis <- function(data_dir = "data", results_dir = "results") {
  e <- read_inputs(data_dir)
  syn <- list(); sub <- list(); sens <- list(); sse <- list(); models <- list()

  for (i in seq_len(nrow(SYNTHESES))) {
    s <- SYNTHESES[i, ]; sid <- s$id
    dd <- synthesis_rows(e, sid)
    k <- nrow(dd); kc <- length(unique(dd$study_id))
    base <- data.frame(synthesis_id = sid, biomarker = s$biomarker, outcome = s$outcome, analysis = s$analysis,
                       k = kc, patients = sum(dd$n, na.rm = TRUE), events = sum(dd$events, na.rm = TRUE),
                       cohorts_events_nr = sum(is.na(dd$events)), studies = paste(dd$study, collapse = "; "))
    pooled <- s$pool && kc >= MIN_K && kc == k
    if (!pooled) {
      syn[[sid]] <- cbind(base, hr_95ci = "Not pooled", hr = NA, ci_lower = NA, ci_upper = NA, p = NA, i2 = NA, tau2 = NA,
                          q_p = NA, pi_lower = NA, pi_upper = NA, hr_per5_95ci = NA,
                          reason = if (s$pool) paste0("k = ", kc, " (fewer than ", MIN_K, " cohorts)") else "Display only (incompatible scales, platforms or outcomes)")
      next
    }
    m <- fit_re(dd); pr <- predict(m); models[[sid]] <- m
    syn[[sid]] <- cbind(base, hr_95ci = hr_ci(m$b, m$ci.lb, m$ci.ub), hr = exp(as.numeric(m$b)), ci_lower = exp(m$ci.lb),
                        ci_upper = exp(m$ci.ub), p = m$pval, i2 = m$I2, tau2 = m$tau2, q_p = m$QEp,
                        pi_lower = exp(pr$pi.lb), pi_upper = exp(pr$pi.ub),
                        hr_per5_95ci = if (grepl("CONT", sid)) hr_ci(5 * m$b, 5 * m$ci.lb, 5 * m$ci.ub) else NA, reason = "")

    ## small-study effects: contour-enhanced funnel plot and Egger's test (SE as predictor)
    if (kc >= MIN_K_EGGER) {
      eg <- regtest(m, model = "rma", predictor = "sei")
      sse[[sid]] <- data.frame(synthesis_id = sid, k = kc, egger_z = eg$zval, egger_p = eg$pval,
                               limit_estimate_95ci = hr_ci(eg$est, eg$ci.lb, eg$ci.ub))
    }

    ## sensitivity analyses
    S <- list()
    for (meth in c("DL", "PM")) {
      mm <- fit_re(dd, meth)
      S[[length(S) + 1]] <- data.frame(synthesis_id = sid,
        analysis = if (meth == "DL") "DerSimonian–Laird estimator of τ²" else "Paule–Mandel estimator of τ²",
        k = k, hr_95ci = hr_ci(mm$b, mm$ci.lb, mm$ci.ub), i2 = round(mm$I2), studies = "")
    }
    l1 <- leave1out(m)
    for (j in seq_len(k)) S[[length(S) + 1]] <- data.frame(synthesis_id = sid, analysis = paste("Leave-one-out: without", dd$study[j]),
      k = k - 1, hr_95ci = hr_ci(l1$estimate[j], l1$ci.lb[j], l1$ci.ub[j]), i2 = round(l1$I2[j]), studies = "")
    small <- (!is.na(dd$n) & dd$n < 20) | (!is.na(dd$events) & dd$events < 10)
    S[[length(S) + 1]] <- sens_row(sid, "Excluding cohorts with <20 participants or <10 events", dd[!small, ])
    S[[length(S) + 1]] <- sens_row(sid, "Directly reported estimates only", dd[dd$direct, ])
    if (grepl("BIN", sid)) S[[length(S) + 1]] <- sens_row(sid, "Excluding data-derived cut-offs (mean, median or ROC)", dd[dd$data_derived_cutoff != "yes", ])
    S[[length(S) + 1]] <- sens_row(sid, "Index-operation specimens only", dd[dd$specimen == "index_primary", ])
    S[[length(S) + 1]] <- sens_row(sid, "Excluding paediatric-only cohorts", dd[dd$paediatric_only != "yes", ])
    if (any(dd$tiantan_cluster == "yes")) S[[length(S) + 1]] <- sens_row(sid, "Excluding the Beijing Tiantan cohort cluster", dd[dd$tiantan_cluster != "yes", ])
    if (any(!is.na(dd$quips_overall))) S[[length(S) + 1]] <- sens_row(sid, "Excluding estimates at high risk of bias", dd[!is.na(dd$quips_overall) & dd$quips_overall != "high", ])
    ## alternative estimates for the same cohort, substituted one at a time
    alt <- e[e$synthesis_id == sid & e$role == "sensitivity" & !is.na(e$yi), ]
    for (j in seq_len(nrow(alt))) {
      dd2 <- dd; w <- which(dd2$study_id == alt$study_id[j])
      if (length(w) == 1) { dd2[w, c("yi", "sei")] <- alt[j, c("yi", "sei")]
        S[[length(S) + 1]] <- sens_row(sid, alt$sensitivity_label[j], dd2) }
    }
    sens[[sid]] <- do.call(rbind, S)

    ## subgroup analysis by Ki-67 cut-off (around 5% [4-6%] versus other values), binary scale;
    ## done when each level has at least two cohorts; levels with fewer than three cohorts
    ## are described without inference
    if (grepl("BIN", sid)) {
      dd$level <- ifelse(!is.na(dd$ki67_cutoff_percent) & dd$ki67_cutoff_percent >= 4 & dd$ki67_cutoff_percent <= 6,
                         "Cut-off around 5% (4–6%)", "Other cut-offs")
      if (all(table(dd$level) >= 2) && length(unique(dd$level)) == 2) {
        mr <- rma(yi = dd$yi, sei = dd$sei, mods = ~ factor(dd$level), method = "REML", test = "adhoc")
        for (lv in sort(unique(dd$level))) {
          sg <- dd[dd$level == lv, ]; ms <- fit_re(sg)
          sub[[paste(sid, lv)]] <- data.frame(synthesis_id = sid, level = lv, k = nrow(sg), hr_95ci = hr_ci(ms$b, ms$ci.lb, ms$ci.ub),
                                             i2 = round(ms$I2), p_subgroup_difference = mr$QMp,
                                             note = ifelse(nrow(sg) < 3, "Fewer than three cohorts: described without inference", ""),
                                             studies = paste(sg$study, collapse = "; "))
        }
      }
    }
  }

  out <- list(estimates = e, models = models,
              syntheses = do.call(rbind, syn), subgroups = do.call(rbind, sub),
              sensitivity = do.call(rbind, sens), small_study_effects = do.call(rbind, sse))
  if (!is.null(results_dir)) {
    dir.create(results_dir, showWarnings = FALSE)
    for (nm in c("syntheses", "subgroups", "sensitivity", "small_study_effects"))
      write.csv(out[[nm]], file.path(results_dir, paste0(nm, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
  }
  out
}

if (sys.nframe() == 0) {
  res <- run_meta_analysis()
  print(res$syntheses[, c("synthesis_id", "k", "hr_95ci", "i2")], row.names = FALSE)
}
