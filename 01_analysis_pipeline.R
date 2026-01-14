# ============================================================
# PROJECT: EMBODY-LATAM (MOSAIC)
# SCRIPT: 01_analysis_pipeline.R
# DESCRIPTION: 
#   Performs the Main Two-Step Mendelian Randomization Analysis:
#   1. Step 1: Social Adversity (Education) -> Epigenetic Aging (Clocks)
#   2. Step 2: Epigenetic Aging (Clocks) -> Oral Cancer
#
# 
#
# AUTHOR: IRMA VALDEBENITO / EMBODY-LATAM Team
# DATE: January 2026
# ============================================================

rm(list = ls())

# --------------------------
# PACKAGES
# --------------------------
pkgs <- c("TwoSampleMR","ieugwasr","dplyr","readr","ggplot2","beepr","openxlsx","tibble","stringr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

# Optional MR-PRESSO
has_presso <- requireNamespace("MRPRESSO", quietly = TRUE)
if (!has_presso) message("MRPRESSO no instalado (opcional). Si lo quieres: install.packages('MRPRESSO')")

# --------------------------
# FOLDERS (SAVE HERE)
# --------------------------
base_dir <- "results"
dir.create(base_dir, showWarnings = FALSE)

dir.create(file.path(base_dir, "00_mapping"), showWarnings = FALSE)
dir.create(file.path(base_dir, "01_step1_Edu_to_Clock"), showWarnings = FALSE)
dir.create(file.path(base_dir, "02_step2_Clock_to_Cancer"), showWarnings = FALSE)
dir.create(file.path(base_dir, "03_summary"), showWarnings = FALSE)
dir.create(file.path(base_dir, "logs"), showWarnings = FALSE)

log_file <- file.path(base_dir, "logs", paste0("run_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
sink(log_file, split = TRUE)

cat("Run started:", as.character(Sys.time()), "\n")
cat("TwoSampleMR:", as.character(packageVersion("TwoSampleMR")), "\n")
cat("ieugwasr:", as.character(packageVersion("ieugwasr")), "\n\n")

# --------------------------
# IDS
# --------------------------
X_id <- "ebi-a-GCST90029013"  # Education (years)
Y_id <- "ieu-b-4961"          # Oral cavity cancer

# --------------------------
# NOTE ABOUT TOKENS
# --------------------------
# If you use a token:
# Sys.setenv(OPENGWAS_JWT="PASTE_YOUR_TOKEN_HERE")
# STEP 1 — AUTH (OpenGWAS JWT)
# --------------------------
# WHY: OpenGWAS requires JWT since May 2024.
# HOW: read OPENGWAS_JWT environment variable; if missing, ask once.

cat("STEP 1: OpenGWAS authentication...\n")

jwt <- trimws(Sys.getenv("OpenGWAS_TOKEN"))

if (jwt == "") {
  cat("\nNo OPENGWAS_JWT found.\nPaste your OpenGWAS JWT token and press ENTER:\n")
  jwt <- trimws(readline(prompt = "> "))
  if (jwt == "") stop("No token provided. Stopping.")
  Sys.setenv(OPENGWAS_JWT = jwt)
  cat("Token loaded into session (not saved to disk).\n")
} else {
  cat("Token found in environment.\n")
}

# Sanity check: JWT has 3 dot-separated parts
if (length(strsplit(jwt, "\\.")[[1]]) != 3) stop("Token does not look like a JWT (expected 3 dot-separated parts).")

# Live API check
cat("Testing API with ieugwasr::user()...\n")
u <- try(ieugwasr::user(), silent = TRUE)
if (inherits(u, "try-error")) stop("Token present but API call failed (401/403). Generate a NEW token and retry.")

cat("✅ Auth confirmed.\n\n")
beepr::beep(2)

# --------------------------
# HELPERS
# --------------------------
safe_write_csv <- function(df, path){
  tryCatch(readr::write_csv(df, path), error = function(e) message("write_csv error: ", e$message))
}

safe_save_plot <- function(p, path, w=8, h=6){
  tryCatch(ggsave(filename = path, plot = p, width = w, height = h, dpi = 300),
           error = function(e) message("ggsave error: ", e$message))
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

# --------------------------
# 0) MAP: FIND EPIGENETIC CLOCK GWAS
# --------------------------
ao <- TwoSampleMR::available_outcomes()

kw <- c(
  "grimage", "dunedinpace", "phenoage", "hannum", "horvath",
  "epigenetic", "dnam", "methylation age", "age acceleration", "ageaccel",
  "pace", "clock"
)
pattern <- paste0("(", paste(kw, collapse="|"), ")")

cand_clocks <- ao %>%
  filter(grepl(pattern, trait, ignore.case = TRUE)) %>%
  select(id, trait, population, sample_size, year, mr) %>%
  arrange(desc(sample_size))

safe_write_csv(cand_clocks, file.path(base_dir, "00_mapping", "candidate_epigenetic_clocks_from_OpenGWAS.csv"))

cat("Candidate clock GWAS found:", nrow(cand_clocks), "\n")
cat("Top 25 by sample size:\n")
print(head(cand_clocks, 25))

if (nrow(cand_clocks) == 0) {
  cat("\nNO clock-like GWAS were found with current keywords.\n")
  sink()
  stop("No candidate clock GWAS found.")
}

# --------------------------
# 0b) RANK BY INSTRUMENT COUNT (REAL SNPs @ p<5e-8)
# --------------------------
rank_by_instruments <- function(M_id){
  out <- tryCatch({
    inst <- TwoSampleMR::extract_instruments(outcomes = M_id, p1 = 5e-8, clump = TRUE)
    tibble(id = M_id, n_snps = nrow(inst))
  }, error = function(e){
    tibble(id = M_id, n_snps = NA_integer_)
  })
  out
}

topN <- min(40, nrow(cand_clocks))
cand_top <- cand_clocks %>% slice(1:topN)

cat("\nRanking top", topN, "candidates by number of instruments...\n")
rank_tbl <- bind_rows(lapply(cand_top$id, rank_by_instruments)) %>%
  left_join(cand_top, by = "id") %>%
  arrange(desc(n_snps), desc(sample_size))

safe_write_csv(rank_tbl, file.path(base_dir, "00_mapping", "ranked_clocks_by_instruments.csv"))

cat("\nTop 20 by instruments:\n")
print(head(rank_tbl, 20))

min_snps <- 10
selected <- rank_tbl %>%
  filter(!is.na(n_snps), n_snps >= min_snps) %>%
  arrange(desc(n_snps), desc(sample_size))

safe_write_csv(selected, file.path(base_dir, "00_mapping", paste0("selected_clocks_nsnps_ge_", min_snps, ".csv")))

cat("\nSelected clocks with n_snps >=", min_snps, ":", nrow(selected), "\n")
if (nrow(selected) == 0) {
  cat("None passed instrument threshold. Lower min_snps if needed.\n")
  sink()
  stop("No clocks passed instrument threshold.")
}

# --------------------------
# MR RUNNER
# --------------------------
run_mr_block <- function(exposure_id, outcome_id, label, outdir){
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cat("\n=============================\n")
  cat("RUN:", label, "\nExposure:", exposure_id, "\nOutcome:", outcome_id, "\n")
  cat("=============================\n")
  
  exp_dat <- tryCatch({
    TwoSampleMR::extract_instruments(outcomes = exposure_id, p1 = 5e-8, clump = TRUE)
  }, error = function(e){
    cat("extract_instruments error:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(exp_dat) || nrow(exp_dat) < 3) {
    cat("Not enough instruments (<3). Skipping.\n")
    return(list(ok=FALSE, reason="too_few_instruments", res=NULL))
  }
  
  out_dat <- tryCatch({
    TwoSampleMR::extract_outcome_data(snps = exp_dat$SNP, outcomes = outcome_id)
  }, error = function(e){
    cat("extract_outcome_data error:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(out_dat) || nrow(out_dat) == 0) {
    cat("No outcome data returned. Skipping.\n")
    return(list(ok=FALSE, reason="no_outcome_data", res=NULL))
  }
  
  dat <- tryCatch({
    TwoSampleMR::harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat)
  }, error = function(e){
    cat("harmonise_data error:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(dat) || nrow(dat) < 3) {
    cat("Not enough harmonised SNPs (<3). Skipping.\n")
    return(list(ok=FALSE, reason="too_few_harmonised", res=NULL))
  }
  
  mr_res <- TwoSampleMR::mr(dat)
  het    <- tryCatch(TwoSampleMR::mr_heterogeneity(dat), error=function(e) NULL)
  pleio  <- tryCatch(TwoSampleMR::mr_pleiotropy_test(dat), error=function(e) NULL)
  single <- tryCatch(TwoSampleMR::mr_singlesnp(dat), error=function(e) NULL)
  loo    <- tryCatch(TwoSampleMR::mr_leaveoneout(dat), error=function(e) NULL)
  
  safe_write_csv(dat,     file.path(outdir, paste0(label, "_harmonised_data.csv")))
  safe_write_csv(mr_res,  file.path(outdir, paste0(label, "_mr_results.csv")))
  if (!is.null(het))   safe_write_csv(het,   file.path(outdir, paste0(label, "_heterogeneity.csv")))
  if (!is.null(pleio)) safe_write_csv(pleio, file.path(outdir, paste0(label, "_pleiotropy_egger_intercept.csv")))
  if (!is.null(single))safe_write_csv(single,file.path(outdir, paste0(label, "_singleSNP.csv")))
  if (!is.null(loo))   safe_write_csv(loo,   file.path(outdir, paste0(label, "_leaveoneout.csv")))
  
  p_scatter <- tryCatch({
    TwoSampleMR::mr_scatter_plot(mr_res, dat)[[1]]
  }, error = function(e) NULL)
  if (!is.null(p_scatter)) safe_save_plot(p_scatter, file.path(outdir, paste0(label, "_scatter.png")))
  
  p_forest <- tryCatch({
    if (!is.null(single) && nrow(single) > 0) TwoSampleMR::mr_forest_plot(single)[[1]] else NULL
  }, error = function(e) NULL)
  if (!is.null(p_forest)) safe_save_plot(p_forest, file.path(outdir, paste0(label, "_forest.png")), w=9, h=7)
  
  p_funnel <- tryCatch({
    if (!is.null(single) && nrow(single) > 0) TwoSampleMR::mr_funnel_plot(single)[[1]] else NULL
  }, error = function(e) NULL)
  if (!is.null(p_funnel)) safe_save_plot(p_funnel, file.path(outdir, paste0(label, "_funnel.png")), w=8, h=7)
  
  p_loo <- tryCatch({
    if (!is.null(loo) && nrow(loo) > 0) TwoSampleMR::mr_leaveoneout_plot(loo)[[1]] else NULL
  }, error = function(e) NULL)
  if (!is.null(p_loo)) safe_save_plot(p_loo, file.path(outdir, paste0(label, "_leaveoneout.png")), w=9, h=7)
  
  return(list(ok=TRUE, reason="ok", mr=mr_res, het=het, pleio=pleio, n_snps=nrow(dat)))
}

# --------------------------
# MAIN LOOP: FOR EACH CLOCK
# --------------------------
all_summary <- list()

for (i in seq_len(nrow(selected))) {
  M_id <- selected$id[i]
  M_trait <- selected$trait[i] %||% "NA"
  
  clock_tag <- paste0("M", sprintf("%02d", i), "_", gsub("[^A-Za-z0-9]+", "_", substr(M_trait, 1, 60)))
  
  cat("\n\n#############################################\n")
  cat("CLOCK:", i, "/", nrow(selected), "\nID:", M_id, "\nTrait:", M_trait, "\n")
  cat("#############################################\n")
  
  # Step 1: Education -> Clock
  outdir1 <- file.path(base_dir, "01_step1_Edu_to_Clock", clock_tag)
  step1 <- run_mr_block(exposure_id = X_id, outcome_id = M_id,
                        label = paste0("STEP1_Edu_to_Clock_", clock_tag),
                        outdir = outdir1)
  
  # Step 2: Clock -> Cancer
  outdir2 <- file.path(base_dir, "02_step2_Clock_to_Cancer", clock_tag)
  step2 <- run_mr_block(exposure_id = M_id, outcome_id = Y_id,
                        label = paste0("STEP2_Clock_to_Cancer_", clock_tag),
                        outdir = outdir2)
  
  sum_row <- tibble(
    clock_id = M_id,
    clock_trait = M_trait,
    step1_ok = step1$ok,
    step1_nsnps = ifelse(step1$ok, step1$n_snps, NA_integer_),
    step1_ivw_p = if (step1$ok) step1$mr %>% filter(method=="Inverse variance weighted") %>% pull(pval) %>% dplyr::first() else NA_real_,
    step1_ivw_b = if (step1$ok) step1$mr %>% filter(method=="Inverse variance weighted") %>% pull(b) %>% dplyr::first() else NA_real_,
    step2_ok = step2$ok,
    step2_nsnps = ifelse(step2$ok, step2$n_snps, NA_integer_),
    step2_ivw_p = if (step2$ok) step2$mr %>% filter(method=="Inverse variance weighted") %>% pull(pval) %>% dplyr::first() else NA_real_,
    step2_ivw_b = if (step2$ok) step2$mr %>% filter(method=="Inverse variance weighted") %>% pull(b) %>% dplyr::first() else NA_real_
  )
  
  all_summary[[i]] <- sum_row
  try(beepr::beep(2), silent = TRUE)
}

summary_df <- bind_rows(all_summary) %>%
  arrange(step2_ivw_p, step1_ivw_p)

safe_write_csv(summary_df, file.path(base_dir, "03_summary", "two_step_summary_all_clocks.csv"))

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Selected_Clocks")
openxlsx::writeData(wb, "Selected_Clocks", selected)

openxlsx::addWorksheet(wb, "TwoStep_Summary")
openxlsx::writeData(wb, "TwoStep_Summary", summary_df)

openxlsx::saveWorkbook(wb, file.path(base_dir, "03_summary", "two_step_summary_all_clocks.xlsx"), overwrite = TRUE)

cat("\n\nDONE ✅\n")
cat("Saved in:\n", base_dir, "\n")
cat("Run finished:", as.character(Sys.time()), "\n")

sink()

# End

