# Run the full pipeline from the raw PhysioNet files.
for (f in c("R/01_extract.R", "R/02_cohort.R", "R/03_analysis.R", "R/04_tables_figures.R")) {
  message("Running ", f)
  source(f, local = new.env())
}
