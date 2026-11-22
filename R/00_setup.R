# Paths to the raw PhysioNet files (csv.gz) and output folders.
# Set MIMIC_DIR to the mimic-iv-3.1 folder (containing hosp/ and icu/)
# and EICU_DIR to the eicu-collaborative-research-database-2.0 folder.
MIMIC_DIR <- Sys.getenv("MIMIC_DIR", "data/mimic-iv-3.1")
EICU_DIR  <- Sys.getenv("EICU_DIR",  "data/eicu-crd-2.0")
WORK_DIR  <- "work"
OUT_DIR   <- "results"
dir.create(WORK_DIR, showWarnings = FALSE)
dir.create(OUT_DIR, showWarnings = FALSE)

pkgs <- c("DBI", "duckdb", "data.table", "survival", "cmprsk", "splines", "mice",
          "EValue", "sandwich", "WeightIt", "cobalt", "ggplot2", "ragg", "scales", "jsonlite")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install packages: ", paste(miss, collapse = ", "))
suppressMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
