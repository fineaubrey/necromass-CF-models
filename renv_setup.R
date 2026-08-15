# ============================================================
# Optional: initialize renv for this repository
#
# Run ONCE from the RStudio project root after the analysis
# scripts run successfully on your machine.
#
# This creates a project-specific renv.lock containing the exact
# R package versions used locally.
# ============================================================

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Initialize renv if the project is not already using it.
if (!file.exists("renv.lock")) {
  renv::init(bare = TRUE)
}

# Discover package dependencies used by project scripts.
deps <- renv::dependencies()

print(
  unique(
    deps$Package
  )
)

# Install any packages that are referenced by the project but
# missing from the renv library.
renv::hydrate()

# Snapshot exact package versions into renv.lock.
renv::snapshot()

message(
  "\nrenv setup complete.\n",
  "Commit renv.lock and renv/activate.R to Git.\n",
  "Do not commit renv/library/."
)
