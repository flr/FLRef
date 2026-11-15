# Render the FLRef guideline documents to docs/ (GitHub Pages) -----------------
# Sources live in inst/rmd/, rendered HTML goes to docs/ (served by GitHub
# Pages from the main branch, /docs folder). Run from the package folder:
#   source("inst/scripts/render_docs.R")
#
# The slow MSE runs are cached as .rds in rdata/ (not pushed; listed in
# .gitignore and .Rbuildignore). Delete a cached file to re-run that MSE.

dir.create("rdata", showWarnings = FALSE)
options(FLRef.rdata   = normalizePath("rdata", mustWork = TRUE),
        FLRef.run_mse = TRUE)

docs <- c(
  lfd_single  = "inst/rmd/lfd-01-single-sex.Rmd"
  # lfd_sex     = "inst/rmd/lfd-02-sex-structured.Rmd",   # unit.type = "gender"
  # lfd_spatial = "inst/rmd/lfd-03-spatial.Rmd"           # unit.type = "spatial"
)

for (d in docs) {
  rmarkdown::render(d, output_dir = "docs", envir = new.env(),
                    output_options = list(self_contained = TRUE))
}

# GitHub Pages: serve files as-is (no Jekyll processing)
if (!file.exists("docs/.nojekyll")) file.create("docs/.nojekyll")
