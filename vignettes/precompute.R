# Pre-compute the methodology vignette: knit the runnable .Rmd.orig into the
# static, shipped methodology.Rmd (executable chunks become inert markdown, so
# R CMD build re-runs nothing), then render to HTML to confirm it builds.
# Run from the package root:  Rscript vignettes/precompute.R
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(devtools::load_all("."))

old <- setwd("vignettes")
on.exit(setwd(old), add = TRUE)

knitr::knit("methodology.Rmd.orig", output = "methodology.Rmd")
rmarkdown::render("methodology.Rmd", quiet = TRUE)  # writes methodology.html

writeLines("PRECOMPUTE_DONE", "precompute.marker")
