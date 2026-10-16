# The line printed when the package is attached. It goes through packageStartupMessage(), so that
# suppressPackageStartupMessages(library(rotationGWAS)) silences it.

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(sprintf("rotationGWAS %s. The means are rotated, the errors are independent, and the intercept is not a column of ones. Start with vignette('rotationGWAS').",
                                utils::packageVersion("rotationGWAS")))
}
