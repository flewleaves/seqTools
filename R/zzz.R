.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Thank you for using seqTools")
  packageStartupMessage("Version ", utils::packageVersion(pkgname))
}