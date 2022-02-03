.onLoad <- function(libname = find.package("scexpr"), pkgname = "scexpr") {

  # https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when

  # CRAN Note avoidance
  if(getRversion() >= "2.15.1")
    utils::globalVariables(

      #.calc_vd
      c("an")
    )
}
