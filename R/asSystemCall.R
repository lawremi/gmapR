##written to get the string of the system call rather than having the
##system call actually executed. Initial motivation was for unit
##testing.
asSystemCall <- function(x) {
  expr <- substitute(x)
  options(systemCallMode = TRUE)
  on.exit(options(systemCallMode = FALSE))
  error <- tryCatch(eval(expr, parent.frame()),
                         error = function(e) return(e))
  error$systemCall
}
