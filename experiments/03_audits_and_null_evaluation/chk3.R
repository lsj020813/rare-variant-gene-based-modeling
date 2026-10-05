configured <- function(value) {
  matches <- regmatches(value, gregexpr("\\$\\{[A-Z][A-Z0-9_]*\\}", value))[[1]]
  for (token in matches) {
    name <- substring(token, 3, nchar(token)-1)
    replacement <- Sys.getenv(name)
    if (!nzchar(replacement)) stop(paste("Set", name, "before running this script"))
    value <- gsub(token, replacement, value, fixed=TRUE)
  }
  value
}

for (t in c("htn","dm","lip","tchl")) {
  p <- paste0(configured("${PROJECT_ROOT}/work/ref/saige_step1_v3/"), t, "_v3.rda")
  o <- get(load(p))
  xn <- NULL
  if (!is.null(o$X) && !is.null(colnames(o$X))) xn <- colnames(o$X)
  cat(sprintf("%-5s N=%-7d trait=%-13s theta=%s\n  covars: %s\n",
      t, length(o$sampleID), o$traitType,
      paste(round(o$theta,5), collapse=","),
      if (is.null(xn)) "(none)" else paste(xn, collapse=",")))
}
