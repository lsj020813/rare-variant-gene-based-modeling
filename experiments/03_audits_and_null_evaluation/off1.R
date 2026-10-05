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

o <- get(load(configured("${PROJECT_ROOT}/work/ref/saige_step1_v4/tchl_v4.rda")))
cat("FIELDS:", paste(names(o), collapse=","), "\n")
for (f in c("linear.predictors","fitted.values","offset","eta","mu","y","residuals","sampleID")) {
  if (f %in% names(o)) {
    v <- o[[f]]
    cat(f, ": length", length(v), " head:", paste(head(v,3), collapse=" "), "\n")
  }
}
