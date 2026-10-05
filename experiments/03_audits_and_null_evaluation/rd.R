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
o <- get(load(configured("${ANALYSIS_ROOT}/work/saige_gene/grch38_chr22/step1/tchl_primary.rda")))
cat("FIELDS:", paste(names(o), collapse=","), "\n")
if (!is.null(o$traitType)) cat("traitType:", o$traitType, "\n")
if (!is.null(o$sampleID))  cat("N:", length(o$sampleID), "\n")
xn <- NULL
for (k in c("X","X1","Xmat")) if (!is.null(o[[k]]) && !is.null(colnames(o[[k]]))) xn <- colnames(o[[k]])
if (is.null(xn) && !is.null(o$obj.noK$X1)) xn <- colnames(o$obj.noK$X1)
cat("COVARIATES:", if (is.null(xn)) "(not stored)" else paste(xn, collapse=","), "\n")
if (!is.null(o$theta)) cat("theta:", paste(round(o$theta,5), collapse=","), "\n")
