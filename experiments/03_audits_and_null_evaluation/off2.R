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

outdir <- configured("${PROJECT_ROOT}/work/ref/l0")
for (t in c("tchl","htn","dm","lip")) {
  o <- get(load(sprintf(configured("${PROJECT_ROOT}/work/ref/saige_step1_v4/%s_v4.rda"), t)))
  df <- data.frame(sampleID=o$sampleID, eta=o$linear.predictors, y=o$y)
  fp <- sprintf("%s/%s.offset.tsv", outdir, t)
  write.table(df, fp, sep="\t", quote=FALSE, row.names=FALSE)
  cat(t, "n=", nrow(df), " traitType=", o$traitType, "\n")
}
cat("OFFSETS_DONE\n")
