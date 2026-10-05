
suppressMessages(library(SAIGE))
src <- unlist(strsplit(paste(deparse(get("extractVarianceRatio", envir=getNamespace("SAIGE"))), collapse="\n"), "\n"))
idx <- grep("ratio|Ratio", src)
cat("=== how the ratio is computed ===\n")
for (i in idx) {
  l <- src[i]
  if (grepl("var1|var2|ratioVec|=.*\\/|mean\\(", l)) cat(sprintf("%4d| %s\n", i, substr(l,1,150)))
}
