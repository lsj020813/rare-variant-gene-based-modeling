
suppressMessages(library(SAIGE))
src <- unlist(strsplit(paste(deparse(get("extractVarianceRatio", envir=getNamespace("SAIGE"))), collapse="\n"), "\n"))
cat("=== lines 200-360 (the estimation loop) ===\n")
for (i in 200:min(360,length(src))) {
  l <- src[i]
  if (grepl("varRatio|Var|MAC|G0|scale|AC|sum\\(|if *\\(", l)) cat(sprintf("%4d| %s\n", i, substr(l,1,140)))
}
