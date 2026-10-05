
suppressMessages(library(SAIGE))
src <- deparse(get("extractVarianceRatio", envir=getNamespace("SAIGE")))
lines <- unlist(strsplit(paste(src, collapse="\n"), "\n"))
idx <- grep("MACdata|cateVarRatio|ratioVec|subMAC|nMarker|sample\\(|nrow\\(", lines)
cat("total lines:", length(lines), " matched:", length(idx), "\n\n")
sel <- sort(unique(unlist(lapply(idx, function(i) max(1,i-1):min(length(lines), i+2)))))
prev <- 0
for (i in sel) {
  if (i > prev + 1) cat("   ...\n")
  cat(sprintf("%4d| %s\n", i, substr(lines[i], 1, 150)))
  prev <- i
}
