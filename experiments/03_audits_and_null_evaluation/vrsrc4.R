
suppressMessages(library(SAIGE))
for (fn in c("setgeno","getIsVarRatioGeno")) {
  f <- tryCatch(get(fn, envir=getNamespace("SAIGE")), error=function(e) NULL)
  if (is.null(f)) { cat(fn, ": absent\n"); next }
  s <- unlist(strsplit(paste(deparse(f), collapse="\n"), "\n"))
  cat("=====", fn, "(", length(s), "lines ) =====\n")
  idx <- grep("VarRatio|minMAF|MAC|isVarRatio|memoryChunk", s)
  for (i in head(idx, 18)) cat(sprintf("%4d| %s\n", i, substr(s[i],1,140)))
  cat("\n")
}
fa <- formals(get("fitNULLGLMM", envir=getNamespace("SAIGE")))
cat("=== args mentioning VarRatio/MAC ===\n")
for (k in names(fa)) if (grepl("VarRatio|MAC|minMAF", k))
  cat(sprintf("  %-36s %s\n", k, paste(deparse(fa[[k]]), collapse="")))
