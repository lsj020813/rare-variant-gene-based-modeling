
suppressMessages(library(SAIGE))
src <- unlist(strsplit(paste(deparse(get("extractVarianceRatio", envir=getNamespace("SAIGE"))), collapse="\n"), "\n"))
idx <- grep("setgeno|minMAFforGRM|isVarRatio|getIsVarRatio|memoryChunk|Geno_forVarRatio", src)
for (i in idx) cat(sprintf("%4d| %s\n", i, substr(src[i],1,150)))
cat("\n=== fitNULLGLMM default args of interest ===\n")
fa <- formals(get("fitNULLGLMM", envir=getNamespace("SAIGE")))
for (k in c("minMAFforGRM","maxMissingRateforGRM","numRandomMarkerforVarianceRatio",
            "ratioCVcutoff","isCateVarianceRatio","cateVarRatioMinMACVecExclude",
            "cateVarRatioMaxMACVecInclude","includeNonautoMarkersforVarRatio")) {
  if (k %in% names(fa)) cat(sprintf("  %-34s %s\n", k, paste(deparse(fa[[k]]), collapse="")))
}
