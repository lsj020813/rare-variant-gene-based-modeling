
suppressMessages(library(SAIGE))
fns <- ls(getNamespace("SAIGE"))
cand <- grep("VarRatio|VarianceRatio|varRatio|Cate", fns, value=TRUE, ignore.case=TRUE)
cat("candidate functions:", paste(cand, collapse=", "), "\n\n")
for (f in cand) {
  src <- tryCatch(paste(deparse(get(f, envir=getNamespace("SAIGE"))), collapse="\n"),
                  error=function(e) "")
  hits <- grep("200|MAC|categ|at least", strsplit(src,"\n")[[1]], value=TRUE)
  if (length(hits)) {
    cat("=====", f, "=====\n")
    cat(paste(head(hits, 14), collapse="\n"), "\n\n")
  }
}
