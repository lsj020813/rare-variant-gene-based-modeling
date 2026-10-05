
model_file <- Sys.getenv("NULL_MODEL_FILE", unset = "")
if (!nzchar(model_file)) stop("Required environment variable missing: NULL_MODEL_FILE")
if (!file.exists(model_file)) stop("NULL_MODEL_FILE does not exist")
obj <- get(load(model_file))
cat("names:", paste(names(obj), collapse=", "), "\n")
for (k in c("traitType","X.colnames","covariateNames","obj.noK","sampleID","y")) {
  if (!is.null(obj[[k]])) {
    v <- obj[[k]]
    if (is.character(v) && length(v) < 40) cat(k, ":", paste(v, collapse=","), "\n")
    else cat(k, ": length", length(v), "class", class(v)[1], "\n")
  }
}
if (!is.null(obj$X)) cat("X colnames:", paste(colnames(obj$X), collapse=","), " nrow:", nrow(obj$X), "\n")
if (!is.null(obj$sampleID)) cat("N samples:", length(obj$sampleID), "\n")
if (!is.null(obj$y)) cat("y summary: mean", mean(obj$y), "sd", sd(obj$y), "\n")
