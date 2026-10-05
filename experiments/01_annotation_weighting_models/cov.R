.required_env <- function(name) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) stop(name, " must be set", call. = FALSE)
  value
}

for (t in c("arth_cur","gastro_cur","tchl_v4")) {
  o <- get(load(sprintf(paste0(.required_env("PROJECT_ROOT"), "/work/ref/saige_step1/%s.rda"), t)))
  cn <- colnames(o$X)
  cat(t, "| N:", length(o$sampleID), "| traitType:", o$traitType, "| covars:", paste(cn, collapse=","), "\n")
}
