.required_env <- function(name) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) stop(name, " must be set", call. = FALSE)
  value
}
for (t in c("tchl","htn","dm","lip")) {
  e <- new.env(); load(sprintf(paste0(.required_env("PROJECT_ROOT"), "/work/ref/saige_step1_v4/%s_v4.rda"), t), envir = e)
  o <- get(ls(e)[1], envir = e); ids <- o$sampleID
  stopifnot(length(ids) > 1000, !any(duplicated(ids)))
  writeLines(ids, sprintf(paste0(.required_env("PROJECT_ROOT"), "/work/run_ourfm/keep/%s.keep"), t))
  cat(t, length(ids), o$traitType, "\n")
}
