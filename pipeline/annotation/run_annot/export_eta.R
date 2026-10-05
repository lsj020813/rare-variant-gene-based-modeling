.required_env <- function(name) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) stop(name, " must be set", call. = FALSE)
  value
}
for (ph in c('tchl','htn','dm','lip')) {
  load(sprintf(paste0(.required_env("PROJECT_ROOT"), '/work/ref/saige_step1_v4/%s_v4.rda'), ph)); m <- modglmm
  n <- length(m$sampleID)
  stopifnot(length(m$y)==n, nrow(m$linear.predictors)==n, nrow(m$residuals)==n)
  d <- data.frame(IID=m$sampleID, y=as.numeric(m$y), eta=as.numeric(m$linear.predictors), mu=as.numeric(m$fitted.values),
                  resid=as.numeric(m$residuals), offset=as.numeric(m$offset), stringsAsFactors=FALSE)
  write.table(d, sprintf(paste0(.required_env("PROJECT_ROOT"), '/work/ref/annot/offset/%s.eta.tsv'), ph), sep='\t', quote=FALSE, row.names=FALSE)
  cat(sprintf('%s traitType=%s n=%d isCovOffset=%s tau=%s ymean=%.4f eta_range=[%.3f,%.3f]\n', ph, m$traitType, n, m$isCovariateOffset,
      paste(round(m$theta,4), collapse='/'), mean(m$y), min(m$linear.predictors), max(m$linear.predictors)))
}
