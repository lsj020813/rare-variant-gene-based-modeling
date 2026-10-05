.required_env <- function(name) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) stop(name, " must be set", call. = FALSE)
  value
}
load(paste0(.required_env("PROJECT_ROOT"), '/work/ref/saige_step1_v4/dm_v4.rda')); cat('objects:', ls(), '\n'); m <- get(ls()[1]); cat('n fields:', length(names(m)), '\n'); print(names(m))
for (f in c('linear.predictors','fitted.values','residuals','offset','y','mu','eta','tau','sampleID','traitType','obj.noK','X','coefficients')) if (!is.null(m[[f]])) cat(sprintf('%-18s %s len=%d\n', f, class(m[[f]])[1], length(m[[f]])))
