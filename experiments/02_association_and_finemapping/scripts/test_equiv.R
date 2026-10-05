.required_env <- function(name) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) stop(name, " must be set", call. = FALSE)
  value
}
.libPaths("~/R/lib"); suppressMessages(library(susieR))
D <- paste0(.required_env("PROJECT_ROOT"), "/work/run_ourfm/smoke_chr22/fm/tchl/tchl_chr22_r1")
s <- read.delim(paste0(D, ".step2.txt"), check.names=FALSE, colClasses=c(CHR="character"))
key_s <- paste(s$CHR, s$POS, s$Allele1, s$Allele2, sep=":")
key_p <- readLines(paste0(D, ".ld.vars")); m <- length(key_p)
con <- file(paste0(D, ".ld.bin"), "rb"); Rm <- readBin(con, "numeric", n=m*m, size=4); close(con); dim(Rm) <- c(m,m)
idx <- match(key_p, key_s); keep <- !is.na(idx); R <- Rm[keep,keep]; rm(Rm); s2 <- s[idx[keep],]; z <- s2$BETA/s2$SE; n <- as.integer(round(median(s2$N)))
cat("m_fm", nrow(R), "n", n, "\n")
set.seed(1); t1 <- system.time(f1 <- susie_rss(z=z, R=R, n=n, L=10, coverage=0.95, min_abs_corr=0.5, max_iter=1000))
adj <- (n-1)/(z^2+n-2); z_adj <- sqrt(adj)*z; R2 <- R*(n-1)
set.seed(1); t2 <- system.time(f2 <- susie_suff_stat(XtX=R2, Xty=sqrt(n-1)*z_adj, n=n, yty=n-1, L=10, standardize=FALSE, estimate_residual_variance=FALSE, check_prior=TRUE, coverage=0.95, min_abs_corr=0.5, max_iter=1000))
cat("max|pip diff|", max(abs(f1$pip - f2$pip)), " max pip", max(f1$pip), max(f2$pip), " ncs", length(f1$sets$cs), length(f2$sets$cs), " niter", f1$niter, f2$niter, "\n")
cat("elbo diff", max(abs(f1$elbo - f2$elbo)), " time", t1[3], t2[3], "\n")
cat("EQUIV_DONE\n")
