#!/usr/bin/env Rscript
.libPaths("~/R/lib"); suppressMessages(library(susieR))
a <- commandArgs(TRUE); ST <- a[1]; LD <- a[2]
s <- read.delim(ST, check.names = FALSE, colClasses = c(CHR = "character"), stringsAsFactors = FALSE)
key_s <- paste(sub("^chr", "", s$CHR), s$POS, s$Allele1, s$Allele2, sep = ":")
key_p <- readLines(paste0(LD, ".vars")); m <- length(key_p)
idx <- match(key_p, key_s); stopifnot(!any(is.na(idx)))
s2 <- s[idx, ]; z <- s2$BETA / s2$SE
n <- if ("N" %in% names(s2)) as.integer(round(median(s2$N))) else as.integer(round(median(s2$N_case + s2$N_ctrl)))
con <- file(paste0(LD, ".bin"), "rb"); R0 <- readBin(con, "numeric", n = m * m, size = 4); close(con); dim(R0) <- c(m, m)
R0 <- (R0 + t(R0)) / 2; diag(R0) <- 1
cat(sprintf("EQ_INPUT m=%d n=%d max|z|=%.2f\n", m, n, max(abs(z))))

set.seed(1)
tA <- system.time(fA <- susie_rss(z = z, R = R0, n = n, L = 10, coverage = 0.95, min_abs_corr = 0.5, max_iter = 100))
pA <- susie_get_pip(fA, prune_by_cs = FALSE, prior_tol = 1e-9)
csA <- if (is.null(fA$sets$cs)) 0L else length(fA$sets$cs)

utils::assignInNamespace("is_symmetric_matrix", function(x) TRUE, ns = "susieR")
set.seed(1)
adj <- (n - 1) / (z^2 + n - 2); z_adj <- sqrt(adj) * z
R <- R0 * (n - 1)
tB <- system.time(fB <- susie_suff_stat(XtX = R, Xty = sqrt(n - 1) * z_adj, n = n, yty = n - 1, L = 10,
        standardize = FALSE, estimate_residual_variance = FALSE, check_prior = TRUE, coverage = NULL,
        min_abs_corr = 0.5, max_iter = 100))
rm(R); invisible(gc())
fB$sets <- susie_get_cs(fB, Xcorr = R0, coverage = 0.95, min_abs_corr = 0.5, check_symmetric = FALSE)
pB <- susie_get_pip(fB, prune_by_cs = FALSE, prior_tol = 1e-9)
csB <- if (is.null(fB$sets$cs)) 0L else length(fB$sets$cs)

d <- max(abs(pA - pB)); sp <- suppressWarnings(cor(pA, pB, method = "spearman"))
setsA <- if (csA) sort(sapply(fA$sets$cs, function(v) paste(sort(v), collapse = ","))) else character(0)
setsB <- if (csB) sort(sapply(fB$sets$cs, function(v) paste(sort(v), collapse = ","))) else character(0)
same <- identical(setsA, setsB)
cat(sprintf("EQ_RESULT maxPIPdiff=%.3g spearman=%.6f nCS_A=%d nCS_B=%d CS_identical=%s conv_A=%s conv_B=%s niter_A=%d niter_B=%d sec_A=%.0f sec_B=%.0f maxPIP_A=%.4f maxPIP_B=%.4f\n",
    d, sp, csA, csB, same, fA$converged, fB$converged, fA$niter, fB$niter, tA[["elapsed"]], tB[["elapsed"]], max(pA), max(pB)))
cat(if (d < 1e-6 && same) "EQ_PASS\n" else "EQ_DIFF\n")
