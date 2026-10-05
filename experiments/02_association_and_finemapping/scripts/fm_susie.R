#!/usr/bin/env Rscript
.libPaths("~/R/lib"); suppressMessages(library(susieR))
utils::assignInNamespace("is_symmetric_matrix", function(x) TRUE, ns = "susieR")
a <- commandArgs(TRUE); T <- a[1]; N <- a[2]; RID <- a[3]; ST <- a[4]; LD <- a[5]; OUT <- a[6]
s <- read.delim(ST, check.names = FALSE, colClasses = c(CHR = "character"), stringsAsFactors = FALSE)
key_s <- paste(sub("^chr", "", s$CHR), s$POS, s$Allele1, s$Allele2, sep = ":")
stopifnot(!any(duplicated(key_s)))
key_p <- readLines(paste0(LD, ".vars"))
pv <- data.frame(k = key_p, stringsAsFactors = FALSE)
m <- nrow(pv)
MMAX <- as.integer(Sys.getenv("FM_MAXVAR_MEM", "27000"))
if (m > MMAX) { cat(sprintf("SUSIE_SKIP mem m=%d > %d\n", m, MMAX)); quit(save = "no", status = 0) }
idx <- match(key_p, key_s)
if (any(is.na(idx))) stop(sprintf("GATE FAIL: LD 변이 %d 개가 step2 에 없음", sum(is.na(idx))))
n_ld <- m; n_step2 <- nrow(s); n_fm <- m
if (n_fm < 2) stop("GATE FAIL: 변이 <2")
s2 <- s[idx, ]
con <- file(paste0(LD, ".bin"), "rb"); R <- readBin(con, "numeric", n = m * m, size = 4); close(con)
if (length(R) != m * m) stop(sprintf("GATE FAIL: LD bin 크기 %d != %d^2", length(R), m))
dim(R) <- c(m, m)
z <- s2$BETA / s2$SE
if (any(!is.finite(z))) stop(sprintf("GATE FAIL: z 결측 %d", sum(!is.finite(z))))
if (any(!is.finite(R))) stop(sprintf("GATE FAIL: R 결측 %d", sum(!is.finite(R))))
if (any(abs(diag(R) - 1) > 1e-4)) stop("GATE FAIL: LD 대각 != 1")
k <- seq_len(min(300L, m)); if (!isTRUE(all.equal(R[k, ], t(R[, k]), tolerance = 1e-6))) stop("GATE FAIL: R 비대칭")
invisible(gc())
n <- if ("N" %in% names(s2)) as.integer(round(median(s2$N))) else as.integer(round(median(s2$N_case + s2$N_ctrl)))
lam <- if (nzchar(Sys.getenv("FM_LAMBDA"))) tryCatch(estimate_s_rss(z, R, n = n), error = function(e) NA_real_) else NA_real_
set.seed(1)
adj <- (n - 1) / (z^2 + n - 2); z_adj <- sqrt(adj) * z
R <- R * (n - 1); invisible(gc())
fit <- susie_suff_stat(XtX = R, Xty = sqrt(n - 1) * z_adj, n = n, yty = n - 1, L = 10, standardize = FALSE,
                       estimate_residual_variance = FALSE, check_prior = TRUE, coverage = NULL, min_abs_corr = 0.5, max_iter = as.integer(Sys.getenv("FM_MAX_ITER", "100")))
rm(R); invisible(gc())
con <- file(paste0(LD, ".bin"), "rb"); R <- readBin(con, "numeric", n = m * m, size = 4); close(con); dim(R) <- c(m, m)
fit$sets <- susie_get_cs(fit, Xcorr = R, coverage = 0.95, min_abs_corr = 0.5, check_symmetric = FALSE)
fit$pip <- susie_get_pip(fit, prune_by_cs = FALSE, prior_tol = 1e-9)
rm(R); invisible(gc())
pip <- fit$pip
cs_id <- rep(NA_integer_, n_fm)
ncs <- 0L; cs_sizes <- integer(0); cs_purity <- numeric(0)
if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 0) {
  ncs <- length(fit$sets$cs)
  for (k in seq_along(fit$sets$cs)) { cs_id[fit$sets$cs[[k]]] <- k }
  cs_sizes <- sapply(fit$sets$cs, length); cs_purity <- fit$sets$purity$min.abs.corr
}
af <- s2$AF_Allele2; maf <- pmin(af, 1 - af)
band <- ifelse(maf >= 0.05, "ge5", ifelse(maf >= 0.01, "1to5", ifelse(maf >= 0.001, "0.1to1", "lt0.1")))
out <- data.frame(variant_hg19 = key_p, chr = sub("^chr", "", s2$CHR), pos = s2$POS, a1 = s2$Allele1, a2 = s2$Allele2,
                  maf = signif(maf, 6), z = signif(z, 6), pip = signif(pip, 6), cs_id = cs_id, region = RID, stringsAsFactors = FALSE)
write.table(out, paste0(OUT, ".pip.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
mass <- tapply(pip, factor(band, levels = c("ge5", "1to5", "0.1to1", "lt0.1")), sum); mass[is.na(mass)] <- 0
nvar <- table(factor(band, levels = c("ge5", "1to5", "0.1to1", "lt0.1")))
summ <- data.frame(trait = T, chr = N, region = RID, n = n, n_ld = n_ld, n_step2 = n_step2, n_fm = n_fm,
                   n_cs = ncs, cs_sizes = paste(cs_sizes, collapse = ";"), cs_min_abs_corr = paste(signif(cs_purity, 3), collapse = ";"),
                   n_pip_ge0.1 = sum(pip >= 0.1), n_pip_ge0.5 = sum(pip >= 0.5), max_pip = signif(max(pip), 4),
                   nvar_ge5 = nvar[["ge5"]], nvar_1to5 = nvar[["1to5"]], nvar_0.1to1 = nvar[["0.1to1"]], nvar_lt0.1 = nvar[["lt0.1"]],
                   mass_ge5 = signif(mass[["ge5"]], 4), mass_1to5 = signif(mass[["1to5"]], 4), mass_0.1to1 = signif(mass[["0.1to1"]], 4), mass_lt0.1 = signif(mass[["lt0.1"]], 4),
                   converged = fit$converged, niter = fit$niter, lambda_s = signif(lam, 4), max_abs_z = signif(max(abs(z)), 4))
write.table(summ, paste0(OUT, ".summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("SUSIE_OK %s n_fm=%d n_cs=%d n_pip01=%d converged=%s lambda=%.3g\n", RID, n_fm, ncs, sum(pip >= 0.1), fit$converged, lam))
