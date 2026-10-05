#!/usr/bin/env Rscript
.libPaths("~/R/lib"); suppressMessages(library(susieR))
utils::assignInNamespace("is_symmetric_matrix", function(x) TRUE, ns = "susieR")
a <- commandArgs(TRUE); T <- a[1]; N <- a[2]; RID <- a[3]; ST <- a[4]; LD <- a[5]; OUT <- a[6]; PHI <- a[7]; NPERM <- as.integer(a[8])
MI <- as.integer(Sys.getenv("FM_MAX_ITER", "100")); UREF <- Sys.getenv("FM_UNIFORM_PIP", "")
s <- read.delim(ST, check.names = FALSE, colClasses = c(CHR = "character"), stringsAsFactors = FALSE)
key_s <- paste(sub("^chr", "", s$CHR), s$POS, s$Allele1, s$Allele2, sep = ":"); stopifnot(!any(duplicated(key_s)))
key_p <- readLines(paste0(LD, ".vars")); m <- length(key_p); M_full <- m
MMAX <- as.integer(Sys.getenv("FM_MAXVAR_MEM", "45000"))
if (m > MMAX) { cat(sprintf("RW_SKIP mem m=%d > %d\n", m, MMAX)); quit(save = "no", status = 0) }
idx <- match(key_p, key_s); keep <- !is.na(idx); n_drop <- sum(!keep)
if (n_drop > 0) { cat(sprintf("NOTE %s: LD 변이 %d 개가 step2 에 없음 -> 제외(부분행렬)\n", RID, n_drop)); key_p <- key_p[keep]; idx <- idx[keep]; m <- length(key_p) }
s2 <- s[idx, ]; z <- s2$BETA / s2$SE; if (any(!is.finite(z))) stop("GATE FAIL: z 결측")
n <- if ("N" %in% names(s2)) as.integer(round(median(s2$N))) else as.integer(round(median(s2$N_case + s2$N_ctrl)))
adj <- (n - 1) / (z^2 + n - 2); z_adj <- sqrt(adj) * z; Xty <- sqrt(n - 1) * z_adj
ph <- read.delim(PHI, stringsAsFactors = FALSE); pi_ <- ph$phi[match(key_p, ph$variant_hg19)]
if (any(is.na(pi_)) || any(pi_ <= 0)) stop("GATE FAIL: φ 결측/비양수")
readR <- function() { con <- file(paste0(LD, ".bin"), "rb"); R <- readBin(con, "numeric", n = M_full * M_full, size = 4); close(con)
  if (length(R) != M_full * M_full) stop("GATE FAIL: LD bin 크기"); dim(R) <- c(M_full, M_full); if (n_drop > 0) { R <- R[keep, keep, drop = FALSE]; invisible(gc()) }; R }
t0 <- Sys.time()
R <- readR(); if (any(abs(diag(R) - 1) > 1e-4)) stop("GATE FAIL: LD 대각"); R <- R * (n - 1); invisible(gc())
set.seed(1)
fit <- susie_suff_stat(XtX = R, Xty = Xty, n = n, yty = n - 1, L = 10, standardize = FALSE, estimate_residual_variance = FALSE,
                       check_prior = TRUE, coverage = NULL, min_abs_corr = 0.5, max_iter = MI)
rm(R); invisible(gc()); R <- readR()
getcs <- function(f, mac = 0.5) { set.seed(1); susie_get_cs(f, Xcorr = R, coverage = 0.95, min_abs_corr = mac, check_symmetric = FALSE) }
fit$sets <- getcs(fit); fit$pip <- susie_get_pip(fit, prune_by_cs = FALSE, prior_tol = 1e-9)
sets_all <- getcs(fit, 0)
fit_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cs_index <- if (is.null(fit$sets$cs_index)) integer(0) else fit$sets$cs_index
eff_all <- which(fit$V > 1e-9)
saveRDS(list(region = RID, key = key_p, alpha = fit$alpha, V = fit$V, pip = fit$pip, sets = fit$sets, n_cs_all = length(sets_all$cs), cs_index = cs_index, niter = fit$niter, converged = fit$converged), paste0(OUT, ".alpha.rds"))
ent <- function(p) { q <- p / sum(p); q <- q[q > 0]; -sum(q * log(q)) }
metr <- function(pip, sets) { sz <- if (length(sets$cs)) sapply(sets$cs, length) else integer(0)
  list(n_cs = length(sz), cs_sizes = paste(sz, collapse = ";"), mean_cs = if (length(sz)) mean(sz) else NA, sum_cs = sum(sz), max_pip = max(pip), entropy = ent(pip),
       n01 = sum(pip >= 0.1), n05 = sum(pip >= 0.5), purity = paste(signif(sets$purity$min.abs.corr, 3), collapse = ";")) }
rw <- function(w, effects) { a2 <- fit$alpha
  for (l in effects) { v <- w * fit$alpha[l, ]; a2[l, ] <- v / sum(v) }
  f2 <- fit; f2$alpha <- a2; sets <- getcs(f2); pip <- susie_get_pip(f2, prune_by_cs = FALSE, prior_tol = 1e-9); list(pip = pip, sets = sets, alpha = a2) }
row <- function(arm, seed, mt, extra = list()) data.frame(trait = T, chr = N, region = RID, arm = arm, seed = seed, n_fm = m, n_cs = mt$n_cs, cs_sizes = mt$cs_sizes,
  mean_cs_size = signif(mt$mean_cs, 5), sum_cs_size = mt$sum_cs, cs_min_abs_corr = mt$purity, max_pip = signif(mt$max_pip, 5), entropy = signif(mt$entropy, 5),
  n_pip_ge0.1 = mt$n01, n_pip_ge0.5 = mt$n05, n_cs_all = if (is.null(extra$n_cs_all)) NA else extra$n_cs_all, n_eff_reweighted = if (is.null(extra$neff)) NA else extra$neff,
  converged = fit$converged, niter = fit$niter, max_iter = MI, fit_s = round(fit_s, 1),
  gate_max_abs_dpip = if (is.null(extra$dpip)) NA else signif(extra$dpip, 4), stringsAsFactors = FALSE)
rows <- list()
mu <- metr(fit$pip, fit$sets)
ex <- list(n_cs_all = length(sets_all$cs), neff = 0L)
if (nzchar(UREF) && file.exists(UREF)) { u <- read.delim(UREF, stringsAsFactors = FALSE); if (!identical(u$variant_hg19, key_p)) stop("GATE FAIL: uniform pip 순서"); ex$dpip <- max(abs(fit$pip - u$pip)) }
rows[[1]] <- row("uniform_refit", NA, mu, ex)
cat(sprintf("GATE %s n_cs=%d n_cs_all=%d cs_index=%s max_abs_dpip=%s fit=%.0fs\n", RID, mu$n_cs, length(sets_all$cs), paste(cs_index, collapse = ","), format(ex$dpip), fit_s))
A <- rw(pi_, cs_index); B <- rw(pi_, eff_all)
rows[[2]] <- row("phi_rwA", NA, metr(A$pip, A$sets), list(neff = length(cs_index), n_cs_all = length(getcs(list(alpha = A$alpha, V = fit$V, null_index = fit$null_index), 0)$cs)))
rows[[3]] <- row("phi_rwB", NA, metr(B$pip, B$sets), list(neff = length(eff_all), n_cs_all = length(getcs(list(alpha = B$alpha, V = fit$V, null_index = fit$null_index), 0)$cs)))
wp <- function(x, suf) write.table(data.frame(variant_hg19 = key_p, pip = signif(x$pip, 6), cs_id = { ci <- rep(NA_integer_, m); if (length(x$sets$cs)) for (j in seq_along(x$sets$cs)) ci[x$sets$cs[[j]]] <- j; ci }, region = RID),
  paste0(OUT, ".", suf, ".pip.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
wp(A, "rwA"); wp(B, "rwB"); wp(list(pip = fit$pip, sets = fit$sets), "rwU")
cat(sprintf("PHI_RW %s A: n_cs=%d max_pip=%.3f | B: n_cs=%d max_pip=%.3f\n", RID, length(A$sets$cs), max(A$pip), length(B$sets$cs), max(B$pip)))
t1 <- Sys.time()
pm <- vector("list", NPERM)
for (k in seq_len(NPERM)) { set.seed(20260910 + k); w <- sample(pi_)
  a_ <- rw(w, cs_index); b_ <- rw(w, eff_all); ma <- metr(a_$pip, a_$sets); mb <- metr(b_$pip, b_$sets)
  pm[[k]] <- data.frame(region = RID, k = k, seed = 20260910 + k, A_n_cs = ma$n_cs, A_mean_cs = signif(ma$mean_cs, 5), A_sum_cs = ma$sum_cs, A_max_pip = signif(ma$max_pip, 5), A_entropy = signif(ma$entropy, 5),
                        B_n_cs = mb$n_cs, B_mean_cs = signif(mb$mean_cs, 5), B_sum_cs = mb$sum_cs, B_max_pip = signif(mb$max_pip, 5), B_entropy = signif(mb$entropy, 5)) }
perm_s <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
write.table(do.call(rbind, pm), paste0(OUT, ".rw.perms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, rows), paste0(OUT, ".rw.arms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("RW_OK %s m=%d nperm=%d fit=%.0fs perms=%.0fs\n", RID, m, NPERM, fit_s, perm_s))
