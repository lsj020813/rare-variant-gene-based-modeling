#!/usr/bin/env Rscript
.libPaths("~/R/lib"); suppressMessages(library(susieR))
utils::assignInNamespace("is_symmetric_matrix", function(x) TRUE, ns = "susieR")
a <- commandArgs(TRUE); T <- a[1]; N <- a[2]; RID <- a[3]; ST <- a[4]; LD <- a[5]; OUT <- a[6]; PHI <- a[7]
ARMS <- strsplit(a[8], ",")[[1]]; NPERM <- as.integer(a[9])
MI <- as.integer(Sys.getenv("FM_MAX_ITER", "100"))
UREF <- Sys.getenv("FM_UNIFORM_PIP", "")
s <- read.delim(ST, check.names = FALSE, colClasses = c(CHR = "character"), stringsAsFactors = FALSE)
key_s <- paste(sub("^chr", "", s$CHR), s$POS, s$Allele1, s$Allele2, sep = ":")
stopifnot(!any(duplicated(key_s)))
key_p <- readLines(paste0(LD, ".vars")); m <- length(key_p); M_full <- m
MMAX <- as.integer(Sys.getenv("FM_MAXVAR_MEM", "45000"))
if (m > MMAX) { cat(sprintf("PRIOR_SKIP mem m=%d > %d\n", m, MMAX)); quit(save = "no", status = 0) }
idx <- match(key_p, key_s); keep <- !is.na(idx); n_drop <- sum(!keep)
if (n_drop > 0) { cat(sprintf("NOTE %s: LD 변이 %d 개가 step2 에 없음 -> 제외(부분행렬)\n", RID, n_drop)); key_p <- key_p[keep]; idx <- idx[keep]; m <- length(key_p) }
s2 <- s[idx, ]
z <- s2$BETA / s2$SE; if (any(!is.finite(z))) stop("GATE FAIL: z 결측")
n <- if ("N" %in% names(s2)) as.integer(round(median(s2$N))) else as.integer(round(median(s2$N_case + s2$N_ctrl)))
adj <- (n - 1) / (z^2 + n - 2); z_adj <- sqrt(adj) * z; Xty <- sqrt(n - 1) * z_adj
pi_ <- NULL
if (any(c("phi", "perm") %in% ARMS)) {
  ph <- read.delim(PHI, stringsAsFactors = FALSE)
  stopifnot(all(c("variant_hg19", "phi") %in% names(ph)))
  pi_ <- ph$phi[match(key_p, ph$variant_hg19)]
  if (any(is.na(pi_))) stop(sprintf("GATE FAIL: φ 결측 %d", sum(is.na(pi_))))
  if (any(pi_ <= 0) || any(!is.finite(pi_))) stop("GATE FAIL: φ <= 0 또는 비유한")
}
readR <- function() { con <- file(paste0(LD, ".bin"), "rb"); R <- readBin(con, "numeric", n = M_full * M_full, size = 4); close(con)
  if (length(R) != M_full * M_full) stop("GATE FAIL: LD bin 크기"); dim(R) <- c(M_full, M_full)
  if (n_drop > 0) { R <- R[keep, keep, drop = FALSE]; invisible(gc()) }
  R }
R <- readR()
if (any(!is.finite(R))) stop("GATE FAIL: R 결측"); if (any(abs(diag(R) - 1) > 1e-4)) stop("GATE FAIL: LD 대각 != 1")
k <- seq_len(min(300L, m)); if (!isTRUE(all.equal(R[k, ], t(R[, k]), tolerance = 1e-6))) stop("GATE FAIL: R 비대칭")
rm(R); invisible(gc())
ent <- function(p) { q <- p / sum(p); q <- q[q > 0]; -sum(q * log(q)) }
fit_arm <- function(w) {
  t0 <- Sys.time()
  R <- readR(); R <- R * (n - 1); invisible(gc())
  set.seed(1)
  fit <- susie_suff_stat(XtX = R, Xty = Xty, n = n, yty = n - 1, L = 10, standardize = FALSE, estimate_residual_variance = FALSE,
                         check_prior = TRUE, coverage = NULL, min_abs_corr = 0.5, max_iter = MI, prior_weights = w)
  rm(R); invisible(gc())
  R <- readR()
  fit$sets <- susie_get_cs(fit, Xcorr = R, coverage = 0.95, min_abs_corr = 0.5, check_symmetric = FALSE)
  fit$pip <- susie_get_pip(fit, prune_by_cs = FALSE, prior_tol = 1e-9)
  rm(R); invisible(gc())
  cs_id <- rep(NA_integer_, m); sizes <- integer(0); pur <- numeric(0)
  if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 0) { for (j in seq_along(fit$sets$cs)) cs_id[fit$sets$cs[[j]]] <- j
    sizes <- sapply(fit$sets$cs, length); pur <- fit$sets$purity$min.abs.corr }
  list(pip = fit$pip, cs_id = cs_id, sizes = sizes, pur = pur, converged = fit$converged, niter = fit$niter,
       elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
summ_row <- function(arm, seed, r, extra = list()) {
  data.frame(trait = T, chr = N, region = RID, arm = arm, seed = seed, n_fm = m, n_cs = length(r$sizes),
             cs_sizes = paste(r$sizes, collapse = ";"), mean_cs_size = if (length(r$sizes)) signif(mean(r$sizes), 5) else NA,
             sum_cs_size = sum(r$sizes), cs_min_abs_corr = paste(signif(r$pur, 3), collapse = ";"),
             max_pip = signif(max(r$pip), 5), entropy = signif(ent(r$pip), 5), n_pip_ge0.1 = sum(r$pip >= 0.1), n_pip_ge0.5 = sum(r$pip >= 0.5),
             sum_pip = signif(sum(r$pip), 5), converged = r$converged, niter = r$niter, max_iter = MI, elapsed_s = round(r$elapsed, 1),
             gate_max_abs_dpip = if (is.null(extra$dpip)) NA else signif(extra$dpip, 4), gate_cs_equal = if (is.null(extra$cseq)) NA else extra$cseq,
             stringsAsFactors = FALSE)
}
rows <- list()
uref <- NULL
if (nzchar(UREF) && file.exists(UREF)) {
  uref <- read.delim(UREF, stringsAsFactors = FALSE)
  if (!identical(uref$variant_hg19, key_p)) stop("GATE FAIL: uniform pip.tsv 변이 순서 불일치")
  us <- as.integer(table(uref$cs_id)); 
  rows[[length(rows) + 1]] <- summ_row("uniform_ref", NA, list(pip = uref$pip, sizes = us, pur = numeric(0), converged = NA, niter = NA, elapsed = 0))
}
writePip <- function(r, suffix) {
  out <- data.frame(variant_hg19 = key_p, pip = signif(r$pip, 6), cs_id = r$cs_id, region = RID, stringsAsFactors = FALSE)
  write.table(out, paste0(OUT, ".", suffix, ".pip.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
if ("gate" %in% ARMS) {
  r <- fit_arm(rep(1 / m, m))
  ex <- list()
  if (!is.null(uref)) { ex$dpip <- max(abs(r$pip - uref$pip)); ex$cseq <- identical(sort(as.integer(table(uref$cs_id))), sort(as.integer(r$sizes))) }
  rows[[length(rows) + 1]] <- summ_row("gate_uniform", NA, r, ex); writePip(r, "gate")
  cat(sprintf("GATE %s n_cs=%d max_abs_dpip=%s cs_equal=%s\n", RID, length(r$sizes), format(ex$dpip), format(ex$cseq)))
}
if ("phi" %in% ARMS) {
  r <- fit_arm(pi_ / sum(pi_)); rows[[length(rows) + 1]] <- summ_row("phi", NA, r); writePip(r, "phi")
  cat(sprintf("PHI %s n_cs=%d max_pip=%.3f niter=%d %.0fs\n", RID, length(r$sizes), max(r$pip), r$niter, r$elapsed))
}
if ("perm" %in% ARMS) for (kk in seq_len(NPERM)) {
  set.seed(20260910 + kk); w <- sample(pi_); w <- w / sum(w)
  r <- fit_arm(w); rows[[length(rows) + 1]] <- summ_row(sprintf("perm%02d", kk), 20260910 + kk, r)
  cat(sprintf("PERM %s k=%d n_cs=%d max_pip=%.3f niter=%d %.0fs\n", RID, kk, length(r$sizes), max(r$pip), r$niter, r$elapsed))
  write.table(do.call(rbind, rows), paste0(OUT, ".arms.tsv.partial"), sep = "\t", quote = FALSE, row.names = FALSE)
}
write.table(do.call(rbind, rows), paste0(OUT, ".arms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
unlink(paste0(OUT, ".arms.tsv.partial"))
cat(sprintf("PRIOR_OK %s arms=%d m=%d\n", RID, length(rows), m))
