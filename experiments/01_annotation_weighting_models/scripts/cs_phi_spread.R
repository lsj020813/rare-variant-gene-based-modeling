#!/usr/bin/env Rscript
.required_env <- function(name) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) stop(name, " must be set", call. = FALSE)
  value
}
.libPaths("~/R/lib"); suppressMessages(library(susieR))
W <- paste0(.required_env("PROJECT_ROOT"), "/work/run_trackA"); F0 <- paste0(.required_env("PROJECT_ROOT"), "/work/run_ourfm/fm")
regs <- read.delim(paste0(W, "/regions_by_size_v2.tsv"), header = FALSE, stringsAsFactors = FALSE)
rows <- list()
stat_phi <- function(p) { p <- p[is.finite(p)]; n <- length(p)
  c(n = n, ratio_max_min = max(p) / min(p), p99_p01 = if (n >= 10) as.numeric(quantile(p, .99) / quantile(p, .01)) else NA, cv = sd(p) / mean(p), max_share = max(p) / sum(p)) }
for (i in seq_len(nrow(regs))) { T <- regs[i, 1]; RID <- regs[i, 2]
  ph <- read.delim(sprintf("%s/phi/%s.phi.tsv", W, RID), stringsAsFactors = FALSE)
  key_full <- readLines(sprintf("%s/%s/%s.ld.vars", F0, T, RID)); M <- length(key_full)
  ar <- sprintf("%s/rw/%s.alpha.rds", W, RID)
  if (file.exists(ar)) { A <- readRDS(ar); key <- A$key; src <- "alpha"
    allcs <- susie_get_cs(list(alpha = A$alpha, V = A$V), coverage = 0.95)$cs
    pure_sets <- lapply(A$sets$cs, sort); pure_pur <- A$sets$purity$min.abs.corr
  } else { u <- read.delim(sprintf("%s/%s/%s.pip.tsv", F0, T, RID), stringsAsFactors = FALSE); key <- u$variant_hg19; src <- "pipfile"
    ids <- sort(unique(u$cs_id[!is.na(u$cs_id)])); allcs <- lapply(ids, function(j) which(u$cs_id == j)); pure_sets <- lapply(allcs, sort); pure_pur <- rep(NA_real_, length(allcs)) }
  if (length(allcs) == 0) next
  pi_ <- ph$phi[match(key, ph$variant_hg19)]; stopifnot(!any(is.na(pi_)))
  reg <- stat_phi(pi_)
  keep <- key_full %in% key; con <- file(sprintf("%s/%s/%s.ld.bin", F0, T, RID), "rb"); R <- readBin(con, "numeric", n = M * M, size = 4); close(con); dim(R) <- c(M, M)
  if (sum(keep) != M) R <- R[keep, keep]
  for (j in seq_along(allcs)) { cs <- sort(allcs[[j]]); k <- length(cs)
    pm <- if (length(pure_sets)) which(vapply(pure_sets, function(s) identical(s, cs), logical(1))) else integer(0); pure <- length(pm) > 0
    if (k > 1) { sub <- if (k > 100) { set.seed(1); sample(cs, 100) } else cs; mac <- min(abs(R[sub, sub][upper.tri(matrix(0, length(sub), length(sub)))])) } else mac <- 1
    sp <- stat_phi(pi_[cs]); pos <- as.integer(sapply(strsplit(key[cs], ":"), `[`, 2))
    rows[[length(rows) + 1]] <- data.frame(trait = T, region = RID, src = src, cs_idx = j, pure = pure, cs_size = k, size_bin = cut(k, c(0, 1, 4, 9, 29, Inf), labels = c("1", "2-4", "5-9", "10-29", "30+")),
      min_abs_corr = signif(mac, 4), min_abs_corr_saved = if (pure) signif(pure_pur[pm[1]], 4) else NA, span_bp = diff(range(pos)),
      cs_phi_ratio_max_min = signif(sp["ratio_max_min"], 4), cs_phi_p99_p01 = signif(sp["p99_p01"], 4), cs_phi_cv = signif(sp["cv"], 4), cs_phi_max_share = signif(sp["max_share"], 4),
      region_phi_ratio_max_min = signif(reg["ratio_max_min"], 4), region_phi_p99_p01 = signif(reg["p99_p01"], 4), region_phi_cv = signif(reg["cv"], 4), n_region = M, stringsAsFactors = FALSE) }
  rm(R); invisible(gc()); cat(sprintf("[%s] src=%s n_cs=%d pure=%d\n", RID, src, length(allcs), length(pure_sets)))
}
out <- do.call(rbind, rows); write.table(out, paste0(W, "/trackA_cs_phi_spread.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("CS_SPREAD_OK rows=%d regions=%d\n", nrow(out), length(unique(out$region))))
