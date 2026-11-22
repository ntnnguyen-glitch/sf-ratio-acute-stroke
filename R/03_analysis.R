# Statistical analysis. Missing covariates are handled by multiple imputation (20 datasets, chained equations)
# within each database; estimates are combined with Rubin's rules. All contrasts are against no S/F dysfunction.
source("R/00_setup.R"); library(WeightIt); set.seed(2026)
D0 <- as.data.table(readRDS(file.path(WORK_DIR, "analytic.rds")))
M2 <- c("age", "male", "sub", "htn", "af", "dm", "chf", "ckd", "copd", "cancer")
M3x <- c("gcs", "map_min", "hr", "temp", "wbc", "hb", "plt", "creat", "bun", "glu", "na", "vaso")
G2 <- c("grpTransient", "grpPersistent"); M <- 20
D0[, `:=`(sfd = -sf_min / 50, h6 = hrs315 / 6, sfmd = -sfmean / 50, gnum = as.integer(grp) - 1)]

## multiple imputation within each database ---------------------------------
impute <- function(d) {
  v <- c(M2, M3x, "rr", if (d$db[1] == "eICU-CRD") "apache", "grp", "sf_min", "dhosp", "imv7")
  imp <- mice(d[, ..v], m = M, maxit = 10, printFlag = FALSE, seed = 2026)
  lapply(seq_len(M), function(i) { x <- copy(d); ci <- complete(imp, i); for (k in c(M3x, "rr", if (d$db[1] == "eICU-CRD") "apache")) x[[k]] <- ci[[k]]; x })
}
IMP <- list(`MIMIC-IV` = impute(D0[db == "MIMIC-IV"]), `eICU-CRD` = impute(D0[db == "eICU-CRD"]))
IMP$Pooled <- lapply(seq_len(M), function(i) rbind(IMP$`MIMIC-IV`[[i]], IMP$`eICU-CRD`[[i]]))
med <- copy(D0); for (v in c(M3x, "rr")) med[, (v) := fifelse(is.na(get(v)), median(get(v), na.rm = TRUE), as.numeric(get(v))), by = db]

## helpers -------------------------------------------------------------------
covs <- function(pooled, extra = NULL, drop = NULL) setdiff(c(M2, M3x, extra, if (pooled) "db"), drop)
fml <- function(y, x, cv) as.formula(paste(y, "~", paste(c(x, cv), collapse = "+")))
vc <- function(m, d) if (uniqueN(d$cl) < nrow(d)) vcovCL(m, cluster = d$cl) else vcov(m)
# Rubin's rules for one or more coefficients; returns estimate, 95% CI and P on the exponentiated scale
rubin <- function(fits, k) {
  Q <- sapply(fits, function(f) f$b[k]); U <- sapply(fits, function(f) f$v[k]); m <- length(fits)
  Q <- matrix(Q, nrow = length(k)); U <- matrix(U, nrow = length(k))
  qb <- rowMeans(Q); ub <- rowMeans(U); bb <- apply(Q, 1, var); tt <- ub + (1 + 1 / m) * bb
  df <- (m - 1) * (1 + ub / ((1 + 1 / m) * bb))^2; se <- sqrt(tt)
  cbind(exp(qb), exp(qb - qt(.975, df) * se), exp(qb + qt(.975, df) * se), 2 * pt(-abs(qb / se), df))
}
# multivariate Wald test combined across imputations (D1 statistic)
wald_mi <- function(fits, k) {
  m <- length(fits); Q <- sapply(fits, function(f) f$b[k]); Qb <- rowMeans(Q)
  Ub <- Reduce(`+`, lapply(fits, function(f) f$V[k, k])) / m; B <- cov(t(Q)); r <- (1 + 1 / m) * sum(diag(B %*% solve(Ub))) / length(k)
  W <- as.numeric(t(Qb) %*% solve(Ub) %*% Qb) / (length(k) * (1 + r)); pf(W, length(k), 1e6, lower.tail = FALSE)
}
gfit <- function(d, y, x, cv) { m <- glm(fml(y, x, cv), binomial, d); V <- vc(m, d); list(b = coef(m), v = diag(V), V = V) }
single <- function(d, y, x, cv, k) { f <- gfit(d, y, x, cv); se <- sqrt(f$v[k]); b <- f$b[k]; cbind(exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se), 2 * pnorm(-abs(b / se))) }
mi <- function(L, FUN, k) rubin(lapply(L, FUN), k)
coxf <- function(d, tt, ev, cv) { m <- coxph(as.formula(paste0("Surv(", tt, ",", ev, ") ~ ", paste(c("grp", cv, "cluster(cl)"), collapse = "+"))), d)
  list(b = coef(m), v = diag(vcov(m))) }
fgf <- function(d, cv) { mm <- model.matrix(as.formula(paste("~", paste(c("grp", cv), collapse = "+"))), d)[, -1]; mm <- mm[, apply(mm, 2, var) > 0]
  f <- crr(d$ft, d$fs, mm, failcode = 1, cencode = 0); v <- diag(f$var); names(v) <- names(f$coef); list(b = f$coef, v = v) }
owf <- function(d, cv) { w <- weightit(as.formula(paste("grp ~", paste(cv, collapse = "+"))), d, method = "glm", estimand = "ATO")
  m <- glm(dhosp ~ grp, quasibinomial, d, weights = w$weights); V <- vcovCL(m, cluster = d$cl, type = "HC0")
  list(b = coef(m), v = diag(V), smd = max(abs(cobalt::bal.tab(w, stats = "m")$Balance.Across.Pairs$Max.Diff.Adj), na.rm = TRUE)) }

run <- function(key, pooled) {
  L <- IMP[[key]]; d0 <- if (pooled) D0 else D0[db == key]; dm <- if (pooled) med else med[db == key]
  cv <- covs(pooled); o <- list(n = nrow(d0))
  o$crude <- mi(L, function(d) gfit(d, "dhosp", "grp", if (pooled) "db"), G2)
  o$m2 <- mi(L, function(d) gfit(d, "dhosp", "grp", c(M2, if (pooled) "db")), G2)
  o$m3 <- mi(L, function(d) gfit(d, "dhosp", "grp", cv), G2)
  o$trend <- mi(L, function(d) gfit(d, "dhosp", "gnum", cv), "gnum")
  o$sfd <- mi(L, function(d) gfit(d, "dhosp", "sfd", cv), "sfd")
  o$h6 <- mi(L, function(d) gfit(d, "dhosp", "h6", cv), "h6")
  sp <- lapply(L, function(d) gfit(d, "dhosp", "ns(sf_min, df = 4)", cv)); o$spline_p <- wald_mi(sp, grep("^ns\\(", names(sp[[1]]$b), value = TRUE))
  o$group_p <- wald_mi(lapply(L, function(d) gfit(d, "dhosp", "grp", cv)), G2)
  # secondary outcomes
  o$fg <- mi(L, function(d) fgf(d, cv), G2)
  o$gray_p <- cuminc(d0$ft, d0$fs, d0$grp)$Tests["1", "pv"]
  o$cif_p <- cuminc(d0$ht, d0$hs, d0$grp)$Tests["1", "pv"]
  o$los <- mi(L, function(d) { dl <- d[icu_los > 0]; m <- lm(fml("log(icu_los)", "grp", cv), dl); list(b = coef(m), v = diag(vc(m, dl))) }, G2)
  if (key == "MIMIC-IV") {
    o$d28 <- mi(L, function(d) gfit(d, "d28", "grp", cv), G2)
    o$d90 <- mi(L, function(d) coxf(d, "t90", "d90", cv), G2); o$d365 <- mi(L, function(d) coxf(d, "t365", "d365", cv), G2) }
  # confounding control
  ow <- lapply(L, function(d) owf(d, cv)); o$ow <- rubin(ow, G2); o$ow_smd <- mean(sapply(ow, `[[`, "smd"))
  e <- o$m3[2, ]; o$evalue <- evalues.OR(e[1], e[2], e[3], rare = mean(d0$dhosp) < 0.15)
  # sensitivity analyses (persistent vs none unless stated)
  S <- list()
  sub_mi <- function(nm, filt, x = "grp", cc = cv, k = G2) {
    S[[nm]] <<- list(r = mi(L, function(d) gfit(d[eval(filt, d)], "dhosp", x, cc), k), n = nrow(L[[1]][eval(filt, L[[1]])])) }
  S$complete_case <- list(r = single(na.omit(d0, cols = c(M2, M3x)), "dhosp", "grp", cv, G2), n = nrow(na.omit(d0, cols = c(M2, M3x))))
  S$median_imputation <- list(r = single(dm, "dhosp", "grp", cv, G2), n = nrow(dm))
  sub_mi("spo2_le97", quote(!is.na(g97)), "g97", k = c("g97Transient", "g97Persistent"))
  sub_mi("flow_le6", quote(!is.na(gf6)), "gf6", k = c("gf6Transient", "gf6Persistent"))
  sub_mi("threshold_235", quote(rep(TRUE, length(id))), "g235", k = c("g235Transient", "g235Persistent"))
  sub_mi("excl_niv_hfnc", quote(niv24 == 0)); sub_mi("excl_comfort", quote(cmo24 == 0)); sub_mi("excl_dnr", quote(dnr == 0))
  sub_mi("excl_early_death", quote(dhosp == 0 | is.na(L) | L > 48)); sub_mi("excl_copd", quote(copd == 0), cc = setdiff(cv, "copd"))
  sub_mi("adj_rr", quote(rep(TRUE, length(id))), cc = c(cv, "rr")); sub_mi("mean_sf", quote(rep(TRUE, length(id))), "sfmd", k = "sfmd")
  if (key == "eICU-CRD") sub_mi("adj_apache", quote(rep(TRUE, length(id))), cc = c(cv, "apache"))
  if (key != "MIMIC-IV") sub_mi("o2_recorded", quote(o2doc == 1))
  o$S <- S
  # exploratory effect modification
  if (pooled) for (v in c("copd", "sub")) {
    cc <- setdiff(cv, v); it <- lapply(L, function(d) gfit(d, "dhosp", paste0("grp*", v), cc))
    lv <- if (v == "copd") c("1", "0") else levels(D0$sub)
    o$em[[v]] <- list(p_int = wald_mi(it, grep("^grp.*:", names(it[[1]]$b), value = TRUE)),
      levels = lapply(setNames(lv, lv), function(l) list(n = sum(as.character(d0[[v]]) == l),
        r = mi(L, function(d) gfit(d[as.character(get(v)) == l], "dhosp", "grp", cc), G2)))) }
  o
}
R <- list(); for (k in c("MIMIC-IV", "eICU-CRD")) R[[k]] <- run(k, FALSE); R$Pooled <- run("Pooled", TRUE)
R$counts <- D0[, .(n = .N, dhosp = sum(dhosp), imv = sum(imv7), niv = sum(niv24)), keyby = .(db, grp)]
saveRDS(R, file.path(OUT_DIR, "results.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "session_info.txt"))
