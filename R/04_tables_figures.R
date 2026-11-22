# Baseline tables (csv) and Figures 1-4 (TIFF, 600 dpi).
source("R/00_setup.R"); library(grid)
D0 <- as.data.table(readRDS(file.path(WORK_DIR, "analytic.rds"))); flow <- read_json(file.path(OUT_DIR, "flow.json"))
M2 <- c("age", "male", "sub", "htn", "af", "dm", "chf", "ckd", "copd", "cancer")
M3x <- c("gcs", "map_min", "hr", "temp", "wbc", "hb", "plt", "creat", "bun", "glu", "na", "vaso")

## baseline tables ---------------------------------------------------------
iqr <- function(x, d = 1) sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), median(x, na.rm = TRUE), quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))
npc <- function(x) sprintf("%d (%.1f)", sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))
fp <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
table1 <- function(d) {
  g <- d$grp; rows <- list()
  num <- function(lab, v, dg = 1) { x <- d[[v]]; rows[[length(rows) + 1]] <<- c(lab, iqr(x, dg), tapply(x, g, iqr, d = dg), fp(kruskal.test(x ~ g)$p.value)) }
  bin <- function(lab, x, p = TRUE) rows[[length(rows) + 1]] <<- c(lab, npc(x), tapply(x, g, npc), if (p) fp(suppressWarnings(chisq.test(table(x, g))$p.value)) else "")
  rows[[1]] <- c("Number of patients", format(nrow(d), big.mark = ","), format(as.integer(table(g)), big.mark = ","), "")
  num("Age, years", "age"); bin("Male sex", d$male)
  rows[[length(rows) + 1]] <- c("Stroke subtype", "", "", "", "", fp(chisq.test(table(d$sub, g))$p.value))
  bin("   Ischaemic stroke", as.integer(d$sub == "AIS"), FALSE); bin("   Intracerebral haemorrhage", as.integer(d$sub == "ICH"), FALSE); bin("   Subarachnoid haemorrhage", as.integer(d$sub == "SAH"), FALSE)
  cm <- c(htn = "Hypertension", af = "Atrial fibrillation", dm = "Diabetes mellitus", chf = "Heart failure", ckd = "Chronic kidney disease", copd = "COPD", cancer = "Malignancy")
  for (v in names(cm)) bin(cm[[v]], d[[v]])
  num("Lowest GCS, 0–24 h", "gcs", 0); num("Lowest MAP, mmHg", "map_min", 0); num("Highest heart rate, beats/min", "hr", 0)
  num("Highest respiratory rate, breaths/min", "rr", 0); num("Highest temperature, °C", "temp"); num("White cell count, ×10⁹/L", "wbc")
  num("Haemoglobin, g/dL", "hb"); num("Platelets, ×10⁹/L", "plt", 0); num("Creatinine, mg/dL", "creat"); num("Blood urea nitrogen, mg/dL", "bun", 0)
  num("Glucose, mg/dL", "glu", 0); num("Sodium, mmol/L", "na", 0); bin("Vasopressor, 0–24 h", d$vaso)
  num("Lowest S/F, 0–24 h", "sf_min", 0); num("Hours with S/F ≤315, 0–24 h", "hrs315", 0)
  bin("Non-invasive ventilation or HFNC, 0–24 h", d$niv24); bin("In-hospital mortality", d$dhosp); bin("Invasive ventilation within 7 days", d$imv7); num("ICU stay after landmark, days", "icu_los")
  if (d$db[1] == "MIMIC-IV") { bin("28-day mortality", d$d28); bin("90-day mortality", d$d90); bin("365-day mortality", d$d365) }
  r <- as.data.table(do.call(rbind, rows))
  setnames(r, c("Characteristic", "Overall", "No S/F dysfunction", "Transient dysfunction", "Persistent dysfunction", "P value")); r
}
T1 <- list(); for (DBX in c("MIMIC-IV", "eICU-CRD")) { T1[[DBX]] <- table1(D0[db == DBX]); fwrite(T1[[DBX]], file.path(OUT_DIR, paste0("table1_", DBX, ".csv"))) }
miss <- D0[, lapply(.SD, function(x) sprintf("%d (%.1f)", sum(is.na(x)), 100 * mean(is.na(x)))), by = db, .SDcols = c(M3x, "rr")]
src <- D0[, .(o2doc = sum(o2doc), niv24 = sum(niv24)), by = db]
fwrite(miss, file.path(OUT_DIR, "missing_values.csv"))
saveRDS(list(T1 = T1, miss = miss, src = src), file.path(OUT_DIR, "tables.rds"))

## figures -----------------------------------------------------------------
LAB <- c(None = "No S/F dysfunction", Transient = "Transient dysfunction", Persistent = "Persistent dysfunction"); COL <- c("black", "#E69F00", "#D55E00")
DBS <- c("MIMIC-IV", "eICU-CRD")
th <- theme_classic(base_size = 10, base_family = "Arial") +
  theme(legend.position = "bottom", strip.background = element_blank(), strip.text = element_text(face = "bold", size = 10))
save_fig <- function(p, f, w = 7.5, h = 3.8) ggsave(file.path(OUT_DIR, f), p, width = w, height = h, dpi = 600, compression = "lzw", bg = "white", device = ragg::agg_tiff)

# Figure 1: flow chart; arrows stop short of the boxes
fmt <- function(x) vapply(x, function(z) format(z, big.mark = ","), character(1))
gn <- D0[, .N, keyby = .(db, grp)]
cnt <- function(DBX) { k <- sub("-", ".", DBX); s <- if (DBX == "MIMIC-IV") flow$m_start else flow$e_start
  a <- flow$alive24[[DBX]]; i <- flow$icu24[[DBX]]; v <- flow$no_imv24[[DBX]]; c0 <- flow$cohort[[DBX]]
  sa <- flow$sf_available[[DBX]]
  c(s, i, v, c0, s - a, a - i, i - v, v - sa, gn[db == DBX, N], sa - c0) }
draw_flow <- function() {
  grid.newpage(); pushViewport(viewport(xscale = c(0, 100), yscale = c(0, 100)))
  box <- function(x, y, w, h, t, fs = 7.2) { grid.rect(unit(x, "native"), unit(y, "native"), unit(w, "native"), unit(h, "native"), gp = gpar(fill = "white", lwd = .8))
    grid.text(t, unit(x, "native"), unit(y, "native"), gp = gpar(fontsize = fs, fontfamily = "Arial")) }
  arr <- function(x0, y0, x1, y1) grid.lines(unit(c(x0, x1), "native"), unit(c(y0, y1), "native"),
                                             arrow = arrow(length = unit(1.6, "mm"), type = "closed"), gp = gpar(fill = "black", lwd = .8))
  column <- function(cx, title, n) {
    ex <- cx + 22; ys <- c(80, 64, 48, 32)
    box(cx, 94, 26, 7, title, 7.8)
    txt <- c(sprintf("First ICU stay for acute stroke,\naged \u226518 years\n(n = %s)", fmt(n[1])), sprintf("Alive and in the ICU at 24 h\n(n = %s)", fmt(n[2])),
             sprintf("Not on invasive ventilation\nduring the first 24 h\n(n = %s)", fmt(n[3])), sprintf("Analytic cohort\n(n = %s)", fmt(n[4])))
    for (i in 1:4) box(cx, ys[i], 26, 9, txt[i])
    arr(cx, 90.5, cx, 85.4); for (i in 1:3) arr(cx, ys[i] - 4.5, cx, ys[i + 1] + 5.4)
    et <- c(sprintf("Died within 24 h (n = %s)\nICU stay <24 h (n = %s)", fmt(n[5]), fmt(n[6])), sprintf("Invasive ventilation\nwithin 24 h (n = %s)", fmt(n[7])),
            if (n[12] > 0) sprintf("S/F missing in either\n12-h window (n = %s)\nHospital discharge status\nnot recorded (n = %s)", fmt(n[8]), fmt(n[12])) else sprintf("S/F missing in either\n12-h window (n = %s)", fmt(n[8])))
    for (i in 1:3) { y <- c(72, 56, 40)[i]; box(ex, y, 17, if (i == 3 && n[12] > 0) 9.5 else 7.5, et[i], 6.6); arr(cx, y, ex - 9.4, y) }
    gx <- cx + c(-10, 0, 10); gt <- sprintf("%s\n(n = %s)", c("No S/F\ndysfunction", "Transient\ndysfunction", "Persistent\ndysfunction"), fmt(n[9:11]))
    for (i in 1:3) { box(gx[i], 13, 9.2, 10, gt[i], 6.6); arr(cx, 27.5, gx[i], 18.9) }
  }
  column(15, "MIMIC-IV v3.1 (single centre)", cnt("MIMIC-IV")); column(62, "eICU-CRD v2.0 (multicentre)", cnt("eICU-CRD"))
  grid.text("Groups defined by the lowest hourly S/F in the 0\u201312 h and 12\u201324 h windows (dysfunction = S/F \u2264315)",
            unit(50, "native"), unit(3, "native"), gp = gpar(fontsize = 6.8, fontfamily = "Arial", fontface = "italic"))
}
ragg::agg_tiff(file.path(OUT_DIR, "Figure1.tiff"), width = 7.5, height = 7, units = "in", res = 600, compression = "lzw"); draw_flow(); dev.off()

# Figure 2: cumulative incidence of in-hospital death (discharge alive as a competing event)
ch <- rbindlist(lapply(DBS, function(DBX) { d <- D0[db == DBX]; ci <- cuminc(d$ht, d$hs, d$grp)
  rbindlist(lapply(levels(d$grp), function(g) { z <- ci[[paste(g, 1)]]; data.table(t = z$time, y = z$est, g = g, db = DBX) })) }))
ch[, `:=`(db = factor(db, DBS), g = factor(LAB[g], LAB))]
lab <- data.table(db = factor(DBS, DBS), t = 1, y = .19, l = "Gray test P < 0.001")
p2 <- ggplot(ch[t <= 30], aes(t, y, colour = g)) + geom_step(linewidth = .7) + facet_wrap(~db) +
  geom_text(data = lab, aes(t, y, label = l), inherit.aes = FALSE, hjust = 0, size = 3, family = "Arial") +
  scale_y_continuous(labels = scales::percent, limits = c(0, .2)) + scale_x_continuous(breaks = seq(0, 30, 5)) +
  scale_colour_manual(values = COL, name = NULL) + labs(x = "Days after the 24-h landmark", y = "Cumulative incidence of in-hospital death") + th
save_fig(p2, "Figure2.tiff")

# Figure 3: adjusted spline of the lowest S/F
D <- copy(D0); for (v in M3x) D[, (v) := fifelse(is.na(get(v)), median(get(v), na.rm = TRUE), as.numeric(get(v))), by = db]
cv <- c(M2, M3x)
rc <- rbindlist(lapply(DBS, function(DBX) { d <- D[db == DBX]; fit <- glm(as.formula(paste("dhosp ~ ns(sf_min, df = 4) +", paste(cv, collapse = "+"))), binomial, d)
  nd <- d[rep(1, 141)]; for (v in cv) if (is.numeric(d[[v]])) nd[[v]] <- median(d[[v]]); nd$sub <- factor("AIS", levels(d$sub))
  nd$sf_min <- seq(quantile(d$sf_min, .01), 470, length.out = 141); p <- predict(fit, nd, se.fit = TRUE); r0 <- predict(fit, transform(nd[1], sf_min = 450))
  data.table(db = DBX, x = nd$sf_min, or = exp(p$fit - r0), lo = exp(p$fit - 1.96 * p$se.fit - r0), hi = exp(p$fit + 1.96 * p$se.fit - r0)) }))
rc[, db := factor(db, DBS)]
p3 <- ggplot(rc, aes(x, or)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80") + geom_line(linewidth = .8) +
  geom_hline(yintercept = 1, linetype = 2) + geom_vline(xintercept = 315, linetype = 3) + facet_wrap(~db) + scale_y_log10() + scale_x_reverse() +
  labs(x = "Lowest SpO\u2082/FiO\u2082 ratio, 0\u201324 h", y = "Adjusted odds ratio for in-hospital death\n(reference S/F = 450)") + th
save_fig(p3, "Figure3.tiff")

# Figure 4: cumulative incidence of invasive ventilation (death as competing event)
cf <- rbindlist(lapply(DBS, function(DBX) { d <- D0[db == DBX]; ci <- cuminc(d$ft, d$fs, d$grp)
  rbindlist(lapply(levels(d$grp), function(g) { z <- ci[[paste(g, 1)]]; data.table(t = z$time, y = z$est, g = g, db = DBX) })) }))
cf[, `:=`(db = factor(db, DBS), g = factor(LAB[g], LAB))]
lab2 <- data.table(db = factor(DBS, DBS), t = 1, y = .125, l = "Gray test P < 0.001")
p4 <- ggplot(cf, aes(t, y, colour = g)) + geom_step(linewidth = .7) + facet_wrap(~db) +
  geom_text(data = lab2, aes(t, y, label = l), inherit.aes = FALSE, hjust = 0, size = 3, family = "Arial") +
  scale_y_continuous(labels = scales::percent, limits = c(0, .13)) + scale_x_continuous(breaks = 0:7) +
  scale_colour_manual(values = COL, name = NULL) + labs(x = "Days after the 24-h landmark", y = "Cumulative incidence of invasive ventilation") + th
save_fig(p4, "Figure4.tiff")
