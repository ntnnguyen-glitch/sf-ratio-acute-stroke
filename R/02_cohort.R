# Build the analytic cohort: 24-h landmark, hourly SpO2/FiO2, exposure groups, covariates and outcomes.
source("R/00_setup.R")
con <- dbConnect(duckdb::duckdb(), file.path(WORK_DIR, "extract.duckdb"), read_only = TRUE)
Q <- function(sql) setDT(dbGetQuery(con, sql))

# S/F per SpO2 reading: each reading is paired with the latest FiO2 charted within the previous 4 h
# (21% when no oxygen record); readings are then summarised as the hourly median.
# max_flow: oxygen flow above this value (L/min) without a charted FiO2 makes the reading unusable.
hourly_sf <- function(sp, fi, max_flow = 15) {
  fi <- fi[(fio2 >= 21 & fio2 <= 100) | flow > max_flow][order(id, h)]
  fi[flow > max_flow, fio2 := -1]
  # one value per time point: a charted FiO2 takes precedence over a value derived from oxygen flow
  fi <- fi[, .(fio2 = if (any(flow == 0)) max(fio2[flow == 0]) else fio2[which.max(flow)]), by = .(id, h)]
  x <- fi[sp[order(id, h)], on = .(id, h), roll = 4]
  x[is.na(fio2), fio2 := 21]
  x <- x[fio2 > 0][, sf := spo2 / (fio2 / 100)]
  x[, .(sf = median(sf), spo2 = median(spo2), fio2 = median(fio2)), by = .(id, h = floor(h))]
}
add_windows <- function(co, x, pre = "") {
  w <- function(a, b) x[h >= a & h < b, .(m = min(sf), k = sum(sf <= 315),
                                          m97 = if (any(spo2 <= 97)) min(sf[spo2 <= 97]) else NA_real_), by = id]
  w1 <- w(0, 12); w2 <- w(12, 24)
  co[w1, on = .(id), (paste0(pre, c("sf_w1", "s97_w1"))) := .(i.m, i.m97)]
  co[w2, on = .(id), (paste0(pre, c("sf_w2", "s97_w2"))) := .(i.m, i.m97)]
  if (pre == "") {
    w0 <- w(0, 24); co[w0, on = .(id), `:=`(sf_min = i.m, hrs315 = i.k)]
    co[x[h >= 0 & h < 24, .(v = mean(sf)), by = id], on = .(id), sfmean := i.v]
  }
  co
}
icd <- function(dx, hadm, a, b) hadm %in% dx[(icd_version == 10 & grepl(paste0("^(", a, ")"), icd_code)) |
                                               (icd_version == 9 & grepl(paste0("^(", b, ")"), icd_code)), hadm_id]
flow <- list()

## MIMIC-IV ---------------------------------------------------------------
c <- Q("select * from mco")[age >= 18]; c[, `:=`(id = stay_id, T0 = intime)]; flow$m_start <- nrow(c)
ce <- Q("select * from mce")[c[, .(stay_id, T0)], on = .(stay_id), nomatch = 0]
ce[, `:=`(id = stay_id, h = as.numeric(difftime(charttime, T0, units = "hours")))]
sp <- ce[itemid == 220277 & valuenum >= 60 & valuenum <= 100, .(id, h, spo2 = valuenum)]
fi <- rbind(ce[itemid == 223835 & !is.na(valuenum), .(id, h, fio2 = ifelse(valuenum <= 1, valuenum * 100, valuenum), flow = 0)],
            ce[itemid == 223834 & valuenum >= 0 & valuenum <= 80, .(id, h, fio2 = 21 + 4 * valuenum, flow = valuenum)])
c <- add_windows(c, hourly_sf(sp, fi))
c <- add_windows(c, hourly_sf(sp, fi, max_flow = 6), "f6_")
g <- dcast(ce[itemid %in% c(220739, 223900, 223901) & h >= -6 & h <= 24], id + charttime ~ itemid, value.var = "valuenum", fun.aggregate = min)
g <- g[is.finite(`220739`) & is.finite(`223900`) & is.finite(`223901`), .(gcs = min(`220739` + `223900` + `223901`)), by = id]
c[g, on = .(id), gcs := i.gcs]
cv <- function(items, fn, lo, hi) { d <- ce[itemid %in% items & valuenum >= lo & valuenum <= hi & h >= -6 & h <= 24, .(v = fn(valuenum)), by = id]; d[c, on = .(id), v] }
c[, map_min := cv(c(220052, 220181), min, 20, 200)][, hr := cv(220045, max, 20, 250)][, rr := cv(c(220210, 224690), max, 4, 70)]
tc <- rbind(ce[itemid == 223762 & valuenum >= 30 & valuenum <= 45, .(id, h, t = valuenum)],
            ce[itemid == 223761 & valuenum >= 86 & valuenum <= 113, .(id, h, t = (valuenum - 32) / 1.8)])[h >= -6 & h <= 24, .(t = max(t)), by = id]
c[tc, on = .(id), temp := i.t]
cs <- ce[itemid == 223758]
c[, cmo24 := as.integer(id %in% cs[grepl("Comfort", value, ignore.case = TRUE) & h <= 24, id])]
c[, dnr := as.integer(id %in% cs[grepl("DNR|DNI|Comfort", value, ignore.case = TRUE), id])]
lb <- Q("select * from mlab")[c[, .(hadm_id, id, T0)], on = .(hadm_id), nomatch = 0]; lb[, h := as.numeric(difftime(charttime, T0, units = "hours"))]
lv <- function(it, fn, lo, hi) { d <- lb[itemid == it & valuenum >= lo & valuenum <= hi & h >= -6 & h <= 24, .(v = fn(valuenum)), by = id]; d[c, on = .(id), v] }
c[, creat := lv(50912, max, .1, 20)][, bun := lv(51006, max, 1, 250)][, wbc := lv(51301, max, .1, 200)][, hb := lv(51222, min, 3, 20)]
c[, plt := lv(51265, min, 1, 1500)][, glu := lv(50931, max, 20, 2000)][, na := lv(50983, max, 100, 180)]
c[, vaso := as.integer(stay_id %in% Q("select * from minp")$stay_id)]
# invasive ventilation: procedure 225792, endotracheal or tracheostomy tube as oxygen device, or an invasive ventilator mode
mv <- Q("select * from mvent")[c[, .(stay_id, T0)], on = .(stay_id), nomatch = 0][, h := as.numeric(difftime(charttime, T0, units = "hours"))]
non_inv_mode <- c("Standby", "SPONT", "CPAP", "Ambient")
pe <- Q("select * from mimv")[c[, .(stay_id, T0)], on = .(stay_id), nomatch = 0][, h := as.numeric(difftime(vs, T0, units = "hours"))]
iv <- rbind(pe[, .(stay_id, h)],
            mv[itemid == 226732 & value %in% c("Endotracheal tube", "Tracheostomy tube"), .(stay_id, h)],
            mv[itemid == 223849 & !is.na(value) & !value %in% non_inv_mode, .(stay_id, h)])
c[, imv_pre := as.integer(stay_id %in% iv[h <= 24, stay_id])]
c[, imv_pre_proc := as.integer(stay_id %in% pe[h <= 24, stay_id])]
c[iv[h > 24, .(s = min(h)), by = stay_id], on = .(stay_id), imv_h := i.s]
c[, niv24 := as.integer(stay_id %in% mv[h <= 24 & itemid == 226732 & value %in% c("Bipap mask", "CPAP mask", "High flow nasal cannula"), stay_id])]
dx <- Q("select * from mdx")
codes <- list(htn = c("I1[0-6]", "40[1-5]"), af = c("I48", "42731"), dm = c("E1[0-4]", "250"), chf = c("I50", "428"), ckd = c("N18", "585"),
              copd = c("J4[34]", "49[126]"), cancer = c("C[0-8]|C9[0-7]", "1[4-9][0-9]|20[0-8]"))
for (k in names(codes)) c[, (k) := as.integer(icd(dx, hadm_id, codes[[k]][1], codes[[k]][2]))]
c[, male := as.integer(gender == "M")]
c[, dt := fifelse(!is.na(deathtime), deathtime, as.POSIXct(as.Date(dod), tz = "UTC") + 43200)][, td := as.numeric(difftime(dt, T0, units = "hours"))]
c[, `:=`(icu_h = los * 24, hosp_h = as.numeric(difftime(dischtime, T0, units = "hours")), dhosp = hospital_expire_flag,
         db = "MIMIC-IV", hosp = 0, apache = NA_real_, o2doc = 1L)]
cm <- c

## eICU-CRD ---------------------------------------------------------------
c <- Q("select * from eco")[!is.na(age) & age >= 18]; c[, id := s]; flow$e_start <- nrow(c)
dx <- Q("select * from edx"); ph <- Q("select * from eph")
c[, sub := fifelse(grepl("Subarachnoid", adx), "SAH", fifelse(grepl("Hemorrhage/hematoma", adx), "ICH",
         fifelse(id %in% dx[grepl("subarachnoid", ds), s], "SAH",
         fifelse(id %in% dx[grepl("hemorrhagic stroke|intracerebral hemorrhage", ds) & !grepl("ischemic", ds), s], "ICH", "AIS"))))]
vp <- Q("select * from evp")[, `:=`(id = s, h = o / 60)]
rc <- Q("select * from erc")[, `:=`(id = s, h = o / 60, x = suppressWarnings(as.numeric(gsub("[^0-9.]", "", v))))]
el <- Q("select * from elab")[, `:=`(id = s, h = o / 60)]
fi <- rbind(rc[k %in% c("FiO2", "FIO2 (%)", "Set Fraction of Inspired Oxygen (FIO2)", "O2 Percentage") & !is.na(x), .(id, h, fio2 = ifelse(x <= 1, x * 100, x), flow = 0)],
            rc[k %in% c("LPM O2", "Oxygen Flow Rate") & x >= 0 & x <= 80, .(id, h, fio2 = 21 + 4 * x, flow = x)],
            el[ln == "FiO2" & !is.na(v), .(id, h, fio2 = ifelse(v <= 1, v * 100, v), flow = 0)])
c[, o2doc := as.integer(id %in% fi$id)]
spe <- vp[spo2 >= 60 & spo2 <= 100, .(id, h, spo2)]
c <- add_windows(c, hourly_sf(spe, fi))
c <- add_windows(c, hourly_sf(spe, fi, max_flow = 6), "f6_")
v24 <- vp[h >= -6 & h <= 24]
ap <- Q("select * from evap")[, `:=`(id = s, h = o / 60)][h >= -6 & h <= 24 & nimap >= 20 & nimap <= 200]
c[rbind(v24[mbp >= 20 & mbp <= 200, .(id, m = mbp)], ap[, .(id, m = nimap)])[, .(m = min(m)), by = id], on = .(id), map_min := i.m]
c[v24[hr >= 20 & hr <= 250, .(v = max(hr)), by = id], on = .(id), hr := i.v]
c[v24[rr >= 4 & rr <= 70, .(v = max(rr)), by = id], on = .(id), rr := i.v]
tp <- rbind(v24[tmp >= 30 & tmp <= 45, .(id, t = tmp)], Q("select s id, o, v t from etemp")[o >= -360 & o <= 1440 & t >= 30 & t <= 45, .(id, t)])
c[tp[, .(t = max(t)), by = id], on = .(id), temp := i.t]
aps <- Q("select * from eaps")
aps[, gcs := fifelse(meds == 1 | eyes < 1 | motor < 1 | verbal < 1, NA_real_, as.numeric(eyes + motor + verbal))]
c[aps[, .(gcs = gcs[1]), by = s], on = .(s), gcs := i.gcs]
ev <- function(n, fn, lo, hi) { d <- el[ln == n & v >= lo & v <= hi & h >= -6 & h <= 24, .(x = fn(v)), by = id]; d[c, on = .(id), x] }
c[, creat := ev("creatinine", max, .1, 20)][, bun := ev("BUN", max, 1, 250)][, wbc := ev("WBC x 1000", max, .1, 200)][, hb := ev("Hgb", min, 3, 20)]
c[, plt := ev("platelets x 1000", min, 1, 1500)][, glu := ev("glucose", max, 20, 2000)][, na := ev("sodium", max, 100, 180)]
inf <- Q("select * from einf")
c[, vaso := as.integer(id %in% inf[o >= -360 & o <= 1440 & grepl("norepi|levophed|epineph|dopamine|vasopressin|phenyleph", dn), s])]
dev <- Q("select * from edev")
dev_imv <- dev[(k == "O2 Device" & v %in% c("ventilator", "ett")) | (k == "Oxygen Delivery Method" & v == "mechanical ventilator") |
               grepl("mechanical ventilator mode|ventilator support mode", tolower(k)), .(s, o)]
iv <- rbind(Q("select s, vs o from erco"), Q("select s, o from etr"), dev_imv)[!is.na(o)]
c[, imv_pre := as.integer(s %in% iv[o <= 1440, s] | s %in% aps[intubated == 1, s])]
c[, imv_pre_proc := as.integer(s %in% rbind(Q("select s, vs o from erco"), Q("select s, o from etr"))[o <= 1440, s] | s %in% aps[intubated == 1, s])]
c[iv[o > 1440, .(o = min(o)), by = s], on = .(s), imv_h := i.o / 60]
niv_dev <- dev[o <= 1440 & ((k == "O2 Device" & grepl("bi-pap|bipap|cpap|high flow nasal", v)) | grepl("non-invasive ventilation mode|bipap", tolower(k))), s]
c[, niv24 := as.integer(s %in% c(niv_dev, Q("select s from eniv")$s))]
has <- function(rx) as.integer(c$id %in% unique(c(dx[grepl(rx, ds), s], ph[grepl(rx, p), s])))
rx <- list(htn = "hypertension", af = "atrial fibrillation", dm = "diabetes", chf = "congestive heart failure|chf", ckd = "renal insufficiency|chronic kidney",
           copd = "copd|chronic obstructive|emphysema", cancer = "cancer|malignan|carcinoma|lymphoma|leukemia")
for (k in names(rx)) c[, (k) := has(rx[[k]])]
cp <- Q("select * from ecp")
c[, cmo24 := as.integer(id %in% cp[grepl("Comfort", v) & o <= 1440, s])]
c[, dnr := as.integer(id %in% cp[grepl("No CPR|Do not resuscitate|No intubation|Comfort", v), s])]
c[, `:=`(male = as.integer(gender == "Male"), td = fifelse(hst == "Expired", ho / 60, NA_real_), icu_h = uo / 60, hosp_h = ho / 60,
         dhosp = as.integer(hst == "Expired"), db = "eICU-CRD")]
c[Q("select * from eapr"), on = .(s), apache := i.apache]
ce_ <- c
dbDisconnect(con, shutdown = TRUE)

## combine, landmark and outcomes -------------------------------------------
cols <- c("db", "id", "hosp", "age", "male", "sub", "htn", "af", "dm", "chf", "ckd", "copd", "cancer", "gcs", "map_min", "hr", "rr", "temp",
          "creat", "bun", "wbc", "hb", "plt", "glu", "na", "vaso", "apache", "cmo24", "dnr", "o2doc", "niv24", "imv_pre_proc",
          "sf_w1", "sf_w2", "s97_w1", "s97_w2", "f6_sf_w1", "f6_sf_w2", "sf_min", "sfmean", "hrs315", "imv_pre", "imv_h", "td", "icu_h", "hosp_h", "dhosp")
D <- rbind(cm[, ..cols], ce_[, ..cols])
cnt <- function(tag) flow[[tag]] <<- as.list(D[, .N, keyby = db][, setNames(N, db)])
D <- D[is.na(td) | td > 24]; cnt("alive24")
D <- D[icu_h >= 24]; cnt("icu24")
flow$imv_extra <- as.list(D[imv_pre == 1 & imv_pre_proc == 0, .N, keyby = db][, setNames(N, db)])
D <- D[imv_pre == 0]; cnt("no_imv24")
D <- D[is.finite(sf_w1) & is.finite(sf_w2)]; cnt("sf_available")
D <- D[!is.na(dhosp)]; cnt("cohort")   # hospital discharge status recorded
lv3 <- c("None", "Transient", "Persistent")
grp3 <- function(a, b, t = 315, tp = t) factor(fifelse(a <= tp & b <= tp, "Persistent", fifelse(a <= t | b <= t, "Transient", "None")), levels = lv3)
D[, grp := grp3(sf_w1, sf_w2)][, g235 := grp3(sf_w1, sf_w2, 315, 235)]
D[, g97 := grp3(s97_w1, s97_w2)][is.na(s97_w1) | is.na(s97_w2), g97 := NA]   # unclassifiable when no SpO2 <=97% in a window
D[, gf6 := grp3(f6_sf_w1, f6_sf_w2)][!is.finite(f6_sf_w1) | !is.finite(f6_sf_w2), gf6 := NA]
D[, L := td - 24]
# in-hospital death after the landmark (primary outcome); cumulative incidence with discharge alive as competing event
D[, ht := pmax(fifelse(dhosp == 1 & !is.na(L), L, hosp_h - 24), 0.5) / 24][, hs := fifelse(dhosp == 1, 1L, 2L)]
D[db == "MIMIC-IV", `:=`(d28 = as.integer(!is.na(L) & L <= 672), d90 = as.integer(!is.na(L) & L <= 2160), t90 = pmin(fifelse(is.na(L), 2160, L), 2160) / 24,
                         d365 = as.integer(!is.na(L) & L <= 8760), t365 = pmin(fifelse(is.na(L), 8760, L), 8760) / 24)]
# invasive ventilation within 7 days of the landmark; death before intubation is a competing event
D[, ih := imv_h - 24]
D[, fs := fifelse(!is.na(ih) & ih <= 168 & (is.na(L) | ih <= L), 1L, fifelse(!is.na(L) & L <= 168, 2L, 0L))]
D[, ft := pmax(fifelse(fs == 1, ih, fifelse(fs == 2, L, pmin(168, hosp_h - 24, na.rm = TRUE))), 0.5) / 24][, imv7 := as.integer(fs == 1)]
D[, icu_los := (icu_h - 24) / 24]
D[, sub := factor(sub, levels = c("AIS", "ICH", "SAH"))][, cl := fifelse(db == "eICU-CRD", paste0("h", hosp), paste0("m", id))]
flow$s97_unclassifiable <- as.list(D[is.na(g97), .N, keyby = db][, setNames(N, db)])
flow$f6_unclassifiable <- as.list(D[is.na(gf6), .N, keyby = db][, setNames(N, db)])
saveRDS(D, file.path(WORK_DIR, "analytic.rds")); write_json(flow, file.path(OUT_DIR, "flow.json"), auto_unbox = TRUE, pretty = TRUE)
print(D[, .(n = .N, hosp_death = sum(dhosp), imv7 = sum(imv7), niv = sum(niv24)), keyby = .(db, grp)]); print(flow)
