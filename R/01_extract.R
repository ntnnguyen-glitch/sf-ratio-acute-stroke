# Extract the data needed for the study from the raw MIMIC-IV v3.1 and eICU-CRD v2.0 files.
source("R/00_setup.R")
con <- dbConnect(duckdb::duckdb(), file.path(WORK_DIR, "extract.duckdb"))
X <- function(sql) { t0 <- Sys.time(); dbExecute(con, sql)
  rej <- tryCatch(dbGetQuery(con, "select count(*) n from reject_errors")$n, error = function(e) NA)
  message(round(difftime(Sys.time(), t0, units = "mins"), 1), " min (rejected rows: ", rej, "): ", substr(gsub("\\s+", " ", sql), 1, 60)) }
X("PRAGMA threads=8")
rm_ <- function(p) sprintf("read_csv('%s/%s.csv.gz', header=true, auto_detect=true, sample_size=200000, ignore_errors=true, store_rejects=true)", MIMIC_DIR, p)
re_ <- function(p) sprintf("read_csv('%s/%s.csv.gz', all_varchar=true, strict_mode=false, ignore_errors=true, store_rejects=true)", EICU_DIR, p)

## MIMIC-IV ---------------------------------------------------------------
X(paste0("create or replace table mdx_all as select hadm_id, seq_num, trim(icd_code) icd_code, icd_version from ", rm_("hosp/diagnoses_icd")))
X("create or replace table mst as with d as (select hadm_id,
  min(case when (icd_version=10 and icd_code like 'I60%') or (icd_version=9 and icd_code like '430%') then seq_num end) sah,
  min(case when (icd_version=10 and icd_code like 'I61%') or (icd_version=9 and icd_code like '431%') then seq_num end) ich,
  min(case when (icd_version=10 and icd_code like 'I63%') or (icd_version=9 and (regexp_matches(icd_code,'^43[34][0-9]1$') or icd_code='436')) then seq_num end) ais
  from mdx_all group by 1) select * from d where least(coalesce(sah,99),coalesce(ich,99),coalesce(ais,99))<=3")
X(paste0("create or replace table mco as with i as (select i.*, row_number() over(partition by i.subject_id order by i.intime) rn
  from ", rm_("icu/icustays"), " i where hadm_id in (select hadm_id from mst))
  select i.subject_id, i.hadm_id, i.stay_id, i.intime, i.outtime, i.los, p.gender, p.anchor_age+(year(i.intime)-p.anchor_year) age, p.dod,
  a.deathtime, a.dischtime, a.hospital_expire_flag,
  case when coalesce(d.sah,99)<=least(coalesce(d.ich,99),coalesce(d.ais,99)) then 'SAH' when coalesce(d.ich,99)<=coalesce(d.ais,99) then 'ICH' else 'AIS' end sub
  from i join ", rm_("hosp/patients"), " p using(subject_id) join ", rm_("hosp/admissions"), " a using(hadm_id) join mst d using(hadm_id) where rn=1"))
X(paste0("create or replace table mce as select c.stay_id, c.itemid, c.charttime, c.valuenum, c.value from ", rm_("icu/chartevents"),
  " c join mco m using(stay_id) where c.itemid in (220277,223835,223834,220739,223900,223901,220045,220210,224690,220052,220181,223761,223762,226512,223758)
  and c.charttime between m.intime - interval 6 hour and m.intime + interval 24 hour"))
X(paste0("create or replace table mlab as select l.hadm_id, l.itemid, l.charttime, l.valuenum from ", rm_("hosp/labevents"),
  " l join mco m using(hadm_id) where l.itemid in (50912,51006,51301,51222,51265,50931,50983,50882,50813)
  and l.charttime between m.intime - interval 6 hour and m.intime + interval 24 hour"))
X(paste0("create or replace table minp as select n.stay_id, n.starttime from ", rm_("icu/inputevents"),
  " n join mco m using(stay_id) where n.itemid in (221906,221289,221662,222315,221749)
  and n.starttime between m.intime - interval 6 hour and m.intime + interval 24 hour"))
X(paste0("create or replace table mimv as select p.stay_id, p.starttime vs, p.endtime ve from ", rm_("icu/procedureevents"),
  " p join mco using(stay_id) where p.itemid=225792"))
X("create or replace table mdx as select d.hadm_id, d.icd_code, d.icd_version from mdx_all d join mco using(hadm_id)")
X(paste0("create or replace table mvent as select c.stay_id, c.itemid, c.charttime, c.value from ", rm_("icu/chartevents"),
  " c join mco m using(stay_id) where c.itemid in (223849,226732) and c.charttime between m.intime - interval 6 hour and m.intime + interval 192 hour"))

## eICU-CRD ---------------------------------------------------------------
X(paste0("create or replace table ept as select * from ", re_("patient")))
X("create or replace table eco as with p as (select t.*, row_number() over(partition by uniquepid order by try_cast(hospitaldischargeyear as int),
  try_cast(patienthealthsystemstayid as int), try_cast(unitvisitnumber as int)) rn from ept t
  where regexp_matches(apacheadmissiondx,'^(CVA, cerebrovascular accident/stroke|Hemorrhage/hematoma, intracranial|Hemorrhage/hematoma-intracranial, surgery for|Subarachnoid hemorrhage/)'))
  select try_cast(patientunitstayid as int) s, try_cast(hospitalid as int) hosp, gender,
  case when age like '>%' then 90 else try_cast(age as int) end age, apacheadmissiondx adx, try_cast(admissionweight as double) w,
  try_cast(hospitaldischargeoffset as int) ho, hospitaldischargestatus hst, try_cast(unitdischargeoffset as int) uo from p where rn=1")
sid <- "try_cast(patientunitstayid as int)"
X(paste0("create or replace table evp as select ", sid, " s, try_cast(observationoffset as int) o, try_cast(sao2 as double) spo2, try_cast(heartrate as double) hr,
  try_cast(respiration as double) rr, try_cast(systemicmean as double) mbp, try_cast(temperature as double) tmp from ", re_("vitalPeriodic"),
  " where ", sid, " in (select s from eco) and try_cast(observationoffset as int) between -360 and 1440"))
X(paste0("create or replace table evap as select ", sid, " s, try_cast(observationoffset as int) o, try_cast(noninvasivemean as double) nimap from ", re_("vitalAperiodic"),
  " where ", sid, " in (select s from eco) and try_cast(observationoffset as int) between -360 and 1440"))
X(paste0("create or replace table erc as select ", sid, " s, try_cast(respchartoffset as int) o, respchartvaluelabel k, respchartvalue v from ", re_("respiratoryCharting"),
  " where ", sid, " in (select s from eco) and try_cast(respchartoffset as int) between -360 and 1440
  and respchartvaluelabel in ('FiO2','FIO2 (%)','Set Fraction of Inspired Oxygen (FIO2)','LPM O2','O2 Percentage','Oxygen Flow Rate')"))
X(paste0("create or replace table elab as select ", sid, " s, try_cast(labresultoffset as int) o, labname ln, try_cast(labresult as double) v from ", re_("lab"),
  " where ", sid, " in (select s from eco) and try_cast(labresultoffset as int) between -360 and 1440
  and labname in ('FiO2','creatinine','BUN','WBC x 1000','Hgb','platelets x 1000','glucose','sodium','bicarbonate','lactate')"))
X(paste0("create or replace table erco as select ", sid, " s, try_cast(ventstartoffset as int) vs from ", re_("respiratoryCare"),
  " where ", sid, " in (select s from eco) and try_cast(ventstartoffset as int)<>0"))
X(paste0("create or replace table etr as select ", sid, " s, try_cast(treatmentoffset as int) o from ", re_("treatment"),
  " where ", sid, " in (select s from eco) and regexp_matches(lower(treatmentstring),'mechanical ventilation|endotracheal|intubation')
  and not regexp_matches(lower(treatmentstring),'non-invasive|noninvasive')"))
X(paste0("create or replace table einf as select ", sid, " s, try_cast(infusionoffset as int) o, lower(drugname) dn from ", re_("infusionDrug"),
  " where ", sid, " in (select s from eco) and try_cast(infusionoffset as int) between -360 and 1440"))
X(paste0("create or replace table edx as select ", sid, " s, lower(diagnosisstring) ds from ", re_("diagnosis"), " where ", sid, " in (select s from eco)"))
X(paste0("create or replace table eph as select ", sid, " s, lower(pasthistorypath) p from ", re_("pastHistory"), " where ", sid, " in (select s from eco)"))
X(paste0("create or replace table eapr as select ", sid, " s, max(try_cast(apachescore as int)) apache from ", re_("apachePatientResult"),
  " where apacheversion='IVa' and ", sid, " in (select s from eco) group by 1"))
X(paste0("create or replace table eaps as select ", sid, " s, try_cast(intubated as int) intubated, try_cast(eyes as int) eyes, try_cast(motor as int) motor,
  try_cast(verbal as int) verbal, try_cast(meds as int) meds, try_cast(temperature as double) temperature, try_cast(sodium as double) sodium,
  try_cast(glucose as double) glucose from ", re_("apacheApsVar"), " where ", sid, " in (select s from eco)"))
X(paste0("create or replace table ecp as select ", sid, " s, try_cast(cplitemoffset as int) o, cplitemvalue v from ", re_("carePlanGeneral"),
  " where ", sid, " in (select s from eco) and cplgroup='Care Limitation'"))
X(paste0("create or replace table etemp as select ", sid, " s, try_cast(nursingchartoffset as int) o, try_cast(nursingchartvalue as double) v from ", re_("nurseCharting"),
  " where ", sid, " in (select s from eco) and nursingchartcelltypevallabel='Temperature' and nursingchartcelltypevalname='Temperature (C)' and try_cast(nursingchartoffset as int) between -360 and 1440"))
X(paste0("create or replace table edev as select ", sid, " s, try_cast(respchartoffset as int) o, respchartvaluelabel k, lower(respchartvalue) v from ", re_("respiratoryCharting"),
  " where ", sid, " in (select s from eco) and try_cast(respchartoffset as int) between -360 and 10080 and regexp_matches(lower(respchartvaluelabel), 'device|delivery method|ventilator mode|ventilation mode|bipap')"))
X(paste0("create or replace table eniv as select ", sid, " s, try_cast(treatmentoffset as int) o from ", re_("treatment"),
  " where ", sid, " in (select s from eco) and regexp_matches(lower(treatmentstring), 'non-invasive ventilation') and try_cast(treatmentoffset as int) between -360 and 1440"))
dbDisconnect(con, shutdown = TRUE)
