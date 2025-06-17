install.packages("dplyr")
install.packages("survival")
install.packages("DiagrammeR")
library(DiagrammeR)
library(dplyr)
library(survival)
library(survminer)
NHANES <- readRDS("/Users/theo_hung/Library/CloudStorage/OneDrive-EmoryUniversity/Varghese, Jithin Sam's files - Theo Hung/working/rth01_selected df.RDS")
nrow(NHANES)
sum(is.na(NHANES$fasting_glucose))
sum(is.na(NHANES$glycohemoglobin))
#1. Prediabetes prediction 
NHANES <- NHANES %>%
  mutate(prediabetes = (dm_doc_told == 2 & is.na(dm_age)) &
           ((fasting_glucose >= 100 & fasting_glucose <= 125) |
              (glycohemoglobin >= 5.6 & glycohemoglobin <= 6.4)))  


NHANES %>%
  dplyr::filter(prediabetes == TRUE) %>%
  nrow()

prediabetics <- NHANES %>%
  dplyr::filter(prediabetes == TRUE)

sum(is.na(prediabetics$glycohemoglobin))
sum(is.na(prediabetics$fasting_glucose))

#2. Classifying prediabetic individuals into subclasses 

#inputting pooled mean and sd for scaling 
pooled_mean <- c(
  age = 56.06255, hba1c = 5.765866, bmi = 29.6381, 
  homa2b = 117.6387, homa2ir = 2.091171, egfr = 79.98367,
  sbp = 132.0694, dbp = 85.29373, ldlc = 126.8684, hdlc = 45.2273
)

pooled_sd <- c(
  age = 10.49422, hba1c = 0.444204, bmi = 6.187436, 
  homa2b = 63.87755, homa2ir = 1.473068, egfr = 21.0366,
  sbp = 22.10103, dbp = 15.94767, ldlc = 36.89973, hdlc = 19.10283
)

#renaming pooled mean and sd varibles to match prediabetics table 
names(pooled_mean)[names(pooled_mean) == "hba1c"] <- "glycohemoglobin"
names(pooled_mean)[names(pooled_mean) == "ldlc"] <- "ldl"
names(pooled_mean)[names(pooled_mean) == "hdlc"] <- "hdl"

names(pooled_sd)[names(pooled_sd) == "hba1c"] <- "glycohemoglobin"
names(pooled_sd)[names(pooled_sd) == "ldlc"] <- "ldl"
names(pooled_sd)[names(pooled_sd) == "hdlc"] <- "hdl"
vars <- names(pooled_mean)

#standardizing individual values 
prediabetics_scaled <- prediabetics %>%
  mutate(across(all_of(vars),
                ~ (. - pooled_mean[cur_column()]) / pooled_sd[cur_column()],
                .names = "scaled_{.col}"))



#standardizing the subtypes 
subtype_centroids <- data.frame(
  subtype = 1:7,
  age = c(59.64709118,59.78944708, 54.89094371, 44.08253801, 59.49204648, 65.64438561, 55.36587772),
  hba1c = c(5.503029412, 6.204511533, 5.63236755, 5.76603968, 5.753035033, 5.766283716, 5.821352478),
  bmi = c(24.98430676, 29.06433963, 26.27414524, 33.06441474, 33.17694144, 29.98786514, 32.18790268),
  homa2b = c(74.63511765, 92.78880597, 79.0343819, 120.8730482, 113.3552203, 252.0642857, 159.3466883),
  homa2ir = c(1.062294118, 1.627645862, 1.262444812, 2.032620459, 1.940062435, 4.853146853, 3.481333951),
  egfr = c(74.51230988, 72.67762125, 72.81202589, 107.5359492, 80.17968411, 71.51735816, 68.67514723),
  sbp = c(134.7091678, 126.4644328, 145.4261477, 115.1221316, 136.6411015, 123.0394516, 145.8778439),
  dbp = c(85.86461382, 78.06768554, 98.86771362, 76.87933071, 88.07331892, 70.14724069, 96.93581813),
  ldlc = c(113.5063959, 126.7215984, 147.7487464, 121.9864682, 122.3806313, 106.6826828, 146.5486526),
  hdlc = c(65.66297647, 38.74463704, 32.6306239, 47.26272481, 51.53480923, 51.8545989, 24.81419592)
)
subtype_centroids <- subtype_centroids %>%
  mutate(subtype = case_when(
    subtype == 1 ~ "overweight",
    subtype == 2 ~ "dysglycemic",
    subtype == 3 ~ "overweight hypertensive",
    subtype == 4 ~ "obese early onset",
    subtype == 5 ~ "obese moderate IR",
    subtype == 6 ~ "older severe IR",
    subtype == 7 ~ "obese severe IR",
    TRUE ~ NA_character_   
  ))

colnames(subtype_centroids)
scaled_centroids <- subtype_centroids

names(scaled_centroids)
vars %in% names(scaled_centroids)

scaled_centroids <- scaled_centroids %>%
  rename(
    glycohemoglobin = hba1c,
    ldl = ldlc,
    hdl = hdlc
  )
vars %in% names(scaled_centroids)

scaled_centroids[vars] <- sweep(scaled_centroids[vars], 2, pooled_mean[vars], FUN = "-")
scaled_centroids[vars] <- sweep(scaled_centroids[vars], 2, pooled_sd[vars], FUN = "/")

assign_subtype <- function(row) {
  distances <- apply(scaled_centroids[vars], 1, function(centroid) {
    sqrt(sum((row - centroid)^2, na.rm = TRUE))
  })
  which.min(distances)
}
prediabetics_scaled$subtype <- apply(
  prediabetics_scaled %>%
    select(starts_with("scaled_")), 1, assign_subtype
)

prediabetics_scaled <- prediabetics_scaled %>%
  mutate(subtype = case_when(
    subtype == 1 ~ "overweight",
    subtype == 2 ~ "dysglycemic",
    subtype == 3 ~ "overweight hypertensive",
    subtype == 4 ~ "obese early onset",
    subtype == 5 ~ "obese moderate IR",
    subtype == 6 ~ "older severe IR",
    subtype == 7 ~ "obese severe IR",
    TRUE ~ NA_character_   
  ))


#3.Study hazards of all-cause mortality using Cox Proportional Hazards regressions for each subtype, relative to Overweight Normotensive subtype



#a. no other covariaties 

sum(NHANES$mortstat == 1, na.rm = TRUE)

prediabetics_scaled$subtype <- relevel(factor(prediabetics_scaled$subtype), ref = "overweight")
surv_obj <- Surv(time = prediabetics_scaled$permth_exm, event = prediabetics_scaled$mortstat)
km_fit <- survfit(surv_obj ~ subtype, data = prediabetics_scaled)
ggsurvplot(
  km_fit,
  data = prediabetics_scaled,
  xlab = "Time (months)",
  ylab = "Survival Probability",
  surv.median.line = "hv",       # adds median survival line
  legend.title = "Subtype",
  pval = TRUE,                   # show log-rank p-value
  conf.int = TRUE                # show confidence intervals
)

#b.	Adjusted for age and sex
coxph(Surv(permth_exm, mortstat) ~ subtype + age + factor(gender), data = prediabetics_scaled)

