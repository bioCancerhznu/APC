#===============================================================

#===============================================================

rm(list = ls())
library(dplyr)
library(MatchIt)
library(survival)
library(survminer)
library(officer)
library(flextable)
library(compareGroups)
library(broom)
library(cobalt)


#===============================================================

#===============================================================

jco_colors <- c('#FF7F5B', '#393939', '#4EA0AE', '#640D5F', '#FFF48F')


data <- read.csv("data.csv", header = T, row.names = 1)

data$Year <- NULL

table(data$AJCC_T_6th)
table(data$AJCC_N_6th)
table(data$AJCC_M_6th)

str(data)

table(data$Surgery)


data$Surgery[data$Surgery == "CPR"] <- 1
data$Surgery[data$Surgery == "nonCPR"] <- 0
data$Surgery <- as.numeric(as.character(data$Surgery))

table(data$Surgery)
data$Surgery[1:6]

data$Surgery <- factor(data$Surgery, levels = c(0, 1), labels = c("nonCPR", "CPR"))

table(data$Surgery)
data$Surgery[1:6]

str(data)

#=================================================================

#===============================================================


surv_obj <- Surv(data$OSTime, data$OS)

fit <- survfit(surv_obj ~ Surgery, data = data)

ggsurvplot(fit, data = data, pval = TRUE)

p <- ggsurvplot(
  fit,
  data = data,
  size = 1,
  conf.int = FALSE,
  pval = TRUE,
  xlim = c(0, 10),
  xlab = "Time in years",
  break.time.by = 2,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  ggtheme = theme_light(),
  palette = jco_colors  
)

print(p)


pdf(file = "orig_surg.pdf", height = 5, width = 5)
print(p)
dev.off()

fit



#===============================================================

#===============================================================

base_tab <- descrTable(Surgery ~ Age + Gender + Race + Marital + Laterality + PrimarySite + Type + 
                         T6th + N6th + M6th + Chemotherapy + Radiation,
                       data = data,  method = c("chi2"))

export2word(base_tab, file='surg_before.docx')

#===============================================================

#===============================================================

table(data$Surgery)
head(data)

str(data)

match_surgery <- matchit(Surgery ~ Age + Gender + Race + Marital + Laterality + PrimarySite + Type + 
                           T6th + N6th + M6th + Chemotherapy + Radiation, 
                         data = data, method = "nearest", distance = "logit", 
                         caliper = 0.2, ratio = 1)


summary(match_surgery)

matched_data <- match.data(match_surgery)

head(matched_data)

#===============================================================

#===============================================================

pdf("DiifMatching.pdf", height = 5, width = 7)

love.plot(match_surgery, stats = "mean.diffs",
          abs = TRUE, threshold = 0.1, 
          stars = "std")

dev.off()


#===============================================================

#===============================================================

surv_obj <- Surv(matched_data$OSTime, matched_data$OS)

fit <- survfit(surv_obj ~ Surgery, data = matched_data)

fit

p <- ggsurvplot(
  fit,
  data = matched_data,
  size = 1,
  conf.int = FALSE,
  pval = TRUE,
  xlim = c(0, 10),
  xlab = "Time in years",
  break.time.by = 2,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  ggtheme = theme_light(),
  palette = jco_colors  
)

print(p)

pdf(file = "spm_surg.pdf", height = 5, width = 5)
print(p)
dev.off()

fit

#===============================================================

#===============================================================

head(matched_data)

base_tab <- descrTable(Surgery ~ Age + Gender + Race + Marital + Laterality + PrimarySite + Type + 
                         T6th + N6th + M6th + Chemotherapy + Radiation,
                       data = matched_data,  method = c("chi2"))

export2word(base_tab, file='surg_after.docx')



#===============================================================

#===============================================================

table(matched_data$T6th, matched_data$N6th, matched_data$M6th)


table(matched_data$Age)

matched_data1 <- subset(matched_data, matched_data$Age == "66above")

surv_obj <- Surv(matched_data1$OSTime, matched_data1$OS)

fit <- survfit(surv_obj ~ Surgery, data = matched_data1)


p <- ggsurvplot(
  fit,
  data = matched_data1,
  size = 1,
  conf.int = FALSE,
  pval = TRUE,
  xlim = c(0, 10),
  xlab = "Time in years",
  break.time.by = 2,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  ggtheme = theme_light(),
  palette = jco_colors  
)

print(p)


pdf(file = "elder_surg.pdf", height = 5, width = 5)
print(p)
dev.off()


#===============================================================

#===============================================================


matched_data1 <- subset(matched_data, matched_data$Age == "56to65")

surv_obj <- Surv(matched_data1$OSTime, matched_data1$OS)

fit <- survfit(surv_obj ~ Surgery, data = matched_data1)

p <- ggsurvplot(
  fit,
  data = matched_data1,
  size = 1,
  conf.int = FALSE,
  pval = TRUE,
  xlim = c(0, 10),
  xlab = "Time in years",
  break.time.by = 2,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  ggtheme = theme_light(),
  palette = jco_colors  
)

print(p)


pdf(file = "median_surg.pdf", height = 5, width = 5)
print(p)
dev.off()



#===============================================================

#===============================================================

# 
matched_data1 <- subset(matched_data, matched_data$Age == "55to")

surv_obj <- Surv(matched_data1$OSTime, matched_data1$OS)

fit <- survfit(surv_obj ~ Surgery, data = matched_data1)


p <- ggsurvplot(
  fit,
  data = matched_data1,
  size = 1,
  conf.int = FALSE,
  pval = TRUE,
  xlim = c(0, 10),
  xlab = "Time in years",
  break.time.by = 2,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  ggtheme = theme_light(),
  palette = jco_colors  
)

print(p)


pdf(file = "young_surg.pdf", height = 5, width = 5)
print(p)
dev.off()


