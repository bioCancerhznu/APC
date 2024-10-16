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

jco_colors <- c('#3bcf99', '#e6e261', '#da71ee', '#7c38cc', '#ebcea0')

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
  xlim = c(0, 5),
  xlab = "Time in years",
  break.time.by = 1,
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

base_tab <- descrTable(Surgery ~ Age + Gender + Race + Marital + Laterality + PrimarySite + Type + 
                         T6th + N6th + M6th + Chemotherapy + Radiation,
                       data = data,  method = c("chi2"))

export2word(base_tab, file='surg_before.docx')

#===============================================================

ggplot(data, aes(x = T6th, fill = Surgery)) +
  geom_bar(position = "fill") +
  labs(title = "Proportion of Surgery by AJCC T Stage",
       x = "AJCC T Stage", y = "Proportion") +
  theme_minimal()

ggplot(data, aes(x = N6th, fill = Surgery)) +
  geom_bar(position = "fill") +
  labs(title = "Proportion of Surgery by AJCC N Stage",
       x = "AJCC N Stage", y = "Proportion") +
  theme_minimal()

ggplot(data, aes(x = M6th, fill = Surgery)) +
  geom_bar(position = "fill") +
  labs(title = "Proportion of Surgery by AJCC M Stage",
       x = "AJCC M Stage", y = "Proportion") +
  theme_minimal()

#===============================================================

# T stage
t_stage_summary <- data %>%
  group_by(T6th, Surgery) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

print("Surgery Percentage by AJCC T Stage:")
print(t_stage_summary)

# N stage
n_stage_summary <- data %>%
  group_by(N6th, Surgery) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

print("Surgery Percentage by AJCC N Stage:")
print(n_stage_summary)

# M stage
m_stage_summary <- data %>%
  group_by(M6th, Surgery) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

print("Surgery Percentage by AJCC M Stage:")
print(m_stage_summary)

#===============================================================

table(data$Surgery)
head(data)

str(data)

# Calculate propensity scores and match
match_surgery <- matchit(Surgery ~ Age + Gender + Race + Marital + Laterality + PrimarySite + Type + 
                           T6th + N6th + M6th + Chemotherapy + Radiation, 
                         data = data, method = "nearest", distance = "logit", 
                         caliper = 0.2, ratio = 1)

summary(match_surgery)

# Extract matched data
matched_data <- match.data(match_surgery)

# Check matched data
head(matched_data)

#===============================================================

# Create a love plot to examine standardized mean differences
pdf("DiifMatching.pdf", height = 5, width = 7)

love.plot(match_surgery, stats = "mean.diffs",
          abs = TRUE, threshold = 0.1, 
          stars = "std")

dev.off()

#===============================================================

# Create survival object
surv_obj <- Surv(matched_data$OSTime, matched_data$OS)

fit <- survfit(surv_obj ~ Surgery, data = matched_data)

fit

# Plot survival curves
p <- ggsurvplot(
  fit,
  data = matched_data,
  size = 1,
  conf.int = FALSE,
  pval = TRUE,
  xlim = c(0, 5),
  xlab = "Time in years",
  break.time.by = 1,
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

# Export chi-square test results as a Word document
head(matched_data)

# Generate descriptive statistics and compute p-values
base_tab <- descrTable(Surgery ~ Age + Gender + Race + Marital + Laterality + PrimarySite + Type + 
                         T6th + N6th + M6th + Chemotherapy + Radiation,
                       data = matched_data,  method = c("chi2"))

export2word(base_tab, file='surg_after.docx')

#===============================================================

# View the number of people in different TNM combinations
table(matched_data$T6th, matched_data$N6th, matched_data$M6th)

table(matched_data$Age)

#===============================================================

# Subset data for ages 66 and above
matched_data1 <- subset(matched_data, matched_data$Age == "66above")

# Create survival object
surv_obj <- Surv(matched_data1$OSTime, matched_data1$OS)

# Survival curve based on surgery
fit <- survfit(surv_obj ~ Surgery, data = matched_data1)

# Plot survival curve
p <- ggsurvplot(
  fit,
  data = matched_data1,
  size = 1,
  conf.int = FALSE,
  pval = TRUE,
  xlim = c(0, 5),
  xlab = "Time in years",
  break.time.by = 1,
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

# Subset data for ages 56 to 65
matched_data1 <- subset(matched_data, matched_data$Age == "56to65")

# Create survival object
surv_obj <- Surv(matched_data1$OSTime, matched_data1$OS)

# Survival curve based on surgery
fit <- survfit(surv_obj ~ Surgery, data = matched_data1)

# Plot survival curve
p <- ggsurvplot(
  fit,
  data = matched_data1,
  size = 1,
  conf.int = FALSE,
  pval = TRUE,
  xlim = c(0, 5),
  xlab = "Time in years",
  break.time.by = 1,
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

# Subset data for ages 55 and below
matched_data1 <- subset(matched_data, matched_data$Age == "55to")

# Create survival object
surv_obj <- Surv(matched_data1$OSTime, matched_data1$OS)

# Survival curve based on surgery
fit <- survfit(surv_obj ~ Surgery, data = matched_data1)

# Plot survival curve
p <- ggsurvplot(
  fit,
  data = matched_data1,
  size = 1,
  conf.int = FALSE,
  pval = TRUE,
  xlim = c(0, 5),
  xlab = "Time in years",
  break.time.by = 1,
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
