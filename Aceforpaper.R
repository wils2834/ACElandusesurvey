# DRACE 2024  ACE (Autoclaved Citrate Extractable) Protein
library(readxl)
library(dplyr)
library(lme4)
library(lmerTest)   
library(emmeans)
library(multcomp) 
library(ggplot2)
# PART 1 Import and clean
drace24 <- read_excel("C:\\Users\\swils\\Downloads\\01 Projects\\2024 CRP\\2024DRACE_03042025\\AllplatesDRACE2024.xlsx", sheet = "raw_data")
colnames(drace24) <- make.names(colnames(drace24), unique = TRUE)

# Drop plate-calibration standards- not field samples.
drace24 <- drace24 %>%
  filter(!Type %in% c("Standard", "Standard 2", "Standard 3", "Standard 4"))

drace24 <- drace24 %>% filter(!is.na(Avg_runs))

drace24 <- drace24 %>%
  mutate(Type = recode(Type, "SH" = "Soil_Health", "RG" = "Rotational_Grazing"))

drace24$Type <- factor(
  drace24$Type,
  levels = c("Conventional", "Soil_Health", "Rotational_Grazing", "New_CRP", "Old_CRP")
)
# True field ID (see header note) use this for all random effects.
drace24 <- drace24 %>%
  mutate(Field = paste(Farmer, Type, sep = "_"))

# Sanity checks
cat("Rows after cleaning:", nrow(drace24), "\n")
cat("Distinct fields:", n_distinct(drace24$Field), "(expect 31)\n")
cat("Distinct sample locations (UID):", n_distinct(drace24$UID), "\n\n")
stopifnot(n_distinct(drace24$Field) == 31)

# PART 2  Mixed model
#Fit raw- and log-scale models, keep whichever gives better-behaved
# (more normal) residuals rather than assuming log is always right.
model_ACE     <- lmer(Avg_runs ~ Type + (1 | Field), data = drace24)
model_ACE_log <- lmer(log(Avg_runs) ~ Type + (1 | Field), data = drace24)

shapiro_raw <- shapiro.test(residuals(model_ACE))
shapiro_log <- shapiro.test(residuals(model_ACE_log))
cat("Shapiro-Wilk, raw-scale residuals:  p =", signif(shapiro_raw$p.value, 3), "\n")
cat("Shapiro-Wilk, log-scale residuals:  p =", signif(shapiro_log$p.value, 3), "\n")

use_log <- shapiro_log$p.value > shapiro_raw$p.value
final_model <- if (use_log) model_ACE_log else model_ACE
cat("Using", if (use_log) "LOG-transformed" else "RAW-scale", "model for inference.\n\n")

summary(final_model)

# Estimated marginal means (back-transformed if log model) + Tukey CLD
emm_ACE       <- emmeans(final_model, "Type")
emm_ACE_table <- summary(emm_ACE, type = "response")
cld_ACE       <- cld(emm_ACE, Letters = letters, adjust = "tukey")

print(emm_ACE_table)
print(cld_ACE)

# Tidy summary table used for plotting
mean_col <- if ("response" %in% names(emm_ACE_table)) "response" else "emmean"

ACE_summary <- as.data.frame(emm_ACE_table) %>%
  dplyr::select(Type, Mean = !!mean_col) %>%
  dplyr::left_join(
    as.data.frame(cld_ACE) %>% dplyr::select(Type, CLD = .group),
    by = "Type"
  )

print(ACE_summary)

# PART 3 Plots

type_labels <- c(
  "Conventional"        = "Conventional Cropping & Tillage",
  "Soil_Health"         = "No Till and/or Cover Crops",
  "Rotational_Grazing"  = "Rotational Grazing",
  "New_CRP"             = "0-10 years CRP",
  "Old_CRP"             = "10+ years CRP"
)

fill_colors <- c(
  "Conventional"        = "#E69F00",
  "Soil_Health"         = "#56B4E9",
  "Rotational_Grazing"  = "#009E73",
  "New_CRP"             = "#F0E442",
  "Old_CRP"             = "#0072B2"
)

#  Overall boxplot
label_positions <- drace24 %>%
  group_by(Type) %>%
  summarise(y_pos = max(Avg_runs, na.rm = TRUE) * 1.05, .groups = "drop") %>%
  dplyr::left_join(ACE_summary %>% dplyr::select(Type, CLD), by = "Type")

p_overall <- ggplot(drace24, aes(x = Type, y = Avg_runs, fill = Type)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  geom_text(
    data = label_positions,
    aes(x = Type, y = y_pos, label = CLD),
    inherit.aes = FALSE, size = 6, fontface = "bold", color = "black"
  ) +
  scale_x_discrete(labels = type_labels) +
  scale_fill_manual(values = fill_colors, guide = "none") +
  labs(
    title = "Autoclaved Citrate Extractable Protein by Management",
    x = NULL,
    y = expression(paste("g Protein kg"^-1, " soil"))
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 13, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 13, face = "bold", color = "black"),
    axis.title.y = element_text(size = 21, face = "bold", color = "black"),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

print(p_overall)
ggsave("ACE_overall_boxplot.png", p_overall, width = 9, height = 6, dpi = 300)

#  Faceted by county 
label_positions_county <- drace24 %>%
  group_by(County, Type) %>%
  summarise(y_pos = max(Avg_runs, na.rm = TRUE) * 1.05, .groups = "drop") %>%
  dplyr::left_join(ACE_summary %>% dplyr::select(Type, CLD), by = "Type")

p_county <- ggplot(drace24, aes(x = Type, y = Avg_runs, fill = Type)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  geom_text(
    data = label_positions_county,
    aes(x = Type, y = y_pos, label = CLD),
    inherit.aes = FALSE, size = 5, fontface = "bold", color = "black"
  ) +
  facet_wrap(~ County) +
  scale_fill_manual(values = fill_colors, labels = type_labels, name = "Management") +
  labs(
    title = "Autoclaved Citrate Extractable Protein by Management and County",
    x = NULL,
    y = expression(paste("g Protein kg"^-1, " soil"))
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 17, face = "bold"),
    strip.text   = element_text(size = 17, face = "bold"),
    plot.title   = element_text(size = 17, face = "bold", hjust = 0.5),
    legend.text  = element_text(size = 12, face = "bold", color = "black"),
    legend.title = element_text(size = 13, face = "bold", color = "black")
  )

print(p_county)
ggsave("ACE_by_county_boxplot.png", p_county, width = 11, height = 6, dpi = 300)

cat("\nSaved ACE_overall_boxplot.png and ACE_by_county_boxplot.png\n")

# ---- Metadata export ----------------------------------------------------
pairs(emm_ACE, adjust = "sidak")
anova_ACE <- anova(final_model)

ACE_meta <- as.data.frame(emm_ACE_table) %>%
  dplyr::rename(Mean = if ("response" %in% names(emm_ACE_table)) "response" else "emmean") %>%
  dplyr::left_join(
    as.data.frame(cld_ACE) %>% dplyr::select(Type, CLD = .group),
    by = "Type"
  ) %>%
  mutate(
    Metric         = "ACE",
    Scale          = if (use_log) "log" else "raw",
    Shapiro_raw_p  = signif(shapiro_raw$p.value, 3),
    Shapiro_log_p  = signif(shapiro_log$p.value, 3),
    LMM_F          = round(anova_ACE$`F value`[1], 3),
    LMM_p          = signif(anova_ACE$`Pr(>F)`[1], 3)
  ) %>%
  dplyr::select(Metric, Scale, Shapiro_raw_p, Shapiro_log_p,
                LMM_F, LMM_p, Type, Mean, SE, CLD)

write.csv(ACE_meta, "ACE_metadata.csv", row.names = FALSE)
cat("Saved ACE_metadata.csv\n")

#  PART 4  Indexed (continuous management score) version

library(MuMIn)

# Load pre-computed scores and correct stale RG/CRP values

scores_full <- read.csv(
  "C:\\Users\\swils\\Downloads\\01 Projects\\2024 CRP\\Landusesurveygithub\\managementoflandusesurvey2024\\soil_health_scores_2024_all_types.csv",
  stringsAsFactors = FALSE
)

rg_crp_farmer_map <- c(
  "5G_Farms"           = "5G_Farms_LLC",
  "Cliff"              = "Cliff",
  "Cliff_2"            = "Cliff_2",
  "Collins"            = "Steve_Collins",
  "Derouin"            = "Michael_Derouin",
  "Eiken"              = "Kay_Eiken",
  "Eiken_2"            = "Kay_Eiken_2",
  "Haugen"             = "Bonnie_Haugen",
  "Lee"                = "Carter_Lee",
  "Lee_2"              = "Carter_Lee_2",
  "McLean"             = "Maurice_McLean",
  "Menne"              = "Kevin_Menne",
  "Mitchell"           = "Thomas_Mitchell",
  "Reichling_Brothers" = "Reichling_Brothers",
  "RLH_Inc"            = "RLH_Inc",
  "RLH_Inc_2"          = "RLH_Inc_2",
  "Schmidt"            = "Dennis_Schmidt",
  "UW_AG"              = "UW_AG",
  "UW_AG_2"            = "UW_AG_2"
)

# CRP contract-age split, from the CRP Fields sheet's "treatment" column
crp_type_map <- c(
  "Schmidt"             = "Old_CRP",
  "Reichling_Brothers"  = "Old_CRP",
  "McLean"              = "Old_CRP",
  "Eiken"               = "Old_CRP",
  "Eiken_2"             = "Old_CRP",
  "Derouin"             = "Old_CRP",
  "Collins"             = "New_CRP",
  "Menne"               = "New_CRP",
  "5G_Farms"            = "New_CRP",
  "Mitchell"            = "New_CRP",
  "RLH_Inc"             = "New_CRP",
  "RLH_Inc_2"           = "New_CRP"
)

scores_full <- scores_full %>%
  mutate(
    Farmer = ifelse(farm_name %in% names(rg_crp_farmer_map),
                    rg_crp_farmer_map[farm_name],
                    Farmer),
    Type = ifelse(Type == "CRP" & farm_name %in% names(crp_type_map),
                  crp_type_map[farm_name],
                  Type)
  )
scores_full <- scores_full %>%
  mutate(Type = recode(Type,
                       "CRP_New" = "New_CRP",
                       "CRP_Old" = "Old_CRP"))
cat("Corrected Farmer/Type for RG/CRP rows:\n")
print(scores_full %>% filter(Type %in% c("Rotational_Grazing","Old_CRP","New_CRP")) %>%
        dplyr::select(farm_name, Farmer, Type))

# Join management score onto the ACE data

drace24_indexed <- drace24 %>%
  mutate(Farmer = as.character(Farmer), Type = as.character(Type)) %>%
  left_join(
    scores_full %>% dplyr::select(Farmer, Type, soil_health_score),
    by = c("Farmer", "Type")
  )

cat("\nRows with missing soil_health_score after join:",
    sum(is.na(drace24_indexed$soil_health_score)), "\n")

# Mixed model with continuous management score

model_ACE_idx     <- lmer(Avg_runs      ~ soil_health_score + (1 | Field), data = drace24_indexed)
model_ACE_idx_log <- lmer(log(Avg_runs) ~ soil_health_score + (1 | Field), data = drace24_indexed)

shapiro_idx_raw <- shapiro.test(residuals(model_ACE_idx))
shapiro_idx_log <- shapiro.test(residuals(model_ACE_idx_log))
cat("Indexed model - Shapiro-Wilk, raw-scale residuals:  p =", signif(shapiro_idx_raw$p.value, 3), "\n")
cat("Indexed model - Shapiro-Wilk, log-scale residuals:  p =", signif(shapiro_idx_log$p.value, 3), "\n")

use_log_idx <- shapiro_idx_log$p.value > shapiro_idx_raw$p.value
final_model_idx <- if (use_log_idx) model_ACE_idx_log else model_ACE_idx
cat("Using", if (use_log_idx) "LOG-transformed" else "RAW-scale", "model for indexed inference.\n\n")

summary(final_model_idx)
anova_ACE_idx <- anova(final_model_idx)
print(anova_ACE_idx)

#  Visual diagnostic: QQ plot for the selected model 
par(mfrow = c(1, 2))
qqnorm(residuals(final_model_idx), main = "QQ Plot: Indexed ACE Model")
qqline(residuals(final_model_idx), col = "red", lwd = 2)
plot(fitted(final_model_idx), residuals(final_model_idx),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs. Fitted")
abline(h = 0, col = "red", lwd = 2)
par(mfrow = c(1, 1))

# Plot: scatter + fitted line

fixed_est   <- summary(final_model_idx)$coefficients["soil_health_score", "Estimate"]
fixed_p     <- anova_ACE_idx$`Pr(>F)`[1]
r2_marginal <- MuMIn::r.squaredGLMM(final_model_idx)[1, "R2m"]

sig_lab <- paste0("R\u00b2(marginal) = ", round(r2_marginal, 2),
                  "  (", ifelse(fixed_p < 0.001, "p < 0.001", paste0("p = ", round(fixed_p, 3))), ")")

plot_data <- drace24_indexed %>%
  mutate(y_plot = if (use_log_idx) log(Avg_runs) else Avg_runs)

p_indexed <- ggplot(plot_data, aes(x = soil_health_score, y = y_plot)) +
  geom_smooth(method = "lm", color = "gray50", fill = "gray85",
              alpha = 0.3, linewidth = 0.8, se = TRUE) +
  geom_point(aes(color = Type, shape = County), size = 3.5, alpha = 0.9) +
  annotate("text", x = Inf, y = Inf, label = sig_lab,
           hjust = 1.05, vjust = 1.8, size = 4.2, color = "gray30", fontface = "italic") +
  scale_color_manual(values = fill_colors, labels = type_labels, name = "Management") +
  scale_shape_manual(values = c("Grant" = 16, "Fillmore" = 17), name = "County") +
  labs(
    title = "ACE Protein vs. Soil Health Management Score",
    x = "Soil Health Management Score  (higher = more favorable)",
    y = if (use_log_idx) expression(paste("log(g Protein kg"^-1, " soil)")) else expression(paste("g Protein kg"^-1, " soil"))
  ) +
  theme_bw() +
  theme(
    axis.text       = element_text(size = 12, face = "bold", color = "black"),
    axis.title      = element_text(size = 13, face = "bold"),
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title    = element_text(size = 12, face = "bold"),
    legend.text     = element_text(size = 11),
    legend.position = "bottom",
    legend.key.size = unit(0.8, "lines")
  )

print(p_indexed)
ggsave("ACE_vs_soil_health_score.png", p_indexed, width = 8.5, height = 6.5, dpi = 300)

cat("\nSaved ACE_vs_soil_health_score.png\n")

# ---- Metadata export (indexed version) ------------------------------------

ACE_meta_indexed <- data.frame(
  Metric         = "ACE",
  Model          = "indexed (soil_health_score continuous)",
  Scale          = if (use_log_idx) "log" else "raw",
  Shapiro_raw_p  = signif(shapiro_idx_raw$p.value, 3),
  Shapiro_log_p  = signif(shapiro_idx_log$p.value, 3),
  Fixed_Estimate = round(fixed_est, 4),
  LMM_F          = round(anova_ACE_idx$`F value`[1], 3),
  LMM_p          = signif(fixed_p, 3),
  R2_marginal    = round(r2_marginal, 3)
)

write.csv(ACE_meta_indexed, "ACE_metadata_indexed.csv", row.names = FALSE)
cat("Saved ACE_metadata_indexed.csv\n")

# PART 5 — Field-averaged scatter (n=31, simple lm, ordinary R^2)

ace_field_avg <- drace24_indexed %>%
  group_by(Farmer, Type, County, soil_health_score) %>%
  summarise(ACE = mean(Avg_runs, na.rm = TRUE), .groups = "drop")

lm_fit_avg  <- lm(ACE ~ soil_health_score, data = ace_field_avg)
r2_avg      <- summary(lm_fit_avg)$r.squared
pval_avg    <- summary(lm_fit_avg)$coefficients[2, 4]
sig_lab_avg <- paste0("R\u00b2 = ", round(r2_avg, 2),
                      "  (", ifelse(pval_avg < 0.001, "p < 0.001", paste0("p = ", round(pval_avg, 3))), ")")

cat("Field-averaged ACE model: n =", nrow(ace_field_avg),
    ", R2 =", round(r2_avg, 3), ", p =", signif(pval_avg, 3), "\n")

p_field_avg <- ggplot(ace_field_avg, aes(x = soil_health_score, y = ACE)) +
  geom_smooth(method = "lm", color = "gray50", fill = "gray85",
              alpha = 0.3, linewidth = 0.8, se = TRUE) +
  geom_point(aes(color = Type, shape = County), size = 4, alpha = 0.9) +
  annotate("text", x = Inf, y = Inf, label = sig_lab_avg,
           hjust = 1.05, vjust = 1.8, size = 4.2, color = "gray30", fontface = "italic") +
  scale_color_manual(values = fill_colors, labels = type_labels, name = "Management") +
  scale_shape_manual(values = c("Grant" = 16, "Fillmore" = 17), name = "County") +
  labs(
    title = "ACE Protein by Management",
    x = "Soil Health Management Score  (higher = more favorable)",
    y = expression(bold(paste("g Protein kg"^-1, " soil")))
  ) +
  theme_bw() +
  theme(
    axis.text       = element_text(size = 12, face = "bold", color = "black"),
    axis.title      = element_text(size = 13, face = "bold"),
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title    = element_text(size = 12, face = "bold"),
    legend.text     = element_text(size = 11),
    legend.position = "bottom",
    legend.key.size = unit(0.8, "lines")
  )

print(p_field_avg)
ggsave("ACE_field_avg_vs_score.png", p_field_avg, width = 8.5, height = 6.5, dpi = 300)
cat("Saved ACE_field_avg_vs_score.png\n")