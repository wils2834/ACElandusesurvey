# DRACE 2024 — ACE (Autoclaved Citrate Extractable) Protein
library(readxl)
library(dplyr)
library(lme4)
library(lmerTest)   
library(emmeans)
library(multcomp) 
library(ggplot2)
# PART 1 — Import and clean
drace24 <- read_excel("C:\\Users\\swils\\Downloads\\01 Projects\\2024 CRP\\2024DRACE_03042025\\AllplatesDRACE2024.xlsx", sheet = "raw_data")
colnames(drace24) <- make.names(colnames(drace24), unique = TRUE)

# Drop plate-calibration standards — not field samples.
drace24 <- drace24 %>%
  filter(!Type %in% c("Standard", "Standard 2", "Standard 3", "Standard 4"))

drace24 <- drace24 %>% filter(!is.na(Avg_runs))

drace24 <- drace24 %>%
  mutate(Type = recode(Type, "SH" = "Soil_Health", "RG" = "Rotational_Grazing"))

drace24$Type <- factor(
  drace24$Type,
  levels = c("Conventional", "Soil_Health", "Rotational_Grazing", "New_CRP", "Old_CRP")
)
# True field ID (see header note) — use this for all random effects.
drace24 <- drace24 %>%
  mutate(Field = paste(Farmer, Type, sep = "_"))

# Sanity checks
cat("Rows after cleaning:", nrow(drace24), "\n")
cat("Distinct fields:", n_distinct(drace24$Field), "(expect 31)\n")
cat("Distinct sample locations (UID):", n_distinct(drace24$UID), "\n\n")
stopifnot(n_distinct(drace24$Field) == 31)

# PART 2 — Mixed model
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

## ----------------------------------------------------------------------
## PART 3 — Plots
## CLD letters come from the mixed model above, so the plot and the
## stats behind it actually agree (the original script's boxplots used
## a separate, pseudoreplicated ANOVA for the letters).
## ----------------------------------------------------------------------

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

# --- Overall boxplot ---
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
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13, face = "bold", color = "black"),
    axis.text.y = element_text(size = 13, face = "bold", color = "black"),
    axis.title.y = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

print(p_overall)
ggsave("ACE_overall_boxplot.png", p_overall, width = 9, height = 6, dpi = 300)

# --- Faceted by county ---
# Letters are identical across both facets because they come from the
# single overall model (pooled across counties), not a per-county test.
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
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
  )

print(p_county)
ggsave("ACE_by_county_boxplot.png", p_county, width = 11, height = 6, dpi = 300)

cat("\nSaved ACE_overall_boxplot.png and ACE_by_county_boxplot.png\n")

# ---- Metadata export ----------------------------------------------------
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
                LMM_F, LMM_p, Type, Mean, CLD)

write.csv(ACE_meta, "ACE_metadata.csv", row.names = FALSE)
cat("Saved ACE_metadata.csv\n")