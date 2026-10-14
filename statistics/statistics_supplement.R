# Statistical analyses - SI Appendix (Fig. S4A-G)
# Kasap et al., PNAS
# Input: long-format CSVs in data/ (one row per device; values as in Source Data)

library(dplyr)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)
library(car)

data_dir <- "data"
rd <- function(f) read.csv(file.path(data_dir, f), fileEncoding = "UTF-8")


# Fig. S4A-B - Claudin-5 / VE-cadherin MFI, HC (2B) vs MS (9B), untreated vs VEGFA-Ab
# Single experiment per donor -> descriptive only

s4ab <- rd("FigS4AB_long.csv")
s4ab %>% group_by(marker, disease, condition) %>%
  summarise(mean = mean(value), sd = sd(value), n = n(), .groups = "drop")


# Fig. S4C - uSiM, BMEC_tri / BPB / NVU / NVU (SCT-Ast) / NVU (hAST)
# Three conditions are single experiments -> descriptive only

s4c <- rd("FigS4C_long.csv")
s4c %>% group_by(condition) %>%
  summarise(mean = mean(value), sd = sd(value), n_devices = n(), n_exp = n_distinct(experiment))


# Fig. S4D - Transwell, BMEC_tri / BPB / NVU / NVU with BLECs
# Linear mixed model, independent experiment as random intercept.
# Welch's ANOVA not possible (NVU+BLECs is a single experiment); raw Levene's test is
# significant, residual Levene's test after accounting for experiment is not.

s4d <- rd("FigS4D_long.csv")
s4d %>% count(condition, experiment)

leveneTest(value ~ factor(condition), data = s4d)

m_s4d <- lmer(value ~ condition + (1 | experiment), data = s4d)
leveneTest(resid(m_s4d) ~ factor(s4d$condition))
anova(m_s4d, ddf = "Kenward-Roger")
emmeans(m_s4d, pairwise ~ condition, adjust = "tukey", lmer.df = "kenward-roger")$contrasts


# Fig. S4E - mean LY permeability per donor line (data of Fig. 4D), descriptive

fig4d <- rd("Fig4D_long.csv")
fig4d %>% group_by(platform, disease, donor) %>%
  summarise(mean = mean(value), sd = sd(value), n_devices = n(), .groups = "drop") %>%
  arrange(platform, mean)


# Fig. S4F-G - Claudin-5 / VE-cadherin MFI, HC (1B, 2B) vs MS (9B, 15A)
# Replicates averaged to one value per donor -> Welch's t-test on donor means (N = 2 per group)

s4fg <- rd("FigS4FG_long.csv")
s4fg_donor <- s4fg %>% group_by(marker, disease, donor) %>% summarise(value = mean(value), .groups = "drop")
s4fg_donor

t.test(value ~ disease, data = subset(s4fg_donor, marker == "claudin5"))
t.test(value ~ disease, data = subset(s4fg_donor, marker == "VE_cadherin"))


sessionInfo()
