# Statistical analyses - main figures (Fig. 3D, 4B, 4D, 4F)
# Kasap et al., PNAS
# Input: long-format CSVs in data/ (one row per device; values as in Source Data)

library(dplyr)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)
library(car)
library(rstatix)

data_dir <- "data"
rd <- function(f) read.csv(file.path(data_dir, f), fileEncoding = "UTF-8")


# Fig. 3D - EECM-BMEC LY permeability, mono- vs triculture conditions, uSiM vs Transwell
# Linear mixed model, independent experiment as random intercept

fig3d <- rd("Fig3D_long.csv")
fig3d %>% count(condition, platform, experiment)

shapiro.test(fig3d$value)
leveneTest(value ~ interaction(condition, platform), data = fig3d)

m3d <- lmer(value ~ condition * platform + (1 | experiment), data = fig3d)
anova(m3d, ddf = "Kenward-Roger")
emmeans(m3d, pairwise ~ condition * platform, adjust = "tukey", lmer.df = "kenward-roger")$contrasts

fig3d %>% group_by(condition, platform) %>%
  summarise(mean = mean(value), sd = sd(value), n_devices = n(), n_exp = n_distinct(experiment), .groups = "drop")


# Fig. 4B - uSiM-NVU vs +dbcAMP vs +VEGFA-Ab
# BMEC_tri (same data as Fig. 3D) is a reference line only, excluded from testing.
# Unequal variances (Levene) -> Welch's ANOVA on experiment means + Holm-corrected pairwise Welch t-tests

fig4b <- rd("Fig4B_long.csv")
fig4b %>% count(condition, experiment)

fig4b_test <- filter(fig4b, condition != "BMEC_tri")
leveneTest(value ~ factor(condition), data = fig4b_test)

agg4b <- fig4b_test %>% group_by(experiment, condition) %>% summarise(value = mean(value), .groups = "drop")
welch_anova_test(agg4b, value ~ condition)
pairwise.t.test(agg4b$value, agg4b$condition, pool.sd = FALSE, p.adjust.method = "holm")

fig4b %>% group_by(condition) %>%
  summarise(mean = mean(value), sd = sd(value), n_devices = n(), n_exp = n_distinct(experiment))


# Fig. 4D - HC vs MS donor lines, uSiM vs Transwell
# Linear mixed model, donor line as random intercept

fig4d <- rd("Fig4D_long.csv")
fig4d %>% count(disease, platform, donor)

shapiro.test(fig4d$value)
leveneTest(value ~ interaction(disease, platform), data = fig4d)

m4d <- lmer(value ~ disease * platform + (1 | donor), data = fig4d)
anova(m4d, ddf = "Kenward-Roger")
emmeans(m4d, pairwise ~ disease * platform, adjust = "tukey", lmer.df = "kenward-roger")$contrasts


# Fig. 4F - PBMC transmigration, HC (2B) vs MS (9B) uSiM-NVU
# Linear mixed model, independent experiment as random intercept

fig4f <- rd("Fig4F_long.csv")
fig4f$disease <- factor(fig4f$disease, levels = c("HC", "MS"))
fig4f %>% count(disease, experiment)

by(fig4f$value, fig4f$disease, shapiro.test)
leveneTest(value ~ disease, data = fig4f)

m4f <- lmer(value ~ disease + (1 | experiment), data = fig4f)
anova(m4f, ddf = "Kenward-Roger")
emm4f <- emmeans(m4f, ~ disease, lmer.df = "kenward-roger")
emm4f
confint(pairs(emm4f))   # mean difference HC - MS with 95% CI (Kenward-Roger df)

# sensitivity: log-transformed values (MS group deviates from normality on raw scale)
by(log(fig4f$value), fig4f$disease, shapiro.test)
m4f_log <- lmer(log(value) ~ disease + (1 | experiment), data = fig4f)
anova(m4f_log, ddf = "Kenward-Roger")


sessionInfo()
