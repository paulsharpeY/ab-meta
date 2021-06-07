# Attentional Blink (AB) meta-analysis

library(tidyverse)
library(effsize)
library(brms)
library(tidybayes)
library(forester)

rm(list=ls())
source('functions.R')

## data

data_dir <- 'data'
studies <- read_csv(paste0(data_dir,'/ab_data.csv'))

### main ###

## BOOKMARK: Sharpe et al. (2021)
# FA+OM (80), relaxation (40)
sharpe <- read_csv(paste0(data_dir,'/sharpe_et_al_2021.csv')) %>%
  pivot_wider(names_from = group, values_from = c(n, mean_ab, sd_ab)) %>%
  mutate(effect_sd = sd_pooled(n_meditation, n_control, sd_ab_meditation, sd_ab_meditation)) %>% # FIXME: wrong sd?
           select(study, mean_ab_control, mean_ab_meditation, effect_sd, n_meditation, n_control, experience)
smd <- set_effect(sharpe)

## BOOKMARK: May et al. experiment 1
# Meditation - Control > 0 means post-test AB change > in Meditation group
may1 <- read_csv(paste0(data_dir, '/may_et_al_2011_1.csv')) %>%
  mutate(mean_ab_diff = mean_ab_diff * 100, sd_ab_diff = sd_ab_diff * 100) %>% # convert to %
  pivot_wider(names_from = group, values_from = c(n, mean_ab_diff, sd_ab_diff)) %>%
  mutate(effect_sd = sd_pooled(n_Meditation, n_Control, sd_ab_diff_Meditation, sd_ab_diff_Control)) %>% # FIXME: wrong sd?
  select(study, mean_ab_diff_Meditation, mean_ab_diff_Control, effect_sd, n_Meditation, n_Control, experience)
smd <- bind_rows(smd, set_effect(may1))

## May et al. experiment 2
# Meditation - Control > 0 means post-test AB change > in Meditation group
# Two participants were excluded from the meditation condition
may2 <- read_csv(paste0(data_dir, '/may_et_al_2011_2.csv')) %>%
  mutate(mean_ab = mean_ab * 100, sd_ab = sd_ab * 100) %>% # convert to %
  pivot_wider(names_from = group, values_from = c(n, mean_ab, sd_ab)) %>%
  mutate(effect_sd = sd_pooled(n_Meditation, n_Control, sd_ab_Meditation, sd_ab_Control)) %>% # FIXME: wrong sd?
  select(study, mean_ab_Meditation, mean_ab_Control, effect_sd, n_Meditation, n_Control, experience)
smd <- bind_rows(smd, set_effect(may2))

## BOOKMARK: Slagter et al. (2007)
# Meditation - Control > 0 means post-test AB change > in Meditation group
# slagter <- studies %>% filter(grepl('Slagter et al', publication))
# foo <- slagter %>%
#   select(treatment1.n:full.fa.short.sd) %>%
#   rename(treatment.n = treatment1.n) %>%
#   pivot_longer(everything(), names_to = 'var') %>%
#   drop_na() %>%
#   pivot_wider(names_from = var)
# write_csv(foo, paste0(data_dir, '/slagter_et_al_2007.csv'))

slagter <- read_csv(paste0(data_dir, '/slagter_et_al_2007.csv')) %>%
  # SDs are pooled from long and short lags for treatment and control conditions
  mutate(
    treatment.pre.m   = treatment.pre.long.m - treatment.pre.short.m,
    treatment.pre.sd  = sd_pooled(treatment.n, treatment.n, treatment.pre.short.sd, treatment.pre.long.sd),
    control.pre.m     = control.pre.long.m - control.pre.short.m,
    control.pre.sd    = sd_pooled(control.n, control.n, control.pre.short.sd, control.pre.long.sd),
    treatment.post.m  = treatment.post.long.m - treatment.post.short.m,
    treatment.post.sd = sd_pooled(treatment.n, treatment.n, treatment.post.short.sd, treatment.post.long.sd),
    control.post.m    = control.post.long.m - control.post.short.m,
    control.post.sd   = sd_pooled(control.n, control.n, control.post.short.sd, control.post.long.sd)
  ) %>%
  select(study, treatment.n, control.n, treatment.pre.m, treatment.pre.sd, treatment.post.m, treatment.post.sd,
         control.pre.m, control.pre.sd, control.post.m, control.post.sd, experience) %>%
  ab_set_diff() %>%
  mutate(effect_sd = sd_pooled(treatment.n, control.n, treatment.diff.sd, control.diff.sd)) %>%
  select(study, treatment.diff.m, control.diff.m, effect_sd, treatment.n, control.n, experience)
smd <- bind_rows(smd, set_effect(slagter))

## BOOKMARK: Colzato et al.
# FAM vs. OMM
# N=60, group=3, lag=4
# Denominator df: (60-1)*(4-1)-(4-1)*(3-1) = 171
# May not be reliable (https://stats.stackexchange.com/questions/139996/computing-the-true-effect-size-from-f-values-and-sample-sizes)
# especially for interactions. Note also that this is across all (4) lags, not 3 vs. 8
colzato <- read_csv(paste0(data_dir, '/colzato_et_al_2015.csv'))

# FAM * OMM (no control): F(3,114) = 3.07, p =.03, MSE = 0.012, partial eta squared = 0.08
# lag3-lag8 slope is steeper from FAM than OMM, indicating that +ve d returned by fes() is FAM-OMM
# i.e. FAM.AB > OMM.AB
library(compute.es) # F -> d
colzato <- bind_cols(colzato,
                     fes(f=colzato$f, n.1=colzato$treatment1.n, n.2=colzato$treatment2.n) %>%
                       select(g, l.d, u.d))
smd <- add_row(smd, study=colzato$study, d=colzato$g, l=colzato$l.d, u=colzato$u.d, group=colzato$experience)

## BOOKMARK: van Vugt & Slagter (2014)
## van Vugt et al.
# Prediction in all samples is that OM.AB < FA.AB, therefore FA - OM > 0 supports prediction
# short=lag4=336ms (see data in email sent by Marieke van Vugt)
# van_vugt <- studies %>% filter(grepl('van Vugt', publication)) %>%
#   select(treatment1.n:last_col()) %>%
#   pivot_longer(cols = everything(), names_to = 'var') %>%
#   drop_na() %>%
#   pivot_wider(names_from = var) %>%
#   mutate(across(full.fa.short.m:low.om.long.m, function(x) { x * 100}))
# write_csv(van_vugt, paste0(data_dir, '/van_vugt_slagter_2014.csv'))

van_vugt <- read_csv(paste0(data_dir, '/van_vugt_slagter_2014.csv')) %>%
  mutate(
    high.fa.m  = high.fa.long.m - high.fa.short.m,
    high.fa.sd = sd_pooled(treatment2.n, treatment2.n, high.fa.short.sd, high.fa.long.sd),
    high.om.m  = high.om.long.m - high.om.short.m,
    high.om.sd = sd_pooled(treatment2.n, treatment2.n, high.om.short.sd, high.om.long.sd),
) %>%
  mutate(effect_sd = sd_pooled(treatment2.n, treatment2.n, high.om.sd, high.fa.sd),
         treatment1.n = treatment2.n) %>% # high experience subsample; n=FA=OM=15
  select(study, high.fa.m, high.om.m, effect_sd, treatment1.n, treatment2.n, experience)
smd <- bind_rows(smd, set_effect(van_vugt))

## BOOKMARK: Braboszcz et al. (see process.R)
brab <- read_csv(paste0(data_dir, '/braboszcz_et_al_2013.csv')) %>%
  pivot_wider(names_from = time, values_from = c(m, sd)) %>%
  mutate(effect_sd = sd_pooled(n, n, sd_1, sd_2),
         n2 = n) %>%
  select(study, m_1, m_2, effect_sd, n, n2, experience)
smd <- bind_rows(smd, set_effect(brab))

## meta-analysis

effects <- smd %>%
  mutate(se = ci95_to_se(u, l)) %>%
  rename(Study = study) %>%
  mutate(Study = as.character(Study))
write.csv(effects, file = paste0(data_dir, '/ab_meta.csv'), row.names = FALSE)

estimate_col_name <- 'AB SMD estimate'

# calculate SEM from CI
model_data <- effects %>%
  select(Study, d, se, ci, l, u, mean1, mean2, group)

# brm config
iter        <- '50e4' # increase to 100e4 for final results
adapt_delta <- 0.99 # Should normally be 0.8 (default) < adapt_delta < 1
sd_prior    <- "cauchy(0, .3)"

## TBC: model all studies for funnel plots
# model_data <- model_data %>%
#   mutate(study_number = as.numeric(rownames(model_data)))
# study_names <- model_data %>% select(study_number, Study)
#
# rem <- brm(
#   d | se(se) ~ 1 + (1 | study_number), # random effects meta-analyses model (see brmsformula)
#   data = model_data,
#   chains=8, iter=iter,
#   prior = c(prior_string("normal(0,1)", class = "Intercept"),
#             prior_string(sd_prior, class = "sd")),
#   control = list(adapt_delta = adapt_delta),
#   file = paste('ab-brms', sep = '-')
# )

# save data for funnel plots
# effect by study
# study <- rem %>%
#   spread_draws(b_Intercept, r_study_number[study_number,]) %>%
#   mean_hdci(fitted_smd = r_study_number + b_Intercept, .width = c(.95)) %>%  # "Study-specific effects are deviations + average"
#   left_join(study_names, by = 'study_number') %>%
#   left_join(model_data, by = 'Study') %>%
#   select(Study, fitted_smd, se, group)
# # average effect
# pooled <- rem %>%
#   spread_draws(b_Intercept) %>%
#   mean_hdci(fitted_smd = b_Intercept, .width = c(.95)) %>%
#   mutate(Study = "pooled", se_max = max(study$se))
# funnel_data <- bind_rows(study, pooled) # SMD is fitted value, but SE is from original study
# saveRDS(funnel_data, paste0(data_dir, '/ab_funnel_data.Rd'))

# model all studies
rem_all <- brm(
  d | se(se) ~ 1 + (1 | Study),
  data = model_data,
  chains=8, iter=iter,
  prior = c(prior_string("normal(0,1)", class = "Intercept"),
            prior_string(sd_prior, class = "sd")),
  control = list(adapt_delta = adapt_delta),
  file = 'cache/all-ab-brms'
)

all_data <- brms_object_to_table(rem_all, effects %>% select(Study, group), cache_label = 'cache/all-ab', subset_col = 'group',
                             iter = iter, sd_prior = sd_prior, adapt_delta = adapt_delta)

# forest plot for all studies
forester(
  select(all_data, Study),
  all_data$est,
  all_data$ci_low,
  all_data$ci_high,
  estimate_col_name = estimate_col_name,
  null_line_at = 0,
  font_family = "serif",
  x_scale_linear = TRUE,
  xlim = c(-1.3, 1.3),
  xbreaks = c(-1.3, -1, -.8, -.5, -.3, 0, .3, .5, .8, 1, 1.3),
  arrows = FALSE,
  arrow_labels = c("Low", "High"),
  nudge_y = -0.2,
  estimate_precision = 2,
  display = FALSE,
  file_path = here::here(paste0("figures/all_ab_forest.png"))
)

# model all studies except ours
rem_not_ours <- brm(
  d | se(se) ~ 1 + (1 | Study),
  data = model_data %>% filter(Study != 'Sharpe (2021, Experiment 6)'),
  chains=8, iter=iter,
  prior = c(prior_string("normal(0,1)", class = "Intercept"),
            prior_string(sd_prior, class = "sd")),
  control = list(adapt_delta = adapt_delta),
  file = 'cache/not-ours-ab-brms'
)

all_except_our_data <- brms_object_to_table(rem_not_ours, effects %>% select(Study, group),
                                            cache_label = 'cache/not-ours-ab', subset_col = 'group',
                                            iter = iter, sd_prior = sd_prior, adapt_delta = adapt_delta)

# forest plot for all studies
forester(
  select(all_except_our_data, Study),
  all_except_our_data$est,
  all_except_our_data$ci_low,
  all_except_our_data$ci_high,
  estimate_col_name = estimate_col_name,
  null_line_at = 0,
  font_family = "serif",
  x_scale_linear = TRUE,
  xlim = c(-1.3, 1.3),
  xbreaks = c(-1.3, -1, -.8, -.5, -.3, 0, .3, .5, .8, 1, 1.3),
  arrows = FALSE,
  arrow_labels = c("Low", "High"),
  nudge_y = -0.2,
  estimate_precision = 2,
  display = FALSE,
  file_path = here::here(paste0("figures/not_ours_ab_forest.png"))
)

## check rhat
# for (type in c('rct_ab', 'nrct_ab')) {
#   for (experience in c('Novice', 'Experienced')) {
#     name <- paste(type, experience, 'brms', sep = '-')
#     if (name %in% c('rct_ab-Experienced-brms')) next
#     model <- read_brms_model(name)
#     print(name)
#     print(summary(model)) # check rhat
#   }
# }
