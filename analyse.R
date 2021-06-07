# Attentional Blink (AB) meta-analysis

library(tidyverse)
library(effsize)
library(BSDA)
library(psy.phd)
library(brms)
library(tidybayes)
library(forester)

rm(list=ls())

## data

data_dir <- 'data'
studies <- read_csv(paste0(data_dir,'/ab_data.csv'))

## functions

# If SMD, calculate pooled SD as per Cochrane Handbook, which in 6.3.1, says to use Hedge's g
# Standardized Mean Difference (SMD, Cohen's d)
# https://training.cochrane.org/handbook/current/chapter-06#section-6-5-1
# I don't think https://training.cochrane.org/handbook/current/chapter-23#section-23-2-7
# is relevant, as these are parallel group trials, not crossover trials

# set pre-post attentional blink differences for treatment and control groups
ab_set_diff <- function(df) {
  # SDs are pooled from pre and post AB
  # May need "Imputing a change-from-baseline standard deviation using a correlation coefficient" (Higgins et al., 2019)
  df %>% mutate(
               treatment.diff.m  = treatment.pre.m - treatment.post.m,
               treatment.diff.sd = sd_pooled(treatment.n, treatment.n, treatment.pre.sd, treatment.post.sd),
               control.diff.m    = control.pre.m - control.post.m,
               control.diff.sd   = sd_pooled(control.n, control.n, control.pre.sd, control.post.sd)
               )
}

#' Mean difference
#'
#' Returns (m1 - m2) / sd. If sd is the pooled standard deviation, then this is Hedge's g.
#'
#' @param m1 mean 1
#' @param m2 mean 2
#' @param sd standard deviation
#' @return numeric
d <- function(m1, m2, sd) {
  # https://www.statisticshowto.com/hedges-g/
  (m1 - m2) / sd
}

#' Pooled standard deviation
#'
#' Returns the pooled standard deviation for two groups with same or different n.
#'
#' @param n1 n for group 1
#' @param n2 n for group 2
#' @param sd1 standard deviation for group 1
#' @param sd2 standard deviation for group 2
#' @return numeric
sd_pooled <- function(n1, n2, sd1, sd2) {
  # https://www.statisticshowto.com/pooled-standard-deviation/
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

#' Convert 95% confidence interval to standard error
#'
#' https://training.cochrane.org/handbook/current/chapter-06#section-6-3-1
#'
#' @param u upper interval
#' @param l lower interval
#' @return numeric
ci95_to_se <- function(u, l) { (u - l) / 3.92 }

#' 95% confidence interval for Cohen's d (Rosnow & Rosenthal, 2009)
#'
#' @param d Cohen's d
#' @param n1 n in group1
#' @param n2 n in group2
#' @return numeric
d_ci95 <- function(d, n1, n2) {
  df <- n1 + n2 - 2
  sqrt((((n1 + n2) / (n1 * n2)) + (d^2 / (2 * df))) * ((n1 + n2) / df))
}

#' Set effect size for a study
#'
#' Calculates mean difference and 95% confidence interval.
#'
#' @param study data frame with publication column
#' @param m1 mean for group 1
#' @param m2 mean for group 2
#' @param sd standard deviation
#' @param n1 n in group 1
#' @param n2 n in group 2
#' @return tibble
#'
#' Assumes data = (study, mean1, mean2, sd, n1, n2, group)
set_effect <- function(data) {
  v <- select(data, 2:7)
  oldnames <- names(v)
  newnames <- c('mean1', 'mean2', 'sd', 'n1', 'n2', 'group')
  v <- v %>%
    rename_with(~ newnames[which(oldnames == .x)], .cols = oldnames)
  df <- tibble(
    study  = data$study,
    mean1  = v$mean1,
    mean2  = v$mean2,
    d      = d(v$mean1, v$mean2, v$sd),
    ci = 0, l = 0, u = 0
  ) %>% mutate(
    ci    = d_ci95(d, v$n1, v$n2),
    l     = d - .data$ci,
    u     = d + .data$ci,
    group = v$group
  )
}

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
iter        <- '100e4'
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
rem <- brm(
  d | se(se) ~ 1 + (1 | Study),
  data = model_data,
  chains=8, iter=iter,
  prior = c(prior_string("normal(0,1)", class = "Intercept"),
            prior_string(sd_prior, class = "sd")),
  control = list(adapt_delta = adapt_delta),
  file = 'cache/all-ab-brms'
)

brms_function <- function(model, data = dat, average_effect_label = 'Pooled effect',
                          iter = '50e4', sd_prior = "cauchy(0, .3)", adapt_delta = 0.8, cache_file = 'foo') {
  ## https://github.com/mvuorre/brmstools
  ## https://vuorre.netlify.com/post/2016/09/29/meta-analysis-is-a-special-case-of-bayesian-multilevel-modeling/

  # store study names to avoid clashes with commas in spread_draws() column specifications
  data <- data %>%
    mutate(study_number = as.numeric(rownames(data)))
  study_names <- data %>% select(study_number, Study)

  model <- brm(
    d | se(se) ~ 1 + (1 | study_number),
    data = data,
    chains=8, iter=iter,
    prior = c(prior_string("normal(0,1)", class = "Intercept"),
              prior_string(sd_prior, class = "sd")),
    control = list(adapt_delta = adapt_delta),
    file = paste(cache_file, 'brms', sep = '-')
  )

  # For an explanation of tidybayes::spread_draws(), refer to http://mjskay.github.io/tidybayes/articles/tidy-brms.html
  # Study-specific effects are deviations + average
  draws <- spread_draws(model, r_study_number[study_number, term], b_Intercept) %>%
    rename(b = b_Intercept) %>%
    mutate(b = r_study_number + b) %>%
    left_join(study_names, by = 'study_number')

  # Average effect
  draws_overall <- spread_draws(model, b_Intercept) %>%
    rename(b = b_Intercept) %>%
    mutate(Study = average_effect_label)

  # Combine average and study-specific effects' data frames
  combined_draws <- bind_rows(draws, draws_overall) %>%
    ungroup() %>%
    mutate(Study = fct_relevel(Study, average_effect_label, after = Inf)) # put overall effect after individual studies

  # summarise in metafor format
  metafor <- group_by(combined_draws, Study) %>%
    mean_hdci(b) %>% # FIXME: parameterise interval
    rename(est = b, ci_low = .lower, ci_high = .upper)

  return(metafor)
}

brms_object_to_table <- function(model, table, overall_estimate = FALSE, subset_col = "Overall",
                                 subset_col_order = NULL, iter = '1e4', sd_prior = "cauchy(0, .3)",
                                 adapt_delta = 0.8, cache_label = 'foo') {

  table <- left_join(data.frame(model$data %>% select(Study, d, se)), table, by = 'Study')

  # Reorder data
  table <- select(table, Study, everything())

  # Clean level names so that they look nice in the table
  table[[subset_col]] <- str_to_sentence(table[[subset_col]])
  levels <- unique(table[[subset_col]])
  if(!(is.null(subset_col_order))){
    levels <- intersect(subset_col_order, levels)
  }

  # Work out if only one level is present. Passed to create_subtotal_row(), so
  # that if only one group, no subtotal is created.
  single_group <- ifelse(length(levels)==1, TRUE, FALSE)

  # Subset data by levels, run user-defined metafor function on them, and
  # recombine along with Overall rma output
  subset <- lapply(levels, function(level){filter(table, !!as.symbol(subset_col) == level)})
  names(subset) <- levels

  # model each data subset
  subset_res <- lapply(levels, function(level){brms_function(model, data = subset[[level]],
                                                             iter = iter, sd_prior = sd_prior, adapt_delta = adapt_delta,
                                                             cache_file = paste(cache_label, level, sep = '-'))})
  names(subset_res) <- levels

  # This binds the table together
  subset_tables <-
    lapply(levels, function(level){
      rbind(
        create_title_row(level),
        dplyr::select(subset_res[[level]], Study, .data$est, .data$ci_low, .data$ci_high)
      )
    })

  subset_table <- do.call("rbind", lapply(subset_tables, function(x) x))

  ordered_table <- rbind(subset_table,
                         if (overall_estimate) {
                           create_subtotal_row(rma, "Overall", add_blank = FALSE)
                         })

  # Indent the studies for formatting purposes
  ordered_table$Study <- as.character(ordered_table$Study)
  ordered_table$Study <- ifelse(!(ordered_table$Study %in% levels) & ordered_table$Study != "Overall",
                                paste0("  ", ordered_table$Study),
                                ordered_table$Study)

  return(ordered_table)
}

# TODO forester plot
library(forester)
# plunder robvis for rem -> forester-friendly
data <- brms_object_to_table(rem, effects %>% select(Study, group), cache_label = 'cache/all-ab', subset_col = 'group',
                             iter = iter, sd_prior = sd_prior, adapt_delta = adapt_delta)

forester(data,
               subset_col = 'group',
               estimate_col_name = estimate_col_name,
               cache_label = 'nrct_ab',
               null_line_at = 0,
               font_family = "serif",
               x_scale_linear = TRUE,
               xlim = c(-.8, 1.3),
               xbreaks = c(-.8, -.5, -.3, 0, .3, .5, .8, 1, 1.3),
               arrows = FALSE,
               arrow_labels = c("Low", "High"),
               nudge_y = -0.2,
               estimate_precision = 2,
               display = FALSE,
               add_tests = FALSE,
               file_path = here::here(paste0("figures/all_ab_forest.png"))
)

# RCTs
rct <- model_data %>% slice(c(1, 5))
rem <- brm(
  d | se(se) ~ 1 + (1 | Study),
  data = rct,
  chains=8, iter=iter,
  prior = c(prior_string("normal(0,1)", class = "Intercept"),
            prior_string(sd_prior, class = "sd")),
  control = list(adapt_delta = adapt_delta),
  file = 'rct-ab-brms'
)
# variance (vi) = sei^2
# yi must be named 'yi' for slab to work
#rct.escalc <- metafor::escalc(measure = "SMD", yi = yi, sei = sei, data = rct, slab = Study)
#rct.rma <- metafor::rma(yi, vi, data = rct.escalc)
rob_blobbogram(rem, rob2,
               iter = iter, sd_prior = sd_prior, adapt_delta = adapt_delta,
               subset_col = 'group',
               estimate_col_name = estimate_col_name,
               cache_label = 'rct_ab',
               null_line_at = 0,
               font_family = "serif",
               rob_colour = "cochrane",
               rob_tool = "ROB2",
               x_scale_linear = TRUE,
               xlim = c(-.3, 1),
               xbreaks = c(-.3, 0, .3, .5, .8, 1),
               arrows = FALSE,
               arrow_labels = c("Low", "High"),
               nudge_y = -0.2,
               estimate_precision = 2,
               display = FALSE,
               add_tests = FALSE,
               file_path = here::here(paste0("figures/rct_ab_forest.png"))
)

read_brms_model <- function(name) {
  readRDS(paste0(name, '.rds'))
}

p_effect <- tibble(type = '', experience = '', small = 0, medium = 0, large = 0, neg_small = 0, neg_medium = 0)
for (type in c('rct_ab', 'nrct_ab')) {
    for (experience in c('Novice', 'Experienced')) {
      name <- paste(type, experience, 'brms', sep = '-')
      if (name %in% c('rct_ab-Experienced-brms')) next
      model <- read_brms_model(name)
      # pooled effect
      draws <- spread_draws(model, b_Intercept) %>% rename(b = b_Intercept)
      df <- draws %>% summarise(small = mean(b > 0.2), medium = mean(b > 0.5), large = mean(b > 0.8),
                                neg_small = mean(b < -.2), neg_medium = mean(b < -.5))
      p_effect <- add_row(p_effect, type = type, experience = experience,
                          small = df$small, medium = df$medium, large = df$large,
                          neg_small = df$neg_small, neg_medium = df$neg_medium)
    }
}
p_effect <- slice(p_effect, 2:n())
saveRDS(p_effect, paste0(data_dir, '/p_effect.Rd'))

## check rhat
for (type in c('rct_ab', 'nrct_ab')) {
  for (experience in c('Novice', 'Experienced')) {
    name <- paste(type, experience, 'brms', sep = '-')
    if (name %in% c('rct_ab-Experienced-brms')) next
    model <- read_brms_model(name)
    print(name)
    print(summary(model)) # check rhat
  }
}

# for (var in c('md', 'smd')) {
#   for (subgroup in c('novice', 'experienced')) {
#     if (var == 'md') {
#       xlab <- 'Mean Difference'
#     } else {
#       xlab <- 'Standardised Mean Difference'
#     }
#     # calculate SEM from CI
#     model_data <- effects %>%
#       filter(score == var & group == subgroup) %>%
#       mutate(se = ci95_to_se(u, l), study = factor(study))
#
#     # this shows you the default priors
#     get_prior(d | se(se) ~ 1 + (1 | study), data=model_data)
#
#     model_data$study <- str_replace(model_data$study, ",", "")  # remove commas in study names
#     # to rebuild model delete cached file var'-rem'
#     rem <- brm(
#       d | se(se) ~ 1 + (1 | study),
#       data = model_data,
#       chains=8, iter=10e4,
#       prior = c(prior_string("normal(0,1)", class = "Intercept"),
#         prior_string("cauchy(0, .5)", class = "sd")),
#       file = paste(subgroup, var, 'rem', sep = '-')
#     )
#
#     # Vuorre also used control=list(adapt_delta = .99). Something to do with non-convergence.
#
#     # extract the posterior samples...
#     posterior_samples <- rem %>% as.data.frame()
#
#     # this dataframe has one column per parameter in the model
#     # (inluding the random effects, so you get each study's divergence from the mean too)
#     posterior_samples %>% names
#
#     # posterior credible interval
#     posterior_samples %>% select(b_Intercept) %>% mean_qi()
#
#     # posterior plot
#     posterior_samples %>%
#       ggplot(aes(b_Intercept)) +
#       geom_density() + xlab("Pooled effect size (posterior density)")
#
#     # test of the hypothesis that the pooled effect is larger than .3 or smaller than -.3
#     # the Evid.Ratio here is a BayesFactor so we have reasonable evidence against.
#     hypothesis(rem, "abs(Intercept)>.3")
#
#     # or calculate the other way. BF=4 that the effect is smaller than < .3, even with so few studies
#     hypothesis(rem, "abs(Intercept)<.3")
#
#     forest <- forest_bayes(rem)
#     # save plot
#     ggsave(filename = paste0('figures/', subgroup, '_', var, '_forest.pdf'), plot = forest,
#            units = "in", width = 8.27, height = 11.69)
#
#     # Rhat = 1.0 is good
#     # Effective Sample Size (ESS):chains*iter-warmp should be above some threshold
#     # Group-Level, Population-Level or both?
#     print(paste(subgroup, var, sep = ':'))
#     print(rem)
#     # prior_summary(rem)
#   }
# }

