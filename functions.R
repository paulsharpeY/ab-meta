# meta-analysis functions

# If SMD, calculate pooled SD as per Cochrane Handbook, which in 6.3.1, says to use Hedge's g
# Standardized Mean Difference (SMD, Cohen's d)
# https://training.cochrane.org/handbook/current/chapter-06#section-6-5-1
# I don't think https://training.cochrane.org/handbook/current/chapter-23#section-23-2-7
# is relevant, as these are parallel group trials, not crossover trials

## BOOKMARK generic functions

get_hostname <- function(){
  return(as.character(Sys.info()["nodename"]))
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

  # free memory?
  rm(model, draws, draws_overall, combined_draws)
  gc()

  return(metafor)
}

create_title_row <- function(title){
  return(data.frame(Study = title, est = NA, ci_low = NA, ci_high = NA))
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

## BOOKMARK: Attentional Blink (AB) functions

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

## BOOKMARK: Attention Network Test (ANT) functions

#' Compute ANT scores for alerting, orienting and conflict.
#'
#' @param df Data frame
#' @importFrom dplyr group_by left_join mutate select
#' @importFrom forcats fct_relevel
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer pivot_wider unite
#' @export
#' @return Data frame
#'
ant_scores <- function(df) {
  alerting_orienting <- df %>%
    pivot_wider(id_cols = c(.data$p,.data$group,.data$t), names_from = .data$cue,
                values_from = .data$rt, values_fn = list(rt = mean)) %>%
    mutate(alerting = .data$nocue - .data$double, orienting = .data$center - .data$spatial) %>%
    select(.data$p, .data$group, .data$t, .data$alerting, .data$orienting)
  conflict <- df %>%
    pivot_wider(id_cols = c(.data$p,.data$group,.data$t),
                names_from = .data$flanker_type, values_from = .data$rt,
                values_fn = list(rt = mean)) %>%
    mutate(conflict = .data$incongruent - .data$congruent) %>%
    select(.data$p, .data$group, .data$t, .data$conflict)
  result <- left_join(alerting_orienting, conflict, by=c('p', 'group', 't')) %>%
    pivot_longer(cols = c(.data$alerting, .data$orienting, .data$conflict),
                 names_to = 'var', values_to = 'rt') %>%
    mutate(var = factor(var))
  # arrange for plot facets to be LtR: Alerting, Orienting, Conflict
  result$var <- fct_relevel(result$var, 'conflict', after = Inf)
  return(result)
}
