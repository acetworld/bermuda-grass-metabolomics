#########################################
#                                       #
#       Bermuda Grass Metabolomics      #
#                                       #
#########################################


# Analysis of Nematode Count and Root Biomass Data

# ACET Lab
# March 2019

source('./analysis/packages.R')


# Load in data -------------------------------------

# Read in Counts of Sting nematode and biomass of roots
nematode_counts <- read_csv('./data/nematode_counts.csv', skip = 5)

# Nematode Counts Analysis-------------------------------------

# Visualize nematode counts
set.seed(24)
dodge_width = 0.48
jitter_width = 0.06
fig1_nematode_counts <- ggplot(nematode_counts, aes(x = line, y = nematode_count, color = treatment)) + 
    geom_point(position = position_jitterdodge(jitter.width = jitter_width,
                                               dodge.width = dodge_width),
               alpha = 0.72) + 
    stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar',
                 position = position_dodge(dodge_width), 
                 width = dodge_width/2) + 
    stat_summary(fun.data = 'mean_cl_boot', geom = 'point',
                 position = position_dodge(dodge_width)) + 
    scale_color_grey(labels = c('Inoculated', 'Uninoculated')) +
    theme_bw() + 
    theme(legend.title = element_blank(),
          legend.position = c(0.78, 0.124),
          legend.direction = 'horizontal',
          legend.background = element_rect(color = 'grey20'),
          panel.border = element_blank(),
          axis.line.x = element_line(size = 0.48, linetype = 'solid', 
                                     color = 'black'),
          axis.line.y = element_line(size = 0.48, linetype = 'solid', 
                                     color = 'black'),
          legend.margin = ggplot2::margin(1, 1, 1, 1)) + 
    labs(x = 'Line', y = 'Sting Nematode Count')

ggsave(plot = fig1_nematode_counts, filename = './figures/raw-figures/fig1a-nematode-counts.pdf',
       width = 5, height = 3)


# Model Nematode Counts

# baseline model for comparison
baseline_count_mod = lm(nematode_count ~ 1,
                        data = nematode_counts)

# model effects of line, treatment, and their interaction on nematode counts
nematode_count_mod <- lm(nematode_count ~ line * treatment, 
                         data = nematode_counts)

# Model summary and diagnosis
summary(nematode_count_mod)
Anova(nematode_count_mod, type = 'III')
AIC(nematode_count_mod)
BIC(nematode_count_mod)
lrtest(nematode_count_mod, baseline_count_mod)
qqPlot(nematode_count_mod$residuals)
plot(nematode_count_mod, which = 1:6)

# Post-hoc analysis of nematode counts
nematode_count_emm <- emmeans(nematode_count_mod, ~ line * treatment)
nematode_count_emm %>% contrast() %>% cld()

contrast(nematode_count_emm, by = 'treatment', method = 'pairwise')

# Root Biomass Analysis --------------------------------


# Visualize Nematode Weights
set.seed(24)
dodge_width = 0.58
jitter_width = 0.06
# note: outliers removed - see modeling analysis below
fig2_root_weights <- ggplot(nematode_counts[-c(31,30),], aes(x = line, y = weight, color = treatment)) +
    geom_point(position = position_jitterdodge(jitter.width = jitter_width,
                                               dodge.width = dodge_width),
               alpha = 0.72) +
    stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar',
                 position = position_dodge(dodge_width),
                 width = dodge_width/2) + 
    stat_summary(fun.data = 'mean_cl_boot', geom = 'point',
                 position = position_dodge(dodge_width)) + 
    scale_color_grey(labels = c('Inoculated', 'Uninoculated')) + 
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = c(0.9, 0.124),
          legend.background = element_rect(color = 'grey30'),
          panel.border = element_blank(),
          axis.line.x = element_line(size = 0.48, linetype = 'solid', 
                                     color = 'black'),
          axis.line.y = element_line(size = 0.48, linetype = 'solid', 
                                     color = 'black'),
          legend.margin = ggplot2::margin(1, 1, 1, 1)) + 
    labs(x = 'Line', y = 'Root Biomass (g)')

ggsave(plot = fig2_root_weights, filename = './figures/raw-figures/fig1b-root-weights.pdf',
       width = 5, height = 3)

# Model root biomass effects 

# baseline mass model for comparison
baseline_mass_mod <- lm(weight ~ 1,
                        data = nematode_counts)

# Model effect of line, treatment, and nematode counts on root biomass
mass_mod <- lm(weight ~ line * treatment + nematode_count, 
               data = nematode_counts)

# Model summary and diagnostics
summary(mass_mod)
Anova(mass_mod, type = 'III')
AIC(mass_mod)
BIC(mass_mod)
lrtest(mass_mod, baseline_mass_mod)
qqPlot(mass_mod$residuals)
plot(mass_mod, which = 1:6)

# it looks like we have two pretty large outliers. Let's double check that.
outlierTest(mass_mod)


# Test removing observation 31
mass_mod <- lm(weight ~ line * treatment + nematode_count, 
               data = nematode_counts[-31,])

summary(mass_mod)
Anova(mass_mod, type = 'III')
AIC(mass_mod)
BIC(mass_mod)
qqPlot(mass_mod$residuals)
plot(mass_mod, which = 1:6)
outlierTest(mass_mod)

# test removing both observations
counts_no_outlier <- nematode_counts %>% slice(-c(30,31))

# new baseline model
baseline_mass_no_outlier <- lm(weight ~ 1, data = counts_no_outlier)

# Model without outliers
mass_mod_no_outlier <- lm(weight ~ line + treatment, 
                          data = counts_no_outlier)

# Model summary and diagnostics
summary(mass_mod_no_outlier)
Anova(mass_mod_no_outlier, type = 'III')
AIC(mass_mod_no_outlier)
BIC(mass_mod_no_outlier)
lrtest(mass_mod_no_outlier, baseline_mass_no_outlier)
qqPlot(mass_mod_no_outlier$residuals)
plot(mass_mod_no_outlier, which = 1:6)
outlierTest(mass_mod_no_outlier)

# post-hoc
mass_emm <- emmeans(mass_mod_no_outlier, ~ line + treatment)
mass_emm %>% contrast() %>% cld()

contrast(mass_emm, by = 'line', method = 'pairwise')
contrast(mass_emm, by = 'treatment', method = 'pairwise')


# Root biomass loss analysis -------------------------------

# We'll use bootstrapping to get at estimates of loss 
# since we do not have exact values

# number of bootstrap reps
bootstrap_reps = 1000

# select only those plants with nematodes
inoculated_plants <- filter(counts_no_outlier, treatment == 'inoc')

# select only those plants without nematodes for baseline
uninoculated_plants <- filter(counts_no_outlier, treatment == 'uninoc')

# holder for bootstrapped data:
bootstrapped_loss <- tibble()
bootstrapped_loss_model <- tibble()

# function for fitting model to results
loss_model <- function(df) {
    lm(loss ~ nematode_count, data = df)
}

# Set Seed for reproducibility
set.seed(72)

# Bootstrap loss by sampling with replacement uninoculated plants,
# calculating baseline mean, then subtracting that mean from observed 
# plants inoculated with nematodes
for (i in seq(1:bootstrap_reps)) {
    
    # sample with replacement and calculate means
    uninoculated_means <- uninoculated_plants %>%
        group_by(line) %>% 
        sample_frac(replace = TRUE) %>%
        summarize(uninoc_mn = mean(weight))
    
    # match means by line and calculate loss
    loss_of_biomass <- left_join(inoculated_plants, uninoculated_means) %>%
        mutate(loss = uninoc_mn - weight,
               bootstrap_rep = i)
    
    # record samples for analysis
    bootstrapped_loss <- bind_rows(bootstrapped_loss, loss_of_biomass)
    
    # model relationship between nematode counts and weights
    biomass_loss_model <- loss_of_biomass %>%
        mutate(line = as.factor(line)) %>%
        select(line, nematode_count, loss) %>% 
        group_by(line) %>% 
        nest() %>%
        mutate(model = purrr::map(data, loss_model),
               model_coefs = purrr::map(model, coef),
        ) %>%
        select(line, model_coefs) %>%
        unnest() %>%
        mutate(bootstrap_rep = i,
               coefficient = rep(c('intercept', 'nematode_count'),3))
    
    # store samples
    bootstrapped_loss_model <- bind_rows(bootstrapped_loss_model, biomass_loss_model)
    
}

# Examine medians of bootstrapped loss for comparison
bootstrapped_loss_plot <- bootstrapped_loss %>%
    select(line, bootstrap_rep, loss) %>%
    group_by(line, bootstrap_rep) %>%
    summarize(median = median(loss)) %>%
    group_by(line) %>%
    summarize(med_med = median(median),
              low = quantile(median, 0.025),
              high = quantile(median, 0.975))

# Permutation Test on estimated root loss
bootstrap_test <- bootstrapped_loss %>% 
    select(line, bootstrap_rep, loss) %>% 
    group_by(line, bootstrap_rep) %>% 
    summarize(mn = median(loss)) %>% 
    spread(line, mn) %>%
    # one-sided comparisons here:
    mutate(AB39_AB03 = AB39 > AB03,
           AB39_AB33 = AB39 > AB33,
           AB33_AB03 = AB33 > AB03) %>%
    select(AB39_AB03, AB39_AB33, AB33_AB03) %>%
    gather(key = 'test', value = 'comparison_result') %>%
    group_by(test) %>%
    summarize(val = mean(comparison_result)) %>%
    # conservatively adjust for multiple comparisons
    mutate(bonferroni_correction = val * 3)

# Display results
error_bar_width = 0.48

fig3_root_loss <- ggplot(bootstrapped_loss_plot, 
                         aes(x = line, y = med_med, color = line)) + 
    geom_point() + 
    geom_errorbar(aes(x = line, ymin = low, ymax = high),
                  width = error_bar_width) + 
    theme_bw() +
    scale_color_grey() + 
    labs(x = 'Line', y = 'Biomass Loss (g)') + 
    theme(legend.position = 'none',
          panel.border = element_blank(),
          axis.line.x = element_line(size = 0.48, linetype = 'solid', 
                                     color = 'black'),
          axis.line.y = element_line(size = 0.48, linetype = 'solid', 
                                     color = 'black'))

ggsave(plot = fig3_root_loss, filename = './figures/raw-figures/fig1c-root-loss.pdf',
       width = 5, height = 3)
