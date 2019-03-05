#########################################
#                                       #
#       Bermuda Grass Metabolomics      #
#                                       #
#########################################


# ACET Lab
# Denis Willett
# March 2019

source('./analysis/packages.R')


# Load in data -------------------------------------

# Read in Counts of Sting nematode and biomass of roots
nematode_counts <- read_csv('./data/nematode_counts.csv', skip = 5)

# Raw Data Visualization -------------------------------------

# Visualize nematode counts
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
          legend.position = c(0.75, 0.124),
          legend.direction = 'horizontal',
          legend.background = element_rect(color = 'grey30')) + 
    labs(x = 'Line', y = 'Sting Nematode Count')

ggsave(plot = fig1_nematode_counts, filename = './figures/fig1-nematode-counts.pdf',
       width = 8, height = 5)

# Visualize Nematode Weights
dodge_width = 0.48
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
          legend.background = element_rect(color = 'grey30')) + 
    labs(x = 'Line', y = 'Root Biomass (g)')

ggsave(plot = fig2_root_weights, filename = './figures/fig2-root-weights.pdf',
       width = 8, height = 5)

# Modeling ----------------------------------

# Model effect of line and treatment on nematode counts
baseline_count_mod = lm(nematode_count ~ 1,
                  data = nematode_counts)

nematode_count_mod <- lm(nematode_count ~ line * treatment, 
                         data = nematode_counts)
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

# Model effect of line, treatment, and nematode counts on root biomass

baseline_mass_mod <- lm(weight ~ 1,
                        data = nematode_counts)

mass_mod <- lm(weight ~ line * treatment + nematode_count, 
               data = nematode_counts)

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
lrtest(mass_mod, baseline_mass_mod)
qqPlot(mass_mod$residuals)
plot(mass_mod, which = 1:6)
outlierTest(mass_mod)

# test removing both observations
counts_no_outlier <- nematode_counts %>% slice(-c(30,31))

baseline_mass_no_outlier <- lm(weight ~ 1, data = counts_no_outlier)

mass_mod_no_outlier <- lm(weight ~ line + treatment + nematode_count, 
               data = counts_no_outlier)

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

# Bootstrapping estimates of loss -------------------------
bootstrap_reps = 1000

# select only those plants with nematodes
inoculated_plants <- filter(counts_no_outlier, treatment == 'inoc')

# select only those plants without nematodes for baseline
uninoculated_plants <- filter(counts_no_outlier, treatment == 'uninoc')

# holder for bootstrapped data:
bootstrapped_loss <- tibble()

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
    
}

ggplot(bootstrapped_loss, aes(x = nematode_count, y = loss, 
                              color = line, shape = line)) +
    geom_point() + 
    stat_smooth(method = 'lm') + 
    scale_color_grey() + 
    theme_bw()



---------------------------------------






ggplot(nematode_counts, 
       aes(x = nematode_count, y = weight, 
           color = treatment, shape = line)) + 
    geom_point() + 
    scale_color_grey() + 
    theme_bw()

ggplot(filter(nematode_counts, treatment == 'inoc'), 
       aes(x = nematode_count, y = weight, color = line)) + 
    geom_point() + 
    scale_color_grey() + 
#    stat_smooth(method = 'lm') + 
    theme_bw()


baseline <- lm(weight ~ line, data = filter(nematode_counts, treatment == 'uninoc'))
summary(baseline)
Anova(baseline)
plot(baseline)        

predict(baseline, newdata = filter(nematode_counts, treatment == 'inoc'))

# Bootstrapping estimates of loss -------------------------
bootstrap_reps = 1000

# select only those plants with nematodes
inoculated_plants <- filter(nematode_counts[-c(30, 31),], treatment == 'inoc')

# select only those plants without nematodes for baseline
uninoculated_plants <- filter(nematode_counts[-c(30, 31),], treatment == 'uninoc')

# holder for bootstrapped data:
bootstrapped_loss <- tibble()

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
    
}

ggplot(bootstrapped_loss, aes(x = nematode_count, y = loss, 
                              color = line, shape = line)) +
    geom_point() + 
    stat_smooth(method = 'lm') + 
    scale_color_grey() + 
    theme_bw()

weight_mod <- lm(weight ~ line + treatment + nematode_count, data = nematode_counts[-c(30, 31),])
summary(weight_mod)
Anova(weight_mod, type = 'III')
plot(weight_mod, which = 1:6)
qqPlot(weight_mod$residuals)
outlierTest(weight_mod)

test_diff <- emmeans(weight_mod, ~ line + treatment)
contrast(test_diff, method = 'pairwise') %>% cld()
