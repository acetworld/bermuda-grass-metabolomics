#########################################
#                                       #
#       Bermuda Grass Metabolomics      #
#                                       #
#########################################


# Analysis of Metabolite Data

# ACET Lab
# March 2019


# Prepare ----------------------------

# Preprocess metabolite data
source('./analysis/preprocess_metabolite_data.R')

# Read in Counts of Sting nematode and biomass of roots
nematode_counts <- read_csv('./data/nematode_counts.csv', skip = 5)

# Reformat data for analysis

# Rearrange data
nematode_treatments <- tibble(
    local_sample_id = colnames(processed_metabolites[-c(1,2)])) %>%
    separate(local_sample_id, 
             into = c('line', 'treatment', 'plant_number'),
             sep = '_',
             remove = FALSE)

# Combine Treatment Data with Metabolites
# and revised names
metabolite_treatment <- processed_metabolites %>% 
    gather(local_sample_id, value, -compounds, -ion_type) %>% 
    left_join(nematode_treatments) %>%
    left_join(compound_names)

# Clean environment
rm(metabolite_imputation, clean_metabolites, metabolite_combined, metabolite_data, negative, positive,
   vsn_model)

# Ordination -------------------------------
#    we want to examine relationship between treatments and lines
#    with observed metabolites.  Canonical correspondence analysis
#    could be a good fit for this but depends upon distributions:

# examine metablite distributions. VSN should help us out here 
ggplot(filter(metabolite_treatment, compounds == '100.1122'), aes(x = value)) + 
    stat_density() + 
    facet_grid(~line)

ggplot(filter(metabolite_treatment, compounds == '245.1386'), aes(x = value)) + 
    stat_density() + 
    facet_grid(~line)

ggplot(filter(metabolite_treatment, compounds == '224.0174'), aes(x = value)) + 
    stat_density() + 
    facet_grid(~line)

# Prepare data for canonical correspondence analysis (CCA)
to_cca <- processed_metabolites %>% 
    gather(local_sample_id, value, -compounds, -ion_type) %>% 
    left_join(nematode_treatments) %>% 
    unite(ion_rt, ion_type, compounds) %>%
    spread(ion_rt, value)

# In this case, metabolites will be our dependent variable, our response
metabolite_response <- to_cca %>% select(-local_sample_id, -line,
                                    -treatment, -plant_number)

# Our independent variables will be line and treatment
line_treatment <- to_cca %>% select(line, treatment)

# Run CCA with interaction 
ccaone <- cca(metabolite_response ~ line * treatment, data = line_treatment)

# Examine Model
# Set Seed for reproducibility
set.seed(776)

# Run Analyses of Variance
anova(ccaone)
anova(ccaone, by = 'term', perm = 1000)
anova(ccaone, by = 'margin', perm = 1000)
anova(ccaone, by = 'axis', perm = 1000)

# Look at explained variance
ccaone$tot.chi
ccaone$CCA$tot.chi

ccaone$CCA$tot.chi/ccaone$tot.chi
ccaone$CCA$eig/ccaone$CCA$tot.chi*100

# Other model diagnostics
summary(ccaone)
vegan::scores(ccaone)
RsquareAdj(ccaone)
goodness(ccaone, summarise = TRUE)
as.mlm(ccaone) 

# visual inspection
plot(ccaone,
     display = c("sp","wa","cn"))
plot(ccaone, display = 'lc')
plot(ccaone, display = 'bp')
plot(ccaone,
     display = c("sp","bp"))

# Get CCA projections
line_treatment_cca <- vegan::scores(ccaone)$sites %>%
    cbind(line_treatment) %>%
    unite(line_treatment, 
          line, treatment, 
          remove = FALSE)

# Determine 95% confidence ellipses for plant projections
out <- dataEllipse(line_treatment_cca$CCA1, line_treatment_cca$CCA2, 
                   groups = as.factor(line_treatment_cca$line_treatment), 
                   levels = 0.95) %>% data.frame() %>%
    gather(H, Values) %>%
    mutate(Coord = str_sub(H, -1, -1),
           Humid = str_sub(H, 1, -3)) %>%
    select(-H) %>% mutate(Coord = as.factor(Coord)) 

# Obtain coordinates for plotting
outx <- filter(out, Coord == 'x') %>% select(-Coord)
outy <- filter(out, Coord == 'y')
ells <- cbind(outx, y =  outy$Values)

# Determine effect vectors and scale for plotting
vectors <- data.frame(ccaone$CCA$biplot * 
                          vecscale(ccaone$CCA$biplot, 
                                   bbox = matrix(c(-5, 3, -5, 5), 2, 2)))
vectors$Name <- rownames(vectors)

# remove TreatmentU because of negligible contribution
# Relabel vectors for display
vectors <- vectors %>%
    mutate(labels = c('AB33', 'AB39', 'Uninoculated',
                      'AB33:Treatment',
                      'AB39:Treatment'),
           new_y = if_else(CCA2 < 0, CCA2 - 0.24, CCA2 + 0.24))

# Plot sites (plants) with vectors
figure2_ordination <- ggplot() + 
    geom_point(data = line_treatment_cca, 
               aes(x = CCA1, y = CCA2, color = treatment, shape = line)) +
    geom_polygon(data = ells, aes(x = Values, y = y, group = Humid), 
                 alpha = 0.1) +
    geom_segment(data = vectors, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), 
                 arrow = arrow(length = unit(0.03, 'npc')),
                 color = 'blue') +
    geom_text(data = vectors, aes(x = CCA1, y = new_y, label = labels),
              color = 'blue') + 
    theme_bw() +
    scale_color_grey(name = 'Treatment', 
                     labels = c('Inoculated', 'Not Inoculated')) + 
    scale_shape_manual(name = 'Line',
                       values = c(15, 17, 19)) +
    theme(legend.position = c(0.75, 0.12),
          legend.background = element_rect(color = 'grey30'),
          legend.direction = 'horizontal',
          panel.border = element_blank(),
          axis.line.x = element_line(size = 0.48, linetype = 'solid', 
                                     color = 'black'),
          axis.line.y = element_line(size = 0.48, linetype = 'solid', 
                                     color = 'black'),
          legend.margin = ggplot2::margin(2, 2, 2, 2)) + 
    labs(x = 'CCA1 (45.3%)', y = 'CCA2 (26.6%)')

ggsave(plot = figure2_ordination, filename = './figures/raw-figures/fig2-ordination.pdf',
       width = 8, height = 5)

# Heirarchical Cluster Analysis -------------------

# reformat for clustering of metabolites
to_clust <- metabolite_treatment %>% 
    filter(!is.na(revised_names)) %>% 
    group_by(revised_names, local_sample_id) %>% 
    summarize(value = median(value)) %>% 
    spread(revised_names, value) 

set.seed(24)
# Run cluster analysis with bootstrapped p-values
cl <- pvclust(to_clust[,-1], parallel = TRUE)

# Grab the dendrogram
dend <- as.dendrogram(cl$hclust)

plot(dend)

# We'll highlight 8 groups
num_clades <- 8

# We want to color the labels and groups
dend <- dend %>% 
    color_branches(k=num_clades, col=brewer.pal(8, 'Dark2')) %>% 
    color_labels(k=num_clades, col=brewer.pal(8, 'Dark2'))

# Extract dendrogram data for plotting
plot_dend <- dendro_data(dend)

# Plot dendrogram
fig3_dendrogram <- ggplot(segment(plot_dend)) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    coord_flip() + 
    theme_dendro()

ggsave(plot = fig3_dendrogram, 
       filename = './figures/raw-figures/fig3-dendrogram.pdf',
       width = 3, height = 12)

# Indicator Species Analysis/Multilevel pattern analysis ---------------

# Select named compounds for examination
named_compounds <- metabolite_treatment %>% 
    filter(!is.na(revised_names)) %>%
    group_by(local_sample_id, line, treatment, revised_names) %>% 
    summarize(value = mean(value)) %>%
    ungroup() %>%
    spread(revised_names, value)

# Separate metabolites and labels for analysis
metabolites <- named_compounds %>% select(-line, -treatment, -local_sample_id)
metabolite_groups <- named_compounds %>% select(line, treatment)

# Examine associations with treatment groups using Phi association index 
# Tichy and Chytry 2006
treatment_indicators <- multipatt(metabolites, metabolite_groups$treatment,
                                  func = 'r.g',
                                  control = how(nperm = 999))

# Examine association
summary(treatment_indicators, indvalcomp = TRUE)

# Examine associations with grass lines
line_indicators <- multipatt(metabolites, metabolite_groups$line,
                             func = 'r.g', restcomb = c(1,2,3),
                             control = how(nperm = 999))

# Examine associations
summary(line_indicators)

# Extract assaciation values for visualization
line_indicator_values <- line_indicators$str %>% 
    data.frame() %>%
    rownames_to_column('revised_names') %>%
    mutate(revised_names = fct_relevel(revised_names, labels(dend))) %>%
    gather(line, rg, -revised_names)

# Plot heatmap of associations
plot_heatmap <- function(legend = FALSE) {
    ggplot() + 
    geom_tile(data = line_indicator_values, aes(x = line, y = revised_names, 
                                         fill = rg)) + 
    scale_fill_gradient(name = expression(phi), 
                        low = 'blue', high = 'red',
                        breaks = seq(-0.75, 0.75, by = 0.25)) + 
    geom_tile(colour="white",size=0.25)+
    labs(x="",y="")+
    scale_y_discrete(expand=c(0,0))+
    coord_fixed() +
    theme(plot.background = element_blank(),
          panel.border = element_blank(),
          legend.position = if_else(legend, 'right', 'none'),
          legend.text = element_text(size = 10, angle = 90, hjust = 1),
          legend.title = element_text(size = 18),
          axis.text.x = element_text(angle = 68, hjust = 0.6,
                                     vjust = 0.5),
          text = element_text(size = 8))
}

# Plot heatmap without legend
fig3_associations <- plot_heatmap()

# Grab legend for editing later
fig3_legend <- cowplot::get_legend(plot_heatmap(legend = TRUE))

# Save without legend
ggsave(plot = fig3_associations, 
       filename = './figures/raw-figures/fig3-associations.pdf',
       width = 4, height = 12)
    
# Save legend separately
cowplot::save_plot(plot = fig3_legend, 
                   filename = './figures/raw-figures/fig3-legend.pdf')

# Examine treatment differences by line -----------------

# Set up model function for wilcoxon test
model_differences <- function(df) {
    wilcox.test(value ~ treatment, data = df, conf.int = TRUE)
}

# Set up function to extract model statistics
extract_stats <- function(wilcox_mod) {
    data.frame(estimate = wilcox_mod$estimate,
               p_value = wilcox_mod$p.value)
}

# Testing for differences...
multiple_tests <- to_clust %>%
    # separate out line, treatment, and rep
    separate(local_sample_id, c('line', 'treatment', 'rep')) %>%
    # do not need rep
    select(-rep) %>%
    ungroup() %>%
    # collapse into long form
    gather(key = 'revised_names', value = 'value', -line, -treatment) %>%
    # regroup for modeling
    group_by(revised_names, line) %>%
    # nested data frames are awesome!!!
    nest() %>%
    # run models and extract stats
    mutate(wilcox_model = purrr::map(data, model_differences),
           stats = purrr::map(wilcox_model, extract_stats)) %>%
    # clean up and unnest
    select(revised_names, line, stats) %>%
    unnest() %>%
    # correct for false discovery rate
   mutate(adj_p = p.adjust(p_value, method = 'BH'))

# Check distribution of location parameter 
# (median of difference, not difference of median)
ggplot(multiple_tests, aes(x = estimate)) + stat_density()

sig <- multiple_tests %>%
    filter(adj_p <= 0.05) %>% 
    arrange(adj_p)

data.frame(sig)

# Examine differences in amino acids -------------------

# Amino acids identified in analysis
amino_acids <- c('L-Tyrosine', 
                 'Threonine/Homoserine',
                 'Leucene',
                 'Citrulline',
                 'Phenylalanine',
                 'Alanine',
                 'L-Glutamic Acid',
                 'Tryptophan',
                 'N-Methyl-D-Aspartic Acid',
                 'Alanine/Sarcosine',
                 'Phenylalanine-HCOOH',
                 'L-Methionine', 
                 'L-Proline',
                 'L-Isoleucine',
                 'N-Alpha-Acetyl-L-Lysine')

# Order amino acids by significance of tests for AB03
amino_acid_order <- multiple_tests %>% 
    filter(revised_names %in% amino_acids,
           line == 'AB03') %>% 
    mutate(revised_names = as.character(revised_names)) %>% 
    arrange(line, adj_p) %>% 
    pull(revised_names)

# For displaying differences, we'll want to normalize for comparison across 
# compounds and lines for uninoculated.  To do that, we'll calculate 
# the mean abundance by compound and line, then subtract it from the raw values
normalization_mean <- metabolite_treatment %>% 
    # filter only amino acids
    filter(revised_names %in% amino_acids) %>%
    # group over positive/negation ions values
    group_by(revised_names, line, treatment, plant_number) %>%
    summarize(value = mean(value)) %>%
    ungroup() %>%
    # only uninoculated treatments for normalization 
    filter(treatment == 'U') %>% 
    # group by line and compound for normalization
    group_by(line, revised_names) %>% 
    # find mean
    summarize(mn = mean(value)) 

# normalize abundances with previously calculated means
normalized_abundances <- metabolite_treatment %>% 
    filter(revised_names %in% amino_acids) %>%
    group_by(revised_names, line, treatment, plant_number) %>%
    summarize(value = mean(value)) %>%
    ungroup() %>%
    # join back previously calculated means
    left_join(normalization_mean) %>% 
    mutate(normalized_value = value - mn,
           revised_names = fct_relevel(revised_names, amino_acid_order))

# determine significances to display in plot
amino_acid_significance <- multiple_tests %>% 
    filter(revised_names %in% amino_acids,
           adj_p <= 0.05) %>% 
    mutate(sig = '*',
           location = 1.5) %>%
    select(revised_names, line, sig, location)

dodge_width = 0.48

# plot amino acid abundance by line
fig4_amino_acids <- ggplot() + 
    geom_point(data = normalized_abundances, 
               aes(x = revised_names, y = normalized_value, 
                   color = treatment),
               position = position_jitterdodge(jitter.width = dodge_width/2,
                                               dodge.width = dodge_width),
               alpha = 0.48) + 
    stat_summary(data = normalized_abundances, 
                 aes(x = revised_names, y = normalized_value, 
                     color = treatment),
                 fun.data = 'mean_cl_boot', geom = 'point',
                 position = position_dodge(dodge_width)) +
    stat_summary(data = normalized_abundances, 
                 aes(x = revised_names, y = normalized_value, 
                     color = treatment),
                 fun.data = 'mean_cl_boot', geom = 'errorbar',
                 position = position_dodge(dodge_width),
                 width = dodge_width) + 
    geom_text(data = amino_acid_significance,
              aes(x = revised_names, y = location,
                  label = '*')) +
    facet_wrap(line ~ ., ncol = 1) + 
    theme_bw() + 
    scale_color_grey(labels = c('Inoculated', 'Not Inoculated')) + 
    scale_y_continuous(limits = c(-3,3)) + 
    labs(x = '', y = 'Normalized Abundance') + 
    theme(text = element_text(size = 10),
          axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          legend.position = 'right',
          legend.background = element_rect(color = 'grey30'),
          legend.title = element_blank())

ggsave(plot = fig4_amino_acids, 
       filename = './figures/raw-figures/fig4-amino-acids.pdf',
       width = 7, height = 8)

# Nematode Counts + Metabolites ------------------

# select compounds of interest based on:
# preliminary screening of differences related to AB39 
# and knowledge of plant defense pathways
compounds_of_interest <- c('D-Glucuronic Acid', 
                           'Phenylalanine', 
                           'Glycolate')

# Metabolite abundances for compound of interest
# only for inoculated plants with nematodes
coi_abundance <- metabolite_treatment %>% 
    filter(revised_names %in% compounds_of_interest,
           treatment == 'I') %>% 
        group_by(revised_names, line, plant_number) %>% 
    summarize(value = mean(value))

# Prepare nematode counts for joining with metabolite data
nematodes <- nematode_counts %>% 
    filter(treatment !='uninoc') %>% 
    mutate(plant_number = as.character(plant_number))

# Normalize nematode counts based on a per line basis
nem_norm <- nematodes %>% 
    group_by(line) %>% 
    mutate(norm_count = nematode_count - mean(nematode_count)) %>% 
    select(line, plant_number, norm_count)

# Normalize metabolite levels by line and compound
metab_norm <- coi_abundance %>% 
    group_by(revised_names, line) %>% 
    mutate(norm_value = value -mean(value))

# Join normalized nematode and metabolite data
nematode_metabolites <- left_join(metab_norm, nem_norm) %>% 
    drop_na()

# Plot nematode metabolite relationship
fig5_compounds_of_interest <- ggplot(nematode_metabolites, aes(x = norm_value, y = norm_count)) + 
    geom_point() + 
    facet_wrap( ~ revised_names, ncol = 2, scales = 'free') + 
    stat_smooth(method = 'lm') + 
    theme_bw() + 
    labs(x = 'Normalized Abundance', y = 'Normalized Nematode Count') + 
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.48, linetype = 'solid', 
                                     color = 'black'),
          axis.line.y = element_line(size = 0.48, linetype = 'solid', 
                                     color = 'black'))

ggsave(plot = fig5_compounds_of_interest, 
       filename = './figures/raw-figures/fig5_compounds_of_interest.pdf',
       width = 12, height = 7)

# Helper function to model counts
lm_mod <- function(df) {
    lm(norm_count ~ norm_value, data = df)
}

# Determine model statistics for compounds of interest
nematode_metabolites %>% 
    group_by(revised_names) %>% 
    nest() %>% 
    mutate(lm_mod = purrr::map(data, lm_mod),
           glance = purrr::map(lm_mod, broom::glance)) %>% 
    select(revised_names, glance) %>% 
    unnest()

# Examine Phenylalanine more closely
phenylalanine <- left_join(coi_abundance, nematodes) %>% 
    ungroup() %>% drop_na() %>% 
    filter(revised_names == 'Phenylalanine') %>% 
    select(line, plant_number, value, nematode_count)

# Relate observed levels of phenylalanine to 
# nematode counts
phenylalanine_model <- lm(nematode_count ~ line + value, 
                          data = phenylalanine)

# double check line-nematode count relationship
emmeans(phenylalanine_model, ~ line) %>% contrast() %>% cld()

# examine model diagnostics
summary(phenylalanine_model)
Anova(phenylalanine_model, type = 'III')
#plot(phenylalanine_model, which = 1:6)
qqPlot(phenylalanine_model)
lrtest(lm(nematode_count ~ 1, data = phenylalanine),
       phenylalanine_model)

# Create sequence from ranges
expand_inputs <- function(df) {
    seq(df$min, df$max, by = 0.01)
}

# Determine ranges of phenylalanine by line
# then create sequences used for plotting
model_inputs <- phenylalanine %>% 
    group_by(line) %>% 
    summarize(min = min(value),
              max = max(value)) %>% 
    group_by(line) %>% 
    nest() %>% 
    mutate(value = purrr::map(data, expand_inputs)) %>% 
    select(line, value) %>% 
    unnest()

# Generate predictions and confidence intervals
# to display results of the model
phenylalanine_model_predictions <- predict(phenylalanine_model, 
                                           newdata = model_inputs, 
                                           interval = 'confidence') %>% 
    data.frame()

# Combine model predictions with inputs
palanine_model_results <- bind_cols(model_inputs, 
                                    phenylalanine_model_predictions)

# Display model results for phenylalanine production by line
fig5_phenylalanine <- ggplot(palanine_model_results, aes(x = value, y = fit, 
                                   fill = line,
                                   linetype = line)) + 
    geom_line() + 
    geom_ribbon(aes(x = value, ymin = lwr, ymax = upr), alpha = 0.5) +
    scale_fill_grey(name = 'Line') + 
    scale_linetype_discrete(name = 'Line') + 
    theme_bw() + 
    labs(x = 'Phenylalanine Abundance', y = 'Nematode Count') + 
    theme(legend.text = element_text(size = 10),
          legend.position = c(0.9, 0.25), 
          legend.direction = 'vertical',
          legend.background = element_rect(color = 'grey30'),
          panel.border = element_blank(),
          axis.line.x = element_line(size = 0.48, linetype = 'solid', 
                                     color = 'black'),
          axis.line.y = element_line(size = 0.48, linetype = 'solid', 
                                     color = 'black'),
          legend.margin = ggplot2::margin(1, 1, 1, 1)) 

ggsave(plot = fig5_phenylalanine, 
       filename = './figures/raw-figures/fig5_phenylalanine.pdf',
       width = 5, height = 3)

# Pipecolic Acid Analysis ----------------------------

# Prep nematode data
nematode_bootstrap <- nematodes %>% 
    select(line, nematode_count) %>% 
    group_by(line) %>% 
    summarize(mean = smean.cl.boot(nematode_count)[1],
              lower = smean.cl.boot(nematode_count)[2],
              upper = smean.cl.boot(nematode_count)[3])

# Pull metabolite data on Pipecolic Acid
pipecolic_acid <- metabolite_treatment %>% 
    filter(revised_names == 'Pipecolate/L-Pipecolic Acid')

# Visualize differences it Pipecolic Acid Production
set.seed(72)

dodge_width = 0.48
jitter_width = 0.12

fig6_pipecolic_acid <- ggplot(pipecolic_acid, 
       aes(x = line, y = value, color = treatment)) + 
    geom_point(position = position_jitterdodge(jitter.width = jitter_width,
                                               dodge.width = dodge_width),
               alpha = 0.48) + 
    stat_summary(fun.data = 'mean_cl_boot', 
                 geom = 'point',
                 position = position_dodge(dodge_width)) + 
    stat_summary(fun.data = 'mean_cl_boot', 
                 geom = 'errorbar',
                 position = position_dodge(dodge_width)) + 
    theme_bw() + 
    scale_color_grey(name = 'Treatment', 
                     labels = c('Inoculated', 'Not Inoculated')) +
    labs(x = 'Line', y = 'L-Pipecolic Acid Abundance') + 
    theme(legend.position = 'bottom',
          legend.background = element_rect(color = 'grey30'),
          panel.border = element_blank(),
          axis.line.x = element_line(size = 0.48, linetype = 'solid', 
                                     color = 'black'),
          axis.line.y = element_line(size = 0.48, linetype = 'solid', 
                                     color = 'black'),
          legend.margin = ggplot2::margin(1, 1, 1, 1)) 

ggsave(plot = fig6_pipecolic_acid, 
       filename = './figures/raw-figures/fig6_pipecolic_acid.pdf',
       width = 5, height = 3.48)

# Test results
t_test <- function(df) {
    ttest <- t.test(value ~ treatment, data = df)
    data.frame(mean = ttest$estimate[1] - ttest$estimate[2],
               lower = ttest$conf.int[1],
               upper = ttest$conf.int[2],
               t_statistic = ttest$statistic,
               df = ttest$parameter,
               pvalue = ttest$p.value)
}

# Look at differences in pipecolic acid production
pipecolic_acid %>% group_by(line) %>% 
    nest() %>%
    mutate(ttest = purrr::map(data, t_test)) %>% 
    select(line, ttest) %>% 
    unnest() %>% 
    mutate(p_adj = p.adjust(pvalue, method = 'bonferroni'))


# Bootstrap differences in pipecolic acid production
B = 1000

bootstrapped_means <- tibble()

for (i in 1:B) {
    # Calculate and store means by line and treatment
    calculate_means <- pipecolic_acid %>% 
        group_by(line, treatment) %>% 
        sample_frac(size = 1, replace = TRUE) %>% 
        summarize(mean = mean(value)) %>% 
        spread(treatment, mean) %>% 
        ungroup() %>% 
        mutate(bootstrap_rep = i)
    
    bootstrapped_means <- bind_rows(bootstrapped_means, calculate_means)

}

# Look and differences for display
pipecolic_summary <- bootstrapped_means %>% 
    mutate(difference = I - U) %>%
    group_by(line) %>% 
    summarize(mean = mean(difference),
              lower = quantile(difference, 0.025),
              upper = quantile(difference, 0.975))

joined_pipecolic <- left_join(pipecolic_summary, 
                              nematode_bootstrap, 
                              by = 'line',
                              suffix = c('_pipecolic', '_nematode'))

# Look at differences associated with nematode production
fig6_pipecolic_differences <- ggplot(joined_pipecolic, 
       aes(x = mean_pipecolic, 
           y = mean_nematode,
           color = line,
           shape = line)) + 
    geom_point() + 
    geom_errorbar(aes(x = mean_pipecolic, 
                      ymin = lower_nematode, 
                      ymax = upper_nematode),
                  width = 0.12) + 
    geom_errorbarh(aes(y = mean_nematode, 
                       xmin = lower_pipecolic, 
                       xmax = upper_pipecolic),
                   height = 10) + 
    geom_vline(xintercept = 0, 
               color = 'black', 
               alpha = 0.72, 
               linetype = 'dashed') +
    scale_color_grey(name = 'Line') + 
    scale_shape_manual(name = 'Line',
                       values = c(15, 17, 19)) + 
    labs(x = 'Increase in Pipecolic Acid', y = 'Sting Nematode Count') + 
    theme_bw() + 
    theme(legend.position = 'bottom',
          legend.background = element_rect(color = 'grey30'),
          panel.border = element_blank(),
          axis.line.x = element_line(size = 0.48, linetype = 'solid', 
                                     color = 'black'),
          axis.line.y = element_line(size = 0.48, linetype = 'solid', 
                                     color = 'black'),
          legend.margin = ggplot2::margin(1, 1, 1, 1)) 

ggsave(plot = fig6_pipecolic_differences, 
       filename = './figures/raw-figures/fig6_pipecolic_differences.pdf',
       width = 5, height = 3.48)