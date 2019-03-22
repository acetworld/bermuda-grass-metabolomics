source('./analysis/preprocess_metabolite_data.R')

# Prepare ----------------------------

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
ggplot() + 
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
                     labels = c('Inoculated', 'Uninoculated')) + 
    scale_shape_manual(name = 'Line',
                       values = c(15, 17, 19)) +
    theme(legend.position = c(0.75, 0.12),
          legend.background = element_rect(color = 'grey30'),
          legend.direction = 'horizontal') + 
    labs(x = 'CCA1 (45.3%)', y = 'CCA2 (26.6%)')


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
ggplot(segment(plot_dend)) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    coord_flip() + 
    theme_dendro()

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
ggplot() + 
    geom_tile(data = line_indicator_values, aes(x = line, y = revised_names, 
                                         fill = rg)) + 
    scale_fill_gradient(low = 'blue', high = 'red') + 
    geom_tile(colour="white",size=0.25)+
    labs(x="",y="")+
    scale_y_discrete(expand=c(0,0))+
    coord_fixed() +
    theme(plot.background = element_blank(),
          panel.border = element_blank())

    
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
           location = 2) %>%
    select(revised_names, line, sig, location)

dodge_width = 0.48

# plot amino acid abundance by line
ggplot() + 
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
    facet_wrap(~line, ncol = 1) + 
    theme_bw() + 
    scale_color_grey(labels = c('Inoculated', 'Uninoculated')) + 
    scale_y_continuous(limits = c(-3,3)) + 
    labs(x = '', y = 'Normalized Abundance') + 
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          legend.position = c(0.75, 0.12),
          legend.background = element_rect(color = 'grey30'),
          legend.title = element_blank())

# Nematode Counts + Metabolites ------------------

# select compounds of interest based on:
# preliminary screening of differences related to AB39 
# and knowledge of plant defense pathways
compounds_of_interest <- c('D-Glucaronic Acid', 
                           'Phenylalanine', 
                           'Glycolate',
                           'Phenylalanine-HCOOH')

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
ggplot(nematode_metabolites, aes(x = norm_value, y = norm_count)) + 
    geom_point() + 
    facet_wrap( ~ revised_names, ncol = 2, scales = 'free') + 
    stat_smooth(method = 'lm') + 
    theme_bw()

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

Phenylalanine Chart...
pala <- left_join(coi_abundance, nematodes) %>%
    filter(revised_names == 'Phenylalanine')

ggplot(pala, aes(x = line, y = value)) + 
    geom_point()

ggplot(pala, aes(x = line, y = nematode_count)) + 
    geom_point()


b_rep = 1000

hold = tibble()

for(i in 1:b_rep) {
    means <- pala %>% 
        group_by(line) %>%
        sample_frac(size = 1, replace = TRUE) %>% 
        summarize(mn_nematode = mean(nematode_count, na.rm = TRUE),
                  mn_value = mean(value, na.rm = TRUE))
    
    hold <- bind_rows(hold, means)
    
    
}

test <- hold %>% 
    group_by(line) %>% 
    summarize(mn_nem = mean(mn_nematode),
              low_nematode = quantile(mn_nematode, 0.025),
              high_nematode = quantile(mn_nematode, 0.975),
              mn_val = mean(mn_value),
              low_val = quantile(mn_value, 0.025),
              high_val = quantile(mn_value, 0.975))

ggplot(test, aes(x = mn_val, y = mn_nem)) + 
    geom_point() + 
    geom_errorbar(aes(x = mn_val, ymin = low_nematode, ymax = high_nematode)) + 
    geom_errorbarh(aes(xmin = low_val, xmax = high_val, y = mn_nem))

ggplot(pala, aes(x = value, y = nematode_count, color = line)) + 
    geom_point() + 
    theme_bw(14) + 
    stat_smooth(method =  'lm') +
    scale_color_grey()

ggplot(pala, aes(x = value, y = nematode_count)) + 
    stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar') 


test_mod <- lm(nematode_count ~ line + value, data = pala)

Anova(test_mod, type = 'III')

summary(test_mod)
plot(test_mod, which = 1:6)
qqPlot(test_mod$residuals)

conts <- emmeans(test_mod, ~ line + value)

conts %>% cld()
emtrends(conts, var = 'value')

emtrends(test_mod, ~ line + value, var = 'value')

# ------------------------


h

mean_test <- function(df) {
    t.test(value ~ treatment, data = df)
}

mean_estimates <- function(test) {
    data.frame(estimate = t_test$estimate[1] - t_test$estimate[2],
               lower = t_test$conf.int[1],
               upper = t_test$conf.int[2])
}

mean_diffs <- metabolite_treatment %>% 
    filter(revised_names %in% 
               c('L-Leucine',
                 'L-Tyrosine', 'Threonine/Homoserine',
                 'Leucene', 'Phenylalanine',
                 'Alanine','L-Glutamic Acid',
                 'Tryptophan-2,3,3-D3',
                 'N-Methyl-D-Aspartic Acid',
                 'Alanine/Sarcosine',
                 'Phenylalanine-HCOOH',
                 'L-Methionine', 'L-Proline',
                 'Tryptophan')) %>%
    group_by(revised_names, line, treatment, plant_number) %>%
    summarize(value = mean(value)) %>%
    ungroup() %>%
    group_by(line, revised_names) %>% 
    mutate(mn = mean(value)) %>% 
    filter(treatment == 'U') %>% 
    select(revised_names, line, mn)

diffs <- metabolite_treatment %>% 
    filter(revised_names %in% 
               c('L-Leucine',
                 'L-Tyrosine', 'Threonine/Homoserine',
                 'Leucene', 'Phenylalanine',
                 'Alanine','L-Glutamic Acid',
                 'Tryptophan-2,3,3-D3',
                 'N-Methyl-D-Aspartic Acid',
                 'Alanine/Sarcosine',
                 'Phenylalanine-HCOOH',
                 'L-Methionine', 'L-Proline',
                 'Tryptophan')) %>%
    group_by(revised_names, line, treatment, plant_number) %>%
    summarize(value = mean(value)) %>%
    ungroup() %>%
    left_join(mean_diffs) %>% 
    mutate(normalized_value = value - mn)
    
    
    
    
    select(-plant_number) %>%
    group_by(revised_names, line) %>%
    nest() %>%
    mutate(t_test = purrr::map(data, mean_test),
           estimates = purrr::map(t_test, mean_estimates)) %>% 
    select(revised_names, line, estimates) %>% 
    unnest()

ggplot(diffs, aes(x = revised_names, y = normalized_value, color = treatment)) + 
    geom_point(position = position_dodge(0.24)) + 
    stat_summary(fun.data = 'mean_cl_boot', geom = 'point',
                 position = position_dodge(0.24)) +
    stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar',
                 position = position_dodge(0.24)) + 
    facet_wrap(~line, ncol = 1) + 
    theme_bw() + 
    scale_color_grey()




#---------------








ggplot(diffs, aes(x = revised_names, y = estimate, shape = line)) + 
    geom_point(position = position_dodge(0.24)) + 
    geom_errorbar(aes(x = revised_names, ymin = lower, ymax = upper),
                  position = position_dodge(0.24)) + 
    theme_bw() +
    coord_flip()
    


ggplot(diffs, aes(x = line, y = value, color = treatment)) + 
    geom_point(position = position_dodge(0.24)) + 
    theme_bw() + 
    facet_wrap(~ revised_names, ncol = 3) +
    scale_color_grey()



ggplot() + 
    geom_tile(data = multiple_tests, aes(x = line, y = compounds, 
                                         fill = estimate)) + 
    scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white') + 
    #add border white colour of line thickness 0.25
    geom_tile(colour="white",size=0.25)+
    #remove x and y axis labels
    labs(x="",y="")+
    #remove extra space
    scale_y_discrete(expand=c(0,0))+
    coord_fixed() +
    theme(plot.background = element_blank(),
          panel.border = element_blank())




-------------------------------

vectors <- data.frame(ccaone$CCA$biplot * 
                          vecscale(ccaone$CCA$biplot, 
                                   bbox = matrix(c(-5, 3, -5, 5), 2, 2)))
vectors$Name <- rownames(vectors)


waplot <- cbind(metab_cause, ccaone$CCA$wa) 

eplot <- cbind(grps, ccaone$CCA$wa)

out <- dataEllipse(eplot$CCA1, eplot$CCA2, groups = eplot$Group, 
                   levels = 0.95) %>% data.frame() %>%
    gather(H, Values) %>%
    mutate(Coord = str_sub(H, -1, -1),
           Humid = str_sub(H, 1, -3)) %>%
    select(-H) %>% mutate(Coord = as.factor(Coord)) 

dataEllipse(vis_nmds$X, vis_nmds$Y, 
            groups = as.factor(vis_nmds$line_treatment), 
            levels = 0.95) %>% data.frame()
# Nonmetric Multidimensional Scaling ------------------------------

# reformat for distance matrix calculations
to_nmds <- metabolite_treatment %>% 
    unite(line_treatment, line, treatment) %>%
    unite(ion_rt, ion_type, compounds) %>% 
    spread(ion_rt, value)

# calculate distances between samples
nema_dist <- dist(to_nmds[, 4:2694])

# non-metric multidimensional scaling
nema_nmds <- MASS::isoMDS(nema_dist, k = 2)

# add treatment labels back
vis_nmds <- data.frame(line_treatment = to_nmds$line_treatment,
                       X = nema_nmds$points[,1],
                       Y = nema_nmds$points[,2]) %>% 
    separate(line_treatment, into = c('line', 'treatment'), remove = FALSE)

# plot on scaling axes
nmds_plot <- ggplot(vis_nmds, aes(x = X, y = Y, 
                                  color = treatment, shape = line)) + 
    geom_point() + 
    labs(x = 'Scaling Axis 1', y = 'Scaling Axis 2') + 
    scale_color_grey(name = 'Treatment',
                     labels = c('Nematodes', 'No Nematodes')) +
    scale_shape_manual(name = 'Genotype',
                       values = c(16, 17, 15)) 

ggsave(plot = nmds_plot, filename = 'NMDS_Plot.pdf',
       width = 8, height = 5)

# Heirarchical Cluster Analysis -------------------

# reformat for clustering of metabolites
to_clust <- metabolite_treatment %>% 
    filter(!is.na(revised_names)) %>% 
    group_by(revised_names, local_sample_id) %>% 
    summarize(value = median(value)) %>% 
    spread(revised_names, value) 

# Run cluster analysis with bootstrapped p-values
cl <- pvclust(to_clust[,-1], parallel = TRUE)

# Grab the dendrogram
dend <- as.dendrogram(cl$hclust)

# We'll highlight 8 groups
num_clades <- 8

# We want to color the labels and groups
dend <- dend %>% 
    color_branches(k=num_clades, col=brewer.pal(8, 'Dark2')) %>% 
    color_labels(k=num_clades, col=brewer.pal(8, 'Dark2'))

# Plot
pdf('dendrogram.pdf', width = 12, height = 12)
par(cex = 0.5, mar = rep(0,4))
circlize_dendrogram(dend, dend_track_height = 0.4,
                    labels_track_height = 0.5) 
dev.off()

# Random Forests ------------------------

# register parallel backend
registerDoMC(cores = 4)

# All compounds
u <- processed_nema %>% unite(c_i, compounds, ion_type) 
ut <- base::t(u[, -1]) %>% data.frame() 
colnames(ut) <- u$c_i
ut$local_sample_id <- names(processed_nema)[3:43]

# Join with treatment data to get infection status
to_rf <- left_join(ut, metabolite_treatment) %>% 
    mutate(infection = `Nematode inoculation`) %>% 
    select(-local_sample_id, -mb_sample_id, - `Nematode inoculation`)

to_rf <- metabolite_treatment %>%
    filter(!is.na(revised_names)) %>%
    group_by(local_sample_id, revised_names) %>%
    summarize(value = mean(value)) %>%
    ungroup() %>% 
    left_join(nematode_counts, by = c('local_sample_id' = 'identifier')) %>%
    filter(treatment == 'inoc') %>%
    select(revised_names, line, nematode_count, value, local_sample_id) %>%
    spread(revised_names, value) %>%
    select(-local_sample_id) %>%
    mutate(line = as.factor(line))


to_rf <- metabolite_treatment %>%
    group_by(local_sample_id, compounds) %>%
    summarize(value = mean(value)) %>%
    ungroup() %>% 
    left_join(nematode_counts, by = c('local_sample_id' = 'identifier')) %>%
    filter(treatment == 'inoc') %>%
    select(local_sample_id, line, nematode_count, compounds, value) %>% 
    spread(compounds, value)

# Split into training and test set
train_index <- createDataPartition(to_rf$nematode_count, 
                                       p = 0.7, list = FALSE)

train_data <- to_rf[train_index,]
test_data <- to_rf[-train_index,]

# Repeated cross validation
rf_control <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           savePredictions = T)

# Run RF model
rf_model <- train(nematode_count ~ .,
                  data = to_rf,
                  trControl = rf_control,
                  method = 'rf',
                  importance = TRUE,
                  verbose = TRUE)

rf_model
varImpPlot(rf_model$finalModel)

plot(predict(rf_model, to_rf), to_rf$nematode_count)

preds <- predict(rf_model, test_data)

plot(preds, test_data$nematode_count)

# Just named compounds
named <- left_join(to_clust, treatment_data) %>% 
    mutate(infection = `Nematode inoculation`) %>% 
    select(-local_sample_id, -mb_sample_id, - `Nematode inoculation`)

# set random seed for reproducibility
set.seed(24)

# Split into training and test set
train_index <- createDataPartition(to_rf$infection, 
                                   p = 0.7, list = FALSE)
train_index_named <- createDataPartition(named$infection, 
                                         p = 0.7, list = FALSE)

train_data <- to_rf[train_index,]
test_data <- to_rf[-train_index,]

train_data_named <- to_rf[train_index_named,]
test_data_named <- to_rf[-train_index_named,]

# Check for imbalances
table(train_data$infection)/nrow(train_data)
table(test_data$infection)/nrow(test_data)

# Repeated cross validation
rf_control <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           summaryFunction = twoClassSummary, 
                           classProbs = T,
                           savePredictions = T)

# Run RF model
rf_model <- train(infection ~ .,
                  data = train_data,
                  trControl = rf_control,
                  method = 'rf',
                  importance = TRUE,
                  verbose = TRUE)

# Run RF model just on named compounds

rf_model_named <- train(infection ~ .,
                        data = train_data_named,
                        trControl = rf_control,
                        method = 'rf',
                        importance = TRUE,
                        verbose = TRUE)

# Variable importance measures for unnamed compounds
var_importance <- varImp(rf_model, scale = TRUE)$importance
var_importance$compound <- rownames(var_importance) 

# Take top 40 to plot
var_plot <- var_importance %>% top_n(40, inoc) %>% 
    arrange(desc(inoc)) %>% 
    mutate(compound = fct_reorder(str_sub(compound, 2, 9), inoc))

varimp_plot <- ggplot(var_plot, aes(x = compound, y = inoc)) + 
    geom_point() + 
    geom_segment(aes(x = compound, xend = compound,
                     y = 0, yend = inoc)) + 
    theme_bw(14) + 
    labs(x = 'Compound or Retention Time', y = 'Relative Variable Importance') +
    coord_flip(ylim = c(70, 102)) + 
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = 'solid',
                                     color = 'black'),
          axis.line.y = element_line(size = 0.5, linetype = 'solid',
                                     color = 'black'),
          panel.background = element_blank(),
          axis.text.y = element_text(size = 8))

ggsave(plot = varimp_plot, filename = 'VarImpPlot.pdf',
       width = 8, height = 5)

# Confusion matrix for all compounds
confusionMatrix(data = predict(rf_model, test_data), 
                reference = test_data$infection)

# Confusion matrix for named compounds
confusionMatrix(data = predict(rf_model_named, test_data_named), 
                reference = test_data_named$infection)


# ROC curves
si <- rf_model$pred$mtry == 2
si_named <- rf_model_named$pred$mtry == 2

roc_all <- cbind(rf_model$pred[si, ], Metabolites = 'All')
roc_named <-  cbind(rf_model_named$pred[si_named, ], Metabolites = 'Named')
proc <- rbind(roc_all, roc_named)

auc <- data.frame(Metabolites = c('All', 'Named'),
                  AUC = round((calc_auc(p_init))$AUC, 2),
                  x = c(0.25, 0.25),
                  y = c(0.96, 0.55))

roc_plot <- ggplot() + 
    geom_roc(data = proc, 
             aes(m=inoc, d=factor(obs, levels = c("uninoc", "inoc")),
                 color = Metabolites), n.cuts=0) + 
    geom_text(data = auc, 
              aes(x = x, y = y, color = Metabolites, 
                  label = paste('AUC =', AUC)),
              show.legend = FALSE) +
    coord_equal() +
    style_roc() +
    theme(text = element_text(size = 14),
          panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = 'solid',
                                     color = 'black'),
          axis.line.y = element_line(size = 0.5, linetype = 'solid',
                                     color = 'black'),
          panel.background = element_blank(),
          legend.background = element_rect(color = 'grey30'),
          legend.position = c(0.75, 0.25))

ggsave(plot = roc_plot, filename = 'ROC_Plot.pdf',
       width = 5, height = 5)

# Inference -----------------------------------

inf_plot_prep <- to_rf %>% 
    select(infection, `Plant genotype`, `401.1206_positive`)

inf_fig <- ggplot(inf_plot_prep, aes(x = `Plant genotype`, y = `401.1206_positive`,
                                     color = infection)) + 
    geom_point(position = position_jitterdodge(jitter.width = 0.12,
                                               dodge.width = 0.48)) + 
    stat_summary(fun.data = 'mean_cl_boot', geom = 'point',
                 position = position_dodge(0.48),
                 shape = 17, size = 2) +
    stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar',
                 position = position_dodge(0.48),
                 width = 0.24) +
    scale_color_grey(name = 'Nematode', labels = c('Inoculated', 'Uninoculated')) +
    labs(y = 'Normalized Area of \nCompound at RT 401.1206 +ESI') +
    theme_bw(14) + 
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = 'solid',
                                     color = 'black'),
          axis.line.y = element_line(size = 0.5, linetype = 'solid',
                                     color = 'black'),
          panel.background = element_blank(),
          legend.background = element_rect(color = 'grey30'),
          legend.position = c(0.24, 0.12))

ggsave(plot = inf_fig,
       filename = 'Inference_Figure.pdf',
       width = 8, height = 5)

aovp(`401.1206_positive` ~ infection + `Plant genotype`,
     data = inf_plot_prep) %>% summary()

