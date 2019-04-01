source('./analysis/load_metabolite_data.R')

# Combine common compounds ------------------------------

# We want to combine results for compunds that have more than 
# one retention time.  To do this, we'll collapse all labeled compounds
metabolite_combined <- metabolite_data %>% 
    mutate(compounds = str_replace(Samples, '_R[12345678]', '')) %>% 
    select(-Samples) %>%
    group_by(compounds, ion_type) %>% 
    summarize_all(funs(sum)) %>%
    ungroup()


# Remove contaminants and adjust for internal standards ------------------------------

# We'll remove phthalates
# Also remove Dinoterb (applied herbicide) and ibuprofen

no_contaminants <- metabolite_combined %>% 
    filter(!(compounds %in% c('Dinoterb (herbicide)',
                              'Ibuprofen')),
           !(grepl('phthal', .$compounds, ignore.case = TRUE)))

# Remove internal standards referenced in procedures:
standard_blend <- tibble(
    compound = c('Tryptophan-2,3,3-D3', 
                 'L-Leucine-D10',
                 'Creatine-D3',
                 'CREATINE',
                 'Caffeine-D3',
                 'Caffeine 13C3',
                 'Salicylic Acid-D4',
                 'Succinic-2,2,3,3-D4 Acid',
                 'CITRIC ACID-13C6'),
    concentration_ug_ml = c(40, rep(4, 8))
)

clean_metabolites <- no_contaminants %>% 
    filter(!(compounds %in% standard_blend$compound))


# Imputation -----------------------------------

# We have missing data in our files 
clean_metabolites %>% gather(sample, value, -compounds, -ion_type) %>% 
    group_by(ion_type, sample) %>% 
    summarize(num_nas = sum(is.na(value)),
              perc_nas = num_nas/n()) %>% 
    arrange(perc_nas) %>%
    data.frame()

# We'll impute (Fill-in) the missing data using a K-nearest neighbors approach

metabolite_imputation <- kNN(clean_metabolites[,3:43], imp_var = FALSE, 
                   useImputedDist = FALSE)


# Normalize ----------------------------------

# For normalization we'll use variance stabilizing normalization
# which tends to have high performance (Li et al 2016) and has the nice
# benefit of tending to remove heteroscedasticity
es_nema <- ExpressionSet(as.matrix(metabolite_imputation))

vsn_model <- vsn2(es_nema)
meanSdPlot(vsn_model)

# Processed nematode grass metabolite data
processed_metabolites <- cbind(clean_metabolites[,1:2],
                        vsn_model@hx)

# Refined compound names for display
compound_names <- processed_metabolites %>%
    select(compounds) %>%
    mutate(numeric_names = as.numeric(compounds)) %>%
    filter(is.na(numeric_names)) %>%
    select(-numeric_names) %>% 
    mutate(revised_names = c('(R)-Malate',
                             'ACC',
                             '2,3-Dihydroxybenzoate',
                             '2,5-Dihydroxybenzoate',
                             '2,4-Dihydroxyacetophenone',
                             '3-Hydroxyphenyl Lactate',
                             '3-Hydroxyphenyl Acetate',
                             '3-Methyl Adipic Acid',
                             '3-Tert-Butyl Adipic Acid',
                             '3,4-Dihydroxybenzoate',
                             '3,4-Dihydroxyphenyl Acetate',
                             '4-Aminobutanoate',
                             '4-Guanidinobutanoate',
                             '4-Hydroxy-L-Phenylglycine',
                             'Oxoproline',
                             '5-Deoxyadenosine',
                             '6-Hydroxycaproic Acid',
                             '6-Phosphogluconic Acid',
                             '6C Sugar Alcohol',
                             'a-Ketoglutaric Acid',
                             'Adenine',
                             'ADP',
                             'Alanine/Sarcosine',
                             'Aldo/Keto-Hexose',
                             'Aldopentose',
                             'Allopurinol',
                             'a-Aminoadipate',
                             'a-Ketoglutaric Acid',
                             'Anthranilate',
                             'Arachidonic Acid',
                             'Arachidonic Acid',
                             'Asparagine',
                             'Aspartate',
                             'Aspartate',
                             'Benzoate',
                             'Betaine',
                             'C12-Disaccharide',
                             'C5 Sugar Alcohol',
                             'Hexose',
                             'Caffeate',
                             'Camphor',
                             'Citrate',
                             'Citrate',
                             'Citrulline',
                             'Cytosine',
                             'D-(+)-Raffinose',
                             'D-Glucuronic Acid',
                             'D-Saccharic Acid',
                             'Diethanolamine',
                             'Dipropylene Glycol',
                             'C6 Disaccharide',
                             'Erucamide',
                             'Erythritol',
                             'Eugenol/Isoeugenol',
                             'Ferulate',
                             'Fumaric Acid',
                             'Gluconic Acid',
                             'Glucosamine',
                             'Glucose/Fructose',
                             'Glucose/Fructose',
                             'Glutathione Disulfide',
                             'Glyceraldehyde',
                             'Glyceraldehyde/Lactate',
                             'Glyceric Acid',
                             'Glycerol',
                             'Glycolate',
                             'Guanine',
                             'Guanosine',
                             'Hexose Alcohol',
                             'Hexose-6-Phosphate',
                             'Hexose-Disaccharide',
                             'Homovanillate',
                             'Inosine',
                             'Isocitric Acid',
                             'Isocytosine',
                             'L-(+)-Lactic Acid', 
                             'L-2-Phosphoglyceric Acid',
                             'L-Carnitine',
                             'L-Cysteic Acid',                                              
                             'L-Glutamic Acid',                                             
                             'L-Glutamine',                                                
                             'L-Isoleucine',                                            
                             'L-Methionine',                                               
                             'L-Proline',                                                   
                             'L-Tyrosine',                                                  
                             'Leucine',
                             'Linoleic Acid',
                             '2-6-Diaminoheptanedioate',
                             'LysoPC',
                             'LysoPE',
                             'Malate',
                             'Malanate',
                             'Methyl-Beta-D-Galactoside',
                             'N-Acetylneuraminate',
                             'N-Alpha-Acetyl-L-Lysine',
                             'N-Butylbenzenesulfonamide',
                             'N-Methyl-D-Aspartic Acid',
                             'Nonanoate',
                             'Orthophosphate',
                             'Palmitoleic Acid',
                             'Pantothenic Acid',
                             'PEG N5',
                             'Phenylethylamine',
                             'Phenylalanine',
                             'Phenylalanine-HCOOH',
                             'Phosphocholine',
                             'Picolinic Acid',
                             'Pipecolate/L-Pipecolic Acid',
                             'Pyridoxal',
                             'Pyroglutamic Acid',
                             'Quinate',
                             'Quinate',
                             'Spermidine',
                             'Stachyose',
                             'Stachyose',
                             'Succinate',
                             'Succinate',
                             'Sulcatol',
                             'Tartaric Acid',
                             'Theophylline',
                             'Threonine/Homoserine',
                             'Trans-Cinnamate',
                             'Triethyl-Phosphate',
                             'Trigonelline',
                             'Tryptophan',
                             'Tryptophan-NH3',
                             'Urate',
                             'Uridine',
                             'Urocanate'))

