source('packages.R')

# Load Data --------------------------

# Text files with metabolomics workbench text files
nematode_files <- c('./data/ST000353_AN000575.txt', './data/ST000353_AN000576.txt')

# For these, the columns in the positive and negative ion sets are not the 
# same and not in order.  

read_mw_file <- function(file, ion_type) {
    # First read in the text and find where the MS data are
    all_data <- read_lines(file)
    start <- which(all_data == 'MS_METABOLITE_DATA_START')
    end <- which(all_data == 'MS_METABOLITE_DATA_END')
    
    # Grab the column names for future use
    column_names <- str_split(all_data[start + 1], 
                              '\t')[[1]] %>%
        str_replace_all(., ' ', '_')
    
    # Read in the MS data (no column headers)
    raw_data <- read_tsv(file, 
                         col_types = paste(c('c', 
                                             rep('d', length(column_names)-1)),
                                           collapse = ''),
                         na = c('\\N', 'NA', 'na'),
                         skip = start+2, 
                         n_max = end-start-3, 
                         col_names = column_names)
    
    # Label positive or negative ions
    table_data <- raw_data %>% mutate(ion_type = ion_type)
    
    # return output to foreach
    table_data
}

# Read in the files
positive <- read_mw_file(nematode_files[1], 'positive')
negative <- read_mw_file(nematode_files[2], 'negative')

# Double check to make sure we did not read in end of file marker
# This will throw an error if it is wrong
test_that('We did not read in end of file designator', {
    expect_length(which(positive$Samples == 'MS_METABOLITE_DATA_END'), 0)
    expect_length(which(negative$Samples == 'MS_METABOLITE_DATA_END'), 0)
})

# Let's combine the data by finding the common columns
common_columns <- intersect(colnames(positive), colnames(negative))

# Put the data together
nematode_data <- rbind(
    select_(positive, .dots = common_columns),
    select_(negative, .dots = common_columns)
)
