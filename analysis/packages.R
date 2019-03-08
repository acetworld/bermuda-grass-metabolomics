

if (!require("devtools")) install.packages("devtools")
if (!require("pacman")) install.packages("pacman")
if (!require('MASS')) pacman::p_install('MASS')
if (!require('Hmisc')) pacman::p_install('Hmisc')
pacman::p_load('tidyverse', 'stringr', 'foreach', 'testthat', 'VIM', 'vsn',
               'pvclust', 'dendextend', 'circlize', 'lettercase', 
               'RColorBrewer', 'kohonen', 'caret', 'vegan', 'lmPerm',
               'randomForest', 'doMC', 'forcats', 'plotROC','car','emmeans',
               'multcompView', 'lmtest', 'vegan', 'candisc', 'indicspecies')
