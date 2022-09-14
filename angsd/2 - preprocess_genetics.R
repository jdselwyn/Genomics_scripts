#### Libraries ####
library(tidyverse)
library(magrittr)
library(adegenet)
library(poppr)
library(SNPRelate)

#### Settings - and order of filtering ####
confidence_cutoff <- 0.99 #require probability of genotype to be greater than this
#Identify Clone Groups 
#Remove Clones based on which has least missing loci
individual_missingness <- 0.3 #Maximum % of loci missing in an individual to keep individual
#Split into two datasets one to compare A. cerv and A. palm & one to compare within A. cerv
#Remove if entirely missing within one species/location
#Remove if monomorphic
locus_missingness <- 0.1 #Maximum % of individuals missing in a locus to keep locus
minor_allele_frequency <- 0.05

#### Functions ####
process_angsd <- function(genotype_file, probability_file, sample_names, cutoff){
  filter_below <- function(genotype, probabilites, cutoff){
    to_swap_to_na <- probabilites < cutoff
    probabilites[to_swap_to_na] <- NA
    
    message('Swapped ', scales::percent(sum(to_swap_to_na, na.rm = TRUE) / prod(dim(probabilites)), accuracy = 1), 
            ' genotypes to NA due to having probability less than ', cutoff)
    
    genotype * round(probabilites, 0)
  }
  
  genotypes <- read_csv(genotype_file, 
                        col_types = cols(.default = col_integer(),
                                         contig = col_character(),
                                         ref = col_character(),
                                         alt = col_character())) %>%
    mutate(locus_number = row_number(),
           genotype = cbind(!!!syms(sample_names)),
           .keep = 'unused')
  
  genotype_probabilities <- read_csv(probability_file, 
                                     col_types = cols(.default = col_double(),
                                                      contig = col_character(),
                                                      ref = col_character(),
                                                      alt = col_character())) %>%
    mutate(locus_number = row_number(),
           genotype_probability = cbind(!!!syms(sample_names)),
           .keep = 'unused')
  
  full_join(genotypes, genotype_probabilities,
            by = c('contig', 'position', 'ref', 'alt', 'locus_number')) %>%
    
    #Filter out genotypes below confidence threshold
    mutate(genotype = filter_below(genotype = genotype, 
                                   probabilites = genotype_probability, 
                                   cutoff = cutoff), 
           .keep = 'unused') %>%
    
    #Remove Loci fixed after confidence cutoff
    mutate(alt_freq = rowSums(genotype, na.rm = TRUE) / (2 * (rowSums(!is.na(genotype)))),
           ref_freq = 1 - alt_freq, 
           .before = genotype) %>%
    filter(ref_freq != 1, alt_freq != 1) %>%
    arrange(locus_number) %>%
    select(-locus_number, -ref_freq, -alt_freq) %>%
    
    #Convert to genlight to more easily use in adegenet etc.
    mutate(snp_loc = str_c(ref, alt, sep = '/'),
           loc_name = str_c(contig, position, sep = '-'),
           .before = genotype) %$%
    
    new('genlight', gen = t(genotype), ploidy = 2,
        loc.all = snp_loc, loc.names = loc_name, 
        chromosome = contig, position = position)
}

remove_monomorphic <- function(snp_clone){
  snp_mat <- as.matrix(snp_clone)
  #n_ind <- dim(snp_mat)[1]
  
  all_0 <- colSums(snp_mat == 0, na.rm = TRUE) / colSums(!is.na(snp_mat))
  # all_2 <- colSums(snp_clone == 2, na.rm = TRUE) / n_ind
  loci_to_remove <- unname(which(all_0 == 1 | all_0 == 0))
  message(scales::comma(length(loci_to_remove)), ' monomorphic loci removed')
  snp_clone[,-loci_to_remove]
}

remove_missing_pop <- function(snp_clone, pop){
  setPop(snp_clone) <- pop
  
  missing_loci_in_pops <- seppop(snp_clone) %>%
    map(~as.matrix(.x) %>%
          is.na(.) %>%
          colSums(.) %>%
          equals(nInd(.x)) %>%
          which()) %>%
    unlist %>%
    names %>%
    unique %>%
    str_remove(str_c(levels(pop(snp_clone)), collapse = '|')) %>%
    str_remove('^\\.') %>%
    unique
  
  message('Removed ', scales::comma(length(missing_loci_in_pops)), 
          ' loci which were all NAs in one of the populations')
  
  snp_clone[,!snp_clone@loc.names %in% missing_loci_in_pops]
}

genlight_to_gds <- function(gen_light, out_name){
  SNPRelate::snpgdsCreateGeno(gds.fn = out_name, 
                              genmat = as.matrix(gen_light), 
                              sample.id = gen_light@ind.names, 
                              snp.id = gen_light@loc.names, 
                              snp.chromosome = gen_light@chromosome, 
                              snp.position = gen_light@position, 
                              snp.allele = gen_light@loc.all, 
                              snpfirstdim = FALSE)
  
  genofile <- SNPRelate::snpgdsOpen(out_name)
  SNPRelate::snpgdsSummary(genofile)
  SNPRelate::snpgdsClose(genofile)
  out_name
}

get_gds <- function(gds, column){
  z <- 0
  if(is.character(gds)){
    gds <- snpgdsOpen(gds)
    z <- 1
  }
  out <- gdsfmt::index.gdsn(gds, index = column) %>%
    gdsfmt::read.gdsn()
  
  if(z == 1){
    snpgdsClose(gds)
  }
  out
}

choose_snps <- function(genomic_data, samples_use = NULL, ...){
  genofile <- snpgdsOpen(genomic_data)
  if(is.null(samples_use)){
    samples_use <- get_gds(genofile, 'sample.id')
  }
  
  snps_used <- snpgdsSelectSNP(genofile, autosome.only = FALSE,
                               sample.id = samples_use, ...)
  snpgdsClose(genofile)
  snps_used
}

#### Read in data ####
sample_data <- read_csv('../intermediate_files/collected_sample_metadata.csv', 
                        show_col_types = FALSE)

genomic_data <- process_angsd(genotype_file = '../Data/genotypes.csv',
                              probability_file = '../Data/genotype_probabilities.csv',
                              sample_names = sample_data$ID, 
                              cutoff = confidence_cutoff)
strata(genomic_data) <- sample_data

#### Identify Clones ####
clone_data <- as.snpclone(genomic_data)

png('../Results/clone_removal.png', width = 7, height = 7, units = 'in', res = 125)
filter_stat_data <- filter_stats(clone_data, distance = 'bitwise.dist', plot = TRUE)
cutoff_clones <- cutoff_predictor(filter_stat_data$farthest$THRESHOLDS) 
abline(v = cutoff_clones, col = 'red') 
dev.off()

mlg.filter(clone_data, distance = 'bitwise.dist', algorithm = "farthest") <- cutoff_clones
# mlg.table(clone_data)

updated_metadata <- mutate(sample_data, 
                           clone_group = str_c('clone', 
                                               mll(clone_data, type = 'contracted'), 
                                               sep = '_'))

#### Initial Missing Data per Individual ####
# choose clone to keep least missing data initially
metadata_clones_removed <- as.matrix(genomic_data) %>%
  is.na %>%
  rowSums() %>%
  divide_by(nLoc(genomic_data)) %>%
  enframe(name = 'ID', value = 'pct_missing') %>%
  left_join(updated_metadata, ., by = 'ID') %>%
  group_by(clone_group) %>%
  filter(pct_missing == min(pct_missing)) %>%
  ungroup %>%
  filter(pct_missing < individual_missingness)

#### Write GDS files for each of two options ####
#Apalm v Acerv with Baum data
apalm_acerv <- genomic_data[genomic_data@ind.names %in% metadata_clones_removed$ID[metadata_clones_removed$species != 'Apr']] %>%
  remove_monomorphic() %>%
  remove_missing_pop(~species)

apalm_acerv_gds <- genlight_to_gds(apalm_acerv, '../intermediate_files/initial_apalm_acerv.gds')

#Acerv just our data
acerv <- genomic_data[genomic_data@ind.names %in% metadata_clones_removed$ID[metadata_clones_removed$data_origin == 'vollmer']] %>%
  remove_monomorphic() %>%
  remove_missing_pop(~location)

acerv_gds <- genlight_to_gds(acerv, '../intermediate_files/initial_acerv.gds')

#### Filter by Locus missingness & MAF ####
species_locus_keep <- choose_snps(apalm_acerv_gds, 
                                  maf = minor_allele_frequency, 
                                  missing.rate = locus_missingness)
filtered_apalm_acerv <- apalm_acerv[,species_locus_keep]


acerv_locus_keep <- choose_snps(acerv_gds, 
                                maf = minor_allele_frequency, 
                                missing.rate = locus_missingness)
filtered_acerv <- acerv[,acerv_locus_keep]

#### Write metadata, RDS of genlight and gds files for future analysis ####
metadata_out <- metadata_clones_removed %>%
  select(-pct_missing) %>%
  left_join(as.matrix(filtered_apalm_acerv) %>%
              is.na %>%
              rowSums() %>%
              divide_by(nLoc(filtered_apalm_acerv)) %>%
              enframe(name = 'ID', value = 'pct_missing_species'),
            by = 'ID') %>%
  left_join(as.matrix(filtered_acerv) %>%
              is.na %>%
              rowSums() %>%
              divide_by(nLoc(filtered_acerv)) %>%
              enframe(name = 'ID', value = 'pct_missing_location'),
            by = 'ID') 
write_csv(metadata_out, '../intermediate_files/preprocessed_metadata.csv')

write_rds(filtered_apalm_acerv, '../intermediate_files/preprocessed_apalm_acerv.rds', compress = 'xz')
write_rds(filtered_acerv, '../intermediate_files/preprocessed_acerv.rds', compress = 'xz')

genlight_to_gds(filtered_apalm_acerv, '../intermediate_files/preprocessed_apalm_acerv.gds')
genlight_to_gds(filtered_apalm_acerv, '../intermediate_files/preprocessed_acerv.gds')