library(tidyverse)
library(magrittr)
library(xml2)
library(adegenet)
library(poppr)
library(reticulate)
# conda_install(channel = 'bioconda', packages = 'pysradb')

main_dir <- '23August2022'
loci_missing_max <- 0.3 #If individuals are missing more than this much data remove initially

#### Functions ####
srr_to_srs <- function(data){
  pysradb <- import('pysradb')
  db = pysradb$SRAweb()
  
  srr_ids <- filter(data, !is.na(SRR_id)) %>%
    pull(SRR_id)
  
  out <- db$srr_to_srs(srr_ids, detailed = FALSE, 
                       sample_attribute = FALSE, 
                       expand_sample_attributes = FALSE) %>%
    as_tibble() %>%
    select(run_accession, sample_accession) %>% 
    rename(SRS_id = sample_accession, SRR_id = run_accession)
  full_join(data, out, by = 'SRR_id')
}


set_pop <- function(gen_obj, pop_form){
  gen_obj <- setPop(gen_obj, pop_form)
  gen_obj
}

convert_to_decimal <- function(x){
  out <- str_split(x, 'W|N') %>%
    unlist(recursive = FALSE) %>%
    as.numeric() %>%
    magrittr::divide_by(c(1, 60)) %>%
    sum
  
  if_else(str_detect(x, 'W'), -1 * out, out)
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


#### Read in Data ####
#get metadata
baum_locations <- read_xml('metadata/biosample_result.xml') %>%
  as_list() %>%
  as_tibble %>%
  unnest_wider(BioSampleSet) %>%
  rowwise %>%
  mutate(SRS_id = Ids[[3]][[1]]) %>%
  mutate(location = Attributes[[3]][[1]] %>% str_remove(' \\(.*\\)')) %>%
  select(SRS_id, location)

florida_locations <- read_csv('metadata/SP.pop.data.csv', show_col_types = FALSE) %>%
  select(Genotype, subregion, reefName, x.coord, y.coord) %>%
  rename(gen_id = Genotype,
         lat = x.coord,
         lon = y.coord) %>%
  mutate(reef = str_c(subregion, reefName, sep = '-'), .keep = 'unused')

panama_locations <- tribble(
  ~'reef', ~'lat_garmin', ~'lon_garmin', 
  'tetas', '82W6.068', '9N16.579',
  'sebastians', '82W7.631', '9N15.274',
  'holy shit', '82W6.929', '9N16.773',
  'CK14', '82W7.557', '9N15.238',
  'CK4', '82W7.625', '9N15.517'
) %>%
  rowwise %>% 
  mutate(lat = convert_to_decimal(lat_garmin),
         lon = convert_to_decimal(lon_garmin)) %>%
  ungroup %>%
  select(-ends_with('garmin')) %>%
  mutate(reef = str_replace_all(reef, c('tetas' = 'Tet', 'sebastians' = 'SR', 'holy shit' = 'HS')))

all_individual_data <- read_delim(str_c(main_dir, 'bam.list', sep = '/'), delim = '\t', col_names = 'file', show_col_types = FALSE) %>%
  mutate(ID = str_remove(file, dirname(file)) %>% str_remove('/') %>% str_remove('-fp.*bam$')) %>%
  select(-file) %>%
  left_join(read_delim(str_c(main_dir, 'individual_missing_data.csv', sep = '/'), show_col_types = FALSE),
            by = 'ID') %>%
  mutate(SRR_id = str_extract(ID, 'SRR[0-9]+')) %>% 
  srr_to_srs() %>%
  left_join(baum_locations, by = 'SRS_id') %>%
  select(-SRR_id, -SRS_id) %>%
  mutate(gen_id = str_extract(ID, '[A-Za-z0-9]+$')) %>%
  left_join(florida_locations, by = 'gen_id') %>%
  mutate(data_origin = if_else(str_detect(ID, 'baum'), 'baum', 'vollmer'),
         species = str_extract(ID, 'A[cp]|ac|apr|apa') %>% str_to_sentence(),
         location = if_else(is.na(location), str_extract(ID, 'FL|PA'), location) %>% str_replace_all(c('^FL$' = 'Florida', '^PA$' = 'Panama')),
         reef = if_else(data_origin == 'vollmer' & location == 'Panama', str_extract(ID, 'CK4|CK14|HS|SR|Tet'), reef)) %>%
  select(ID, gen_id, data_origin, species, percent_missing_loci, location, reef, lat, lon) %>%
  left_join(panama_locations, by = 'reef') %>%
  mutate(lat = coalesce(lat.x, lat.y),
         lon = coalesce(lon.x, lon.y),
         .keep = 'unused')

#Read in Genotypes to genlight object
read_delim('../genome_assembly/contig_length.dat', delim = '\t', col_names = c('chrom', 'length'))


genotypes <- read_csv(str_c(main_dir, 'genotypes.csv', sep = '/'), 
                      col_types = cols(.default = col_integer(),
                                       contig = col_character(),
                                       ref = col_character(),
                                       alt = col_character())) %>%
  mutate(locus_number = row_number(),
         genotypes = cbind(!!!syms(all_individual_data$ID)),
         .keep = 'unused') %>%
  mutate(snp_loc = str_c(ref, alt, sep = '/'),
         loc_name = str_c(contig, position, sep = '-'),
         alt_freq = rowSums(genotypes, na.rm = TRUE) / (2 * (rowSums(!is.na(genotypes)))),
         ref_freq = 1 - alt_freq, 
         .before = genotypes) %>%
  filter(ref_freq != 1, alt_freq != 1) %$%
  
  new('genlight', gen = t(genotypes), ploidy = 2,
      loc.all = snp_loc, loc.names = loc_name, 
      chromosome = contig, position = position,
      strata = all_individual_data) %>%
  set_pop(~species/location/reef) %>%
  as.snpclone()

#### Identify Clones ####
png(str_c(main_dir, '/clone_removal.png'), width = 7, height = 7, units = 'in', res = 125)
filter_stat_data <- filter_stats(genotypes, distance = 'bitwise.dist', plot = TRUE)
abline(v = 0.03, col = 'red') 
dev.off()
# A genetic distance of 0.03 is right between the initial peak and the main body of the histogram
# use UPGMA or Farthest Neighbor
cutoff_predictor(filter_stat_data$farthest$THRESHOLDS) #matches what my eye says

mlg.filter(genotypes, distance = 'bitwise.dist', algorithm = "farthest") <- 0.03
mlg.table(genotypes)

updated_metadata <- mutate(all_individual_data, clone_group = str_c('clone', mll(genotypes, type = 'contracted'), sep = '_'))

#### Output files ####
#Subset individuals
updated_metadata %>%
  group_by(clone_group) %>%
  # filter(reef == 'HS') %>%
  mutate(keep = percent_missing_loci == min(percent_missing_loci)) %>%
  ungroup %>%
  mutate(keep = as.integer(keep & (percent_missing_loci < loci_missing_max))) %>%
  select(keep) %>%
  write_csv(str_c(main_dir, '/subset_noClones/subset_data.txt', sep = '/'), col_names = FALSE)


#Full MetaData
write_csv(updated_metadata, '../../Combo_Tank_Survival/Data/snp_metadata.csv')

#Genetics File
write_rds(genotypes, '../../Combo_Tank_Survival/Data/snp_cloneGenLight.rds', compress = 'xz')
