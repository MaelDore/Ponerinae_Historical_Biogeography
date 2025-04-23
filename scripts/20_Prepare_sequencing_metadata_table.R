##### Script 20: Prepare Sequencing metadata table #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Prepare Sequencing metadata table for Supplementary Data 2 (includes all taxa in the phylogeny)
# Provides summary stats from sequencing and biological attributes of the extractions
# Update Voucher table needed for raw sequence upload in NCBI SRA - BioSample with proper updated biomaterial information

###

### Inputs

# SD1: Voucher specimens metadata
# Extractions metadata from Global Taxon File
# Parsing for biological material


###

### Outputs

# SD2: Sequencing metadata table
# Voucher table following template provided by NCBI - BioSample (with proper updated biomaterial information)

###

# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### 1.1/ Load packages ####

library(tidyverse)
library(readxl)
library(xlsx)      # Need the Java Development Kit (JDK) installed
library(openxlsx)  # Use Rccp. No need of Java

### 1.2/ Load data ####

## Load SD1 table

# SD1_Voucher_specimens_metadata <- read.xlsx(xlsxFile = "./input_data/Molecular_data/SD1_Voucher_specimens_metadata.xlsx")
SD1_Voucher_specimens_metadata <- readRDS(file = "./input_data/Molecular_data/SD1_Voucher_specimens_metadata.rds")

## Load sequencing metadata
sequencing_metadata_df <- read.xlsx(xlsxFile = "./input_data/Molecular_data/ponerinae-792-uce-stats.xlsx")

## Load Global Taxon File
GTL_df <- read.xlsx(xlsxFile = "./input_data/Molecular_data/Global_Taxon_List_2025_01_13.xlsx", sheet = "taxa")

## Load parsing table for biological material
NCBI_biological_material_parsing <- read.xlsx(xlsxFile = "./input_data/Molecular_data/ncbi specimen parsing.xlsx")

## Load file for sample name conversion
sample_name_updates_df <- read.table(file = "./input_data/Molecular_data/taxonomy-updates.txt", header = T)
sample_name_updates_df$original_name <- str_split(string = sample_name_updates_df$original_name.updated_name, pattern = ",", simplify = T)[,1]
sample_name_updates_df$updated_name <- str_split(string = sample_name_updates_df$original_name.updated_name, pattern = ",", simplify = T)[,2]


##### 2/ Generate template and populate with already available metadata from SD1 ####

SD2_Sequencing_metadata <- SD1_Voucher_specimens_metadata %>% 
  select(Extraction_code, Taxa_name, Voucher_specimen_code, Outgroup, Status_data, BioProject_ID, BioSample_ID, Read_sequences_SRA_SRR_ID, COX1_GenBank_ID, Reference, Reference_DOI)

##### 3/ Merge with sequencing metadata ####

# Rename fields
sequencing_metadata_df_to_merge <- sequencing_metadata_df %>% 
  rename("Nb_UCEs" = `#uces`) %>% 
  rename("Total_length_bp" = `length.(bp)`) %>% 
  rename("Mean_length_bp" = `mean.length.(bp)`) %>% 
  rename("SD_length_bp" = `stderr.length.(bp)`) %>% 
  rename("Min_length_bp" = `min.length.(bp)`) %>% 
  rename("Median_length_bp" = `median.length.(bp)`) %>% 
  rename("Nb_contigs_longer_than_1Kbp" = `contigs.>.1kb`) %>%
  mutate(Nb_reads = NA, Mean_depth = NA, Perc_coverage = NA, Extraction_source = NA) %>%
  dplyr::select(samples, Extraction_source, Nb_UCEs, Nb_reads, Nb_contigs_longer_than_1Kbp, Mean_depth, Perc_coverage, Total_length_bp, Mean_length_bp, SD_length_bp, Min_length_bp, Median_length_bp)

# Update sample names
sequencing_metadata_df_to_merge$Sample_name <- sample_name_updates_df$updated_name[match(x = sequencing_metadata_df_to_merge$samples, table = sample_name_updates_df$original_name)]

# Extract Taxa names
Sample_names_list <- str_split(string = sequencing_metadata_df_to_merge$Sample_name, pattern = "_")
nb_items <- unlist(lapply(X = Sample_names_list, FUN = length))
Voucher_specimen_codes <- c()
Taxa_names_list <- list()
for (i in seq_along(Sample_names_list))
{
  nb_items_i <- nb_items[i]
  Voucher_specimen_codes[i] <- Sample_names_list[[i]][nb_items_i]
  Taxa_names_list[[i]] <- Sample_names_list[[i]][1:(nb_items_i-2)]
}

Taxa_names <- unlist(lapply(X = Taxa_names_list, FUN = paste0, collapse = "_"))
Taxa_names <- str_remove(string = Taxa_names, pattern = "_EX.*")
Taxa_names[Taxa_names == "NewGenus_bucki"] <- "Neoponera_bucki"

sequencing_metadata_df_to_merge$Taxa_name <- Taxa_names

table(sequencing_metadata_df_to_merge$Taxa_name %in% SD2_Sequencing_metadata$Taxa_name)

# Merge sequencing metadata
SD2_Sequencing_metadata <- SD2_Sequencing_metadata %>% 
  left_join(y = sequencing_metadata_df_to_merge) %>% 
  dplyr::select(-samples, -Sample_name)

##### 4/ Merge information on biological material ####

## Load Global Taxon File
GTL_df <- read.xlsx(xlsxFile = "./input_data/Molecular_data/Global_Taxon_List_2025_01_13.xlsx", sheet = "taxa")

# Deal with variations of extraction codes
EXcode_variants <- SD2_Sequencing_metadata$Extraction_code[!(SD2_Sequencing_metadata$Extraction_code %in% GTL_df$ExtractionCode)]
EXcode_initial <- str_remove(string = EXcode_variants, pattern = "a$")

GTL_df$ExtractionCode[GTL_df$ExtractionCode %in% EXcode_initial] <- EXcode_variants

table(GTL_df$ExtractionCode %in% SD2_Sequencing_metadata$Extraction_code)
table(SD2_Sequencing_metadata$Extraction_code %in% GTL_df$ExtractionCode)

GTL_df_to_merge <- GTL_df %>% 
  dplyr::select(ExtractionCode, material) %>%
  filter(ExtractionCode %in% SD2_Sequencing_metadata$Extraction_code)
         
# Deal with duplicates
EXcode_duplicates <- GTL_df_to_merge$ExtractionCode[duplicated(GTL_df_to_merge$ExtractionCode)]
View(GTL_df_to_merge[GTL_df_to_merge$ExtractionCode == EXcode_duplicates, ])

GTL_df_to_merge <- GTL_df_to_merge %>% 
  dplyr::distinct(ExtractionCode, .keep_all = T)

## Parse into proper fields

GTL_df_to_merge <- GTL_df_to_merge %>% 
  left_join(y = NCBI_biological_material_parsing, by = join_by(material == original)) %>% 
  dplyr::select(-material) %>% 
  rename(destructive_or_not = `D/ND`)

# Merge biological information
SD2_Sequencing_metadata <- SD2_Sequencing_metadata %>% 
  left_join(y = GTL_df_to_merge, by = join_by(Extraction_code == ExtractionCode))

##### 5/ Reorder and export #####

SD2_Sequencing_metadata <- SD2_Sequencing_metadata %>%
  dplyr::select(Extraction_code, Taxa_name, Voucher_specimen_code, Outgroup, caste, dev_stage, sex, tissue, destructive_or_not,
                Status_data, BioProject_ID, BioSample_ID, Read_sequences_SRA_SRR_ID, COX1_GenBank_ID, Reference, Reference_DOI, Extraction_source,
                Nb_UCEs, Nb_reads, Nb_contigs_longer_than_1Kbp, Mean_depth, Perc_coverage, Total_length_bp, Mean_length_bp, SD_length_bp, Min_length_bp, Median_length_bp)

# Save updated table
saveRDS(object = SD2_Sequencing_metadata, file = "./input_data/Molecular_data/SD2_Sequencing_metadata.rds")
write.xlsx(x = SD2_Sequencing_metadata, file = "./input_data/Molecular_data/SD2_Sequencing_metadata.xlsx")


##### 6/ Update the NCBI Biosample voucher table with proper biological information #####

NCBI_BioSample_Voucher_table <- readRDS(file = "./input_data/Molecular_data/NCBI_BioSample_Voucher_table.rds")
# NCBI_BioSample_Voucher_table <- read.xlsx(xlsxFile = "./input_data/Molecular_data/NCBI_BioSample_Voucher_table.xlsx")

Biomaterial_info_to_merge <- SD2_Sequencing_metadata %>% 
  dplyr::select(Voucher_specimen_code, caste, dev_stage, sex, tissue, destructive_or_not)
  
# Update tissue with N/BD info
Biomaterial_info_to_merge$tissue[!is.na(Biomaterial_info_to_merge$destructive_or_not)] <- paste0(Biomaterial_info_to_merge$tissue, " (",Biomaterial_info_to_merge$destructive_or_not," extraction)")[!is.na(Biomaterial_info_to_merge$destructive_or_not)]

# Fill NA with "missing"
Biomaterial_info_to_merge$tissue[is.na(Biomaterial_info_to_merge$tissue)] <- "missing"
Biomaterial_info_to_merge$caste[is.na(Biomaterial_info_to_merge$caste)] <- "missing"
Biomaterial_info_to_merge$dev_stage[is.na(Biomaterial_info_to_merge$dev_stage)] <- "missing"
Biomaterial_info_to_merge$sex[is.na(Biomaterial_info_to_merge$sex)] <- "missing"

Biomaterial_info_to_merge <- Biomaterial_info_to_merge %>% 
  rename(life_stage = caste) %>% 
  dplyr::select(-destructive_or_not)
  
table(NCBI_BioSample_Voucher_table$specimen_voucher %in% SD2_Sequencing_metadata$Voucher_specimen_code)

# Remove previous fields
NCBI_BioSample_Voucher_table <- NCBI_BioSample_Voucher_table %>% 
  dplyr::select(-tissue, -dev_stage, -sex)

# Merge biological information
NCBI_BioSample_Voucher_table <- NCBI_BioSample_Voucher_table %>% 
  left_join(y = Biomaterial_info_to_merge, by = join_by(specimen_voucher == Voucher_specimen_code)) %>% 
  dplyr::select(sample_name, bioproject_accession, organism, isolate, breed, host, isolation_source, collection_date, geo_loc_name, tissue, dev_stage, life_stage, sex, 
                altitude, biomaterial_provider, collected_by, identified_by, lat_lon, specimen_voucher, infra.specific.rank, infra.specific.name, description)

### Save/Export the Voucher table
saveRDS(object = NCBI_BioSample_Voucher_table, file = "./input_data/Molecular_data/NCBI_BioSample_Voucher_table.rds")
write.xlsx(x = NCBI_BioSample_Voucher_table, file = "./input_data/Molecular_data/NCBI_BioSample_Voucher_table.xlsx")
  
