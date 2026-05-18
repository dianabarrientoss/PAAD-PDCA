##############
# Script: Mutational Signatures Analysis for PAAD
# Used by: Yuliana Sosa
# Date: 14/03/26
# Description: Analysis of Mutational Signatures for PAAD dataset, with parameters 
# nmf_replicates=100 and maximum_signaures=25
#############

### STEP 1: LOAD ENVIRONMENT ###
# Indicates libraries position
.libPaths('/mnt/data/bioinfo-estadistica-2/ysosa/env/libs_R')

# Load libraries
library(reticulate)
library(tidyverse)
library(SigProfilerMatrixGeneratorR)
library(SigProfilerPlottingR)
library(SigProfilerExtractorR)
library(SigProfilerAssignmentR)

# Base directory
setwd('/mnt/data/bioinfo-estadistica-2/ysosa/funcional/PAAD/')

# Load Conda Environment
use_condaenv('/mnt/data/bioinfo-estadistica-2/ysosa/env/mutational_signatures', required = TRUE)

### STEP 2: LOADING CANCER SEQUENCING DATA ###
# Read cBioPortal style MAF file
maf_cbioportal = read.delim('/mnt/data/bioinfo-estadistica-2/ysosa/funcional/PAAD/data_mutations.txt')

# Selection of specific columns needed by SigProfiler
maf_sp = maf_cbioportal %>%
  select(Hugo_Symbol, Entrez_Gene_Id, Center, NCBI_Build, Chromosome,
         Start_Position, End_Position, Strand, Variant_Classification,
         Variant_Type, Reference_Allele, Tumor_Seq_Allele1,
         Tumor_Seq_Allele2, dbSNP_RS, dbSNP_Val_Status, Tumor_Sample_Barcode)

# Filter for only considering single base substitutions -> for SBS96 analysis
maf_sp = maf_sp %>%
  filter(Variant_Type == 'SNP')

# Create new folder for signature analysis results
# dir.create('signatures')
# Create new folder for updated MAF file (needed for SigProfilerMatrixGenerator)
# dir.create('signatures/SPMG/')

# Write updated MAF file
write.table(maf_sp, 
            '/mnt/data/bioinfo-estadistica-2/ysosa/funcional/PAAD/signatures/SPMG/data_mutations.maf',
            quote = FALSE, 
            row.names = FALSE, 
            sep = '\t')

cat("Ready! Rows:", nrow(maf_sp), "\n")
# We got 30,836 SNPs

### STEP 3: GENERATING MUTATIONAL MATRICES ###

# Install reference genome (only required once, previously done)
# message("Installing GRCh37 reference (skipped if already present)...")
# install('GRCh37', rsync=FALSE, bash=TRUE)

# Generates matrix of 96 possible mutations, according to reference genome (bases before and after)
message("Generating SBS96 matrix...")
matrices <- SigProfilerMatrixGeneratorR(project = "PAAD",
                                        genome = "GRCh37",
                                        matrix_path = "./signatures/SPMG",
                                        plot = F,
                                        exome = T)

### STEP 4: VISUALIZING MUTATIONAL PROFILES ###
# MUTATIONAL PROFILES FOR EACH SAMPLE
plotSBS(matrix_path = 'signatures/SPMG/output/SBS/PAAD.SBS96.exome',
        output_path = 'signatures/SPMG/output/SBS/',
        project = 'PAAD',
        plot_type = '96',
        percentage = FALSE)

# AVERAGE MUTATIONAL PROFILES
# Generate average mutational profiles
mut_matrix = matrices[['96']]

# Get relative mutational matrix, to avoid bias from hypermutatios
relative_mut_matrix = apply(mut_matrix, 2, prop.table)

# Get average mutational matrix
# Calculates mean of each row
average_mut_matrix = rowMeans(relative_mut_matrix)
average_mut_matrix = data.frame(Average_PAAD = average_mut_matrix)

# Add row names as column and print
average_mut_matrix_to_print = cbind(rownames(average_mut_matrix),
                                    average_mut_matrix)
colnames(average_mut_matrix_to_print)[1] = 'MutationType'
write.table(average_mut_matrix_to_print, 'signatures/avg_PAAD.SBS96.all',
            quote = F, row.names = F, sep = '\t')

# Plot average mutational profiles with percentages
plotSBS(matrix_path = 'signatures/avg_PAAD.SBS96.all',
        output_path = 'signatures/',
        project = 'avg_PAAD',
        plot_type = '96',
        percentage = TRUE)

# AVERAGE MUTATIONAL PROFILES PER SUBGROUP
# Read clinical file with metadata
metadata = read.delim('paad_tcga_pan_can_atlas_2018_clinical_data.tsv')

# Filtering metadata file to use only samples where we have mutation information
metadata = metadata %>%
  filter(Sample.ID %in% maf_sp$Tumor_Sample_Barcode)

# Get subtypes in our data
message("Subtypes found in metadata:")
print(table(metadata$Subtype))

# We'll use the PAAD subgroup
TARGET_SUBTYPE <- "PAAD"

# Get samples from group
samples_group = metadata %>%
  filter(Subtype == 'PAAD') %>%
  pull(Sample.ID)

# If there aren't samples, we skip this subtype
if (length(samples_group) == 0) {
  warning(paste("No samples found for subtype:", TARGET_SUBTYPE,
                "– skipping subtype plot. Check the table printed above."))
} else {
  # If we have PAAD samples...
  #select group samples from main matrix and get average
  mm_group <- rowMeans(relative_mut_matrix[, samples_group, drop = FALSE])
  mm_group <- data.frame(mm_group)
  
  #add row names as column and print
  mm_group_to_print <- cbind(MutationType = rownames(mm_group), mm_group)
  colnames(mm_group_to_print) <- c('MutationType', TARGET_SUBTYPE)
  
  # Writes table only with PAAD info
  out_file <- paste0('signatures/avg_', TARGET_SUBTYPE, '.SBS96.all')
  write.table(mm_group_to_print, out_file,
              quote = FALSE, row.names = FALSE, sep = '\t')
  
  # Create average mutational profiles from that particular type
  plotSBS(matrix_path = out_file,
          output_path = 'signatures/',
          project     = paste0('avg_', TARGET_SUBTYPE),
          plot_type   = '96',
          percentage  = TRUE)
}

### STEP 5: EXTRACTING DE NOVO MUTATIONAL SIGNATURES ###
sigprofilerextractor(input_type = 'matrix',
                     output = 'signatures/SPE/',
                     input_data = 'signatures/SPMG/output/SBS/PAAD.SBS96.exome',
                     nmf_replicates = 100,
                     minimum_signatures = 1,
                     maximum_signatures = 25,
                     exome = T)

### STEP 6: COSMIC REFERENCE MUTATIONAL SIGNATURES
# Done with portal

### STEP 7: ASIGNING REFERENCE MUTATIONAL SIGNATURES ###
# Run assignment analysis
cosmic_fit(samples = 'signatures/SPMG/output/SBS/PAAD.SBS96.exome',
           output = 'signatures/SPA',
           input_type='matrix',
           exome = T)

### STEP 8: DOWNSTREAM ANALYSIS OF SIGNATURE ASSIGNMENT RESULTS ###
stats = read.delim('signatures/SPA/Assignment_Solution/Solution_Stats/Assignment_Solution_Samples_Stats.txt')

## COSINE SIMILARITY
ggplot(stats) +
  aes(x=Cosine.Similarity) +
  labs(x='')+
  geom_histogram(aes(y = after_stat(density))) +
  geom_density(col = 4, lwd = 1.5) +
  geom_vline(aes(xintercept = 0.9),
             col = 2, lwd = 1.5) +
  labs(x = 'Cosine Similarity') +
  theme_bw()
ggsave('signatures/PAAD_cosine_similarity.pdf',
       width = 6, height = 4)

## SIGNATURE ACTIVITIES
# Read activities matrix
acts = read.delim('signatures/SPE/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt')

# Calculate average activities per signature
avg_acts = colMeans(acts[,-1])
message("Average activities per signature:")
print(avg_acts)
# Visualize average activities per signature

# Transform avg_acts into compatible dataframe with ggplot
df_avg_acts <- data.frame(
  Signature = names(avg_acts),
  Average_Mutations = as.numeric(avg_acts)
)

# Generate graph with ggplot
ggplot(df_avg_acts) +
  aes(x = Signature, y = Average_Mutations, fill = Signature) +
  geom_bar(stat = 'identity', show.legend = FALSE) +  # show.legend = FALSE porque los nombres ya están en el eje X
  theme_bw() +
  labs(
    title = "PAAD - Average signature activities",
    x = "Signatures",
    y = "Average number of mutations"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
    plot.title = element_text(face = "bold", hjust = 0.5)           
  )

# Save graph
ggsave('PAAD_avg_activities_ggplot.pdf', width = 6, height = 5)

# Reformat dataframe to use ggplot
acts_tidy = acts %>%
  pivot_longer(cols = !Samples,
               names_to = 'Signature',
               values_to = 'Mutations')

## STACKED BARPLOT
# Generate stacked barplot (percent stacked)
ggplot(acts_tidy) +
  aes(x = Samples, y = Mutations, fill = Signature) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  labs(title = "PAAD - Activities per sample") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave('signatures/PAAD_activities_all_samples.pdf',
       width = 10, height = 5)

## TOP 10 MUTATED CASES
# Calculate number of mutations per sample
number_of_mutations = rowSums(acts[,-1])

# Selecting the activities of only the top 10 mutated cases
top_10_mutated_samples = acts[order(number_of_mutations,
                                    decreasing = T)[1:10],]

# Reformatting and plotting
top_10_mutated_samples %>%
  pivot_longer(cols = !Samples,
               names_to = 'Signature',
               values_to = 'Mutations') %>%
  ggplot() +
  aes(x = reorder(Samples, Mutations), y = Mutations, fill = Signature) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  labs(x = 'Samples', title="PAAD - Top 10 most mutated samples")  +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave('signatures/PAAD_TOP10.pdf',
       width = 8, height = 5)


### Associating signature with specific metadata.
## SIGNATURE ACTIVITIES PER SUBTYPE
# Merge activities and metadata tables
acts_and_metadata = acts_tidy %>%
  rename(Sample.ID = Samples) %>%
  left_join(metadata,  by = 'Sample.ID')

# Calculate average activities per subtype
acts_per_subgroup = acts_and_metadata %>%
  group_by(Subtype, Signature) %>%
  summarise(Avg_mutations = mean(Mutations), .groups = 'drop') %>%
  filter(grepl('PAAD', Subtype))

head(acts_per_subgroup)

# Selecting only PAAD subtypes
acts_per_subgroup = acts_per_subgroup %>%
  filter(grepl('PAAD', Subtype))

# Plotting stacked barplot per subtype
ggplot(acts_per_subgroup) +
  aes(x = reorder(Subtype, Avg_mutations), y = Avg_mutations, fill = Signature) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  labs(x = 'PAAD Subtype',
       y = 'Average number of mutations',
       title = 'PAAD - Signature activities per subtype')

ggsave('signatures/PAAD_persubgroup.pdf',
       width = 8, height = 5)
