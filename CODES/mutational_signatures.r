.libPaths('/mnt/data/bioinfo-estadistica-2/dbarrientos/env/Rlibs')
#cargar librerías
library(tidyverse)
library(reticulate)
library(SigProfilerMatrixGeneratorR)
library(SigProfilerPlottingR)
library(SigProfilerExtractorR)
library(SigProfilerAssignmentR)

#ambiente de conda
use_condaenv('/mnt/data/bioinfo-estadistica-2/dbarrientos/env/mutational_signatures', required = TRUE)

#CARGAR DATA
message("Loading MAF file...")
maf_cbioportal <- read.delim('data_mutations.txt')

# Selection of specific columns needed by SigProfiler
maf_sp = maf_cbioportal %>%
  select(Hugo_Symbol, Entrez_Gene_Id, Center, NCBI_Build, Chromosome,
         Start_Position, End_Position, Strand, Variant_Classification,
         Variant_Type, Reference_Allele, Tumor_Seq_Allele1,
         Tumor_Seq_Allele2, dbSNP_RS, dbSNP_Val_Status, Tumor_Sample_Barcode)

# only considering single base substitutions
maf_sp = maf_sp %>%
    filter(Variant_Type == 'SNP')
message(paste("SNPs retained:", nrow(maf_sp)))

#ESTRUCTURA DE LOS DATOS
# generate folders
dir.create('signatures')
dir.create('signatures/SPMG/')

# Write updated MAF file
write.table(maf_sp, 'signatures/SPMG/data_mutations.maf', quote = F,
            row.names = F, sep = '\t')

#CREAR LA MATRIZ MUTACIONAL
message("Installing GRCh37 reference (skipped if already present)...")
install('GRCh37', rsync=FALSE, bash=TRUE)

message("Generating SBS96 matrix...")
matrices <- SigProfilerMatrixGeneratorR(project = "PAAD",
                                        genome = "GRCh37",
                                        matrix_path = "./signatures/SPMG",
                                        plot = F,
                                        exome = T)
#HACER PLOTS POR SAMPLE MUTATIONAL
plotSBS(matrix_path = 'signatures/SPMG/output/SBS/PAAD.SBS96.exome',
        output_path = 'signatures/SPMG/output/SBS/',
        project = 'PAAD',
        plot_type = '96',
        percentage = FALSE)

#PROMEDIO EN PERFIL MUTACIONAL DE TODAS LAS MUESTRAS
# Generate average mutational profiles
mut_matrix = matrices[['96']]

# Get relative mutational matrix
relative_mut_matrix = apply(mut_matrix, 2, prop.table)

# Get average mutational matrix
average_mut_matrix = rowMeans(relative_mut_matrix)
average_mut_matrix = data.frame(Average_PAAD = average_mut_matrix)

# Add row names as column and print
average_mut_matrix_to_print = cbind(rownames(average_mut_matrix),
                                    average_mut_matrix)

colnames(average_mut_matrix_to_print)[1] = 'MutationType'

write.table(average_mut_matrix_to_print, 'signatures/avg_PAAD.SBS96.all',
            quote = F, row.names = F, sep = '\t')

# Plot average mutational profiles (note the percentage parameter now)
plotSBS(matrix_path = 'signatures/avg_PAAD.SBS96.all',
        output_path = 'signatures/',
        project = 'avg_PAAD',
        plot_type = '96',
        percentage = TRUE)


#CLINICAL METADATA
message("Loading clinical metadata...")
# Read clinical file with metadata
metadata = read.delim('paad_tcga_pan_can_atlas_2018_clinical_data.tsv')

# Filtering metadata file to use only samples where we have mutation information
metadata = metadata %>%
    filter(Sample.ID %in% maf_sp$Tumor_Sample_Barcode)

message(paste("Samples with both mutation and clinical data:", nrow(metadata)))

# Subtype-level average profile
message("Subtypes found in metadata:")
print(table(metadata$Subtype))

# PAAD BASAL BELOW
# Get samples from group
# MAL (puede tener caracteres invisibles)
TARGET_SUBTYPE <- "PAAD_Basal"

# BIEN — usa el subtype que realmente existe en tus datos
TARGET_SUBTYPE <- "PAAD"

samples_group = metadata %>%
    filter(Subtype == TARGET_SUBTYPE) %>%
    pull(Sample.ID)

if (length(samples_group) == 0) {
  warning(paste("No samples found for subtype:", TARGET_SUBTYPE,
                "– skipping subtype plot. Check the table printed above."))
} else {
#select group samples from main matrix and get average
  mm_group <- rowMeans(relative_mut_matrix[, samples_group, drop = FALSE])
  mm_group <- data.frame(mm_group)
 #add row names as column and print
  mm_group_to_print <- cbind(MutationType = rownames(mm_group), mm_group)
  colnames(mm_group_to_print) <- c('MutationType', TARGET_SUBTYPE)
 
  out_file <- paste0('signatures/avg_', TARGET_SUBTYPE, '.SBS96.all')
  write.table(mm_group_to_print, out_file,
              quote = FALSE, row.names = FALSE, sep = '\t')
 
  plotSBS(matrix_path = out_file,
          output_path = 'signatures/',
          project     = paste0('avg_', TARGET_SUBTYPE),
          plot_type   = '96',
          percentage  = TRUE)
}


# DE NOVO SIGNATURE EXTRACTION

sigprofilerextractor(
  input_type         = 'matrix',
  output             = 'signatures/SPE/',
  input_data         = 'signatures/SPMG/output/SBS/PAAD.SBS96.exome',
  nmf_replicates     = 3,
  minimum_signatures = 1,
  maximum_signatures = 3,
  exome              = TRUE
)
#COSMIC SIGNATUR ASSIGNMENT 
# Run assignment analysis
cosmic_fit(samples = 'signatures/SPMG/output/SBS/PAAD.SBS96.exome',
           output = 'signatures/SPA',
           input_type='matrix',
           exome = T)
#Cosine-similarity QC plot
stats = read.delim('signatures/SPA/Assignment_Solution/Solution_Stats/Assignment_Solution_Samples_Stats.txt')

ggplot(stats) +
  aes(x = Cosine.Similarity) +
  geom_histogram(aes(y = after_stat(density)), bins = 30) +
  geom_density(col = 4, lwd = 1.5) +
  geom_vline(xintercept = 0.9, col = 2, lwd = 1.5) +
  labs(x = 'Cosine Similarity', title = 'PAAD - Assignment quality') +
  theme_bw()

ggsave('signatures/PAAD_cosine_similarity.pdf',
       width = 6, height = 4)

#SIGNATURES CTIVITIES
# Read activities matrix
acts <- read.delim(
  'signatures/SPE/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt')

SBS96_Activities <- acts
# Average activities per signature
avg_acts <- colMeans(acts[, -1])
message("Average activities per signature:")
print(avg_acts)


# Visualize average activities per signature
# Save barplot of average activities
# Barplot average activities
pdf('signatures/PAAD_avg_activities_barplot.pdf')
barplot(avg_acts, main = 'PAAD - Average signature activities', col="pink")
dev.off()

# Tidy format
acts_tidy <- acts %>%
  pivot_longer(cols      = !Samples,
               names_to  = 'Signature',
               values_to = 'Mutations')

# Generate stacked barplot (percent stacked)
ggplot(acts_tidy) +
  aes(x = Samples, y = Mutations, fill = Signature) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  labs(title = 'PAAD - Signature activities (all samples)') +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

ggsave('signatures/PAAD_activities_all_samples.pdf',
       width = 10, height = 5)


#TOP 10 MOST MUTATED SAMPLES
# Calculate number of mutations per sample
number_of_mutations = rowSums(acts[,-1])
# Selecting the activities of only the top 10 mutated cases
top_10_mutated_samples = acts[order(number_of_mutations,
                                    decreasing = T)[1:10],]




# Reformatting and plotting
top_10_mutated_samples %>%
  pivot_longer(cols = !Samples, names_to = 'Signature', values_to = 'Mutations') %>%
  ggplot() +
  aes(x = reorder(Samples, Mutations), y = Mutations, fill = Signature) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  labs(title = 'PAAD - Top 10 most mutated samples') +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

ggsave('signatures/PAAD_activities_top10.pdf', width=8, height=5)

#activities per ubtype


# Merge activities and metadata tables
# (The samples column needs to be renamed in one of them)
acts_and_metadata = acts_tidy %>%
    rename(Sample.ID = Samples) %>%
    left_join(metadata,  by = 'Sample.ID')

# Calculate average activities per subtype
acts_per_subgroup = acts_and_metadata %>%
    group_by(Subtype, Signature) %>%
    summarise(Avg_mutations = mean(Mutations), .groups = 'drop') %>%
  filter(grepl('PAAD', Subtype))

head(acts_per_subgroup)

# Selecting only COAD subtypes
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

ggsave('signatures/PAAD_activities_per_subtype.pdf',
       width = 8, height = 5)
# SIGMATURE ACTIVITY FOR SPECIFIC METADATA

# BRCA1/2 Homologous Recombination Deficiency (SBS3)
message("Checking BRCA1/2 variants and SBS3 activity...")

# Muestras con variantes en BRCA1 o BRCA2
maf_BRCA <- maf_cbioportal %>%
  filter(Hugo_Symbol %in% c('BRCA1', 'BRCA2'))

Samples_with_BRCA_var <- maf_BRCA %>%
  distinct(Tumor_Sample_Barcode)

message(paste("Samples with BRCA1/2 variants:", nrow(Samples_with_BRCA_var)))

# Actividad de SBS3 en esas muestras
if ('SBS3' %in% colnames(SBS96_Activities)) {
  SBS3_activities <- SBS96_Activities %>%
    select(Samples, SBS3)

  message("SBS3 activity in BRCA1/2 variant samples:")
  print(
    SBS3_activities %>%
      filter(Samples %in% Samples_with_BRCA_var$Tumor_Sample_Barcode) %>%
      arrange(desc(SBS3))
  )

  # Comparar SBS3 entre muestras con y sin variante BRCA
  SBS3_activities <- SBS3_activities %>%
    mutate(BRCA_variant = ifelse(Samples %in% Samples_with_BRCA_var$Tumor_Sample_Barcode,
                                  'BRCA1/2 variant', 'No variant'))

  p_brca <- ggplot(SBS3_activities) +
    aes(x = BRCA_variant, y = SBS3, fill = BRCA_variant) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_bw() +
    labs(title = 'PAAD – SBS3 activity by BRCA1/2 status',
         x = '', y = 'SBS3 mutations assigned') +
    theme(legend.position = 'none')

  ggsave('signatures/PAAD_BRCA_SBS3.pdf', p_brca, width = 6, height = 5)

} else {
  message("SBS3 not detected in this cohort.")
}

# Carga mutacional de muestras con variante BRCA
if (nrow(Samples_with_BRCA_var) > 0) {
  maf_BRCA_sign <- maf_cbioportal %>%
    filter(Tumor_Sample_Barcode %in% Samples_with_BRCA_var$Tumor_Sample_Barcode)

  p_brca_burden <- ggplot(data = maf_BRCA_sign) +
    geom_bar(mapping = aes(x = Tumor_Sample_Barcode), fill = "#AED6F1") +
    theme_bw() +
    labs(title = 'PAAD – Mutational burden in BRCA1/2 variant samples',
         x = 'Sample', y = 'Number of mutations') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  ggsave('signatures/PAAD_BRCA_burden.pdf', p_brca_burden, width = 8, height = 5)
}


