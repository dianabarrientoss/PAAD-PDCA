# Diana Barrientos
# 26/03/2026

library(dndscv)
maf = read.table("data_mutations.txt", header=T, sep="\t", 
                 stringsAsFactors=F, comment.char="#")

# Verificar columnas disponibles, como ya se realizó pues ya no se volverá a hacer apra la nueva corrida
#head(maf)
#colnames(maf)
#de maf a mut para dndscv
muts = data.frame(
  sampleID = maf$Tumor_Sample_Barcode,
  chr      = maf$Chromosome,
  pos      = maf$Start_Position,
  ref      = maf$Reference_Allele,
  mut      = maf$Tumor_Seq_Allele2
)
#exploracion del archivo mut
head(muts)

length(unique(muts$sampleID))

nrow(muts)
#mutaicones por samples
number_muts_IDs <- table(muts$sampleID) #We count the number of mutations of each donor
#hipermutadores
hypermuts <- muts[muts$sampleID %in% names(number_muts_IDs[number_muts_IDs >= 500]), ] #reduce the original dataset to those that are hypermutators
print(number_muts_IDs[unique(hypermuts$sampleID)])#We print the names and number of mutations of each hypermutator
#barplot con carga mutacional
pdf("barplot_PAAD.pdf", width=20, height=10)
barplot(sort(table(muts$sampleID)),
        ylab="Number of mutations",
        xlab="Donors",
        names.arg=NA,   
        col="pink")

dev.off()
#correr dndnscv con filtros
dout = dndscv(muts,
max_muts_per_gene_per_sample=3,max_coding_muts_per_sample=500,outmats=T)
#nombres d eoutput
names(dout)
#driver genes significativos

dout$sel_cv[which(dout$sel_cv$qglobal_cv<0.1),]
#dnds global
print(dout$globaldnds)
# Load the the Cancer Gene Census (v81) genes
data("cancergenes_cgc81", package="dndscv")
dout_cancergenes = dndscv(muts, outmats=T, gene_list=known_cancergenes)
#Hotspots con sitednds
sout = sitednds(dout_cancergenes)
names(sout)
sout_recursites <- sout$recursites[which(sout$recursites$qval<0.1),]
freq_sout_gene <- matrix(sout_recursites$freq, nrow = 1)
colnames(freq_sout_gene) <- sout_recursites$gene
cat("\nHotspot genes and their frequencies:\n")
print(freq_sout_gene)
#codones recurrentes
data("refcds_hg19", package = "dndscv")
RefCDS_codon = buildcodon(RefCDS)
codon_dnds = codondnds(dout_cancergenes, RefCDS_codon,
theta_option="conservative", min_recurr=2)
codon_dnds$recurcodons[which(codon_dnds$recurcodons$qval<0.1),]
