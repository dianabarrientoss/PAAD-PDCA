.libPaths(c("~/R/library", .libPaths()))
library(maftools)
library(maftools)

PAAD = read.maf(maf = "data_mutations.txt")

#nota los codones se sacaron del primer Rscript

# KRAS - hotspots G12 y Q61
pdf("lollipop_KRAS_codones_PAAD.pdf", width=10, height=5)
lollipopPlot(maf=PAAD, gene="KRAS", refSeqID="NM_004985",
             showMutationRate=TRUE, labelPos=c(12, 61))
dev.off()

# CDKN2A p16INK4a - hotspot H83
pdf("lollipop_CDKN2A_p16_codones_PAAD.pdf", width=10, height=5)
lollipopPlot(maf=PAAD, gene="CDKN2A", refSeqID="NM_000077",
             showMutationRate=TRUE, labelPos=83)
dev.off()

# CDKN2A p14ARF - hotspots A97 y P94
pdf("lollipop_CDKN2A_codones_p14_PAAD.pdf", width=10, height=5)
lollipopPlot(maf=PAAD, gene="CDKN2A", refSeqID="NM_058195",
             showMutationRate=TRUE, labelPos=c(94, 97))
dev.off()
