.libPaths(c("~/R/library", .libPaths()))
library(maftools)
library(maftools)

PAAD = read.maf(maf = "data_mutations.txt")

#oncoplots de dif estilos
pdf("oncoplot_sencillo_PAAD.pdf", width=30, height=10)
oncoplot(maf = PAAD, genes = c("TP53","KRAS","CDKN2A","SMAD4", "RNF43", "ARID1A"))
dev.off()

pdf("oncoplot_con_titv_PAAD.pdf", width=30, height=10)
oncoplot(maf = PAAD, genes = c("TP53","KRAS","CDKN2A","SMAD4", "RNF43", "ARID1A"), draw_titv= TRUE)
dev.off()

pdf("oncoplot_vias_PAAD.pdf", width=30, height=10)
oncoplot(maf = PAAD, pathways= 'sigpw', gene_mar =5, fontSize =0.3)
dev.off()

#dashboard de resumen
pdf("maf_summary_PAAD.pdf", width=15, height=10)
plotmafSummary(maf=PAAD, addStat = 'median', dashboard = TRUE)
dev.off()

#lolliplot plots de hotspots
# KRAS - hotspot en codón 12
pdf("lollipop_KRAS_PAAD.pdf", width=10, height=5)
lollipopPlot(maf=PAAD, gene="KRAS", refSeqID="NM_004985", showMutationRate=TRUE, labelPos=12)
dev.off()

# p16INK4a
pdf("lollipop_CDKN2A_p16_PAAD.pdf", width=10, height=5)
lollipopPlot(maf=PAAD, gene="CDKN2A", refSeqID="NM_000077", showMutationRate=TRUE)
dev.off()

# p14ARF
pdf("lollipop_CDKN2A_p14_PAAD.pdf", width=10, height=5)
lollipopPlot(maf=PAAD, gene="CDKN2A", refSeqID="NM_058195", showMutationRate=TRUE)
dev.off()

#oncoplot normal
pdf("oncoplot_normal_PAAD.pdf")
oncoplot(maf=PAAD)
dev.off()
#
pdf("plotTivTv.pdf")
PAAD_titv = titv(maf = PAAD, plot=FALSE)
plotTiTv(res=PAAD_titv)
dev.off()
