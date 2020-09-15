library(fgsea)

##create rank and gmt file
rnk.file <- system.file("extdata", "dmso_upper.rnk", package="fgsea")
gmt.file <- system.file("extdata", "ALL.gmt", package="fgsea")

##loading ranks
ranks <- read.table(rnk.file,
                    header=TRUE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$t, ranks$ID)
str(ranks)

##loading pathways
pathways <- gmtPathways(gmt.file)
str(head(pathways))

##run fgsea
fgseaRes <- fgsea(pathways = pathways,
                  stats = ranks,
                  minSize=5,
                  maxSize=500,
                  nperm=10000)

sum(fgseaRes[, padj < 0.01])
sum(fgseaRes[, padj < 0.05])

##plotting
library(ggplot2)
library(ggpubr)
a <- plotEnrichment(pathways[["tissue_morphogenesis(4)"]],ranks) + labs(title="Tissue Morphogenesis", caption="padj=0.003806")

b <- plotEnrichment(pathways[["tube_morphogenesis(4)"]],ranks) + labs(title="Tube Morphogenesis", caption="padj=0.003806")

c <- plotEnrichment(pathways[["tissue_remodeling(4)"]],ranks) + labs(title="Tissue Remodeling", caption="padj=0.012072")

d <- plotEnrichment(pathways[["developmental_maturation(4)"]],ranks) + labs(title="Developmental Maturation", caption="padj=0.003806")

e <- plotEnrichment(pathways[["cell_fate_commitment(5)"]],ranks) + labs(title="Cell Fate Commitment", caption="padj=0.003806")

f <- plotEnrichment(pathways[["mesenchymal_to_epithelial_transition(7)"]],ranks) + labs(title="MET", caption="padj=0.007929")

g <- plotEnrichment(pathways[["spinal_cord_development(4)"]],ranks) + labs(title="Spinal Cord Development", caption="padj=0.003806")

h <- plotEnrichment(pathways[["skeletal_system_development(5)"]],ranks) + labs(title="Skeletal System Development", caption="padj=0.004896")

i <- plotEnrichment(pathways[["angiogenesis(4)"]],ranks) + labs(title="Angiogenesis", caption="padj=0.003806")



plot <- ggarrange(a,b,c,d,e,f,g,h,i, ncol=3, nrow=3) + theme_bw(base_size=16)


##output table
data <- as.matrix(fgseaRes)
write.table(data, sep="\t", file="ALL_enrichment_dmso_beccari.xls")
