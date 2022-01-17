##### GO analysis #####
lista_geni <- names(dati.geni.log.tolist)
g.list <- unique(lista_geni)
tot.geni <- length(g.list)

# hy-test
geni.selezionati.our <- unique(names(idPvals05Bonf))
geneList.our <- factor(as.integer(g.list %in% geni.selezionati.our))
names(geneList.our) <- g.list

sampleGOdata.our <- new("topGOdata",
                        description = "Simple session", 
                        ontology = "BP",
                        allGenes = geneList.our,
                        annot = annFUN.org, 
                        mapping = "org.Hs.eg.db", 
                        ID = "ensembl")
sampleGOdata.our

resultFisher.our <- runTest(sampleGOdata.our, algorithm = "classic", statistic = "fisher")
total.score.our <- score(resultFisher.our)

allRes.our <- GenTable(sampleGOdata.our, classicFisher = resultFisher.our, 
                       orderBy = "classicFisher", ranksOf = "classicFisher", 
                       topNodes = length(total.score.our), numChar = 999)
allRes.our$classicFisher[allRes.our$classicFisher == "< 1e-30"] <- 1E-30
allRes.our$classicFisher <- as.numeric(allRes.our$classicFisher)

pval.our <- data.frame("terms" = names(total.score.our), "pval.our" = total.score.our)

allRes.our <- merge.data.frame(allRes.our, pval.our, by.x = "GO.ID", by.y = "terms", all = "TRUE")

# moderated t-test
geni.selezionati.t <- unique(rownames(topBonf))
geneList.t <- factor(as.integer(g.list %in% geni.selezionati.t))
names(geneList.t) <- g.list

sampleGOdata.t <- new("topGOdata",
                      description = "Simple session", 
                      ontology = "BP",
                      allGenes = geneList.t,
                      annot = annFUN.org, 
                      mapping = "org.Hs.eg.db", 
                      ID = "ensembl")
sampleGOdata.t

resultFisher.t <- runTest(sampleGOdata.t, algorithm = "classic", statistic = "fisher")
total.score.t <- score(resultFisher.t)

allRes.t <- GenTable(sampleGOdata.t, classicFisher = resultFisher.t,
                     orderBy = "classicFisher", ranksOf = "classicFisher", 
                     topNodes = length(total.score.t), numChar = 999)
allRes.t$classicFisher[allRes.t$classicFisher == "< 1e-30"] <- 1E-30
allRes.t$classicFisher <- as.numeric(allRes.t$classicFisher)

pval.t <- data.frame("terms" = names(total.score.t), "pval.t" = total.score.t)

allRes.t <- merge.data.frame(allRes.t, pval.t, by.x = "GO.ID", by.y = "terms", all = "TRUE")

# sam 
geni.selezionati.sam <- unique(names(sam.obj.pv.Bonf[sam.obj.pv.Bonf <= .05]))
geneList.sam <- factor(as.integer(g.list %in% geni.selezionati.sam))
names(geneList.sam) <- g.list

sampleGOdata.sam <- new("topGOdata",
                        description = "Simple session", 
                        ontology = "BP",
                        allGenes = geneList.sam,
                        annot = annFUN.org, 
                        mapping = "org.Hs.eg.db", 
                        ID = "ensembl")
sampleGOdata.sam

resultFisher.sam <- runTest(sampleGOdata.sam, algorithm = "classic", statistic = "fisher")
total.score.sam <- score(resultFisher.sam)

allRes.sam <- GenTable(sampleGOdata.sam, classicFisher = resultFisher.sam,
                       orderBy = "classicFisher", ranksOf = "classicFisher", 
                       topNodes = length(total.score.sam), numChar = 999)
allRes.sam$classicFisher[allRes.sam$classicFisher == "< 1e-30"] <- 1E-30
allRes.sam$classicFisher <- as.numeric(allRes.sam$classicFisher)

pval.sam <- data.frame("terms" = names(total.score.sam), "pval.sam" = total.score.sam)

allRes.sam <- merge.data.frame(allRes.sam, pval.sam, by.x = "GO.ID", by.y = "terms", all = "TRUE")

# ebam 
# geni.selezionati.ebam <- unique(rownames(summary(ebam.obj, delta = .95)@mat.sig))
# geneList.ebam <- factor(as.integer(g.list %in% geni.selezionati.ebam))
# names(geneList.ebam) <- g.list
# 
# sampleGOdata.ebam <- new("topGOdata",
#                         description = "Simple session", 
#                         ontology = "BP",
#                         allGenes = geneList.ebam,
#                         annot = annFUN.org, 
#                         mapping = "org.Hs.eg.db", 
#                         ID = "ensembl")
# sampleGOdata.ebam
# 
# resultFisher.ebam <- runTest(sampleGOdata.ebam, algorithm = "classic", statistic = "fisher")
# total.score.ebam <- score(resultFisher.ebam)
# 
# allRes.ebam <- GenTable(sampleGOdata.ebam, classicFisher = resultFisher.ebam,
#                        orderBy = "classicFisher", ranksOf = "classicFisher", 
#                        topNodes = length(total.score.ebam), numChar = 999)
# allRes.ebam$classicFisher[allRes.ebam$classicFisher == "< 1e-30"] <- 1E-30
# allRes.ebam$classicFisher <- as.numeric(allRes.ebam$classicFisher)
# 
# pval.ebam <- data.frame("terms" = names(total.score.ebam), "pval.ebam" = total.score.ebam)
# 
# allRes.ebam <- merge.data.frame(allRes.ebam, pval.ebam, by.x = "GO.ID", by.y = "terms", all = "TRUE")

# merge of all results
risultati <-  merge.data.frame(
  # merge.data.frame(
    merge.data.frame(allRes.our, allRes.t, by = "GO.ID", all = "TRUE", suffixes = c(".hy", ".t")), 
    allRes.sam, by = "GO.ID", all = TRUE, suffixes = c(".t", ".sam"))#, 
  # allRes.ebam, by = "GO.ID", all = TRUE, suffixes = c(".sam", ".ebam"))
risultati <- risultati[order(risultati$GO.ID), ]

# let consider the terms with at least one annotated gene 
risultati <- subset(risultati, Annotated.hy >= 1 & Annotated.t >= 1 & Annotated >= 1)# & Annotated.ebam >= 1)
str(risultati)
summary(risultati)

# complete previous data.frame with adjusted p-values and indicator of significant gene
bonf.our <- p.adjust(risultati$pval.our, method = "bonferroni")
BH.our <- p.adjust(risultati$pval.our, method = "BH")

bonf.t <- p.adjust(risultati$pval.t, method = "bonferroni")
BH.t <- p.adjust(risultati$pval.t, method = "BH")

bonf.sam <- p.adjust(risultati$pval.sam, method = "bonferroni")
BH.sam <- p.adjust(risultati$pval.sam, method = "BH")

# bonf.ebam <- p.adjust(risultati$pval.ebam, method = "bonferroni")
# BH.ebam <- p.adjust(risultati$pval.ebam, method = "BH")

risult.all <- cbind.data.frame(risultati, 
                               bonf.our, BH.our, bonf.t, BH.t, bonf.sam, BH.sam, #bonf.ebam, BH.ebam,
                               rif.bonf.our = 1*(bonf.our <= 0.05),
                               rif.BH.our = 1*(BH.our <= 0.05),
                               rif.bonf.t = 1*(bonf.t <= 0.05),
                               rif.BH.t = 1*(BH.t <= 0.05),
                               rif.bonf.sam = 1*(bonf.sam <= 0.05),
                               rif.BH.sam = 1*(BH.sam <= 0.05))#,
                               # rif.bonf.ebam = 1*(bonf.ebam <= 0.05),
                               # rif.BH.ebam = 1*(BH.ebam <= 0.05))

table(risult.all$rif.bonf.our)
table(risult.all$rif.bonf.t)
table(risult.all$rif.bonf.sam)
# table(risult.all$rif.bonf.ebam)
table(risult.all$rif.BH.our)
table(risult.all$rif.BH.t)
table(risult.all$rif.BH.sam)
# table(risult.all$rif.BH.ebam)

write.table(risult.all, "risultati_all_new.txt" , row.names = F)

vpj <- list("hy-test" = subset(risult.all, rif.bonf.our == 1)$Term.hy, 
            "mod t-test" = subset(risult.all, rif.bonf.t == 1)$Term.t,
            "sam" = subset(risult.all, rif.bonf.sam == 1)$Term)#,
            # "ebam" = subset(risult.all, rif.bonf.ebam == 1)$Term.ebam)
# colori <- c("mediumseagreen", "brown1", "blue")#, "yellow")
# colori_tr <- adjustcolor(colori, alpha.f = 0.2)
# colori_bordo <- adjustcolor(colori, alpha.f = 0.5)

venn.diagram(vpj, filename = paste0(TISSUE, "_EnrichmentAnalysis.png"),
             imagetype = "png" ,
             height = 640 , 
             width = 640 , 
             resolution = 300,
             compression = "lzw",
             main.cex = 0.5,
             main = "",#paste0(ifelse(TISSUE == "brca", "A)","B)"), " Enrichment Analysis"),
             lwd = 1,
             col = colori ,
             fill = colori_tr,
             cex = 0.5,  # grandezza of area lab
             fontfamily = "sans",
             cat.cex = 0.5,    ## cat si riferisce ai nomi dei cerchi
             cat.default.pos = "outer",
             cat.fontfamily = "sans",
             cat.col = colori ,
             cat.dist = c(0.05, 0.05, -.45),
             cat.pos = 0,
             margin = 0.0 ,
             scaled = FALSE ,
             force.unique = FALSE)

# vpj <- list("hy-test" = subset(risult.all, rif.BH.our == 1)$Term.hy, 
#             "moderated t-test" = subset(risult.all, rif.BH.t == 1)$Term.t,
#             "sam" = subset(risult.all, rif.BH.sam == 1)$Term)#,
# "ebam" = subset(risult.all, rif.bonf.ebam == 1)$Term.ebam)

# venn.diagram(vpj, filename = "VennGOBH.png",
#              imagetype = "png" ,
#              height = 640 , 
#              width = 640 , 
#              resolution = 300,
#              compression = "lzw",
#              main.cex = 0.5,
#              main = "Venn of GO terms",
#              lwd = 1,
#              col = colori ,
#              fill = colori_tr,
#              cex = 0.5,  # grandezza of area lab
#              fontfamily = "sans",
#              cat.cex = 0.5,    ## cat si riferisce ai nomi dei cerchi
#              cat.default.pos = "outer",
#              cat.fontfamily = "sans",
#              cat.col = colori ,
#              cat.dist = c(0.05, 0.05, -.45),
#              cat.pos = 0,
#              margin = 0.1 ,
#              scaled = FALSE ,
#              force.unique = FALSE)


# all(subset(risult.all, rif.bonf.our == 1 | rif.bonf.t == 1 | rif.bonf.sam == 1 | rif.bonf.ebam == 1)$Term.hy %in% 
#       subset(risult.all, rif.BH.our == 1 | rif.BH.t == 1 | rif.BH.sam == 1 | rif.BH.ebam == 1)$Term.hy)
risult.all2 <- subset(risult.all, rif.BH.our == 1 | rif.BH.t == 1 | rif.BH.sam == 1)# | rif.BH.ebam == 1)
dim(risult.all2)


save.image(paste0(TISSUE, "_go_enrichment_", format(Sys.Date(), "%d%m%y"), ".RData"))
