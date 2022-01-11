n <- ncol(dati.geni.log) / 2
tissue <- rep(1:2, each = n)
design <- cbind(Grp1 = 1, Grp2 = tissue - 1)

##### moderated t-test #####
fit <- lmFit(dati.geni.log, design)
fit <- eBayes(fit)
topBonf <- topTable(fit, coef = 2, p.value = .05, number = 30000, lfc = 1, adjust.method = "bonferroni")
topBH <- topTable(fit, coef = 2, p.value = .05, number = 30000, lfc = 1, adjust.method = "BH")
dim(topBonf)
summary(topBonf)
dim(topBH)
summary(topBH)
write.csv2(topBonf, file = "mod_t_test_bonf.csv")
write.csv2(topBH, file = "mod_t_test_bh.csv")

# Volcano plot
volcanoplot(fit, coef = 2, highlight = 0)

# Mean-difference plot
plotMD(fit, column = 2)

# Q-Q plot of moderated t-statistics
qqt(fit$t[,2], df = fit$df.residual + fit$df.prior)
abline(0,1)


##### Significance analysis of microarrays #####
sam.obj <- sam(dati.geni.log, c(-(1:n), 1:n))
str(summary(sam.obj, delta = 2.4))
sam.obj.pv.Bonf <- p.adjust(sam.obj@p.value, method = "bonferroni")
summary(sam.obj.pv.Bonf)
sam.obj.pv.BH <- p.adjust(sam.obj@p.value, method = "BH")
summary(sam.obj.pv.BH)

plot(sam.obj, 2.4)

write.csv2(sam.obj.pv.Bonf, file = "sam_bonf.csv")
write.csv2(sam.obj.pv.BH, file = "sam_bh.csv")


##### Empirical Bayes Analysis of Microarrays #####
find.obj <- find.a0(dati.geni.log, c(-(1:n), 1:n), rand = 123)
ebam.obj <- ebam(find.obj, delta = c(.90, .95, .99))
summary(ebam.obj, delta = .95)

write.csv2(summary(ebam.obj, delta = .95)@mat.sig, file = "ebam_bh.csv")

##### hy-test #####
dati.geni.log.tolist <- as.list(data.frame(t(dati.geni.log)))

out <- pbmclapply(dati.geni.log.tolist, 
                  function(gene, tissue) hy.test(gene ~ tissue, data = data.frame(gene, tissue)), 
                  tissue = tissue, mc.cores = 4L)

pvals <- sapply(out, function(x) x$p.value)

hist(pvals)
table(pvals <= .05)

hist(p.adjust(pvals, method = "bonferroni"))
table(p.adjust(pvals, method = "bonferroni") <= .05)

hist(p.adjust(pvals, method = "BH"))
table(p.adjust(pvals, method = "BH") <= .05)

idPvals05Bonf <- which(p.adjust(pvals, method = "bonferroni") <= .05)
idPvals05BH <- which(p.adjust(pvals, method = "BH") <= .05)

hist(pvals[idPvals05Bonf])
hist(pvals[idPvals05BH])
outPvals05Bonf <- out[idPvals05Bonf] 
outPvals05BH <- out[idPvals05BH] 

plot(outPvals05Bonf[[4]]) 

vpj <- list("hy-test" = names(outPvals05Bonf),
            "mod t-test" = rownames(topBonf),
            "sam" = names(sam.obj.pv.Bonf[sam.obj.pv.Bonf <= .05]))

colori<- gray.colors(n = 3, start = .3, end = .7)
colori_tr <- adjustcolor(colori, alpha.f = 0.2)
colori_bordo<- adjustcolor(colori, alpha.f = 0.5)

venn.diagram(vpj, filename = paste0(TISSUE, "_DEA.png"),
             imagetype = "png" ,
             height = 640 , 
             width = 640 , 
             resolution = 300,
             compression = "lzw",
             main.cex = 0.5,
             main = "",#paste0(ifelse(TISSUE == "brca", "A)","B)"), " DEA"),
             lwd = 1,
             col = colori ,
             fill = colori_tr,
             cex = 0.5,  # grandezza of area lab
             fontfamily = "sans",
             cat.cex = 0.5,    ## cat si riferisce ai nomi dei cerchi
             cat.default.pos = "outer",
             cat.fontfamily = "sans",
             cat.col = colori ,
             cat.dist = c(0.1, .1, .05),
             cat.pos = 0,
             margin = 0.0,
             scaled = FALSE ,
             force.unique = FALSE)

vpj <- list("hy-test" = names(outPvals05BH),
            "moderated t-test" = rownames(topBH),
            "sam" = names(sam.obj.pv.Bonf[sam.obj.pv.BH <= .05]))#,
            # "ebam" = names(summary(ebam.obj, delta = .95)@row.sig.genes))

colori<- c("mediumseagreen", "brown1", "blue")#, "yellow")
colori_tr <- adjustcolor(colori, alpha.f = 0.2)
colori_bordo<- adjustcolor(colori, alpha.f = 0.5)

venn.diagram(vpj, filename = "VennBH.png" ,
             imagetype = "png" ,
             height = 640, 
             width = 640, 
             resolution = 300,
             compression = "lzw",
             main.cex=0.5,
             main = "",
             lwd = 1,
             col = colori ,
             fill = colori_tr,
             cex = 0.5,  # grandezza of area lab
             fontfamily = "sans",
             cat.cex = 0.5,    ## cat si riferisce ai nomi dei cerchi
             cat.default.pos = "outer",
             cat.fontfamily = "sans",
             cat.col = colori ,
             cat.dist = c(0.05, .05, .05),
             cat.pos = 1,
             margin = 0.0,
             scaled = FALSE ,
             force.unique = FALSE)


save.image(paste0(TISSUE, "_deg_analysis_", format(Sys.Date(), "%d%m%y"), ".RData"))
