##### PUBMED research after GO enrichment #####

# tissue <- ifelse(TISSUE == "kirc", "\"kidney\"", "\"breast\"")
# tumor <- "\"cancer\""
# termini <- paste0("\"", risult.all2$Term.hy, "\"")
universe_text <- ""
tissue <- ifelse(TISSUE == "kirc", "kidney", "breast")
tumor <- paste(tissue, "cancer")
terms <- paste0(risult.all2$Term.hy)
tumor_and_terms <- paste(tumor, terms, sep = " AND ")
GO_ID <- risult.all2$GO.ID

# numer of the urn
universe <- EUtilsSummary(universe_text, type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2020)
# num_N <- QueryCount(universe)
# number of white balls
white_balls <- EUtilsSummary(tumor, type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2020)
# num_m <- QueryCount(white_balls)
# number of black balls
# num_n <- num_N - num_m
# number of balls drawn
# balls_drawn <- try(EUtilsSummary(terms[i], type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2020), silent = TRUE)
# while(inherits(balls_drawn, "try-error")) balls_drawn <- try(EUtilsSummary(terms[i], type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2020), silent = TRUE)
# num_k <- QueryCount(balls_drawn)
# # number of white balls drawn
# number_white_balls_drawn <- try(EUtilsSummary(tumor_and_terms[i], type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2020), silent = TRUE)
# while(inherits(number_white_balls_drawn, "try-error")) number_white_balls_drawn <- try(EUtilsSummary(tumor_and_terms[i], type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2020), silent = TRUE)
# num_x <- QueryCount(number_white_balls_drawn)
# # compute p-value
# phyper(num_x - 1, m = num_m, n = num_n, k = num_k, lower.tail = FALSE)

# hypergeometric test for pubmed research of significant terms
pvalHyper <- function(universe, tumor, terms, tumor_and_terms){
  # numer of the urn
  if(!inherits(universe, "EUtilsSummary")) 
    universe <- EUtilsSummary(universe, type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2020)
  num_N <- QueryCount(universe)
  text_N <- universe
  res_N <- universe@querytranslation
  # number of white balls
  if(!inherits(tumor, "EUtilsSummary")) 
    white_balls <- EUtilsSummary(tumor, type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2020)
  else white_balls <- tumor
  num_m <- QueryCount(white_balls)
  text_m <- tumor
  res_m <- white_balls@querytranslation
  # number of black balls
  num_n <- num_N - num_m
  # number of balls drawn
  if(!inherits(terms, "EUtilsSummary")) {
    balls_drawn <- try(EUtilsSummary(terms, type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2020), silent = TRUE)
    while(inherits(balls_drawn, "try-error")) balls_drawn <- try(EUtilsSummary(terms, type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2020), silent = TRUE)
  } else balls_drawn <- terms
  num_k <- QueryCount(balls_drawn)
  text_k <- terms
  res_k <- balls_drawn@querytranslation
  # number of white balls drawn
  if(!inherits(tumor_and_terms, "EUtilsSummary")) {
    number_white_balls_drawn <- try(EUtilsSummary(tumor_and_terms, type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2020), silent = TRUE)
    while(inherits(number_white_balls_drawn, "try-error")) number_white_balls_drawn <- try(EUtilsSummary(tumor_and_terms, type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2020), silent = TRUE)
  } else number_white_balls_drawn <- tumor_and_terms
  num_x <- QueryCount(number_white_balls_drawn)
  text_x <- tumor_and_terms
  res_x <- number_white_balls_drawn@querytranslation
  # calcolo del p-val
  pval <- phyper(num_x - 1, num_m, num_n, num_k, lower.tail = FALSE)
  
  out <- data.frame(num_N = num_N, num_m = num_m, num_n = num_n, num_k = num_k, num_x = num_x, pval = pval)
  attr(out, "details") <- list(
    text_N = text_N, res_N = res_N, 
    text_m = text_m, res_m = res_m, 
    text_k = text_k, res_k = res_k, 
    text_x = text_x, res_x = res_x)
  out
}

out2 <- NULL
pb <- txtProgressBar(min = 0, max = length(terms), style = 3)
for(i in 1:length(terms)){  #
  setTxtProgressBar(pb, i)
  
  temp <- pvalHyper(universe = universe, tumor = white_balls, terms = terms[i], tumor_and_terms = tumor_and_terms[i])
  out2 <- rbind(out2, temp)
  print(paste0("Iter ", i, " pval ", out2[i, 6L]))
}
close(pb)

out2 <- data.frame(termini = risult.all2$Term.hy, GO_ID = GO_ID, out2)
out2$pval.bonf <- p.adjust(out2$pval, method = "bonferroni")
out2$pval.BH <- p.adjust(out2$pval, method = "BH")

sum(out2$pval.bonf <= .05)
sum(out2$pval.BH <= .05)

dim(subset(risult.all, rif.BH.our == 1))
sum(subset(risult.all, rif.BH.our == 1)$Term.hy %in% out2[out2$pval.BH <= .05, "termini"])
dim(subset(risult.all, rif.BH.t == 1))
sum(subset(risult.all, rif.BH.t == 1)$Term.t %in% out2[out2$pval.BH <= .05, "termini"])
dim(subset(risult.all, rif.BH.sam == 1))
sum(subset(risult.all, rif.BH.sam == 1)$Term %in% out2[out2$pval.BH <= .05, "termini"])
# dim(subset(risult.all, rif.BH.ebam == 1))
# sum(subset(risult.all, rif.BH.ebam == 1)$Term.ebam %in% out2[out2$pval.BH <= .05, "termini"])

dim(subset(risult.all, rif.bonf.our == 1))
sum(subset(risult.all, rif.bonf.our == 1)$Term.hy %in% out2[out2$pval.bonf <= .05, "termini"])
dim(subset(risult.all, rif.bonf.t == 1))
sum(subset(risult.all, rif.bonf.t == 1)$Term.t %in% out2[out2$pval.bonf <= .05, "termini"])
dim(subset(risult.all, rif.bonf.sam == 1))
sum(subset(risult.all, rif.bonf.sam == 1)$Term %in% out2[out2$pval.bonf <= .05, "termini"])
# dim(subset(risult.all, rif.bonf.ebam == 1))
# sum(subset(risult.all, rif.bonf.ebam == 1)$Term.ebam %in% out2[out2$pval.bonf <= .05, "termini"])

# vpj <- list("hy-test" = subset(risult.all, rif.BH.our == 1)$Term.hy[subset(risult.all, rif.BH.our == 1)$Term.hy %in% out2[out2$pval.BH <= .05, "termini"]], 
#             "moderated t-test" = subset(risult.all, rif.BH.t == 1)$Term.t[subset(risult.all, rif.BH.t == 1)$Term.t %in% out2[out2$pval.BH <= .05, "termini"]],
#             "sam" = subset(risult.all, rif.BH.sam == 1)$Term[subset(risult.all, rif.BH.sam == 1)$Term %in% out2[out2$pval.BH <= .05, "termini"]])#,
# # "ebam" = subset(risult.all, rif.BH.ebam == 1)$Term.ebam[subset(risult.all, rif.BH.ebam == 1)$Term.ebam %in% out2[out2$pval.BH <= .05, "termini"]])
# colori <- c("mediumseagreen", "brown1", "blue")#, "yellow")
# colori_tr <- adjustcolor(colori, alpha.f = 0.2)
# colori_bordo <- adjustcolor(colori, alpha.f = 0.5)
# 
# venn.diagram(vpj, filename = "VennTermsBH.png",
#              imagetype = "png" ,
#              height = 640 , 
#              width = 640 , 
#              resolution = 300,
#              compression = "lzw",
#              main.cex = 0.5,
#              main = "Venn of terms",
#              lwd = 1,
#              col = colori ,
#              fill = colori_tr,
#              cex = 0.5,  # grandezza of area lab
#              fontfamily = "sans",
#              cat.cex = 0.5,    ## cat si riferisce ai nomi dei cerchi
#              cat.default.pos = "outer",
#              cat.fontfamily = "sans",
#              cat.col = colori ,
#              cat.dist = c(0.05, 0.05, .05),
#              cat.pos = 0,
#              margin = 0.05 ,
#              scaled = FALSE ,
#              force.unique = FALSE)
# 
# vpj <- list("hy-test" = subset(risult.all, rif.bonf.our == 1)$Term.hy[subset(risult.all, rif.bonf.our == 1)$Term.hy %in% out2[out2$pval.bonf <= .05, "termini"]], 
#             "moderated t-test" = subset(risult.all, rif.bonf.t == 1)$Term.t[subset(risult.all, rif.bonf.t == 1)$Term.t %in% out2[out2$pval.bonf <= .05, "termini"]],
#             "sam" = subset(risult.all, rif.bonf.sam == 1)$Term[subset(risult.all, rif.bonf.sam == 1)$Term %in% out2[out2$pval.bonf <= .05, "termini"]])#,
# # "ebam" = subset(risult.all, rif.bonf.ebam == 1)$Term.ebam[subset(risult.all, rif.bonf.ebam == 1)$Term.ebam %in% out2[out2$pval.bonf <= .05, "termini"]])
# 
# venn.diagram(vpj, filename = "VennTermsBonf.png",
#              imagetype = "png" ,
#              height = 640 , 
#              width = 640 , 
#              resolution = 300,
#              compression = "lzw",
#              main.cex = 0.5,
#              main = "Venn of terms",
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
#              margin = 0.05 ,
#              scaled = FALSE ,
#              force.unique = FALSE)


risult.all3 <- subset(risult.all, rif.bonf.our == 1 | rif.bonf.t == 1 | rif.bonf.sam == 1)# | rif.BH.ebam == 1)
out3 <- subset(out2, termini %in% risult.all3$Term.hy)
out3$pval.bonf <- p.adjust(out3$pval, method = "bonferroni")
out3$pval.BH <- p.adjust(out3$pval, method = "BH")

sum(out3$pval.bonf <= .05)
sum(out3$pval.BH <= .05)

dim(subset(risult.all, rif.bonf.our == 1))
sum(subset(risult.all, rif.bonf.our == 1)$Term.hy %in% out3[out3$pval.bonf <= .05, "termini"])
dim(subset(risult.all, rif.bonf.t == 1))
sum(subset(risult.all, rif.bonf.t == 1)$Term.t %in% out3[out3$pval.bonf <= .05, "termini"])
dim(subset(risult.all, rif.bonf.sam == 1))
sum(subset(risult.all, rif.bonf.sam == 1)$Term %in% out3[out3$pval.bonf <= .05, "termini"])


vpj <- list("hy-test" = subset(risult.all, rif.bonf.our == 1)$Term.hy[subset(risult.all, rif.bonf.our == 1)$Term.hy %in% out3[out3$pval.bonf <= .05, "termini"]], 
            "mod t-test" = subset(risult.all, rif.bonf.t == 1)$Term.t[subset(risult.all, rif.bonf.t == 1)$Term.t %in% out3[out3$pval.bonf <= .05, "termini"]],
            "sam" = subset(risult.all, rif.bonf.sam == 1)$Term[subset(risult.all, rif.bonf.sam == 1)$Term %in% out3[out3$pval.bonf <= .05, "termini"]])#,
# "ebam" = subset(risult.all, rif.bonf.ebam == 1)$Term.ebam[subset(risult.all, rif.bonf.ebam == 1)$Term.ebam %in% out2[out2$pval.bonf <= .05, "termini"]])

venn.diagram(vpj, filename = paste0(TISSUE, "_PubMedResearch.png"),
             imagetype = "png" ,
             height = 640 , 
             width = 640 , 
             resolution = 300,
             compression = "lzw",
             main.cex = 0.5,
             main = "",#paste0(ifelse(TISSUE == "brca", "A)","B)"), " PubMed research"),
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

risult.all4 <- subset(risult.all, Term %in% out3[out3$pval.bonf <= .05, "termini"])[, c("Term", "Annotated", "Significant.hy", "Significant.t", "Significant")]
risult.all4 <- merge.data.frame(risult.all4, out3, by.x = "Term", by.y = "termini") 
risult.all4 <- merge.data.frame(risult.all4, risult.all3[, c("Term", "rif.bonf.our", "rif.bonf.t", "rif.bonf.sam")], by.x = "Term", by.y = "Term")
risult.all4$PVAL <- ifelse(risult.all4$pval.bonf < 1.11E-16, "<1.11E-16", formatC(risult.all4$pval.bonf, digits = 2, width = 3, format = "E"))
risult.all4$group <- apply(risult.all4[,15:17], 1, function(x) ifelse(all(x == 1), "all", 
                                                                      ifelse(x[1] == 1 & x[2] == 1, "hy-test/mod t-test", 
                                                                             ifelse(x[1] == 1 & x[3] == 1, "hy-test/sam", 
                                                                                    ifelse(x[2] == 1 & x[3] == 1, "mod t-test/sam", 
                                                                                           ifelse(x[1] == 1, "hy-test", 
                                                                                                  ifelse(x[2] == 1, "mod t-test", "sam")))))))
risult.all4$group <- factor(risult.all4$group, 
                            levels = c("hy-test", "mod t-test", "sam", "hy-test/mod t-test", 
                                       "hy-test/sam", "mod t-test/sam", "all"))
table(risult.all4$group)
write.csv2(risult.all4[order(risult.all4$group, risult.all4$Term),], "risultati_all4_new.csv" , row.names = F)

save.image(paste0(TISSUE, "_pubmed_research_", format(Sys.Date(), "%d%m%y"), "_new.RData"))
