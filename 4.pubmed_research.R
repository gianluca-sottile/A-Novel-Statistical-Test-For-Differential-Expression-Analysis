##### PUBMED research after GO enrichment #####

tissue <- ifelse(TISSUE == "kirc", "\"kidney\"", "\"breast\"")
tumor <- "\"cancer\""
termini <- paste0("\"", risult.all2$Term.hy, "\"")
GO_ID <- risult.all2$GO.ID

tissue_and_termini <- paste(tissue, termini, sep = " AND ")
tumor_and_termini <- paste(tumor, termini, sep = " AND ")
tissue_and_tumor_and_terms <- paste(tissue, tumor, termini, sep = " AND ")

# hypergeometric test for pubmed research of significant terms
pvalHyper2 <- function(termini, tissue_and_termini, tumor_and_termini, tissue_and_tumor_and_terms){
  res_N <- EUtilsSummary(termini, type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2020)
  num_N <- QueryCount(res_N)
  
  res_k <- res_m <- res_x <- NA
  if(num_N != 0){ 
    # per trovare k - numero di palle estratte
    # text_k <- tissue_and_termini[i]
    res_k <- EUtilsSummary(tissue_and_termini, type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2020)
    num_k <- QueryCount(res_k)
    
    #per trovare x - intersezione
    # text_x <- tissue_and_tumor_and_terms[i]
    res_m <- EUtilsSummary(tumor_and_termini, type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2020)
    num_m <- QueryCount(res_m)
    
    ### per trovare i casi sfavorevoli fare N-m
    num_b = num_N - num_m
    
    if(num_k != 0 & num_m != 0){ 
      res_x <- EUtilsSummary(tissue_and_tumor_and_terms, type="esearch", db="pubmed", datetype='pdat', mindate=2000, maxdate=2020)
      num_x <- QueryCount(res_x)
      
      # calcolo del p-val
      pval <- phyper(num_x - 1, num_m, num_b, num_k, lower.tail = FALSE)
    }
    else{
      num_x <- 0
      pval <- 1
    }
  }
  else{
    num_k <- 0
    num_m <- 0
    num_x <- 0
    num_b <- num_N
    pval <- 1
  }
  
  out <- data.frame(num_N = num_N, num_k = num_k, num_m = num_m, num_b = num_b, num_x = num_x, pval = pval)
  attr(out, "details") <- list(
    text_N = termini, res_N = res_N@querytranslation, 
    text_k = tissue_and_termini, res_k = ifelse(isS4(res_k), res_k@querytranslation, res_k), 
    text_m = tumor_and_termini, res_m = ifelse(isS4(res_m), res_m@querytranslation, res_m), 
    text_x = tissue_and_tumor_and_terms, res_x = ifelse(isS4(res_x), res_x@querytranslation, res_x))
  out
}

out2 <- NULL
pb <- txtProgressBar(min = 0, max = length(termini), style = 3)
for(i in 1:length(termini)){  #
  setTxtProgressBar(pb, i)
  
  temp <- try(pvalHyper2(termini = termini[i], 
                         tissue_and_termini = tissue_and_termini[i],
                         tumor_and_termini = tumor_and_termini[i], 
                         tissue_and_tumor_and_terms = tissue_and_tumor_and_terms[i]), 
              silent = TRUE)
  while(inherits(temp, "try-error")) temp <- try(pvalHyper2(termini = termini[i],
                                                            tissue_and_termini = tissue_and_termini[i],
                                                            tumor_and_termini = tumor_and_termini[i],
                                                            tissue_and_tumor_and_terms = tissue_and_tumor_and_terms[i]),
                                                 silent = TRUE)
  out2 <- rbind(out2, temp)
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
sum(subset(risult.all, rif.BH.sam == 1)$Term.sam %in% out2[out2$pval.BH <= .05, "termini"])
dim(subset(risult.all, rif.BH.ebam == 1))
sum(subset(risult.all, rif.BH.ebam == 1)$Term.ebam %in% out2[out2$pval.BH <= .05, "termini"])

dim(subset(risult.all, rif.bonf.our == 1))
sum(subset(risult.all, rif.bonf.our == 1)$Term.hy %in% out2[out2$pval.bonf <= .05, "termini"])
dim(subset(risult.all, rif.bonf.t == 1))
sum(subset(risult.all, rif.bonf.t == 1)$Term.t %in% out2[out2$pval.bonf <= .05, "termini"])
dim(subset(risult.all, rif.bonf.sam == 1))
sum(subset(risult.all, rif.bonf.sam == 1)$Term.sam %in% out2[out2$pval.bonf <= .05, "termini"])
dim(subset(risult.all, rif.bonf.ebam == 1))
sum(subset(risult.all, rif.bonf.ebam == 1)$Term.ebam %in% out2[out2$pval.bonf <= .05, "termini"])


vpj <- list("hy-test" = subset(risult.all, rif.BH.our == 1)$Term.hy[subset(risult.all, rif.BH.our == 1)$Term.hy %in% out2[out2$pval.BH <= .05, "termini"]], 
            "moderated t-test" = subset(risult.all, rif.BH.t == 1)$Term.t[subset(risult.all, rif.BH.t == 1)$Term.t %in% out2[out2$pval.BH <= .05, "termini"]],
            "sam" = subset(risult.all, rif.BH.sam == 1)$Term.sam[subset(risult.all, rif.BH.sam == 1)$Term.sam %in% out2[out2$pval.BH <= .05, "termini"]],
            "ebam" = subset(risult.all, rif.BH.ebam == 1)$Term.ebam[subset(risult.all, rif.BH.ebam == 1)$Term.ebam %in% out2[out2$pval.BH <= .05, "termini"]])
colori <- c("mediumseagreen", "brown1", "blue", "yellow")
colori_tr <- adjustcolor(colori, alpha.f = 0.2)
colori_bordo <- adjustcolor(colori, alpha.f = 0.5)

venn.diagram(vpj, filename = "VennTermsBH.png",
             imagetype = "png" ,
             height = 640 , 
             width = 640 , 
             resolution = 300,
             compression = "lzw",
             main.cex = 0.5,
             main = "Venn of terms",
             lwd = 1,
             col = colori ,
             fill = colori_tr,
             cex = 0.5,  # grandezza of area lab
             fontfamily = "sans",
             cat.cex = 0.5,    ## cat si riferisce ai nomi dei cerchi
             cat.default.pos = "outer",
             cat.fontfamily = "sans",
             cat.col = colori ,
             cat.dist = c(0.2, 0.2, .1, .1),
             cat.pos = 0,
             margin = 0.0 ,
             scaled = FALSE ,
             force.unique = FALSE)

vpj <- list("hy-test" = subset(risult.all, rif.bonf.our == 1)$Term.hy[subset(risult.all, rif.bonf.our == 1)$Term.hy %in% out2[out2$pval.bonf <= .05, "termini"]], 
            "moderated t-test" = subset(risult.all, rif.bonf.t == 1)$Term.t[subset(risult.all, rif.bonf.t == 1)$Term.t %in% out2[out2$pval.bonf <= .05, "termini"]],
            "sam" = subset(risult.all, rif.bonf.sam == 1)$Term.sam[subset(risult.all, rif.bonf.sam == 1)$Term.sam %in% out2[out2$pval.bonf <= .05, "termini"]],
            "ebam" = subset(risult.all, rif.bonf.ebam == 1)$Term.ebam[subset(risult.all, rif.bonf.ebam == 1)$Term.ebam %in% out2[out2$pval.bonf <= .05, "termini"]])

venn.diagram(vpj, filename = "VennTermsBonf.png",
             imagetype = "png" ,
             height = 640 , 
             width = 640 , 
             resolution = 300,
             compression = "lzw",
             main.cex = 0.5,
             main = "Venn of terms",
             lwd = 1,
             col = colori ,
             fill = colori_tr,
             cex = 0.5,  # grandezza of area lab
             fontfamily = "sans",
             cat.cex = 0.5,    ## cat si riferisce ai nomi dei cerchi
             cat.default.pos = "outer",
             cat.fontfamily = "sans",
             cat.col = colori ,
             cat.dist = c(0.2, 0.2, .1, .1),
             cat.pos = 0,
             margin = 0.0 ,
             scaled = FALSE ,
             force.unique = FALSE)


save.image(paste0(TISSUE, "_pubmed_research_", format(Sys.Date(), "%d%m%y"), ".RData"))
