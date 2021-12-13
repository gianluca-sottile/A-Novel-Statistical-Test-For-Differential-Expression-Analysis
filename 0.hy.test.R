hy.test <- function (x, ...) UseMethod("hy.test")

hy.test.default <- function (x, y, nthresh = 50L, paired = TRUE, ...) {
  method <- ifelse(paired, "Paired Hy-test", "Two Sample Hy-test")
  
  if(missing(x)) stop("'x' is missing for hyper test")
  if(missing(y)) stop("'y' is missing for hyper test")
  dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  xok <- yok <- complete.cases(x, y)
  x <- x[xok]
  y <- y[yok]
  
  nx <- length(x)
  ny <- length(y)
  if(nx != ny) stop("'x' and 'y' have different lengths")
  
  min.val <- max(min(x), min(y))
  minimum <- floor(min(min(x), min(y)))
  max.val <- min(max(x), max(y))
  maximum <- ceiling(max(max(x), max(y)))
  thresh <- seq(min.val, max.val, length.out = nthresh)
  thresh <- CJ("V1" = thresh, "V2" = thresh) 
  thresh <- thresh[which(thresh$V1 > thresh$V2)]
  nthresh <- nrow(thresh)
  cutoff <- data.table("geni.sani" = 1:nthresh, "geni.malati" = 1:nthresh)
  n.cutoff <- nrow(cutoff)
  
  tempp <- sapply(seq_len(nthresh), function(i) {
    x <- cut(x, breaks = c(minimum, thresh[i,], maximum), labels = c(-1, 0, 1), include.lowest = TRUE)
    y <- cut(y, breaks = c(minimum, thresh[i,], maximum), labels = c(-1, 0, 1), include.lowest = TRUE)
    list(x = x, y = y)
  }, simplify = FALSE)
           
  stat <- sapply(seq_len(nthresh), function(i) {
    Npm <- sum(tempp[[i]]$x == 1 & tempp[[i]]$y == -1)
    Nmp <- sum(tempp[[i]]$x == -1 & tempp[[i]]$y == 1)
    Npm + Nmp
  })
  id <- which.max(stat)[1]
  
  up <- thresh$V1[id]
  low <- thresh$V2[id]
  xdisc <- tempp[[id]]$x
  ydisc <- tempp[[id]]$y
  par <- data.frame(N = nx + ny, 
                    Kp = sum(xdisc == 1) + sum(ydisc == 1), 
                    Km = sum(xdisc == -1) + sum(ydisc == -1),
                    K0 = sum(xdisc == 0) + sum(ydisc == 0), 
                    Kgp = sum(xdisc == 1), 
                    Kgm = sum(xdisc == -1), 
                    Kg0 = sum(xdisc == 0),
                    Krp = sum(ydisc == 1), 
                    Krm = sum(ydisc == -1), 
                    Kr0 = sum(ydisc == 0),
                    E = (sum(xdisc == 1) * sum(ydisc == -1) + sum(xdisc == -1) * sum(ydisc == 1)) / nx,
                    x.hat = stat[id])
  
  pval <- pval_hy.test(par, paired = paired)
  
  rval <- list(thresh = c("up" = up, "low" = low), 
               x0 = x, y0 = y, xdisc = xdisc, ydisc = ydisc,
               tab = table(x = xdisc, y = ydisc),
               par = par, p.value = pval, 
               data.name = dname, method = method)
  class(rval) <- "hy.test"
  rval
}

hy.test.formula <- function (formula, data, subset, na.action, ...) {
  if (missing(formula) || (length(formula) != 3L)) 
    stop("'formula' missing or incorrect")
  oneSampleOrPaired <- FALSE
  if (length(attr(terms(formula[-2L]), "term.labels")) != 1L) 
    if (formula[[3]] == 1L) 
      oneSampleOrPaired <- TRUE
  else stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame()))) 
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  if (!oneSampleOrPaired) {
    g <- factor(mf[[-response]])
    if (nlevels(g) != 2L) 
      stop("grouping factor must have exactly 2 levels")
    DATA <- setNames(split(mf[[response]], g), c("x", "y"))
    y <- do.call("hy.test", c(DATA, list(...)))
  }
  else {
    respVar <- mf[[response]]
    if (inherits(respVar, "Pair")) {
      DATA <- list(x = respVar[, 1], y = respVar[, 2], 
                   paired = TRUE)
      y <- do.call("hy.test", c(DATA, list(...)))
    }
    else {
      stop("'formula' missing or incorrect")
    }
  }
  y$data.name <- DNAME
  y
}

pval_hy.test <- function(obj, paired = TRUE){
  
  N <- obj$N
  N2 <- N / 2
  Kp <- obj$Kp
  Km <- obj$Km
  K0 <- obj$K0
  x.hat <- obj$x.hat
  
  if(!paired){
    Kgp <- obj$Kgp; Kgm <- obj$Kgm; Kg0 <- obj$Kg0
    Krp <- obj$Krp; Krm <- obj$Krm; Kr0 <- obj$Kr0
    
    # n++, n+-, n--, n-+
    vett.n.pos.pos <- seq.int(0, min(Kgp, Krp))
    vett.n.neg.neg <- seq.int(0, min(Kgm, Krm))
    vett.n.pos.neg <- seq.int(0, min(Kgp, Krm))
    vett.n.neg.pos <- seq.int(0, min(Kgm, Krp))
    vett <- CJ("+-" = vett.n.pos.neg, "-+" = vett.n.neg.pos, "++" = vett.n.pos.pos, "--" = vett.n.neg.neg)
    vett$`x` <- (vett$`+-` + vett$`-+`)
    id.vett <- (vett$x >= x.hat)
    vett <- vett[id.vett, ]
    vett$`+0` <- Kgp - vett$`+-` - vett$`++`
    vett$`-0` <- Kgm - vett$`--` - vett$`-+`
    id.vett <- vett$`+0` >= 0 & 
      vett$`-0` >= 0 & 
      Kr0 >= vett$`+0` & 
      Krp - vett$`++` >= vett$`-+` & 
      Krm - vett$`+-` >= vett$`--` &
      Kr0 - vett$`+0` >= vett$`-0`
    vett <- vett[id.vett, ]
    rm(list = c("id.vett"))
    # all(vett$x >= x.hat) # ok
    
    due <- lchoose(Krp, vett$`++`) + lchoose(Krm, vett$`+-`) + lchoose(Kr0, vett$`+0`) - lchoose(N2, Kgp)
    tre <- lchoose(Krp - vett$`++`, vett$`-+`) + lchoose(Krm - vett$`+-`, vett$`--`) + lchoose(Kr0 - vett$`+0`, vett$`-0`) - lchoose(N2 - Kgp, Kgm)
    pval <- sum(exp(due + tre))
  }
  else{
    # Kr-, Kr+, Kr0, Kg+, Kg-, Kg0
    # max(0, Kp - N2) <= Krp <= min(Kp, N2)
    vett.Krp <- seq.int(max(0, Kp - N2), min(Kp, N2))
    # max(0, Km - N2) <= Krp <= min(Km, N2)
    vett.Krm <- seq.int(max(0, Km - N2), min(Km, N2))
    vett2 <- CJ("r+" = vett.Krp, "r-" = vett.Krm)
    vett2$`r0` <- - vett2$`r-` - vett2$`r+` + N2
    vett2$`g+` <- Kp - vett2$`r+`
    vett2$`g-` <- Km - vett2$`r-`
    vett2$`g0` <- - Km - Kp + vett2$`r-` + vett2$`r+` + N2
    id.vett <- vett2$`r0` >= 0 &  vett2$`g0` >= 0
    vett2 <- vett2[id.vett,]
    
    pval <- sapply(seq_len(nrow(vett2)), function(i){
      temp.vett <- vett2[i, ]
      
      # n++, n+-, n--, n-+
      vett.n.pos.pos <- seq.int(0, min(temp.vett$`g+`, temp.vett$`r+`))
      vett.n.neg.neg <- seq.int(0, min(temp.vett$`g-`, temp.vett$`r-`))
      vett.n.pos.neg <- seq.int(0, min(temp.vett$`g+`, temp.vett$`r-`))
      vett.n.neg.pos <- seq.int(0, min(temp.vett$`g-`, temp.vett$`r+`))
      vett <- CJ("+-" = vett.n.pos.neg, "-+" = vett.n.neg.pos, 
                 "++" = vett.n.pos.pos, "--" = vett.n.neg.neg)
      
      vett$`x` <- (vett$`+-` + vett$`-+`)
      id.vett <- (vett$x >= x.hat)
      vett <- vett[id.vett, ]
      
      vett$`+0` <- temp.vett$`g+` - vett$`+-` - vett$`++`
      vett$`-0` <- temp.vett$`g-` - vett$`--` - vett$`-+`
      id.vett <- vett$`+0` >= 0 & vett$`-0` >= 0
      vett <- vett[id.vett, ]
      
      id.vett <- (temp.vett$`r+` >= vett$`++`) & 
        (temp.vett$`r-` >= vett$`+-`) &
        (temp.vett$`r0` >= vett$`+0`) &
        (temp.vett$`r+` - vett$`++` >= vett$`-+`) &
        (temp.vett$`r-` - vett$`+-` >= vett$`--`) &
        (temp.vett$`r0` - vett$`+0` >= vett$`-0`)
      vett <- vett[id.vett, ]
      
      uno <- lchoose(Kp, temp.vett$`r+`) + 
        lchoose(Km, temp.vett$`r-`) + 
        lchoose(K0, temp.vett$`r0`) - 
        lchoose(N, N2)
      # uno <- test[cbind(Kp, temp.vett$`r+`)+1] + test[cbind(Km, temp.vett$`r-`)+1] + test[cbind(K0, temp.vett$`r0`)+1] - test[cbind(N, N2)+1]
      
      due <- lchoose(temp.vett$`r+`, vett$`++`) + 
        lchoose(temp.vett$`r-`, vett$`+-`) + 
        lchoose(temp.vett$`r0`, vett$`+0`) - 
        lchoose(N2, temp.vett$`g+`)
      # due <- test[cbind(temp.vett$`r+`, vett$`++`)+1] + test[cbind(temp.vett$`r-`, vett$`+-`)+1] + test[cbind(temp.vett$`r0`, vett$`+0`)+1] - test[cbind(N2, temp.vett$`g+`)+1]
      
      tre <- lchoose(temp.vett$`r+` - vett$`++`, vett$`-+`) + 
        lchoose(temp.vett$`r-` - vett$`+-`, vett$`--`) + 
        lchoose(temp.vett$`r0` - vett$`+0`, vett$`-0`) - 
        lchoose(N2 - temp.vett$`g+`, temp.vett$`g-`)
      # tre <- test[cbind(temp.vett$`r+` - vett$`++`, vett$`-+`)+1] + test[cbind(temp.vett$`r-` - vett$`+-`, vett$`--`)+1] + test[cbind(temp.vett$`r0` - vett$`+0`, vett$`-0`)+1] - test[cbind(N2 - temp.vett$`g+`, temp.vett$`g-`)+1]
      
      sum(exp(uno + due + tre))
    })
    pval <- sum(pval)
  }
  pval
}
  
print.hy.test <- function(x, digits = getOption("digits"), prefix = "\t", ...) {
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  cat("data:  ", x$data.name, "\n", sep = "")
  out <- character()
  if (!is.null(x$par$x.hat)) 
    out <- c(out, paste("hat(x)", "=", format(x$par$x.hat, digits = max(1L, digits - 2L))))
  if (!is.null(x$par$E)) 
    out <- c(out, paste("E[x]", "=", format(x$par$E, digits = max(1L, digits - 2L))))
  if (!is.null(x$p.value)) {
    fp <- format.pval(x$p.value, digits = max(1L, digits - 3L))
    out <- c(out, paste("p-value", if (startsWith(fp, "<")) fp else paste("=", fp)))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")

  if (!is.null(x$tab)) {
    cat("\n")
    print(x$tab, digits = digits, ...)
  }
  cat("\n")
  invisible(x)
}

plot.hy.test <- function(x, select = c("both", "x", "y"), ...){
  select <- match.arg(select)
  x0 <- x$x0
  y0 <- x$y0
  up <- x$thresh["up"]
  low <- x$thresh["low"]
  minimum <- floor(min(min(x0), min(y0)))
  maximum <- ceiling(max(max(x0), max(y0)))
  switch(select,
         "x" = plot(density(x0), ...),
         "y" = plot(density(y0), ...),
         "both" = {
           plot(density(x0), xlim = c(minimum, maximum), ...)
           lines(density(y0), lty = 2, ...)
         }
  )
  abline(v = c(up, low), col = 2, lty = 3, lwd = 1.2)
  
  invisible(x)
}




