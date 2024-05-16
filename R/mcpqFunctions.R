
#################################################
###  Point estimate of quantile and variances ###
#################################################
mcpqest <- function(y, f, event = NULL, Right.Censored = FALSE, p = 0.5, ...){
  aargs <- list(...)
  if(is.null(aargs$bw.selec)) {
    bw.selec <- "plug-in"
  } else {
    bw.selec <- aargs$bw.selec
  }
  if(length(y) != length(f)){
    stop("y and f must have the same length")
  }
  if(!is.numeric(y)){
    stop("y must be numeric")
  }
  ff <- droplevels(as.factor(f))
  ni <- tapply(ff, ff, length)
  ngroups <- length(ni)
  groupLabel <- names(ni)
  if(Right.Censored && !is.null(event)){
    if(length(base::unique(event)) != 2){
      stop("event must be a binary variable")
    }
    if(length(y) != length(event)){
      stop("y and event must have the same length")
    }
    Dat <- data.frame(f = f, time = y, event = event)
    sm <- summary(survfit(Surv(time, event) ~ f, data = Dat))
    sdata <- data.frame(strata = sm$strata, time = sm$time,
                        prob = sm$surv, std.err = sm$std.err)
    survquant <- sapply(1:ngroups, function(i){
      stratai <- (levels(sdata$strata))[i]
      datai <- sdata[sdata$strata == stratai, ]
      quanti <- datai[datai$prob <= (1-p), "time"]
      if(length(quanti) == 0){
        stop(paste(p), "th quantile of group ", paste(stratai),
             " does not exist \n")
      }else if(sum(is.na(quanti)) == length(quanti)){
        stop(paste(p), "th quantile of group ", paste(stratai),
             " does not exist \n")
      }
      quant <- head(quanti[!(is.na(quanti))],1)
      std.err <- datai[datai$time == quant, "std.err"]
      if(is.nan(std.err) | is.na(std.err)){
        stop("In group ", paste(stratai), " std.err of survival 
             function at ", paste(p), "th quantile does not exist \n")
      }
      return(c(quant, std.err))
    })
    quantileEST <- survquant[1,]
    survden <- sapply(1:ngroups, function(i){
      datai <- Dat[Dat$f == groupLabel[i], ]
      deni <- presmooth(times = time, status = event, dataset = datai, 
                        estimand = "f", bw.selec = bw.selec, 
                        x.est = quantileEST[i])$estimate
      return(deni)
    })
    varEST <- (survquant[2,]/survden)^2
  }
  if(is.null(event) && !Right.Censored ){
    Dat <- data.frame(f = f, y = y)
    quant <- tapply(y, f, quantile, p, na.rm=TRUE,names=FALSE)
    kd <- sapply(1:ngroups, function(i){
      den.grp <- density(Dat[f == groupLabel[i],2,drop = TRUE], na.rm = TRUE)
      datfram.grp <- data.frame(x = den.grp$x, y = den.grp$y)
      ix <- which.min(abs(den.grp$x - quant[i]))
      return(datfram.grp[ix, "y"])
    })
    quantileEST <- quant
    varEST <- (p*(1 - p))/(ni*(kd^2))
  }
  return(list(quantileEST = quantileEST, varEST = varEST, n = ni))
}

###################################
##  mcpq function for difference ##
###################################
mcpqdci <- function(y, f, event = NULL, Right.Censored = FALSE, p = 0.5, 
                     conf.level = 0.95, type = "Dunnett", base = 1, 
                     cmat = NULL, ...){
  if(!is.numeric(p) || (length(p) != 1 | p >= 1 | p < 0)) { 
    stop("p should be a single numeric value between 0 and 1")
  } 
  if(!is.numeric(conf.level)||(length(conf.level)!=1|
                               conf.level>=1|conf.level<0.5)){
    stop("conf.level should be a single numeric value between 0.5 and 1")
  }
  aargs <- list(...)
  if(is.null(aargs$dist)) {
    dist <- "MVT"
  } else {
    dist <- aargs$dist
  }
  if(is.null(aargs$bw.selec)) {
    bw.selec <- "plug-in"
  } else {
    bw.selec <- aargs$bw.selec
  }
  if(Right.Censored){
    if(is.null(event)){
      stop("If Right.Censored = TRUE, event argument must be provided")
    }
  }
  if(!is.null(event)){
    if(!Right.Censored){
      stop("If event argument is provided, Right.Censored 
           argument must be TRUE")
    }
    if(any(y < 0)){
      stop("y must be non-negative for right-censored data")
    }
  }
  if(length(f) != length(y))
    stop("Argument y and f must have the same length")
  if(!is.numeric(y)){
    stop("y must be numeric")
  }
  ff <- droplevels(as.factor(f))
  ni <- tapply(ff, ff, length)
  ngroups <- length(ni)
  groupLabel <- names(ni)
  if(!(ngroups > 1)){
    stop("Number of groups should be at least 2")
  }
  if (any(ni  < 2)) {
    stop("the number of observations in each group should be at least 2")
  }
  if (is.null(cmat)) {
    if (type == "Dunnett") {
      cmat   <- contrMat(n = ni, type = type, base = base)
    } else {
      cmat   <- contrMat(n = ni, type = type)
    }
  } else {
    if (!is.matrix(cmat)|| (ncol(cmat) != ngroups)) {
      stop("cmat must be a matrix with number of columns = number of groups")
    }
    if(sum(apply(cmat, 1, sum))!=0){
      stop("each row of cmat must give a sum of 0")
    }
    colnames(cmat) <- groupLabel
    tmp <- apply(abs(cmat), 1, sum)
    m <- dim(cmat)[1]
    if(all(tmp == 2)){
      mlabel <- sapply(1:m, function(x){
        i1 <-  which(cmat[x,] == 1)
        i2 <-  which(cmat[x,] == -1)
        groupNames <- groupLabel[c(i1,i2)]
        paste0(groupNames[1], " - ", groupNames[2])
      })
      rownames(cmat) <- mlabel
    } else {
      rownames(cmat) <- paste0("C", 1:m)
    }
  }
  mcpq.estimates <- mcpqest(y, f, event, Right.Censored, p, ...)
  quantileEST <- mcpq.estimates$quantileEST
  varEST   <- mcpq.estimates$varEST
  SIGMA <- diag(varEST)
  covMat <- cmat%*%SIGMA%*%t(cmat)
  invD <- solve(sqrt(diag(diag(covMat))))
  corMat <- invD%*%covMat%*%invD
  contrastEST <- cmat%*%quantileEST  
  std.err <- as.matrix(sqrt(diag(covMat)))
  if (dist == "MVN"){
    eQ <- qmvnorm(conf.level, corr = corMat, tail = "both.tails")$quantile
  } else if (dist == "MVT"){
    eQ <- qmvt(conf.level, corr = corMat, 
               df = sum(ni)-ngroups, tail = "both.tails")$quantile
  }
  lCI <- contrastEST - eQ*std.err
  uCI <- contrastEST + eQ*std.err
  CI  <- cbind(lCI, uCI)
  colnames(CI) <- c("lower", "upper")
  attr(cmat, "type") <- NULL
  attr(cmat, "class") <- NULL
  return(list(cmat = cmat, conf.level = conf.level, estimate = contrastEST, 
              std.err = std.err, conf.int = CI))
}



###################################
##  mcpq.rci function for ratio  ##
###################################

mcpqrci <- function(y, f, event = NULL, Right.Censored = FALSE, p = 0.5, 
                     conf.level = 0.95, type = "Dunnett", base = 1, 
                     Num.cmat = NULL, Den.cmat = NULL, 
                     method = c("Wald", "Fieller"), ...){
  if(!is.numeric(p) || (length(p) != 1 | p >= 1 | p <= 0)) { 
    stop("p should be a single numeric value between 0 and 1")
  } 
  if(!is.numeric(conf.level)||(length(conf.level)!=1|
                               conf.level>=1|conf.level<0.5)){
    stop("conf.level should be a single numeric value between 0.5 and 1")
  }
  method <- match.arg(method)
  aargs <- list(...)
  if(is.null(aargs$dist)) {
    dist <- "MVT"
  } else {
    dist <- aargs$dist
  }
  if(is.null(aargs$bw.selec)) {
    bw.selec <- "plug-in"
  } else {
    bw.selec <- aargs$bw.selec
  }
  if(Right.Censored){
    if(is.null(event)){
      stop("If Right.Censored = TRUE event argument must be provided")
    }
  }
  if(!is.null(event)){
    if(!Right.Censored){
      stop("If Event argument is provided, Right.Censored 
           argument must be TRUE")
    }
    if(any(y < 0)){
      stop("y must be non-negative for right-censored data")
    }
  }
  if(length(y) != length(f))
    stop("Argument y and f must have the same length")
  if(!is.numeric(y)){
    stop("y must be numeric")
  }
  ff <- droplevels(as.factor(f))
  ni <- tapply(ff, ff, length)
  ngroups <- length(ni)
  groupLabel <- names(ni)
  if(!(ngroups > 1)){
    stop("Number of groups should be at least 2")
  }
  if (any(ni < 2)) {
    stop("the number of observations in each group should be at least 2")
  }
  if (!is.null(Num.cmat) || !is.null(Num.cmat)){
    if (!is.null(Num.cmat) && is.null(Den.cmat)) {
      stop("Num.cmat is specified, but Den.cmat is missing")
    }
    if (is.null(Num.cmat) && !is.null(Den.cmat)) {
      stop("Den.cmat is specified, but Num.cmat is missing")
    }
    if (!is.null(Num.cmat) && !is.null(Den.cmat)) {
      if(!is.matrix(Num.cmat) || !is.matrix(Den.cmat)){
        stop("Num.cmat and Den.cmat must be matrices")
      }
      if (dim(Den.cmat)[1] != dim(Num.cmat)[1]) {
        stop("number of rows in Num.cmat and Den.cmat are not the same")
      }
      if (ncol(Den.cmat) != ngroups && ncol(Num.cmat) != ngroups) {
        stop("number of columns in Num.cmat or Den.cmat 
             is not the same as number of groups")
      }
      if(dim(Den.cmat)[1] <= 1){
        stop("number of rows in Num.cmat or Den.cmat must be atleast 2")
      }
      NC0 <- apply(Num.cmat, 1, function(x){
        all(x == 0)
      })
      DC0 <- apply(Den.cmat, 1, function(x){
        all(x == 0)
      })
      if (any(c(NC0, DC0))) {
        message("At least one row of the numerator or denominator contrast 
                matrices is a vector with all components equal to zero")
      }
      Num.C <- Num.cmat
      Den.C <- Den.cmat
      Userdefined <- TRUE
      if (is.null(rownames(Num.C)) && is.null(rownames(Den.C))) {
        compnames <- paste("C", 1:nrow(Num.C), sep = "")
      } else {
        if (any(rownames(Num.C) != rownames(Den.C))) {
          compnames <- paste(rownames(Num.C), rownames(Den.C), sep = "/")
        } else {
          compnames <- rownames(Num.C)
        }
      }
      }
  } else {
    Userdefined <- FALSE
    Cmat <- contrMatRatio(n = ni, type = type, base = base)
    Num.C <- Cmat$numC
    Den.C <- Cmat$denC
    compnames <- Cmat$rnames
  }
  mcpq.estimates <- mcpqest(y, f, event, Right.Censored = Right.Censored, 
                             p, ...)
  quantileEST <- mcpq.estimates$quantileEST
  varEST   <- mcpq.estimates$varEST
  ratioEST <- (Num.C%*%quantileEST)/(Den.C%*%quantileEST)
  Dmatplugin <- apply(Den.C, 2, function(x){ x*ratioEST }) - Num.C
  SIGMA <- diag(varEST)
  covMat <- (Dmatplugin)%*%SIGMA%*%(t(Dmatplugin))
  invD <- solve(sqrt(diag(diag(covMat))))
  corMat <- invD%*%covMat%*%invD
  if (dist == "MVN"){
    eQ <- qmvnorm(conf.level, corr = corMat, tail = "both.tails")$quantile
  } else if(dist == "MVT"){
    eQ <- qmvt(conf.level, corr = corMat,
               df = sum(ni)-ngroups, tail = "both.tails")$quantile
  }
  switch(method,
         Wald ={
           m <- dim(Num.C)[1]
           std.err  <- sapply(1:m, function(j){
             aa <- t(Num.C[j,])%*%SIGMA%*%((Num.C[j,]))
             bb <- (t(Den.C[j,])%*%SIGMA%*%((Den.C[j,])))*(ratioEST[j])^2
             ab <- (t(Num.C[j,])%*%SIGMA%*%((Den.C[j,])))*(ratioEST[j])*2
             return((1/(Den.C[j,]%*%quantileEST))*sqrt((aa+bb-ab)))
           })
           std.err <- as.matrix(std.err)
           rownames(std.err) <- compnames
           lCI <- ratioEST - eQ*std.err
           uCI <- ratioEST + eQ*std.err
           CI <- cbind(lCI, uCI)
         }, 
         Fieller = {
           m <- dim(Num.C)[1]
           G <- sapply(1:m, function(j){
             (eQ)^2*(t(Den.C[j,])%*%SIGMA%*%((Den.C[j,]
                                              )))/((Den.C[j,]%*%quantileEST)^2)
           })
           if(any(G >= 1)){
             message("Fieller method resulted in unbounded intervals,\n",
                     "consider changing the method to Wald")
             return(list(conf.level = conf.level, estimate = ratioEST ))
           }
           CI <- matrix(unlist(lapply(1:m, function(j){
             V11 <- t(Num.C[j,])%*%SIGMA%*%((Num.C[j,]))
             V22 <- t(Den.C[j,])%*%SIGMA%*%((Den.C[j,]))
             V12 <- t(Num.C[j,])%*%SIGMA%*%((Den.C[j,]))
             VV12 <- ratioEST[j]   - (V12/V22)*G[j]
             VV11 <- ratioEST[j]^2 - (V11/V22)*G[j]
             lCI <- (1/(1-G[j]))*(VV12 - sqrt((VV12)^2-(1-G[j])*(VV11)))
             uCI <- (1/(1-G[j]))*(VV12 + sqrt((VV12)^2-(1-G[j])*(VV11)))
             return(c(lCI, uCI))
           })), nrow = m, byrow = TRUE, ncol = 2)
           rownames(CI) <- compnames
         })
  if (Userdefined) {
    colnames(Num.C) <- groupLabel
    colnames(Den.C) <- groupLabel
    rownames(Num.C) <- compnames
    rownames(Den.C) <- compnames
    cmat <- Num.C+Den.C
    tmp <- apply(cmat, 1, sum)
    if(all(tmp == 2)){
      m <- dim(Den.C)[1]
      mlabel <- sapply(1:m, function(x){
        i1 <-  which(Num.cmat[x,] == 1)
        i2 <-  which(Den.cmat[x,] == 1)
        groupNames <- groupLabel[c(i1,i2)]
        paste0(groupNames[1], "/", groupNames[2])
      })
      if(method == "Wald"){
        rownames(std.err) <- mlabel
      }
      rownames(ratioEST) <- mlabel
      rownames(CI) <- mlabel
    } else {
      rownames(ratioEST) <- compnames
      rownames(CI) <- compnames
    }
  }
  colnames(CI) <- c("lower", "upper")
  if(method == "Wald"){
    return(list(Num.Constrast = Num.C, Den.Contrast = Den.C,
                conf.level = conf.level, estimate = ratioEST,
                std.err = std.err, conf.int = CI))
  } else {
    return(list(Num.Constrast = Num.C, Den.Contrast = Den.C,
                conf.level = conf.level,
                estimate = ratioEST, conf.int = CI))
  }
}


