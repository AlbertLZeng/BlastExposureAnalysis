smahal = function(z,X) {
  X<-as.matrix(X)
  n<-dim(X)[1]
  rownames(X)<-1:n
  k<-dim(X)[2]
  m<-sum(z)
  for(j in 1:k) {
    X[,j]<- rank(X[,j])
  }
  cv<-cov(X)
  vuntied<-var(1:n)
  rat<-sqrt(vuntied/diag(cv))
  cv<-diag(rat)%*%cv%*%diag(rat)
  out<-matrix(NA,m,n-m)
  Xc<-X[z==0,]
  Xt<-X[z==1,]
  rownames(out)<-rownames(X)[z==1]
  colnames(out)<-rownames(X)[z==0]
  library(MASS)
  icov<-ginv(cv)
  if(ncol(X)>1) {
    for(i in 1:m) {
      out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
    }
  }
  if(ncol(X)==1) {
    for(i in 1:m) {
      out[i,]=(Xc-Xt[i])^2
    }
  }
  out
}

addcaliper = function(dmat,z,logitp,calipersd=.5,penalty=1000) {
  sd.logitp=sqrt((sd(logitp[z==1])^2+sd(logitp[z==0])^2)/2)
  adif=abs(outer(logitp[z==1],logitp[z==0],"-"))
  adif=(adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
  dmat=dmat+adif*penalty
  dmat
}

fitmatchedregression <- function(outcome, propscore.model, matchvec, data, treatment = "surveyed") 
  {
  
  reg_formula <- update(propscore.model$formula, as.formula(paste0(outcome, " ~ + ", treatment, " + matchvec + .")))
  matched_reg_model <- lm(reg_formula, data=data)
  return(matched_reg_model)
}

fitpropensitymodel <- function(data, treatment, covariates) {
  covariate_formula <- paste(covariates, collapse = "+")
  formula <- as.formula(paste(treatment, "~", covariate_formula))
  
  propscore.model <- glm(formula, family = binomial, x=TRUE,y=TRUE, data=data)
  
  return(propscore.model)
}

preparematch <- function(data, propscore.model, treatment) {
  if(!is.null(propscore.model$na.action)) {
    data <- data[-propscore.model$na.action, ]
  }
  
  dmy <- dummyVars(propscore.model$formula, data=data)
  Xmat <- data.frame(predict(dmy, newdata=data))
  Xmatmahal <- Xmat
  
  data$logit.ps <- predict(propscore.model)
  logit.propscore <- predict(propscore.model)
  
  treated <- data[[treatment]]
  pooled.sd <- sqrt(var(logit.propscore[treated==1])/2+var(logit.propscore[treated==0])/2)
  
  min.treated <- min(logit.propscore[treated == 1])
  max.control <- max(logit.propscore[treated==0])
  
  which.remove <- which((logit.propscore>(max.control+0.5*pooled.sd)) | (logit.propscore<(min.treated-0.5*pooled.sd)))
  
  data.original <- data
  Xmat.original <- Xmat
  
  if(length(which.remove)>0) {
    data.trimmed <- data[-which.remove,]
    Xmat.trimmed <- Xmat[-which.remove,]
    Xmatmahal.trimmed <- Xmatmahal[-which.remove,]
    
    data.full <- rbind(data.trimmed, data.original[which.remove,])
    Xmat.full <- rbind(Xmat.trimmed, Xmat.original[which.remove,])
  }
  else {
    data.trimmed <- data
    Xmat.trimmed <- Xmat
    Xmatmahal.trimmed <- Xmatmahal
    data.full <- data
    Xmat.full <- Xmat
  }
  
  rownames(data.trimmed) <- seq(1, nrow(data.trimmed), 1)
  rownames(Xmat.trimmed) <- seq(1, nrow(data.trimmed), 1)
  rownames(Xmatmahal.trimmed) <- seq(1, nrow(data.trimmed), 1)
  
  rownames(data.full) <- seq(1, nrow(data.full), 1)
  rownames(Xmat.full) <- seq(1, nrow(data.full), 1)
  
  return(
    list(
      data = data.trimmed,
      Xmat = Xmat.trimmed,
      Xmatmahal = Xmatmahal.trimmed,
      data.full = data.full,
      Xmat.full = Xmat.full
    ))
}

performmatch <- function(data, Xmatmahal, treatment, caliper.sd = 0.5, controls.per.match = 1, add.penalty = FALSE) {
  Xmatmahal <- Xmatmahal[, colSums(Xmatmahal) > 0]
  
  distmat <- smahal(data[[treatment]], Xmatmahal)
  distmat <- addcaliper(distmat, data[[treatment]], data$logit.ps, calipersd = caliper.sd)
  
  if(add.penalty) {
  ## asymmetric penalty for Race, take out for other types of matching
  var = rank(data$Race)
  dif =outer(var[data[[treatment]]==1], var[data[[treatment]]==0],"-")/sd(var)
  distmat = distmat +(dif^2)*(dif<0)
  }
  rownames(distmat) <- rownames(data)[data[[treatment]]==1]
  colnames(distmat) <- rownames(data)[data[[treatment]]==0]
  
  matchvec <- pairmatch(distmat, controls = controls.per.match, data=data)
  
  treated.index <- rep(0, sum(data[[treatment]]==1))
  control.index.mat <- matrix(rep(0, controls.per.match*length(treated.index)), ncol=controls.per.match)
  matchedset.index <- substr(matchvec, start = 3, stop = 10)
  matchedset.index.numeric <- as.numeric(matchedset.index)
  
  valid.matchsets <- sort(unique(matchedset.index.numeric[!is.na(matchedset.index.numeric)]))
  
  treated.index <- rep(0, length(valid.matchsets))
  control.index.mat <- matrix(rep(0, controls.per.match * length(valid.matchsets)), ncol = controls.per.match)
  
  for(i in seq_along(valid.matchsets)) {
    matched.set.temp <- which(matchedset.index.numeric == valid.matchsets[i])
    treated.temp.index <- which(data[[treatment]][matched.set.temp] == 1)
    treated.index[i] <- matched.set.temp[treated.temp.index]
    control.index.mat[i,] <- matched.set.temp[-treated.temp.index]
  }
  
  return(list(
    matchvec = matchvec,
    treated.index = treated.index,
    control.index = control.index.mat
  ))
}

calculatebalance <- function(Xmat.full, data.full, treatment, treated.index, control.index, missing.mat = NULL) {
  if(!is.null(missing.mat)) {
    Xmat.without.missing <- Xmat.full
    for(i in 1:ncol(Xmat.full)) {
      Xmat.without.missing[missing.mat[,i]==1, i] <- NA
    }
  }
  else {
    Xmat.without.missing <- Xmat.full
  }
  
  treatedmat <- Xmat.without.missing[data.full[[treatment]]==1,]
  controlmat.before <- Xmat.without.missing[data.full[[treatment]]==0,]
  controlmean.before <- apply(controlmat.before, 2, mean, na.rm = TRUE)
  treatmean <- apply(treatedmat, 2, mean, na.rm = TRUE)
  treatvar <- apply(treatedmat, 2, var, na.rm=TRUE)
  controlvar <- apply(controlmat.before, 2, var, na.rm = TRUE)
  
  stand.diff.before <- (treatmean - controlmean.before) / sqrt((treatvar+controlvar)/2)
  
  treatmat.after <- Xmat.without.missing[treated.index,]
  controlmat.after <- Xmat.without.missing[control.index,]
  
  controlmean.after <- apply(controlmat.after, 2, mean, na.rm=TRUE)
  treatmean.after <- apply(treatmat.after, 2, mean, na.rm=TRUE)
  
  stand.diff.after <- (treatmean.after-controlmean.after) / sqrt((treatvar+controlvar)/2)
  
  return(
    list(
    stand.diff.before = stand.diff.before,
    stand.diff.after = stand.diff.after
  ))
}

graphloveplots <- function(stand.diff.before, stand.diff.after) {
  library(ggplot2)
  abs.stand.diff.before <- abs(stand.diff.before[-1])
  abs.stand.diff.after <- abs(stand.diff.after[-1])
  covariates <- names(stand.diff.before[-1])
  
  plot.dataframe <- data.frame(
    abs.stand.diff = c(abs.stand.diff.before, abs.stand.diff.after),
    covariates = rep(covariates,2),
    type = c(rep("Before", length(covariates)), rep("After", length(covariates)))
  )
  
  ggplot(plot.dataframe, aes(x = abs.stand.diff, y=covariates)) +
    geom_point(size = 5, aes(shape = factor(type))) +
    scale_shape_manual(values = c(4,1)) +
    geom_vline(xintercept = c(0.1,0.2), lty = 2) + 
    labs(shpae = "Type", x = "Absolute Standardized Difference", y = "Covariates")
}

  matchanalysis <- function(data, treatment, covariates, controls.per.match = 1, caliper.sd = 0.5, add.penalty=FALSE) {
  
  propscore.model <- fitpropensitymodel(data, treatment, covariates)
  prepared.data <- preparematch(data, propscore.model, treatment)
  
  matching.results <- performmatch(
    prepared.data$data,
    prepared.data$Xmatmahal,
    treatment,
    caliper.sd,
    controls.per.match,
    add.penalty
  )
  
  balance.stats <- calculatebalance(
    prepared.data$Xmat.full,
    prepared.data$data.full,
    treatment,
    matching.results$treated.index,
    matching.results$control.index
  )
  
  prepared.data$data$matchvec <- matching.results$matchvec
  
  return(
  list(
    propscore.model = propscore.model,
    data = prepared.data$data,
    matchvec = matching.results$matchvec,
    treated.index = matching.results$treated.index,
    control.index = matching.results$control.index,
    balance = balance.stats
  ))
  
}
  



