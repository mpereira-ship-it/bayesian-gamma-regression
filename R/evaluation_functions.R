logmeanexp <- function(v) {
  if (anyNA(v) || any(is.nan(v))) stop('Densidad predictiva con NA/NaN')
  m <- max(v)
  if (is.infinite(m) && m<0) return(-Inf)
  m + log(mean(exp(v-m)))
}

eval_mixture_gamma <- function(y,X,B,V,max_quantile_draws=4000L) {
  stopifnot(ncol(X)==ncol(B),ncol(X)==ncol(V),nrow(B)==nrow(V),
            length(y)==nrow(X),all(y>0),nrow(B)>0L)
  S <- nrow(B);n <- length(y)
  qprobs <- c(.025,.10,.25,.75,.90,.975)
  qnames <- c('q025','q10','q25','q75','q90','q975')
  out <- matrix(NA_real_,n,3L+length(qprobs))
  colnames(out) <- c('log_pred','mean','pit',qnames)
  ii <- unique(as.integer(round(seq.int(1L,S,length.out=min(S,max_quantile_draws)))))
  for(i in seq_len(n)) {
    mu <- exp(pmin(drop(B %*% X[i,]),700))
    phi <- exp(pmin(drop(V %*% X[i,]),700))
    rate <- phi/mu
    if(any(!is.finite(mu)) || any(!is.finite(phi)) || any(!is.finite(rate)) ||
       any(mu<=0) || any(phi<=0) || any(rate<=0))
      stop('Parametros posteriores predictivos no finitos, obs=',i)
    lp <- stats::dgamma(y[i],shape=phi,rate=rate,log=TRUE)
    ff <- stats::pgamma(y[i],shape=phi,rate=rate)
    if (anyNA(lp) || anyNA(ff) || any(!is.finite(ff)))
      stop('Densidad/PIT no evaluable, obs=',i)
    pa <- phi[ii];ra <- rate[ii]
    mix_cdf <- function(z) mean(stats::pgamma(z,shape=pa,rate=ra))
    mix_q <- function(prob) {
      hi <- max(stats::qgamma(prob,shape=pa,rate=ra))
      if(!is.finite(hi)||hi<=0) stop('Cuantil predictivo no finito, obs=',i)
      stats::uniroot(function(z) mix_cdf(z)-prob,
          lower=0,upper=hi*(1+1e-8),tol=1e-7)$root
    }
    out[i,] <- c(logmeanexp(lp),mean(mu),mean(ff),
                 vapply(qprobs,mix_q,numeric(1)))
  }
  df <- as.data.frame(out)
  df$n_draws_density <- S
  df$n_draws_quantiles <- length(ii)
  df
}

eval_gamma_plugin <- function(y,mu,phi) {
  stopifnot(length(y)==length(mu),length(y)==length(phi),
            all(is.finite(mu)),all(is.finite(phi)),all(mu>0),all(phi>0))
  r <- phi/mu
  ret <- data.frame(log_pred=stats::dgamma(y,shape=phi,rate=r,log=TRUE),
                    mean=mu,pit=stats::pgamma(y,shape=phi,rate=r))
  pp <- c(.025,.10,.25,.75,.90,.975)
  for(j in seq_along(pp)) ret[[c('q025','q10','q25','q75','q90','q975')[j]]] <-
    stats::qgamma(pp[j],shape=phi,rate=r)
  ret$n_draws_density <- NA_integer_;ret$n_draws_quantiles <- NA_integer_
  ret
}

prediction_metrics <- function(y,p,method,group='All') {
  stopifnot(length(y)==nrow(p),all(is.finite(p$mean)),all(is.finite(p$pit)),
            all(p$pit>=0 & p$pit<=1),all(is.finite(p$log_pred)))
  data.frame(method=method,Private=group,n=length(y),NLPD=-mean(p$log_pred),
   RMSE=sqrt(mean((y-p$mean)^2)),MAE=mean(abs(y-p$mean)),
   coverage50=mean(y>=p$q25 & y<=p$q75),
   coverage80=mean(y>=p$q10 & y<=p$q90),
   coverage95=mean(y>=p$q025 & y<=p$q975),
   width50=mean(p$q75-p$q25),width80=mean(p$q90-p$q10),
   width95=mean(p$q975-p$q025),PIT_mean=mean(p$pit),
   PIT_KS_D=unname(suppressWarnings(stats::ks.test(p$pit,'punif')$statistic)))
}

make_predictive_table <- function(preds,y,Private) {
  do.call(rbind,lapply(names(preds),function(m){
    P <- preds[[m]]
    do.call(rbind,c(list(prediction_metrics(y,P,m,'All')),
      lapply(c('No','Yes'),function(g){ix<-which(Private==g)
        prediction_metrics(y[ix],P[ix,,drop=FALSE],m,g)})))
  }))
}

all_draws <- function(chains,slot) do.call(rbind,lapply(chains,function(x)x[[slot]]))

selection_table <- function(fits,terms) {
  stopifnot(length(terms)==17L)
  o <- data.frame(term=terms,stringsAsFactors=FALSE)
  for(m in names(fits)) {
    sl <- if(m=='HB-Lasso')'HB_Lasso' else if(m=='HB-Ridge')'HB_Ridge' else 'SS_BEN'
    B <- all_draws(fits[[m]],'B')[,-1,drop=FALSE]
    V <- all_draws(fits[[m]],'V')[,-1,drop=FALSE]
    o[[paste0(sl,'_magnitude_mu')]] <- colMeans(abs(B)>.05)
    o[[paste0(sl,'_magnitude_phi')]] <- colMeans(abs(V)>.05)
    o[[paste0(sl,'_active_mu')]] <- colMeans(abs(B)>.05)>.5
    o[[paste0(sl,'_active_phi')]] <- colMeans(abs(V)>.05)>.5
    if(m=='SS-BEN') {
      pB <- colMeans(all_draws(fits[[m]],'GB'))
      pV <- colMeans(all_draws(fits[[m]],'GV'))
      o$SS_BEN_PIP_mu <- pB;o$SS_BEN_PIP_phi <- pV
      o$SS_BEN_included_mu <- pB>.5;o$SS_BEN_included_phi <- pV>.5
    }
  }
  o
}

coefficient_table <- function(fits,centers,scales,terms) {
  summarize_slot <- function(chains,slot,method,pr) {
    M <- all_draws(chains,slot)
    raw <- M
    raw[,-1L] <- sweep(M[,-1L,drop=FALSE],2L,scales,'/')
    raw[,1L] <- M[,1L]-drop(raw[,-1L,drop=FALSE]%*%centers)
    stat <- function(a) c(mean=mean(a),SD=stats::sd(a),
       CrI025=unname(stats::quantile(a,.025)),CrI975=unname(stats::quantile(a,.975)))
    out <- do.call(rbind,lapply(seq_len(ncol(M)),function(j){
      a<-stat(M[,j]);b<-stat(raw[,j])
      interp <- if (j==1L) rep(NA_real_,nrow(M)) else if(pr=='mu')
        exp(pmin(raw[,j],700)) else exp(pmin(-0.5*raw[,j],700))
      c(a,b,transformed_mean=if(j==1L)NA_real_ else mean(interp),
        transformed_CrI025=if(j==1L)NA_real_ else unname(quantile(interp,.025)),
        transformed_CrI975=if(j==1L)NA_real_ else unname(quantile(interp,.975)))
    }))
    data.frame(method=method,predictor=pr,term=terms,
       standardized_mean=out[,1],standardized_SD=out[,2],
       standardized_CrI025=out[,3],standardized_CrI975=out[,4],
       original_mean=out[,5],original_SD=out[,6],
       original_CrI025=out[,7],original_CrI975=out[,8],
       transformed_mean=out[,9],transformed_CrI025=out[,10],
       transformed_CrI975=out[,11],
       transformed_description=if(pr=='mu')'Mean ratio per original unit' else
          'CV ratio per original unit',stringsAsFactors=FALSE)
  }
  do.call(rbind,lapply(names(fits),function(m) {
    rbind(summarize_slot(fits[[m]],'B',m,'mu'),
          summarize_slot(fits[[m]],'V',m,'phi'))
  }))
}

mcmc_diagnostics <- function(chains,method) {
  if(!requireNamespace('posterior',quietly=TRUE))stop('Instalar posterior')
  get_stat <- function(slot) {
    arr <- simplify2array(lapply(chains,`[[`,slot))
    if(length(dim(arr))!=3L)stop('Dimension inesperada para ',slot)
    arr <- aperm(arr,c(1L,3L,2L))
    dimnames(arr) <- list(iteration=NULL,chain=paste0('chain',seq_along(chains)),
        variable=colnames(chains[[1]][[slot]]))
    draws <- posterior::as_draws_array(arr)
    stat <- posterior::summarise_draws(draws,posterior::rhat,
       posterior::ess_bulk,posterior::ess_tail,posterior::mcse_mean)
    dd <- as.data.frame(stat,check.names=FALSE)
    names(dd) <- sub('^posterior::','',names(dd))
    for(nm in c('rhat','ess_bulk','ess_tail','mcse_mean'))
      if(!nm%in%names(dd))stop('Falta diagnostico ',nm)
    data.frame(method=method,block=slot,dd,check.names=FALSE)
  }
  do.call(rbind,lapply(if(method=='SS-BEN')c('B','V','GB','GV') else c('B','V'),get_stat))
}

indicator_stability <- function(chains,terms) {
  one <- function(slot,pred) {
    ms <- lapply(chains,`[[`,slot)
    pip <- sapply(ms,colMeans)
    switching <- sapply(ms,function(A)colMeans(A[-1L,,drop=FALSE]!=A[-nrow(A),,drop=FALSE]))
    data.frame(term=rep(terms,times=length(ms)),predictor=pred,
      chain=rep(seq_along(ms),each=length(terms)),
      PIP=as.vector(pip),switch_fraction=as.vector(switching))
  }
  rbind(one('GB','mu'),one('GV','phi'))
}

bootstrap_nlpd <- function(preds,Private,B=2000L,seed=20261013L) {
  if(length(preds)<2L)return(data.frame())
  set.seed(seed)
  pairs <- combn(names(preds),2L,simplify=FALSE)
  idx <- split(seq_along(Private),Private)
  mat <- sapply(preds,function(p)-p$log_pred)
  do.call(rbind,lapply(pairs,function(pair){
    difference <- mat[,pair[1]]-mat[,pair[2]]
    boot <- replicate(B,{
      bs <- unlist(lapply(idx,function(u)sample(u,length(u),replace=TRUE)),use.names=FALSE)
      mean(difference[bs])
    })
    data.frame(method_A=pair[1],method_B=pair[2],
       delta_NLPD_A_minus_B=mean(difference),
       bootstrap_low=unname(quantile(boot,.025)),
       bootstrap_high=unname(quantile(boot,.975)))
  }))
}

cv_profile_table <- function(fits,F) {
  tr <- F$train_df
  raw <- matrix(1,2L,length(F$column_names),dimnames=list(c('Public','Private'),F$column_names))
  raw[,2L] <- c(0,1)
  for(j in 3L:ncol(raw))raw[,j] <- stats::median(tr[[colnames(raw)[j]]])
  X <- raw
  X[,-1L] <- sweep(sweep(raw[,-1L,drop=FALSE],2L,F$center,'-'),2L,F$scale,'/')
  do.call(rbind,lapply(names(fits),function(m){
    V <- all_draws(fits[[m]],'V')
    ratio <- exp(pmin(-.5*drop(V%*%(X[2,]-X[1,])),700))
    data.frame(method=m,contrast='Private Yes / No; other covariates at training medians',
     CV_ratio=mean(ratio),CrI025=unname(quantile(ratio,.025)),
     CrI975=unname(quantile(ratio,.975)))
  }))
}
