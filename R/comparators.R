
nll_gamma <- function(par,X,Z,y) {
  nb <- ncol(X);nz <- ncol(Z)
  eta <- drop(X %*% par[seq_len(nb)])
  xi <- drop(Z %*% par[nb+seq_len(nz)])
  if(any(!is.finite(eta)) || any(!is.finite(xi)) ||
      any(eta > 700) || any(eta < -700) || any(xi > 700) || any(xi < -700))
    return(1e100)
  mu <- exp(eta); phi <- exp(xi)
  ll <- suppressWarnings(stats::dgamma(y,shape=phi,rate=phi/mu,log=TRUE))
  if(any(!is.finite(ll)))return(1e100)
  -mean(ll)
}

fit_constant_gamma <- function(F) {
  df <- data.frame(y=F$y_train,F$X_train[,-1L,drop=FALSE],check.names=FALSE)
  fm <- stats::reformulate(colnames(F$X_train)[-1L],response='y')
  mod <- stats::glm(fm,data=df,family=stats::Gamma(link='log'),
                    control=stats::glm.control(maxit=100L))
  if(!isTRUE(mod$converged) || any(!is.finite(coef(mod))))stop('GLM Gamma constante no converge')
  mu <- as.numeric(stats::predict(mod,newdata=df,type='response'))
  fn <- function(log_phi) -mean(stats::dgamma(F$y_train,shape=exp(log_phi),
                            rate=exp(log_phi)/mu,log=TRUE))
  opt <- stats::optimize(fn,c(log(1e-3),log(1e5)))
  if(!is.finite(opt$objective))stop('Gamma constante: forma ML invalida')
  list(beta=unname(stats::coef(mod)),phi=exp(opt$minimum),
       train_nlpd=opt$objective,glm=mod,convergence='OK')
}

grad_gamma <- function(par,X,Z,y) {
  nb <- ncol(X); nz <- ncol(Z)
  eta <- drop(X %*% par[seq_len(nb)])
  xi <- drop(Z %*% par[nb+seq_len(nz)])
  mu <- exp(eta); phi <- exp(xi)
  if (any(!is.finite(mu)) || any(!is.finite(phi)) ||
      any(mu<=0) || any(phi<=0)) stop('Gradiente Gamma: parametros no finitos')
  a <- phi * (y/mu-1)
  b <- phi * (log(phi)+1-log(mu)+log(y)-y/mu-digamma(phi))
  ans <- -c(drop(crossprod(X,a)),drop(crossprod(Z,b)))/length(y)
  if (any(!is.finite(ans))) stop('Gradiente Gamma: resultado no finito')
  ans
}

fit_variable_gamma <- function(F,constant) {
  X<-F$X_train;Z<-F$Z_train;y<-F$y_train
  nb<-ncol(X); nz<-ncol(Z)
  p0<-c(constant$beta,log(constant$phi),rep(0,nz-1L))
  low<-c(rep(-10,nb),rep(-8,nz));up<-c(rep(15,nb),rep(8,nz))
  if(any(p0<=low)|any(p0>=up))stop('MLE variable: inicio en cota')
  obj<-function(x)nll_gamma(x,X,Z,y)
  score<-function(x)grad_gamma(x,X,Z,y)
  o0 <- obj(p0)

  ans<-stats::optim(p0,fn=obj,method='L-BFGS-B',lower=low,upper=up,
                   control=list(maxit=4000L,factr=1e7,pgtol=1e-7,
                                ndeps=rep(1e-5,length(p0))))
  stage1_nll <- obj(ans$par)
  stage1_grad <- max(abs(score(ans$par)))
  refine_status <- 'NOT_NEEDED'

  if (is.finite(stage1_nll) && is.finite(stage1_grad) &&
      stage1_nll <= o0+1e-7 && stage1_grad>1e-3) {
    refine_status <- 'ATTEMPTED'
    refined <- tryCatch(stats::optim(ans$par,fn=obj,gr=score,
                  method='L-BFGS-B',lower=low,upper=up,
                  control=list(maxit=8000L,factr=1e4,pgtol=1e-9)),
                  error=function(e)e)
    if (inherits(refined,'error')) {
      refine_status <- paste0('ERROR: ',conditionMessage(refined))
    } else {
      new_score <- tryCatch(max(abs(score(refined$par))),
                             error=function(e) Inf)
      if (is.finite(refined$value) && is.finite(new_score) &&
          refined$value <= stage1_nll+1e-10 &&
          new_score < stage1_grad && refined$convergence==0L) {
        ans <- refined
        refine_status <- 'ACCEPTED'
      } else refine_status <- sprintf('NO_IMPROVEMENT: code=%s; score=%.6g',
                                      refined$convergence,new_score)
    }
  }
  final <- obj(ans$par)
  if (!is.finite(final) || !is.finite(o0) || final>o0+1e-7)
    stop('MLE Gamma variable: no se conserva el maximo del modelo anidado')
  beta<-ans$par[seq_len(nb)];gamma<-ans$par[nb+seq_len(nz)]
  grad <- score(ans$par)
  max_grad <- max(abs(grad))
  at_bound<-any(abs(ans$par-low)<1e-5|abs(ans$par-up)<1e-5)
  status <- if(ans$convergence==0L && !at_bound && max_grad<=1e-3)
    'OK_GRADIENTE_Y_CONVERGENCIA' else 'REVISAR_OPTIMIZACION'
  list(beta=beta,gamma=gamma,nll_initial=o0,nll_train=final,
       convergence=ans$convergence,max_abs_gradient=max_grad,
       at_bound=at_bound,status=status,message=ans$message,
       reduction_NLL_total=length(y)*(o0-final),
       stage1_nll=stage1_nll,stage1_max_gradient=stage1_grad,
       refinement=refine_status)
}


fit_boosted_gamma <- function(F, phi_initial) {
  if (!requireNamespace('gamboostLSS', quietly=TRUE) ||
      !requireNamespace('mboost', quietly=TRUE)) stop('Packages gamboostLSS and mboost are required.')
  df <- data.frame(y=F$y_train, F$X_train[,-1L,drop=FALSE], check.names=FALSE)
  fml <- stats::reformulate(colnames(F$X_train)[-1L], response='y')
  mu_initial <- mean(F$y_train)
  stopifnot(is.finite(mu_initial), mu_initial>0, is.finite(phi_initial), phi_initial>0)

  literal_family <- substitute(
    gamboostLSS::GammaLSS(mu=MU_CONST, sigma=PHI_CONST, stabilization='L2'),
    list(MU_CONST=as.numeric(mu_initial), PHI_CONST=as.numeric(phi_initial)))
  literal_call <- substitute(
    gamboostLSS::glmboostLSS(formula=fml, data=df, families=FAMILY_LITERAL,
      method='noncyclic', control=mboost::boost_control(mstop=2000L)),
    list(FAMILY_LITERAL=literal_family))
  captured <- character()
  capture_warning <- function(expr) withCallingHandlers(expr, warning=function(w) {
    captured <<- c(captured, conditionMessage(w))
    invokeRestart('muffleWarning')
  })
  fit <- capture_warning(eval(literal_call, envir=environment()))
  saved_call <- attr(fit,'call')
  if (is.null(saved_call[['families']]) ||
      !inherits(eval(saved_call[['families']], envir=baseenv()),'families'))
    stop('GammaLSS family cannot be reconstructed during cross-validation.')

  set.seed(20261015L)
  foldid <- integer(nrow(df))
  for (g in c('No','Yes')) {
    u <- which(F$Private_train==g)
    foldid[u] <- sample(rep(seq_len(5L),length.out=length(u)))
  }
  folds <- matrix(1,nrow=nrow(df),ncol=5L)
  for (k in seq_len(5L)) folds[foldid==k,k] <- 0
  storage.mode(folds) <- 'double'
  grid <- seq.int(20L,2000L,by=20L)
  cv <- capture_warning(mboost::cvrisk(fit, folds=folds, grid=grid, papply=lapply))
  risk <- unclass(as.matrix(cv))
  native_grid <- suppressWarnings(as.integer(attr(cv,'mstop')))
  if (!identical(native_grid,grid)) stop('Unexpected GammaLSS cross-validation grid.')
  if (identical(as.integer(dim(risk)),c(length(grid),5L))) {
    risk <- risk
  } else if (identical(as.integer(dim(risk)),c(5L,length(grid)))) {
    risk <- t(risk)
  } else stop('Unexpected GammaLSS cross-validation dimensions.')
  if (any(!is.finite(risk)) ||
      any(grepl('fold(s) encountered an error',captured,fixed=TRUE)))
    stop('GammaLSS cross-validation failed.')
  best <- as.integer(mboost::mstop(cv))
  if (length(best)!=1L || !is.finite(best) || !(best %in% grid))
    stop('GammaLSS selected an invalid boosting iteration.')
  if (best==max(grid)) stop('GammaLSS optimum is on the upper grid boundary.')
  fit <- fit[best]
  newdata <- data.frame(y=F$y_test, F$X_test[,-1L,drop=FALSE],check.names=FALSE)
  pred <- stats::predict(fit,newdata=newdata,type='response')
  if (!all(c('mu','sigma') %in% names(pred)))
    stop('GammaLSS prediction must contain mu and sigma.')
  cv_table <- data.frame(mstop=grid,risk,check.names=FALSE)
  names(cv_table)[-1L] <- paste0('fold_',seq_len(5L))
  cv_table$mean_cv_risk <- rowMeans(risk)
  list(prediction=eval_gamma_plugin(F$y_test,as.vector(pred$mu),as.vector(pred$sigma)),
       mstop=best,cv_table=cv_table)
}

run_comparators <- function(F) {
  constant <- fit_constant_gamma(F)
  varying <- fit_variable_gamma(F,constant)
  if (!identical(varying$status,'OK_GRADIENTE_Y_CONVERGENCIA'))
    stop('Gamma varying-precision MLE did not meet the original convergence criteria.')
  mu_const <- exp(drop(F$X_test %*% constant$beta))
  mu_vary <- exp(drop(F$X_test %*% varying$beta))
  phi_vary <- exp(drop(F$Z_test %*% varying$gamma))
  boosted <- fit_boosted_gamma(F,constant$phi)
  list(predictions=list(
    'Gamma constant precision'=eval_gamma_plugin(F$y_test,mu_const,
                                                 rep(constant$phi,length(mu_const))),
    'Gamma varying precision'=eval_gamma_plugin(F$y_test,mu_vary,phi_vary),
    'GammaLSS'=boosted$prediction),
    gammalss_mstop=boosted$mstop, gammalss_cv=boosted$cv_table)
}
