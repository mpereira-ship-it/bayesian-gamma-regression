#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
packages <- c('ISLR','posterior','ggplot2','gamboostLSS','mboost')
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop('Missing R packages: ',paste(missing,collapse=', '))

for (d in c('frozen','outputs/chains','outputs/predictions','outputs/tables',
            'outputs/diagnostics','outputs/figures'))
  dir.create(d,recursive=TRUE,showWarnings=FALSE)
for (file in c('prepare_college','hb_lasso_engine','hb_ridge_engine','ss_ben_engine',
               'evaluation_functions','comparators','figures'))
  source(file.path('R',paste0(file,'.R')))

if (!file.exists('frozen/college_design.rds'))
  prepare_college('frozen/college_design.rds',seed=20261012L)
F <- readRDS('frozen/college_design.rds')
stopifnot(identical(F$data,'ISLR::College'),identical(as.integer(F$seed),20261012L),
          identical(F$X_train,F$Z_train),identical(F$X_test,F$Z_test),
          nrow(F$X_train)==621L,nrow(F$X_test)==156L,
          ncol(F$X_train)==18L,ncol(F$X_test)==18L)

engines <- list('HB-Lasso'=engine_lasso,'HB-Ridge'=engine_ridge,'SS-BEN'=engine_ssben)
paths <- c('HB-Lasso'='R/hb_lasso_engine.R',
           'HB-Ridge'='R/hb_ridge_engine.R','SS-BEN'='R/ss_ben_engine.R')
slugs <- c('HB-Lasso'='hblasso','HB-Ridge'='hbridge','SS-BEN'='ssben')
seeds <- list('HB-Lasso'=20262101:20262104,
              'HB-Ridge'=20262201:20262204,
              'SS-BEN'=20262301:20262304)
burns <- c('HB-Lasso'=10000L,'HB-Ridge'=10000L,'SS-BEN'=20000L)
iterations <- c('HB-Lasso'=60000L,'HB-Ridge'=60000L,'SS-BEN'=180000L)
design_md5 <- unname(tools::md5sum('frozen/college_design.rds'))

run_chain <- function(method,k) {
  iterations_m <- unname(iterations[[method]])
  burn_m <- unname(burns[[method]])
  path <- sprintf('outputs/chains/%s_chain%02d.rds',slugs[[method]],k)
  signature <- list(method=method,chain=k,seed=as.integer(seeds[[method]][k]),
     n_iter=iterations_m,burn=burn_m,thin=1L,tau='AM',
     design_md5=design_md5,engine_md5=unname(tools::md5sum(paths[[method]])))
  if (file.exists(path)) {
    obj <- readRDS(path)
    if (!identical(obj$signature,signature)) stop('Incompatible saved chain: ',path)
    return(obj$fit)
  }
  set.seed(signature$seed)
  args <- list(y=F$y_train,X=F$X_train,Z=F$Z_train,
               n_iter=iterations_m,n_burn=burn_m,n_thin=1L,
               s_b0=10,s_v0=10,param_gamma='AM')
  if (method=='SS-BEN')
    args <- c(args,list(lam1_b=1,lam2_b=.25,lam1_v=1,lam2_v=.25,
                        a_pi_b=1,b_pi_b=1,a_pi_v=1,b_pi_v=1))
  else
    args <- c(args,list(a_b=1,b_b=1,r_b=1,a_v=1,b_v=1,r_v=1))
  fit <- do.call(engines[[method]],args)
  stopifnot(nrow(fit$B)==iterations_m-burn_m,nrow(fit$V)==iterations_m-burn_m,
            ncol(fit$B)==18L,ncol(fit$V)==18L,
            all(is.finite(fit$B)),all(is.finite(fit$V)))
  if (method=='SS-BEN')
    stopifnot(ncol(fit$GB)==17L,ncol(fit$GV)==17L,
              all(fit$GB %in% c(0,1)),all(fit$GV %in% c(0,1)))
  saveRDS(list(signature=signature,fit=fit),path,compress='xz')
  fit
}

fits <- lapply(names(engines),function(m) lapply(seq_len(4L),function(k) run_chain(m,k)))
names(fits) <- names(engines)

preds <- list()
for (m in names(fits)) {
  cache <- file.path('outputs/predictions',paste0(slugs[[m]],'.rds'))
  chain_files <- sprintf('outputs/chains/%s_chain%02d.rds',slugs[[m]],1:4)
  sig <- list(design_md5=design_md5,chain_md5=unname(tools::md5sum(chain_files)))
  if (file.exists(cache)) {
    saved <- readRDS(cache)
    if (!identical(saved$signature,sig)) stop('Incompatible predictive cache: ',cache)
    preds[[m]] <- saved$prediction
  } else {
    B <- all_draws(fits[[m]],'B')
    V <- all_draws(fits[[m]],'V')
    preds[[m]] <- eval_mixture_gamma(F$y_test,F$X_test,B,V,max_quantile_draws=4000L)
    saveRDS(list(signature=sig,prediction=preds[[m]]),cache,compress='xz')
    rm(B,V)
  }
}

comparators <- run_comparators(F)
preds <- c(preds,comparators$predictions)
expected_methods <- c('HB-Lasso','HB-Ridge','SS-BEN','Gamma constant precision',
                      'Gamma varying precision','GammaLSS')
stopifnot(identical(names(preds),expected_methods))

selection <- selection_table(fits,F$column_names[-1L])
effects <- coefficient_table(fits,F$center,F$scale,F$column_names)
performance <- make_predictive_table(preds,F$y_test,F$Private_test)
diag <- do.call(rbind,lapply(names(fits),function(m) mcmc_diagnostics(fits[[m]],m)))
bootstrap <- bootstrap_nlpd(preds,F$Private_test,B=2000L,seed=20261013L)
utils::write.csv(selection,'outputs/tables/selection.csv',row.names=FALSE)
utils::write.csv(effects,'outputs/tables/effects.csv',row.names=FALSE)
utils::write.csv(performance,'outputs/tables/predictive_performance.csv',row.names=FALSE)
utils::write.csv(bootstrap,'outputs/tables/bootstrap_nlpd.csv',row.names=FALSE)
utils::write.csv(diag,'outputs/diagnostics/mcmc.csv',row.names=FALSE)
utils::write.csv(comparators$gammalss_cv,'outputs/diagnostics/gammalss_cv.csv',row.names=FALSE)
utils::write.csv(data.frame(mstop=comparators$gammalss_mstop),
                 'outputs/diagnostics/gammalss_mstop.csv',row.names=FALSE)
utils::write.csv(indicator_stability(fits[['SS-BEN']],F$column_names[-1L]),
                 'outputs/diagnostics/ssben_pip_by_chain.csv',row.names=FALSE)
for (m in names(preds)) {
  slug <- gsub('[^a-z0-9]+','_',tolower(m))
  utils::write.csv(data.frame(row_id=F$index_test,Private=F$Private_test,
                              y=F$y_test,preds[[m]]),
                   file.path('outputs/predictions',paste0(slug,'.csv')),row.names=FALSE)
}
make_college_figures(F,fits,selection,preds)
