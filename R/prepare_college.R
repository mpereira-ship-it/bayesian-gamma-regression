prepare_college <- function(path = 'frozen/college_design.rds', seed = 20261012L) {
  if (file.exists(path)) stop('El diseno College ya existe: no sobrescribir.')
  if (!requireNamespace('ISLR', quietly=TRUE)) stop('Instalar paquete ISLR.')
  data('College', package='ISLR', envir=environment())
  dd <- as.data.frame(College)
  predictors <- c('Private','Apps','Accept','Enroll','Top10perc','Top25perc',
    'F.Undergrad','P.Undergrad','Room.Board','Books','Personal','PhD',
    'Terminal','S.F.Ratio','perc.alumni','Expend','Grad.Rate')
  stopifnot(nrow(dd)==777L, ncol(dd)==18L,
            identical(setdiff(names(dd),'Outstate'),predictors),
            is.numeric(dd$Outstate), all(is.finite(dd$Outstate)), all(dd$Outstate>0),
            !anyNA(dd[,c('Outstate',predictors)]))
  dd$Private <- factor(as.character(dd$Private), levels=c('No','Yes'))
  stopifnot(!anyNA(dd$Private), all(table(dd$Private)>1L))
  set.seed(seed)
  index_train <- sort(unlist(lapply(split(seq_len(nrow(dd)),dd$Private),
    function(ids) sample(ids, size=floor(.80*length(ids)))),use.names=FALSE))
  index_test <- setdiff(seq_len(nrow(dd)),index_train)
  fml <- stats::reformulate(predictors)
  raw <- stats::model.matrix(fml,data=dd,contrasts.arg=list(Private='contr.treatment'))
  expected <- c('(Intercept)','PrivateYes',predictors[-1L])
  stopifnot(identical(colnames(raw),expected),ncol(raw)==18L,
            all(raw[,1L]==1),length(index_train)+length(index_test)==777L)
  means <- colMeans(raw[index_train,-1L,drop=FALSE])
  sds <- apply(raw[index_train,-1L,drop=FALSE],2L,stats::sd)
  if (any(!is.finite(sds)) || any(sds<=0)) stop('Columna con SD cero en train.')
  standardize <- function(m) {
    x <- m
    x[,-1L] <- sweep(sweep(m[,-1L,drop=FALSE],2L,means,'-'),2L,sds,'/')
    storage.mode(x) <- 'double';x
  }
  Xtr <- standardize(raw[index_train,,drop=FALSE]);Xte <- standardize(raw[index_test,,drop=FALSE])
  stopifnot(identical(colnames(Xtr),expected),
            isTRUE(all.equal(Xtr, Xtr, tolerance=0)),
            qr(Xtr)$rank==ncol(Xtr), all(is.finite(Xtr)),all(is.finite(Xte)))
  obj <- list(data='ISLR::College', package_version=as.character(utils::packageVersion('ISLR')),
     seed=seed, predictors=predictors, column_names=expected,
     index_train=index_train,index_test=index_test,
     y_train=as.numeric(dd$Outstate[index_train]),y_test=as.numeric(dd$Outstate[index_test]),
     Private_train=as.character(dd$Private[index_train]),
     Private_test=as.character(dd$Private[index_test]),
     X_train=Xtr,Z_train=Xtr,X_test=Xte,Z_test=Xte,
     raw_train=raw[index_train,,drop=FALSE], raw_test=raw[index_test,,drop=FALSE],
     center=means,scale=sds,
     train_df=dd[index_train,c('Outstate',predictors),drop=FALSE],
     test_df=dd[index_test,c('Outstate',predictors),drop=FALSE])
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  saveRDS(obj,path,version=3)
  utils::write.csv(data.frame(row_id=seq_len(777L),
      set=ifelse(seq_len(777L)%in%index_train,'train','test'),
      Private=as.character(dd$Private)), 'frozen/split_manifest.csv',row.names=FALSE)
  utils::write.csv(data.frame(term=names(means),training_mean=means,training_sd=sds),
       'frozen/scaling_constants.csv',row.names=FALSE)
  writeLines(c(paste('source:',obj$data),paste('ISLR_version:',obj$package_version),
    paste('seed:',seed),paste('n_train:',length(index_train)),
    paste('n_test:',length(index_test)),capture.output(table(dd$Private))),
    'frozen/provenance.txt')
  invisible(obj)
}
