# Figures for manuscript: captions, axes and legend in English.
make_college_figures <- function(F,fit,selected,preds) {
 if(!requireNamespace('ggplot2',quietly=TRUE))stop('Install ggplot2')
 terms <- F$column_names[-1L]
 rows <- list()
 for(m in names(fit)) {
   prefix <- switch(m,'HB-Lasso'='HB_Lasso','HB-Ridge'='HB_Ridge','SS-BEN'='SS_BEN')
   for(pred in c('mu','phi')) {
      rows[[length(rows)+1L]] <- data.frame(term=terms,method=m,predictor=pred,
       criterion='Magnitude',probability=selected[[paste0(prefix,'_magnitude_',pred)]])
      if(m=='SS-BEN')rows[[length(rows)+1L]] <- data.frame(term=terms,method=m,
       predictor=pred,criterion='PIP',probability=selected[[paste0(prefix,'_PIP_',pred)]])
   }
 }
 dt<-do.call(rbind,rows)
 dt$term<-factor(dt$term,levels=rev(terms))
 dt$predictor <- factor(dt$predictor,levels=c('mu','phi'),
    labels=c('Central tendency (mu)','Precision (phi)'))
 fig <- ggplot2::ggplot(dt,ggplot2::aes(x=term,y=probability,shape=criterion,
                                     group=criterion)) +
   ggplot2::geom_point(position=ggplot2::position_dodge(width=.4),size=1.7) +
   ggplot2::geom_hline(yintercept=.5,linetype='dashed',linewidth=.3) +
   ggplot2::coord_flip() +
   ggplot2::facet_grid(predictor~method) +
   ggplot2::scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25))+
   ggplot2::labs(x=NULL,y='Posterior probability',shape='Criterion') +
   ggplot2::theme_bw(base_size=11)
 ggplot2::ggsave('outputs/figures/figure1_selection.pdf',fig,width=12,height=8)

 pitdf <- do.call(rbind,lapply(names(preds),function(m) {
   data.frame(model=m,pit=preds[[m]]$pit,Private=F$Private_test)
 }))
 f2 <- ggplot2::ggplot(pitdf,ggplot2::aes(x=pit,colour=Private)) +
   ggplot2::stat_ecdf(geom='step',linewidth=.65)+
   ggplot2::geom_abline(intercept=0,slope=1,linetype='dashed',colour='grey35') +
   ggplot2::coord_equal(xlim=c(0,1),ylim=c(0,1))+
   ggplot2::facet_wrap(~model,ncol=3)+
   ggplot2::labs(x='Probability integral transform (PIT)',
                 y='Empirical cumulative distribution',colour='Private')+
   ggplot2::theme_bw(base_size=10)
 ggplot2::ggsave('outputs/figures/figure2_pit_calibration.pdf',f2,width=11,height=7)

  invisible(NULL)
}
