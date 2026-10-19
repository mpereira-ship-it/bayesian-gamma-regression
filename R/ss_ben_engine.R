engine_ssben <- local({

  tam_bloque <- 5L

  ridge_precond <- 0.10

tau_gamma <- function(param) {
  param <- match.arg(param, c("AM", "MD", "GM", "MO", "HM"))
  if (param == "AM") return(0)
  if (param == "MD") return(1 / 3)
  if (param == "GM") return(1 / 2)
  if (param %in% c("MO", "HM")) return(1)
  stop("Parametrizacion gamma no reconocida: ", param)
}

loglik_gamma_tau <- function(y, mu, fi, param_gamma = "AM") {
  if (any(y <= 0) || any(!is.finite(y))) {
    return(-Inf)
  }
  if (any(mu <= 0) || any(fi <= 0) ||
      any(!is.finite(mu)) || any(!is.finite(fi))) {
    return(-Inf)
  }
  tau <- tau_gamma(param_gamma)
  delta <- (sqrt(fi * (fi + 4 * tau)) + fi) / 2
  shape <- delta + tau
  rate  <- delta / mu
  if (any(!is.finite(delta)) || any(delta <= 0) ||
      any(!is.finite(shape)) || any(shape <= 0) ||
      any(!is.finite(rate))  || any(rate <= 0)) {
    return(-Inf)
  }
  sum(dgamma(y, shape = shape, rate = rate, log = TRUE))
}

lp_en_ss <- function(th_b, th_v, gb, gv,pi_b, pi_v,lam1_b, lam2_b,lam1_v, lam2_v,X, Z, y,s_b0, s_v0,param_gamma = "AM") {
  if (any(!is.finite(th_b)) || any(!is.finite(th_v))) return(-Inf)
  if (any(!(gb %in% c(0, 1))) || any(!(gv %in% c(0, 1)))) return(-Inf)
  if (!is.finite(pi_b) || !is.finite(pi_v)) return(-Inf)
  if (pi_b <= 0 || pi_b >= 1 || pi_v <= 0 || pi_v >= 1) return(-Inf)
  
  b <- c(th_b[1], gb * th_b[-1])
  v <- c(th_v[1], gv * th_v[-1])
  
  eta <- drop(X %*% b)
  del <- drop(Z %*% v)
  
  if (any(!is.finite(eta)) || any(!is.finite(del))) {
    return(-Inf)
  }
  
  mu <- exp(pmin(eta, 700))
  fi <- exp(pmin(del, 700))
  
  if (any(!is.finite(mu)) || any(!is.finite(fi)) || any(mu <= 0) || any(fi <= 0)) {
    return(-Inf)
  }
  
  
  
  
  
  if (any(y <= 0) || any(!is.finite(y))) {
    return(-Inf)
  }
  
  
  
  
  ll <- loglik_gamma_tau(y = y,mu = mu,fi = fi,param_gamma = param_gamma)
  
  if (!is.finite(ll)) {
    return(-Inf)
  }
  
  pb0 <- dnorm(th_b[1], 0, s_b0, log = TRUE)
  pv0 <- dnorm(th_v[1], 0, s_v0, log = TRUE)
  
  pb_slab <- -sum(lam1_b * abs(th_b[-1]) + lam2_b * th_b[-1]^2)
  pv_slab <- -sum(lam1_v * abs(th_v[-1]) + lam2_v * th_v[-1]^2)
  
  pg_b <- sum(gb * log(pi_b) + (1 - gb) * log(1 - pi_b))
  pg_v <- sum(gv * log(pi_v) + (1 - gv) * log(1 - pi_v))
  
  val <- ll + pb0 + pv0 + pb_slab + pv_slab + pg_b + pg_v
  if (!is.finite(val)) {
    return(-Inf)
  }
  
  val
}

crear_precondicionadores <- function(M, bloques, ridge = 0.10) {
  lapply(bloques, function(idx) {
    d <- length(idx)
    if (d == 1L) return(matrix(1, 1, 1))
    R <- suppressWarnings(stats::cor(M[, idx, drop = FALSE]))
    if (any(!is.finite(R))) R <- diag(d)
    Q <- tryCatch(
      solve(R + ridge * diag(d)),
      error = function(e) diag(d)
    )
    esc <- sqrt(pmax(diag(Q), 1e-10))
    C <- Q / outer(esc, esc)
    C <- (C + t(C)) / 2
    tryCatch(
      chol(C + 1e-8 * diag(d)),
      error = function(e) diag(d)
    )
  })
}

mwg_opt <- function(y, X, Z, n_iter, n_burn, n_thin,
                    s_b0, s_v0,
                    lam1_b, lam2_b,
                    lam1_v, lam2_v,
                    a_pi_b = 1, b_pi_b = 1,
                    a_pi_v = 1, b_pi_v = 1,
                    param_gamma = "AM") {
  p <- ncol(X); q <- ncol(Z); n_shared <- min(p-1L,q-1L)
  if(n_shared<1L)stop("Se requiere al menos un predictor compartido.")

  adapt_interval<-100L;adapt_cov_start<-1000L;adapt_shrink<-0.15;adapt_ridge<-0.05
  target_multi<-0.234;target_intercept<-0.35
  group_gibbs_every <- 20L
  inactive_refresh_every <- 5L

  block_diag2<-function(A,B){out<-matrix(0,nrow(A)+nrow(B),ncol(A)+ncol(B));out[seq_len(nrow(A)),seq_len(ncol(A))]<-A;ii<-nrow(A)+seq_len(nrow(B));jj<-ncol(A)+seq_len(ncol(B));out[ii,jj]<-B;out}
  safe_chol<-function(S){S<-(S+t(S))/2;d<-nrow(S);if(d==1L)return(matrix(sqrt(max(S[1,1],1e-8)),1,1));esc<-max(1e-8,stats::median(pmax(diag(S),1e-8)));for(jj in 0:8){rr<-(10^jj)*1e-10*esc;C<-tryCatch(chol(S+rr*diag(d)),error=function(e)NULL);if(!is.null(C))return(C)};diag(sqrt(pmax(diag(S),1e-6)),d)}
  init_adapt<-function(d,cov0=diag(d),init_scale=0.12)list(d=d,n=0L,mean=rep(0,d),M2=matrix(0,d,d),chol=safe_chol(cov0),log_scale=log(init_scale),acc_w=0L,prop_w=0L,upd=0L,acc_total=0L,prop_total=0L,acc_post=0L,prop_post=0L)
  observe_adapt<-function(a,x){x<-as.numeric(x);a$n<-a$n+1L;delta<-x-a$mean;a$mean<-a$mean+delta/a$n;delta2<-x-a$mean;a$M2<-a$M2+tcrossprod(delta,delta2);a}
  record_adapt<-function(a,ok,it){a$prop_total<-a$prop_total+1L;if(ok)a$acc_total<-a$acc_total+1L;if(it>n_burn){a$prop_post<-a$prop_post+1L;if(ok)a$acc_post<-a$acc_post+1L}else{a$prop_w<-a$prop_w+1L;if(ok)a$acc_w<-a$acc_w+1L};a}
  refresh_adapt<-function(a,target){a$upd<-a$upd+1L;if(a$prop_w>0L){gamma_ad<-0.40/sqrt(a$upd);ar<-a$acc_w/a$prop_w;a$log_scale<-pmin(log(5),pmax(log(1e-3),a$log_scale+gamma_ad*(ar-target)))};if(a$n>=adapt_cov_start&&a$n>2L){S<-a$M2/(a$n-1L);vv<-pmax(diag(S),1e-6);D<-sqrt(vv);C<-S/outer(D,D);C[!is.finite(C)]<-0;diag(C)<-1;C<-(1-adapt_shrink)*C+adapt_shrink*diag(a$d);Sreg<-(D%o%D)*C+adapt_ridge*stats::median(vv)*diag(a$d);a$chol<-safe_chol(Sreg)};a$acc_w<-0L;a$prop_w<-0L;a}
  rw_prop<-function(x,a){d<-length(x);as.numeric(x+exp(a$log_scale)*(2.38/sqrt(d))*drop(t(a$chol)%*%rnorm(d)))}
  acc_rate<-function(lst,post=TRUE){pr<-sum(vapply(lst,function(a)if(post)a$prop_post else a$prop_total,numeric(1)));ac<-sum(vapply(lst,function(a)if(post)a$acc_post else a$acc_total,numeric(1)));if(pr==0)NA_real_ else ac/pr}

  r_slab_en1 <- function(l1,l2) {
    l1 <- max(as.numeric(l1),1e-10); l2 <- max(as.numeric(l2),0)
    sg <- if (runif(1) < 0.5) -1 else 1
    if (l2 < 1e-10) return(sg * rexp(1, rate=l1))
    sd0 <- 1/sqrt(2*l2); mu0 <- -l1/(2*l2)
    p0 <- pnorm(0, mean=mu0, sd=sd0)
    u <- runif(1, min=p0, max=1)
    x <- qnorm(u, mean=mu0, sd=sd0)
    if (!is.finite(x) || x < 0) x <- 0
    sg*x
  }
  refresh_inactive <- function(theta,gamma,l1,l2) {
    idx <- which(gamma==0L)
    if(length(idx)) for(j in idx) theta[j+1L] <- r_slab_en1(l1[j],l2[j])
    theta
  }
  sample_logweights <- function(logw) {
    m <- max(logw); w <- exp(logw-m); w[!is.finite(w)] <- 0
    if(sum(w)<=0) return(which.max(logw))
    sample.int(length(w),1,prob=w)
  }

  shared0<-split(seq_len(n_shared),ceiling(seq_along(seq_len(n_shared))/tam_bloque))
  joint<-lapply(shared0,function(ii)list(base=ii,b=ii+1L,v=ii+1L))
  extra0<-if((p-1L)>n_shared)(n_shared+1L):(p-1L) else integer(0)
  extra0<-if(length(extra0))split(extra0,ceiling(seq_along(extra0)/tam_bloque)) else list()
  extra<-lapply(extra0,function(ii)list(base=ii,b=ii+1L))
  groups_b <- c(lapply(joint,`[[`,"base"),lapply(extra,`[[`,"base"))
  groups_v <- lapply(joint,`[[`,"base")
  cov_joint0<-lapply(joint,function(g){cb<-crossprod(crear_precondicionadores(X,list(g$b),ridge_precond)[[1]]);cv<-crossprod(crear_precondicionadores(Z,list(g$v),ridge_precond)[[1]]);block_diag2(cb,cv)})
  cov_extra0<-lapply(extra,function(g)crossprod(crear_precondicionadores(X,list(g$b),ridge_precond)[[1]]))

  th_b<-c(log(mean(y))+rnorm(1,0,0.10),rnorm(p-1L,0,0.05));th_v<-c(rnorm(1,0,0.10),rnorm(q-1L,0,0.05))
  gb<-rbinom(p-1L,1L,0.5);gv<-rbinom(q-1L,1L,0.5);pi_b<-0.5;pi_v<-0.5
  if(length(lam1_b)==1L)lam1_b<-rep(lam1_b,p-1L);if(length(lam2_b)==1L)lam2_b<-rep(lam2_b,p-1L)
  if(length(lam1_v)==1L)lam1_v<-rep(lam1_v,q-1L);if(length(lam2_v)==1L)lam2_v<-rep(lam2_v,q-1L)
  init_ok<-FALSE
  for(ii in 0:100){if(is.finite(lp_en_ss(th_b,th_v,gb,gv,pi_b,pi_v,lam1_b,lam2_b,lam1_v,lam2_v,X,Z,y,s_b0,s_v0,param_gamma))){init_ok<-TRUE;break};th_b[1]<-th_b[1]-0.25}
  if(!init_ok)stop("No fue posible encontrar un estado inicial finito.")

  ad_int<-init_adapt(2L,diag(2),0.10);ad_j<-Map(function(C)init_adapt(nrow(C),C,0.10),cov_joint0);ad_e<-Map(function(C)init_adapt(nrow(C),C,0.10),cov_extra0)
  group_change_b<-group_prop_b<-group_change_v<-group_prop_v<-0L
  group_change_b_post<-group_prop_b_post<-group_change_v_post<-group_prop_v_post<-0L
  n_g<-floor((n_iter-n_burn)/n_thin);B<-matrix(NA_real_,n_g,p);V<-matrix(NA_real_,n_g,q);GB<-matrix(NA_real_,n_g,p-1L);GV<-matrix(NA_real_,n_g,q-1L)
  colnames(B)<-paste0("b",0:(p-1L));colnames(V)<-paste0("v",0:(q-1L));colnames(GB)<-paste0("gb",seq_len(p-1L));colnames(GV)<-paste0("gv",seq_len(q-1L));h<-1L

  update_group_b <- function(gamma, idx, lp_current) {
    d <- length(idx); states <- as.matrix(expand.grid(rep(list(0:1),d)))
    logw <- rep(-Inf,nrow(states))
    for(rr in seq_len(nrow(states))){gg<-gamma;gg[idx]<-as.integer(states[rr,]);logw[rr]<-lp_en_ss(th_b,th_v,gg,gv,pi_b,pi_v,lam1_b,lam2_b,lam1_v,lam2_v,X,Z,y,s_b0,s_v0,param_gamma)}
    kk<-sample_logweights(logw);new<-gamma;new[idx]<-as.integer(states[kk,]);list(gamma=new,lp=logw[kk],changed=any(new[idx]!=gamma[idx]))
  }
  update_group_v <- function(gamma, idx, lp_current) {
    d <- length(idx); states <- as.matrix(expand.grid(rep(list(0:1),d)))
    logw <- rep(-Inf,nrow(states))
    for(rr in seq_len(nrow(states))){gg<-gamma;gg[idx]<-as.integer(states[rr,]);logw[rr]<-lp_en_ss(th_b,th_v,gb,gg,pi_b,pi_v,lam1_b,lam2_b,lam1_v,lam2_v,X,Z,y,s_b0,s_v0,param_gamma)}
    kk<-sample_logweights(logw);new<-gamma;new[idx]<-as.integer(states[kk,]);list(gamma=new,lp=logw[kk],changed=any(new[idx]!=gamma[idx]))
  }

  for(it in seq_len(n_iter)){
    lp_actual<-lp_en_ss(th_b,th_v,gb,gv,pi_b,pi_v,lam1_b,lam2_b,lam1_v,lam2_v,X,Z,y,s_b0,s_v0,param_gamma)
    if(!is.finite(lp_actual))stop("Estado MCMC actual no finito.")

    prop<-rw_prop(c(th_b[1],th_v[1]),ad_int);bp<-th_b;vp<-th_v;bp[1]<-prop[1];vp[1]<-prop[2]
    lpp<-lp_en_ss(bp,vp,gb,gv,pi_b,pi_v,lam1_b,lam2_b,lam1_v,lam2_v,X,Z,y,s_b0,s_v0,param_gamma)
    ok<-is.finite(lpp-lp_actual)&&log(runif(1))<(lpp-lp_actual);if(ok){th_b<-bp;th_v<-vp;lp_actual<-lpp};ad_int<-record_adapt(ad_int,ok,it)

    for(g in seq_along(joint)){ib<-joint[[g]]$b;iv<-joint[[g]]$v;prop<-rw_prop(c(th_b[ib],th_v[iv]),ad_j[[g]]);bp<-th_b;vp<-th_v;bp[ib]<-prop[seq_along(ib)];vp[iv]<-prop[length(ib)+seq_along(iv)]
      lpp<-lp_en_ss(bp,vp,gb,gv,pi_b,pi_v,lam1_b,lam2_b,lam1_v,lam2_v,X,Z,y,s_b0,s_v0,param_gamma);ok<-is.finite(lpp-lp_actual)&&log(runif(1))<(lpp-lp_actual);if(ok){th_b<-bp;th_v<-vp;lp_actual<-lpp};ad_j[[g]]<-record_adapt(ad_j[[g]],ok,it)}
    for(g in seq_along(extra)){ib<-extra[[g]]$b;bp<-th_b;bp[ib]<-rw_prop(th_b[ib],ad_e[[g]]);lpp<-lp_en_ss(bp,th_v,gb,gv,pi_b,pi_v,lam1_b,lam2_b,lam1_v,lam2_v,X,Z,y,s_b0,s_v0,param_gamma);ok<-is.finite(lpp-lp_actual)&&log(runif(1))<(lpp-lp_actual);if(ok){th_b<-bp;lp_actual<-lpp};ad_e[[g]]<-record_adapt(ad_e[[g]],ok,it)}

    for(j in seq_len(p-1L)){alt<-gb;alt[j]<-1L-gb[j];lalt<-lp_en_ss(th_b,th_v,alt,gv,pi_b,pi_v,lam1_b,lam2_b,lam1_v,lam2_v,X,Z,y,s_b0,s_v0,param_gamma);if(!is.finite(lalt))next
      if(gb[j]==0L){m<-max(lp_actual,lalt);pr1<-exp(lalt-m)/(exp(lp_actual-m)+exp(lalt-m))}else{m<-max(lalt,lp_actual);pr1<-exp(lp_actual-m)/(exp(lalt-m)+exp(lp_actual-m))};new<-rbinom(1,1,pr1);if(new!=gb[j]){gb[j]<-new;lp_actual<-lalt}}
    for(k in seq_len(q-1L)){alt<-gv;alt[k]<-1L-gv[k];lalt<-lp_en_ss(th_b,th_v,gb,alt,pi_b,pi_v,lam1_b,lam2_b,lam1_v,lam2_v,X,Z,y,s_b0,s_v0,param_gamma);if(!is.finite(lalt))next
      if(gv[k]==0L){m<-max(lp_actual,lalt);pr1<-exp(lalt-m)/(exp(lp_actual-m)+exp(lalt-m))}else{m<-max(lalt,lp_actual);pr1<-exp(lp_actual-m)/(exp(lalt-m)+exp(lp_actual-m))};new<-rbinom(1,1,pr1);if(new!=gv[k]){gv[k]<-new;lp_actual<-lalt}}

    if(it %% group_gibbs_every == 0L){
      for(idx in groups_b){old<-gb[idx];uu<-update_group_b(gb,idx,lp_actual);gb<-uu$gamma;lp_actual<-uu$lp;group_prop_b<-group_prop_b+1L;if(uu$changed)group_change_b<-group_change_b+1L;if(it>n_burn){group_prop_b_post<-group_prop_b_post+1L;if(uu$changed)group_change_b_post<-group_change_b_post+1L}}
      for(idx in groups_v){old<-gv[idx];uu<-update_group_v(gv,idx,lp_actual);gv<-uu$gamma;lp_actual<-uu$lp;group_prop_v<-group_prop_v+1L;if(uu$changed)group_change_v<-group_change_v+1L;if(it>n_burn){group_prop_v_post<-group_prop_v_post+1L;if(uu$changed)group_change_v_post<-group_change_v_post+1L}}
    }

    pi_b<-rbeta(1,a_pi_b+sum(gb),b_pi_b+length(gb)-sum(gb));pi_v<-rbeta(1,a_pi_v+sum(gv),b_pi_v+length(gv)-sum(gv))

    if(it %% inactive_refresh_every == 0L){th_b<-refresh_inactive(th_b,gb,lam1_b,lam2_b);th_v<-refresh_inactive(th_v,gv,lam1_v,lam2_v);lp_actual<-lp_en_ss(th_b,th_v,gb,gv,pi_b,pi_v,lam1_b,lam2_b,lam1_v,lam2_v,X,Z,y,s_b0,s_v0,param_gamma);if(!is.finite(lp_actual))stop("Refresco exacto de theta inactivo produjo estado no finito.")}

    if(it<=n_burn){ad_int<-observe_adapt(ad_int,c(th_b[1],th_v[1]));for(g in seq_along(joint)){ib<-joint[[g]]$b;iv<-joint[[g]]$v;ad_j[[g]]<-observe_adapt(ad_j[[g]],c(th_b[ib],th_v[iv]))};for(g in seq_along(extra)){ib<-extra[[g]]$b;ad_e[[g]]<-observe_adapt(ad_e[[g]],th_b[ib])};if(it%%adapt_interval==0L){ad_int<-refresh_adapt(ad_int,target_intercept);ad_j<-lapply(ad_j,refresh_adapt,target=target_multi);ad_e<-lapply(ad_e,refresh_adapt,target=target_multi)}}
    if(it>n_burn&&((it-n_burn)%%n_thin==0L)){B[h,]<-c(th_b[1],gb*th_b[-1]);V[h,]<-c(th_v[1],gv*th_v[-1]);GB[h,]<-gb;GV[h,]<-gv;h<-h+1L}
  }
  b_list<-c(ad_j,ad_e);v_list<-ad_j
  chb<-if(group_prop_b_post>0)group_change_b_post/group_prop_b_post else NA_real_;chv<-if(group_prop_v_post>0)group_change_v_post/group_prop_v_post else NA_real_
  list(B=B,V=V,GB=GB,GV=GV,acc_b=acc_rate(b_list,FALSE),acc_v=acc_rate(v_list,FALSE),acc_b_post=acc_rate(b_list,TRUE),acc_v_post=acc_rate(v_list,TRUE),
       acc_b0_post=if(ad_int$prop_post>0)ad_int$acc_post/ad_int$prop_post else NA_real_,acc_v0_post=if(ad_int$prop_post>0)ad_int$acc_post/ad_int$prop_post else NA_real_,
       acc_group_b_post=chb,acc_group_v_post=chv,group_change_b_post=chb,group_change_v_post=chv,
       algoritmo_v3="adaptive_joint_exact_group_gibbs")
}
  mwg_opt
})

