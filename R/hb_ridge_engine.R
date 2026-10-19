engine_ridge <- local({

  tam_bloque <- 5L

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

mwg_opt <- function(y, X, Z, n_iter, n_burn, n_thin,
                    s_b0, s_v0,
                    a_b, b_b, r_b,
                    a_v, b_v, r_v,
                    param_gamma = "AM") {

  p <- ncol(X)
  q <- ncol(Z)

  loglik_eta_delta <- function(eta, delta) {
    if (any(!is.finite(eta)) || any(!is.finite(delta))) return(-Inf)
    mu <- exp(pmin(eta, 700))
    fi <- exp(pmin(delta, 700))
    if (any(!is.finite(mu)) || any(!is.finite(fi)) ||
        any(mu <= 0) || any(fi <= 0)) return(-Inf)
    ans <- loglik_gamma_tau(y=y, mu=mu, fi=fi, param_gamma=param_gamma)
    if (!is.finite(ans)) -Inf else ans
  }

  ess_block <- function(x_cur, prior_sd, ll_cur, eval_ll,
                        max_bracket = 2000L) {
    if (length(x_cur) != length(prior_sd)) {
      stop("Dimension incompatible en ESS por bloques.")
    }
    if (any(!is.finite(prior_sd)) || any(prior_sd <= 0)) {
      stop("Escala gaussiana no valida en ESS por bloques.")
    }
    if (!is.finite(ll_cur)) stop("Likelihood actual no finita en ESS.")

    nu <- stats::rnorm(length(x_cur), 0, prior_sd)
    log_y <- ll_cur + log(stats::runif(1))

    theta <- stats::runif(1, 0, 2*pi)
    theta_min <- theta - 2*pi
    theta_max <- theta
    n_eval <- 0L

    repeat {
      n_eval <- n_eval + 1L
      x_prop <- x_cur * cos(theta) + nu * sin(theta)
      ev <- eval_ll(x_prop)

      if (is.finite(ev$ll) && ev$ll >= log_y) {
        return(list(x=x_prop, ll=ev$ll, eta=ev$eta,
                    delta=ev$delta, n_eval=n_eval))
      }

      if (theta < 0) theta_min <- theta else theta_max <- theta
      if (n_eval >= max_bracket) {
        stop("ESS por bloques excedio max_bracket.")
      }
      theta <- stats::runif(1, theta_min, theta_max)
    }
  }

  b <- c(log(mean(y)) + rnorm(1,0,0.10), rnorm(p-1L,0,0.05))
  v <- c(rnorm(1,0,0.10), rnorm(q-1L,0,0.05))
  k_b <- rep(1,p-1L)
  k_v <- rep(1,q-1L)
  ro_b <- 1
  ro_v <- 1

  eta <- drop(X %*% b)
  delta <- drop(Z %*% v)
  ll_cur <- loglik_eta_delta(eta,delta)
  for (ii in 0:100) {
    if (is.finite(ll_cur)) break
    b[1] <- b[1] - 0.25
    eta <- drop(X %*% b)
    ll_cur <- loglik_eta_delta(eta,delta)
  }
  if (!is.finite(ll_cur)) stop("No fue posible inicializar ESS por bloques.")

  shared_b <- if (q > 1L) 2:q else integer(0)
  shared_v <- if (q > 1L) 2:q else integer(0)
  shared_groups <- if (length(shared_b)) {
    split(seq_along(shared_b), ceiling(seq_along(shared_b)/tam_bloque))
  } else list()

  extra_b <- if (p > q) (q+1L):p else integer(0)
  extra_groups <- if (length(extra_b)) {
    split(extra_b, ceiling(seq_along(extra_b)/tam_bloque))
  } else list()

  bloques <- list(list(tipo="intercepto"))
  if (length(shared_groups)) {
    for (gg in shared_groups) {
      bloques[[length(bloques)+1L]] <- list(
        tipo="compartido",
        ib=shared_b[gg],
        iv=shared_v[gg]
      )
    }
  }
  if (length(extra_groups)) {
    for (ib in extra_groups) {
      bloques[[length(bloques)+1L]] <- list(tipo="solo_b",ib=ib)
    }
  }

  n_g <- floor((n_iter-n_burn)/n_thin)
  B <- matrix(NA_real_,n_g,p)
  V <- matrix(NA_real_,n_g,q)
  colnames(B) <- paste0("b",0:(p-1L))
  colnames(V) <- paste0("v",0:(q-1L))
  h <- 1L

  n_eval_total <- 0
  n_eval_max <- 0L
  n_updates <- 0L

  for (it in seq_len(n_iter)) {

    orden <- sample.int(length(bloques))

    for (kk in orden) {
      bl <- bloques[[kk]]

      if (bl$tipo == "intercepto") {
        x0 <- c(b[1],v[1])
        sd0 <- c(s_b0,s_v0)

        evfun <- function(xx) {
          eta_p <- eta + X[,1] * (xx[1]-b[1])
          del_p <- delta + Z[,1] * (xx[2]-v[1])
          list(ll=loglik_eta_delta(eta_p,del_p),
               eta=eta_p,delta=del_p)
        }

        es <- ess_block(x0,sd0,ll_cur,evfun)
        b[1] <- es$x[1]
        v[1] <- es$x[2]
        eta <- es$eta
        delta <- es$delta
        ll_cur <- es$ll

      } else if (bl$tipo == "compartido") {
        ib <- bl$ib
        iv <- bl$iv
        x0 <- c(b[ib],v[iv])
        sd0 <- c(
          1/sqrt(pmax(k_b[ib-1L],1e-12)),
          1/sqrt(pmax(k_v[iv-1L],1e-12))
        )
        nb <- length(ib)

        evfun <- function(xx) {
          db <- xx[seq_len(nb)] - b[ib]
          dv <- xx[nb + seq_along(iv)] - v[iv]
          eta_p <- eta + drop(X[,ib,drop=FALSE] %*% db)
          del_p <- delta + drop(Z[,iv,drop=FALSE] %*% dv)
          list(ll=loglik_eta_delta(eta_p,del_p),
               eta=eta_p,delta=del_p)
        }

        es <- ess_block(x0,sd0,ll_cur,evfun)
        b[ib] <- es$x[seq_len(nb)]
        v[iv] <- es$x[nb + seq_along(iv)]
        eta <- es$eta
        delta <- es$delta
        ll_cur <- es$ll

      } else if (bl$tipo == "solo_b") {
        ib <- bl$ib
        x0 <- b[ib]
        sd0 <- 1/sqrt(pmax(k_b[ib-1L],1e-12))

        evfun <- function(xx) {
          db <- xx - b[ib]
          eta_p <- eta + drop(X[,ib,drop=FALSE] %*% db)
          list(ll=loglik_eta_delta(eta_p,delta),
               eta=eta_p,delta=delta)
        }

        es <- ess_block(x0,sd0,ll_cur,evfun)
        b[ib] <- es$x
        eta <- es$eta
        ll_cur <- es$ll
      }

      n_eval_total <- n_eval_total + es$n_eval
      n_eval_max <- max(n_eval_max,es$n_eval)
      n_updates <- n_updates + 1L
    }

    k_b <- pmax(
      rgamma(p-1L,shape=r_b+0.5,
             rate=pmax(ro_b+b[-1]^2/2,1e-10)),
      1e-12
    )
    k_v <- pmax(
      rgamma(q-1L,shape=r_v+0.5,
             rate=pmax(ro_v+v[-1]^2/2,1e-10)),
      1e-12
    )
    ro_b <- max(
      rgamma(1,shape=a_b+(p-1L)*r_b,rate=b_b+sum(k_b)),
      1e-12
    )
    ro_v <- max(
      rgamma(1,shape=a_v+(q-1L)*r_v,rate=b_v+sum(k_v)),
      1e-12
    )

    if (it > n_burn && ((it-n_burn) %% n_thin == 0L)) {
      B[h,] <- b
      V[h,] <- v
      h <- h+1L
    }
  }

  list(
    B=B,V=V,
    acc_b=NA_real_,acc_v=NA_real_,
    acc_b_post=NA_real_,acc_v_post=NA_real_,
    acc_b0_post=NA_real_,acc_v0_post=NA_real_,
    acc_joint_post=NA_real_,
    acc_standardized_post=NA_real_,
    ess_lik_evals_mean=if(n_updates>0) n_eval_total/n_updates else NA_real_,
    ess_lik_evals_max=n_eval_max,
    n_block_updates=n_updates,
    algoritmo_v5="blocked_elliptical_slice_gibbs"
  )
}
  mwg_opt
})
