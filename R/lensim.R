#' Construct an inverse age-length key
#'
#' Constructs the probability distribution of length conditional on age from
#' a growth model and a coefficient of variation in length-at-age.
#'
#' @param params Growth parameters coercible to an `FLPar`. Parameters required
#'   by `model` must be named; the default von Bertalanffy model uses `linf`,
#'   `k`, and `t0`.
#' @param model Growth function or model object used to predict mean
#'   length-at-age. Defaults to `vonbert`.
#' @param age Numeric vector of ages.
#' @param cv Numeric coefficient of variation in length-at-age. Defaults to
#'   `0.1`.
#' @param lmin Numeric lower limit of the length grid. Defaults to `5`.
#' @param lmax Numeric multiplier applied to `linf` when calculating the upper
#'   length limit. Defaults to `1.2`.
#' @param bin Numeric length-bin width. Defaults to `1`.
#' @param max Numeric upper limit of the length grid. Defaults to
#'   `ceiling(linf * lmax)`.
#' @param reflen Optional reference length used to define a constant standard
#'   deviation as `cv * reflen`. When `NULL`, the standard deviation is
#'   `cv` times predicted length-at-age.
#'
#' @return An `FLPar` containing an age-by-length inverse age-length key. Rows
#'   represent ages, columns represent length classes, and each age row sums
#'   to one.
#'
#' @examples
#' \dontrun{
#' ialk <- iALK(
#'   params = FLPar(linf = 45, k = 0.4, t0 = -0.3),
#'   age = 0:5,
#'   lmin = 5,
#'   lmax = 1.2,
#'   bin = 1
#' )
#' }
#'
#' @export

iALK <- function(params, model=vonbert, age, cv=0.1,lmin=5, lmax=1.2, bin=1,
                 max=ceiling(linf * lmax), reflen=NULL) {
  
  linf <- c(params['linf'])
  
  # FOR each age
  bins <- seq(lmin, max, bin)
  
  # METHOD
  if(isS4(model))
    len <- do.call(model, list(age=age,params=params))
  else {
    lparams <- as(FLPar(params), "list")
    len <- do.call(model, c(list(age=age),
                            lparams[names(lparams) %in% names(formals(model))]))
  }
  
  if(is.null(reflen)) {
    sd <- abs(len * cv)
  } else {
    sd <- reflen * cv
  }
  
  probs <- Map(function(x, y) {
    p <- c(pnorm(1, x, y),
           dnorm(bins[-c(1, length(bins))], x, y),
           pnorm(bins[length(bins)], x, y, lower.tail=FALSE))
    return(p / sum(p))
  }, x=len, y=sd)
  
  res <- do.call(rbind, probs)
  
  alk <- FLPar(array(res, dim=c(length(age), length(bins), 1)),
               dimnames=list(age=age, len=bins, iter=1), units="")
  
  return(alk)
} 


# Internal null-default operator (base R >= 4.4.0 provides its own).
# Not exported.
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Length grid shared by iALK(), build_gear() and all downstream helpers
#'
#' Lower bin limits `seq(lmin, ceiling(linf * lmax_mult), bin)`, exactly as
#' [iALK()] builds them; the last bin is the plus group
#' \eqn{[l_{max}, \infty)}.
#'
#' Background: before this helper, `sel_la()` (exclusive upper limit) and
#' `iALK()` (plus-group upper limit) built their grids differently, and
#' every helper that rebuilt an `iALK()` (`plot_condALK()`,
#' `plot_condALK_bias()`, `lfd.sim(timing = )`, `sel_a_seasonal_avg()`,
#' `wt_a_seasonal()`, `plot_sel_a_seasonal()`) used its own default
#' `lmax_mult = 1.2`. With `build_gear(lmax_mult = 1.3)` those helpers built
#' a shorter grid than the gear's (e.g. 63 vs 69 bins for Linf = 55.7), and
#' `ialk %*% len_mid` / `condition_alk()` no longer matched. They now all
#' take the grid from `gear$grid` (via [gear_ialk()] or an explicit `max`).
#'
#' @param linf Asymptotic length.
#' @param lmin,lmax_mult,bin Grid settings.
#' @return Numeric vector of lower bin limits.
#' @export
len_grid <- function(linf, lmin = 5, lmax_mult = 1.2, bin = 1) {
  seq(lmin, ceiling(linf * lmax_mult), bin)
}


#' Biological inverse ALK on a gear's own length grid
#'
#' Builds \eqn{P(l \mid a, t)} with exactly the grid, CV and timing stored in
#' `gear$grid`, so its columns always match `gear$sel_len` and
#' `gear$condALK`. The upper limit is passed to [iALK()] explicitly via
#' `max`, so `iALK()`'s own `lmax` default can no longer diverge.
#'
#' @param gear A gear object from [build_gear()] (needs `lhpar`, `grid`).
#' @param age Ages. Default: rows of `gear$condALK`.
#' @param timing Within-year time; default `gear$grid$timing`.
#' @param params `FLPar` with `linf`, `k`, `t0`; default `gear$lhpar`.
#' @param model,reflen Passed to [iALK()].
#' @return Age x length matrix (rows sum to one).
#' @export
gear_ialk <- function(gear, age = as.numeric(rownames(gear$condALK)),
                      timing = gear$grid$timing %||% 0,
                      params = gear$lhpar, model = vonbert, reflen = NULL) {
  g <- gear$grid
  if (is.null(g)) stop("Gear has no $grid; rebuild it with build_gear().")
  linf <- c(params["linf"])
  lens <- len_grid(linf, g$lmin, g$lmax_mult, g$bin)
  ialk <- iALK(
    params = c(linf = linf, k = c(params["k"]), t0 = c(params["t0"]) - timing),
    model = model, age = age, cv = g$cv, lmin = g$lmin, bin = g$bin,
    max = max(lens), reflen = reflen
  )
  m <- c(ialk); dim(m) <- dim(ialk)[1:2]
  dimnames(m) <- list(age = dimnames(ialk)$age, len = dimnames(ialk)$len)
  if (!is.null(gear$condALK) && !identical(colnames(m), colnames(gear$condALK)))
    stop("Rebuilt iALK grid does not match gear$condALK.")
  m
}


# }}}

#' Convert an inverse age-length key to an age-length key
#'
#' Combines numbers-at-age with an inverse age-length key and conditions the
#' result on length to obtain age proportions within each length class.
#'
#' @param N_a Numeric vector of numbers-at-age for one sampling event.
#' @param iALK Inverse age-length key returned by [iALK()].
#'
#' @return An `FLPar` containing age proportions by length class.
#'
#' @examples
#' \dontrun{
#' alk <- ALK(N_a = c(100, 80, 50, 20), iALK = ialk)
#' }
#'
#' @export
ALK <- function(N_a,iALK){
  alk = iALK
  alk[] = N_a
  alk = alk*iALK
  alksum =apply(alk,2,sum)
  for(i in 1:dim(alk)[1]){
    alk[i,]=alk[i,]/an(alksum)
  }   
  return(alk)
}


# }}}

#' Sample annual age-length keys by length class
#'
#' Generates multinomial age samples within length classes from annual
#' age-length keys. The number sampled in each length class is limited by
#' `nbin` and by the corresponding length-frequency sample.
#'
#' @param lfds An `FLQuant` containing length-frequency data by year and
#'   iteration.
#' @param alks An `FLPars` object containing annual age-length keys, typically
#'   returned by [ALKs()].
#' @param nbin Maximum number of age observations sampled per length class.
#'   Defaults to `20`.
#' @param n.sample Optional total length-sample size used to rescale `lfds`.
#'   A value of `1` uses the supplied frequencies without rescaling.
#'
#' @return An `FLPars` object containing sampled annual age-length keys.
#'
#' @examples
#' \dontrun{
#' sampled_alks <- alk.sample(lfds, alks, nbin = 20)
#' }
#'
#' @export

alk.sample <- function(lfds,alks,nbin = 20,n.sample=1){
  res = alks
  if(n.sample>1){
    lfdn = lfds%/%apply(lfds,2:6,sum)*n.sample
  } else {
    lfdn = lfds
  }
  for(i in seq(dims(lfds)$iter)){
    for(y in seq(dims(lfds)$year)){
      nL = pmin(lfdn[,y],nbin )  # check len samples
      for(l in seq(dim(alks[[1]])[2])){
        if(nL[l]>1){  
          res[[y]][,l][] = apply(rmultinom(nL[l], 1, prob = c(alks[[y]][,l])),1,sum)
        } else {
          res[[y]][,l][] = 0}
        
      }
    }
  }
  res
}  


# }}}

#' Construct annual age-length keys
#'
#' Combines annual numbers-at-age with an inverse age-length key and returns
#' age proportions within each length class for every year and iteration.
#'
#' @param object An `FLQuant` containing numbers-at-age by year and iteration.
#' @param iALK Inverse age-length key returned by [iALK()].
#'
#' @return An `FLPars` object with one age-length key per year.
#'
#' @examples
#' \dontrun{
#' alks <- ALKs(stock.n(stk), ialk)
#' }
#'
#' @export
ALKs <- function(object,iALK){
  it = dim(object)[6]
  nyr= dim(object)[2]
  year = (dimnames(object)$year)
  alks = FLPars(lapply(year,function(x){
    alk = propagate(iALK,it)
    alk[] = object[,ac(x)]
    for(i in seq(it)){
      iter(alk,i) = iter(alk,i) *iALK
      iter(alk,i) =  sweep(iter(alk,i),2,apply(iter(alk,i),2,sum),"/")
    }
    
    return(alk)
  }))
  names(alks) = ac(year)
  return(alks)
}


# }}}

#' Convert length frequencies to numbers-at-age
#'
#' Applies annual age-length keys to length-frequency data to estimate
#' numbers-at-age by year and iteration.
#'
#' @param lfds An `FLQuant` containing numbers-at-length.
#' @param alks An annual `FLPars` collection returned by [ALKs()], or a single
#'   `FLPar` applied to every year.
#'
#' @return An `FLQuant` containing estimated numbers-at-age.
#'
#' @examples
#' \dontrun{
#' numbers_at_age <- applyALK(lfds, alks)
#' }
#'
#' @export

applyALK <- function(lfds,alks){
  yr = dimnames(lfds)$year
  year = an(yr)
  if(class(alks)=="FLPar"){
    alks = FLPars(lapply(yr,function(x){
      alks
    }))
    names(alks) = yr
  }
  age = ac(dimnames(alks[[1]])$age)
  
  if(any(dimnames(lfds)$year !=  names(alks)))
    stop("ALKs list must have the same years as length data")
  its = dims(lfds)$iter
  
  res <- FLQuant(NA, units="1000",
                 dimnames=list(age=age, year=year,iter=1:its))
  
  for(i in seq(its)){
    for(y in 1:length(year)){
      mf = (model.frame(iter(alks[[y]],i))[,seq(dim(alks[[y]])[1])])
      mf = mf/pmax(apply(mf,1,sum),0.01)
      iter(res,i)[,y] = apply(t(mf[,seq(dim(alks[[y]])[1])])%*%c(iter(lfds[,y],i)),1,sum)
    }
  }
  return(res)
}


#' Fold length selectivity into an inverse ALK
#'
#' Reweights an inverse age-length key by selectivity-at-length and normalises
#' each age row to obtain the conditional distribution of length among fish
#' caught by a gear.
#'
#' @param ialk_mat Numeric matrix with ages in rows and length classes in
#'   columns. Each age row should sum to one.
#' @param sel_len_vec Named numeric vector of selectivity-at-length. Values are
#'   reordered to match `colnames(ialk_mat)`.
#'
#' @return A numeric matrix with the same dimensions as `ialk_mat`. Rows sum
#'   to one, except ages with zero selected probability, which are returned as
#'   zero rows.
#'
#' @export
condition_alk <- function(ialk_mat, sel_len_vec) {
  
  sel_len_vec <- sel_len_vec[colnames(ialk_mat)]
  
  weighted <- sweep(ialk_mat, 2, sel_len_vec, "*")
  row_tot <- rowSums(weighted)
  cond <- weighted
  valid <- row_tot > 0
  cond[valid, ] <- weighted[valid, , drop = FALSE] / row_tot[valid]
  cond[!valid, ] <- 0
  cond
}


#' Simulate length-frequency samples from numbers-at-age
#'
#' Projects numbers-at-age through an inverse age-length key at one or more
#' within-year sampling times and draws a multinomial length sample. When a
#' gear is supplied, the inverse key is conditioned on its
#' selectivity-at-length.
#'
#' @param N_a An `FLQuant` containing numbers-at-age or catch-at-age.
#' @param params Growth parameters coercible to an `FLPar`; the default growth
#'   model requires `linf`, `k`, and `t0`.
#' @param model Growth function or model object. Defaults to `vonbert`.
#' @param ess Integer effective sample size. It is divided equally, after
#'   rounding, among the values in `timing`. Defaults to `250`.
#' @param timing Numeric vector of sampling times expressed as fractions of a
#'   year. Growth at time `t` is evaluated using `t0 - t`. Defaults to monthly
#'   sampling times from `0` to `11/12`.
#' @param unit Character label for the length unit. Defaults to `"cm"`.
#' @param scale Logical; if `TRUE`, rescale each year and iteration to the
#'   corresponding total in `N_a`.
#' @param reflen,bin,cv,lmin,lmax Arguments passed to [iALK()] (`lmax` is
#'   the multiplier on `linf`). When `gear` is supplied, `bin`, `cv`, `lmin`
#'   and `lmax` default to `gear$grid`.
#' @param gear Optional gear object as returned by [build_gear()]. If
#'   supplied, `gear$sel_len` is incorporated through [condition_alk()].
#' @param u Optional list of common-random-number arrays, one per timing
#'   event, as returned by [rUnif_len()]. If `NULL`, multinomial samples are
#'   drawn from the current random-number stream.
#'
#' @return An `FLQuant` containing sampled length frequencies by year and
#'   iteration.
#'
#' @details
#' This is a single-stage length sampler. In contrast, [lfd.sim()] first
#' samples catch-at-age at low effective sample size and then conditions the
#' length draw on that sampled age composition.
#'
#' @examples
#' \dontrun{
#' lfd_survey <- len.sim(stock.n(om)[, "2024"], params = lhpar, ess = 300,
#'                        timing = 0.5)
#' u_len <- rUnif_len(trawl, timing = seq(0, 11/12, 1/12),
#'                     years = an(dimnames(N_a)$year), nsim = dims(N_a)$iter,
#'                     seed = 456)
#' lfd_survey_g <- len.sim(N_a, params = lhpar, ess = trawl$ess_len,
#'                          timing = seq(0, 11/12, 1/12),
#'                          gear = trawl, u = u_len)
#' }
#'
#' @export
len.sim <- function(N_a, params, model = vonbert, ess = 250,
                    timing = seq(0, 11/12, 1/12), unit = "cm", scale = TRUE,
                    reflen = NULL, bin = gear$grid$bin %||% 1,
                    cv = gear$grid$cv %||% 0.1, lmin = gear$grid$lmin %||% 5,
                    lmax = gear$grid$lmax_mult %||% 1.2,
                    gear = NULL, u = NULL) {
  
  gp <- c(params)
  age <- an(dimnames(N_a)$age)
  its <- dims(N_a)$iter
  yrs <- dimnames(N_a)$year
  ess_t <- round(ess / length(timing))
  
  out <- NULL   # built on the first timing event, once we know the len grid
  
  for (t in seq_along(timing)) {
    
    ialk <- iALK(
      params = c(linf = gp[["linf"]], k = gp[["k"]], t0 = gp[["t0"]] - timing[t]),
      model = model, age = age, cv = cv, reflen = reflen, bin = bin, lmin = lmin,
      max = max(len_grid(gp[["linf"]], lmin, lmax, bin))
    )
    ialk_mat <- c(ialk)
    dim(ialk_mat) <- dim(ialk)[1:2]
    dimnames(ialk_mat) <- list(age = dimnames(ialk)$age, len = dimnames(ialk)$len)
    
    if (is.null(gear)) {
      cond <- ialk_mat                      # no selectivity: original behaviour
    } else {
      sel_len_vec <- as.numeric(gear$sel_len)
      names(sel_len_vec) <- dimnames(gear$sel_len)$len
      cond <- condition_alk(ialk_mat, sel_len_vec)
    }
    
    if (is.null(out)) {
      out <- FLQuant(
        0,
        dimnames = list(
          len = colnames(cond), year = yrs, unit = "unique",
          season = "all", area = "unique", iter = seq_len(its)
        )
      )
    }
    
    for (y in seq_along(yrs)) {
      for (i in seq_len(its)) {
        
        n_a <- as.numeric(N_a[, y, , , , i])
        if (sum(n_a) == 0) next
        
        p_len <- as.numeric((n_a / sum(n_a)) %*% cond)
        
        if (is.null(u)) {
          len_n <- apply(rmultinom(ess_t, 1, prob = p_len), 1, sum)
        } else {
          len_n <- rmultinom_crn(u[[t]][, y, i], p_len)
        }
        
        out[, y, , , , i] <- out[, y, , , , i] + len_n
      }
    }
  }
  
  units(out) <- unit
  
  if (scale) {
    fac <- apply(N_a, 2:6, sum) / apply(out, 2:6, sum)
    fac[!is.finite(fac)] <- 0
    out <- out %*% fac
  }
  
  out
}


#' Generate common random numbers for length sampling
#'
#' Generates one array of uniform deviates for each within-year timing event
#' for use by [len.sim()].
#'
#' @param gear Gear object returned by [build_gear()]. `gear$ess_len` defines
#'   the total length-sample size.
#' @param timing Numeric vector of within-year sampling times.
#' @param years Vector of years.
#' @param nsim Integer number of simulation iterations.
#' @param seed Optional integer random seed.
#'
#' @return A list with one `[draw, year, iter]` array per timing event.
#'
#' @examples
#' \dontrun{
#' u_len <- rUnif_len(trawl, timing = seq(0, 11/12, 1/12),
#'                     years = 2025:2054, nsim = 100, seed = 456)
#' length(u_len)          # one array per timing event
#' dim(u_len[[1]])        # ess_t x years x nsim
#' }
#'
#' @export
rUnif_len <- function(gear, timing, years, nsim, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  ess_t <- round(gear$ess_len / length(timing))
  yrs <- as.character(years)
  
  lapply(seq_along(timing), function(t) {
    array(
      stats::runif(ess_t * length(yrs) * nsim),
      dim = c(ess_t, length(yrs), nsim),
      dimnames = list(draw = seq_len(ess_t), year = yrs, iter = seq_len(nsim))
    )
  })
}


#' Simulate a two-stage gear-specific length-frequency sample
#'
#' Generates fishery-dependent length-frequency data from gear-specific
#' catch-at-age using a two-stage sampling design.
#'
#' @param object An `FLQuant` containing true catch-at-age or numbers-at-age
#'   for one gear.
#' @param gear Gear object returned by [build_gear()].
#' @param ess_age Integer effective sample size for the catch-at-age draw.
#' @param ess_len Integer effective sample size for the length draw.
#' @param u_age Optional common-random-number array with dimensions
#'   `[ess_age, year, iter]` for the age draw.
#' @param u_len Optional common-random-number array for one sampling time, or
#'   a list of arrays when `timing` contains multiple values.
#' @param scale Logical; if `TRUE`, rescale each sampled length distribution
#'   to the corresponding total in `object`.
#' @param timing Optional numeric vector of within-year sampling times. If
#'   `NULL`, use the conditional ALK stored in `gear`. Otherwise, rebuild the
#'   conditional ALK at each time using `t0 - timing`.
#' @param params Growth parameters coercible to an `FLPar`. Required when
#'   `timing` is supplied.
#' @param model Growth function or model object. Defaults to `vonbert`.
#' @param reflen,cv,lmin,lmax_mult,bin Arguments for rebuilding the
#'   conditional ALK when `timing` is supplied; `cv`, `lmin`, `lmax_mult` and
#'   `bin` default to `gear$grid`, so the rebuilt key always matches the gear.
#'
#' @return An `FLQuant` containing sampled length frequencies by year and
#'   iteration.
#'
#' @details
#' Sampling proceeds as follows:
#' \enumerate{
#'   \item A catch-at-age composition is sampled with effective sample size
#'     `ess_age` using [ca.sim()].
#'   \item The sampled age composition is projected through
#'     \eqn{P(l \mid a, g)} and sampled at effective sample size `ess_len`.
#' }
#'
#' @examples
#' \dontrun{
#' lfd_trawl <- lfd.sim(catch_n_gear$Trawl, trawl)
#' devs <- rUnif_lfd(om_gears, years = an(dimnames(catch_n_gear$Trawl)$year),
#'                    nsim = dim(catch_n_gear$Trawl)[6], seed = 456)
#' lfd_trawl_crn <- lfd.sim(
#'   catch_n_gear$Trawl,
#'   trawl,
#'   u_age = devs$age$Trawl,
#'   u_len = devs$len$Trawl
#' )
#' }
#'
#' @export
lfd.sim <- function(object, gear,
                    ess_age = gear$ess_age, ess_len = gear$ess_len,
                    u_age = NULL, u_len = NULL, scale = TRUE,
                    timing = NULL, params = gear$lhpar, model = vonbert,
                    reflen = NULL, cv = gear$grid$cv %||% 0.1,
                    lmin = gear$grid$lmin %||% 5,
                    lmax_mult = gear$grid$lmax_mult %||% 1.2,
                    bin = gear$grid$bin %||% 1) {
  
  age <- an(dimnames(object)$age)
  years <- dimnames(object)$year
  its <- dims(object)$iter
  
  rebuild_alk <- !is.null(timing)
  multi_timing <- rebuild_alk && length(timing) > 1L

  ## 'rebuild_alk': any timing supplied means gear$condALK (built at
  ## build_gear()'s own timing) is not reused; the key is rebuilt at this
  ## call's timing on the gear's grid. 'multi_timing': several timing
  ## events within one year, whose length draws are summed.
  if (rebuild_alk && is.null(params)) {
    stop("'params' (FLPar(linf, k, t0)) is required when 'timing' is supplied ",
         "-- gear$lhpar was also NULL. Pass params explicitly, or rebuild gear ",
         "via build_gear() so it carries its own lhpar.")
  }

  build_condALK_at <- function(t) {
    ## gear_ialk() builds the key on the gear's grid and stops if it does not
    ## match gear$condALK (e.g. an inconsistent lmax_mult override)
    g <- gear
    g$grid <- list(lmin = lmin, lmax_mult = lmax_mult, bin = bin, cv = cv,
                   timing = t)
    ialk_mat <- gear_ialk(g, age = age, timing = t, params = params,
                          model = model, reflen = reflen)
    condition_alk(ialk_mat, sel_len_vec)
  }
  
  ## per-timing-event conditioned ALKs, built once up front (not per
  ## year/iter -- growth timing doesn't depend on year or iteration)
  if (rebuild_alk) {
    sel_len_vec <- as.numeric(gear$sel_len)
    names(sel_len_vec) <- dimnames(gear$sel_len)$len
  }
  
  if (multi_timing) {
    
    cond_list <- lapply(timing, build_condALK_at)
    
    ## split_ess(): largest-remainder split so the per-timing ESS values
    ## sum exactly to ess_len, rather than rep(round(ess_len/S), S), which
    ## can silently over/under-count (e.g. 50/4 -> 4*12=48, not 50)
    split_ess <- function(n, S) {
      rep(n %/% S, S) + (seq_len(S) <= n %% S)
    }
    ess_len_t <- split_ess(ess_len, length(timing))
    len_lower <- colnames(cond_list[[1]])
    
  } else if (rebuild_alk) {
    ## single scalar timing: rebuild once, then behave exactly like the
    ## no-rebuild branch below (same ess_len, same condALK-shaped object)
    condALK <- build_condALK_at(timing)
    len_lower <- colnames(condALK)
    
  } else {
    condALK <- gear$condALK
    len_lower <- colnames(condALK)
  }
  
  if (!identical(dimnames(gear$condALK)$age, as.character(age))) {
    stop(
      "'gear$condALK' age dimension does not match 'object' age dimension. ",
      "Check 'object' and 'gear' were built for the same age range."
    )
  }
  
  out <- FLQuant(
    NA_real_,
    dimnames = list(
      len = len_lower, year = years, unit = "unique",
      season = "all", area = "unique", iter = seq_len(its)
    )
  )
  
  ## --- stage 1: low-ESS catch-at-age sample, single annual draw ---------
  age_n_flq <- ca.sim(object, ess = ess_age, u = u_age)
  
  ## --- stage 2: larger length sample, optionally summed across timing --
  for (y in seq_along(years)) {
    for (i in seq_len(its)) {
      
      age_n <- as.numeric(age_n_flq[, y, , , , i])
      
      if (sum(age_n) == 0) {
        age_n <- as.numeric(object[, y, , , , i])
      }
      
      age_p <- age_n / sum(age_n)
      
      if (multi_timing) {
        
        len_n <- rep(0, length(len_lower))
        for (t in seq_along(timing)) {
          p_len_t <- as.numeric(age_p %*% cond_list[[t]])
          
          if (is.null(u_len)) {
            draw_t <- apply(rmultinom(ess_len_t[t], 1, prob = p_len_t), 1, sum)
          } else {
            ## u_len must be a list of length(timing), each element a
            ## [ess_len_t, year, iter] array -- see rUnif_lfd()'s
            ## per-timing extension
            draw_t <- rmultinom_crn(u_len[[t]][, y, i], p_len_t)
          }
          len_n <- len_n + draw_t
        }
        
      } else {
        p_len <- as.numeric(age_p %*% condALK)
        
        if (is.null(u_len)) {
          len_n <- apply(rmultinom(ess_len, 1, prob = p_len), 1, sum)
        } else {
          len_n <- rmultinom_crn(u_len[, y, i], p_len)
        }
      }
      
      out[, y, , , , i] <- len_n
    }
  }
  
  units(out) <- "cm"
  
  if (scale) {
    fac <- apply(object, 2:6, sum) / apply(out, 2:6, sum)
    fac[!is.finite(fac)] <- 0
    out <- out %*% fac
  }
  
  out
}


#' Simulate seasonal length-frequency samples
#'
#' Applies [lfd.sim()] separately to each season of a season-structured
#' catch-at-age object. Seasonal midpoint timing is
#' `(season - 0.5) / n_seasons`.
#'
#' @param object An `FLQuant` containing catch-at-age with a populated season
#'   dimension.
#' @param gear Gear object returned by [build_gear()].
#' @param n_seasons Integer number of seasons. Defaults to the size of the
#'   season dimension in `object`.
#' @param ess_age,ess_len Effective sample sizes. A scalar is divided equally,
#'   after rounding, among seasons; a vector of length `n_seasons` specifies
#'   season-specific values.
#' @param params Growth parameters coercible to an `FLPar`; must contain
#'   `linf`, `k`, and `t0`.
#' @param ... Additional arguments passed to [lfd.sim()].
#'
#' @return An `FLQuant` containing sampled length frequencies by season, year,
#'   and iteration.
#'
#' @details
#' A separate stage-1 age sample is drawn for each season. Scalar sample sizes
#' are divided using `round(ess / n_seasons)` in every season.
#'
#' @examples
#' \dontrun{
#' lfd.sim.season(obj, trawl, n_seasons = 4, params = lhpars)
#' lfd.sim.season(obj, trawl, n_seasons = 4, ess_len = c(0, 200, 200, 100),
#'                 params = lhpars)
#' }
#'
#' @export
lfd.sim.season <- function(object, gear, n_seasons = dim(object)[4],
                           ess_age = gear$ess_age, ess_len = gear$ess_len,
                           params = NULL, ...) {
  
  if (is.null(params)) {
    stop("'params' (FLPar(linf, k, t0)) is required for seasonal timing adjustment.")
  }
  
  ## resolve ess_age/ess_len to length-n_seasons vectors
  resolve_ess <- function(ess, label) {
    if (length(ess) == 1) {
      rep(round(ess / n_seasons), n_seasons)
    } else if (length(ess) == n_seasons) {
      ess
    } else {
      stop("'", label, "' must be a single value or a vector of length n_seasons (",
           n_seasons, "), got length ", length(ess), ".")
    }
  }
  
  ess_age_s <- resolve_ess(ess_age, "ess_age")
  ess_len_s <- resolve_ess(ess_len, "ess_len")
  
  out_list <- lapply(seq_len(n_seasons), function(s) {
    obj_s <- object[, , , s, , ]
    
    lfd.sim(
      obj_s, gear,
      ess_age = ess_age_s[s], ess_len = ess_len_s[s],
      timing = (s - 0.5) / n_seasons,
      params = params,
      ...
    )
  })
  
  ## assemble into a genuine season-dimensioned FLQuant
  template <- out_list[[1]]
  out <- FLQuant(
    NA_real_,
    dimnames = list(
      len = dimnames(template)$len,
      year = dimnames(template)$year,
      unit = "unique",
      season = ac(seq_len(n_seasons)),
      area = "unique",
      iter = dimnames(template)$iter
    )
  )
  
  for (s in seq_len(n_seasons)) {
    out[, , , s, , ] <- out_list[[s]]
  }
  
  units(out) <- "cm"
  
  out
}


#' Generate common random numbers for two-stage length sampling
#'
#' Generates uniform-deviate arrays for the age and length stages of
#' [lfd.sim()] for every gear.
#'
#' @param om_gears Named list of gear objects returned by [build_gear()].
#' @param years Vector of years.
#' @param nsim Integer number of simulation iterations.
#' @param seed Optional integer random seed.
#'
#' @return A list with components `age` and `len`. Each component is a named
#'   list of `[draw, year, iter]` arrays, one per gear.
#'
#' @examples
#' \dontrun{
#' om_gears <- list(Trawl = trawl, Gillnet = gillnet)
#' devs <- rUnif_lfd(om_gears, years = 2025:2054, nsim = 100, seed = 456)
#' }
#'
#' @export
rUnif_lfd <- function(om_gears, years, nsim, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  gears <- names(om_gears)
  yrs <- as.character(years)
  
  mk <- function(field) {
    stats::setNames(lapply(gears, function(g) {
      n <- om_gears[[g]][[field]]
      array(
        stats::runif(n * length(yrs) * nsim),
        dim = c(n, length(yrs), nsim),
        dimnames = list(draw = seq_len(n), year = yrs, iter = seq_len(nsim))
      )
    }), gears)
  }
  
  list(age = mk("ess_age"), len = mk("ess_len"))
}


#' Disaggregate gear catch-at-age by season
#'
#' Disaggregates annual gear-specific catch-at-age into equal within-year
#' seasons under constant instantaneous fishing and natural mortality.
#'
#' @param stock An `FLStock` or `FLStockR` containing annual stock numbers,
#'   fishing mortality, and natural mortality.
#' @param f_age_gear Named list of gear-specific fishing-mortality-at-age
#'   `FLQuant` objects.
#' @param n_seasons Integer number of equal within-year seasons.
#'
#' @return A named list of gear-specific `FLQuant` objects with a populated
#'   season dimension.
#'
#' @details
#' Population numbers at the start of season `s` are
#' \eqn{N_s=N_0\exp[-Z(s-1)/S]}, where `S` is `n_seasons`. Gear catch in
#' that season is calculated with the corresponding partial fishing mortality.
#'
#' @examples
#' \dontrun{
#' seasonal_catch <- seasonalize_catch_n_gear(om_hist, f_age_gear, n_seasons = 4)
#' }
#'
#' @export
seasonalize_catch_n_gear <- function(stock, f_age_gear, n_seasons) {
  gears <- names(f_age_gear)
  stock_n <- stock.n(stock)
  z_total <- harvest(stock) + m(stock)   # OM's own aggregate Z
  
  ## f_age_gear may span the full projection (e.g. 1972:2024) while `stock`
  ## has already been subset to a shorter window (e.g. the last 5 data
  ## years) -- align years before any arithmetic against z_total/stock_n,
  ## or the division below is non-conformable
  yrs <- dimnames(stock)$year
  f_age_gear <- lapply(f_age_gear, function(x) x[, yrs])
  
  mk_season_template <- function(x) {
    d <- dimnames(x); d$season <- ac(seq_len(n_seasons))
    FLQuant(NA_real_, dimnames = d)
  }
  
  n_season <- mk_season_template(stock_n)
  for (s in seq_len(n_seasons)) {
    n_season[, , , s, , ] <- stock_n * exp(-z_total * (s - 1) / n_seasons)
  }
  
  catch_season <- setNames(lapply(gears, function(g) {
    out <- mk_season_template(stock_n)
    for (s in seq_len(n_seasons)) {
      out[, , , s, , ] <- n_season[, , , s, , ] *
        (1 - exp(-z_total / n_seasons)) * (f_age_gear[[g]] / z_total)
    }
    out
  }), gears)
  
  catch_season
}


#' Calculate population numbers-at-age by season
#'
#' Propagates annual beginning-of-year numbers-at-age through equal
#' within-year seasons under constant total instantaneous mortality.
#'
#' @param stock An `FLStock` or `FLStockR` containing annual stock numbers,
#'   fishing mortality, and natural mortality.
#' @param n_seasons Integer number of equal within-year seasons.
#'
#' @return An `FLQuant` containing numbers-at-age at the start of each season.
#'
#' @examples
#' \dontrun{
#' seasonal_numbers <- pop_n_season(om_hist, n_seasons = 4)
#' }
#'
#' @export
pop_n_season <- function(stock, n_seasons) {
  stock_n <- stock.n(stock)
  z_total <- harvest(stock) + m(stock)

  d <- dimnames(stock_n); d$season <- ac(seq_len(n_seasons))
  n_season <- FLQuant(NA_real_, dimnames = d)

  for (s in seq_len(n_seasons)) {
    n_season[, , , s, , ] <- stock_n * exp(-z_total * (s - 1) / n_seasons)
  }

  n_season
}


#' Plot seasonal length-frequency distributions
#'
#' Plots seasonal length-frequency distributions for a selected year, with
#' seasons distinguished by colour and gears shown in separate panels.
#'
#' @param lfd_season Named list of seasonal length-frequency `FLQuant` objects.
#' @param year Character or numeric year to plot.
#' @param iter Integer simulation iteration. Defaults to `1`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' plot_lfd_season(list(Trawl = lfd_season_trawl, Gillnet = lfd_season_gillnet),
#'                  year = "2020")
#' }
#'
#' @export
plot_lfd_season <- function(lfd_season, year, iter = 1) {
  
  df <- do.call(rbind, lapply(names(lfd_season), function(g) {
    x <- lfd_season[[g]][, ac(year), , , , iter]
    d <- as.data.frame(x, cohort = FALSE)
    d$gear <- g
    d
  }))
  
  df$len <- as.numeric(as.character(df$len))
  df$season <- factor(df$season)
  
  ggplot2::ggplot(df, ggplot2::aes(x = len, y = data, colour = season, group = season)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~gear, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      x = "Length (cm)", y = "Sampled count", colour = "Season",
      title = paste("Seasonal length-frequency samples,", year)
    ) +
    ggplot2::theme_bw()
}


#' Plot length-frequency distributions as temporal ridgelines
#'
#' Displays annual or seasonal length-frequency distributions along a
#' continuous time axis, with one panel per gear. Optional von Bertalanffy
#' cohort trajectories can be overlaid for reference.
#'
#' @param panels Named list of length-frequency `FLQuant` objects. List names
#'   determine panel order and labels.
#' @param lhpar Optional growth parameters containing `linf`, `k`, and `t0`
#'   for cohort reference lines.
#' @param scale Numeric ridge-height multiplier. Defaults to `0.9`.
#' @param cohort_pad Integer number of additional cohorts generated beyond
#'   each end of the plotted time range. Defaults to `15`.
#' @param annual Logical; if `TRUE`, sum the season dimension before plotting.
#'
#' @return A `ggplot` object.
#'
#' @details
#' For seasonal data, ridge positions are calculated as
#' `year + (season - 0.5) / n_seasons`. When `annual = TRUE`, seasonal data
#' are aggregated with `seasonSums()`.
#'
#' @examples
#' \dontrun{
#' plot_lfd_ridgeline(lfds, lhpar = lhpars)
#' plot_lfd_ridgeline(lfd_season, annual = TRUE, lhpar = lhpars)
#' }
#'
#' @export
plot_lfd_ridgeline <- function(panels, lhpar = NULL, scale = 0.9, cohort_pad = 15,
                               annual = FALSE) {
  
  if (annual) panels <- lapply(panels, seasonSums)
  
  panel_levels <- names(panels)
  
  df <- do.call(rbind, lapply(panel_levels, function(p) {
    d <- as.data.frame(panels[[p]], cohort = FALSE)
    d$len <- as.numeric(as.character(d$len))
    d$year <- as.numeric(as.character(d$year))
    d$panel <- p
    d
  }))
  df$panel <- factor(df$panel, levels = panel_levels)
  
  has_season <- !annual && length(unique(df$season)) > 1
  
  if (has_season) {
    n_seasons <- length(unique(df$season))
    df$season <- factor(as.numeric(as.character(df$season)))
    df$time <- df$year + (as.numeric(as.character(df$season)) - 0.5) / n_seasons
  } else {
    df$time <- df$year
  }
  
  group_cols <- c("panel", "time")
  df <- df[do.call(order, df[group_cols]), ]
  
  yr_step <- min(diff(sort(unique(df$year))))
  slot_step <- if (has_season) yr_step / n_seasons else yr_step
  
  df <- do.call(rbind, lapply(split(df, df[group_cols], drop = TRUE), function(g) {
    if (!nrow(g) || all(g$data == 0)) { g$dens <- 0; return(g) }
    g$dens <- g$data / max(g$data) * scale * slot_step
    g
  }))
  
  df$grp <- interaction(df$panel, df$time)
  
  p <- ggplot(df)
  
  if (has_season) {
    p <- p + geom_ribbon(
      aes(x = len, ymin = time, ymax = time + dens, group = grp, fill = season),
      alpha = 0.8, colour = NA
    ) + scale_fill_viridis_d(name = "Month", option = "D")
  } else {
    p <- p + geom_ribbon(
      aes(x = len, ymin = time, ymax = time + dens, group = grp),
      fill = "#5B9BD5", alpha = 0.75, colour = NA
    )
  }
  
  p <- p +
    coord_flip() +
    facet_wrap(~panel, ncol = 1) +
    labs(x = "Length (cm)", y = "Year") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  if (!is.null(lhpar)) {
    yr_range <- range(df$time)
    len_range <- range(df$len)
    linf <- c(lhpar["linf"]); k <- c(lhpar["k"]); t0 <- c(lhpar["t0"])
    cohorts <- seq(floor(yr_range[1]) - cohort_pad, ceiling(yr_range[2]) + cohort_pad)
    age_grid <- seq(0, -log(1 - 0.98) / k, by = 0.1)
    
    ref <- do.call(rbind, lapply(cohorts, function(cy) {
      L <- linf * (1 - exp(-k * (age_grid - t0)))
      data.frame(t = age_grid + cy, L = L, cohort = cy)
    }))
    ref <- ref[ref$t >= yr_range[1] & ref$t <= yr_range[2] &
                 ref$L >= len_range[1] & ref$L <= len_range[2], ]
    
    p <- p + geom_line(
      data = ref, aes(x = L, y = t, group = cohort),
      colour = "grey50", linetype = 2, linewidth = 0.3, inherit.aes = FALSE
    )
  }
  
  p
}


#' Partial (Baranov) catch-at-age by gear
#'
#' Splits the catch implied by the stock's realised total F among gears:
#' \deqn{C_{g,a,y}=N_{a,y}\frac{\phi_{g,a,y}F_{a,y}}{Z_{a,y}}
#'       \left(1-e^{-Z_{a,y}}\right),\qquad Z_{a,y}=F_{a,y}+M_{a,y}.}
#' Because shares sum to one, \eqn{\sum_g C_{g,a,y}} equals the Baranov total
#' catch exactly. This is the canonical definition used by the OM and by the
#' OEM (`sampling.lfd.oem()`), which should call this function rather than
#' reproduce it.
#'
#' @param stock `FLStock` with realised `stock.n`, `harvest` and `m`.
#' @param shares Either a list returned by [build_gear()] (static
#'   `f_share`, recycled over years and iterations) or `FLQuants` of shares
#'   from [f_share_gear()].
#' @param years Years to compute. Default all years of `stock`.
#' @param check Logical; verify that gear catches sum to the Baranov total
#'   (relative tolerance `1e-8`). Default `TRUE`.
#' @return `FLQuants` of catch numbers-at-age by gear. The Baranov total is
#'   attached as `attr(, "total")`.
#' @examples
#' \dontrun{
#' cn <- catch_n_gear(om_hist, om_gears, years = 1972:2024)  # static shares
#' cn <- catch_n_gear(om_hist, f_share_gear(f_age))          # time-varying
#' }
#' @export
catch_n_gear <- function(stock, shares, years = dimnames(stock)$year,
                         check = TRUE) {
  years <- as.character(years)
  stk <- stock[, years]
  Ft  <- harvest(stk)
  Z   <- Ft + m(stk)
  tot <- stock.n(stk) * Ft / Z * (1 - exp(-Z))
  tot[Z == 0] <- 0

  out <- if (inherits(shares, "FLQuants")) {
    FLQuants(lapply(shares, function(s) tot %*% s[, years]))
  } else {
    ## static age-only share: recycle over year/unit/season/area/iter of the
    ## stock (age varies fastest), so unit/area dimnames always match `tot`
    FLQuants(lapply(shares, function(g) {
      if (is.null(g$f_share))
        stop("Rebuild om_gears with build_gear(..., fbar_range = ).")
      if (length(c(g$f_share)) != dim(tot)[1])
        stop("f_share ages do not match the stock's age dimension.")
      s_full <- tot
      s_full[] <- c(g$f_share)
      tot * s_full
    }))
  }

  if (check) {
    err <- max(abs(Reduce(`+`, out) - tot), na.rm = TRUE) /
      max(abs(tot), na.rm = TRUE)
    if (err > 1e-8) stop("Gear catches do not sum to total (rel. error ",
                         signif(err, 3), "). Do the shares sum to one?")
  }
  attr(out, "total") <- tot
  out
}


#' Catch by gear and year
#'
#' @param catch_n `FLQuants` of catch numbers-at-age by gear
#'   (from [catch_n_gear()]).
#' @param wt `FLStock` (its `catch.wt` is used) or an `FLQuant` of
#'   weight-at-age. Ignored for `type = "numbers"`.
#' @param type `"weight"` (default) or `"numbers"`.
#' @return `FLQuants` of annual catch by gear (age summed).
#' @export
catch_gear <- function(catch_n, wt = NULL, type = c("weight", "numbers")) {
  type <- match.arg(type)
  if (type == "weight") {
    if (is.null(wt)) stop("'wt' required for type = 'weight'.")
    if (is(wt, "FLStock")) wt <- catch.wt(wt)
  }
  FLQuants(lapply(catch_n, function(cn) {
    if (type == "numbers") return(quantSums(cn))
    quantSums(cn * wt[, dimnames(cn)$year])
  }))
}


#' Relative catch by gear (fiticc() `catch_by_gear` input)
#'
#' Sums gear catches over a window of years and normalises across gears:
#' \deqn{\pi_g=\frac{\sum_{y\in\mathcal Y}C_{g,y}}
#'                  {\sum_h\sum_{y\in\mathcal Y}C_{h,y}}.}
#' Supplying the realised catch split over the LFD window (rather than the
#' nominal `f_mult`) mirrors what an assessment would actually be given, and
#' reflects that a gear's catch share differs from its Fbar share whenever
#' gears select different ages.
#'
#' @inheritParams catch_gear
#' @param years Years to sum over (e.g. `lfd_years`). Default all.
#' @param iter Optional iteration(s). With a single `iter` a named vector is
#'   returned; otherwise an `iter x gear` matrix (one row per fit).
#' @return Named numeric vector or matrix of catch shares (rows sum to one).
#' @examples
#' \dontrun{
#' catch_by_gear(catch_n, om_hist, years = lfd_years, iter = 1)
#' cbg <- catch_by_gear(catch_n, om_hist, years = lfd_years)  # all iters
#' fiticc(lfd_i, stklen_i, sel_fun = sel_fun, catch_by_gear = cbg[i, ])
#' }
#' @export
catch_by_gear <- function(catch_n, wt = NULL, years = NULL, iter = NULL,
                          type = c("weight", "numbers")) {
  cg <- catch_gear(catch_n, wt, type = type)
  if (!is.null(years)) cg <- FLQuants(lapply(cg, function(x) x[, as.character(years)]))
  its <- dim(cg[[1]])[6]
  m <- vapply(cg, function(x) apply(x@.Data, 6, sum, na.rm = TRUE),
              numeric(its))
  m <- matrix(m, nrow = its, dimnames = list(iter = seq_len(its),
                                             gear = names(cg)))
  m <- m / rowSums(m)
  if (!is.null(iter)) {
    m <- m[iter, , drop = FALSE]
    if (nrow(m) == 1) m <- setNames(m[1, ], colnames(m))
  }
  m
}


#' Raise a length-frequency sample to catch weight (sum of products)
#'
#' Length analogue of `FLRef::sops()` (which raises catch-at-age with
#' `catch.wt`): scales each year/unit/iteration of a length sample so that
#' its sum of products with weight-at-length equals the catch weight,
#' \deqn{\tilde N_{l,y,u}=N_{l,y,u}\,
#'       \frac{C_{y,u}}{\sum_{l'}N_{l',y,u}\,W(l')}.}
#' Years/units with zero sample or zero catch return zero.
#'
#' @param lfd `FLQuant` length sample (len x year x unit x ... x iter).
#' @param wt_len Named numeric weight-at-length (names = length bins of
#'   `lfd`), e.g. `gear$wt_len` from [build_gear()], or a gear object.
#' @param catch `FLQuant` catch weight (year x unit x ... x iter) in the same
#'   units as `wt_len`.
#' @return `FLQuant` of raised numbers-at-length, with `attr(, "raising")`
#'   the raising factor per year/unit/iter.
#' @seealso [sampling.lfd.oem()] (`lfd.raising = TRUE`)
#' @export
lfd_sop <- function(lfd, wt_len, catch) {
  if (is.list(wt_len)) wt_len <- wt_len$wt_len
  if (is.null(wt_len))
    stop("No weight-at-length: build the gear with lhpar containing 'a' and 'b'.")
  lens <- dimnames(lfd)$len
  if (!all(lens %in% names(wt_len)))
    stop("Length bins of 'lfd' are not all in names(wt_len).")
  w <- lfd
  w[] <- wt_len[lens]                       # recycled over year/unit/iter
  f <- catch / quantSums(lfd * w)
  f[!is.finite(f)] <- 0
  out <- lfd %*% f
  attr(out, "raising") <- f
  out
}


#' Length-frequency sampling observation error model
#'
#' An `FLoem` method for `mse::mp()` that generates fishery-dependent
#' length-frequency (LFD) samples and a reported catch split by gear (and,
#' optionally, by area), as a structural analogue of `mse::sampling.oem()`.
#' The abundance-index block of `sampling.oem()` is replaced by an LFD block:
#' true catch-at-age is partitioned by gear with [catch_n_gear()] and passed
#' through the two-stage sampler [lfd.sim()].
#'
#' A single-area stock (one `unit`, typically `"unique"`) is the default
#' case. Multiple areas are represented by the `unit` dimension of the
#' operating-model stock, each area with its own gear set.
#'
#' @param stk `FLStock` passed by `mp()` (the OM stock).
#' @param deviances List with optional elements
#'   \describe{
#'     \item{`stk`}{`FLQuants` of multiplicative deviances on `FLStock`
#'       slots (e.g. `catch.n`), applied exactly as in `sampling.oem()`.}
#'     \item{`lfd`}{Uniform common-random-number streams for the two
#'       sampling stages. Single area: the list returned by [rUnif_lfd()],
#'       i.e. `list(age = list(<gear> = array), len = list(<gear> = array))`,
#'       each array `[draw, year, iter]`. Multiple areas: a list named by
#'       area (unit) of such lists. If `NULL`, draws use the live RNG.}
#'     \item{`catch_gear`}{`FLQuants` named by gear of multiplicative
#'       reporting deviances on catch weight (year x unit x iter). If `NULL`,
#'       the catch split is reported without error.}
#'   }
#' @param observations List with elements
#'   \describe{
#'     \item{`stk`}{Observed `FLStock`, as in `sampling.oem()`.}
#'     \item{`lfd`}{`FLQuants` (or list) named by gear: len x year x unit x
#'       iter, spanning at least `y0:fy`, seeded with historical samples up
#'       to the first data year. Areas where a gear is absent stay 0.}
#'     \item{`catch_gear`}{`FLQuants` (or list) named by gear: catch weight
#'       by year x unit x iter, seeded like `lfd`.}
#'   }
#' @param args `mp()` arguments; uses `y0`, `dy`, `dys`. `y0` must be set in
#'   `mp()`'s own `args` (the first year of the LFD series), not in
#'   `FLoem(args = )`.
#' @param tracking `mp()` tracking object, returned unchanged.
#' @param om_gears Gear objects from [build_gear()] (built with
#'   `fbar_range`, so each gear carries `f_share`). Single area: the named
#'   list of gears. Multiple areas: a list named by area (matching
#'   `dimnames(stk)$unit`) of named gear lists; a gear may be absent from an
#'   area. A plain gear list with a multi-unit stock is applied to every area.
#' @param f_share Optional `FLQuants` named by gear of time-varying shares of
#'   total F-at-age (age x year x unit x iter, covering the projection
#'   years), e.g. from [f_share_gear()]. Default `NULL` uses each gear's
#'   static `f_share`, which is exact whenever gear Fbar trajectories are
#'   proportional and needs no extension into the projection period.
#' @param lfd.raising Logical. If `TRUE` (intended for more than one unit),
#'   the per-unit samples passed to `est` are raised to each unit's reported
#'   gear catch by sum of products with the gear's weight-at-length
#'   (`gear$wt_len`, from `a` and `b` in `lhpar`; see [lfd_sop()]) and
#'   summed over units, so the pooled length composition is weighted by
#'   catch per unit rather than by sample size per unit. The reported catch
#'   by gear is summed over units as well. `observations$lfd` always keeps
#'   the raw samples by unit. Default `FALSE` passes the samples by unit
#'   unchanged.
#' @param lfd.scale With `lfd.raising = TRUE`: `"sample"` (default) rescales
#'   the pooled, catch-weighted LFD to the total number of fish measured
#'   across units, so the sample size seen by the likelihood is the real one;
#'   `"catch"` keeps the raised numbers caught (requires `wt_len` and catch
#'   in the same units).
#'
#' @return A list with `stk` (the observed stock, window `y0:dy`, with
#'   attributes `"lfd"` and `"catch_gear"`: `FLQuants` by gear holding the
#'   accumulated series `y0:dy`, by unit or, with `lfd.raising = TRUE`,
#'   pooled over units), `idx` (empty `FLIndices`), and the updated
#'   `observations` and `tracking`. The attributes are how the samples reach
#'   the `est` module (e.g. `flicc.sa()`), because `mp()` passes only `stk`
#'   and `idx` to it.
#'
#' @details
#' Each call is incremental, as in `sampling.oem()`: only years `dys` are
#' simulated and written into `observations`, so each CRN deviate is used
#' once and the same `FLoem` gives identical sampling noise to every MP.
#'
#' For area \eqn{u}, gear \eqn{g} and year \eqn{y \in} `dys`:
#' \enumerate{
#'   \item gear catch-at-age by the partial Baranov equation,
#'     \eqn{C_{g,a,y,u}=N_{a,y,u}\,\phi_{g,a,y,u}F_{a,y,u}
#'     (1-e^{-Z_{a,y,u}})/Z_{a,y,u}}, using the OM's realised F and the
#'     (possibly noised) stock;
#'   \item two-stage LFD sample with [lfd.sim()]: a multinomial age sample of
#'     size `ess_age`, projected through the gear's \eqn{P(l \mid a, g)} and
#'     sampled with size `ess_len`;
#'   \item reported catch weight
#'     \eqn{\sum_a C_{g,a,y,u}w_{a,y,u}\,\varepsilon_{g,y,u}}.
#' }
#' \strong{Raising across units} (`lfd.raising = TRUE`). For gear \eqn{g},
#' year \eqn{y} and unit \eqn{u} with sampled counts \eqn{n_{l,y,u}} and
#' reported catch weight \eqn{\hat C_{g,y,u}}:
#' \deqn{N_{l,y}=\sum_u n_{l,y,u}\,
#'   \frac{\hat C_{g,y,u}}{\sum_{l'} n_{l',y,u}W(l')},\qquad
#'   N^{(sample)}_{l,y}=N_{l,y}\,
#'   \frac{\sum_u\sum_l n_{l,y,u}}{\sum_l N_{l,y}}.}
#' Without raising, a unit that is heavily sampled relative to its catch
#' would dominate the pooled composition. With one unit, raising followed by
#' `lfd.scale = "sample"` returns the sample unchanged.
#'
#' Zero catch-at-age entries are set to `sqrt(.Machine$double.eps)` to avoid
#' 0/0 in [lfd.sim()]; a gear with no catch at all in an area-year leaves a
#' zero LFD.
#'
#' @seealso [build_gear()], [catch_n_gear()], [lfd.sim()], [rUnif_lfd()],
#'   `mse::sampling.oem()`
#'
#' @examples
#' \dontrun{
#' lfd.devs <- rUnif_lfd(om_gears, years = 2010:2034, nsim = 20, seed = 456)
#' oem <- FLoem(
#'   method = sampling.lfd.oem,
#'   observations = list(stk = stock(om), lfd = lfd_obs, catch_gear = cg_obs),
#'   deviances = list(stk = stk.devs, lfd = lfd.devs, catch_gear = cg.devs),
#'   args = list(om_gears = om_gears)
#' )
#' run <- mp(om, oem = oem, ctrl = ctrl,
#'           args = list(iy = 2024, fy = 2034, y0 = 2012, it = 20))
#'
#' ## two areas (units): gear sets by area, independent CRN streams by area,
#' ## LFDs raised to catch by area and pooled for the assessment
#' oem2 <- FLoem(
#'   method = sampling.lfd.oem,
#'   observations = list(stk = stock(om2), lfd = lfd_obs2, catch_gear = cg_obs2),
#'   deviances = list(stk = stk.devs2,
#'                    lfd = list(North = rUnif_lfd(gears_N, 2010:2034, 20, 1),
#'                               South = rUnif_lfd(gears_S, 2010:2034, 20, 2)),
#'                    catch_gear = cg.devs2),
#'   args = list(om_gears = list(North = gears_N, South = gears_S),
#'               lfd.raising = TRUE)
#' )
#' }
#' @export
sampling.lfd.oem <- function(stk, deviances, observations, args, tracking,
                             om_gears, f_share = NULL,
                             lfd.raising = FALSE,
                             lfd.scale = c("sample", "catch")) {

  spread(args)                                 # y0, dy, dys, ay, it, ...
  dys <- ac(dys)
  lfd.scale <- match.arg(lfd.scale)

  stk <- window(stk, start = y0, end = dy, extend = FALSE)
  obs <- window(observations$stk, start = y0, end = dy, extend = FALSE)
  areas <- dimnames(stk)$unit

  ## --- normalise om_gears to area -> gear --------------------------------
  is_gear <- function(x) is.list(x) && !is.null(x$condALK)
  if (all(vapply(om_gears, is_gear, logical(1))))
    om_gears <- setNames(rep(list(om_gears), length(areas)), areas)
  if (!all(areas %in% names(om_gears)))
    stop("sampling.lfd.oem(): names(om_gears) must match dimnames(stk)$unit: ",
         paste(areas, collapse = ", "))
  om_gears <- om_gears[areas]
  gears <- unique(unlist(lapply(om_gears, names)))

  ## --- normalise CRN streams to area -> list(age, len) -> gear ----------
  u_lfd <- deviances$lfd
  if (!is.null(u_lfd) && all(c("age", "len") %in% names(u_lfd))) {
    if (length(areas) > 1)
      stop("sampling.lfd.oem(): with several areas, deviances$lfd must be a ",
           "list by area of rUnif_lfd() outputs (independent streams per area).")
    u_lfd <- setNames(list(u_lfd), areas)
  }

  ## --- STK: observation error on stock slots (as sampling.oem()) --------
  if (!is.null(deviances$stk)) {
    for (i in names(deviances$stk)) {
      slot(stk, i)[, dys] <-
        do.call(i, list(object = stk))[, dys] %*% deviances$stk[[i]][, dys] + 1e-8
    }
    landings(stk)[, dys] <- computeLandings(stk[, dys])
    discards(stk)[, dys] <- computeDiscards(stk[, dys])
    catch(stk)[, dys]    <- computeCatch(stk[, dys])
  }

  ## --- gear catch-at-age by area, years dys only ------------------------
  cn <- lapply(setNames(nm = areas), function(a) {
    ga <- om_gears[[a]]
    sh <- if (is.null(f_share)) ga else
      FLQuants(lapply(f_share[names(ga)], function(x) x[, dys, a]))
    catch_n_gear(stk[, dys, a], sh, check = is.null(f_share))
  })

  ## --- LFD: two-stage sample per gear x area ----------------------------
  for (g in gears) {
    res <- observations$lfd[[g]][, dys]
    res[] <- 0
    for (a in areas) {
      gear_ga <- om_gears[[a]][[g]]
      if (is.null(gear_ga)) next                       # gear absent in area
      x <- cn[[a]][[g]]
      if (sum(x, na.rm = TRUE) == 0) next              # no catch: zero LFD
      x[x <= 0] <- sqrt(.Machine$double.eps)           # guard 0/0
      u_age <- u_lfd[[a]]$age[[g]]
      u_len <- u_lfd[[a]]$len[[g]]
      res[, , a] <- lfd.sim(
        x, gear = gear_ga,
        u_age = if (!is.null(u_age)) u_age[, dys, , drop = FALSE],
        u_len = if (!is.null(u_len)) u_len[, dys, , drop = FALSE],
        scale = FALSE
      )
    }
    observations$lfd[[g]][, dys] <- res
  }

  ## --- CATCH_GEAR: reported catch weight by gear x area -----------------
  for (g in gears) {
    tru <- observations$catch_gear[[g]][, dys]
    tru[] <- 0
    for (a in areas) {
      if (is.null(om_gears[[a]][[g]])) next
      tru[, , a] <- quantSums(cn[[a]][[g]] * catch.wt(stk)[, dys, a])
    }
    dev <- deviances$catch_gear[[g]]
    observations$catch_gear[[g]][, dys] <-
      if (is.null(dev)) tru else tru %*% dev[, dys]
  }

  ## --- store the observed stock (as sampling.oem()) ---------------------
  slots <- c("landings", "discards", "catch", "landings.n", "discards.n",
             "catch.n", "stock.n", "harvest",
             "landings.wt", "discards.wt", "catch.wt", "stock.wt")
  for (i in slots) slot(obs, i)[, dys] <- slot(stk, i)[, dys]
  units(obs) <- units(stk)
  observations$stk[, dys] <- obs[, dys]

  ## --- channels to est ---------------------------------------------------
  lfd_out <- lapply(setNames(nm = gears), function(g)
    window(observations$lfd[[g]], start = y0, end = dy))
  cg_out  <- lapply(setNames(nm = gears), function(g)
    window(observations$catch_gear[[g]], start = y0, end = dy))

  if (lfd.raising) {
    ## raise each unit's sample to that unit's reported gear catch (SOP),
    ## pool over units, then (default) rescale to the total number sampled
    lfd_out <- lapply(setNames(nm = gears), function(g) {
      x <- lfd_out[[g]]
      raised <- x
      raised[] <- 0
      for (a in areas) {
        gear_ga <- om_gears[[a]][[g]]
        if (is.null(gear_ga)) next
        raised[, , a] <- lfd_sop(x[, , a], gear_ga$wt_len, cg_out[[g]][, , a])
      }
      pooled <- unitSums(raised)
      if (lfd.scale == "sample") {
        f <- unitSums(quantSums(x)) / quantSums(pooled)
        f[!is.finite(f)] <- 0
        pooled <- pooled %*% f
      }
      pooled
    })
    cg_out <- lapply(cg_out, unitSums)
  }

  attr(obs, "lfd")        <- FLQuants(lfd_out)
  attr(obs, "catch_gear") <- FLQuants(cg_out)

  list(stk = obs, idx = FLIndices(), observations = observations,
       tracking = tracking)
}
