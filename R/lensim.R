# iALK() {{{
#
#' inverse ALK function with lmin added to FLCore::invALK 
#' @param params growth parameter, default FLPar(linf,k,t0)
#' @param model growth model, only option currently vonbert
#' @param age age vector
#' @param cv of length-at-age
#' @param lmax maximum upper length specified lmax*linf
#' @param max maximum size value
#' @param lmin milonimum length
#' @param reflen evokes fixed sd for L_a at sd = cv*reflen
#' @param bin length bin size, dafault 1
#' @param timing t0 assumed 1st January, default seq(0,11/12,1/12), but can be single event 0.5
#' @param unit default is "cm"
#' @return FLPar age-length matrix
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
# }}}

# ALK() {{{
#
#' ALK function
#' @param N_a numbers at age sample for single event
#' @param iALK from iALK() outout
#' @return FLPar of ALK
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

# alk.sample() {{{
#
#' generates annual ALK sample with length stratified sampling
#' @param lfds length frequency *FLQuant*
#' @param alks annual ALK proportions at age output form ALKs() *FLPars*
#' @param nbin number of samples per length bin
#' @param n.sample sample size of lfd 
#' @return FLPars of sampled ALK
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

# ALKs() {{{
#
#' annual ALK function
#' @param object FLQuant with numbers at age
#' @param iALK from iALK() outout
#' @return FLPars of ALK
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

# applyALK() {{{
#
#' applyALK function to length to age
#' @param lfd *FLQuant* with numbers at length
#' @param alks *FLPars* annual ALKs
#' @return FLQuant for numbers at age
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



#{{{
# condition
#' Fold length selectivity into an inverse ALK
#'
#' Shared helper: re-weights an inverse age-length key (age x len, rows
#' sum to 1) by a length-selectivity vector on the same length grid, and
#' renormalises each age row back to sum to 1. Ages with zero selected
#' probability get an all-zero row rather than a divide-by-zero. Used by
#' both `build_gear()` and the gear-aware `len.sim()` so the two stay
#' consistent rather than each carrying its own copy of this logic.
#'
#' @param ialk_mat Numeric matrix, age (rows) x len (cols), rows sum to 1.
#' @param sel_len_vec Numeric vector of length-selectivity, named by
#'   length. Reordered internally to match `colnames(ialk_mat)`.
#'
#' @return Numeric matrix, same dimensions as `ialk_mat`, rows sum to 1
#'   (or all-zero where no length was selected for that age).
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
#}}}

#' Generate survey (pulse) or continuous length-frequency samples
#'
#' Samples a length-frequency distribution from `N_a` by projecting it
#' through the biological inverse age-length key at one or more
#' within-year timing events (`timing`), optionally re-weighted by a
#' gear's length selectivity.
#'
#' Unlike [lfd.sim()]'s two-stage fishery-dependent design (a noisy
#' catch-at-age sample feeding a length draw), `len.sim()` samples length
#' directly from `N_a` at each timing event and sums across events -- the
#' single-stage design appropriate for a survey or continuous sampling
#' programme, where `N_a` is already the relevant numbers/catch-at-age
#' for that observation process rather than something to be re-sampled
#' for age first.
#'
#' Because growth (`t0 + timing[t]`) shifts within a year, the ALK -- and
#' therefore any gear-selectivity weighting applied to it via
#' [condition_alk()] -- is rebuilt fresh at each timing event, rather
#' than reusing a single pre-built `gear$condALK` as `lfd.sim()` does.
#'
#' @param N_a Numbers (or catch)-at-age `FLQuant`.
#' @param params Growth parameters, `FLPar(linf, k, t0)`.
#' @param model Growth model. Default `vonbert`.
#' @param ess Effective sample size, split evenly across `timing` events
#'   (`round(ess / length(timing))` per event, as in the original).
#'   Default `250`. If using `gear`, consider passing `gear$ess_len`
#'   explicitly.
#' @param timing Within-year timing of sampling events (as fractions of a
#'   year added to `t0`). Default `seq(0, 11/12, 1/12)` (roughly monthly,
#'   continuous sampling); pass a single value (e.g. `0.5`) for a single
#'   survey pulse.
#' @param unit Length unit label. Default `"cm"`.
#' @param scale Logical. If `TRUE` (default), rescale the output so its
#'   per-year/iter total matches `N_a`'s true total.
#' @param reflen,bin,cv,lmin,lmax Passed to `iALK()`, as in the original.
#' @param gear Optional gear object as returned by [build_gear()]. If
#'   supplied, its `sel_len` is folded into the ALK at every timing event
#'   via [condition_alk()]. If `NULL` (default), no selectivity weighting
#'   is applied -- reproduces the original `len.sim()` behaviour exactly.
#' @param u Optional list of length `length(timing)`, each element a
#'   `[ess_t, year, iter]` array of pre-generated `Uniform(0,1)` deviates
#'   (`ess_t = round(ess / length(timing))`) for CRN-safe sampling via
#'   [rmultinom_crn()] -- see [rUnif_len()]. If `NULL` (default), falls
#'   back to `stats::rmultinom()` against the live RNG, i.e. the
#'   original behaviour.
#'
#' @return `FLQuant` with a `len` dimension, the sampled length-frequency
#'   by year and iter.
#'
#' @examples
#' \dontrun{
#' ## original behaviour, unchanged: no gear, no CRN
#' lfd_survey <- len.sim(stock.n(om)[, "2024"], params = lhpar, ess = 300,
#'                        timing = 0.5)
#'
#' ## gear-aware, CRN-safe
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
                    reflen = NULL, bin = 1, cv = 0.1, lmin = 5, lmax = 1.2,
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
      model = model, age = age, lmax = lmax, reflen = reflen, bin = bin, lmin = lmin
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

#' Pre-generate CRN uniform-deviate streams for len.sim()
#'
#' Generates one `Uniform(0,1)` array per within-year timing event, each
#' shaped `[ess_t, year, iter]`, for use with [len.sim()]'s `u` argument
#' via [rmultinom_crn()]. Analogous to [rUnif_lfd()], but shaped for
#' `len.sim()`'s timing-event loop rather than the two-stage age/length
#' design.
#'
#' @param gear A gear object as returned by [build_gear()]; used for
#'   `ess_len` (per [len.sim()]'s `ess` argument, split across timing
#'   events the same way `len.sim()` does internally).
#' @param timing Within-year timing vector, as passed to [len.sim()].
#' @param years Integer vector of years to generate deviates for.
#' @param nsim Integer. Number of OM iterations.
#' @param seed Optional integer seed. If `NULL` (default), the current
#'   RNG state is used as-is.
#'
#' @return A list of length `length(timing)`, each element an array of
#'   dimension `[round(gear$ess_len / length(timing)), length(years), nsim]`.
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



# lfd_sim.R
# Draft, not yet run against a live R session.
#
# Completes FLRef::lfd.sim() (R/OMsim.R), which as currently committed
# takes a `sel` argument but never uses it, and references an undefined
# `N_a` instead of its own `object` argument. This version:
#   1. actually folds gear selectivity into the ALK (via a pre-built
#      build_gear() condALK, so the reweighting only happens once, at
#      build_gear() time, not on every sampling call), and
#   2. implements the two-stage design explicitly: a low-ESS catch-at-age
#      sample (ca.sim()) feeds the length draw, rather than sampling
#      length straight from the true age composition.
#
# Depends on ca.sim() (ca_sim_crn.R) and rmultinom_crn() (rmultinom_crn.R)
# already being loaded.

#' Two-stage gear-selective length-frequency sample from an operating model
#'
#' Generates a length-frequency sample for one gear in two stages:
#' \enumerate{
#'   \item a low-effective-sample-size catch-at-age sample is drawn from
#'     the gear's true catch-at-age (via [ca.sim()]);
#'   \item that sampled age composition is projected through the gear's
#'     conditional inverse age-length key,
#'     \eqn{P(l \mid a, \text{caught by gear})} (built once by
#'     [build_gear()], so selectivity is already folded in), and a larger
#'     length sample is drawn from the resulting expected length
#'     distribution.
#' }
#' This mimics a sampling design where a small, possibly more costly,
#' age-reading sample (e.g. observer-collected otoliths) informs the age
#' structure applied to a larger, cheaper length sample (e.g. market or
#' landings-site length measurements) -- rather than treating the two as
#' independent draws from the true population.
#'
#' @param object `FLQuant` of the gear's true catch (or numbers)-at-age,
#'   e.g. one element of `catch_n_gear` from the operating model.
#' @param gear A gear object as returned by [build_gear()]; must contain
#'   `condALK` (age x len matrix, rows sum to 1) and, unless overridden
#'   below, `ess_age`/`ess_len`.
#' @param ess_age,ess_len Integer. Effective sample sizes for the two
#'   stages. Default to `gear$ess_age`/`gear$ess_len`.
#' @param u_age Optional `[ess_age, year, iter]` array of pre-generated
#'   `Uniform(0,1)` deviates for the stage-1 draw (see [rUnif_lfd()]). If
#'   `NULL` (default), stage 1 falls back to `ca.sim()`'s own
#'   `stats::rmultinom()` behaviour (no CRN).
#' @param u_len Optional `[ess_len, year, iter]` array of pre-generated
#'   `Uniform(0,1)` deviates for the stage-2 draw. If `NULL` (default),
#'   stage 2 falls back to a plain `stats::rmultinom()` draw (no CRN).
#' @param scale Logical. If `TRUE` (default), rescale the returned length
#'   sample so its per-year/iter total matches the true total in `object`,
#'   rather than returning raw counts on the `ess_len` scale.
#'
#' @return `FLQuant` with a `len` dimension (matching `gear$len_bins`),
#'   the sampled length frequencies by year and iter.
#'
#' @examples
#' \dontrun{
#' ## trawl <- build_gear(...) from build_gear.R examples
#' ## catch_n_gear$Trawl : FLQuant, true catch-at-age for Trawl, from the OM
#'
#' lfd_trawl <- lfd.sim(catch_n_gear$Trawl, trawl)
#'
#' ## sanity check: shape should show the dome, not just decline with length
#' plot(lfd_trawl)
#'
#' ## CRN-safe version, reproducible across calls
#' devs <- rUnif_lfd(om_gears, years = an(dimnames(catch_n_gear$Trawl)$year),
#'                    nsim = dim(catch_n_gear$Trawl)[6], seed = 456)
#' lfd_trawl_a <- lfd.sim(catch_n_gear$Trawl, trawl,
#'                         u_age = devs$age$Trawl, u_len = devs$len$Trawl)
#' lfd_trawl_b <- lfd.sim(catch_n_gear$Trawl, trawl,
#'                         u_age = devs$age$Trawl, u_len = devs$len$Trawl)
#' stopifnot(identical(lfd_trawl_a, lfd_trawl_b))
#' }
#'
#' @param timing Optional. Within-year timing of the length-sampling stage
#'   only (as fractions of a year added to `t0`), e.g.
#'   `seq(0, 11/12, 1/12)` for roughly monthly, continuous sampling. The
#'   age-sampling stage (stage 1) is NOT split by timing -- it remains a
#'   single annual draw regardless -- on the assumption that age-reading
#'   (e.g. a limited observer programme) and length measurement (e.g.
#'   ongoing market/landings sampling) are typically different sampling
#'   programmes with different temporal coverage. If that assumption
#'   doesn't match your fishery, this needs a different design, not just
#'   a different `timing` value. If `NULL` (default) or length 1, `gear`'s
#'   pre-built `condALK` is used directly, exactly as before -- no ALK
#'   rebuilding, and `params`/`model`/etc. below are unused.
#' @param params Growth parameters, `FLPar(linf, k, t0)`. Required only
#'   if `timing` has length > 1 (the ALK must be rebuilt fresh at each
#'   within-year timing event, since growth advances within the year --
#'   `gear$condALK` alone, built at one fixed timing, cannot be reused
#'   across several).
#' @param model,reflen,cv,lmin,lmax_mult,bin Passed to `iALK()` when
#'   rebuilding the ALK per timing event; ignored if `timing` is `NULL`
#'   or length 1.
#'
#' @export


# lfd_sim_season.R

#' Seasonal length-frequency samples from a season-structured operating model
#'
#' Thin wrapper around [lfd.sim()] that loops over each season of a
#' season-structured `object` (e.g. from [seasonalize_catch_n_gear()]),
#' calling `lfd.sim()` once per season with a within-year `timing` offset
#' of `(s - 0.5) / n_seasons`, and assembling the results into a genuine
#' `season`-dimensioned `FLQuant` (not a list keyed by season number), so
#' `seasonSums()` gives a free round-trip check against the annual total.
#'
#' `ess_age` and `ess_len` each accept either:
#' \itemize{
#'   \item a single value -- interpreted as a fixed annual total, split
#'     evenly across seasons (`rep(round(ess / n_seasons), n_seasons)`),
#'     so the season-summed total matches `gear$ess_age`/`gear$ess_len`
#'     exactly (same invariant as the non-seasonal round-trip check); or
#'   \item a vector of length `n_seasons` -- interpreted as already being
#'     per-season effective sample sizes, used as-is (e.g. `c(20, 10, 50,
#'     30)` for uneven seasonal sampling effort, such as a closed season
#'     or a fleet that's mostly active in one quarter).
#' }
#'
#' Each season's stage-1 age sample is drawn independently (not shared
#' across seasons) -- a real modeling choice (an annual observer
#' programme would typically NOT be redrawn each quarter), not just an
#' implementation detail. If that doesn't match your sampling design,
#' this function needs a different approach (e.g. drawing one annual age
#' sample up front and reusing it across seasons).
#'
#' @param object `FLQuant` with a real `season` dimension (`dim 4 > 1`),
#'   e.g. one element of the list returned by [seasonalize_catch_n_gear()].
#' @param gear A gear object as returned by [build_gear()].
#' @param n_seasons Integer. Number of seasons; defaults to `dim(object)[4]`.
#' @param ess_age,ess_len Single value or length-`n_seasons` vector; see
#'   Details. Default to `gear$ess_age`/`gear$ess_len` (annual total,
#'   split evenly).
#' @param params Growth parameters, `FLPar(linf, k, t0)`. Passed to
#'   [lfd.sim()]; required for the within-season `timing` adjustment.
#' @param ... Passed through to [lfd.sim()] (e.g. `u_age`, `u_len`,
#'   `scale`, `model`, `reflen`, `cv`, `lmin`, `lmax_mult`, `bin`).
#'
#' @return `FLQuant` with a real `season` dimension (`1:n_seasons`).
#'
#' @examples
#' \dontrun{
#' ## even split (default): total ess_len matches gear$ess_len
#' lfd.sim.season(obj, trawl, n_seasons = 4, params = lhpars)
#'
#' ## uneven seasonal effort, e.g. closed season in Q1
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


# {{{
# rUnif_lfd 
#
#' Pre-generate CRN uniform-deviate streams for a two-stage length OEM
#'
#' Generates one `Uniform(0,1)` array per gear per sampling stage (age,
#' length), each shaped `[draw, year, iter]`, for use with
#' [rmultinom_crn()]. Intended to be called once per MSE run/scenario set
#' and reused everywhere the length OEM draws a sample, so the same
#' underlying randomness is shared across MP variants (CRN).
#'
#' @param om_gears Named list of gear objects as returned by
#'   [build_gear()]; must each contain `ess_age` and `ess_len`.
#' @param years Integer vector of years to generate deviates for
#'   (typically the full projection horizon).
#' @param nsim Integer. Number of OM iterations.
#' @param seed Optional integer seed, set via `set.seed()` before
#'   generation for reproducibility across sessions. If `NULL` (default),
#'   the current RNG state is used as-is.
#'
#' @return A list with two elements, `age` and `len`, each a named list
#'   (one array per gear) of dimension `[ess, length(years), nsim]`.
#'
#' @examples
#' \dontrun{
#' om_gears <- list(Trawl = trawl, Gillnet = gillnet)   # from build_gear()
#' devs <- rUnif_lfd(om_gears, years = 2025:2054, nsim = 100, seed = 456)
#'
#' dim(devs$age$Trawl)     # 50 x 30 x 100  (ess_age x years x nsim)
#' dim(devs$len$Gillnet)   # 300 x 30 x 100 (ess_len x years x nsim)
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
#}}}

# seasonalize_catch_n_gear.R
# generalised to multiple gears sharing one total Z, and to use the OM's
# own harvest(stock)+m(stock) rather than re-summing f_age_gear (which
# only shapes each gear's share of catch, not the total decay rate).
#
# Stores output as a genuine FLQuant with the 'season' dimension populated
# (not a list keyed by "season1", "season2", ...), so seasonSums() gives a
# free round-trip check against the annual truth, and any FLQuant-aware
# plotting/aggregation works without custom list-handling.

#' Disaggregate an annual OM's gear catch-at-age into within-year seasons
#'
#' Splits `stock`'s total instantaneous mortality (`harvest(stock) +
#' m(stock)`) evenly across `n_seasons`, propagates numbers-at-age through
#' each season via standard exponential decay, and allocates catch to
#' each gear proportional to that gear's share of total F
#' (`f_age_gear[[g]] / z_total`) -- `f_age_gear` shapes *which* catch goes
#' to which gear; it does not drive the survival/decay rate, which comes
#' from the OM's own already-validated total Z.
#'
#' @param stock `FLStock` (or `FLStockR`), the OM's true annual stock --
#'   `stock.n(stock)`, `harvest(stock)`, `m(stock)` are used directly.
#' @param f_age_gear Named list of `FLQuant`s, one per gear, the OM's true
#'   annual fishing-mortality-at-age by gear (same dimensions as
#'   `harvest(stock)`; `Reduce("+", f_age_gear)` should be close to
#'   `harvest(stock)` -- worth checking directly if this function's
#'   output doesn't round-trip cleanly).
#' @param n_seasons Integer. Number of within-year seasons.
#'
#' @return Named list, one `FLQuant` per gear, each with a real `season`
#'   dimension of length `n_seasons` (all other dimensions matching
#'   `stock.n(stock)`). `seasonSums(out[[g]])` should closely match the
#'   true annual `catch.n` for that gear.
#'
#' @examples
#' \dontrun{
#' seasonal_catch <- seasonalize_catch_n_gear(om_hist, f_age_gear, n_seasons = 4)
#'
#' ## round-trip check: should closely match the true annual catch per gear
#' plot(seasonSums(seasonal_catch$Trawl) - catch_n_gear$Trawl)
#' }
#'
#' @export
seasonalize_catch_n_gear <- function(stock, f_age_gear, n_seasons) {
  
  gears <- names(f_age_gear)
  stock_n <- stock.n(stock)
  z_total <- harvest(stock) + m(stock)   # OM's own aggregate Z
  
  mk_season_template <- function(x) {
    d <- dimnames(x)
    d$season <- ac(seq_len(n_seasons))
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

# plot_lfd_season_ggplot.R
#' Season-coloured length-frequency diagnostic, one panel per gear
#'
#' @param lfd_season Named list of `FLQuant`s, one per gear, each with a
#'   populated `season` dimension (e.g. from [lfd.sim.season()]).
#' @param year Character or numeric. Which year to plot.
#' @param iter Integer. Which iteration to plot. Default `1`.
#'
#' @return A `ggplot` object: length on x, sampled count on y, one line
#'   per season (coloured), faceted by gear.
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
# plot_lfd_ridgeline.R

#' Ridgeline plot of length-frequency samples over time (annual or seasonal)
#'
#' Plots one or more length-frequency `FLQuant`s (e.g. sampled LFDs from
#' [lfd.sim()]/[lfd.sim.season()], keyed by gear) as ridges along a
#' continuous time axis via `geom_ribbon()` + `coord_flip()`, faceted by
#' panel (gear).
#'
#' If any panel carries a real `season` dimension, each season's ridge is
#' placed at its actual fractional-year position (`year + (season - 0.5)
#' / n_seasons`), so seasons appear in chronological order within each
#' year rather than as separate facet rows, and are distinguished by
#' fill colour rather than by expanding the grid.
#'
#' Cohort-tracking VBGF reference lines are optional: pass `lhpar` to
#' overlay them, or leave `NULL` (default) for ridges only.
#'
#' @param panels Named list of `FLQuant`s (names become facet/panel
#'   labels, typically gear names), each with `len` and `year`
#'   dimensions, single iteration. May carry a real `season` dimension.
#' @param lhpar Optional `FLPar`/named vector with `linf`, `k`, `t0` for
#'   dashed diagonal cohort reference lines. `NULL` (default) skips them.
#' @param scale Numeric. Controls how far each ridge protrudes, relative
#'   to the spacing between time points (a year, or a season-slot if
#'   seasonal). Default `0.9`.
#' @param cohort_pad Integer. Extra cohorts drawn beyond the plotted time
#'   range at each end (only used when `lhpar` is supplied). Default `15`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' ## annual, one ridge column per gear
#' plot_lfd_ridgeline(lapply(lfds, function(x) iter(x, test_iter)),lhpar=lhpars)
#'
#' ## seasonal, colour-coded by season, same layout
#' plot_lfd_ridgeline(lapply(lfd_season, function(x) iter(x, test_iter)))
#' }
#'
#' @export
plot_lfd_ridgeline <- function(panels, lhpar = NULL, scale = 0.9, cohort_pad = 15,
                               annual = FALSE) {
  
  if (annual) {
    panels <- lapply(panels, seasonSums)
  }
  
  df <- do.call(rbind, lapply(names(panels), function(p) {
    d <- as.data.frame(panels[[p]], cohort = FALSE)
    d$len <- as.numeric(as.character(d$len))
    d$year <- as.numeric(as.character(d$year))
    d$panel <- p
    d
  }))
  
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
  
  p <- ggplot2::ggplot(df)
  
  if (has_season) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(x = len, ymin = time, ymax = time + dens, group = grp, fill = season),
      alpha = 0.8, colour = NA
    ) + ggplot2::scale_fill_brewer(palette = "Set2", name = "Season")
  } else {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(x = len, ymin = time, ymax = time + dens, group = grp),
      fill = "#5B9BD5", alpha = 0.75, colour = NA
    )
  }
  
  p <- p +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~panel, ncol = 1) +
    ggplot2::labs(x = "Length (cm)", y = "Year") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  
  if (!is.null(lhpar)) {
    yr_range <- range(df$time)
    len_range <- range(df$len)
    linf <- c(lhpar["linf"]); k <- c(lhpar["k"]); t0 <- c(lhpar["t0"])
    cohorts <- seq(floor(yr_range[1]) - cohort_pad, ceiling(yr_range[2]) + cohort_pad)
    age_grid <- seq(0, -log(1 - 0.98) / k, by = 0.1)
    
    ref <- do.call(rbind, lapply(cohorts, function(cy) {
      L <- linf * (1 - exp(-k * (age_grid + t0)))
      data.frame(t = age_grid + cy, L = L, cohort = cy)
    }))
    ref <- ref[ref$t >= yr_range[1] & ref$t <= yr_range[2] &
                 ref$L >= len_range[1] & ref$L <= len_range[2], ]
    
    p <- p + ggplot2::geom_line(
      data = ref, ggplot2::aes(x = L, y = t, group = cohort),
      colour = "grey50", linetype = 2, linewidth = 0.3, inherit.aes = FALSE
    )
  }
  
  p
}