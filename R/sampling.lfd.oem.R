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
