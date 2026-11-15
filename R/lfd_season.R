# lfd_season.R - within-year length sampling designs (continuous fishery
# sampling, pulse surveys) from an annual operating model.
# FLRef/R/lfd_season.R


#' Split a sample size exactly among seasons
#'
#' Largest-remainder split: the seasonal values differ by at most one and sum
#' exactly to `n` (unlike `rep(round(n / S), S)`, which gives 48 for 50 / 12).
#'
#' @param n Integer total sample size.
#' @param S Integer number of seasons.
#' @return Integer vector of length `S` summing to `n`.
#' @examples
#' split_ess(50, 12)
#' sum(split_ess(50, 12))
#' @export
split_ess <- function(n, S) {
  n <- as.integer(round(n))
  rep(n %/% S, S) + as.integer(seq_len(S) <= n %% S)
}


#' Within-year distribution of the annual catch over equal seasons
#'
#' Under constant instantaneous F and M within the year, the fraction of the
#' annual catch-at-age taken in season \eqn{s} of \eqn{S} equal seasons is
#' \deqn{\pi_{a,y,s}=\frac{e^{-Z_{a,y}(s-1)/S}\left(1-e^{-Z_{a,y}/S}\right)}
#'                        {1-e^{-Z_{a,y}}},\qquad \sum_s\pi_{a,y,s}=1,}
#' (\eqn{1/S} where \eqn{Z=0}). Multiplying annual gear catch-at-age (from
#' [catch_n_gear()]) by \eqn{\pi} gives seasonal catch-at-age that sums
#' exactly to the annual catch, with older and heavily fished ages
#' concentrated early in the year. This is the observation-side analogue of
#' [seasonalize_catch_n_gear()] and needs no seasonal operating model.
#'
#' @param stock `FLStock` with annual `harvest` and `m` (one season).
#' @param n_seasons Integer number of equal seasons (e.g. 4 or 12).
#' @return `FLQuant` age x year x unit x season x area x iter of fractions.
#' @examples
#' \dontrun{
#' data(ple4)
#' fr <- season_fraction(ple4[, ac(2015:2017)], 12)
#' range(seasonSums(fr))  # 1
#' }
#' @export
season_fraction <- function(stock, n_seasons) {
  Z <- harvest(stock) + m(stock)
  if (dim(Z)[4] != 1)
    stop("season_fraction(): 'stock' must have a single (annual) season.")
  d <- dimnames(Z)
  d$season <- ac(seq_len(n_seasons))
  out <- FLQuant(NA_real_, dimnames = d, units = "")
  den <- 1 - exp(-Z)
  for (s in seq_len(n_seasons)) {
    fr <- exp(-Z * (s - 1) / n_seasons) * (1 - exp(-Z / n_seasons)) / den
    fr[!is.finite(fr)] <- 1 / n_seasons
    out[, , , s] <- fr
  }
  out
}


#' Survey numbers-at-age at the time of a pulse survey
#'
#' \deqn{I_{a,y}=N_{a,y}\,e^{-Z_{a,y}t_g}\,s_{g,a},}
#' the at-age composition a survey gear samples at its `timing` \eqn{t_g}.
#' Passed to [lfd.sim()] with the same gear, whose conditional ALK
#' \eqn{P(l\mid a,g)} is built at the same timing.
#'
#' @param stock `FLStock`.
#' @param gear Gear object from [build_gear()] (uses `sel_a`, `grid$timing`).
#' @param years Years. Default all years of `stock`.
#' @param timing Within-year timing. Default `gear$grid$timing`.
#' @return `FLQuant` of selected numbers-at-age.
#' @export
survey_n_gear <- function(stock, gear, years = dimnames(stock)$year,
                          timing = gear$grid$timing %||% 0) {
  stk <- stock[, ac(years)]
  n <- stock.n(stk) * exp(-(harvest(stk) + m(stk)) * timing)
  s <- n
  s[] <- c(gear$sel_a)                      # recycled over year/unit/iter
  n * s
}


#' Simulate length samples for a set of gears (one area or sex)
#'
#' The sampling step of [sampling.lfd.oem()] for one `unit`, exported so that
#' sampling designs can be built and checked outside an MSE with exactly the
#' code the OEM runs. For each gear in `om_gears`:
#' \describe{
#'   \item{fishery, `n_seasons = 1`}{gear catch-at-age by [catch_n_gear()],
#'     one two-stage sample with [lfd.sim()] using the gear's conditional
#'     ALK (built at `gear$grid$timing`).}
#'   \item{fishery, `n_seasons = S > 1`}{gear catch-at-age split into
#'     \eqn{S} seasons by [season_fraction()]; an independent two-stage
#'     sample per season with [lfd.sim.season()] (ALK at each season
#'     mid-point, `ess_age` and `ess_len` split among seasons). This emulates
#'     continuous sampling of the fishery.}
#'   \item{survey}{survey numbers-at-age at the gear's timing,
#'     [survey_n_gear()], one two-stage sample with [lfd.sim()].}
#' }
#'
#' @param stock `FLStock` with a single `unit` (subset areas or sexes first).
#' @param om_gears Named list of gears from [build_gear()].
#' @param years Years to sample.
#' @param u_lfd Optional common random numbers for this unit, as returned by
#'   [rUnif_lfd()] (spanning at least `years`); `NULL` uses the live RNG.
#' @param f_share Optional `FLQuants` of time-varying F shares for this unit
#'   (see [sampling.lfd.oem()]).
#' @param season Logical. For gears with `n_seasons > 1`: `TRUE` keeps the
#'   seasonal samples in the season dimension; `FALSE` (default) returns the
#'   annual sample, pooled with `seasonSums()`. Gears sampled once a year
#'   always have one season.
#' @return A list with `lfd` (`FLQuants` of raw sampled counts by gear,
#'   `scale = FALSE`) and `catch_n` (`FLQuants` of annual catch-at-age for the
#'   fishery gears).
#' @examples
#' \dontrun{
#' gears <- build_gear(list(
#'   Commercial = list(type = "dnormal", Lpeak = 30, sd_left = 4, sd_right = 15,
#'                     f_mult = 1, ess_age = 120, ess_len = 1200, n_seasons = 12),
#'   Survey = list(type = "logistic", L50 = 15, L95 = 20, obs = "survey",
#'                 ess_age = 30, ess_len = 200, timing = 1/24)),
#'   lhpar = lhpar, age = 1:10, timing = 0.5, fbar_range = c(2, 6))
#' s <- sample_lfd_gears(ple4, gears, years = 2008:2017, season = TRUE)
#' plot_lfd_ridgeline(s$lfd, lhpar = lhpar)
#' }
#' @export
sample_lfd_gears <- function(stock, om_gears, years = dimnames(stock)$year,
                             u_lfd = NULL, f_share = NULL, season = FALSE) {

  years <- ac(years)
  if (dim(stock)[3] != 1)
    stop("sample_lfd_gears(): 'stock' must have one unit; loop over units.")
  stk <- stock[, years]
  its <- dims(stk)$iter

  obs  <- vapply(om_gears, function(g) g$obs %||% "fishery", character(1))
  fish <- names(om_gears)[obs == "fishery"]

  ## annual gear catch-at-age (shares of the fishery gears sum to one)
  cn <- if (length(fish)) {
    sh <- if (is.null(f_share)) om_gears[fish] else
      FLQuants(lapply(f_share[fish], function(x) x[, years]))
    catch_n_gear(stk, sh, years = years, check = is.null(f_share))
  } else FLQuants()

  sl <- function(u) {
    if (is.null(u)) return(NULL)
    if (is.list(u)) return(lapply(u, sl))
    if (!all(years %in% dimnames(u)$year))
      stop("Common random numbers do not cover years ",
           paste(range(years), collapse = "-"), ".")
    u[, years, , drop = FALSE]
  }
  guard0 <- function(x) { x[x <= 0 | is.na(x)] <- sqrt(.Machine$double.eps); x }

  lfd <- lapply(setNames(nm = names(om_gears)), function(g) {
    gr    <- om_gears[[g]]
    S     <- gr$n_seasons %||% 1
    u_age <- sl(u_lfd$age[[g]])
    u_len <- sl(u_lfd$len[[g]])

    if (obs[[g]] == "survey") {
      x <- survey_n_gear(stk, gr, years)
      return(lfd.sim(guard0(x), gr, u_age = u_age, u_len = u_len, scale = FALSE))
    }

    x <- cn[[g]]
    if (sum(x, na.rm = TRUE) == 0)                     # no catch: zero sample
      return(lfd_template(gr, years, its, dimnames(stk)$unit,
                          season = season))
    if (S <= 1)
      return(lfd.sim(guard0(x), gr, u_age = u_age, u_len = u_len, scale = FALSE))

    xs <- x %*% season_fraction(stk, S)                # seasonal catch-at-age
    lfd.sim.season(guard0(xs), gr, n_seasons = S, u_age = u_age, u_len = u_len,
                   season = season, scale = FALSE)
  })

  list(lfd = FLQuants(lfd), catch_n = cn)
}


#' Empty length-sample template for a gear
#'
#' @param gear Gear object from [build_gear()].
#' @param years,nsim Years and number of iterations.
#' @param units `unit` dimnames (e.g. `"unique"`, areas or sexes).
#' @param season Logical; give gears with `n_seasons > 1` a season
#'   dimension (as [sample_lfd_gears()] with `season = TRUE`).
#' @return `FLQuant` of zeros, len x year x unit x season x area x iter.
#' @export
lfd_template <- function(gear, years, nsim, units = "unique", season = FALSE) {
  S <- gear$n_seasons %||% 1
  FLQuant(0, dimnames = list(
    len = colnames(gear$condALK), year = ac(years), unit = units,
    season = if (season && S > 1) ac(seq_len(S)) else "all",
    area = "unique", iter = seq_len(nsim)), units = "cm")
}


#' Initialise the LFD observations of a sampling.lfd.oem() FLoem
#'
#' Builds `observations$lfd` and `observations$catch_gear` for
#' [sampling.lfd.oem()] over `years` (typically `y0:fy`) and fills the
#' historical years with samples drawn by [sample_lfd_gears()], i.e. with the
#' same code and, when `u_lfd` is given, the same common random numbers the
#' OEM will use. Templates have the season dimension implied by `season`, so
#' they always match the OEM's output.
#'
#' @param stock OM `FLStock` (all units).
#' @param om_gears Gears: a named gear list (applied to every unit) or a list
#'   by unit of gear lists, as for [sampling.lfd.oem()].
#' @param years All years of the observation series.
#' @param hist_years Years to fill with historical samples.
#' @param u_lfd Common random numbers as for `deviances$lfd` of
#'   [sampling.lfd.oem()].
#' @param season,f_share As in [sampling.lfd.oem()].
#' @return A list with `lfd` and `catch_gear` (`FLQuants` by gear).
#' @export
lfd_obs_init <- function(stock, om_gears, years, hist_years, u_lfd = NULL,
                         season = FALSE, f_share = NULL) {

  units <- dimnames(stock)$unit
  om_gears <- .lfd_gears_by_unit(om_gears, units)
  u_lfd <- .lfd_crn_by_unit(u_lfd, units)
  flat <- .lfd_flat_gears(om_gears)
  nsim <- dims(stock)$iter
  hy <- ac(hist_years)

  lfd <- FLQuants(lapply(flat, lfd_template, years = years, nsim = nsim,
                         units = units, season = season))
  fish <- names(flat)[vapply(flat, function(g)
    (g$obs %||% "fishery") == "fishery", logical(1))]
  cg <- FLQuants(lapply(setNames(nm = fish), function(g)
    FLQuant(0, dimnames = list(year = ac(years), unit = units,
                               iter = seq_len(nsim)))))

  for (u in units) {
    fs <- if (is.null(f_share)) NULL else
      FLQuants(lapply(f_share, function(x) x[, , u]))
    sm <- sample_lfd_gears(stock[, hy, u], om_gears[[u]], years = hy,
                           u_lfd = u_lfd[[u]], f_share = fs, season = season)
    for (g in names(sm$lfd)) lfd[[g]][, hy, u] <- sm$lfd[[g]]
    for (g in names(sm$catch_n))
      cg[[g]][, hy, u] <- quantSums(sm$catch_n[[g]] * catch.wt(stock)[, hy, u])
  }
  list(lfd = lfd, catch_gear = cg)
}


#' Compare length-sampling designs for one year
#'
#' Overlays the annual length composition (proportions) obtained under
#' different sampling designs, e.g. one mid-year sample against monthly
#' sampling pooled over the year, optionally with the expected composition.
#' Seasonal samples are pooled with `seasonSums()`.
#'
#' @param designs Named list of length-sample `FLQuant`s (same length bins).
#' @param year Year to plot.
#' @param iter Iteration. Default `1`.
#' @param expected Optional named numeric vector of expected proportions
#'   (names = length bins), drawn as a dashed line.
#' @return A `ggplot` object.
#' @export
plot_lfd_designs <- function(designs, year, iter = 1, expected = NULL) {
  df <- do.call(rbind, lapply(names(designs), function(d) {
    x <- designs[[d]]
    if (dim(x)[4] > 1) x <- seasonSums(x)
    x <- x[, ac(year), , , , iter]
    v <- c(x) / sum(x)
    data.frame(design = d, len = an(dimnames(x)$len), p = v)
  }))
  df$design <- factor(df$design, levels = names(designs))
  p <- ggplot2::ggplot(df, ggplot2::aes(len, p, colour = design)) +
    ggplot2::geom_step(linewidth = 0.7, direction = "mid") +
    ggplot2::labs(x = "Length (cm)", y = "Proportion", colour = NULL,
                  title = paste("Annual length composition,", year)) +
    ggplot2::theme_bw()
  if (!is.null(expected)) {
    p <- p + ggplot2::geom_line(
      data = data.frame(len = an(names(expected)), p = expected),
      ggplot2::aes(len, p), inherit.aes = FALSE, linetype = 2)
  }
  p
}


# Internal helpers shared by sampling.lfd.oem() and lfd_obs_init() -----------

.lfd_is_gear <- function(x) is.list(x) && !is.null(x$condALK)

# gear list or unit -> gear list  ->  unit -> gear list
.lfd_gears_by_unit <- function(om_gears, units) {
  if (all(vapply(om_gears, .lfd_is_gear, logical(1))))
    om_gears <- setNames(rep(list(om_gears), length(units)), units)
  if (!all(units %in% names(om_gears)))
    stop("names(om_gears) must match the stock's unit dimnames: ",
         paste(units, collapse = ", "))
  om_gears[units]
}

# rUnif_lfd() output or unit -> rUnif_lfd() output  ->  unit -> output
.lfd_crn_by_unit <- function(u_lfd, units) {
  if (is.null(u_lfd)) return(setNames(vector("list", length(units)), units))
  if (all(c("age", "len") %in% names(u_lfd))) {
    if (length(units) > 1)
      stop("With several units, deviances$lfd must be a list by unit of ",
           "rUnif_lfd() outputs (independent streams per unit).")
    u_lfd <- setNames(list(u_lfd), units)
  }
  u_lfd[units]
}

# unit -> gear list  ->  one gear list with every gear name once
.lfd_flat_gears <- function(om_gears) {
  gl <- unlist(unname(om_gears), recursive = FALSE)
  gl[!duplicated(names(gl))]
}
