# Quick checks for the within-year LFD sampling code ----------------------------
# Run after installing FLRef:  source("inst/scripts/test_lfd_season.R")
# Each check stops with an error if it fails.

library(FLCore); library(FLRef)
data(ple4)

ages  <- an(dimnames(ple4)$age)
yrs   <- ac(2013:2017)
lhpar <- FLPar(linf = 55, k = 0.15, t0 = -0.5, a = 8.9e-06, b = 3.05)

mk <- function(n_seasons, ess_age = 120, ess_len = 1200)
  build_gear(list(
    Commercial = list(type = "dnormal", Lpeak = 28, sd_left = 4, sd_right = 15,
                      f_mult = 1, ess_age = ess_age, ess_len = ess_len,
                      n_seasons = n_seasons),
    Survey = list(type = "logistic", L50 = 12, L95 = 18, obs = "survey",
                  ess_age = 30, ess_len = 200, timing = 1/24)),
    lhpar = lhpar, age = ages, timing = 0.5, fbar_range = c(2, 6))

g1  <- mk(1)
g12 <- mk(12)

## build_gear: design entries and per-gear timing --------------------------
stopifnot(g12$Commercial$n_seasons == 12, g12$Survey$obs == "survey",
          g12$Survey$f_mult == 0, g12$Survey$grid$timing == 1/24,
          g12$Commercial$grid$timing == 0.5)
stopifnot(all(abs(c(g12$Commercial$f_share) - 1) < 1e-12 |
              c(g12$Commercial$sel_a) == 0))
## old-style gear_cfg (no design entries) still works
g_old <- build_gear(list(T = list(type = "logistic", L50 = 20, L95 = 25,
                                  f_mult = 1, ess_age = 50, ess_len = 150)),
                    lhpar = lhpar, age = ages, timing = 0.5, fbar_range = c(2, 6))
stopifnot(g_old$T$n_seasons == 1, g_old$T$obs == "fishery")

## split_ess / season_fraction ---------------------------------------------
stopifnot(sum(split_ess(50, 12)) == 50, max(split_ess(50, 12)) - min(split_ess(50, 12)) <= 1)
fr <- season_fraction(ple4[, yrs], 12)
stopifnot(max(abs(seasonSums(fr) - 1)) < 1e-12)

## rUnif_lfd: seasonal streams sized by split_ess ---------------------------
u12 <- rUnif_lfd(g12, years = yrs, nsim = 2, seed = 1)
stopifnot(is.list(u12$len$Commercial), length(u12$len$Commercial) == 12,
          sum(sapply(u12$len$Commercial, nrow)) == 1200,
          is.array(u12$len$Survey))
## annual gears: identical deviates to the previous implementation
u1 <- rUnif_lfd(g1, years = yrs, nsim = 2, seed = 1)
set.seed(1); a1 <- runif(120 * length(yrs) * 2)
stopifnot(identical(c(u1$age$Commercial), a1))

## sample_lfd_gears: totals, season switch, CRN reproducibility -------------
stk2 <- propagate(ple4, 2)
sF <- sample_lfd_gears(stk2, g12, yrs, u_lfd = u12, season = FALSE)
sT <- sample_lfd_gears(stk2, g12, yrs, u_lfd = u12, season = TRUE)
stopifnot(dim(sF$lfd$Commercial)[4] == 1, dim(sT$lfd$Commercial)[4] == 12,
          dim(sT$lfd$Survey)[4] == 1)
stopifnot(all(c(quantSums(sF$lfd$Commercial)) == 1200),
          all(c(quantSums(sF$lfd$Survey)) == 200))
stopifnot(identical(c(seasonSums(sT$lfd$Commercial)), c(sF$lfd$Commercial)))
stopifnot(all(c(quantSums(sT$lfd$Commercial)[, , , 1]) == split_ess(1200, 12)[1]))
sF2 <- sample_lfd_gears(stk2, g12, yrs, u_lfd = u12, season = FALSE)
stopifnot(identical(sF, sF2))

## lfd_sop with seasonal samples: SOP of all seasons equals annual catch -----
cw  <- catch_gear(sT$catch_n, stk2)$Commercial
rs  <- lfd_sop(sT$lfd$Commercial, g12$Commercial, cw)
w   <- rs; w[] <- g12$Commercial$wt_len[dimnames(rs)$len]
stopifnot(max(abs(c(seasonSums(quantSums(rs * w)) / cw) - 1)) < 1e-8)

## lfd_obs_init: templates match the season setting -------------------------
ini <- lfd_obs_init(stk2, g12, years = 2010:2020, hist_years = yrs,
                    u_lfd = rUnif_lfd(g12, 2010:2020, 2, seed = 2), season = TRUE)
stopifnot(dim(ini$lfd$Commercial)[4] == 12, dim(ini$lfd$Survey)[4] == 1,
          names(ini$catch_gear) == "Commercial",
          all(c(quantSums(seasonSums(ini$lfd$Commercial[, yrs]))) == 1200),
          all(c(ini$lfd$Commercial[, "2010"]) == 0))

message("All within-year LFD checks passed.")
