# Developer notes: length-frequency simulation in FLRef

Companion to `inst/rmd/lfd-01-single-sex.Rmd`. Covers the documentation set-up,
issues found in the current code while writing the guideline, and the design
for sex-structured (`unit.type = "gender"`) and spatial (`unit.type = "spatial"`)
stocks.

## 1. Documentation set-up (as in cpueTargeting)

```
FLRef/
  inst/rmd/lfd-01-single-sex.Rmd     # sources
  inst/rmd/lfd-02-sex-structured.Rmd # (next)
  inst/rmd/lfd-03-spatial.Rmd        # (next, ICES deep-sea WG)
  inst/scripts/render_docs.R         # renders inst/rmd -> docs/
  docs/index.html                    # landing page (hand-written)
  docs/.nojekyll
  rdata/                             # cached MSE runs (not pushed)
```

* `.gitignore`: add `rdata/`, `inst/rmd/*.html`, `inst/rmd/rdata/`.
* `.Rbuildignore`: add `^docs$`, `^rdata$`.
* GitHub: Settings > Pages > Deploy from a branch > `main` / `/docs`.
* `DESCRIPTION` Suggests: `rmarkdown`, `knitr` (rendering only).
* Render from the package root: `source("inst/scripts/render_docs.R")`.
* The Rmd needs FLRef installed with the current `R/` (it calls the exported
  functions), so `make install` (or `devtools::install()`) first.

## 2. Issues found in the current code

1. **`sampling.lfd.oem()` is defined twice**: in `R/lensim.R` (line ~1370)
   and in `R/sampling.lfd.oem.R`. They are identical today; with the default
   alphabetical collation the second wins, so edits to the copy in
   `lensim.R` would silently have no effect. Keep `R/sampling.lfd.oem.R`
   and delete the block from `lensim.R` (roxygen + function).

2. **Common random numbers are not split by `mp()` in parallel mode.**
   `mp()` with several workers calls `iter(oem, j)`, which subsets only
   `FLQuant` deviances (`rapply(..., classes = "FLQuant")`). The arrays from
   `rUnif_lfd()` are plain arrays, so every worker would index iterations
   `1:length(j)` of the full array, i.e. all chunks reuse the same deviates.
   The observation `lfd`/`catch_gear` FLQuants are also passed through
   ``lapply(..., `[`, iter = j)``, which needs checking for `FLQuants`.
   Until fixed, run with `parallel = FALSE` (the guideline does).
   Fix: let `rUnif_lfd()` return `FLQuant`s with quant `"draw"`
   (`draw x year x unit x season x area x iter`); in `ca.sim()`, `lfd.sim()`
   and `len.sim()` use `c(u[, y, , , , i])` in `rmultinom_crn()`.

3. **`ALKs()` with more than one iteration.** `alk[] = object[, ac(x)]`
   fills an `age x len x iter` FLPar from an `age x iter` vector by
   recycling, which places iteration 2 in length bin 2, etc. Correct for one
   iteration (as used in the guideline); for several, fill per iteration:
   `iter(alk, i)[] <- c(iter(object[, ac(x)], i))`.

4. **Survey gears** (resolved in this version): `build_gear()` accepts
   `obs = "survey"` (implies `f_mult = 0`); `sample_lfd_gears()`, and hence
   `sampling.lfd.oem()`, samples surveys from
   \(N_{a,y}e^{-Z_{a,y}t_g}s_{g,a}\) at the gear's timing
   (`survey_n_gear()`). No `catch_gear` is reported for surveys.

5. **One `timing` per `build_gear()` call** (resolved): a `timing` entry in
   a `gear_cfg` element overrides the call-level default.

6. **`mp()` simplifies the OM stock to the observation's dimensions**
   (`simplify(stock(om), ...)` where `observations$stk` has one unit). For
   sex- or area-structured OEMs, `observations$stk` must keep its `unit`
   dimension or `sampling.lfd.oem()` never sees the units.

7. Candidate helper for FLRef: `fit_gear_sel()` from the guideline (fits
   length selectivity so that the ALK-integrated `sel_a` matches a target
   age selectivity, e.g. `catch.sel(stk)` or survey catchability).

## 3. Within-year sampling designs (this version)

Code changes (files `R/lfd_season.R` new; `R/lensim.R`, `R/sel_la.R`,
`R/sampling.lfd.oem.R` changed; `NAMESPACE`):

| function | change |
|----------|--------|
| `build_gear()` | optional per-gear `obs`, `timing`, `n_seasons`; stored in the gear (`gear$obs`, `gear$n_seasons`, `gear$grid$timing`). Old `gear_cfg` lists work unchanged. |
| `split_ess()` | new, exact largest-remainder split of a sample size among seasons |
| `season_fraction()` | new, \(\pi_{a,y,s}\): share of the annual catch-at-age taken in each of \(S\) equal seasons under constant F and M |
| `survey_n_gear()` | new, survey numbers-at-age at the gear timing |
| `lfd.sim.season()` | per-season CRN (`u_age`, `u_len` lists), exact ESS split, `params = gear$lhpar` default, `season = TRUE/FALSE`, keeps unit/area dimnames |
| `rUnif_lfd()` | gears with `n_seasons > 1` get a list of seasonal arrays; annual gears get the same deviates as before |
| `lfd_sop()` | seasonal samples raised to annual catch with one factor per year |
| `sample_lfd_gears()` | new, the per-unit sampling step of the OEM (fishery annual / within-year, survey), exported so designs can be checked outside an MSE with the OEM's own code |
| `lfd_template()`, `lfd_obs_init()` | new, observation templates and historical seeding with the right season dimension and the OEM's CRN |
| `plot_lfd_designs()` | new, compares pooled annual compositions of several designs |
| `sampling.lfd.oem()` | calls `sample_lfd_gears()` per unit; new `season = FALSE`; surveys; seasonal raising; informative error if `observations$lfd` has the wrong season dimension |

`inst/scripts/test_lfd_season.R` runs quick checks of all of the above.

Design notes:

* The OM stays annual. Seasons exist only on the observation side, so any
  annual `FLom` can be sampled monthly; a seasonal OM is not required.
* `season = TRUE` is meant for a seasonal assessment (FLicc with seasonal
  transitions): season \(s\) of the sample corresponds to timing
  \((s-0.5)/S\). Survey samples keep one season and carry their timing in
  the gear (`gear$grid$timing`), so the estimator should read timings from
  the gears rather than infer them from the season index.
* Different gears can have different `n_seasons`; `observations$lfd` holds
  one `FLQuant` per gear, so their season dimensions need not agree.
* Sampling effort is split equally among seasons (or by a user vector in
  `lfd.sim.season()`). Allocation proportional to seasonal catch would need
  year-varying sample sizes; possible later via `ess_alloc = "catch"`.

## 4. Sex-structured stocks: `unit.type = "gender"`

### Life-history parameters: `FLPar` with a `unit` dimension

Prefer one `FLPar` with dimensions `params x unit x iter` over `FLPars`:

* `lhpar[, "F"]` mirrors `stk[, , "F"]`, so the unit names are the link
  between biology and the stock object (one check:
  `identical(dimnames(lhpar)$unit, dimnames(stk)$unit)`);
* iterations propagate in one object (`propagate(lhpar, nsim)`), which
  `FLPars` does not do as a whole;
* one argument passes through `build_gear()`, `lfd.sim()`, `len.sim()`.

```r
lhpar <- FLPar(array(c(60, 0.10, -0.5, 8.9e-6, 3.05,    # F
                       45, 0.15, -0.5, 8.9e-6, 3.05),   # M
                     dim = c(5, 2, 1),
                     dimnames = list(params = c("linf", "k", "t0", "a", "b"),
                                     unit = c("F", "M"), iter = 1)))
lhpar[, "F"]
```

Note that the current helpers extract `c(params["linf"])`, which with a
`unit` dimension returns two values; all of them must be called per unit
(`lhpar[, u]`), never with the full object.

### `build_gear()`

* If `"unit" %in% names(dimnames(lhpar))`: build each unit with
  `lhpar[, u]` and return a list **unit -> gear** (the same nesting as the
  area -> gear list already accepted by `sampling.lfd.oem()`), with
  `attr(, "unit.type") <- "gender"`.
* Selectivity-at-length \(s_g(l)\) is shared; \(s_{g,a,u}\), `condALK` and
  `f_share` differ by sex through growth.
* **Common length grid**: `len_grid()` depends on `linf`, so the sexes would
  get different grids. Build all units on the grid of the largest
  \(L_\infty\) (`gear_ialk()`/`iALK(max = )` already accept an explicit
  upper limit; `build_gear()` needs a `linf_grid` or `lmax` argument).
* `ess_age`, `ess_len` stay per gear (not per sex).

### Sampling: one event per gear on the combined catch

For gear \(g\) the sampled fish come from the combined catch of both sexes:
the age stage draws `ess_age` fish over the (sex, age) cells and the length
stage draws `ess_len` fish over the (sex, length) cells. This can reuse
`lfd.sim()` unchanged by stacking the sexes into one pseudo-age dimension
with a block-diagonal conditional key (sketch, untested):

```r
lfd.sim.gender <- function(object, gears, ess_age, ess_len,
                           u_age = NULL, u_len = NULL, scale = FALSE) {
  un <- dimnames(object)$unit
  na <- dim(object)[1]
  nl <- ncol(gears[[un[1]]]$condALK)
  lens <- colnames(gears[[un[1]]]$condALK)
  nu <- length(un)
  # stacked catch-at-age, pseudo-ages 1..nu*na (lfd.sim() needs numeric ages)
  st <- FLQuant(0, dimnames = list(age = seq_len(nu * na),
                                   year = dimnames(object)$year,
                                   iter = dimnames(object)$iter))
  K <- matrix(0, nu * na, nu * nl,
              dimnames = list(age = seq_len(nu * na), len = seq_len(nu * nl)))
  for (k in seq_len(nu)) {
    ia <- (k - 1) * na + seq_len(na)
    il <- (k - 1) * nl + seq_len(nl)
    st[ia, ] <- object[, , un[k]]
    K[ia, il] <- gears[[un[k]]]$condALK
  }
  x <- lfd.sim(st, list(condALK = K, ess_age = ess_age, ess_len = ess_len),
               u_age = u_age, u_len = u_len, scale = scale)
  out <- FLQuant(0, dimnames = list(len = lens, year = dimnames(object)$year,
                                    unit = un, iter = dimnames(object)$iter),
                 units = "cm")
  for (k in seq_len(nu)) out[, , un[k]] <- x[(k - 1) * nl + seq_len(nl), ]
  out
}
```

### `sampling.lfd.oem(..., unit.type = c("spatial", "gender"), sexed = TRUE)`

| step | `"spatial"` (current behaviour) | `"gender"` |
|------|---------------------------------|------------|
| `om_gears` | area -> gear | sex -> gear (from `build_gear()` with unit `lhpar`) |
| gear catch | `catch_n_gear()` per unit | `catch_n_gear()` per unit (sex-specific `f_share`) |
| LFD draw | `lfd.sim()` per unit, independent | `lfd.sim.gender()` once per gear |
| `deviances$lfd` | list by area of `rUnif_lfd()` | one `rUnif_lfd()` (per gear) |
| `catch_gear` | by unit | `unitSums()` (catch is not reported by sex) |
| to `est` | by unit, or raised and pooled (`lfd.raising = TRUE`) | by sex if `sexed = TRUE`, else `unitSums()`; no raising |

Default `unit.type`: keep `"spatial"` when `length(unit) > 1` for backward
compatibility, or read `attr(om_gears, "unit.type")` set by `build_gear()`.

## 5. Spatial stocks: `unit.type = "spatial"` (Doc 3)

The current `sampling.lfd.oem()` already implements this case (area-specific
gear sets, independent CRN streams per area, `lfd.raising` with
`lfd.scale = "sample"/"catch"`). Doc 3 will document it with a two-area
example and the checks for raising (with one unit, raising followed by
`lfd.scale = "sample"` returns the sample unchanged), as supporting material
for the ICES deep-sea working group report.
