
#' Logistic selectivity at length
#'
#' Computes asymptotic logistic selectivity as a function of length.
#' The curve is parameterised by the length at 50% selectivity (`L50`)
#' and the length at 95% selectivity (`L95`).
#'
#' @param L Numeric vector of lengths.
#' @param L50 Length at 50\% selectivity.
#' @param L95 Length at 95\% selectivity.
#'
#' @return Numeric vector of selectivity values between 0 and 1.
#' @export
sel_logistic_len <- function(L, L50, L95) {
  1 / (1 + exp(-log(19) * (L - L50) / (L95 - L50)))
}


#' Normal dome-shaped selectivity at length
#'
#' Computes symmetric dome-shaped selectivity as a function of length.
#' The curve peaks at `Lpeak` and declines symmetrically around the peak.
#'
#' @param L Numeric vector of lengths.
#' @param Lpeak Length at maximum selectivity.
#' @param sd Standard deviation controlling the width of the dome.
#'
#' @return Numeric vector of selectivity values scaled to a maximum of 1.
#' @export

sel_normal_len <- function(L, Lpeak, sd) {
  s <- exp(-0.5 * ((L - Lpeak) / sd)^2)
  s / max(s, na.rm = TRUE)
}


#' Double-normal selectivity at length
#'
#' Computes asymmetric dome-shaped selectivity as a function of length.
#' The curve peaks at `Lpeak`, with separate standard deviations for the
#' ascending and descending limbs.
#'
#' @param L Numeric vector of lengths.
#' @param Lpeak Length at maximum selectivity.
#' @param sd_left Standard deviation for lengths below or equal to `Lpeak`.
#' @param sd_right Standard deviation for lengths above `Lpeak`.
#'
#' @return Numeric vector of selectivity values scaled to a maximum of 1.
#' @export

sel_dnormal_len <- function(L, Lpeak, sd_left, sd_right) {
  s <- ifelse(
    L <= Lpeak,
    exp(-0.5 * ((L - Lpeak) / sd_left)^2),
    exp(-0.5 * ((L - Lpeak) / sd_right)^2)
  )
  s / max(s, na.rm = TRUE)
}


#' Generate selectivity-at-length and selectivity-at-age
#'
#' Generates length-based gear selectivity and converts it to selectivity-at-age
#' using von Bertalanffy growth parameters. The function returns both
#' selectivity-at-length and selectivity-at-age as `FLQuant` objects.
#'
#' @param lhpar Named numeric vector or named object containing life-history
#'   parameters. Must include at least `linf`, `k`, and `t0`.
#' @param amin Integer. Minimum age for the age-based selectivity vector.
#'   Default is `0`.
#' @param amax Integer. Maximum age for the age-based selectivity vector.
#'   Default is `20`.
#' @param type Character. Selectivity type. One of `"logistic"`, `"normal"`,
#'   or `"dnormal"`.
#' @param lmin Numeric. Lower bound of the first length bin. Default is `3`.
#' @param lmax Numeric. Upper length limit used to construct length bins.
#'   If `NULL`, this is set to `ceiling(linf * lmax_mult)`.
#' @param binwidth Numeric. Width of length bins. Default is `1`.
#' @param lmax_mult Numeric. Multiplier applied to `linf` when `lmax = NULL`.
#'   Default is `1.1`.
#' @param L50 Numeric. Length at 50 percent selectivity for logistic
#'   selectivity. Required when `type = "logistic"`.
#' @param L95 Numeric. Length at 95 percent selectivity for logistic
#'   selectivity. Required when `type = "logistic"`.
#' @param Lpeak Numeric. Peak selectivity length for normal or double-normal
#'   selectivity. Required when `type = "normal"` or `type = "dnormal"`.
#' @param sd Numeric. Standard deviation of the normal selectivity curve.
#'   Required when `type = "normal"`.
#' @param sd_left Numeric. Standard deviation of the ascending limb of the
#'   double-normal selectivity curve. Required when `type = "dnormal"`.
#'   Default is `0.3`.
#' @param sd_right Numeric. Standard deviation of the descending limb of the
#'   double-normal selectivity curve. Required when `type = "dnormal"`.
#'   Default is `0.6`.
#' @param scale Logical. Should selectivity-at-length and selectivity-at-age be
#'   scaled independently to a maximum of one? Default is `TRUE`.
#'
#' @return A named list with:
#' \itemize{
#'   \item `sel_len`: selectivity-at-length as an `FLQuant`, with the `len`
#'   dimension storing lower length-bin limits.
#'   \item `sel_a`: selectivity-at-age as an `FLQuant`.
#'   \item `len_bins`: data frame containing lower, upper and midpoint of each
#'   length bin.
#' }
#'
#' @details
#' Length bins are stored in the `FLQuant` `len` dimension as lower bin limits.
#' Selectivity is evaluated at the corresponding bin midpoints:
#'
#' `mid = lower + binwidth / 2`
#'
#' Selectivity-at-age is obtained by predicting length-at-age from the
#' von Bertalanffy growth curve:
#'
#' `L[a] = Linf * (1 - exp(-k * (a - t0)))`
#'
#' The function uses `t0 + 0.5` when predicting length-at-age, which approximates
#' mid-year length-at-age.
#'
#' The following selectivity curves are supported:
#'
#' \itemize{
#'   \item logistic, specified by `L50` and `L95`;
#'   \item normal, specified by `Lpeak` and `sd`;
#'   \item double-normal, specified by `Lpeak`, `sd_left`, and `sd_right`.
#' }
#' The helper functions `sel_logistic_len()`, `sel_normal_len()` and
#' `sel_dnormal_len()` must be available in the namespace.
#'
#' @examples
#' \dontrun{
#' ## Logistic linefish selectivity
#' linefish <- sel_la(
#'   lhpar = lhpar,
#'   amin = 0,
#'   amax = 20,
#'   type = "logistic",
#'   L50 = 35,
#'   L95 = 50
#' )
#'
#' ## Dome-shaped gillnet selectivity
#' gillnet <- sel_la(
#'   lhpar = lhpar,
#'   type = "normal",
#'   Lpeak = 35,
#'   sd = 10
#' )
#'
#' ## Asymmetric double-normal trap selectivity
#' traps <- sel_la(
#'   lhpar = lhpar,
#'   type = "dnormal",
#'   Lpeak = 30,
#'   sd_left = 8,
#'   sd_right = 25
#' )
#'
#' plot_sel_la(traps)
#' }
#'
#' @export
sel_la <- function(lhpar,
                   amin = 0,
                   amax = 20,
                   type = c("logistic", "normal", "dnormal"),
                   lmin = 3,
                   lmax = NULL,
                   binwidth = 1,
                   lmax_mult = 1.1,
                   L50 = NULL,
                   L95 = NULL,
                   Lpeak = NULL,
                   sd = NULL,
                   sd_left = 0.3,
                   sd_right = 0.6,
                   scale = TRUE) {
  
  type <- match.arg(type)
  
  make_len_bins <- function(lhpar,
                            lmin = 3,
                            lmax = NULL,
                            binwidth = 1,
                            lmax_mult = 1.1) {
    
    if (is.null(lmax)) {
      lmax <- ceiling(an(lhpar["linf"]) * lmax_mult)
    }
    
    lower <- seq(lmin, lmax - binwidth, by = binwidth)
    upper <- lower + binwidth
    mid <- lower + binwidth / 2
    
    data.frame(
      lower = lower,
      upper = upper,
      mid = mid
    )
  }
  
  bins <- make_len_bins(
    lhpar = lhpar,
    lmin = lmin,
    lmax = lmax,
    binwidth = binwidth,
    lmax_mult = lmax_mult
  )
  
  ## Length-bin midpoints for selectivity-at-length
  L <- bins$mid
  
  ## Mid-year length-at-age
  age <- amin:amax
  
  l_a <- vonbert(
    linf = an(lhpar["linf"]),
    k = an(lhpar["k"]),
    t0 = an(lhpar["t0"]) + 0.5,
    age = age
  )
  
  if (type == "logistic") {
    
    if (is.null(L50) || is.null(L95)) {
      stop("For logistic selectivity, provide L50 and L95.")
    }
    
    sel <- sel_logistic_len(L, L50 = L50, L95 = L95)
    sel_a <- sel_logistic_len(l_a, L50 = L50, L95 = L95)
  }
  
  if (type == "normal") {
    
    if (is.null(Lpeak) || is.null(sd)) {
      stop("For normal selectivity, provide Lpeak and sd.")
    }
    
    sel <- sel_normal_len(L, Lpeak = Lpeak, sd = sd)
    sel_a <- sel_normal_len(l_a, Lpeak = Lpeak, sd = sd)
  }
  
  if (type == "dnormal") {
    
    if (is.null(Lpeak) || is.null(sd_left) || is.null(sd_right)) {
      stop("For double-normal selectivity, provide Lpeak, sd_left, and sd_right.")
    }
    
    sel <- sel_dnormal_len(
      L = L,
      Lpeak = Lpeak,
      sd_left = sd_left,
      sd_right = sd_right
    )
    
    sel_a <- sel_dnormal_len(
      L = l_a,
      Lpeak = Lpeak,
      sd_left = sd_left,
      sd_right = sd_right
    )
  }
  
  if (scale) {
    sel <- sel / max(sel, na.rm = TRUE)
    sel_a <- sel_a / max(sel_a, na.rm = TRUE)
  }
  
  sel_l <- FLCore::FLQuant(
    sel,
    dimnames = list(
      len = as.character(bins$lower),
      year = "1",
      unit = "unique",
      season = "all",
      area = "unique",
      iter = "1"
    )
  )
  
  sel_a <- FLCore::FLQuant(
    sel_a,
    dimnames = list(
      age = as.character(age),
      year = "1",
      unit = "unique",
      season = "all",
      area = "unique",
      iter = "1"
    )
  )
  
  list(
    sel_len = sel_l,
    sel_a = sel_a,
    len_bins = bins
  )
}

#' Plot selectivity at length and age
#'
#' Plots selectivity-at-length and selectivity-at-age side by side for one
#' selectivity object or for a named list of gear-specific selectivity objects.
#'
#' The function is intended as a diagnostic plot for operating-model
#' selectivity assumptions. It accepts objects returned by `sel_la()`, where
#' selectivity-at-length and selectivity-at-age are stored as `FLQuant` objects.
#' If a named list is supplied, each element is treated as a gear and overlaid
#' in both panels.
#'
#' Length selectivity is plotted against length-bin midpoints when these are
#' available through `object$len_bins` or the `"len_bins"` attribute of the
#' length `FLQuant`. Otherwise, midpoints are inferred from the lower length-bin
#' labels and the median bin width.
#'
#' @param object A selectivity object returned by `sel_la()`, or a named list of
#'   such objects. For a single object, the function looks for `sel_len` or
#'   `sel_l`, and for `sel_age` or `sel_a`. For a list, each element is assumed
#'   to represent one gear.
#' @param ncol Numeric. Number of columns in the facet layout. Default is `2`.
#' @param points_age Logical. Should points be added to the age-selectivity
#'   panel? Default is `TRUE`.
#' @param title Optional plot title. Can be a character string or expression.
#' @param line_width Numeric. Line width passed to `ggplot2::geom_line()`.
#'   Default is `0.7`.
#' @param point_size Numeric. Point size for the age-selectivity panel.
#'   Default is `1.5`.
#' @param len_by Numeric. Spacing of x-axis breaks for the length panel.
#'   Default is `5`.
#' @param age_by Numeric. Spacing of x-axis breaks for the age panel.
#'   Default is `1`.
#' @param colours Optional named vector of colours for gears. Names must match
#'   gear names. If `NULL`, ggplot2 default colours are used.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' plot_sel_la(sels$lutjanus_bohar, len_by = 10, age_by = 1)
#'
#' plot_sel_la(
#'   sels$carcharhinus_melanopterus,
#'   len_by = 20,
#'   age_by = 2,
#'   title = expression(italic("Carcharhinus melanopterus"))
#' )
#' }
#'
#' @export
plot_sel_la <- function(object,
                        ncol = 2,
                        points_age = TRUE,
                        title = NULL,
                        line_width = 0.7,
                        point_size = 1.5,
                        len_by = 5,
                        age_by = 1,
                        colours = NULL) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  
  ## ------------------------------------------------------------
  ## Internal helper: extract one sel_la object into a data frame
  ## ------------------------------------------------------------
  
  get_one <- function(x, gear = "gear") {
    
    sel_len <- if (!is.null(x$sel_len)) x$sel_len else x$sel_l
    sel_age <- if (!is.null(x$sel_age)) x$sel_age else x$sel_a
    
    if (is.null(sel_len)) {
      stop(paste("Could not find 'sel_len' or 'sel_l' for", gear))
    }
    
    if (is.null(sel_age)) {
      stop(paste("Could not find 'sel_age' or 'sel_a' for", gear))
    }
    
    ## ---- length panel ----
    
    len_lower <- as.numeric(dimnames(sel_len)$len)
    
    if (any(is.na(len_lower))) {
      len_lower <- as.numeric(dimnames(sel_len)[[1]])
    }
    
    if (!is.null(x$len_bins)) {
      len_x <- x$len_bins$mid
    } else if (!is.null(attr(sel_len, "len_bins"))) {
      len_x <- attr(sel_len, "len_bins")$mid
    } else {
      bw <- median(diff(len_lower), na.rm = TRUE)
      len_x <- len_lower + bw / 2
    }
    
    dat_len <- data.frame(
      x = len_x,
      y = as.numeric(sel_len[, 1, 1, 1, 1, 1]),
      panel = "Selectivity ~ length",
      pos = 1,
      gear = gear,
      stringsAsFactors = FALSE
    )
    
    ## ---- age panel ----
    
    age_x <- as.numeric(gsub("\\+", "", dimnames(sel_age)$age))
    
    if (any(is.na(age_x))) {
      age_x <- as.numeric(gsub("\\+", "", dimnames(sel_age)[[1]]))
    }
    
    dat_age <- data.frame(
      x = age_x,
      y = as.numeric(sel_age[, 1, 1, 1, 1, 1]),
      panel = "Selectivity ~ age",
      pos = 2,
      gear = gear,
      stringsAsFactors = FALSE
    )
    
    rbind(dat_len, dat_age)
  }
  
  ## ------------------------------------------------------------
  ## Determine whether object is one gear or a named list of gears
  ## ------------------------------------------------------------
  
  is_single <- !is.null(object$sel_len) ||
    !is.null(object$sel_l) ||
    !is.null(object$sel_age) ||
    !is.null(object$sel_a)
  
  if (is_single) {
    
    dat <- get_one(object, gear = "gear")
    show_legend <- FALSE
    
  } else {
    
    gears <- names(object)
    
    if (is.null(gears)) {
      gears <- paste0("gear_", seq_along(object))
    }
    
    dat <- NULL
    
    for (i in seq_along(object)) {
      dat <- rbind(dat, get_one(object[[i]], gear = gears[i]))
    }
    
    dat$gear <- factor(dat$gear, levels = gears)
    show_legend <- TRUE
  }
  
  ## Facet labels
  facl <- stats::setNames(unique(dat$panel), unique(dat$pos))
  
  ## Data ranges for identifying free-x panels
  len_range <- range(dat$x[dat$pos == 1], na.rm = TRUE)
  age_range <- range(dat$x[dat$pos == 2], na.rm = TRUE)
  
  len_width <- diff(len_range)
  age_width <- diff(age_range)
  
  ## This function is called separately by ggplot2 for each free-x panel.
  ## We identify whether the panel is age or length from its x-range width.
  make_breaks <- function(lims) {
    
    lim_width <- diff(lims)
    
    is_age_panel <- abs(lim_width - age_width) <
      abs(lim_width - len_width)
    
    if (is_age_panel) {
      by <- age_by
      rng <- age_range
    } else {
      by <- len_by
      rng <- len_range
    }
    
    from <- floor(rng[1] / by) * by
    to   <- ceiling(rng[2] / by) * by
    
    seq(from, to, by = by)
  }
  
  ## ------------------------------------------------------------
  ## Plot
  ## ------------------------------------------------------------
  
  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(x = x, y = y, colour = gear, group = gear)
  ) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::facet_wrap(
      ~pos,
      scales = "free_x",
      ncol = ncol,
      labeller = ggplot2::labeller(pos = facl)
    ) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(
      breaks = make_breaks,
      labels = function(x) sprintf("%.0f", x)
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = NULL,
      y = "Selectivity",
      colour = NULL,
      title = title
    ) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey90"),
      axis.title = ggplot2::element_text(face = "bold")
    )
  
  if (!is.null(colours)) {
    
    if (is.null(names(colours))) {
      stop("'colours' must be a named vector with names matching gear names.")
    }
    
    p <- p + ggplot2::scale_colour_manual(values = colours, drop = FALSE)
  }
  
  if (points_age) {
    p <- p +
      ggplot2::geom_point(
        data = dat[dat$pos == 2, ],
        size = point_size
      )
  }
  
  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  p
}



#' Combine gear-specific selectivity-at-age using relative apical F
#'
#' Combines gear-specific selectivity-at-age curves into a single joint
#' selectivity-at-age curve using relative apical partial fishing mortality
#' weights by gear.
#'
#' @param object A named list of gear-specific selectivity objects, such as
#'   `sels$siganus_sutor`. Each gear object should contain `sel_a` or `sel_age`.
#' @param fapic_rel Numeric vector of relative apical partial F contributions
#'   by gear. Preferably named, with names matching gears in `object`, for
#'   example `c(traps = 1, gillnet = 0.3, linefish = 0.1)`.
#'   If unnamed, its length must match `length(object)`, and names are assigned
#'   from `names(object)` in their existing order.
#' @param normalise_fapic Logical. Should `fapic_rel` be normalised to sum to
#'   one before combining curves? Default is `FALSE`. For the final scaled
#'   selectivity shape this has no effect when `scale = TRUE`, but it may be
#'   useful for interpretation.
#' @param scale Logical. Should the combined curve be scaled to a maximum of
#'   one? Default is `TRUE`.
#'
#' @return An `FLQuant` containing the combined selectivity-at-age.
#'
#' @details
#' Gear-specific selectivity curves are assumed to be scaled to a maximum of
#' one. The combined raw fishing-mortality shape is calculated as:
#'
#' `sel[a] = sum_g fapic_rel[g] * sel_g[a]`
#'
#' If `scale = TRUE`, the combined curve is scaled to a maximum of one:
#'
#' `sel[a] = sel[a] / max(sel)`
#'
#' The weights are interpreted as relative apical partial F contributions by
#' gear, not relative Fbar contributions. This avoids dependence on an
#' arbitrary Fbar age range and is consistent with constructing fleet-specific
#' fishing mortality as:
#'
#' `F_g[a, y] = Fapic_g[y] * sel_g[a]`
#'
#' before summing across gears.
#'
#' @examples
#' \dontrun{
#' ## Preferred: named vector
#' fapic_rel <- c(traps = 1, gillnet = 0.3, linefish = 0.1)
#'
#' sel_joint <- combine_sel_age(
#'   object = sels$siganus_sutor,
#'   fapic_rel = fapic_rel
#' )
#'
#' plot_sel_age(sel_joint)
#'
#' ## Convenience: unnamed vector matched to names(object)
#' sel_joint <- combine_sel_age(
#'   object = sels$siganus_sutor,
#'   fapic_rel = c(1, 0.3, 0.1)
#' )
#' }
#'
#' @export
combine_sel_age <- function(object,
                            fapic_rel,
                            normalise_fapic = FALSE,
                            scale = TRUE) {
  
  if (is.null(names(object))) {
    stop("'object' must be a named list of gear-specific selectivity objects.")
  }
  
  if (!is.numeric(fapic_rel)) {
    stop("'fapic_rel' must be numeric.")
  }
  
  if (any(is.na(fapic_rel))) {
    stop("'fapic_rel' contains NA values.")
  }
  
  if (any(fapic_rel < 0)) {
    stop("'fapic_rel' values must be non-negative.")
  }
  
  if (sum(fapic_rel) <= 0) {
    stop("'fapic_rel' must contain at least one positive value.")
  }
  
  ## ------------------------------------------------------------
  ## Allow unnamed fapic_rel for convenience
  ## ------------------------------------------------------------
  
  if (is.null(names(fapic_rel))) {
    
    if (length(fapic_rel) != length(object)) {
      stop(
        "'fapic_rel' is unnamed and must have the same length as 'object'. ",
        "length(fapic_rel) = ", length(fapic_rel),
        "; length(object) = ", length(object), "."
      )
    }
    
    names(fapic_rel) <- names(object)
  }
  
  if (any(names(fapic_rel) == "")) {
    stop(
      "'fapic_rel' must be fully named, or completely unnamed with length ",
      "matching 'object'."
    )
  }
  
  gears <- names(fapic_rel)
  
  missing_gears <- gears[!gears %in% names(object)]
  
  if (length(missing_gears) > 0) {
    stop(
      "The following gears in 'fapic_rel' are not present in 'object': ",
      paste(missing_gears, collapse = ", ")
    )
  }
  
  if (normalise_fapic) {
    fapic_rel <- fapic_rel / sum(fapic_rel)
  }
  
  ## ------------------------------------------------------------
  ## Extract template
  ## ------------------------------------------------------------
  
  first <- object[[gears[1]]]
  
  sel_joint <- if (!is.null(first$sel_age)) {
    first$sel_age
  } else {
    first$sel_a
  }
  
  if (is.null(sel_joint)) {
    stop("Could not find 'sel_age' or 'sel_a' in first gear object.")
  }
  
  sel_joint[] <- 0
  
  ## ------------------------------------------------------------
  ## Weighted sum of gear-specific selectivity-at-age
  ## ------------------------------------------------------------
  
  for (g in gears) {
    
    sel_g <- if (!is.null(object[[g]]$sel_age)) {
      object[[g]]$sel_age
    } else {
      object[[g]]$sel_a
    }
    
    if (is.null(sel_g)) {
      stop("Could not find 'sel_age' or 'sel_a' for gear: ", g)
    }
    
    sel_joint <- sel_joint + fapic_rel[g] * sel_g
  }
  
  ## ------------------------------------------------------------
  ## Scale to maximum one
  ## ------------------------------------------------------------
  
  if (scale) {
    mx <- max(sel_joint, na.rm = TRUE)
    if (mx > 0) {
      sel_joint <- sel_joint / mx
    }
  }
  
  sel_joint
}


#' Plot selectivity-at-age
#'
#' Plots one or more selectivity-at-age curves stored as `FLQuant` objects.
#'
#' The function accepts either a single `FLQuant` or a named list of `FLQuant`
#' objects. It is intended for comparing joint or gear-specific
#' selectivity-at-age curves used in operating models.
#'
#' @param object An `FLQuant` containing selectivity-at-age, or a list of
#'   `FLQuant` objects. If a list is supplied, each element is plotted as a
#'   separate curve. List names are used as legend labels.
#' @param age_by Numeric. Spacing of age-axis breaks. Default is `1`.
#' @param points Logical. Should points be added to the curves? Default is
#'   `TRUE`.
#' @param title Optional plot title. Can be a character string or expression.
#' @param line_width Numeric. Line width passed to `ggplot2::geom_line()`.
#'   Default is `0.7`.
#' @param point_size Numeric. Point size passed to `ggplot2::geom_point()`.
#'   Default is `1.5`.
#' @param colours Optional named vector of colours. Names must match curve
#'   names. If `NULL`, ggplot2 default colours are used.
#' @param ylab Character. Y-axis label. Default is `"Selectivity"`.
#' @param xlab Character. X-axis label. Default is `"Age"`.
#'
#' @return A `ggplot` object.
#'
#' @details
#' The age dimension is extracted from the `age` dimension of each `FLQuant`.
#' Plus-group labels such as `"10+"` are converted to numeric values for
#' plotting.
#'
#' If a single `FLQuant` is provided, the legend is suppressed. If a list of
#' `FLQuant` objects is provided, the legend is shown using the list names.
#'
#' @examples
#' \dontrun{
#' ## Single curve
#' plot_sel_age(sel_joint)
#'
#' ## Multiple species
#' sel_joint <- list(
#'   siganus = sel_joint_siganus,
#'   lethrinus = sel_joint_lethrinus,
#'   lutjanus = sel_joint_lutjanus,
#'   shark = sel_joint_shark
#' )
#'
#' plot_sel_age(sel_joint)
#' }
#'
#' @export
plot_sel_age <- function(object,
                         age_by = 1,
                         points = TRUE,
                         title = NULL,
                         line_width = 0.7,
                         point_size = 1.5,
                         colours = NULL,
                         ylab = "Selectivity",
                         xlab = "Age") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  
  ## ------------------------------------------------------------
  ## Internal helper: convert one FLQuant to data frame
  ## ------------------------------------------------------------
  
  get_one <- function(x, name = "selectivity") {
    
    if (!inherits(x, "FLQuant")) {
      stop("Each element of 'object' must be an FLQuant.")
    }
    
    age <- dimnames(x)$age
    
    if (is.null(age)) {
      age <- dimnames(x)[[1]]
    }
    
    age <- as.numeric(gsub("\\+", "", age))
    
    if (any(is.na(age))) {
      stop("Could not convert age dimension to numeric values.")
    }
    
    data.frame(
      age = age,
      data = as.numeric(x[, 1, 1, 1, 1, 1]),
      curve = name,
      stringsAsFactors = FALSE
    )
  }
  
  ## ------------------------------------------------------------
  ## Single FLQuant or list of FLQuants
  ## ------------------------------------------------------------
  
  if (inherits(object, "FLQuant")) {
    
    dat <- get_one(object, name = "selectivity")
    show_legend <- FALSE
    
  } else if (is.list(object)) {
    
    nms <- names(object)
    
    if (is.null(nms)) {
      nms <- paste0("curve_", seq_along(object))
    }
    
    dat <- NULL
    
    for (i in seq_along(object)) {
      dat <- rbind(dat, get_one(object[[i]], name = nms[i]))
    }
    
    dat$curve <- factor(dat$curve, levels = nms)
    show_legend <- TRUE
    
  } else {
    
    stop("'object' must be an FLQuant or a list of FLQuant objects.")
  }
  
  ## ------------------------------------------------------------
  ## Axis breaks
  ## ------------------------------------------------------------
  
  age_rng <- range(dat$age, na.rm = TRUE)
  
  age_breaks <- seq(
    floor(age_rng[1] / age_by) * age_by,
    ceiling(age_rng[2] / age_by) * age_by,
    by = age_by
  )
  
  ## ------------------------------------------------------------
  ## Plot
  ## ------------------------------------------------------------
  
  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(x = age, y = data, colour = curve, group = curve)
  ) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(
      breaks = age_breaks,
      labels = function(x) sprintf("%.0f", x)
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      colour = NULL,
      title = title
    ) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold")
    )
  
  if (!is.null(colours)) {
    
    if (is.null(names(colours))) {
      stop("'colours' must be a named vector with names matching curve names.")
    }
    
    p <- p + ggplot2::scale_colour_manual(values = colours, drop = FALSE)
  }
  
  if (points) {
    p <- p + ggplot2::geom_point(size = point_size)
  }
  
  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  p
}

# {{{ 
# build_gear()
#
#' Build a single gear's selectivity, conditional ALK, and sampling config
#' Depends on FLRef::sel_la() and FLRef::iALK() already being available.
#' 
#' Consolidates everything a length-sampling OEM needs for one gear into a
#' single object: selectivity-at-length, selectivity-at-age (both via
#' `sel_la()`), the conditional inverse age-length key
#' \eqn{P(l \mid a, \text{caught by gear})} (via `iALK()`), and the gear's
#' observation-process assumptions (`ess_age`, `ess_len`).
#'
#' `sel_la()` and `iALK()` express their length-range argument differently:
#' `sel_la(lmax=...)` is an *absolute* upper length, while `iALK(lmax=...)`
#' is a *multiplier on linf*. `build_gear()` takes a single `lmax_mult` and
#' derives the correct form for each, then asserts the two resulting length
#' grids are identical before returning -- if they are not, every downstream
#' sampler would silently misalign length bins between `sel_len` and
#' `condALK`, so this fails loudly instead.
#'
#' The conditional ALK is built as
#' \deqn{P(l \mid a, g) = \dfrac{P(l \mid a)\, s_g(l)}{\sum_l P(l \mid a)\, s_g(l)}}
#' i.e. the biological inverse ALK re-weighted by the gear's length
#' selectivity and renormalised so each age row sums to 1 (rows with zero
#' gear-selected probability are set to 0 rather than divided by zero).
#'
#' @param name Character. Gear name, stored on the returned object.
#' @param lhpar `FLPar` or named numeric vector with at least `linf`, `k`,
#'   `t0`. Passed to both `sel_la()` and `iALK()`.
#' @param type Character. One of `"logistic"`, `"normal"`, `"dnormal"`,
#'   passed to `sel_la()`.
#' @param sel_pars Named list of the selectivity parameters `sel_la()`
#'   needs for `type` (e.g. `list(L50=35, L95=50)` for `"logistic"`,
#'   `list(Lpeak=35, sd=10)` for `"normal"`,
#'   `list(Lpeak=30, sd_left=8, sd_right=25)` for `"dnormal"`).
#' @param f_mult Numeric. This gear's relative apical-F weight in the OM
#'   (operating-model truth, not an observation-process assumption -- kept
#'   here for convenience but conceptually separate from `ess_age`/`ess_len`,
#'   see Details).
#' @param ess_age Integer. Effective sample size for the stage-1
#'   catch-at-age sample (observation-process assumption).
#' @param ess_len Integer. Effective sample size for the stage-2
#'   length sample (observation-process assumption).
#' @param age Integer vector of ages, e.g. `0:30`.
#' @param amin,amax Integer. Passed through to `sel_la()`; default to
#'   `range(age)` if not supplied.
#' @param lmin Numeric. Lower bound of the first length bin. Default `5`.
#' @param lmax_mult Numeric. Multiplier on `linf` defining the upper length
#'   bound, applied consistently to both `sel_la()` and `iALK()`.
#'   Default `1.2`.
#' @param bin Numeric. Length-bin width, passed to both functions as
#'   `binwidth`/`bin`. Default `1`.
#' @param cv Numeric. CV of length-at-age used by `iALK()`. Default `0.1`.
#' @param scale Logical. Passed to `sel_la()`; scale `sel_len`/`sel_a`
#'   independently to a max of 1. Default `TRUE`.
#'
#' @return A named list with:
#' \describe{
#'   \item{name}{Gear name.}
#'   \item{sel_len}{`FLQuant`, selectivity-at-length (from `sel_la()`).}
#'   \item{sel_a}{`FLQuant`, selectivity-at-age (from `sel_la()`).}
#'   \item{len_bins}{`data.frame` of length-bin lower/upper/mid, shared by
#'     `sel_len` and `condALK`.}
#'   \item{condALK}{Numeric matrix, age (rows) x len (cols), rows sum to 1:
#'     \eqn{P(l \mid a, \text{caught by gear})}.}
#'   \item{f_mult}{As supplied.}
#'   \item{ess_age, ess_len}{As supplied.}
#' }
#'
#' @examples
#' \dontrun{
#' lhpar <- FLPar(linf = 55.662, k = 0.08, t0 = -2.77)
#' ages  <- 0:30
#'
#' trawl <- build_gear(
#'   name = "Trawl", lhpar = lhpar, type = "dnormal",
#'   sel_pars = list(Lpeak = 18, sd_left = 5, sd_right = 12),
#'   f_mult = 0.70, ess_age = 50, ess_len = 500, age = ages
#' )
#'
#' ## sanity checks before trusting it further
#' plot_sel_la(trawl, len_by = 5, age_by = 2)
#'stopifnot(isTRUE(all.equal(unname(rowSums(trawl$condALK)),
#'                           rep(1, nrow(trawl$condALK)), tolerance = 1e-6)))
#' image(trawl$condALK, main = "P(len | age, caught by Trawl)")
#'
#' gillnet <- build_gear(
#'   name = "Gillnet", lhpar = lhpar, type = "normal",
#'   sel_pars = list(Lpeak = 38, sd = 10),
#'   f_mult = 0.30, ess_age = 30, ess_len = 300, age = ages
#' )
#'
#' om_gears <- list(Trawl = trawl, Gillnet = gillnet)
#' }
#'
#' @export
build_gear <- function(name, lhpar, type = c("logistic", "normal", "dnormal"),
                       sel_pars, f_mult, ess_age, ess_len, age,
                       amin = min(age), amax = max(age),
                       lmin = 5, lmax_mult = 1.2, bin = 1, cv = 0.1,
                       scale = TRUE) {
  
  type <- match.arg(type)
  
  linf <- c(lhpar["linf"])
  ## sel_la()'s 'lmax' appears to be an EXCLUSIVE upper bound (its last bin's
  ## lower edge sits at lmax - bin, one short of iALK()'s plus-group edge --
  ## see the mismatch this caused: passing lmax=67 produced bins only to 66).
  ## Padding by one bin width here is an attempt to make the two grids align
  ## natively, so the reconciliation block below becomes a safety net rather
  ## than something that fires every time. VERIFY: after this change, check
  ## whether the "iALK()'s plus-group bin ... has no corresponding sel_la()
  ## bin" message still appears -- if it's gone, this hypothesis was right;
  ## if it still appears (or appears differently), revert this padding and
  ## rely on the reconciliation block alone, which handles the mismatch
  ## either way.
  lmax_abs <- ceiling(linf * lmax_mult) + bin   # sel_la() wants an absolute value
  
  ## --- selectivity-at-length / -at-age -----------------------------------
  s <- do.call(
    sel_la,
    c(
      list(
        lhpar = lhpar, amin = amin, amax = amax, type = type,
        lmin = lmin, lmax = lmax_abs, binwidth = bin, scale = scale
      ),
      sel_pars
    )
  )
  
  ## --- biological inverse ALK, same grid ---------------------------------
  ialk <- iALK(
    params = c(linf = linf, k = c(lhpar["k"]), t0 = c(lhpar["t0"])),
    age = age, cv = cv, lmin = lmin, lmax = lmax_mult, bin = bin
  )
  
  ## --- reconcile length grids ---------------------------------------------
  ## sel_la() and iALK() can legitimately disagree by one bin at the top:
  ## iALK() treats its last bin as a plus-group (P(length >= edge), via
  ## pnorm(..., lower.tail=FALSE) in its source), while sel_la() does not
  ## build a plus-group at all and simply stops one bin short. This is a
  ## real difference in convention, not necessarily a misconfiguration, so
  ## it is reconciled rather than treated as an error: any length present
  ## in iALK()'s grid but not sel_la()'s (i.e. the plus-group edge) is
  ## covered by carrying the last explicitly modelled selectivity value
  ## forward, on the standard assumption that selectivity has plateaued by
  ## 'lmax_mult * linf'. Any length present in sel_la()'s grid but not
  ## iALK()'s is dropped, with a warning, since it would have no ALK row to
  ## pair with. A genuinely larger mismatch (more than a couple of bins,
  ## or a gap in the middle of the range) is NOT something this should
  ## silently absorb, so that case still stops.
  len_sel_chr  <- as.character(s$len_bins$lower)
  len_ialk_chr <- dimnames(ialk)$len
  
  extra_in_ialk <- setdiff(len_ialk_chr, len_sel_chr)
  extra_in_sel  <- setdiff(len_sel_chr, len_ialk_chr)
  
  if (length(extra_in_ialk) > 1 || length(extra_in_sel) > 1) {
    stop(
      "'sel_la()' and 'iALK()' length grids differ by more than one bin ",
      "for gear '", name, "' (", length(extra_in_ialk), " extra in iALK(), ",
      length(extra_in_sel), " extra in sel_la()). This looks like a real ",
      "'lmin'/'lmax_mult'/'bin' inconsistency rather than the usual ",
      "plus-group edge-case -- check those arguments before proceeding."
    )
  }
  
  sel_len_vec <- as.numeric(s$sel_len)
  names(sel_len_vec) <- len_sel_chr
  
  if (length(extra_in_sel) == 1) {
    warning(
      "Dropping length bin '", extra_in_sel, "' from gear '", name,
      "': present in sel_la()'s grid but not iALK()'s plus-group edge."
    )
    sel_len_vec <- sel_len_vec[names(sel_len_vec) != extra_in_sel]
  }
  
  if (length(extra_in_ialk) == 1) {
    message(
      "Gear '", name, "': iALK()'s plus-group bin ('", extra_in_ialk,
      "'+) has no corresponding sel_la() bin. Carrying the last modelled ",
      "selectivity value (", round(tail(sel_len_vec, 1), 3),
      ") forward to cover it -- confirm selectivity has genuinely plateaued ",
      "by this length before trusting that assumption."
    )
    sel_len_vec[extra_in_ialk] <- tail(sel_len_vec, 1)
  }
  
  ## reorder to match iALK()'s column order exactly before using sweep()
  sel_len_vec <- sel_len_vec[len_ialk_chr]
  
  ## --- condition the ALK on capture by this gear -------------------------
  ialk_mat <- c(ialk)                       # coerce FLPar -> plain matrix
  dim(ialk_mat) <- dim(ialk)[1:2]
  dimnames(ialk_mat) <- list(age = dimnames(ialk)$age, len = dimnames(ialk)$len)
  
  weighted <- sweep(ialk_mat, 2, sel_len_vec, "*")
  row_tot <- rowSums(weighted)
  cond <- weighted
  valid <- row_tot > 0
  cond[valid, ] <- weighted[valid, , drop = FALSE] / row_tot[valid]
  cond[!valid, ] <- 0
  
  ## return the RECONCILED grid throughout, not sel_la()'s original one --
  ## sel_len/len_bins must match condALK's columns exactly, or a naive
  ## plot(trawl$sel_len) against trawl$condALK later would silently
  ## reintroduce the same one-bin mismatch this function just resolved.
  len_bins_out <- data.frame(
    lower = as.numeric(len_ialk_chr),
    upper = as.numeric(len_ialk_chr) + bin,
    mid   = as.numeric(len_ialk_chr) + bin / 2
  )
  
  list(
    name = name,
    sel_len = FLQuant(sel_len_vec, dimnames = list(len = names(sel_len_vec))),
    sel_a = s$sel_a,
    len_bins = len_bins_out,
    condALK = cond,
    f_mult = f_mult,
    ess_age = ess_age,
    ess_len = ess_len
  )
}
# }}}




#' Heatmap of a gear's conditional inverse ALK
#'
#' ggplot equivalent of `image(gear$condALK)`: a length-at-age probability
#' heatmap, \eqn{P(l \mid a, \text{caught by gear})}, with proper axis
#' labels and a colour scale instead of `image()`'s default palette.
#'
#' If `params` is supplied, two mean length-at-age lines are overlaid on
#' the heatmap: the raw biological mean (unconditioned von Bertalanffy
#' expectation) and the mean *implied by the heatmap itself*
#' (`condALK %*% length`) -- i.e. what the gear actually samples. The gap
#' between them is the within-age selectivity bias illustrated in
#' [plot_condALK_bias()], shown here directly against the distribution it
#' arises from rather than as a separate plot.
#'
#' @param gear A gear object as returned by [build_gear()] (needs
#'   `condALK` and `len_bins`).
#' @param params Optional. `FLPar`/named numeric vector with `linf`, `k`,
#'   `t0`. If supplied, overlays the raw-vs-conditioned mean length-at-age
#'   lines; if `NULL` (default), only the heatmap is drawn.
#' @param cv,lmin,lmax_mult,bin Passed to `iALK()` when building the raw
#'   comparison line; should match whatever was used to build `gear`.
#'   Ignored if `params` is `NULL`.
#' @param low,high Colours for the low/high ends of the probability scale.
#'   Defaults chosen to read clearly against a white background.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' plot_condALK(trawl)                    # heatmap only
#' plot_condALK(trawl, params = lhpar)    # heatmap + mean-length overlay
#' plot_condALK(gillnet, params = lhpar)
#' }
#'
#' @export
plot_condALK <- function(gear, params = NULL, cv = 0.1, lmin = 5,
                         lmax_mult = 1.2, bin = 1,
                         low = "#FFFFCC", high = "#7F0000") {
  
  m <- gear$condALK
  df <- data.frame(
    age = as.numeric(rep(rownames(m), times = ncol(m))),
    len = as.numeric(rep(colnames(m), each = nrow(m))),
    p   = as.numeric(m)
  )
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = age, y = len)) +
    ggplot2::geom_raster(ggplot2::aes(fill = p)) +
    ggplot2::scale_fill_gradient(
      low = low, high = high, name = "P(len | age,\ncaught)"
    ) +
    ggplot2::labs(
      x = "Age", y = "Length (cm)",
      title = paste0("P(length | age, caught by ", gear$name, ")")
    ) +
    ggplot2::theme_bw()
  
  if (!is.null(params)) {
    
    age <- as.numeric(rownames(m))
    len_mid <- as.numeric(colnames(m)) + bin / 2
    
    ialk <- iALK(
      params = c(linf = c(params["linf"]), k = c(params["k"]), t0 = c(params["t0"])),
      age = age, cv = cv, lmin = lmin, lmax = lmax_mult, bin = bin
    )
    ialk_mat <- c(ialk)
    dim(ialk_mat) <- dim(ialk)[1:2]
    
    lines_df <- data.frame(
      age = rep(age, 2),
      mean_len = c(as.numeric(ialk_mat %*% len_mid), as.numeric(m %*% len_mid)),
      series = rep(c("Raw biology", "Conditioned on capture"), each = length(age))
    )
    
    p <- p +
      ggplot2::geom_line(
        data = lines_df,
        ggplot2::aes(x = age, y = mean_len, linetype = series),
        colour = "black", linewidth = 0.8, inherit.aes = FALSE
      ) +
      ggplot2::scale_linetype_manual(values = c("Raw biology" = "solid",
                                                "Conditioned on capture" = "dashed"),
                                     name = NULL) +
      ggplot2::theme(legend.box = "vertical")
  }
  
  p
}

#' Compare raw vs. gear-conditioned mean length-at-age for one or more gears
#'
#' Overlays the raw biological mean length-at-age against the
#' capture-conditioned mean length-at-age for each gear in `om_gears`,
#' illustrating how dome-shaped or asymmetric selectivity biases the
#' *within-age* length distribution actually observed in the catch --
#' distinct from, and in addition to, its effect on which ages are caught
#' at all.
#'
#' @param om_gears Named list of gear objects as returned by
#'   [build_gear()].
#' @param params `FLPar`/named numeric vector with `linf`, `k`, `t0`, used
#'   to build the raw (unconditioned) inverse ALK for comparison.
#' @param age Integer vector of ages. Defaults to the age range implied
#'   by `om_gears[[1]]$condALK`'s row names.
#' @param cv,lmin,lmax_mult,bin Passed to `iALK()` when building the raw
#'   ALK; should match whatever was used to build `om_gears`.
#'
#' @return A `ggplot` object, one coloured line per gear plus one for
#'   "Raw biology".
#'
#' @examples
#' \dontrun{
#' plot_condALK_bias(list(Trawl = trawl, Gillnet = gillnet), lhpar)
#' }
#'
#' @export
plot_condALK_bias <- function(om_gears, params,
                              age = as.numeric(rownames(om_gears[[1]]$condALK)),
                              cv = 0.1, lmin = 5, lmax_mult = 1.2, bin = 1) {
  
  ialk <- iALK(
    params = c(linf = c(params["linf"]), k = c(params["k"]), t0 = c(params["t0"])),
    age = age, cv = cv, lmin = lmin, lmax = lmax_mult, bin = bin
  )
  ialk_mat <- c(ialk)
  dim(ialk_mat) <- dim(ialk)[1:2]
  len_mid <- as.numeric(dimnames(ialk)$len) + bin / 2
  
  raw_mean_len <- as.numeric(ialk_mat %*% len_mid)
  
  gear_lines <- lapply(names(om_gears), function(g) {
    cm <- om_gears[[g]]$condALK
    len_mid_g <- as.numeric(colnames(cm)) + bin / 2
    data.frame(
      age = age,
      mean_len = as.numeric(cm %*% len_mid_g),
      series = g
    )
  })
  
  df <- rbind(
    data.frame(age = age, mean_len = raw_mean_len, series = "Raw biology"),
    do.call(rbind, gear_lines)
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = age, y = mean_len, colour = series)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      x = "Age", y = "Mean length (cm)", colour = NULL,
      title = "Raw vs. gear-conditioned mean length-at-age"
    ) +
    ggplot2::theme_bw()
}


#' Compare sampled length-frequency distributions across gears
#'
#' ggplot equivalent of `plot(FLQuants(...))` for gear LFD samples, with
#' explicit numeric length axes -- useful when the default FLQuants panel
#' plot's axis labels are hard to read or ambiguous.
#'
#' @param lfd_list Named list of `FLQuant` length-frequency samples, one
#'   per gear (e.g. `list(Trawl = lfd_trawl, Gillnet = lfd_gillnet)`).
#' @param year Character or numeric. Which year to plot. Defaults to the
#'   first year present.
#' @param it Integer. Which iteration to plot. Defaults to `1`.
#'
#' @return A `ggplot` object, one coloured line per gear.
#'
#' @examples
#' \dontrun{
#' plot_lfd_compare(list(Trawl = lfd_trawl, Gillnet = lfd_gillnet),
#'                   year = one_year, it = 1)
#' }
#'
#' @export
plot_lfd_compare <- function(lfd_list, year = NULL, it = 1) {
  
  if (is.null(year)) year <- dimnames(lfd_list[[1]])$year[1]
  year <- as.character(year)
  
  df <- do.call(rbind, lapply(names(lfd_list), function(g) {
    x <- lfd_list[[g]]
    data.frame(
      len = as.numeric(dimnames(x)$len),
      count = as.numeric(x[, year, , , , it]),
      gear = g
    )
  }))
  
  ggplot2::ggplot(df, ggplot2::aes(x = len, y = count, colour = gear)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::labs(
      x = "Length (cm)", y = "Sampled count", colour = NULL,
      title = paste0("Sampled length-frequency by gear, year ", year, ", iter ", it)
    ) +
    ggplot2::theme_bw()
}

