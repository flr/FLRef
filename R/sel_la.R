
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


#' Construct selectivity-at-length and selectivity-at-age
#'
#' Evaluates a logistic, normal, or double-normal selectivity curve on a
#' length grid and maps the same curve to predicted length-at-age.
#'
#' @param lhpar Life-history parameters containing `linf`, `k`, and `t0`.
#' @param amin,amax Integer minimum and maximum ages. Defaults to `0` and `20`.
#' @param type Selectivity curve: `"logistic"`, `"normal"`, or `"dnormal"`.
#' @param lmin Numeric lower limit of the length grid. Defaults to `3`.
#' @param lmax Optional numeric upper limit of the length grid. When `NULL`,
#'   it is calculated from `linf * lmax_mult`.
#' @param binwidth Numeric length-bin width. Defaults to `1`.
#' @param lmax_mult Numeric multiplier applied to `linf` when `lmax` is
#'   `NULL`. Defaults to `1.1`.
#' @param L50,L95 Lengths at 50 and 95 percent selectivity for a logistic
#'   curve.
#' @param Lpeak Length at maximum selectivity for normal and double-normal
#'   curves.
#' @param sd Standard deviation of a normal curve.
#' @param sd_left,sd_right Standard deviations of the ascending and descending
#'   limbs of a double-normal curve.
#' @param scale Logical; if `TRUE`, scale the length and age curves separately
#'   to a maximum of one.
#'
#' @return A list with components `sel_len`, an `FLQuant` of
#'   selectivity-at-length; `sel_a`, an `FLQuant` of selectivity-at-age; and
#'   `len_bins`, a data frame of lower, upper, and midpoint lengths.
#'
#' @details
#' Length-dimension labels contain lower bin limits, while selectivity is
#' evaluated at bin midpoints. The returned age curve is a point evaluation at
#' predicted length-at-age. Use [build_gear()] when selectivity-at-age must be
#' integrated over a length-at-age distribution.
#'
#' @examples
#' \dontrun{
#' linefish <- sel_la(
#'   lhpar = lhpar,
#'   amin = 0,
#'   amax = 20,
#'   type = "logistic",
#'   L50 = 35,
#'   L95 = 50
#' )
#' gillnet <- sel_la(
#'   lhpar = lhpar,
#'   type = "normal",
#'   Lpeak = 35,
#'   sd = 10
#' )
#' traps <- sel_la(
#'   lhpar = lhpar,
#'   type = "dnormal",
#'   Lpeak = 30,
#'   sd_left = 8,
#'   sd_right = 25
#' )
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
#' Build gear selectivity and sampling configurations
#'
#' Constructs selectivity-at-length, integrated selectivity-at-age, a
#' capture-conditioned inverse age-length key, and observation sample sizes
#' for each gear in a multigear operating model.
#'
#' @param gear_cfg Named list of gear configurations. Each gear must contain
#'   `type`, `f_mult`, `ess_age`, and `ess_len`; remaining entries are passed
#'   to [sel_la()] as selectivity parameters.
#' @param lhpar Life-history parameters containing `linf`, `k`, and `t0`.
#' @param age Numeric vector of ages.
#' @param amin,amax Integer minimum and maximum ages passed to [sel_la()].
#' @param lmin Numeric lower limit of the length grid. Defaults to `5`.
#' @param lmax_mult Numeric multiplier applied to `linf` for the upper length
#'   limit. Defaults to `1.2`.
#' @param bin Numeric length-bin width. Defaults to `1`.
#' @param cv Numeric coefficient of variation in length-at-age. Defaults to
#'   `0.1`.
#' @param timing Numeric within-year timing expressed as a fraction of a year.
#'   Growth at time `t` is evaluated using `t0 - t`. Defaults to `0`.
#' @param scale Logical; if `TRUE`, scale selectivity curves to a maximum of
#'   one.
#'
#' @return A named list with one gear object per element. Each gear contains:
#' \describe{
#'   \item{name}{Gear name.}
#'   \item{sel_len}{Selectivity-at-length as an `FLQuant`.}
#'   \item{sel_a}{ALK-integrated selectivity-at-age as an `FLQuant`.}
#'   \item{len_bins}{Data frame of length-bin limits and midpoints.}
#'   \item{condALK}{Matrix containing \eqn{P(l \mid a,g)}.}
#'   \item{f_mult, ess_age, ess_len}{Values supplied in `gear_cfg`.}
#'   \item{lhpar}{Life-history parameters supplied in `lhpar`.}
#' }
#'
#' @details
#' Integrated selectivity-at-age is
#' \deqn{s_{g,a}=\sum_l P(l\mid a,t)s_g(l).}
#' The conditional inverse age-length key is
#' \deqn{P(l\mid a,g)=\frac{P(l\mid a,t)s_g(l)}{s_{g,a}}.}
#' The length grids created by [sel_la()] and [iALK()] are reconciled at the
#' upper plus group before these quantities are calculated.
#'
#' @examples
#' \dontrun{
#' lhpars <- FLPar(linf = 45, k = 0.4, t0 = -0.3, a = 0.00000612, b = 3.03)
#' ages   <- 0:5
#'
#' gear_cfg <- list(
#'   Ringnet = list(type = "dnormal", Lpeak = 7.3, sd_left = 0.7,
#'                   sd_right = 11.4, f_mult = 0.70, ess_age = 50, ess_len = 500),
#'   Gillnet = list(type = "normal", Lpeak = 39.7, sd = 10.6,
#'                   f_mult = 0.30, ess_age = 30, ess_len = 300)
#' )
#'
#' om_gears <- build_gear(gear_cfg, lhpar = lhpars, age = ages,
#'                         timing = 0.5)
#' }
#'
#' @export
build_gear <- function(gear_cfg, lhpar, age,
                       amin = min(age), amax = max(age),
                       lmin = 5, lmax_mult = 1.2, bin = 1, cv = 0.1,
                       timing = 0, scale = TRUE) {
  
  build_one <- function(cfg, name) {
    
    reserved <- c("type", "f_mult", "ess_age", "ess_len")
    missing_req <- setdiff(reserved, names(cfg))
    if (length(missing_req)) {
      stop(
        "Gear '", name, "' is missing required entries: ",
        paste(missing_req, collapse = ", "),
        ". Every gear's config list must include type/f_mult/ess_age/ess_len ",
        "plus that type's own sel_la() selectivity parameters."
      )
    }
    
    type     <- match.arg(cfg$type, c("logistic", "normal", "dnormal"))
    sel_pars <- cfg[setdiff(names(cfg), reserved)]
    f_mult   <- cfg$f_mult
    ess_age  <- cfg$ess_age
    ess_len  <- cfg$ess_len
    
    linf <- c(lhpar["linf"])
    ## sel_la()'s 'lmax' appears to be an EXCLUSIVE upper bound (its last
    ## bin's lower edge sits at lmax - bin, one short of iALK()'s
    ## plus-group edge). Padding by one bin width here makes the two grids
    ## align natively in the common case; the reconciliation block below
    ## remains a safety net for cases where they still don't.
    lmax_abs <- ceiling(linf * lmax_mult) + bin
    
    ## --- selectivity-at-length / -at-age ---------------------------------
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
    
    ## --- biological inverse ALK, same grid -------------------------------
    ## t0 - timing: growth evaluated 'timing' years past each integer age
    ## on the unmodified age grid (timing = 0 -> exactly at age, the old
    ## default; see the 'timing' argument doc for the sign convention).
    ialk <- iALK(
      params = c(linf = linf, k = c(lhpar["k"]), t0 = c(lhpar["t0"]) - timing),
      age = age, cv = cv, lmin = lmin, lmax = lmax_mult, bin = bin
    )
    
    ## --- reconcile length grids -------------------------------------------
    ## sel_la() and iALK() can legitimately disagree by one bin at the top:
    ## iALK() treats its last bin as a plus-group (P(length >= edge)),
    ## while sel_la() does not build a plus-group and simply stops one bin
    ## short. This is a real difference in convention, not necessarily a
    ## misconfiguration, so it is reconciled rather than treated as an
    ## error: any length present in iALK()'s grid but not sel_la()'s (the
    ## plus-group edge) is covered by carrying the last explicitly modelled
    ## selectivity value forward, on the standard assumption that
    ## selectivity has plateaued by 'lmax_mult * linf'. Any length present
    ## in sel_la()'s grid but not iALK()'s is dropped, with a warning. A
    ## genuinely larger mismatch (more than a couple of bins, or a gap in
    ## the middle of the range) is NOT something this should silently
    ## absorb, so that case still stops.
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
        ") forward to cover it. Confirm that selectivity has plateaued ",
        "by this length."
      )
      sel_len_vec[extra_in_ialk] <- tail(sel_len_vec, 1)
    }
    
    ## reorder to match iALK()'s column order exactly before using sweep()
    sel_len_vec <- sel_len_vec[len_ialk_chr]
    
    ## --- condition the ALK on capture by this gear -----------------------
    ialk_mat <- c(ialk)                       # coerce FLPar -> plain matrix
    dim(ialk_mat) <- dim(ialk)[1:2]
    dimnames(ialk_mat) <- list(age = dimnames(ialk)$age, len = dimnames(ialk)$len)
    
    weighted <- sweep(ialk_mat, 2, sel_len_vec, "*")
    row_tot <- rowSums(weighted)
    cond <- weighted
    valid <- row_tot > 0
    cond[valid, ] <- weighted[valid, , drop = FALSE] / row_tot[valid]
    cond[!valid, ] <- 0
    
    ## --- integrated selectivity-at-age, NOT sel_la()'s own sel_a ---------
    ## row_tot (pre-normalisation) is exactly s_{g,a} = sum_l P(l|a,t)*s_g(l)
    ## -- the length-selectivity curve integrated over the length-at-age
    ## distribution at this gear's own timing, using the same timed,
    ## reconciled ALK as condALK. This is deliberately NOT sel_la()'s own
    ## sel_a output, which evaluates s_g() at a single length-at-age value
    ## (and may use a different, possibly incorrect, within-year timing
    ## internally) -- see the 'timing' argument doc above.
    sel_a_vec <- row_tot
    if (scale && max(sel_a_vec, na.rm = TRUE) > 0) {
      sel_a_vec <- sel_a_vec / max(sel_a_vec, na.rm = TRUE)
    }
    sel_a_out <- FLQuant(sel_a_vec, dimnames = list(age = rownames(ialk_mat)))
    
    ## return the RECONCILED grid throughout, not sel_la()'s original one --
    ## sel_len/len_bins must match condALK's columns exactly, or a naive
    ## plot(gear$sel_len) against gear$condALK later would silently
    ## reintroduce the same one-bin mismatch this function just resolved.
    len_bins_out <- data.frame(
      lower = as.numeric(len_ialk_chr),
      upper = as.numeric(len_ialk_chr) + bin,
      mid   = as.numeric(len_ialk_chr) + bin / 2
    )
    
    list(
      name = name,
      sel_len = FLQuant(sel_len_vec, dimnames = list(len = names(sel_len_vec))),
      sel_a = sel_a_out,
      len_bins = len_bins_out,
      condALK = cond,
      f_mult = f_mult,
      ess_age = ess_age,
      ess_len = ess_len,
      lhpar = lhpar      # needed as lfd.sim()'s default params = gear$lhpar
    )
  }
  
  setNames(
    lapply(names(gear_cfg), function(nm) build_one(gear_cfg[[nm]], nm)),
    names(gear_cfg)
  )
}


#' Calculate seasonally averaged selectivity-at-age
#'
#' Integrates selectivity-at-length over the length-at-age distribution at
#' seasonal midpoint timings and averages the resulting selectivity-at-age
#' curves across seasons.
#'
#' @param lhpar Life-history parameters containing `linf`, `k`, and `t0`.
#' @param sel_len Selectivity-at-length as an `FLQuant` or named numeric
#'   vector.
#' @param age Numeric vector of ages.
#' @param n_seasons Integer number of equal within-year seasons.
#' @param model Growth function or model object. Defaults to `vonbert`.
#' @param reflen Optional reference length passed to [iALK()].
#' @param cv Numeric coefficient of variation in length-at-age.
#' @param lmin Numeric lower limit of the length grid.
#' @param lmax_mult Numeric multiplier applied to `linf` for the upper length
#'   limit.
#' @param bin Numeric length-bin width.
#' @param scale Logical; if `TRUE`, scale the averaged curve to a maximum of
#'   one.
#'
#' @return A named numeric vector of seasonally averaged selectivity-at-age.
#'
#' @details
#' Seasonal midpoint `s` is evaluated at
#' `t = (s - 0.5) / n_seasons`. Within-year growth is represented by replacing
#' `t0` with `t0 - t` in the growth parameters.
#'
#' @examples
#' \dontrun{
#' sel_a <- sel_a_seasonal_avg(
#'   lhpars,
#'   om_gears$Trawl$sel_len,
#'   age = 0:5,
#'   n_seasons = 4
#' )
#' }
#'
#' @export
sel_a_seasonal_avg <- function(lhpar, sel_len, age, n_seasons,
                               model = vonbert, reflen = NULL, cv = 0.1,
                               lmin = 5, lmax_mult = 1.2, bin = 1, scale = TRUE) {
  
  sel_len_vec <- if (inherits(sel_len, "FLQuant")) {
    v <- as.numeric(sel_len); names(v) <- dimnames(sel_len)$len; v
  } else sel_len
  
  linf <- c(lhpar["linf"]); k <- c(lhpar["k"]); t0 <- c(lhpar["t0"])
  
  sel_a_mat <- sapply(seq_len(n_seasons), function(s) {
    t <- (s - 0.5) / n_seasons
    ialk_t <- iALK(
      params = c(linf = linf, k = k, t0 = t0 - t),   # note: t0 - t, see note above
      model = model, age = age, lmax = lmax_mult, reflen = reflen,
      bin = bin, lmin = lmin
    )
    ialk_mat <- c(ialk_t); dim(ialk_mat) <- dim(ialk_t)[1:2]
    dimnames(ialk_mat) <- list(age = dimnames(ialk_t)$age, len = dimnames(ialk_t)$len)
    
    sv <- sel_len_vec[colnames(ialk_mat)]; sv[is.na(sv)] <- 0
    as.numeric(ialk_mat %*% sv)
  })
  
  sel_a <- rowMeans(sel_a_mat)
  names(sel_a) <- as.character(age)
  
  if (scale) sel_a <- sel_a / max(sel_a)
  
  sel_a
}

#' Calculate seasonally weighted weight-at-age
#'
#' Calculates weight-at-age at seasonal midpoint timings and averages across
#' seasons using seasonal numbers-at-age as weights.
#'
#' @param n_season An `FLQuant` containing seasonal population or catch
#'   numbers-at-age.
#' @param params Life-history parameters containing `linf`, `k`, `t0`, `a`,
#'   and `b`.
#' @param gear Optional gear object returned by [build_gear()]. If supplied,
#'   the age-length key is conditioned on the gear's selectivity-at-length.
#' @param n_seasons Integer number of equal within-year seasons.
#' @param model Growth function or model object. Defaults to `vonbert`.
#' @param reflen Optional reference length passed to [iALK()].
#' @param cv Numeric coefficient of variation in length-at-age.
#' @param lmin Numeric lower limit of the length grid.
#' @param lmax_mult Numeric multiplier applied to `linf` for the upper length
#'   limit.
#' @param bin Numeric length-bin width.
#'
#' @return An `FLQuant` containing seasonally weighted weight-at-age by year
#'   and iteration.
#'
#' @details
#' When `gear` is supplied, seasonal catch numbers weight a
#' capture-conditioned age-length key. When `gear = NULL`, the biological key
#' is weighted by seasonal population numbers.
#'
#' @examples
#' \dontrun{
#' catch_wt <- wt_a_seasonal(
#'   seasonal_catch$Trawl,
#'   lhpars,
#'   gear = om_gears$Trawl,
#'   n_seasons = 4
#' )
#' }
#'
#' @export
wt_a_seasonal <- function(n_season, params, gear = NULL, n_seasons,
                          model = vonbert, reflen = NULL, cv = 0.1,
                          lmin = 5, lmax_mult = 1.2, bin = 1) {
  
  age <- an(dimnames(n_season)$age)
  years <- dimnames(n_season)$year
  its <- dims(n_season)$iter
  linf <- c(params["linf"]); k <- c(params["k"]); t0 <- c(params["t0"])
  a_par <- c(params["a"]); b_par <- c(params["b"])
  
  sel_len_vec <- if (!is.null(gear)) {
    v <- as.numeric(gear$sel_len); names(v) <- dimnames(gear$sel_len)$len; v
  } else NULL
  
  W_mat <- sapply(seq_len(n_seasons), function(s) {
    t <- (s - 0.5) / n_seasons
    ialk_t <- iALK(
      params = c(linf = linf, k = k, t0 = t0 - t),   # note: t0 - t, see note above
      model = model, age = age, lmax = lmax_mult, reflen = reflen,
      bin = bin, lmin = lmin
    )
    ialk_mat <- c(ialk_t); dim(ialk_mat) <- dim(ialk_t)[1:2]
    dimnames(ialk_mat) <- list(age = dimnames(ialk_t)$age, len = dimnames(ialk_t)$len)
    
    if (!is.null(sel_len_vec)) {
      sv <- sel_len_vec[colnames(ialk_mat)]; sv[is.na(sv)] <- 0
      ialk_mat <- condition_alk(ialk_mat, sv)
    }
    
    len_mid <- an(colnames(ialk_mat)) + bin / 2
    as.numeric(ialk_mat %*% (a_par * len_mid^b_par))
  })
  
  out <- FLQuant(
    NA_real_,
    dimnames = list(age = age, year = years, unit = "unique",
                    season = "all", area = "unique", iter = seq_len(its))
  )
  
  for (y in seq_along(years)) {
    for (i in seq_len(its)) {
      n_mat <- sapply(seq_len(n_seasons), function(s)
        as.numeric(n_season[, y, , s, , i]))
      
      num <- rowSums(n_mat * W_mat)
      den <- rowSums(n_mat)
      wt <- ifelse(den > 0, num / den, W_mat[, ceiling(n_seasons / 2)])
      
      out[, y, , , , i] <- wt
    }
  }
  
  out
}


# }}}

#' Heatmap of a gear's conditional inverse ALK
#'
#' ggplot equivalent of `image(gear$condALK)`: a length-at-age probability
#' heatmap, \eqn{P(l \mid a,g)}, with labelled axes
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


#' Compare single-timing and seasonal selectivity-at-age
#'
#' Plots the selectivity-at-age stored in a gear object together with the
#' corresponding seasonally averaged curve. Individual seasonal curves can
#' be displayed as reference lines.
#'
#' @param gear Gear object returned by [build_gear()].
#' @param lhpar Life-history parameters containing `linf`, `k`, and `t0`.
#' @param age Numeric vector of ages.
#' @param n_seasons Integer number of equal within-year seasons.
#' @param show_seasons Logical; if `TRUE`, add the individual seasonal curves.
#' @param model Growth function or model object. Defaults to `vonbert`.
#' @param reflen Optional reference length passed to [iALK()].
#' @param cv Numeric coefficient of variation in length-at-age.
#' @param lmin Numeric lower limit of the length grid.
#' @param lmax_mult Numeric multiplier applied to `linf` for the upper length
#'   limit.
#' @param bin Numeric length-bin width.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' plot_sel_a_seasonal(
#'   om_gears$Trawl,
#'   lhpars,
#'   age = 0:5,
#'   n_seasons = 12
#' )
#' }
#'
#' @export
plot_sel_a_seasonal <- function(gear, lhpar, age, n_seasons,
                                show_seasons = TRUE,
                                model = vonbert, reflen = NULL, cv = 0.1,
                                lmin = 5, lmax_mult = 1.2, bin = 1) {
  
  sel_len_vec <- as.numeric(gear$sel_len)
  names(sel_len_vec) <- dimnames(gear$sel_len)$len
  
  linf <- c(lhpar["linf"]); k <- c(lhpar["k"]); t0 <- c(lhpar["t0"])
  
  season_df <- do.call(rbind, lapply(seq_len(n_seasons), function(s) {
    t <- (s - 0.5) / n_seasons
    ialk_t <- iALK(
      params = c(linf = linf, k = k, t0 = t0 - t),
      model = model, age = age, lmax = lmax_mult, reflen = reflen,
      bin = bin, lmin = lmin
    )
    ialk_mat <- c(ialk_t); dim(ialk_mat) <- dim(ialk_t)[1:2]
    dimnames(ialk_mat) <- list(age = dimnames(ialk_t)$age, len = dimnames(ialk_t)$len)
    sv <- sel_len_vec[colnames(ialk_mat)]; sv[is.na(sv)] <- 0
    data.frame(age = age, sel = as.numeric(ialk_mat %*% sv), season = factor(s))
  }))
  
  sel_avg <- sel_a_seasonal_avg(lhpar, gear$sel_len, age, n_seasons,
                                model = model, reflen = reflen, cv = cv,
                                lmin = lmin, lmax_mult = lmax_mult, bin = bin)
  
  main_df <- rbind(
    data.frame(age = age, sel = as.numeric(gear$sel_a), method = "Single timing (naive)"),
    data.frame(age = age, sel = sel_avg, method = "Seasonally averaged")
  )
  
  p <- ggplot()
  
  if (show_seasons) {
    p <- p + geom_line(
      data = season_df, aes(x = age, y = sel, group = season),
      colour = "grey20", linewidth = 0.4, linetype = 3
    )
  }
  
  p +
    geom_line(data = main_df, aes(x = age, y = sel, colour = method), linewidth = 1) +
    geom_point(data = main_df, aes(x = age, y = sel, colour = method), size = 1.5) +
    scale_colour_manual(
      values = c("Single timing (naive)" = "#C0504D", "Seasonally averaged" = "#5B9BD5"),
      name = NULL
    ) +
    labs(
      x = "Age", y = "Selectivity",
      title = paste0(gear$name, ": naive vs seasonally-averaged selectivity-at-age"),
      subtitle = if (show_seasons) "Grey dotted lines: individual within-year season curves" else NULL
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}
