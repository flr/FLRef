#' Logistic selectivity at length
#'
#' Computes asymptotic logistic selectivity as a function of length.
#' The curve is parameterised by the length at 50% selectivity (`L50`)
#' and the length at 95% selectivity (`L95`).
#'
#' @param L Numeric vector of lengths.
#' @param L50 Length at 50% selectivity.
#' @param L95 Length at 95% selectivity.
#'
#' @return Numeric vector of selectivity values between 0 and 1.
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
sel_dnormal_len <- function(L, Lpeak, sd_left, sd_right) {
  s <- ifelse(
    L <= Lpeak,
    exp(-0.5 * ((L - Lpeak) / sd_left)^2),
    exp(-0.5 * ((L - Lpeak) / sd_right)^2)
  )
  s / max(s, na.rm = TRUE)
}





sel_la <- function(lhpar,amin=0,amax=20,
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
    
    if(is.null(lmax)) {
      lmax <- ceiling(an(lhpar["linf"]) * lmax_mult)
    }
    
    lower <- seq(lmin, lmax - binwidth, by = binwidth)
    upper <- lower + binwidth
    mid   <- lower + binwidth / 2
    
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
  
  L <- bins$mid
  # Length-at-age
  age <- amin:amax
  l_a <-vonbert(an(lhpar["linf"]),an(lhpar["k"]),an(lhpar["t0"])+0.5,amin:amax)
  
  if(type == "logistic") {
    if(is.null(L50) || is.null(L95)) {
      stop("For logistic selectivity, provide L50 and L95.")
    }
    sel <- sel_logistic_len(L, L50 = L50, L95 = L95)
    sel_a <- sel_logistic_len(l_a, L50 = L50, L95 = L95)
  }
  
  #browser()
  if(type == "normal") {
    if(is.null(Lpeak) || is.null(sd)) {
      stop("For normal selectivity, provide Lpeak and sd.")
    }
    sel <- sel_normal_len(L, Lpeak = Lpeak, sd = sd)
    sel_a <- sel_normal_len(l_a, Lpeak = Lpeak, sd = sd)
  } 

  if(type == "dnormal") {
    if(is.null(Lpeak) || is.null(sd_left) || is.null(sd_right)) {
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
  
  if(scale) {
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
  
 
  
 list(sel_len=sel_l,sel_a=sel_a)
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
