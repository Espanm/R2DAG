#' Plot (Net) Pairwise Directional Connectedness as a Network
#'
#' Visualizes a square connectedness / spillover matrix as a directed network using
#' **NPDC** (net pairwise directional connectedness): \eqn{\mathrm{NPDC}_{ij} = S_{ij} - S_{ji}}.
#'
#' If `spillover_post` is provided, the function plots either:
#' - **two panels** (pre vs post), or
#' - **three panels** (pre vs post + difference), controlled by `label_diff`.
#'
#' Node size reflects total connectedness (TO + FROM, off-diagonal sums). Node pies show
#' the balance between TO vs FROM (red = more TO, green = more FROM).
#'
#' @param spillover_pre Square numeric matrix (p x p) with column names (node labels).
#' @param spillover_post Optional square numeric matrix (p x p) with identical colnames as `spillover_pre`.
#'   Set to `FALSE` (default) for a single-panel plot.
#' @param label_pre Character title for the first (pre) panel.
#' @param label_post Character title for the second (post) panel. If `FALSE` and `spillover_post` is provided,
#'   defaults to `"Second period"`.
#' @param label_diff Character title for the difference panel. If `FALSE` (default), no difference panel is shown.
#' @param file Optional file path (character) to save the plot as PNG. If `NULL` (default), plots to the active device.
#' @param weights Logical. If `TRUE`, edge weights are printed as edge labels.
#'
#' @return Invisibly returns `NULL`. Called for its side effect (plot).
#'
#' @examples
#' # Example (single panel):
#' # networkplot(spillover_pre = S1, label_pre = "Pre")
#' #
#' # Example (two panels):
#' # networkplot(S1, S2, label_pre = "Pre", label_post = "Post")
#' #
#' # Example (three panels with diff):
#' # networkplot(S1, S2, label_pre = "Pre", label_post = "Post", label_diff = "Diff")
#' #
#' # Save to file:
#' # networkplot(S1, label_pre = "Pre", file = "net_pre.png")
#'
#' @export
networkplot <- function(spillover_pre,
                        spillover_post = FALSE,
                        label_pre,
                        label_post = FALSE,
                        label_diff = FALSE,
                        file = NULL,
                        weights = FALSE) {

  # --- checks ---
  spillover_pre <- as.matrix(spillover_pre)
  if (nrow(spillover_pre) != ncol(spillover_pre)) stop("`spillover_pre` must be a square matrix.")

  names_ <- colnames(spillover_pre)
  if (is.null(names_)) stop("`spillover_pre` must have colnames (used as node labels).")

  has_post <- !isFALSE(spillover_post)
  if (has_post) {
    spillover_post <- as.matrix(spillover_post)
    if (nrow(spillover_post) != ncol(spillover_post)) stop("`spillover_post` must be a square matrix.")
    if (!identical(colnames(spillover_post), names_)) {
      stop("`spillover_post` colnames must match `spillover_pre` colnames.")
    }
  }

  if (isFALSE(label_post) && has_post) label_post <- "Second period"
  # label_diff can remain FALSE: controls whether we show 2 panels or 3 panels

  # --- helpers ---
  thrsh_topmean <- function(data, top_n_avg = FALSE, top = 250, pct = 0.15) {
    if (top_n_avg) {
      data <- as.matrix(data)
      top_n_mean <- mean(sort(as.vector(data), decreasing = TRUE)[1:top])
      data_thrsh <- ifelse(data < top_n_mean, 0, data)
      limit <- top_n_mean
    } else {
      threshold <- 1 - pct
      pos <- ifelse(data == 0, NA, data)
      limit <- as.numeric(stats::quantile(pos, threshold, na.rm = TRUE))
      data_thrsh <- ifelse(data < limit, 0, data)
    }
    list(data_thrsh = data_thrsh, limit = limit)
  }

  thrsh2 <- function(data, limit) {
    list(data_thrsh = ifelse(data < limit, 0, data))
  }

  scale_netspill <- function(x, NAMES, scale_ = FALSE, scaler_ = 5) {
    x_ <- as.matrix(x) * (-1)
    x_ <- ifelse(x_ < 0, 0, x_)
    colnames(x_) <- NAMES
    rownames(x_) <- NAMES
    diag(x_) <- 0
    if (scale_) {
      x_ <- x_ - min(x_, na.rm = TRUE)
      mx <- max(x_, na.rm = TRUE)
      if (mx > 0) x_ <- x_ / mx
      return(x_ * scaler_)
    }
    x_
  }

  size_pie_map <- function(ctable, NAMES, scale_ = TRUE, rescale_ = c(5, 10)) {
    ctable <- as.matrix(ctable)
    colnames(ctable) <- NAMES
    rownames(ctable) <- NAMES

    diag_ <- diag(diag(ctable))
    ctable_offdiag <- ctable - diag_

    vsize_map_TO <- colSums(ctable_offdiag)
    vsize_map_FROM <- rowSums(ctable_offdiag)
    vsize_map <- vsize_map_TO + vsize_map_FROM
    names(vsize_map) <- NAMES

    pie_map <- ifelse(vsize_map == 0, 0.5, vsize_map_TO / vsize_map)

    if (scale_) {
      rng <- range(vsize_map, na.rm = TRUE)
      if (rng[1] == rng[2]) {
        vsize_map_norm <- rep(mean(rescale_), length(vsize_map))
      } else {
        vsize_map_norm <- (vsize_map - rng[1]) / (rng[2] - rng[1])
        vsize_map_norm <- rescale_[1] + vsize_map_norm * (rescale_[2] - rescale_[1])
      }
      return(list(vsize_map = vsize_map_norm, pie_map = pie_map))
    } else {
      return(list(vsize_map = vsize_map, pie_map = pie_map))
    }
  }

  size_pie_map2 <- function(ctable1, ctable2, NAMES, scale_ = TRUE, rescale_ = c(14, 22)) {

    c1 <- as.matrix(ctable1); colnames(c1) <- NAMES; rownames(c1) <- NAMES
    c2 <- as.matrix(ctable2); colnames(c2) <- NAMES; rownames(c2) <- NAMES

    off1 <- c1 - diag(diag(c1))
    off2 <- c2 - diag(diag(c2))

    v1 <- colSums(off1) + rowSums(off1)
    v2 <- colSums(off2) + rowSums(off2)

    v_all <- c(v1, v2)

    if (scale_) {
      rng <- range(v_all, na.rm = TRUE)
      if (rng[1] == rng[2]) {
        v_all_sc <- rep(mean(rescale_), length(v_all))
      } else {
        v01 <- (v_all - rng[1]) / (rng[2] - rng[1])
        v_all_sc <- rescale_[1] + v01 * (rescale_[2] - rescale_[1])
      }
      v1_sc <- v_all_sc[1:length(v1)]
      v2_sc <- v_all_sc[(length(v1) + 1):(2 * length(v1))]
      names(v1_sc) <- NAMES
      names(v2_sc) <- NAMES
      return(list(vsize_map1 = v1_sc, vsize_map2 = v2_sc))
    } else {
      names(v1) <- NAMES; names(v2) <- NAMES
      return(list(vsize_map1 = v1, vsize_map2 = v2))
    }
  }

  # --- edge label args: only pass if weights=TRUE ---
  edge_label_args <- function(M, weights, digits = 3, cex = 1.8) {
    if (!weights) return(list())
    M <- as.matrix(M)
    lab <- matrix("", nrow(M), ncol(M))
    idx <- !is.na(M) & (M != 0)
    if (any(idx)) {
      lab[idx] <- format(round(M[idx], digits), nsmall = digits, trim = TRUE)
    }
    diag(lab) <- ""
    dimnames(lab) <- dimnames(M)

    list(
      edge.labels = lab,
      edge.label.cex = cex,
      edge.label.color = "black"
    )
  }

  # --- NPDC ---
  npdc_pre_raw <- spillover_pre - t(spillover_pre)
  NPDC_pre <- scale_netspill(npdc_pre_raw, NAMES = names_, scale_ = FALSE)

  if (has_post) {
    npdc_post_raw <- spillover_post - t(spillover_post)
    NPDC_post <- scale_netspill(npdc_post_raw, NAMES = names_, scale_ = FALSE)
    NPDC_diff <- scale_netspill((npdc_post_raw - npdc_pre_raw), NAMES = names_, scale_ = FALSE)
  }

  # --- threshold (shared limit) ---
  pct_keep <- 1
  limit1 <- thrsh_topmean(NPDC_pre, top_n_avg = FALSE, pct = pct_keep)$limit
  if (has_post) {
    limit2 <- thrsh_topmean(NPDC_post, top_n_avg = FALSE, pct = pct_keep)$limit
    lim_use <- min(limit1, limit2)
  } else {
    lim_use <- limit1
  }

  NPDC_pre_thrsh <- thrsh2(NPDC_pre, lim_use)$data_thrsh
  if (has_post) NPDC_post_thrsh <- thrsh2(NPDC_post, lim_use)$data_thrsh

  # --- edge size scaling ---
  max_esize <- 7
  if (has_post) {
    erate  <- max(NPDC_post_thrsh, na.rm = TRUE) / max(NPDC_pre_thrsh, na.rm = TRUE) * max_esize
    erate_ <- max(NPDC_pre_thrsh,  na.rm = TRUE) / max(NPDC_post_thrsh, na.rm = TRUE) * max_esize
    post_esize <- ifelse(erate  > max_esize, max_esize, erate)
    pre_esize  <- ifelse(erate_ > max_esize, max_esize, erate_)
  } else {
    pre_esize <- max_esize
  }

  # --- pies + node sizes ---
  pie_map_pre <- size_pie_map(ctable = spillover_pre, NAMES = names_)$pie_map

  if (has_post) {
    vsize_pair <- size_pie_map2(ctable1 = spillover_pre, ctable2 = spillover_post,
                                NAMES = names_, rescale_ = c(14, 22))
    vsize_map_pre  <- vsize_pair$vsize_map1
    vsize_map_post <- vsize_pair$vsize_map2

    pie_map_post <- size_pie_map(ctable = spillover_post, NAMES = names_)$pie_map
  } else {
    vsize_map_pre <- size_pie_map(ctable = spillover_pre, NAMES = names_, rescale_ = c(14, 22))$vsize_map
  }

  # --- colors ---
  col_ <- tryCatch({
    pal_simpsons("springfield")(16)[c(7)]
  }, error = function(e) "grey80")

  groups <- list("First" = seq_along(names_))

  label_font <- rep(2, length(names_))
  label_cex  <- rep(1.1, length(names_))

  label_color2 <- rep("black", length(names_))
  for (i in seq_along(names_)) {
    if (pie_map_pre[i] > 0.55) label_color2[i] <- "red2"
    if (pie_map_pre[i] < 0.45) label_color2[i] <- "green3"
  }

  if (has_post) {
    label_color3 <- rep("black", length(names_))
    for (i in seq_along(names_)) {
      if (pie_map_post[i] > 0.55) label_color3[i] <- "red2"
      if (pie_map_post[i] < 0.45) label_color3[i] <- "green3"
    }
  }

  pieColor_pre  <- rep("red2",   length(names_))
  pieColor2_pre <- rep("green3", length(names_))
  if (has_post) {
    pieColor_post  <- rep("red2",   length(names_))
    pieColor2_post <- rep("green3", length(names_))
  }

  # --- qgraph uniform settings ---
  uniform_settings <- list(
    groups = groups,
    layout = "groups",
    layoutScale = c(1, 1),
    # NOTE: edge.labels are NOT included here; passed conditionally if weights=TRUE
    label.font = label_font,
    label.cex = label_cex,
    shape = "circle",
    labels = names_,
    curveAll = TRUE,
    legend.mode = "groups",
    color = rep("dodgerblue3", length(names_)),
    edge.color = rep("dodgerblue3", length(names_)),
    fade = FALSE,
    pieBorder = 0.25,
    label.fill.vertical = 0.98,
    label.fill.horizontal = 0.98,
    curve = 1.1,
    asize = 6,
    bidirectional = FALSE
  )

  # --- open device if needed ---
  if (!is.null(file)) {
    if (has_post) {
      grDevices::png(file, width = 2800, height = 1400)
    } else {
      grDevices::png(file, width = 1800, height = 1400)
    }
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  graphics::par(mai = rep(0.5, 4), mar = c(1, 1, 1, 1))

  # --- panels ---
  if (!has_post) {
    args_labels <- edge_label_args(NPDC_pre_thrsh, weights)

    do.call(qgraph::qgraph, c(
      uniform_settings,
      list(
        NPDC_pre_thrsh,
        maximum = max(NPDC_pre_thrsh, na.rm = TRUE),
        vsize = vsize_map_pre,
        label.color = label_color2,
        legend = FALSE,
        esize = pre_esize,
        pie = pie_map_pre,
        pieColor = pieColor_pre,
        pieColor2 = pieColor2_pre,
        pieBorder = 1
      ),
      args_labels
    ))
    graphics::text(x = -0.8, y = 1, labels = label_pre, xpd = NA, cex = 5, font = 2)
    return(invisible(NULL))
  }

  show_diff <- !isFALSE(label_diff)

  if (!show_diff) {
    graphics::layout(
      matrix(c(
        1,1,1, 0, 2,2,2,
        1,1,1, 0, 2,2,2,
        1,1,1, 0, 2,2,2
      ), ncol = 7, byrow = TRUE),
      widths = c(1, 1, 1, 0.35, 1, 1, 1)
    )

    # PRE
    args_labels_pre <- edge_label_args(NPDC_pre_thrsh, weights)
    do.call(qgraph::qgraph, c(
      uniform_settings,
      list(
        NPDC_pre_thrsh,
        maximum = max(NPDC_pre_thrsh, na.rm = TRUE),
        vsize = vsize_map_pre,
        label.color = label_color2,
        legend = FALSE,
        esize = pre_esize,
        pie = pie_map_pre,
        pieColor = pieColor_pre,
        pieColor2 = pieColor2_pre,
        pieBorder = 1
      ),
      args_labels_pre
    ))
    graphics::text(x = -0.8, y = 1, labels = label_pre, xpd = NA, cex = 5, font = 2)

    # POST
    args_labels_post <- edge_label_args(NPDC_post_thrsh, weights)
    do.call(qgraph::qgraph, c(
      uniform_settings,
      list(
        NPDC_post_thrsh,
        maximum = max(NPDC_post_thrsh, na.rm = TRUE),
        vsize = vsize_map_post,
        label.color = label_color3,
        legend = FALSE,
        esize = post_esize,
        pie = pie_map_post,
        pieColor = pieColor_post,
        pieColor2 = pieColor2_post,
        pieBorder = 1
      ),
      args_labels_post
    ))
    graphics::text(x = -1, y = 1, labels = label_post, xpd = NA, cex = 5, font = 2)

    return(invisible(NULL))
  }

  # 3 panel layout
  graphics::layout(matrix(c(
    1,1,1, 0, 2,2,2,
    1,1,1, 0, 2,2,2,
    1,1,1, 0, 2,2,2,
    0,0,3,3,3,0,0,
    0,0,3,3,3,0,0
  ), ncol = 7, byrow = TRUE))

  # PRE
  args_labels_pre <- edge_label_args(NPDC_pre_thrsh, weights)
  do.call(qgraph::qgraph, c(
    uniform_settings,
    list(
      NPDC_pre_thrsh,
      maximum = max(NPDC_pre_thrsh, na.rm = TRUE),
      vsize = vsize_map_pre,
      label.color = label_color2,
      legend = FALSE,
      esize = pre_esize,
      pie = pie_map_pre,
      pieColor = pieColor_pre,
      pieColor2 = pieColor2_pre,
      pieBorder = 1
    ),
    args_labels_pre
  ))
  graphics::text(x = -0.8, y = 1, labels = label_pre, xpd = NA, cex = 5, font = 2)

  # POST
  args_labels_post <- edge_label_args(NPDC_post_thrsh, weights)
  do.call(qgraph::qgraph, c(
    uniform_settings,
    list(
      NPDC_post_thrsh,
      maximum = max(NPDC_post_thrsh, na.rm = TRUE),
      vsize = vsize_map_post,
      label.color = label_color3,
      legend = FALSE,
      esize = post_esize,
      pie = pie_map_post,
      pieColor = pieColor_post,
      pieColor2 = pieColor2_post,
      pieBorder = 1
    ),
    args_labels_post
  ))
  graphics::text(x = -1, y = 1, labels = label_post, xpd = NA, cex = 5, font = 2)

  # DIFF
  pie_diff <- pie_map_post - pie_map_pre
  vsize_map_diff <- vsize_map_post - vsize_map_pre

  bordercolor <- rep("black", length(names_))
  bordersize  <- rep(3, length(names_))
  label_color <- rep("black", length(names_))
  vsize_map   <- rep(18, length(names_))

  for (n in seq_along(pie_diff)) {
    bordersize[n] <- ifelse(pie_diff[n] >  0.20, 6, bordersize[n])
    bordersize[n] <- ifelse(pie_diff[n] < -0.20, 6, bordersize[n])

    piecol <- ifelse(pie_diff[n] > 0, "red2", "green3")
    label_color[n] <- ifelse(pie_diff[n] > 0, "red2", "green3")
    bordercolor[n] <- ifelse(pie_diff[n] == 0, "black", piecol)

    vsize_map[n] <- ifelse(vsize_map_diff[n] >  3, 22, vsize_map[n])
    vsize_map[n] <- ifelse(vsize_map_diff[n] < -3, 14, vsize_map[n])
  }

  args_labels_diff <- edge_label_args(NPDC_diff, weights)
  do.call(qgraph::qgraph, c(
    uniform_settings,
    list(
      NPDC_diff,
      maximum = max(NPDC_diff, na.rm = TRUE),
      pieBorder = 0,
      legend = FALSE,
      label.color = label_color,
      vsize = vsize_map,
      esize = 4,
      border.width = bordersize,
      border.color = bordercolor
    ),
    args_labels_diff
  ))
  graphics::text(x = -1, y = 1, labels = label_diff, xpd = NA, cex = 5, font = 2)

  invisible(NULL)
}
