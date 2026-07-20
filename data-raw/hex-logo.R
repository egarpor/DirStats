# -----------------------------------------------------------------------------
# DirStats hexagon logo
#
# A sphere (S^2) whose surface is colored by a *nonparametric* kernel density
# estimate -- computed with DirStats::kde_dir() itself -- from a clustered
# directional sample, with the sample points overlaid. Rendered as a 2D
# orthographic projection with ggplot2 (deterministic, no rgl/OpenGL needed),
# with Lambertian shading so the sphere reads as a lit 3D globe.
#
# Run from the package root:  Rscript data-raw/hex-logo.R
# -----------------------------------------------------------------------------

library(DirStats)
library(ggplot2)

# ---- Palette --------------------------------------------------------------
navy   <- "#0b1b2b"  # hexagon background
border <- "#57c79a"  # hexagon border accent (viridis-teal)
ink    <- "#eef3f5"  # wordmark, points, graticule

# ---- Geometry -------------------------------------------------------------
Rs   <- 0.95   # sphere radius on the canvas
cy   <- 0.14   # sphere vertical offset (leaves room for the wordmark)
Rhex <- 1.5    # hexagon circumradius
res  <- 700    # raster resolution (pixels per side)

# Orthonormal frame (e1, e2, pole) from a pole direction
frame_from_pole <- function(pole) {
  pole <- pole / sqrt(sum(pole^2))
  seed <- if (abs(pole[3]) < 0.9) c(0, 0, 1) else c(1, 0, 0)
  e1 <- seed - sum(seed * pole) * pole
  e1 <- e1 / sqrt(sum(e1^2))
  e2 <- c(pole[2] * e1[3] - pole[3] * e1[2],
          pole[3] * e1[1] - pole[1] * e1[3],
          pole[1] * e1[2] - pole[2] * e1[1])
  list(e1 = e1, e2 = e2, pole = pole)
}

# Project a matrix of unit vectors (view frame) onto the canvas, breaking the
# path where the sphere turns away from the viewer (z <= 0)
project_line <- function(P, grp) {
  hidden <- P[, 3] <= 0.02
  data.frame(X = ifelse(hidden, NA, Rs * P[, 1]),
             Y = ifelse(hidden, NA, cy + Rs * P[, 2]),
             grp = grp)
}

# Build the graticule (meridians + parallels) as broken paths
make_graticule <- function(pole = c(0.13, 0.86, 0.49), n = 220) {
  fr <- frame_from_pole(pole)
  e1 <- fr$e1; e2 <- fr$e2; pl <- fr$pole
  lines <- list()
  k <- 0
  # Parallels
  for (lat in seq(-60, 60, by = 30) * pi / 180) {
    lon <- seq(0, 2 * pi, length.out = n)
    P <- outer(cos(lat) * cos(lon), e1) + outer(cos(lat) * sin(lon), e2) +
      matrix(sin(lat) * pl, nrow = n, ncol = 3, byrow = TRUE)
    k <- k + 1
    lines[[k]] <- project_line(P, paste0("par", k))
  }
  # Meridians
  for (lon in seq(0, 330, by = 30) * pi / 180) {
    lat <- seq(-pi / 2, pi / 2, length.out = n)
    P <- outer(cos(lat) * cos(lon), e1) + outer(cos(lat) * sin(lon), e2) +
      outer(sin(lat), pl)
    k <- k + 1
    lines[[k]] <- project_line(P, paste0("mer", k))
  }
  do.call(rbind, lines)
}

# ---- Logo builder ---------------------------------------------------------
make_logo <- function(seed = 20240720, h = 0.3, graticule = TRUE,
                      palette = viridisLite::viridis, wordmark = TRUE) {

  # Clustered directional sample on the front hemisphere (view frame)
  set.seed(seed)
  mus <- rbind(normalize(c(-0.5,  0.55,  0.85)),
               normalize(c( 0.72, -0.1,   0.7)),
               normalize(c( 0.08,  0.82,  0.5)))
  kappas <- c(32, 24, 42)
  ns     <- c(46, 40, 46)
  samp <- do.call(rbind, lapply(seq_len(3), function(i)
    rotasym::r_vMF(n = ns[i], mu = mus[i, ], kappa = kappas[i])))

  # Front-hemisphere raster over the unit disk
  s <- seq(-1, 1, length.out = res)
  g <- expand.grid(u = s, v = s)
  r2 <- g$u^2 + g$v^2
  inside <- r2 <= 1
  w <- sqrt(pmax(0, 1 - r2[inside]))
  pts <- cbind(g$u[inside], g$v[inside], w)

  # Directional KDE (DirStats eats its own dog food)
  dens <- kde_dir(x = pts, data = samp, h = h)
  dn <- (dens - min(dens)) / (max(dens) - min(dens))
  dn <- 0.06 + 0.94 * dn                      # lift the floor off pure black

  # Base color from the palette
  cols <- palette(256)
  base <- t(grDevices::col2rgb(cols[pmax(1, pmin(256, round(dn * 255) + 1))])) / 255

  # Lambertian shading + specular highlight (light from the upper-left-front)
  light <- c(-0.5, 0.6, 0.75); light <- light / sqrt(sum(light^2))
  diffuse <- pmax(0, as.numeric(pts %*% light))
  shade <- 0.4 + 0.6 * diffuse
  half <- light + c(0, 0, 1); half <- half / sqrt(sum(half^2))
  spec <- pmax(0, as.numeric(pts %*% half))^30
  # Clamp with indexing (pmin() would drop the matrix dimensions)
  lit <- base * shade + spec * 0.28
  lit[lit > 1] <- 1

  g$fill <- NA_character_
  g$fill[inside] <- grDevices::rgb(lit[, 1], lit[, 2], lit[, 3])
  g$X <- Rs * g$u
  g$Y <- cy + Rs * g$v

  # Hexagon (pointy top/bottom, vertical sides: standard R sticker shape)
  ang <- (seq(0, 300, by = 60) + 30) * pi / 180
  hex <- data.frame(X = Rhex * cos(ang), Y = Rhex * sin(ang))

  # Sphere rim and visible sample points
  th <- seq(0, 2 * pi, length.out = 400)
  rim <- data.frame(X = Rs * cos(th), Y = cy + Rs * sin(th))
  front <- samp[, 3] > 0
  pd <- data.frame(X = Rs * samp[front, 1], Y = cy + Rs * samp[front, 2],
                   depth = samp[front, 3])

  # Assemble
  p <- ggplot() +
    geom_polygon(data = hex, aes(X, Y), fill = navy, color = NA) +
    geom_raster(data = g, aes(X, Y, fill = fill)) +
    scale_fill_identity()
  if (graticule) {
    p <- p + geom_path(data = make_graticule(), aes(X, Y, group = grp),
                       color = ink, alpha = 0.12, linewidth = 0.25,
                       na.rm = TRUE)
  }
  p <- p +
    geom_path(data = rim, aes(X, Y), color = ink, alpha = 0.45,
              linewidth = 0.5) +
    geom_point(data = pd, aes(X, Y, alpha = depth), color = ink,
               size = 0.8, show.legend = FALSE) +
    scale_alpha(range = c(0.45, 0.95)) +
    geom_polygon(data = hex, aes(X, Y), fill = NA, color = border,
                 linewidth = 3.2, lineend = "round", linejoin = "round")
  if (wordmark) {
    p <- p + annotate("text", x = 0, y = -Rhex * 0.7, label = "DirStats",
                      color = ink, family = "sans", fontface = "bold",
                      size = 12.5)
  }
  p +
    coord_fixed(xlim = c(-1.36, 1.36), ylim = c(-1.57, 1.57), expand = FALSE) +
    theme_void() +
    theme(plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.margin = margin(0, 0, 0, 0))
}

# ---- Render ---------------------------------------------------------------
# Fixed physical size (hex aspect sqrt(3)/2) so proportions are independent of
# resolution; dpi controls the pixel size. The output is losslessly shrunk with
# zopflipng afterwards, if available.
save_logo <- function(p, file, dpi = 360) {
  ggsave(file, p, width = 2.598, height = 3.0, units = "in", dpi = dpi,
         bg = "transparent")
  if (nzchar(Sys.which("zopflipng"))) {
    system2("zopflipng", c("-y", shQuote(file), shQuote(file)),
            stdout = FALSE, stderr = FALSE)
  }
}

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

# Primary logo
save_logo(make_logo(), "man/figures/logo.png")

# Variants for comparison (kept in data-raw, not shipped)
save_logo(make_logo(graticule = FALSE), "data-raw/logo-variant-nogrid.png")
save_logo(make_logo(palette = viridisLite::magma),
          "data-raw/logo-variant-magma.png")

message("Logo written to man/figures/logo.png")
