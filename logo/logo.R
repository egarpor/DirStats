# -----------------------------------------------------------------------------
# DirStats hexagon logo
#
# A sphere (S^2) whose surface is colored by a *nonparametric* kernel density
# estimate -- computed with DirStats::kde_dir() itself -- from a clustered
# directional sample, with the sample points overlaid. Rendered as a 2D
# orthographic projection with ggplot2 (deterministic, no rgl/OpenGL needed),
# with Lambertian shading so the sphere reads as a lit 3D globe. The globe is
# the sticker subplot; hexSticker::sticker() draws the hexagon, wordmark and URL
# so the font and GitHub-link placement match the other egarpor package logos.
#
# Run from the package root:  Rscript logo/logo.R
# Output: logo/logo.png (master) and man/figures/logo.png (shipped mirror).
# -----------------------------------------------------------------------------

library(DirStats)
library(ggplot2)
library(hexSticker)

## ---- Shared logo standard (identical across egarpor packages) -------------
# Aller_Rg is bundled with (and auto-registered by) hexSticker, so no
# showtext/font_add setup is needed for the sticker() idiom.
FONT    <- "Aller_Rg"   # typeface for the package name and the GitHub URL
P_SIZE  <- 31.2         # package-name size (shared across packages)
U_SIZE  <- 9.0          # GitHub URL size (large enough to read)
U_X     <- 1.00         # GitHub URL position: along the lower-right hex edge
U_Y     <- 0.08
U_ANGLE <- 30
H_SIZE  <- 1.5          # hexagon border thickness
DPI     <- 600

# ---- Palette --------------------------------------------------------------
navy   <- "#0b1b2b"  # hexagon background
border <- "#57c79a"  # hexagon border accent (viridis-teal)
ink    <- "#eef3f5"  # wordmark, points, graticule

# ---- Geometry -------------------------------------------------------------
# The globe is drawn as a unit sphere centred at the origin in the subplot's own
# coordinates; hexSticker insets and scales it via s_x/s_y/s_width.
Rs   <- 1      # sphere radius in the subplot's coordinates
cy   <- 0      # sphere centre (origin)
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

# ---- Globe subplot builder ------------------------------------------------
# Returns the lit KDE globe as a transparent ggplot (no hexagon, no text): the
# subplot passed to hexSticker::sticker().
make_globe <- function(seed = 20240720, h = 0.3, graticule = TRUE,
                       palette = viridisLite::viridis) {

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

  # Sphere rim and visible sample points
  th <- seq(0, 2 * pi, length.out = 400)
  rim <- data.frame(X = Rs * cos(th), Y = cy + Rs * sin(th))
  front <- samp[, 3] > 0
  pd <- data.frame(X = Rs * samp[front, 1], Y = cy + Rs * samp[front, 2],
                   depth = samp[front, 3])

  # Assemble the globe (transparent background -> the navy hex shows through)
  p <- ggplot() +
    geom_raster(data = g, aes(X, Y, fill = fill)) +
    scale_fill_identity()
  if (graticule) {
    p <- p + geom_path(data = make_graticule(), aes(X, Y, group = grp),
                       color = ink, alpha = 0.12, linewidth = 0.25,
                       na.rm = TRUE)
  }
  p +
    geom_path(data = rim, aes(X, Y), color = ink, alpha = 0.45,
              linewidth = 0.5) +
    geom_point(data = pd, aes(X, Y, alpha = depth), color = ink,
               size = 0.8, show.legend = FALSE) +
    scale_alpha(range = c(0.45, 0.95)) +
    coord_fixed(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05), expand = FALSE) +
    theme_void() +
    theme(plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.margin = margin(0, 0, 0, 0))
}

# ---- Render ---------------------------------------------------------------
# hexSticker draws the hexagon, the "DirStats" wordmark and the GitHub URL with
# the shared font/size/position; the globe is the inset subplot.
sticker_logo <- function(globe, file) {
  sticker(
    subplot = globe, s_x = 1, s_y = 1.14, s_width = 1.18, s_height = 1.18,
    package = "DirStats", p_x = 1, p_y = 0.40, p_size = P_SIZE,
    p_color = ink, p_family = FONT,
    url = "github.com/egarpor/DirStats",
    u_x = U_X, u_y = U_Y, u_angle = U_ANGLE, u_size = U_SIZE,
    u_color = "#9fb6c2", u_family = FONT,
    h_fill = navy, h_color = border, h_size = H_SIZE,
    spotlight = FALSE, white_around_sticker = FALSE,
    dpi = DPI, filename = file
  )
  invisible(file)
}

dir.create("logo", showWarnings = FALSE)
dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

# Primary logo (master + shipped mirror)
sticker_logo(make_globe(), "logo/logo.png")
file.copy("logo/logo.png", "man/figures/logo.png", overwrite = TRUE)

# Variants for comparison (kept in logo/, not shipped)
sticker_logo(make_globe(graticule = FALSE), "logo/logo-variant-nogrid.png")
sticker_logo(make_globe(palette = viridisLite::magma), "logo/logo-variant-magma.png")

message("Wrote logo/logo.png and man/figures/logo.png")
