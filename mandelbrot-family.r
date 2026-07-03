# -------------------------------------------------------------------------
# The Mandelbrot Family (Mandelbrot, Julia Set, Burning Ship, Buddhabrot)
#
# All four share the same equation family: z = z^2 + c
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(viridis)  
  library(png)    
  library(parallel)
  library(doParallel)
  library(foreach)
  library(Rcpp)
  library(ambient)
})

# -- Shared output config ---------------------------------------
SiZE     <- 1200L   # pixel dimensions (square)
DPi      <- 150L
MAX_iTER <- 500L    # escape iterations (Buddhabrot overrides to 1500)

# Normalise, log-scale for contrast, then render via image()
render <- function(m, file, palette_fn) {
  m   <- log1p(pmax(m, 0))
  rng <- range(m)
  if (diff(rng) > 0) m <- (m - rng[1]) / diff(rng)
  png(file, width = SiZE, height = SiZE, res = DPi, bg = "black")
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  # t(m): image() expects [x, y] — our matrix is [row=y, col=x]
  image(t(m), col = palette_fn(1024L), axes = FALSE, useRaster = TRUE)
  dev.off()
  message("\u2713  ", file)
}


# -- 1. Mandelbrot ---------------------------------------------
# For each c in the complex plane, iterate z = Z^2 + c from z = 0.
# Track how many steps before |z| > 2 (escape). interior = never escapes.
# Smooth (continuous) colouring removes the ugly iteration bands.

mandelbrot <- function(xmin     = -2.2,
                       xmax     =  0.8,
                       ymin     = -1.5,
                       ymax     =  1.5,
                       width    = SiZE,
                       height   = SiZE,
                       max_iter = MAX_iTER) {
  
  x <- seq(xmin, xmax, length.out = width)
  y <- seq(ymin, ymax, length.out = height)
  
  # Grid of c values; row = imaginary axis, col = real axis
  C <- outer(y, x, function(Im, re) complex(real = re, imaginary = Im))
  
  z      <- matrix(0 + 0i, height, width)
  counts <- matrix(0,      height, width)
  active <- matrix(TRUE,   height, width)   # points not yet escaped
  
  for (k in seq_len(max_iter)) {
    z[active] <- z[active]^2 + C[active]
    
    esc <- active & Mod(z) > 2
    # Smooth escape time: removes banding, gives continuous gradient
    counts[esc] <- k - log2(pmax(log2(pmax(Mod(z[esc]), 1.0)), 1e-10))
    active[esc] <- FALSE
    
    if (!any(active)) break
  }
  counts
}

render(
  mandelbrot(),
  "mandelbrot.png",
  function(n) magma(n, direction = -1)
)


# -- 2. Julia Set ----------------------------------------------
# Same iteration (z = Z^2 + c), but c is now a fixed constant.
# instead of asking "does this c escape?", we ask "does this starting
# point Z0 escape under this fixed rule?"
#
# Changing c by even a tiny amount produces a completely different fractal.
# cx = -0.7269, cy = 0.1889  ->  intricate spiral ("Siegel disk" Julia set)
# Try also: cx = -0.4, cy = -0.6  (the original from your script)

julia <- function(cx       = -0.7269,
                  cy       =  0.1889,
                  xmin     = -1.5,
                  xmax     =  1.5,
                  ymin     = -1.5,
                  ymax     =  1.5,
                  width    = SiZE,
                  height   = SiZE,
                  max_iter = MAX_iTER) {
  
  x <- seq(xmin, xmax, length.out = width)
  y <- seq(ymin, ymax, length.out = height)
  
  K <- complex(real = cx, imaginary = cy)   # fixed parameter
  z <- outer(y, x, function(Im, re) complex(real = re, imaginary = Im))
  
  counts <- matrix(0,    height, width)
  active <- matrix(TRUE, height, width)
  
  for (k in seq_len(max_iter)) {
    z[active] <- z[active]^2 + K
    
    esc <- active & Mod(z) > 2
    counts[esc] <- k - log2(pmax(log2(pmax(Mod(z[esc]), 1.0)), 1e-10))
    active[esc] <- FALSE
    
    if (!any(active)) break
  }
  counts
}

render(
  julia(),
  "julia.png",
  function(n) plasma(n, direction = -1)
)


# -- 3. Burning Ship -------------------------------------------
# One modification: fold z back into the first quadrant before squaring.
#   z_{n+1} = (|Re(z)| + i|Im(z)|)^2 + c
# This breaks the vertical symmetry of the Mandelbrot and produces
# a shape that looks eerily like a fleet of ships on fire.
# The "hull" sits near c ≈ -1.755 + 0.028i in the lower half-plane.

burning_ship <- function(xmin     = -2.5,
                         xmax     =  1.5,
                         ymin     = -2.0,
                         ymax     =  0.5,
                         width    = SiZE,
                         height   = SiZE,
                         max_iter = MAX_iTER) {
  
  x <- seq(xmin, xmax, length.out = width)
  y <- seq(ymin, ymax, length.out = height)
  C <- outer(y, x, function(Im, re) complex(real = re, imaginary = Im))
  
  z      <- matrix(0 + 0i, height, width)
  counts <- matrix(0,      height, width)
  active <- matrix(TRUE,   height, width)
  
  for (k in seq_len(max_iter)) {
    za        <- z[active]
    z[active] <- complex(real      = abs(Re(za)),
                         imaginary = abs(Im(za)))^2 + C[active]
    
    esc <- active & Mod(z) > 2
    counts[esc] <- k - log2(pmax(log2(pmax(Mod(z[esc]), 1.0)), 1e-10))
    active[esc] <- FALSE
    
    if (!any(active)) break
  }
  counts
}

render(
  burning_ship(),
  "burning_ship.png",
  function(n) inferno(n, direction = -1)
)


# -- 4. Buddhabrot ---------------------------------------------
# inverts the Mandelbrot question entirely.
# instead of: "does this point escape?" ->  colour it.
# We ask:     "where did escaping points travel?" ->  record the path.
#
# Sample millions of random c values. Discard those that never escape
# (they are the Mandelbrot interior). For those that do escape, record
# every z the orbit visited and accumulate into a density map.
# High-traffic regions glow bright. Low-traffic regions stay dark.
#
# The result, at enough samples, resembles a seated figure.
# Named "Buddhabrot" by Melinda Green (1993).
#
# C++ (via Rcpp) handles the inner orbit loop for speed.
# Parallel workers distribute the random sampling.

buddhabrot_cpp_src <- '
#include <Rcpp.h>
#include <complex>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<std::complex<double>> orbit(std::complex<double> c, int max_iter) {
  std::complex<double> z(0.0, 0.0);
  std::vector<std::complex<double>> path;
  path.reserve(max_iter);

  for (int i = 0; i < max_iter; ++i) {
    z = z * z + c;
    if (std::abs(z) > 2.0) return path;   // escaped: return the orbit so far
    path.push_back(z);
  }
  return {};   // did not escape: belongs to the interior — discard
}
'

sourceCpp(code = buddhabrot_cpp_src)

buddhabrot <- function(n_points = 2e6,
                       size     = SiZE,
                       max_iter = 1500L) {
  
  density <- matrix(0L, size, size)
  
  # Map complex coordinate -> integer pixel index
  to_px <- function(z_vec) {
    col <- pmin(pmax(floor((Re(z_vec) + 2) / 4 * (size - 1)) + 1L, 1L), size)
    row <- pmin(pmax(floor((Im(z_vec) + 2) / 4 * (size - 1)) + 1L, 1L), size)
    cbind(row, col)
  }
  
  cores <- max(1L, detectCores() - 1L)
  cl    <- makeCluster(cores)
  registerDoParallel(cl)
  on.exit({ stopCluster(cl); closeAllConnections() }, add = TRUE)
  
  clusterExport(cl, c("buddhabrot_cpp_src", "max_iter"), envir = environment())
  clusterEvalQ(cl, { library(Rcpp); sourceCpp(code = buddhabrot_cpp_src) })
  
  batch_size  <- 1e5L
  n_batches   <- ceiling(n_points / batch_size)
  t0          <- proc.time()
  
  for (b in seq_len(n_batches)) {
    
    orbits <- foreach(i        = seq_len(batch_size),
                      .combine = "c",
                      .export  = "max_iter") %dopar% {
                        c_pt <- complex(real = runif(1L, -2, 2), imaginary = runif(1L, -2, 2))
                        list(orbit(c_pt, max_iter))
                      }
    
    valid <- unlist(orbits, recursive = FALSE)
    valid <- valid[lengths(valid) > 0L]
    
    if (length(valid)) {
      px           <- do.call(rbind, lapply(valid, to_px))
      density[px]  <- density[px] + 1L
    }
    gc()
    cat(sprintf("  batch %d / %d\n", b, n_batches))
  }
  
  elapsed <- (proc.time() - t0)["elapsed"]
  cat(sprintf("Done in %.1f min\n", elapsed / 60))
  density
}

# Custom palette: deep black -> violet -> electric pink -> white
buddhabrot_pal <- function(n) {
  colorRampPalette(
    c("black", "#0d0221", "#190c4e", "#4b1571",
      "#9c1de7", "#e040fb", "#ffd6ff", "white")
  )(n)
}

render(buddhabrot(), "buddhabrot.png", buddhabrot_pal)


# -------------------------------------------------------------------------
# Fractal Procedural Terrain (fBm + Ridged Multi-Fractal)


suppressPackageStartupMessages({
  library(ambient)
  library(png)
})

# -- Configuration ----------------------------------------------
SiZE     <- 1200L   # pixel dimensions (square)
DPi      <- 150L

# -- The Fractal Terrain Generator ------------------------------
# 1 Base fBm: Iterates smooth noise to form continents and plains.
# 2 Ridged Fractal: Takes the absolute value of the noise at each step, 
#    folding the waves into sharp, jagged peaks for mountain ranges.
stylegan_terrain <- function(width       = SiZE,
                             height      = SiZE,
                             octaves     = 9,
                             persistence = 0.5,
                             lacunarity  = 2.0) {
  
  # Create the Cartesian grid
  grid <- long_grid(x = seq(0, 4.5, length.out = width),
                    y = seq(0, 4.5, length.out = height))
  
  # Initialize two parallel fractal layers
  base_vec  <- rep(0, width * height)
  ridge_vec <- rep(0, width * height)
  
  amplitude <- 1.0
  frequency <- 1.0
  
  # FRACTAL ITERATION LOOP
  for (k in seq_len(octaves)) {
    
    # Generate the raw Simplex noise wave at the current frequency
    raw_noise <- gen_simplex(grid$x * frequency, grid$y * frequency)
    
    # Accumulate smooth continents
    base_vec <- base_vec + (raw_noise * amplitude)
    
    # Accumulate sharp mountains: (1 - |noise|)^2
    sharp_noise <- (1.0 - abs(raw_noise))^2
    ridge_vec   <- ridge_vec + (sharp_noise * amplitude)
    
    # Scale parameters for the next iteration
    amplitude <- amplitude * persistence
    frequency <- frequency * lacunarity
  }
  
  # Normalize both layers to [0, 1]
  base_norm  <- (base_vec - min(base_vec)) / diff(range(base_vec))
  ridge_norm <- (ridge_vec - min(ridge_vec)) / diff(range(ridge_vec))
  
  # MASKING & SHAPING
  # base^2.5 flattens valleys into oceans.
  # Multiplying ridges by the base ensures mountains only grow inland.
  final_terrain <- (base_norm^2.5) + ((base_norm^1.5) * ridge_norm * 1.8)
  
  # Final normalization for color mapping
  final_terrain <- (final_terrain - min(final_terrain)) / diff(range(final_terrain))
  matrix(final_terrain, nrow = height, ncol = width)
}

# -- Dedicated Terrain Renderer ---------------------------------
render_terrain <- function(m, file, palette_fn) {
  png(file, width = SiZE, height = SiZE, res = DPi, bg = "black")
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  
  image(t(m), col = palette_fn(2048L), axes = FALSE, useRaster = TRUE)
  
  dev.off()
  message("\u2713  ", file)
}

# -- Topographic Palette ------------------------------
stylegan_pal <- function(n) {
  colorRampPalette(
    c("#0b1626", "#172d47", "#294a6e", "#3a6894",   # Deep Ocean / Shallows
      "#baa37b", "#a18a60", "#7a6a4b",              # Plains / Heartlands
      "#4a5c37", "#324022", "#1f2b13",              # Forests
      "#464a4d", "#646b70", "#868e94",              # Rock / Mountain Ridges
      "#dce1e6", "#ffffff")                         # Snowcaps
  )(n)
}

render_terrain(
  stylegan_terrain(), 
  "perlin_landscape.png", 
  stylegan_pal
)
