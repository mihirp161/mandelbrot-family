# Mandelbrot Family Fractals

Renders four fractals from the `z = z^2 + c` family as high-resolution PNGs.

## Requirements

```r
install.packages(c("viridis", "png", "parallel", "doParallel", "foreach", "Rcpp"))
```

A C++ compiler is required for Rcpp (GCC on Linux/macOS, Rtools on Windows).

## Usage

```r
source("fractals.R")
```

Outputs `mandelbrot.png`, `julia.png`, `burning_ship.png`, and `buddhabrot.png` to the working directory. The Buddhabrot takes several minutes; the others finish in under a minute.

## Parameters

Defined at the top of the script:

| Variable | Default | Description |
|---|---|---|
| `SiZE` | `1200` | Image dimensions in pixels |
| `DPi` | `150` | PNG resolution (DPI) |
| `MAX_iTER` | `500` | Max escape iterations (Buddhabrot uses `1500`) |

Each fractal function also accepts its own arguments (bounds, `max_iter`, etc.) for zooming or tweaking.
