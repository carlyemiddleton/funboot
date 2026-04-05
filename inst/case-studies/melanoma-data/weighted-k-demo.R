set.seed(123)

library(ggplot2)

# -----------------------------
# 1. Create a toy B intensity field on a pixel grid
# -----------------------------
img_size <- 100
x_grid <- seq(0, 1, length.out = img_size)
y_grid <- seq(0, 1, length.out = img_size)

grid_df <- expand.grid(x = x_grid, y = y_grid)

B_fun <- function(x, y) {
  2.0 * exp(-((x - 0.72)^2 + (y - 0.68)^2) / 0.015) +
    1.2 * exp(-((x - 0.30)^2 + (y - 0.28)^2) / 0.030)
}

grid_df$B <- with(grid_df, B_fun(x, y))

# pixel area
dx <- diff(x_grid)[1]
dy <- diff(y_grid)[1]
delta_a <- dx * dy

# -----------------------------
# 2. Simulate Type A cells that are more likely in high-B regions
# -----------------------------
n_candidates <- 5000
cand <- data.frame(
  x = runif(n_candidates),
  y = runif(n_candidates)
)

cand$B <- with(cand, B_fun(x, y))
cand$accept_prob <- cand$B / max(cand$B)
cand$keep <- runif(n_candidates) < cand$accept_prob * 0.20

cells_A <- cand[cand$keep, c("x", "y")]
rownames(cells_A) <- NULL

# cap number of cells just for a clean toy figure
n_A <- min(100, nrow(cells_A))
cells_A <- cells_A[seq_len(n_A), ]

# -----------------------------
# 3. Weighted point-field cross-K
#    K_AB(r) = average accumulated B intensity
#              within distance r of A cells
# -----------------------------
weighted_cross_K_field <- function(cells, pixels, r_vals, delta_a) {
  
  nA <- nrow(cells)
  pixel_coords <- as.matrix(pixels[, c("x","y")])
  pixel_weights <- pixels$B * delta_a
  
  K_vals <- numeric(length(r_vals))
  
  for(i in seq_len(nA)) {
    
    # distance from this cell to all pixels
    dx <- pixel_coords[,1] - cells$x[i]
    dy <- pixel_coords[,2] - cells$y[i]
    d  <- sqrt(dx^2 + dy^2)
    
    for(k in seq_along(r_vals)) {
      r <- r_vals[k]
      inside <- d <= r
      K_vals[k] <- K_vals[k] + sum(pixel_weights[inside])
    }
    
  }
  
  K_vals <- K_vals / nA
  
  data.frame(r = r_vals, Kab = K_vals)
}

r_vals <- seq(0.01, 0.25, by = 0.01)
K_obs <- weighted_cross_K_field(
  cells = cells_A,
  pixels = grid_df,
  r_vals = r_vals,
  delta_a = delta_a
)

# -----------------------------
# 4. Null comparison:
#    randomize A-cell locations uniformly on the image
# -----------------------------
n_perm <- 199
K_perm <- matrix(NA_real_, nrow = length(r_vals), ncol = n_perm)

for (b in seq_len(n_perm)) {
  cells_null <- data.frame(
    x = runif(nrow(cells_A)),
    y = runif(nrow(cells_A))
  )
  
  K_perm[, b] <- weighted_cross_K_field(
    cells = cells_null,
    pixels = grid_df,
    r_vals = r_vals,
    delta_a = delta_a
  )$Kab
}

K_obs$perm_mean <- rowMeans(K_perm)
K_obs$perm_lo <- apply(K_perm, 1, quantile, probs = 0.025)
K_obs$perm_hi <- apply(K_perm, 1, quantile, probs = 0.975)
K_obs$excess <- K_obs$Kab - K_obs$perm_mean

# -----------------------------
# 5. Plots
# -----------------------------

# image with A cells overlaid
p1 <- ggplot() +
  geom_raster(data = grid_df, aes(x = x, y = y, fill = B)) +
  geom_point(data = cells_A, aes(x = x, y = y), color = "black", size = 1) +
  coord_equal() +
  theme_bw() +
  labs(
    title = "Toy image: B intensity field and Type A cells",
    x = "x",
    y = "y",
    fill = "B intensity"
  )

# observed curve vs null envelope
p2 <- ggplot(K_obs, aes(x = r)) +
  geom_line(aes(y = Kab, color = "Observed"), linewidth = 1) +
  geom_line(aes(y = perm_mean, color = "Permuted mean"), linewidth = 1) +
  theme_bw() +
  labs(
    title = "Weighted K",
    x = "r",
    y = expression(hat(K)[A*B](r)),
    color = NULL
  )


print(p1)
print(p2)