#' Calculates a smoothed velocity at particular points within an embedding
#'
#' @examples
#' points <- matrix(c(0, 0, 1, 1, 0, 1, 1, 0), nrow = 4)
smooth_velocity <- function(
  arsd,
  emb,
  points,
  sd = max(apply(emb, 2, max) - apply(emb, 2, min)) / 20,
  scale_arrows = 1
) {
  if(ncol(points) != ncol(emb)) stop("points and emb should have the same number of columns (= dimensions)")
  if(nrow(arsd) != nrow(emb)) stop("arsd and emb should have the same number of rows (= cells)")

  distances <- pdist::pdist(emb, points) %>% as.matrix()
  weights <- dnorm(distances, sd = sd)
  weights[weights < 0.01] <- 0.00000001

  density <- colSums(weights)
  weights <- t(t(weights) / colSums(weights))

  arsd_points <- t(t(arsd) %*% weights)

  colnames(arsd_points) <- colnames(arsd)

  arsd_points <- arsd_points / max(sqrt(rowSums(arsd_points**2))) * scale_arrows

  list(
    arsd_points = arsd_points,
    density = density
  )
}



create_grid_points <- function(
  emb,
  grid_n = 10
) {

  points
}



#' Calculates a smoothened velocity as a grid within an embedding
smooth_velocity_grid <- function(
  arsd,
  emb,
  grid_n = 10,
  relative_density_cutoff = 0.2
) {
  # create grid points
  ranges <- apply(emb, 2, range)
  positions <- ranges %>% apply(2, function(x) {seq(x[1], x[2], length.out = grid_n)}) %>% as.data.frame()
  points <- expand.grid(positions)

  # smooth to grid points
  difference <- max(diff(positions[, 1]))
  smoothened <- smooth_velocity(arsd, emb, points, scale_arrows = difference)

  # create arrows from points, arsd of points, and local density
  arrows <- bind_cols(
    as_tibble(points),
    as_tibble(smoothened$arsd_points),
    tibble(density = smoothened$density)
  ) %>%
    mutate(
      length = sqrt(rowSums(smoothened$arsd_points**2)),
      angle = map2_dbl(smoothened$arsd_points[,2], smoothened$arsd_points[,1], atan2) / pi * 180
    )

  # filter based on local density
  arrows <- arrows %>%
    filter(
      density > quantile(density, 0.9) * relative_density_cutoff
    )

  arrows
}



calculate_arrows <- function(

) {

}