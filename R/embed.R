embed_velocity_difference_waypointed <- function(
  emb,
  vel,
  corr_sigma = 0.01
) {
  # select waypoint cells
  cells <- rownames(emb)
  waypoint_cells <- sample(cells, n_waypoints)

  # prepare for correlation calculation
  em <- as.matrix(vel$current)
  ccells <- intersect(rownames(emb), colnames(em))
  em <- em[, ccells]
  emb <- emb[ccells, ]
  nd <- as.matrix(vel$deltaE[, ccells])
  cgenes <- intersect(rownames(em), rownames(nd))
  nd <- nd[cgenes, ]
  em <- em[cgenes, ]

  # calculate correlation
  # this is an adapted version of colDeltaCorLog10 with waypoints
  transfo <- function(x) (log10(abs(x) + 1) * sign(x))
  ndtransfo(nd)
  cc <- colDeltaCorLog10(
    em,
    nd,
    match(colnames(nd2), colnames(em)) - 1
  ) %>% t()
  colnames(cc) <- colnames(em)
  rownames(cc) <- colnames(nd)

  # only look among local waypoints to determine arrow
  knn_result <- RANN::nn2(emb[waypoint_cells,], emb, k = k, eps = 0)
  # knn_result <- FNN::get.knnx(emb[waypoint_cells,], emb, k = k)
  i <- rep(seq_len(nrow(knn_result$nn.idx)), each = k)
  j <- as.integer(t(knn_result$nn.idx))
  knn <- Matrix::sparseMatrix(i, j, dims = c(nrow(emb), nrow(emb)))
  colnames(knn) <- c(waypoint_cells, setdiff(rownames(emb), waypoint_cells))
  rownames(knn) <- rownames(emb)
  knn <- as.matrix(Matrix::t(knn))[rownames(emb), rownames(emb)]

  # calculate transition probabilities (from col to row)
  tp <- exp(cc / corr_sigma) * (cc != 0) * knn
  tp <- Matrix::t(Matrix::t(tp) / Matrix::colSums(tp))
  tp <- as(tp,'dgCMatrix')

  # arrow estimates for each cell
  arsd <- data.frame(t(embArrows(emb, tp, 1, nthreads = 1)))
  rownames(arsd) <- rownames(emb)
  colnames(arsd) <- paste0(colnames(emb), "_velocity")

  if (any(is.na(arsd))) stop()

  arsd
}





#' Calculate the velocity of cells within an existing embedding
#'
#' This is described in page 7-8 Supplementary note
#' To speed things up, we use waypoints cells, and not the full-blown cc calculation as in the original version
#'
#' @param emb The embedding, a matrix with cells in rows and dimensions in columns
#' @param vel The relative velocity output
#' @param n_waypoints Number of waypoints
#' @param corr_sigma Sigma parameter used to translate velocity-(expression delta) correlation into a transition probability. Higher sigma values mean
embed_velocity_difference_waypointed <- function(
  emb,
  vel,
  n_waypoints = 50,
  k = ceiling(n_waypoints / 10),
  corr_sigma = 10
) {
  # select waypoint cells
  cells <- rownames(emb)
  waypoint_cells <- sample(cells, n_waypoints)

  # prepare for correlation calculation
  em <- as.matrix(vel$current)
  ccells <- intersect(rownames(emb), colnames(em))
  em <- em[, ccells]
  emb <- emb[ccells, ]
  nd <- as.matrix(vel$deltaE[, waypoint_cells])
  cgenes <- intersect(rownames(em), rownames(nd))
  nd <- nd[cgenes, ]
  em <- em[cgenes, ]

  # calculate correlation
  # this is an adapted version of colDeltaCorLog10 with waypoints
  transfo <- function(x) (log10(abs(x) + 1) * sign(x))
  nd2 <- transfo(nd)
  cc2 <- colDeltaCorLog10(
    em,
    nd2,
    match(colnames(nd2), colnames(em)) - 1
  ) %>% t()
  colnames(cc2) <- colnames(em)
  rownames(cc2) <- colnames(nd2)

  # add empty rows for non-waypoint cells, necessary be used by embArrows
  cc <- Matrix::rbind2(
    Matrix(cc2, sparse = T),
    Matrix(0, ncol = ncol(cc2), nrow = nrow(emb) - nrow(cc2), sparse = TRUE, dimnames = list(setdiff(rownames(emb), rownames(cc2)), colnames(cc2)))
  )
  cc <- as.matrix(cc)[rownames(emb), rownames(emb)]

  # only look among local waypoints to determine arrow
  knn_result <- RANN::nn2(emb[waypoint_cells,], emb, k = k, eps = 0)
  # knn_result <- FNN::get.knnx(emb[waypoint_cells,], emb, k = k)
  i <- rep(seq_len(nrow(knn_result$nn.idx)), each = k)
  j <- as.integer(t(knn_result$nn.idx))
  knn <- Matrix::sparseMatrix(i, j, dims = c(nrow(emb), nrow(emb)))
  colnames(knn) <- c(waypoint_cells, setdiff(rownames(emb), waypoint_cells))
  rownames(knn) <- rownames(emb)
  knn <- as.matrix(Matrix::t(knn))[rownames(emb), rownames(emb)]

  # calculate transition probabilities (from col to row)
  tp <- exp(cc / corr_sigma) * (cc != 0) * knn
  tp <- Matrix::t(Matrix::t(tp) / Matrix::colSums(tp))
  tp <- as(tp,'dgCMatrix')

  # arrow estimates for each cell
  arsd <- data.frame(t(embArrows(emb, tp, 1, nthreads = 1)))
  rownames(arsd) <- rownames(emb)
  colnames(arsd) <- paste0(colnames(emb), "_velocity")

  if (any(is.na(arsd))) stop()

  arsd
}




# use knn to draw arrows
embed_velocity_knn <- function(
  emb,
  vel,
  k = 10
) {
  knn_result2 <- RANN::nn2(t(log(current+1)), t(log(current+1)), k = k)
  knn_result <- RANN::nn2(t(log(current+1)), t(log(projected+1)), k = k)

  emb_tmp <- apply(knn_result$nn.idx, 1, function(nn) {
    colMeans(emb[nn,])
  }) %>% t() %>% magrittr::set_colnames(paste0(colnames(emb), "_next"))
  emb_notnext <- apply(knn_result2$nn.idx, 1, function(nn) {
    colMeans(emb[nn,])
  }) %>% t() %>% magrittr::set_colnames(paste0(colnames(emb), "_next"))

  emb_next <- emb + (emb_tmp - emb) - (emb_notnext - emb)
  arsd <- emb - emb_next
  arsd
}
