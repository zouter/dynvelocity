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
