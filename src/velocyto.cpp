// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;


// calculates arrows on a provided embedding (emb) given cell-cell transition probabilities (tp)
// returns arrow deltas
// [[Rcpp::export]]
arma::mat  embArrows(const arma::mat& emb, const arma::sp_mat& tp, double arrowScale=1.0, int nthreads=1) {
  arma::mat dm(emb.n_cols,emb.n_rows);

  // calculate tpb (= n), by filling up tprs first and than 1/n
  arma::sp_mat tpb(tp); // make binarized version to give equal weights to each cell in the neighborhood
  arma::vec tprs(tp.n_cols,arma::fill::zeros);
  arma::sp_mat::iterator ei=tpb.end();
  for(arma::sp_mat::iterator ci=tpb.begin(); ci!=ei; ++ci) {
    tprs[ci.col()]++;
  }
  for(arma::sp_mat::iterator ci=tpb.begin(); ci!=ei; ++ci) {
    (*ci)=1.0/tprs[ci.col()];
  }

  arma::colvec zv(emb.n_cols,arma::fill::zeros);
  arma::mat temb=trans(emb);
#pragma omp parallel for shared(dm) num_threads(nthreads)
  for(int i=0;i<emb.n_rows;i++) {
    arma::mat di(temb); // copy of the embedding matrix
    di.each_col()-=di.col(i); // coordinate difference to every cell
    di=arma::normalise(di,2,0) * arrowScale; // normalized and scaled coordinate difference
    di.col(i)=zv; // no distance to itself
    arma::vec ds=di * tp.col(i) - di * tpb.col(i); // left part is p_ij, right part 1/n
    dm.col(i)=ds;
  }
  return(dm);
}




// calculates correlation matrix between e-e_i and d_i with log10 transformation on e-e_i
// [[Rcpp::export]]
arma::mat colDeltaCorLog10(
    const arma::mat& e,
    const arma::mat& d,
    const arma::vec& waypoints,
    double pseudocount=1.0,
    int nthreads=1
) {
  arma::mat rm(e.n_cols,waypoints.n_elem);

  int waypoint=0;

#pragma omp parallel for shared(rm) num_threads(nthreads)
  for(int i=0;i<waypoints.n_elem;i++) {
    waypoint = waypoints[i]; // get the index of the waypoint in e
    arma::mat t(e); t.each_col() -= e.col(waypoint); // calculate difference
    t=log10(abs(t)+pseudocount) % sign(t); // transformation function
    rm.col(i)=cor(t,d.col(i)); // calculate correlation
    rm(waypoint,i) = 0; // set diagonal to 0

    // Rcout << "-->" << std::endl << rm(waypoint,i) << std::endl;
  }

  return(rm);
}
