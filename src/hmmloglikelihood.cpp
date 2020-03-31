#include <Rcpp.h>
using namespace Rcpp;


// This function computes the log-likelihood
// associated with a hidden Markov model. ******* USE AT YOUR OWN RISK!
//
// The hidden variables (IBD_t) are a Markov chain in {0,1},
// where the index refers to a site on the genome.
// The initial probabilities are (1 - r, r), i.e IBD_1 = 1 with probability r.
//
// Then the transition matrix is:
// A(t) =
// (a_00(t) a_01(t))
// (a_10(t) a_11(t)),
// where a_jl(t) is the probability of transitioning from state j at site t-1 to state l at site t
//
// a_01(t) = r(1 - exp(-k * rho * gendist_t))
// a_10(t) = (1 - r) * (1 - exp(-k * rho * gendist_t))
// where k is a parameter, rho is a constant assumed known e.g. 7.4×10^(-7) in M bp^{-1}
// and gendist_t is a distance between sites t and t-1 (possibly + infinity if t is the first site,
// or if the sites t-1 and t are in different chromosomes).
//
// Given IBD_t, some variables G_t^{(i)}, G_t^{(j)} follow
// P(G_t^{(i)} = g^{(i)}, G_t^{(j)} = g^{(j)}| IBD_t = 0) = f_t(g^{(i)}) f_t(g^{(j)})
// P(G_t^{(i)} = g^{(i)}, G_t^{(j)} = g^{(j)}| IBD_t = 1) = f_t(g^{(i)}) 1(g^{(i)} == g^{(j)}).
// where f_t(g) denotes the population frequency of allele g at site t
// and g takes values in 0,1,... (depending on number of different alleles possible).
// For instance if g is in (0,1), then f_t is made of (f_t(0), f_t(1)) (summing to one) and the matrix f
// should have two columns.
//
// The variables G_t^{(i)}, G_t^{(j)} are not observed; instead
// we observe Y_t^{(i)}, Y_t^{(j)}, assumed independent given  G_t^{(i)}, G_t^{(j)}, that follow
// P(Y_t^{(i)} = g^{(i)} |G_t^{(i)} = g)
// = {1 - (#{g} -1) * epsilon  if g^{(i)} = g,
//   { epsilon                 if g^{(i)} != g.
// Above, #{g} is the number of possible allele types g, i.e. if g is in {0,1}, #{g} = 2.
// The above means that there are genotyping errors occurring with probability (#{g}-1) epsilon,
// and if they occur, the observed 'g' is taken uniformly among the (#{g}-1) other types.
//
// The parameters to be estimated are k and r.
// The observations Yt^{(i)},Yt^{(j)} are provided as vectors of equal size 'ndata'.
// The matrix f should have ndata rows and max{#{g}} columns. For instance, if it has
// five columns, then each row represents the frequencies of five different types of alleles.
// If only say 3 allele types are possible for this position,
// then only the first 3 entries of this row should be non-zero.
//
// The genetic distances are provided as a vector of size ndata.
// The first entry should be the distance between position 1 and 2.
// The last entry could be anything, because it is not used; only the first ndata-1 entries are used.
//
// Finally the constants epsilon and rho can be provided, e.g. epsilon = 0.001 and rho = 7.4×10^(-7).
//
// Technically the function implements the forward algorithm, considering (IBD_t) to be the latent
// Markov chain, and G_t^{(i)}, G_t^{(j)} are integrated out for a cost of (#{g})^2 (simply via a double sum
// over all possible values of G_t^{(i)}, G_t^{(j)}).

// [[Rcpp::export]]
double loglikelihood_cpp(double k, double r,
                         const IntegerMatrix & Ys,
                         const NumericMatrix & f, const NumericVector & gendist,
                         double epsilon, double rho){
  // loglikelihood to be computed
  double loglikelihood_value = 0.;
  double l_idata, lk0, lk1, exp_, a01, a11, incr; // temporary quantities used to break down calculations
  if (r < 0 || r > 1 || k < 0){ // if r or k are not in the feasible range, return -infinity.
    return(log(0.));
  }
  // predictive distribution of latent chain given past observations
  // initially set to (1-r, r)
  NumericVector current_predictive(2);
  current_predictive(0) = 1. - r;
  current_predictive(1) = r;
  // filtering distribution of latent chain given past and present
  NumericVector current_filter(2);
  // size of data
  int ndata = Ys.nrow();
  // maximum number of alleles
  int maxnstates = f.ncol();
  int nstates;
  // loop over data to implement the forward algorithm
  for (int idata = 0; idata < ndata; idata ++){
    // compute number of different alleles at that site
    nstates = 0;
    while ((nstates < maxnstates) && (f(idata,nstates) > 1e-20)){
      nstates += 1;
    }
    // update predictive into filter, using observation 'idata'
    // likelihood of observations given IBD_t = 0
    lk0 = 0.;
    incr = 0.;
    // double sum to integrate out G_t^{(i)}, G_t^{(j)}
    for (int g = 0; g < nstates; g ++){
      for (int gprime = 0; gprime < nstates; gprime ++){
        incr = f(idata, g) * f(idata, gprime);
        if (Ys(idata,0) == g){
          incr *= (1 - (nstates - 1) * epsilon);
        } else {
          incr *= epsilon;
        }
        if (Ys(idata,1) == gprime){
          incr *= (1 - (nstates - 1) * epsilon);
        } else {
          incr *= epsilon;
        }
        lk0 += incr;
      }
    }
    // likelihood of observations given IBD_t = 1
    lk1 = 0.;
    incr = 0.;
    // this time, no need to sum over g != g' since this has zero probability given the event {IBD_t = 1}
    for (int g = 0; g < nstates; g ++){
      incr = f(idata, g);
      if (Ys(idata,0) == g){
        incr *= (1 - (nstates - 1) * epsilon);
      } else {
        incr *= epsilon;
      }
      if (Ys(idata,1) == g){
        incr *= (1 - (nstates - 1) * epsilon);
      } else {
        incr *= epsilon;
      }
      lk1 += incr;
    }
    // filtering distribution of latent_t given y_1,...,y_{t}
    // obtained by Bayes formula, given predictive and conditional likelihood
    current_filter(0) = current_predictive(0) * lk0;
    current_filter(1) = current_predictive(1) * lk1;
    // likelihood is p(y_t | y_1, ..., y_{t-1}), the normalizing constant in Bayes formula
    l_idata = current_filter(0) + current_filter(1);
    // add log-likelihood
    loglikelihood_value += log(l_idata);
    // rest is not performed at the last iteration of the for loop
    if (idata < ndata-1){
      // normalize filtering distribution
      current_filter = current_filter / l_idata;
      // next, obtain the next predictive distribution using the formula
      // P(IBD_t = 1 | past data) = P(IBD_{t-1} = 1 | past data) a_11(t) +
      //                            P(IBD_{t-1} = 0 | past data) a_01(t)
      // i.e. using filtering distribution and transition matrix
      // where a_{jl}(t) is proba to transition from j at site t-1 to l at site t
      exp_ = exp(- k * rho * gendist(idata));
      a01 = r * (1 - exp_);
      a11 = r + (1-r) * exp_;
      current_predictive(1) = current_filter(0) * a01 + current_filter(1) * a11;
      // the other probability has to be one minus the above one
      current_predictive(0) = 1 - current_predictive(1);
    }
  }
  // return log-likelihood: log p(y_1,...y_t; parameters).
  return(loglikelihood_value);
}
