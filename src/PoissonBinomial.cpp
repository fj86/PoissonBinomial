// [[Rcpp::interfaces(r, cpp)]]
//'

#include <fftw3.h>
#define STRICT_R_HEADERS
#include <Rcpp.h>
using namespace Rcpp;

// used to make fft computations more readable
#define REAL 0
#define IMAG 1
//#define PI2 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651e+00

/******************************************************************/
/**   Helper Functions                                           **/
/******************************************************************/

// helper function for normalisation of PMFs (i.e. ensure that sum = 1)
void norm_dpb(NumericVector &pmf){
  // sums of PMF
  double new_sum = sum(pmf), old_sum = 0, older_sum = 0, oldest_sum = 0;
  //Rcout << ((new_sum < 1)?"l ":((new_sum == 1)?"e ":"g "));
  while(new_sum != 1){
    oldest_sum = older_sum;
    older_sum = old_sum;
    old_sum = new_sum;
    NumericVector old_pmf = pmf;
    pmf = pmf / new_sum;
    new_sum = sum(pmf);
    //Rcout << ((new_sum < 1)?"l ":((new_sum == 1)?"e ":"g "));
    if(new_sum >= 1 || new_sum == old_sum || new_sum == older_sum || new_sum == oldest_sum) break;
    if(new_sum < 1 && new_sum <= old_sum){
      pmf = old_pmf;
      break;
    }
  }
  //Rcout << "\n";
}

// "generic" function for computing some of the PMFs
NumericVector dpb_generic(const IntegerVector obs, const NumericVector cdf){
  // maximum observed value
  const int max_q = obs.length() ? max(obs) : cdf.length() - 1;
  // results vector
  NumericVector results(max_q + 1);
  
  // compute masses
  results[0] = cdf[0];
  for(int i = 1; i <= max_q; i++)
    results[i] = cdf[i] - cdf[i - 1];
  
  // return final results
  if(obs.length()) return results[obs]; else return results;
}

// "generic" function for computing some of the CDFs
NumericVector ppb_generic(const IntegerVector obs, const NumericVector pmf, bool lower_tail = true){
  // distribution size
  const int size = pmf.length();
  // maximum observed value
  const int max_q = obs.length() ? max(obs) : size - 1;
  // results vector
  NumericVector results = NumericVector(std::min<int>(max_q + 1, size));
  
  // compute cumulative probabilities
  if(lower_tail){
    results[0] = pmf[0];
    for(int i = 1; i <= max_q; i++)
      results[i] = pmf[i] + results[i - 1];
  }else{
    const int min_q = obs.length() ? min(obs) : 0;
    const int len = pmf.length() - 1;
    for(int i = len; i > min_q; i--){
      if(i > max_q) results[max_q] += pmf[i];
      else results[i - 1] = pmf[i] + results[i];
    }
  }
  
  // "correct" numerically too large results
  results[results > 1] = 1;
  
  // return final results
  if(obs.length()) return results[obs]; else return results;
}

IntegerVector order(NumericVector x, bool decreasing = false){
  NumericVector uni = unique(x).sort();
  if(decreasing) uni = NumericVector(rev(uni));
  IntegerVector order(x.length());
  int k = 0;
  for(int i = 0; i < uni.length(); i++){
    for(int j = 0; j < x.length(); j++){
      if(uni[i] == x[j]) order[k++] = j;
    }
  }
  return order;
}

// [[Rcpp::export]]
int vectorGCD(const IntegerVector x){
  // input size
  const int size = x.length();
  
  if(size == 0) return 0;
  
  // make all values positive
  IntegerVector y;
  y = abs(x);
  
  // minimum of 'x' (add 1 to make sure that it is greater than the first value)
  int xmin = y[0] + 1;
  
  // search for minimum, one and zero; return it, if found
  for(int i = 0; i < size; i++){
    if(xmin > y[i]){
      xmin = y[i];
      if(xmin <= 1) return xmin;
    }
  }
  
  int a, b, r, i = 0, gcd = xmin;
  
  while(gcd > 1 && i < size){
    a = std::max<int>(gcd, x[i]);
    b = std::min<int>(gcd, x[i]);
    
    while(b != 0){
      r = a % b;
      a = b;
      b = r;
    }
    gcd = a;
    
    i++;
  }
  
  return gcd;
}

/******************************************************************/
/**   Functions for "ordinary" Poisson binomial distribution     **/
/******************************************************************/

// PMFs
/*NumericVector dpb_conv(IntegerVector obs, NumericVector probs);
NumericVector dpb_dc(IntegerVector obs, NumericVector probs);
NumericVector dpb_dftcf(IntegerVector obs, NumericVector probs);
NumericVector dpb_rf(IntegerVector obs, NumericVector probs);
NumericVector dpb_mean(IntegerVector obs, NumericVector probs);
NumericVector dpb_gmba(IntegerVector obs, NumericVector probs, bool anti = false);
NumericVector dpb_pa(IntegerVector obs, NumericVector probs);
NumericVector dpb_na(IntegerVector obs, NumericVector probs, bool refined = true);

// CDFs
NumericVector ppb_conv(IntegerVector obs, NumericVector probs, bool lower_tail = true);
NumericVector ppb_dc(IntegerVector obs, NumericVector probs, bool lower_tail = true);
NumericVector ppb_dftcf(IntegerVector obs, NumericVector probs, bool lower_tail = true);
NumericVector ppb_rf(IntegerVector obs, NumericVector probs, bool lower_tail = true);
NumericVector ppb_mean(IntegerVector obs, NumericVector probs, bool lower_tail = true);
NumericVector ppb_gmba(IntegerVector obs, NumericVector probs, bool anti = false, bool lower_tail = true);
NumericVector ppb_pa(IntegerVector obs, NumericVector probs, bool lower_tail = true);
NumericVector ppb_na(IntegerVector obs, NumericVector probs, bool refined = true, bool lower_tail = true);*/

// Direct Convolution
// [[Rcpp::export]]
NumericVector dpb_conv(const IntegerVector obs, const NumericVector probs){
  // number of input probabilities
  const int size = probs.length();
  
  // results vector
  NumericVector results(size + 1);
  results[0] = 1 - probs[0];
  results[1] = probs[0];
  
  for(int i = 1; i < size; i++){
    checkUserInterrupt();
    if(probs[i]){
      for(int j = i; j >= 0; j--){
        if(results[j]){
          results[j + 1] += results[j] * probs[i];
          results[j] *= 1 - probs[i];
        }
      }
    }
  }
  // make sure that probability masses sum up to 1
  norm_dpb(results);
  
  // return final results
  if(obs.length()) return results[obs]; else return results;
}

// [[Rcpp::export]]
NumericVector ppb_conv(const IntegerVector obs, const NumericVector probs, const bool lower_tail = true){
  // number of input probabilities
  const int size = probs.length();
  
  // highest observed value
  const int max_q = obs.length() ? max(obs) : size;
  
  // probability masses
  const NumericVector pmf = dpb_conv(IntegerVector(), probs);
  
  // compute CDF
  NumericVector results = ppb_generic(obs, pmf, lower_tail);
  
  // ensure that (for lower tail) sum = 1, if last value = n (the highest observable value)
  if(obs.length()){
    if(max_q == size) results[obs == max_q] = (double)lower_tail;
  }else results[size] = (double)lower_tail;
  
  // return final results
  return results;
}

// Divide & Conquer FFT (DC-FFT)
NumericVector fft_probs(const NumericVector probsA, const NumericVector probsB){
  // sizes of input vectors and the result
  const int sizeA = probsA.length();
  const int sizeB = probsB.length();
  const int sizeResult = sizeA + sizeB - 1;
  
  // results vector
  double *result_vec = new double[sizeResult];
  
  // allocate memory for FFTs of the probs and the convolution result
  fftw_complex *probsA_fft, *probsB_fft, *result_fft;
  
  // 0-padding of probsA vector and perform FFT of it
  NumericVector padded_probsA(sizeResult);
  padded_probsA[Range(0, sizeA - 1)] = probsA;
  probsA_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sizeResult);
  fftw_plan planA = fftw_plan_dft_r2c_1d(sizeResult, padded_probsA.begin(), probsA_fft, FFTW_ESTIMATE);
  fftw_execute(planA);
  fftw_destroy_plan(planA);
  
  // 0-padding of probsB vector and perform FFT of it
  NumericVector padded_probsB(sizeResult);
  padded_probsB[Range(0, sizeB - 1)] = probsB;
  probsB_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sizeResult);
  fftw_plan planB = fftw_plan_dft_r2c_1d(sizeResult, padded_probsB.begin(), probsB_fft, FFTW_ESTIMATE);
  fftw_execute(planB);
  fftw_destroy_plan(planB);
  
  // convolute by complex multiplication of the transformed input probs
  result_fft = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sizeResult);
  for(int i = 0; i < sizeResult; i++){
    result_fft[i][REAL] = (probsA_fft[i][REAL]*probsB_fft[i][REAL] - probsA_fft[i][IMAG]*probsB_fft[i][IMAG])/sizeResult;
    result_fft[i][IMAG] = (probsA_fft[i][REAL]*probsB_fft[i][IMAG] + probsA_fft[i][IMAG]*probsB_fft[i][REAL])/sizeResult;
  }
  
  // inverse tranformation of the above multiplications
  //fftw_plan planResult = fftw_plan_dft_c2r_1d(sizeResult, result_fft, result.begin(), FFTW_ESTIMATE);
  fftw_plan planResult = fftw_plan_dft_c2r_1d(sizeResult, result_fft, result_vec, FFTW_ESTIMATE);
  fftw_execute(planResult);
  fftw_destroy_plan(planResult);
  
  // garbage collection
  fftw_free(probsA_fft);
  fftw_free(probsB_fft);
  fftw_free(result_fft);
  
  // return final results
  NumericVector result(sizeResult);
  for(int i = 0; i < sizeResult; i++) result[i] = result_vec[i];
  delete[] result_vec;
  return result;
}

// [[Rcpp::export]]
NumericVector dpb_dc(const IntegerVector obs, const NumericVector probs){//, const int splits = -1){
  // number of probabilities of success
  const int size = probs.length();
  
  // automatically determine number of splits, if size is above 1950
  //int num_splits = splits < 0 ? std::max<int>(0, (int)std::ceil(std::log(size / 1950) / std::log(2.0))) : splits;
  int num_splits = size > 1950 ? (int)std::ceil(std::log(size / 1950) / std::log(2.0)) : 0;
  // direct convolution is sufficient in case of 0 splits
  if(num_splits == 0) return dpb_conv(obs, probs);
  // number of groups
  int num_groups = std::pow(2, num_splits);
  // reduce number of splits and groups if too large
  while(num_splits > 0 && num_groups > size){
    num_splits -= 1;
    num_groups /= 2;
  }
  // direct convolution is sufficient, if no splits are necessary
  if(num_splits == 0) return dpb_conv(obs, probs);
  
  // range variables
  int start, end;
  
  // compute group sizes with minimum size disparity
  IntegerVector group_sizes(num_groups, size / num_groups);
  const int remainder = size % num_groups;
  for(int i = 0; i < remainder; i++) group_sizes[i]++;
  
  // compute first and last indices of the groups
  IntegerVector starts(num_groups), ends(num_groups);
  starts[0] = 0;
  ends[0] = group_sizes[0] - 1;
  for(int i = 1; i < num_groups; i++){
    starts[i] = starts[i - 1] + group_sizes[i - 1];
    ends[i] = ends[i - 1] + group_sizes[i];
  }
  
  // results vector; direct allocation will increase size of each group by 1
  NumericVector results(size + num_groups);
  
  // compute direct convolutions for each group
  for(int i = 0; i < num_groups; i++){
    checkUserInterrupt();
    // compute new starting and ending indices, because groups grow by 1
    start = starts[i] + i;
    end = ends[i] + i + 1;
    
    // target range
    Range target(start, end);
    
    // direct convolution
    results[target] = dpb_conv(IntegerVector(), probs[Range(starts[i], ends[i])]);
    
    // update starting and ending indices
    starts[i] = start;
    ends[i] = end;
  }
  
  int num_groups_reduced = num_groups / 2;
  while(num_splits > 0){
    for(int i = 0; i < num_groups_reduced; i++){
      checkUserInterrupt();
      // compute new starting and ending indices, because group sizes are
      // reduced by 1, due to FFT convolution
      start = starts[2*i] - i;
      end = ends[2*i + 1] - i - 1;
      
      //convolution
      results[Range(start, end)] = fft_probs(results[Range(starts[2*i], ends[2*i])], results[Range(starts[2*i + 1], ends[2*i + 1])]);
      
      // update starting and ending indices
      starts[i] = start;
      ends[i] = end;
    }
    num_groups_reduced /= 2;
    num_splits -= 1;
  }
  
  // select final results
  results = NumericVector(results[Range(0, size)]);
  
  // "correct" numerically false (and thus useless) results
  results[results < 5.55e-17] = 0;
  results[results > 1] = 1;
  
  // make sure that probability masses sum up to 1
  norm_dpb(results);
  
  // return final results
  if(obs.length()) return results[obs]; else return results;
}

// [[Rcpp::export]]
NumericVector ppb_dc(const IntegerVector obs, const NumericVector probs, const bool lower_tail = true){
  // number of input probabilities
  const int size = probs.length();
  
  // highest observed value
  const int max_q = obs.length() ? max(obs) : size;
  
  // probability masses
  const NumericVector pmf = dpb_dc(IntegerVector(), probs);
  
  // compute CDF
  NumericVector results = ppb_generic(obs, pmf, lower_tail);
  
  // ensure that (for lower tail) sum = 1, if last value = n (the highest observable value)
  if(obs.length()){
    if(max_q == size) results[obs == max_q] = (double)lower_tail;
  }else results[size] = (double)lower_tail;
  
  // return final results
  return results;
}

// Discrete Fourier Transformation of Characteristic Function (DFT-CF)
// [[Rcpp::export]]
NumericVector dpb_dftcf(const IntegerVector obs, const NumericVector probs){
  // number of probabilities of success
  const int sizeIn = probs.length();
  // number of distribution
  const int sizeOut = sizeIn + 1;
  
  // "initialize" DFT input vector
  fftw_complex *input_fft;
  input_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sizeOut);
  input_fft[0][REAL] = 1.0;
  input_fft[0][IMAG] = 0.0;
  
  // initialize complex numbers for "C" and "C to the power of i"
  const std::complex<double> C = exp(std::complex<double>(0.0, 2.0) * M_PI / ((double)sizeOut));
  std::complex<double> C_power = 1.0;
  
  // compute closed-form expression of Hernandez and Williams
  const int mid = sizeIn / 2 + 1;
  for(int i = 1; i <= mid; i++){
    checkUserInterrupt();
    
    C_power *= C;
    
    std::complex<double> product = 1.0;
    for(int j = 0; j < sizeIn; j++) product *= 1.0 + (C_power - 1.0) * probs[j];
    
    input_fft[i][REAL] = product.real();
    input_fft[i][IMAG] = product.imag();
    input_fft[sizeOut - i][REAL] = product.real();
    input_fft[sizeOut - i][IMAG] = -product.imag();
  }
  
  // vector of DFT results
  fftw_complex *result_fft;
  result_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sizeOut);
  
  // perform DFT
  fftw_plan planDFT;
  planDFT = fftw_plan_dft_1d(sizeOut, input_fft, result_fft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(planDFT);
  
  // gather results
  NumericVector results(sizeOut);
  for(int i = 0; i < sizeOut; i++) results[i] = result_fft[i][REAL] / sizeOut;
  
  // garbage collection
  fftw_destroy_plan(planDFT);
  fftw_free(input_fft);
  fftw_free(result_fft);
  
  // "correct" numerically false (and thus useless) results
  results[results < 2.22e-16] = 0;
  results[results > 1] = 1;
  
  // make sure that probability masses sum up to 1
  norm_dpb(results);
  
  // return final results
  if(obs.length()) return results[obs]; else return results;
}

// [[Rcpp::export]]
NumericVector ppb_dftcf(const IntegerVector obs, const NumericVector probs, const bool lower_tail = true){
  // number of input probabilities
  const int size = probs.length();
  
  // highest observed value
  const int max_q = obs.length() ? max(obs) : size;
  
  // probability masses
  const NumericVector pmf = dpb_dftcf(IntegerVector(), probs);
  
  // compute CDF
  NumericVector results = ppb_generic(obs, pmf, lower_tail);
  
  // ensure that (for lower tail) sum = 1, if last value = n (the highest observable value)
  if(obs.length()){
    if(max_q == size) results[obs == max_q] = (double)lower_tail;
  }else results[size] = (double)lower_tail;
  
  // return final results
  return results;
}

// Recursive Formula
// [[Rcpp::export]]
NumericVector dpb_rf(const IntegerVector obs, const NumericVector probs){
  // number of input probabilities
  const int size = probs.length();
  
  NumericMatrix dist(size + 1, 2);
  NumericVector results(size + 1);
  int col_new = 0, col_old = 1;
  
  dist(0, col_new) = 1.0;
  dist(1, col_new) = 1 - probs[0];
  for(int j = 1; j < size; j++) dist(j + 1, col_new) = (1 - probs[j]) * dist(j, col_new);
  results[0] = dist(size, col_new);
  
  for(int i = 1; i <= size; i++){
    checkUserInterrupt();
    col_new -= std::pow(-1, i);
    col_old += std::pow(-1, i);
    
    for(int j = 0; j <= i - 1; j++)
      dist(j, col_new) = 0;
    
    for(int j = i - 1; j < size; j++){
      dist(j + 1, col_new) = (1 - probs[j]) * dist(j, col_new) + probs[j] * dist(j, col_old);
    }
    
    results[i] = dist(size, col_new);
  }
  
  // make sure that probability masses sum up to 1
  norm_dpb(results);
  
  // return final results
  if(obs.length()) return results[obs]; else return results;
}

// [[Rcpp::export]]
NumericVector ppb_rf(const IntegerVector obs, const NumericVector probs, const bool lower_tail = true){
  // number of input probabilities
  int size = probs.length();
  
  // highest observed value
  int max_q = obs.length() ? max(obs) : size;
  
  // probability masses
  const NumericVector pmf = dpb_rf(IntegerVector(), probs);
  
  // compute CDF
  NumericVector results = ppb_generic(obs, pmf, lower_tail);
  
  // make sure that largest observation has probability of 1 (or 0, depending on lower_tail)
  if(obs.length()){
    if(max_q == size) results[obs == max_q] = (double)lower_tail;
  }else results[size] = (double)lower_tail;
  
  // return final results
  return results;
}

// Arithmetic Mean Binomial Approximation
// [[Rcpp::export]]
NumericVector dpb_mean(IntegerVector obs, const NumericVector probs){
  // number of input probabilities
  const int size = probs.length();
  
  // mean of probabilities is the approximate binomial probability
  const double bin_prob = mean(probs);
  
  // compute probability masses and return
  if(obs.length() == 0)
    return dbinom(IntegerVector(Range(0, size)), (double)size, bin_prob);
  else return dbinom(obs, (double)size, bin_prob);
}

// [[Rcpp::export]]
NumericVector ppb_mean(const IntegerVector obs, const NumericVector probs, const bool lower_tail = true){
  // number of input probabilities
  const int size = probs.length();
  
  // mean of probabilities is the approximate binomial probability
  const double bin_prob = mean(probs);
  
  // compute cumulative probabilities and return
  if(obs.length() == 0)
    return pbinom(IntegerVector(Range(0, size)), (double)size, bin_prob, lower_tail);
  else return pbinom(obs, (double)size, bin_prob, lower_tail);
}

// Geometric Mean Binomial Approximations
// [[Rcpp::export]]
NumericVector dpb_gmba(const IntegerVector obs, const NumericVector probs, const bool anti = false){
  // number of probabilities of success
  const int size = probs.length();
  
  // logarithms of 'probs' (sums of logarithms are numerically more stable than
  // products of probabilities, especially when the probabilities are small)
  NumericVector logs;
  double bin_prob;
  
  if(anti){
    logs = NumericVector(log(1 - probs));
    bin_prob = 1 - std::exp(mean(logs));
  }else{
    logs = NumericVector(log(probs));
    bin_prob = std::exp(mean(logs));
  }
  
  // compute probability masses and return
  if(obs.length() == 0)
    return dbinom(IntegerVector(Range(0, size)), (double)size, bin_prob);
  else return dbinom(obs, (double)size, bin_prob);
}

// [[Rcpp::export]]
NumericVector ppb_gmba(const IntegerVector obs, const NumericVector probs, const bool anti = false, const bool lower_tail = true){
  // number of probabilities of success
  const int size = probs.length();
  
  // logarithms of 'probs' (sums of logarithms are numerically more stable than
  // products of probabilities, especially when the probabilities are small)
  NumericVector logs;
  double bin_prob;
  
  if(anti){
    logs = NumericVector(log(1 - probs));
    bin_prob = 1 - std::exp(mean(logs));
  }else{
    logs = NumericVector(log(probs));
    bin_prob = std::exp(mean(logs));
  }
  
  // compute cumulative probabilities and return
  if(obs.length() == 0)
    return pbinom(IntegerVector(Range(0, size)), (double)size, bin_prob, lower_tail);
  else return pbinom(obs, (double)size, bin_prob, lower_tail);
}

// Poisson Approximation
// [[Rcpp::export]]
NumericVector dpb_pa(const IntegerVector obs, const NumericVector probs){
  // number of probabilities of success
  const int size = probs.length();
  // sum of probability is the expectation of the Poisson approximation
  const double lambda = sum(probs);
  
  // compute probability masses
  NumericVector results;
  if(obs.length() == 0){
    results = dpois(IntegerVector(Range(0, size)), lambda);
    results[size] += R::ppois(size, lambda, false, false);
  }else{
    results = dpois(obs, lambda);
    for(int i = 0; i < obs.length(); i++)
      if(obs[i] == size) results[i] += R::ppois(size, lambda, false, false);
  }
  
  // return final results
  return results;
}

// [[Rcpp::export]]
NumericVector ppb_pa(const IntegerVector obs, const NumericVector probs, bool lower_tail = true){
  // sum of probability is the expectation of the Poisson approximation
  const double lambda = sum(probs);
  
  // compute cumulative probabilities
  const IntegerVector observed = obs.length() ? obs : IntegerVector(Range(0, probs.length()));
  
  NumericVector results = ppois(observed, lambda, lower_tail);
  
  // make sure that largest possible observation has probability of 1 (or 0, depending on lower_tail)
  results[observed == probs.length()] = (double)lower_tail;
  
  // return final results
  return results;
}

// [[Rcpp::export]]
NumericVector ppb_na(const IntegerVector obs, const NumericVector probs, const bool refined = true, const bool lower_tail = true){
  // number of probabilities of success
  const int size = probs.length();
  // highest observed value
  const int max_q = obs.length() ? max(obs) : size;
  // mu
  const double mu = sum(probs);
  // p * q
  const NumericVector pq = probs * (1 - probs);
  // sigma
  const double sigma = std::sqrt(sum(pq));
  // standardized observations with continuity correction
  NumericVector obs_std;
  if(obs.length() == 0) obs_std = (NumericVector(IntegerVector(Range(0, size))) + 0.5 - mu)/sigma;
  else obs_std = (NumericVector(obs) + 0.5 - mu)/sigma;
  // vector to store results
  NumericVector results = Rcpp::pnorm(obs_std, 0.0, 1.0, lower_tail);
  // cumulative probabilities
  if(refined){
    // gamma
    const double gamma = sum(pq * (1 - 2 * probs));
    // probabilities
    if(lower_tail)
      results += gamma/(6 * std::pow(sigma, 3.0)) * (1 - pow(obs_std, 2.0)) * dnorm(obs_std);
    else results += -gamma/(6 * std::pow(sigma, 3.0)) * (1 - pow(obs_std, 2.0)) * dnorm(obs_std);
  }
  // make sure that all probabilities do not exceed 1 and are at least 0
  results[results < 0] = 0;
  results[results > 1] = 1;
  
  // make sure largest possible value has cumulative probability of 1 (lower tail) or 0 (upper tail)
  if(obs.length()){
    if(max_q >= size) results[obs >= max_q] = (double)lower_tail;
  }else results[size] = (double)lower_tail;
  
  // return final results
  return results;
}

// Normal Approximations (NA, RNA)
// [[Rcpp::export]]
NumericVector dpb_na(const IntegerVector obs, const NumericVector probs, const bool refined = true){
  // number of probabilities of success
  const int size = probs.length();
  // highest observed value
  const int max_q = obs.length() ? max(obs) : size;
  // rounded down expectation + 0.5 (continuity correction)
  const int mid = (int)floor(sum(probs) + 0.5);
  
  // cumulative probabilities
  const NumericVector cdf_lower = ppb_na(IntegerVector(Range(0, std::min<int>(mid, max_q))), probs, refined, true);
  const NumericVector cdf_upper = ppb_na(IntegerVector(Range(std::min<int>(mid, max_q), max_q)), probs, refined, false);
  
  // vector to store results
  NumericVector results(max_q + 1);
  
  // compute probability masses
  results[0] = cdf_lower[0];
  for(int i = 1; i <= max_q; i++){
    if(i <= mid) results[i] = cdf_lower[i] - cdf_lower[i - 1]; else results[i] = cdf_upper[i - 1 - mid] - cdf_upper[i - mid];
  }
  
  // compute and return results
  if(obs.length()) return results[obs]; else return results;
}

// Bernoulli Random Number Generator
// [[Rcpp::export]]
IntegerVector rpb_bernoulli(const int n, const NumericVector probs){
  // number of probabilities of success
  const int size = probs.length();
  
  // vector to store results
  NumericVector results(n);
  
  // generate random numbers
  for(int i = 0; i < size; i++) 
    for(int j = 0; j < n; j++)
      results[j] += R::rbinom(1.0, probs[i]);
  
  // return results
  return IntegerVector(results);
}


/******************************************************************/
/**   Functions for generalized Poisson binomial distribution    **/
/******************************************************************/

// PMFs
/*NumericVector dgpb_conv(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q);
NumericVector dgpb_dc(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q);
NumericVector dgpb_dftcf(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q);
NumericVector dgpb_na(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q, bool refined = true);

// CDFs
NumericVector pgpb_conv(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q, bool lower_tail = true);
NumericVector pgpb_dc(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q, bool lower_tail = true);
NumericVector pgpb_dftcf(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q, bool lower_tail = true);
NumericVector pgpb_na(IntegerVector obs, NumericVector probs, NumericVector val_p, NumericVector val_q, bool refined = true, bool lower_tail = true);*/


// Generalized Direct Convolution (G-DC)
// [[Rcpp::export]]
NumericVector dgpb_conv(const IntegerVector obs, const NumericVector probs, const IntegerVector val_p, const IntegerVector val_q){
  // number of probabilities of success
  const int sizeIn = probs.length();
  // determine pairwise minimum and maximum
  const IntegerVector v = pmin(val_p, val_q);
  const IntegerVector u = pmax(val_p, val_q);
  // compute differences
  IntegerVector d = u - v;
  // final output size
  const int sizeOut = sum(d) + 1;
  // greatest common divisor of the differences
  const int gcd = vectorGCD(d[d > 0]);
  // rescale differences according to gcd
  if(gcd > 1) d = d / gcd;
  // theoretical rescaled maximum
  const int max_rescaled = sum(d);
  // output size
  const int sizeOut_rescaled = max_rescaled + 1;
  
  // if val_p[i] was not the larger one, the respective probs[i] has to be 'flipped'
  // furthermore: if difference is 0 (i.e. u[i] equals v[i]), a non-zero outcome is impossible
  NumericVector probs_flipped(sizeIn);
  for(int i = 0; i < sizeIn; i++){
    if(!d[i]) probs_flipped[i] = 0.0;
    else if(val_p[i] < u[i]) probs_flipped[i] = 1 - probs[i];
      else probs_flipped[i] = probs[i];
  }
  
  // results vectors
  NumericVector results(sizeOut);
  NumericVector results_rescaled;
  
  // if maximum rescaled difference equals 1, we have an ordinary poisson binomial distribution
  if(max(d) == 1){
    // compute ordinary distribution
    results_rescaled = dpb_conv(IntegerVector(), probs_flipped[d > 0]);
  }else{
    results_rescaled = NumericVector(sizeOut_rescaled);
    // initialize result of first convolution step
    results_rescaled[0] = 1.0;
    // ending position of last computed iteration
    int end = 0;
    // perform convolution
    for(int i = 0; i < sizeIn; i++){
      checkUserInterrupt();
      if(probs_flipped[i]){
        for(int j = end; j >= 0; j--){
          if(d[i] && results_rescaled[j] && probs_flipped[i]){
            results_rescaled[j + d[i]] += results_rescaled[j] * probs_flipped[i];
            results_rescaled[j] *= 1 - probs_flipped[i];
          }
        }
        
        // update ending position
        end += d[i];
      }
    }
  }
  // "correct" numerically false (and thus useless) results
  results_rescaled[results_rescaled > 1] = 1;
  
  // make sure that probability masses sum up to 1
  norm_dpb(results_rescaled);
  
  // map results to generalized distribution (scale-back)
  for(int i = 0; i < sizeOut_rescaled; i++)
    results[i * gcd] = results_rescaled[i];
  
  // return final results
  if(obs.length()) return results[obs - sum(v)]; else return results;
}

// [[Rcpp::export]]
NumericVector pgpb_conv(const IntegerVector obs, const NumericVector probs, const IntegerVector val_p, const IntegerVector val_q, bool lower_tail = true){
  // theoretical minimum
  const int min_v = sum(pmin(val_p, val_q));
  // theoretical maximum
  const int max_v = sum(pmax(val_p, val_q));
  // maximum observed value
  const int max_q = obs.length() ? max(obs) : max_v;
  
  // probability masses
  const NumericVector pmf = dgpb_conv(IntegerVector(), probs, val_p, val_q);
  
  // compute CDF
  NumericVector results = ppb_generic(obs - min_v, pmf, lower_tail);
  
  // ensure that sum = 1 (or 0), if last value equals the highest observable value
  if(obs.length()){
    if(max_q == max_v) results[obs == max_q] = (double)lower_tail;
  }else results[max_v - min_v] = (double)lower_tail;
  
  // return final results
  return results;
}

// Generalized Divide & Conquer FFT Tree Convolution (G-DC-FFT)
// [[Rcpp::export]]
NumericVector dgpb_dc(const IntegerVector obs, const NumericVector probs, const IntegerVector val_p, const IntegerVector val_q){//, const int splits = -1){
  // number of probabilities of success
  const int sizeIn = probs.length();
  // determine pairwise minimum and maximum
  IntegerVector v = pmin(val_p, val_q);
  IntegerVector u = pmax(val_p, val_q);
  // theoretical minimum
  const int min_v = sum(v);
  // compute differences
  IntegerVector d = u - v;
  // final output size
  const int sizeOut = sum(d) + 1;
  // greatest common divisor of the differences
  const int gcd = vectorGCD(d[d > 0]);
  // rescale differences according to gcd
  if(gcd > 1) d = d / gcd;
  // theoretical rescaled maximum
  const int max_rescaled = sum(d);
  // output size
  const int sizeOut_rescaled = max_rescaled + 1;
  
  // if val_p[i] was not the larger one, the respective probs[i] has to be 'flipped'
  // furthermore: if difference is 0 (i.e. u[i] equals v[i]), a non-zero outcome is impossible
  NumericVector probs_flipped(sizeIn);
  for(int i = 0; i < sizeIn; i++){
    if(!d[i]) probs_flipped[i] = 0.0;
    else if(val_p[i] < u[i]) probs_flipped[i] = 1 - probs[i];
      else probs_flipped[i] = probs[i];
  }
  
  // results vectors
  NumericVector results(sizeOut);
  NumericVector results_rescaled;
  
  // if max_rescaled equals input size, we have an ordinary poisson binomial distribution
  if(max(d) == 1){
    // compute ordinary distribution
    results_rescaled = dpb_dc(IntegerVector(), probs_flipped);
  }else{
    // number of tree splits
    //int num_splits = splits < 0 ? std::max<int>(0, (int)std::ceil(std::log(sizeIn / 860) / std::log(2.0))) : splits;
    int num_splits = sizeIn > 860 ? std::max<int>(0, (int)std::ceil(std::log(sizeIn / 860) / std::log(2.0))) : 0;
    // direct convolution is sufficient in case of 0 splits
    if(num_splits == 0) return dgpb_conv(obs, probs_flipped, u, v);
    // number of groups
    int num_groups = std::pow(2, num_splits);
    // fraction of total size per group
    double frac = (double)(sizeOut_rescaled - 1)/num_groups;
    // reduce number of splits and groups if inner-group sizes are too large
    while(num_splits > 0 && (num_groups > sizeIn || frac < max(d))){
      num_splits -= 1;
      num_groups /= 2;
      frac *= 2;
    }
    // direct convolution is sufficient, if no splits are necessary
    if(num_splits == 0) return dgpb_conv(obs, probs_flipped, u, v);
    
    // compute group sizes with minimum size disparity
    IntegerVector group_sizes(num_groups);
    IntegerVector group_indices(sizeIn, -1);
    
    // assign each probability and outcome to a group
    IntegerVector ord = order(NumericVector(d), true);
    IntegerVector d_ordered = d[ord];
    probs_flipped = probs_flipped[ord];
    NumericVector remainder(num_groups, frac);
    int g = 0;
    int inc = 1;
    for(int i = 0; i < sizeIn; i++){
      checkUserInterrupt();
      if(g == num_groups || g == -1){
        inc *= -1;
        g += inc;
      }
      if(d_ordered[i] > remainder[g]){
        g = 0;
        for(int j = 1; j < num_groups; j++){
          if(remainder[j] > remainder[g]) g = j;
        }
      }
      group_sizes[g] += d_ordered[i];
      remainder[g] -= d_ordered[i];
      group_indices[i] = g;
      g += inc;
    }
    
    // compute first and last indices of the groups
    IntegerVector group_starts(num_groups);
    IntegerVector group_ends(num_groups);
    group_starts[0] = 0;
    group_ends[0] = group_sizes[0] - 1;
    for(int i = 1; i < num_groups; i++){
      group_starts[i] = group_starts[i - 1] + group_sizes[i - 1];
      group_ends[i] = group_ends[i - 1] + group_sizes[i];
    }
    
    // results vector; direct convolution will increase size of each group by 1
    results_rescaled = NumericVector(sizeOut_rescaled - 1 + num_groups);
    
    // compute direct convolutions for each group
    int start = 0, end = 0;
    for(int i = 0; i < num_groups; i++){
      checkUserInterrupt();
      // compute new starting and ending indices, because groups grow by 1
      start = group_starts[i] + i;
      end = group_ends[i] + i + 1;
      
      // target range
      Range target(start, end);
      u = d_ordered[group_indices == i];
      v = IntegerVector(u.length());
      
      // direct convolution
      results_rescaled[target] = dgpb_conv(IntegerVector(), probs_flipped[group_indices == i], u, v);
      
      // update starting and ending positions
      group_starts[i] = start;
      group_ends[i] = end;
    }
    
    int num_groups_reduced = num_groups / 2;
    while(num_splits > 0){
      for(int i = 0; i < num_groups_reduced; i++){
        checkUserInterrupt();
        // compute new starting and ending indices, because group sizes are
        // reduced by 1, due to FFT convolution
        start = group_starts[2*i] - i;
        end = group_ends[2*i + 1] - i - 1;
        
        // target range
        Range target(start, end);
        // FFT convolution
        results_rescaled[target] = fft_probs(results_rescaled[Range(group_starts[2*i], group_ends[2*i])], results_rescaled[Range(group_starts[2*i + 1], group_ends[2*i + 1])]);
        
        // update starting and ending indices
        group_starts[i] = start;
        group_ends[i] = end;
      }
      num_groups_reduced /= 2;
      num_splits -= 1;
    }
  }
  // "correct" numerically false (and thus useless) results
  results_rescaled[results_rescaled < 5.55e-17] = 0;
  results_rescaled[results_rescaled > 1] = 1;
  
  // make sure that probability masses sum up to 1
  norm_dpb(results_rescaled);
  
  // map results to generalized distribution (scale-back)
  for(int i = 0; i < sizeOut_rescaled; i++)
    results[i * gcd] = results_rescaled[i];
  
  // return final results
  if(obs.length()) return results[obs - min_v]; else return results;
}

// [[Rcpp::export]]
NumericVector pgpb_dc(const IntegerVector obs, const NumericVector probs, const IntegerVector val_p, const IntegerVector val_q, const bool lower_tail = true){
  // theoretical minimum
  const int min_v = sum(pmin(val_p, val_q));
  // theoretical maximum
  const int max_v = sum(pmax(val_p, val_q));
  // maximum observed value
  const int max_q = obs.length() ? max(obs) : max_v;
  
  // probability masses
  const NumericVector pmf = dgpb_dc(IntegerVector(), probs, val_p, val_q);
  
  // compute CDF
  NumericVector results = ppb_generic(obs - min_v, pmf, lower_tail);
  
  // ensure that sum = 1 (or 0), if last value equals the highest observable value
  if(obs.length()){
    if(max_q == max_v) results[obs == max_v] = (double)lower_tail;
  }else results[max_v - min_v] = (double)lower_tail;
  
  // return final results
  return results;
}

// Generalized Discrete Fourier Transformation of Characteristic Function (G-DFT-CF)
// [[Rcpp::export]]
NumericVector dgpb_dftcf(const IntegerVector obs, const NumericVector probs, const IntegerVector val_p, const IntegerVector val_q){
  // number of probabilities of success
  const int sizeIn = probs.length();
  // determine pairwise minimum and maximum
  const IntegerVector v = pmin(val_p, val_q);
  const IntegerVector u = pmax(val_p, val_q);
  // compute differences
  IntegerVector d = u - v;
  // final output size
  const int sizeOut = sum(d) + 1;
  // greatest common divisor of the differences
  const int gcd = vectorGCD(d[d > 0]);
  // rescale differences according to gcd
  if(gcd > 1) d = d / gcd;
  // theoretical rescaled maximum
  const int max_rescaled = sum(d);
  // output size
  const int sizeOut_rescaled = max_rescaled + 1;
  
  // if val_p[i] was not the larger one, the respective probs[i] has to be 'flipped'
  // furthermore: if difference is 0 (i.e. u[i] equals v[i]), a non-zero outcome is impossible
  NumericVector probs_flipped(sizeIn);
  for(int i = 0; i < sizeIn; i++){
    if(!d[i]) probs_flipped[i] = 0.0;
    else if(val_p[i] < u[i]) probs_flipped[i] = 1 - probs[i];
      else probs_flipped[i] = probs[i];
  }
  
  // results vectors
  NumericVector results(sizeOut);
  NumericVector results_rescaled;
  
  // if max_rescaled equals input size, we have an ordinary poisson binomial distribution
  if(max(d) == 1){
    // compute ordinary distribution
    results_rescaled = dpb_dftcf(IntegerVector(), probs_flipped[d > 0]);
  }else{
    // "initialize" DFT input vector
    fftw_complex *input_fft;
    input_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sizeOut_rescaled);
    input_fft[0][REAL] = 1.0;
    input_fft[0][IMAG] = 0.0;
    
    // initialize complex numbers for "C" and "C to the power of i"
    std::vector< std::complex<double> > C(sizeIn, 1.0);
    std::vector< std::complex<double> > C_power(sizeIn, 1.0);
    for(int i = 0; i < sizeIn; i++){
      if(d[i]) C[i] = exp(std::complex<double>(0.0, d[i] * 2.0) * M_PI / ((double)sizeOut_rescaled));
    }
    
    // compute closed-form expression of Hernandez and Williams
    for(int l = 1; l <= sizeOut_rescaled / 2; l++){
      checkUserInterrupt();
      std::complex<double> product = 1.0;
      for(int k = 0; k < sizeIn; k++){
        if(probs_flipped[k]){
          if(d[k]) C_power[k] *= C[k];////
          product *= 1.0 + probs_flipped[k] * (C_power[k] - 1.0);
        }
      }
      
      input_fft[l][REAL] = product.real();
      input_fft[l][IMAG] = product.imag();
      input_fft[sizeOut_rescaled - l][REAL] = product.real();
      input_fft[sizeOut_rescaled - l][IMAG] = -product.imag();
    }
    
    // vector of DFT results
    fftw_complex *result_fft;
    result_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sizeOut_rescaled);
    
    // perform DFT
    fftw_plan planDFT;
    planDFT = fftw_plan_dft_1d(sizeOut_rescaled, input_fft, result_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(planDFT);
    
    // gather results
    results_rescaled = NumericVector(sizeOut_rescaled);
    for(int i = 0; i < sizeOut_rescaled; i++) results_rescaled[i] = result_fft[i][REAL] / sizeOut_rescaled;
    
    // garbage collection
    fftw_destroy_plan(planDFT);
    fftw_free(input_fft);
    fftw_free(result_fft);
  }
  // "correct" numerically false (and thus useless) results
  results_rescaled[results_rescaled < 2.22e-16] = 0;
  results_rescaled[results_rescaled > 1] = 1;
  
  // make sure that probability masses sum up to 1
  norm_dpb(results_rescaled);
  
  // map results to generalized distribution (scale-back)
  for(int i = 0; i < sizeOut_rescaled; i++)
    results[i * gcd] = results_rescaled[i];
  
  // return final results
  if(obs.length()) return results[obs - sum(v)]; else return results;
}

// [[Rcpp::export]]
NumericVector pgpb_dftcf(const IntegerVector obs, const NumericVector probs, const IntegerVector val_p, const IntegerVector val_q, const bool lower_tail = true){
  // theoretical minimum
  const int min_v = sum(pmin(val_p, val_q));
  // theoretical maximum
  const int max_v = sum(pmax(val_p, val_q));
  // maximum observed value
  const int max_q = obs.length() ? max(obs) : max_v;
  
  // probability masses
  const NumericVector pmf = dgpb_dftcf(IntegerVector(), probs, val_p, val_q);
  
  // compute CDF
  NumericVector results = ppb_generic(obs - min_v, pmf, lower_tail);
  
  // ensure that sum = 1 (or 0), if last value equals the highest observable value
  if(obs.length()){
    if(max_q == max_v) results[obs == max_v] = (double)lower_tail;
  }else results[max_v - min_v] = (double)lower_tail;
  
  // return final results
  return results;
}

// [[Rcpp::export]]
NumericVector pgpb_na(const IntegerVector obs, const NumericVector probs, const IntegerVector val_p, const IntegerVector val_q, const bool refined = true, const bool lower_tail = true){
  // number of probabilities of success
  const int sizeIn = probs.length();
  // determine pairwise minimum and maximum
  const IntegerVector v = pmin(val_p, val_q);
  const IntegerVector u = pmax(val_p, val_q);
  // theoretical (unshifted!) minimum
  const int min_v = sum(v);
  // compute differences
  IntegerVector d = u - v;
  // maximum observable value
  const int max_q = sum(d);
  // greatest common divisor of the differences
  const int gcd = vectorGCD(d[d > 0]);
  // rescale differences according to gcd
  if(gcd > 1) d = d / gcd;
  // theoretical (shifted!) maximum
  const int max_rescaled = max_q / gcd;
  // rescaled observations
  const IntegerVector obs_range = obs.length() ? (obs - min_v) / gcd : IntegerVector(Range(0, max_rescaled));
  
  // if val_p[i] was not the larger one, the respective probs[i] has to be 'flipped'
  // furthermore: if difference is 0 (i.e. u[i] equals v[i]), a non-zero outcome is impossible
  NumericVector probs_flipped(sizeIn);
  for(int i = 0; i < sizeIn; i++){
    if(!d[i]) probs_flipped[i] = 0.0;
    else if(val_p[i] < u[i]) probs_flipped[i] = 1 - probs[i];
    else probs_flipped[i] = probs[i];
  }
  
  // if maximum difference equals 1, we have an ordinary poisson binomial distribution
  if(max(d) == 1){
    // compute ordinary distribution
    return ppb_na(obs_range, probs_flipped[d > 0], refined, lower_tail);
  }else{
    // mu
    const double mu = sum(probs_flipped * NumericVector(d));
    // p * q
    const NumericVector pq = probs_flipped * (1 - probs_flipped);
    // sigma
    const double sigma = std::sqrt(sum(pq * pow(NumericVector(d), 2.0)));
    // standardized observations with continuity correction
    const NumericVector obs_std = (NumericVector(obs_range) + 0.5 - mu)/sigma;
    // cumulative probabilities
    NumericVector results_rescaled = pnorm(obs_std, 0.0, 1.0, lower_tail);
    // refine
    if(refined && sigma){
      // gamma
      const double gamma = sum(pq * (1 - 2 * probs_flipped) * pow(NumericVector(d), 3.0))/std::pow(sigma, 3.0);
      // probabilities
      if(lower_tail)
        results_rescaled += gamma * (1 - pow(obs_std, 2.0)) * dnorm(obs_std) / 6;
      else
        results_rescaled += -gamma * (1 - pow(obs_std, 2.0)) * dnorm(obs_std) / 6;
    }
    
    // make sure that all probabilities do not exceed 1 and are at least 0
    results_rescaled[results_rescaled < 0] = 0;
    results_rescaled[results_rescaled > 1] = 1;
    
    // make sure largest possible value has cumulative probability of 1
    if(max_q >= max_rescaled) results_rescaled[obs_range >= max_rescaled] = (double)lower_tail;
    
    // return final results
    if(obs.length() || gcd == 1) return results_rescaled;
    else{
      // map results to generalized distribution (scale-back)
      NumericVector results(max_q + 1);
      for(int i = 0; i <= max_q; i++)
          results[i] = results_rescaled[i/gcd];
      
      return results;
    }
  }
}

// Generalized Normal Approximations (G-NA, G-RNA)
// [[Rcpp::export]]
NumericVector dgpb_na(const IntegerVector obs, const NumericVector probs, const IntegerVector val_p, const IntegerVector val_q, const bool refined = true){
  // smallest possible value
  const int min_v = sum(pmin(val_p, val_q));
  // highest observed value
  const int max_q = obs.length() ? max(obs) : sum(pmax(val_p, val_q));
  // rounded down expectation + 0.5 (continuity correction)
  const int mid = (int)floor(sum(probs * NumericVector(val_p) + (1 - probs) * NumericVector(val_q)) + 0.5);
  
  // cumulative probabilities
  NumericVector cdf_lower = pgpb_na(IntegerVector(Range(min_v, std::min<int>(mid, max_q))), probs, val_p, val_q, refined, true);
  NumericVector cdf_upper = pgpb_na(IntegerVector(Range(std::min<int>(mid, max_q), max_q)), probs, val_p, val_q, refined, false);
  
  // vector to store results
  NumericVector results(max_q - min_v + 1);
  
  // compute probability masses
  results[0] = cdf_lower[0];
  for(int i = 1; i <= max_q - min_v; i++){
    if(i + min_v <= mid) 
      results[i] = cdf_lower[i] - cdf_lower[i - 1]; 
    else results[i] = cdf_upper[i - 1 - mid + min_v] - cdf_upper[i - mid + min_v];
  }
  
  // compute and return results
  if(obs.length()) return results[obs - min_v]; else return results;
}

// Bernoulli Random Number Generator
// [[Rcpp::export]]
IntegerVector rgpb_bernoulli(const int n, const NumericVector probs, const IntegerVector val_p, const IntegerVector val_q){
  // number of probabilities of success
  const int size = probs.length();
  // sum of values that occur with probability q = 1 - p
  const double sum_v = (double)sum(val_q);
  // differences
  const IntegerVector d = val_p - val_q;
  
  // vector to store results
  NumericVector results(n, sum_v);
  
  // generate random numbers
  for(int i = 0; i < size; i++) 
    for(int j = 0; j < n; j++)
      results[j] += d[i] * R::rbinom(1.0, probs[i]);
  
  // return results
  return IntegerVector(results);
}