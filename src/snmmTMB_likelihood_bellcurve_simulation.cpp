// snmmTMB_likelihood_bellcurve_simulation_corr.cpp
// ---------------------------------------------------------------------------
// SNMM likelihood for the bell-curve simulation of manuscript Section 3.
//
// Model:
//   y_ij = alpha + b1_i + f(t_ij - b2_i) + eps_ij
//   b_i = (b1_i, b2_i) ~ N(0, D),
//        D = [ sd_b1^2            rho*sd_b1*sd_b2 ]
//            [ rho*sd_b1*sd_b2    sd_b2^2         ]
//   eps_ij ~ N(0, sigma^2)
//
// b1_i is a subject vertical shift and b2_i a subject horizontal shift. The
// two are correlated, through rho = tanh(transf_rho), which stays inside
// (-1,1) for any real transf_rho and keeps D positive definite. Holding
// transf_rho fixed with map() in MakeADFun() fixes rho.
//
// The range of t_ij - b2_i is not known a priori, so the knots are fixed on
// [0,1] and the shifted covariate is scaled into that interval
// (Section 2.3.1, second case):
//
//   v_ij = (t_ij - b2_i - a) / s,   a = min(t) - 3*sd_b2,
//                                   s = max(t) + 3*sd_b2 - a.
//
// The bounds depend on sd_b2, so the mapping is re-evaluated as that estimate
// changes, while the basis and the penalty stay fixed.
//
// The centering constants Bmean are computed once in R, as the column means of
// B over a fixed grid spanning the knot range [0,1] (on the scaled v axis), and
// passed in as data.
//
// Data:
//   y, t               response and time values
//   group, nGroup      subject index (1 to nGroup) and number of subjects
//   subj_flag          which subjects to report subject-level curves for
//   knots, degree, K   B-spline basis specification (knots on [0,1])
//   Bmean              fixed-grid column means of the basis (length K-1)
//   Upos, dpos, U0     penalty eigendecomposition (see starting_points.cpp)
//   spline_ci          flag: 1 = evaluate the population curve on t_grid
//   t_grid             grid of t-values at which to report curves
//
// Parameters:
//   alpha                    population intercept
//   b1, b2                   per-subject random effects (vertical shift,
//                            horizontal shift)
//   log_sd_b1, log_sd_b2     log standard deviations of b1, b2
//   transf_rho               atanh of the correlation between b1 and b2
//   log_sigma                log residual standard deviation
//   vpos, gamma0             spline coefficients
//   log_lambda               smoothing parameter
//
// Reported:
//   rho        REPORT and ADREPORT - correlation on its natural scale
//   h_grid     ADREPORT - population curve on t_grid at b = 0, including
//              alpha
//   mu_sel     ADREPORT - subject curves for the subjects flagged in subj_flag
//
// ---------------------------------------------------------------------------
#include <TMB.hpp>

template<class Type>
vector<Type> bspline_basis(Type x, vector<Type> knots, int degree) {

  // AD-safe functions for comparisons
  using CppAD::CondExpGt; // greater than
  using CppAD::CondExpGe; // greater than or equal
  using CppAD::CondExpLt; // less than

  // number of knots and basis functions
  int nKnots = knots.size();
  int nBasis = nKnots - degree - 1;

  // initialize basis vector to zero
  vector<Type> N(nBasis);
  for(int i=0;i<nBasis;i++) N(i) = Type(0.0);


  Type eps = Type(1e-8);

  // Degree 0 basis
  for(int i = 0; i < nBasis; ++i) {

    // indicator function = (x >= knots[i]) && (x < knots[i+1])
    Type ge_left = CondExpGe(x, knots(i), Type(1.0), Type(0.0));
    Type lt_right = CondExpLt(x, knots(i+1), Type(1.0), Type(0.0));

    Type indicator = ge_left * lt_right;

    N(i) = indicator;
  }

  // Cox–de Boor recursion
  for(int p = 1; p <= degree; ++p){

    vector<Type> Np(nBasis);
    for(int i = 0; i < nBasis; ++i) Np(i) = Type(0.0);

    for(int i = 0; i < nBasis; ++i){
      // left term
      Type leftDen = knots(i + p) - knots(i);
      Type leftTerm = Type(0.0);
      leftTerm = CondExpGt(leftDen, eps,
                           ((x - knots(i)) / leftDen) * N(i),
                           Type(0.0));
                           
      // right term
      Type rightDen = knots(i + p + 1) - knots(i + 1);
      Type rightTerm = Type(0.0);

      // null right term for the last basis function
      if((i + 1) < nBasis) {
        rightTerm = CondExpGt(rightDen, eps,
                              ((knots(i + p + 1) - x) / rightDen) * N(i + 1),
                              Type(0.0));
      } else {
        rightTerm = Type(0.0);
      }

      
      Np(i) = leftTerm + rightTerm;
    }
    N = Np;
  }

  return N;
}













template<class Type>
Type objective_function<Type>::operator() ()
{
  // data
  DATA_VECTOR(y);          // response
  DATA_VECTOR(t);            // time values
  DATA_IVECTOR(group);       // group index 1 to nGroup
  DATA_INTEGER(nGroup);

  DATA_VECTOR(knots);        // knot vector
  DATA_INTEGER(degree);      // spline degree
  DATA_INTEGER(K);           // number of basis functions

  DATA_VECTOR(Bmean);        // fixed-grid centering constants (length K-1)

  DATA_INTEGER(spline_ci);
  DATA_VECTOR(t_grid);
  DATA_IVECTOR(subj_flag);

  // penalty eigendecomposition
  DATA_MATRIX(Upos);       
  DATA_VECTOR(dpos);        
  DATA_MATRIX(U0);           

  // parameters 
  PARAMETER(alpha);                    // intercept
  PARAMETER_VECTOR(b1);                // random intercepts (vertical shift)
  PARAMETER_VECTOR(b2);                // random horizontal shift
  PARAMETER(log_sd_b1);                // sd for b1
  PARAMETER(log_sd_b2);                // sd for b2
  PARAMETER(log_sigma);                // residual sd
  PARAMETER(transf_rho);               // atanh(correlation between b1 and b2)

  PARAMETER_VECTOR(vpos);              // penalized spline coordinates
  PARAMETER_VECTOR(gamma0);            // unpenalized spline coordinates
  PARAMETER(log_lambda);               // smoothing parameter

  int N = y.size();
  Type sd_b1 = exp(log_sd_b1);
  Type sd_b2 = exp(log_sd_b2);
  Type sigma = exp(log_sigma);
  Type rho = tanh(transf_rho);         // in (-1, 1) by construction
  Type lambda = exp(log_lambda);







  
  // reconstruct coefficients c = Upos * vpos + U0 * gamma0
  int Kint = K-1;
  vector<Type> c(Kint);
  for(int k=0;k<Kint;k++){
    Type val = Type(0.0);
    if(Upos.cols() > 0){
      for(int j=0;j<Upos.cols(); j++) val += Upos(k,j) * vpos(j);
    }
    if(U0.cols() > 0){
      for(int j=0;j<U0.cols(); j++) val += U0(k,j) * gamma0(j);
    }
    c(k) = val;
  }



  // shift data points according to random horizontal shift b2
  vector<Type> u(N);
  for(int i=0;i<N;i++){
    int g = group(i) - 1;
    u(i) = t(i) - b2(g);
  }



  // bounds of the scaling map (Section 2.3.1): +-3 sd_b2 around the observed
  // range of t, so that t - b2_i lands inside [0,1] for plausible b2_i
  Type t_min = t(0);
  Type t_max = t(0);
  for(int i = 1; i < N; ++i){
  if(t(i) < t_min) t_min = t(i);
  if(t(i) > t_max) t_max = t(i);
  }

  Type umin = t_min - Type(3) * sd_b2; // lower quantile
  Type umax = t_max + Type(3) * sd_b2; // upper quantile

  Type a = umin;
  Type s = (umax - umin);   // range




  // likelihood
  Type nll = Type(0.0);
  for(int i=0;i<N;i++){
    int g = group(i) - 1;
    Type v = (u(i) - a) / s;
    vector<Type> B = bspline_basis(v, knots, (int)degree);
    Type h = Type(0.0);
    for(int k=0;k<Kint;k++){
      B(k) -= Bmean(k);   // center column k
      h += c(k) * B(k);
    }

    Type mu = alpha + b1(g) + h;

    nll -= dnorm(y(i), mu, sigma, true);
  }

  // random effects: bivariate normal with unstructured 2x2 covariance, written
  // out in closed form rather than through density::MVNORM_t
  Type two_pi = Type(6.283185307179586476925286766559);
  Type one_m_rho2 = Type(1.0) - rho * rho;
  Type log_norm_const = log(two_pi) + log_sd_b1 + log_sd_b2
                        + Type(0.5) * log(one_m_rho2);
  for(int j=0;j<nGroup;j++){
    Type z1 = b1(j) / sd_b1;
    Type z2 = b2(j) / sd_b2;
    Type quad = (z1*z1 - Type(2.0)*rho*z1*z2 + z2*z2) / one_m_rho2;
    nll += log_norm_const + Type(0.5) * quad;
  }

  // penalized spline terms
  if(dpos.size() > 0){
    int mpos = dpos.size();
    matrix<Type> Sigma_v(mpos, mpos);
    Sigma_v.setZero();
    for(int j=0;j<mpos;j++) Sigma_v(j,j) = Type(1.0) / ( lambda * dpos(j) );
    density::MVNORM_t<Type> mvn(Sigma_v);
    nll += mvn(vpos);
  }























  // Spline evaluation for population curve

  int ngrid = t_grid.size();
  vector<Type> h_grid(ngrid);

  if(spline_ci == 1){
    for(int ig = 0; ig < ngrid; ++ig){

      Type ui = t_grid(ig);
      Type v = (ui - a) / s;

      Type left = knots(0);
      Type right = knots(knots.size()-1);
      Type tiny_eps = Type(1e-12);
      if(v < left) v = left;
      if(v > right - tiny_eps) v = right - tiny_eps;

      vector<Type> B = bspline_basis(v, knots, (int)degree);
      Type hval = Type(0.0);
      for(int k=0;k<Kint;k++){
        B(k) -= Bmean(k);   // center column k
        hval += c(k) * B(k);
      }

      hval += alpha;
      h_grid(ig) = hval;
    }
  }






  // Spline evaluation for subject curves of flagged subjects
  int nsel = 0;
  for(int is = 0; is < nGroup; ++is) if(subj_flag[is] == 1) ++nsel;

  vector<Type> mu_sel( nsel * ngrid ); 
  int sel_idx = 0;

  for(int is = 0; is < nGroup; ++is){
    if(subj_flag[is] == 0) continue; 

    Type b1_subj = b1(is);
    Type b2_subj = b2(is);


    for(int ig = 0; ig < ngrid; ++ig){

      Type ui = t_grid(ig) - b2_subj;
      Type v = (ui - a) / s;


      vector<Type> B = bspline_basis(v, knots, (int)degree); 
      
      Type hval = Type(0.0);
      for(int k=0; k<Kint; ++k){
        B(k) -= Bmean(k);
        hval += c(k) * B(k); 
      }

      Type mu_ij = alpha + b1_subj + hval;


      mu_sel[ sel_idx * ngrid + ig ] = mu_ij;
    }
    ++sel_idx;
  } 

  // rho on its natural scale, so sdreport() returns a delta-method SE for the
  // correlation rather than for atanh(rho)
  REPORT(rho);
  ADREPORT(rho);

  // ADREPORT for spline evaluations
  if(spline_ci == 1) ADREPORT(h_grid);
  if(spline_ci == 1) ADREPORT(mu_sel);




  return nll;
}
