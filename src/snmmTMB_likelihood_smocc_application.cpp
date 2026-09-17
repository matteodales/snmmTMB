// smocc_data_tmb_model2_scaledknots_fixedgrid_unstructured.cpp
// ---------------------------------------------------------------------------
// SNMM likelihood for the SMOCC height application of manuscript Section 4.
//
// Model:
//   y_ij = beta_intercept + beta_intercept_sex * sex_i + b_intercept_i
//          + exp(beta_amplitude_sex * sex_i)
//            * f( age_ij + ga_i * beta_shift_ga + b_shift_i ) + eps_ij
//
//   ( b_intercept_i )         ( sd1^2        rho*sd1*sd2 )
//   (               ) ~ N2 0,                             ,
//   (   b_shift_i   )         ( rho*sd1*sd2  sd2^2       )
//
//   eps_ij ~ N(0, sigma^2),   sd1 = sd_b_intercept, sd2 = sd_b_shift
//
// b_intercept_i is a subject vertical shift and b_shift_i a subject shift along
// the age axis. The two are correlated through rho = tanh(transf_rho).
//
// The range of the spline argument is not known a priori, so the knots are
// fixed and the shifted age is scaled into their interval (Section 2.3.1,
// second case):
//
//   v_ij = (age_ij + ga_i*beta_shift_ga + b_shift_i - a) / s,
//   a = age_min - 3*sd_b_shift + ga_min*beta_shift_ga,
//   s = age_max + 3*sd_b_shift + ga_max*beta_shift_ga - a.
//
//
// The centering constants Bmean are computed once in R, as the column means of
// B over a fixed grid in v, and passed in as data.
//
// sex_subj and ga_subj have length nGroup and are indexed by subject; sex and
// ga have length N and are indexed by observation.
//
// Data:
//   y, age               response and age values
//   age_min, age_max     range of age, for the knot-scaling map
//   ga_min, ga_max       range of ga, for the knot-scaling map
//   sex, ga              observation-level covariates (length N)
//   sex_subj, ga_subj    subject-level covariates (length nGroup)
//   group, nGroup        subject index (1 to nGroup) and number of subjects
//   subj_flag            which subjects to report subject-level curves for
//   knots, degree, K     B-spline basis specification (on the v scale)
//   Bmean                fixed-grid column means of the basis (length K-1)
//   Upos, dpos, U0       penalty eigendecomposition (see starting_points.cpp)
//   spline_ci            flag: 1 = evaluate the population curve on age_grid
//   age_grid             grid of ages at which to report curves
//
// Parameters:
//   beta_intercept           population intercept
//   beta_intercept_sex       sex effect on the intercept
//   beta_amplitude_sex       sex effect on the curve amplitude (log scale)
//   beta_shift_ga            gestational-age effect on the horizontal shift
//   b_intercept, b_shift     per-subject random effects (vertical shift,
//                            age shift)
//   log_sd_b_intercept       log standard deviation of b_intercept
//   log_sd_b_shift           log standard deviation of b_shift
//   transf_rho               atanh of the correlation between the two
//   log_sigma                log residual standard deviation
//   vpos, gamma0             spline coefficients
//   log_lambda               smoothing parameter
//
// Reported:
//   rho        REPORT and ADREPORT - correlation on its natural scale
//   c          REPORT - reconstructed spline coefficients
//   a, s       REPORT - offset and scale of the knot-scaling map
//   h_grid     ADREPORT - population curve on age_grid at the reference
//              level (sex = 0, ga = 0, no random effects), including
//              beta_intercept
//   mu_sel     ADREPORT - subject curves for the subjects flagged in subj_flag
//
// ---------------------------------------------------------------------------
#include <TMB.hpp>

// de Boor spline basis evaluation
template<class Type>
vector<Type> bspline_basis(Type x_in, vector<Type> knots, int degree) {
  using CppAD::CondExpGt;
  using CppAD::CondExpGe;
  using CppAD::CondExpLt;
  using CppAD::CondExpLe;
  using CppAD::abs;

  int nKnots = knots.size();
  int nBasis = nKnots - degree - 1;
  vector<Type> N(nBasis);
  for(int i=0;i<nBasis;i++) N(i) = Type(0.0);
  if(nBasis <= 0) return N;

  // scale-aware eps
  Type knot_left = knots(0);
  Type knot_right = knots(nKnots - 1);
  Type span = CondExpGt(knot_right - knot_left, Type(0.0), knot_right - knot_left, Type(1.0));
  const double base_eps = 1e-8;
  Type eps = Type(base_eps) * span;

  // clamp into [knot_left, knot_right - eps]
  Type x = x_in;
  x = CondExpLt(x, knot_left, knot_left, x);
  x = CondExpGt(x, knot_right - eps, knot_right - eps, x);

  // Degree 0 basis
  for(int i = 0; i < nBasis; ++i) {
    Type ge_left = CondExpGe(x, knots(i), Type(1.0), Type(0.0));
    Type lt_right = CondExpLt(x, knots(i+1), Type(1.0), Type(0.0));
    Type near_right = CondExpLe(abs(x - knot_right), eps, Type(1.0), Type(0.0));
    Type last_marker = CondExpGe(Type(i), Type(nBasis - 1), near_right, Type(0.0));
    Type indicator = ge_left * lt_right + last_marker;
    N(i) = indicator;
  }

  // Cox-de Boor recursion
  for(int p = 1; p <= degree; ++p){
    vector<Type> Np(nBasis);
    for(int i = 0; i < nBasis; ++i) Np(i) = Type(0.0);
    for(int i = 0; i < nBasis; ++i){
      Type leftDen = knots(i + p) - knots(i);
      Type leftTerm = Type(0.0);
      leftTerm = CondExpGt(leftDen, eps,
                           ((x - knots(i)) / leftDen) * N(i),
                           Type(0.0));
      Type rightDen = knots(i + p + 1) - knots(i + 1);
      Type rightTerm = Type(0.0);
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
  DATA_VECTOR(y);
  DATA_VECTOR(age);
  DATA_SCALAR(age_min);
  DATA_SCALAR(age_max);
  DATA_SCALAR(ga_min);
  DATA_SCALAR(ga_max);

  DATA_VECTOR(knots);
  DATA_INTEGER(degree);
  DATA_INTEGER(K);

  DATA_VECTOR(Bmean);        // fixed-grid centering constants (length K-1)

  DATA_INTEGER(spline_ci);
  DATA_VECTOR(age_grid);
  DATA_IVECTOR(subj_flag);

  DATA_VECTOR(sex);          // observation level, length N
  DATA_VECTOR(ga);           // observation level, length N
  DATA_VECTOR(sex_subj);     // subject level, length nGroup
  DATA_VECTOR(ga_subj);      // subject level, length nGroup

  // penalty eigendecomposition from R
  DATA_MATRIX(Upos);
  DATA_VECTOR(dpos);
  DATA_MATRIX(U0);

  DATA_IVECTOR(group);
  DATA_INTEGER(nGroup);

  // parameters
  PARAMETER(beta_intercept);
  PARAMETER(beta_intercept_sex);
  PARAMETER(beta_amplitude_sex);
  PARAMETER(beta_shift_ga);
  PARAMETER_VECTOR(b_intercept);
  PARAMETER_VECTOR(b_shift);
  PARAMETER(log_sd_b_intercept);
  PARAMETER(log_sd_b_shift);
  PARAMETER(transf_rho);              // atanh(correlation)
  PARAMETER(log_sigma);

  PARAMETER_VECTOR(vpos);
  PARAMETER_VECTOR(gamma0);
  PARAMETER(log_lambda);

  int N = y.size();
  Type sd_b_intercept = exp(log_sd_b_intercept);
  Type sd_b_shift = exp(log_sd_b_shift);
  Type rho = tanh(transf_rho);
  Type sigma = exp(log_sigma);
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

  // horizontal shift of the design points
  vector<Type> u_all(N);
  for(int i=0;i<N;i++){
    int g = group(i) - 1;
    u_all(i) = (age(i) + ga(i)*beta_shift_ga + b_shift(g));
  }

  // knot-scaling map
  Type umin = (age_min - Type(3) * sd_b_shift + ga_min*beta_shift_ga);
  Type umax = (age_max + Type(3) * sd_b_shift + ga_max*beta_shift_ga);
  Type a = umin;
  Type s = (umax - umin);

  // The sum-to-zero centering constants are the data vector Bmean, computed in
  // R on a fixed grid in v, so they are constant here.

  // likelihood
  Type nll = Type(0.0);
  for(int i=0;i<N;i++){
    int g = group(i) - 1;
    Type v = (u_all(i) - a) / s;
    vector<Type> B = bspline_basis(v, knots, (int)degree);
    Type h = Type(0.0);
    for(int k=0;k<Kint;k++){
      B(k) -= Bmean(k);
      h += c(k) * B(k);
    }

    Type mu = beta_intercept + beta_intercept_sex*sex(i) + b_intercept(g)
              + exp(beta_amplitude_sex*sex(i))*h;

    nll -= dnorm(y(i), mu, sigma, true);
  }

  // random effects: bivariate normal with unstructured 2x2 covariance, written
  // out in closed form
  Type two_pi = Type(6.283185307179586476925286766559);
  Type one_m_rho2 = Type(1.0) - rho * rho;
  Type log_norm_const = log(two_pi) + log_sd_b_intercept + log_sd_b_shift
                        + Type(0.5) * log(one_m_rho2);
  for(int j=0;j<nGroup;j++){
    Type z1 = b_intercept(j) / sd_b_intercept;
    Type z2 = b_shift(j) / sd_b_shift;
    Type quad = (z1*z1 - Type(2.0)*rho*z1*z2 + z2*z2) / one_m_rho2;
    nll += log_norm_const + Type(0.5) * quad;
  }

  if(dpos.size() > 0){
    int mpos = dpos.size();
    for(int j=0;j<mpos;j++){
      nll -= dnorm(vpos(j), Type(0.0), Type(1.0) / sqrt( lambda * dpos(j) ), true);
    }
  }

  // ---- population curve (reference level: sex = 0, ga = 0, no random effects)
  int ngrid = age_grid.size();
  vector<Type> h_grid(ngrid);

  if(spline_ci == 1){
    for(int ig = 0; ig < ngrid; ++ig){

      Type ui = (age_grid(ig));
      Type v = (ui - a) / s;

      Type left = knots(0);
      Type right = knots(knots.size()-1);
      Type tiny_eps = Type(1e-12);
      if(v < left) v = left;
      if(v > right - tiny_eps) v = right - tiny_eps;

      vector<Type> B = bspline_basis(v, knots, (int)degree);

      Type hval = Type(0.0);
      for(int k=0;k<Kint;k++){
        B(k) -= Bmean(k);
        hval += c(k) * B(k);
      }

      h_grid(ig) = beta_intercept + hval;
    }
  }

  // ---- subject curves for flagged subjects ---------------------------------
  int nsel = 0;
  for(int is = 0; is < nGroup; ++is) if(subj_flag[is] == 1) ++nsel;

  vector<Type> mu_sel( nsel * ngrid );
  int sel_idx = 0;

  for(int is = 0; is < nGroup; ++is){
    if(subj_flag[is] == 0) continue;

    Type b_intercept_subj = b_intercept(is);
    Type b_shift_subj = b_shift(is);

    for(int ig = 0; ig < ngrid; ++ig){
      Type ui = (age_grid(ig) + ga_subj(is)*beta_shift_ga + b_shift_subj);
      Type v =  (ui - a) / s;

      vector<Type> B = bspline_basis(v, knots, (int)degree);

      Type hval = Type(0.0);
      for(int k=0; k<Kint; ++k){
        B(k) -= Bmean(k);
        hval += c(k) * B(k);
      }

      Type mu_ij = beta_intercept + beta_intercept_sex*sex_subj(is)
                   + b_intercept_subj
                   + exp(beta_amplitude_sex*sex_subj(is))*hval;

      mu_sel[ sel_idx * ngrid + ig ] = mu_ij;
    }
    ++sel_idx;
  }

  REPORT(c);
  REPORT(a);
  REPORT(s);
  REPORT(rho);
  ADREPORT(rho);

  if(spline_ci == 1) ADREPORT(h_grid);
  if(spline_ci == 1) ADREPORT(mu_sel);
  return nll;
}
