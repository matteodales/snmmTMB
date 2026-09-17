// starting_points.cpp
// ---------------------------------------------------------------------------
// Penalized B-spline regression used to generate starting values.
//
// Model (no subject-level random effects):
//   y_i = beta1 + f(x_i) + eps_i,       eps_i ~ N(0, sigma^2)
//   f(x) = [B(x) - m]' c,               m_k = (1/N) sum_i B_k(x_i)
//   c    = Upos %*% vpos + U0 %*% gamma0
//
// The smooth is centered on the column means m of its own design matrix, so
// f is sum-to-zero over the observed x and the intercept beta1 is identified.
//
// This is step one of the two-step starting-value strategy
// of manuscript Section 2.3.2 ("Starting values"): this fitted
// penalized spline is then given to nlme() to obtain starting values
// for the fixed effects and variance components of the full SNMM. It is used by:
//   - Section 3, sine-curve simulation
//   - Section 3, bell-curve simulation
//   - Section 4, SMOCC application
//
//
// Data:
//   y          response vector (length N)
//   x          covariate values, on the domain of the knot vector.
//   knots      B-spline knot vector
//   degree     spline degree (cubic, degree = 3)
//   K          nominal number of basis functions, knots.size() - degree - 1;
//   Upos, dpos, U0   eigendecomposition of the 2nd-difference penalty matrix:
//                    Upos (K-1 x mpos) and dpos (mpos) span the positive
//                    eigenvalue space, U0 (K-1 x r) the null space
//   spline_ci  flag: 1 = evaluate the fitted curve on x_grid
//   x_grid     grid of x-values at which to report the fitted curve
//
// Parameters:
//   beta1      intercept
//   log_sigma  log residual standard deviation
//   vpos       penalized spline coordinates
//   gamma0     unpenalized spline coordinates
//   log_lambda log smoothing parameter
//
// Reported:
//   h_grid     ADREPORT, only if spline_ci == 1 - fitted curve on x_grid
//   c          ADREPORT - reconstructed spline coefficients
//   m          REPORT - the centering constants. The R callers read these
//              back with obj$report()$m to rebuild the centered
//              basis to pass to nlme()
//
// -------------------------------------------------------------------



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

  // scale-aware eps (tunable)
  Type knot_left = knots(0);
  Type knot_right = knots(nKnots - 1);
  Type span = CondExpGt(knot_right - knot_left, Type(0.0), knot_right - knot_left, Type(1.0));
  const double base_eps = 1e-8;
  Type eps = Type(base_eps) * span;

  // clamp into [knot_left, knot_right - eps]
  Type x = x_in;
  x = CondExpLt(x, knot_left, knot_left, x);                       // if x < left -> left else x
  x = CondExpGt(x, knot_right - eps, knot_right - eps, x);         // if x > right-eps -> right-eps else x

  // Degree 0 basis
  for(int i = 0; i < nBasis; ++i) {
    // indicator = (x >= knots[i]) && (x < knots[i+1])
    Type ge_left = CondExpGe(x, knots(i), Type(1.0), Type(0.0));
    Type lt_right = CondExpLt(x, knots(i+1), Type(1.0), Type(0.0));
    Type near_right = CondExpLe(abs(x - knot_right), eps, Type(1.0), Type(0.0));
    Type last_marker = CondExpGe(Type(i), Type(nBasis - 1), near_right, Type(0.0));
    Type indicator = ge_left * lt_right + last_marker;

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
  DATA_VECTOR(y);            // N
  DATA_VECTOR(x);            // on the domain spanned by knots

  DATA_VECTOR(knots);        // knot vector
  DATA_INTEGER(degree);      // spline degree
  DATA_INTEGER(K);           // number of basis functions (knots - degree -1)

  DATA_INTEGER(spline_ci);
  DATA_VECTOR(x_grid);

  // penalty eigendecomposition from R (Upos, dpos, U0)
  DATA_MATRIX(Upos);         // K x mpos
  DATA_VECTOR(dpos);         // length mpos
  DATA_MATRIX(U0);           // K x r

  // parameters 
  PARAMETER(beta1);
  PARAMETER(log_sigma);                // residual sd

  PARAMETER_VECTOR(vpos);              // penalized spline coordinates (random effects)
  PARAMETER_VECTOR(gamma0);            // unpenalized spline coordinates (fixed)
  PARAMETER(log_lambda);               // smoothing parameter

  int N = y.size();
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


  // Compute column-means
  vector<Type> m(Kint);
  for(int k=0;k<Kint;k++) m(k) = Type(0.0);

  for(int i=0;i<N;i++){

    vector<Type> B = bspline_basis(x(i), knots, (int)degree);
    for(int k=0;k<Kint;k++){
    m(k) += B(k);
    }
  }
  for(int k=0;k<Kint;k++) m(k) /= Type(N);

  int ngrid = x_grid.size();
  vector<Type> h_grid(ngrid);


  // compute on x_grid if requested
  if(spline_ci == 1){
    for(int ig = 0; ig < ngrid; ++ig){


      Type left = knots(0);
      Type right = knots(knots.size()-1);
      Type tiny_eps = Type(1e-12);
      if(x_grid(ig) < left) x_grid(ig) = left;
      if(x_grid(ig) > right - tiny_eps) x_grid(ig) = right - tiny_eps;

      vector<Type> B = bspline_basis(x_grid(ig), knots, (int)degree);
      
      Type hval = Type(0.0);
      for(int k=0;k<Kint;k++){
        B(k) -=  m(k);   // center column k
        hval += c(k) * B(k);
      }

      hval += beta1;
      h_grid(ig) = hval;
    }
  }





  // ---------- Likelihood evaluation ----------
  Type nll = Type(0.0);
  for(int i=0;i<N;i++){
    vector<Type> B = bspline_basis(x(i), knots, (int)degree);
    Type h = Type(0.0);
    for(int k=0;k<Kint;k++){
      B(k) -=  m(k);   // center column k
      h += c(k) * B(k);
    }

    Type mu = beta1 + h;

    nll -= dnorm(y(i), mu, sigma, true);
  }


  // ---------- prior for penalized spline coordinates vpos ----------
  if(dpos.size() > 0){
    int mpos = dpos.size();
    matrix<Type> Sigma_v(mpos, mpos);
    Sigma_v.setZero();
    for(int j=0;j<mpos;j++){
        Sigma_v(j,j) = Type(1.0) / ( lambda * dpos(j) );
    }
    density::MVNORM_t<Type> mvn(Sigma_v);
    nll += mvn(vpos);
  }
  
  if(spline_ci == 1) ADREPORT(h_grid);
  ADREPORT(c);
  REPORT(m);

  return nll;
}
