// snmmTMB_derivpenalty.cpp
// ---------------------------------------------------------------------------
// Penalized B-spline regression with a monotonicity penalty, for the
// illustration of Online Appendix D.
//
// Model (no subject-level random effects):
//   y_i = beta1 + f(x_i) + eps_i,       eps_i ~ N(0, sigma^2)
//   f(x) = [B(x) - Bmean]' c
//   c    = Upos %*% vpos + U0 %*% gamma0
//
// A penalty on the negative part of f' promotes a nondecreasing curve without
// constraining the spline coefficients directly. Writing
// f'(x) = sum_k c_k B'_k(x), the penalty is
//
//   lambda_c * sum_i [ min(f'(x_i), 0) ]^2,
//
// summed over the points of x_grid, with the nondifferentiable min replaced by
// the smooth approximation
//
//   min(f'(x), 0) ~ 0.5 * ( sqrt(f'(x)^2 + eps) - f'(x) ),   eps = 1e-8.
//
// The basis derivatives B'_k are evaluated in R with splineDesign(deriv = 1)
// and passed in as X_deriv. lambda_c is the penalty weight with
// larger values make the fitted curve less decreasing. lambda_c = 0 gives the
// unconstrained fit.
//
// The centering constants Bmean are computed once in R, as the column means of
// B over a fixed grid spanning the knot range, and passed in as data.
//
// Data:
//   y, x               response and covariate (x on the domain of the knots)
//   knots, degree, K   B-spline basis specification
//   Bmean              fixed-grid column means of the basis (length K-1)
//   Upos, dpos, U0     penalty eigendecomposition (see starting_points.cpp)
//   X_deriv            ngrid x (K-1) derivatives of the basis on x_grid
//   lambda_c           weight of the monotonicity penalty
//   spline_ci          flag: 1 = evaluate the fitted curve on x_grid
//   x_grid             grid of x-values, used both to report the curve and
//                      as the points at which f' is penalized
//
// Parameters:
//   beta1              intercept
//   log_sigma          log residual standard deviation
//   vpos, gamma0       spline coefficients
//   log_lambda         smoothing parameter
//
// Reported:
//   h_grid     ADREPORT - fitted curve on x_grid, including beta1
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

  // Cox-de Boor recursion
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
  DATA_VECTOR(y);            // response
  DATA_VECTOR(x);            // covariate, on the domain spanned by knots

  DATA_VECTOR(knots);        // knot vector
  DATA_INTEGER(degree);      // spline degree
  DATA_INTEGER(K);           // number of basis functions

  DATA_VECTOR(Bmean);        // fixed-grid centering constants (length K-1)

  DATA_INTEGER(spline_ci);
  DATA_VECTOR(x_grid);

  // penalty eigendecomposition
  DATA_MATRIX(Upos);
  DATA_VECTOR(dpos);
  DATA_MATRIX(U0);

  // monotonicity penalty
  DATA_MATRIX(X_deriv);      // ngrid x (K-1), basis derivatives on x_grid
  DATA_SCALAR(lambda_c);     // weight of the monotonicity penalty

  // parameters
  PARAMETER(beta1);                    // intercept
  PARAMETER(log_sigma);                // residual sd

  PARAMETER_VECTOR(vpos);              // penalized spline coordinates
  PARAMETER_VECTOR(gamma0);            // unpenalized spline coordinates
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


  // likelihood
  Type nll = Type(0.0);
  for(int i=0;i<N;i++){
    vector<Type> B = bspline_basis(x(i), knots, (int)degree);
    Type h = Type(0.0);
    for(int k=0;k<Kint;k++){
      B(k) -= Bmean(k);   // center column k
      h += c(k) * B(k);
    }

    Type mu = beta1 + h;

    nll -= dnorm(y(i), mu, sigma, true);
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

  // monotonicity penalty on the negative part of f', over x_grid
  int ngrid = x_grid.size();
  Type mono_pen = Type(0.0);
  Type eps = Type(1e-8);
  for(int ig=0;ig<ngrid;ig++){
    Type fprime = Type(0.0);
    for(int k=0;k<Kint;k++){
      fprime += X_deriv(ig,k) * c(k);
    }
    // smooth approximation of min(fprime, 0)
    Type neg = Type(0.5) * (sqrt(fprime*fprime + eps) - fprime);
    mono_pen += neg*neg;
  }

  nll += lambda_c * mono_pen;




  // Spline evaluation for the fitted curve
  vector<Type> h_grid(ngrid);

  if(spline_ci == 1){
    for(int ig = 0; ig < ngrid; ++ig){

      vector<Type> B = bspline_basis(x_grid(ig), knots, (int)degree);
      Type hval = Type(0.0);
      for(int k=0;k<Kint;k++){
        B(k) -= Bmean(k);   // center column k
        hval += c(k) * B(k);
      }
      h_grid(ig) = beta1 + hval;
    }
  }

  // ADREPORT for spline evaluation
  if(spline_ci == 1) ADREPORT(h_grid);


  // Return negative log-likelihood
  return nll;
}
