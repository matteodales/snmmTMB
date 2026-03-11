#include <TMB.hpp>







// x: the input value at which to evaluate the B-spline basis functions
// knots: the vector of knot positions for the B-spline basis
// degree: the degree of the B-spline basis functions
// returns: a vector of B-spline basis function values evaluated at x

// de Boor algorithm for evaluation of B-splines compatible with Automatic Differentiation.

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
  DATA_VECTOR(x);            // values in [0,1] (j/m)
  DATA_IVECTOR(group);       // group index 1 to nGroup
  DATA_INTEGER(nGroup);

  DATA_VECTOR(knots);        // knot vector
  DATA_INTEGER(degree);      // spline degree
  DATA_INTEGER(K);           // number of basis functions

  DATA_INTEGER(spline_ci);
  DATA_VECTOR(x_grid);
  DATA_IVECTOR(subj_flag);

  // penalty eigendecomposition
  DATA_MATRIX(Upos);       
  DATA_VECTOR(dpos);        
  DATA_MATRIX(U0);           

  // parameters 
  PARAMETER(alpha);                    // intercept
  PARAMETER_VECTOR(b1);                // random intercepts
  PARAMETER_VECTOR(b3);                // random phase
  PARAMETER(log_sd_b1);                // sd for b1
  PARAMETER(log_sd_b3);                // sd for b3
  PARAMETER(log_sigma);                // residual sd

  PARAMETER_VECTOR(vpos);              // penalized spline coordinates
  PARAMETER_VECTOR(gamma0);            // unpenalized spline coordinates
  PARAMETER(log_lambda);               // smoothing parameter

  int N = y.size();
  Type sd_b1 = exp(log_sd_b1);
  Type sd_b3 = exp(log_sd_b3);
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



  // shift data points according to random phase b3
  vector<Type> u(N);
  for(int i=0;i<N;i++){
    int g = group(i) - 1;
    u(i) = x(i) - b3(g);
  }



  // scale observed data to match the range of knots [0,1]
  Type x_min = x(0);
  Type x_max = x(0);
  for(int i = 1; i < N; ++i){
  if(x(i) < x_min) x_min = x(i);
  if(x(i) > x_max) x_max = x(i);
  }

  Type umin = x_min - Type(3) * sd_b3; // lower quantile
  Type umax = x_max + Type(3) * sd_b3; // upper quantile

  Type a = umin;
  Type s = (umax - umin);   // range




  // Compute column-means
  vector<Type> m(Kint);
  for(int k=0;k<Kint;k++) m(k) = Type(0.0);

  for(int i=0;i<N;i++){
    // map u to v
    Type v = (u(i) - a) / s;

    vector<Type> B = bspline_basis(v, knots, (int)degree);
    for(int k=0;k<Kint;k++) m(k) += B(k);
  }
  for(int k=0;k<Kint;k++) m(k) /= Type(N);


  // likelihood
  Type nll = Type(0.0);
  for(int i=0;i<N;i++){
    int g = group(i) - 1;
    Type v = (u(i) - a) / s;
    vector<Type> B = bspline_basis(v, knots, (int)degree);
    Type h = Type(0.0);
    for(int k=0;k<Kint;k++){
      B(k) -=  m(k);   // center column k
      h += c(k) * B(k);
    }

    Type mu = alpha + b1(g) + h;

    nll -= dnorm(y(i), mu, sigma, true);
  }

  // random-effects (b1, b3)
  for(int j=0;j<nGroup;j++){
    nll -= dnorm(b1(j), Type(0.0), sd_b1, true);
    nll -= dnorm(b3(j), Type(0.0), sd_b3, true);
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

  int ngrid = x_grid.size();
  vector<Type> h_grid(ngrid);

  if(spline_ci == 1){
    for(int ig = 0; ig < ngrid; ++ig){

      Type ui = x_grid(ig);
      Type v = (ui - a) / s;

      Type left = knots(0);
      Type right = knots(knots.size()-1);
      Type tiny_eps = Type(1e-12);
      if(v < left) v = left;
      if(v > right - tiny_eps) v = right - tiny_eps;

      vector<Type> B = bspline_basis(v, knots, (int)degree);
      Type hval = Type(0.0);
      for(int k=0;k<Kint;k++){
        B(k) -=  m(k);   // center column k
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
    Type b3_subj = b3(is);


    for(int ig = 0; ig < ngrid; ++ig){

      Type ui = x_grid(ig) - b3_subj;
      Type v = (ui - a) / s;


      vector<Type> B = bspline_basis(v, knots, (int)degree); 
      
      Type hval = Type(0.0);
      for(int k=0; k<Kint; ++k){
        B(k) -= m(k);            
        hval += c(k) * B(k); 
      }

      Type mu_ij = alpha + b1_subj + hval;


      mu_sel[ sel_idx * ngrid + ig ] = mu_ij;
    }
    ++sel_idx;
  } 

  // ADREPORT for spline evaluations
  if(spline_ci == 1) ADREPORT(h_grid);
  if(spline_ci == 1) ADREPORT(mu_sel);




  return nll;
}
