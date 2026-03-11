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
  DATA_VECTOR(y);            
  DATA_VECTOR(x);            
  DATA_SCALAR(x_min);
  DATA_SCALAR(x_max);
  DATA_SCALAR(ga_min);
  DATA_SCALAR(ga_max);

  DATA_VECTOR(knots);        
  DATA_INTEGER(degree);   
  DATA_INTEGER(K);           

  DATA_INTEGER(spline_ci);
  DATA_VECTOR(x_grid);
  DATA_IVECTOR(subj_flag);

  DATA_VECTOR(sex);
  DATA_VECTOR(ga);

  // penalty eigendecomposition
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
  PARAMETER_VECTOR(b_intercept);                // random intercepts (nGroup)
  PARAMETER_VECTOR(b_shift);                    // random phase (nGroup)
  PARAMETER(log_sd_b_intercept);                // sd for b_intercept
  PARAMETER(log_sd_b_shift);                    // sd for b_shift
  PARAMETER(log_sigma);                         // residual sd

  PARAMETER_VECTOR(vpos);                       // penalized spline coordinates (random effects)
  PARAMETER_VECTOR(gamma0);                     // unpenalized spline coordinates (fixed)
  PARAMETER(log_lambda);                        // smoothing parameter (controls variance of vpos)

  int N = y.size();
  Type sd_b_intercept = exp(log_sd_b_intercept);
  Type sd_b_shift = exp(log_sd_b_shift);
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

  //apply gamma
  vector<Type> u_all(N);
  for(int i=0;i<N;i++){
    int g = group(i) - 1;
    u_all(i) = (x(i) + ga(i)*beta_shift_ga + b_shift(g));
  }

  Type umin = (x_min  - Type(5) * sd_b_shift + ga_min*beta_shift_ga); // lower quantile
  Type umax = (x_max + Type(5) * sd_b_shift + ga_max*beta_shift_ga); // upper quantile

  Type a = umin;
  Type s = (umax - umin);   // range

  // Center columns to enforce sum-to-zero constraint
  vector<Type> m(Kint);
  for(int k=0;k<Kint;k++) m(k) = Type(0.0);

  for(int i=0;i<N;i++){

    Type v = (u_all(i) - a) / s;

    vector<Type> B = bspline_basis(v, knots, (int)degree);
    for(int k=0;k<Kint;k++){
    m(k) += B(k);
    }
  }

  for(int k=0;k<Kint;k++) m(k) /= Type(N);




  // likelihood
  Type nll = Type(0.0);
  for(int i=0;i<N;i++){
    int g = group(i) - 1;
    Type v = (u_all(i) - a) / s;
    vector<Type> B = bspline_basis(v, knots, (int)degree);
    Type h = Type(0.0);
    for(int k=0;k<Kint;k++){
      B(k) -=  m(k);
      h += c(k) * B(k);
    }

    Type mu = beta_intercept + beta_intercept_sex*sex(i) + b_intercept(g) + exp(beta_amplitude_sex*sex(i))*h;

    nll -= dnorm(y(i), mu, sigma, true);
  }


  // random effects
  for(int j=0;j<nGroup;j++){
    nll -= dnorm(b_intercept(j), Type(0.0), sd_b_intercept, true);
    nll -= dnorm(b_shift(j), Type(0.0), sd_b_shift, true);
  }


  // penalized spline random effects
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















  // Spline evaluation for population curve
  int ngrid = x_grid.size();
  vector<Type> h_grid(ngrid);

  if(spline_ci == 1){
    for(int ig = 0; ig < ngrid; ++ig){

      Type ui = (x_grid(ig));
      Type v = (ui- a) / s;

      Type left = knots(0);
      Type right = knots(knots.size()-1);
      Type tiny_eps = Type(1e-12);
      if(v < left) v = left;
      if(v > right - tiny_eps) v = right - tiny_eps;

      vector<Type> B = bspline_basis(v, knots, (int)degree);
      
      
      Type hval = Type(0.0);
      for(int k=0;k<Kint;k++){
        B(k) -=  m(k);
        hval += c(k) * B(k);
      }

      h_grid(ig) = beta_intercept + hval;
    }
  }






  // Spline evaluation for subject curves of flagged subjects
  int nsel = 0;
  for(int is = 0; is < nGroup; ++is) if(subj_flag[is] == 1) ++nsel;

  vector<Type> mu_sel( nsel * ngrid ); 
  int sel_idx = 0;

  for(int is = 0; is < nGroup; ++is){
    if(subj_flag[is] == 0) continue; 


    Type b_intercept_subj = b_intercept(is);
    Type b_shift_subj = b_shift(is);


    for(int ig = 0; ig < ngrid; ++ig){
      
      Type ui = (x_grid(ig) + ga(is)*beta_shift_ga + b_shift_subj);
      Type v =  (ui- a) / s;


      vector<Type> B = bspline_basis(v, knots, (int)degree);
      
      Type hval = Type(0.0);
      for(int k=0; k<Kint; ++k){
        B(k) -= m(k);  
        hval += c(k) * B(k);   
      }

      Type mu_ij = beta_intercept + beta_intercept_sex*sex(is) + b_intercept_subj + exp(beta_amplitude_sex*sex(is))*hval;


      mu_sel[ sel_idx * ngrid + ig ] = mu_ij;
    } 
    ++sel_idx;
  }




  REPORT(m);
  REPORT(c);
  REPORT(a);
  REPORT(s);

  if(spline_ci == 1) ADREPORT(h_grid);
  if(spline_ci == 1) ADREPORT(mu_sel);
  return nll;
}
