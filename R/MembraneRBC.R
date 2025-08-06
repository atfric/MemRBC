# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2025 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, S. (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# Frickenhaus S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#


#
#  global standard model parameters
#' @export
load_param_RBC<-function(msg)
{
  data(M.mu,M.C0,M.mu,M.Ka,M.K_b,M.K_ADE,M.Es,M.rho,M.a2,M.a3,M.a4,M.b0,M.b1,M.b2,M.rho,M.Rcpp,M.Rcpp_ncores,M.scr1,M.scr2,package = "MemRBC",envir = .GlobalEnv)
}

#
# M.xxx parameters to be set in main program
#

#' @export
citation.MemRBC<-function() {cat("when using this software for publications you must cite it as:\n Frickenhaus, S. (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340 \n")}


.onLoad <- function(libname, pkgname) {
  # startup message
  msg="MemRBC, the red blood cell shape modeling R package... \n... with compiled SCM energy, gradient and Hessian by Rcpp.\n Needs time for compiling..."

  if(!interactive()) {
    msg[1] <- paste("Package 'MemRBC'")
  }

  packageStartupMessage(msg)

{
if (!dir.exists(".Rc")) dir.create(".Rc")
Rcpp::sourceCpp(cacheDir = ".Rc", code='\\
#include <Rcpp.h>
#include<RcppEigen.h>
using namespace std;
using namespace Rcpp;

#include <omp.h>

//[[Rcpp::plugins(RcppEigen]]
//[[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]

// integration with Rcpp, see RcppCore,Bates/Eddelbuettel, https://github.com/RcppCore/RcppEigen/blob/master/README.md
using Eigen::Map;                       // maps rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision

void reshape(MatrixXd &x,unsigned int const r, unsigned int const c )	// slow:copies to a temp matrix, but does it Matlab style
{
  MatrixXd temp = x;
  x.resize(r,c);
  int counter = 0;
  int rp, cp;
  for(int cix = 0;cix<(temp.cols());cix++)
    for(int rix = 0;rix<temp.rows();rix++)
    {
      rp = int(counter%r);
      cp = int(counter%c);
      x(rp, cp) = temp(rix,cix);
      counter++;
    }
}

double factorial ( int n )
{
  int i;
  double value;
  value = 1.0;
  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }
return value;
}
double *pm_polynomial ( int mm, int n, int m, double x[] )
{//****************************************************************************80
  //
  //  Purpose:
  //
  //    PM_POLYNOMIAL evaluates the Legendre polynomials Pm(n,m,x).
  //
  //  Differential equation:
  //
  //    (1-X*X) * Y" - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
  //
  //  First terms:
  //
  //    M = 0  ( = Legendre polynomials of first kind P(N,X) )
  //
  //    Pm(0,0,x) =    1
  //    Pm(1,0,x) =    1 X
  //    Pm(2,0,x) = (  3 X^2 -   1)/2
  //    Pm(3,0,x) = (  5 X^3 -   3 X)/2
  //    Pm(4,0,x) = ( 35 X^4 -  30 X^2 +   3)/8
  //    Pm(5,0,x) = ( 63 X^5 -  70 X^3 +  15 X)/8
  //    Pm(6,0,x) = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
  //    Pm(7,0,x) = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
  //
  //    M = 1
  //
  //    Pm(0,1,x) =   0
  //    Pm(1,1,x) =   1 * SQRT(1-X^2)
  //    Pm(2,1,x) =   3 * SQRT(1-X^2) * X
  //    Pm(3,1,x) = 1.5 * SQRT(1-X^2) * (5*X^2-1)
  //    Pm(4,1,x) = 2.5 * SQRT(1-X^2) * (7*X^3-3*X)
  //
  //    M = 2
  //
  //    Pm(0,2,x) =   0
  //    Pm(1,2,x) =   0
  //    Pm(2,2,x) =   3 * (1-X^2)
  //    Pm(3,2,x) =  15 * (1-X^2) * X
  //    Pm(4,2,x) = 7.5 * (1-X^2) * (7*X^2-1)
  //
  //    M = 3
  //
  //    Pm(0,3,x) =   0
  //    Pm(1,3,x) =   0
  //    Pm(2,3,x) =   0
  //    Pm(3,3,x) =  15 * (1-X^2)^1.5
  //    Pm(4,3,x) = 105 * (1-X^2)^1.5 * X
  //
  //    M = 4
  //
  //    Pm(0,4,x) =   0
  //    Pm(1,4,x) =   0
  //    Pm(2,4,x) =   0
  //    Pm(3,4,x) =   0
  //    Pm(4,4,x) = 105 * (1-X^2)^2
  //
  //  Recursion:
  //
  //    if N < M:
  //      Pm(N,M,x) = 0
  //    if N = M:
  //      Pm(N,M,x) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
  //      all the odd integers less than or equal to N.
  //    if N = M+1:
  //      Pm(N,M,x) = X*(2*M+1)*Pm(M,M,x)
  //    if M+1 < N:
  //      Pm(N,M,x) = ( X*(2*N-1)*Pm(N-1,M,x) - (N+M-1)*Pm(N-2,M,x) )/(N-M)
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    14 March 2012
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Reference:
  //
  //    Milton Abramowitz, Irene Stegun,
  //    Handbook of Mathematical Functions,
  //    National Bureau of Standards, 1964,
  //    ISBN: 0-486-61272-4,
  //    LC: QA47.A34.
  //
  //  Parameters:
  //
  //    Input, int MM, the number of evaluation points.
  //
  //    Input, int N, the maximum first index of the Legendre
  //    function, which must be at least 0.
  //
  //    Input, int M, the second index of the Legendre function,
  //    which must be at least 0, and no greater than N.
  //
  //    Input, double X[MM], the point at which the function is to be
  //    evaluated.
  //
  //    Output, double PM_POLYNOMIAL[MM*(N+1)], the function values.
  //
  double fact;
  int i;
  int j;
  int k;
  double *v;
  v = new double[mm*(n+1)];
  for ( j = 0; j < n + 1; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = 0.0;
    }
  }
  //
  //  J = M is the first nonzero function.
  //
  if ( m <= n )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+m*mm] = 1.0;
    }
    fact = 1.0;
    for ( k = 0; k < m; k++ )
    {
      for ( i = 0; i < mm; i++ )
      {
        v[i+m*mm] = - v[i+m*mm] * fact * sqrt ( 1.0 - x[i] * x[i] );
      }
      fact = fact + 2.0;
    }
  }
  //
  //  J = M + 1 is the second nonzero function.
  //
  if ( m + 1 <= n )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+(m+1)*mm] = x[i] * ( double ) ( 2 * m + 1 ) * v[i+m*mm];
    }
  }
  //
  //  Now we use a three term recurrence.
  //
  for ( j = m + 2; j <= n; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = ( ( double ) ( 2 * j     - 1 ) * x[i] * v[i+(j-1)*mm]
                      + ( double ) (   - j - m + 1 ) *        v[i+(j-2)*mm] )
      / ( double ) (     j - m     );
    }
  }

  return v;
}

void legendre(int L, MatrixXd &ct, MatrixXd &pl)	//mimic Matlabs legendre function
{	// accepts a vector ct with elements between -1 and 1 (usually called as cos of the angle)
  // Output: pl of dimensions (Lo + 1) x (rowsxcols)
  // To do: pm_polynomial calculates other L values as well, which we throw away at the moment
  int dim = ct.cols()*ct.rows();
  if(ct.cols()>1 && ct.rows()>1) {reshape(ct,1,ct.rows()*ct.cols());}	// make sure ct and st are reshaped as row vectors
  if(ct.cols()==1){ct.transposeInPlace();}	// make sure ct is a row vector
  double* x = new double[dim];for (int i=0; i<dim; i++) { x[i] = ct(0,i);} // create and initialize to values in ct
  double* v ;//= new double[dim];for (int i=0; i<dim; i++) { v[i] = 0.0;}	// points to the array returned from
  pl.resize(L+1, ct.cols());		// make sure pl is the right size
  // loop over K values and fill in the output matrix pl
  for(int K=L;K>-1;K--)
  {
    int counter = 0;
    v = pm_polynomial (dim, L, K, x);	//call recursion formula
    //for (int i=dim; i<dim*2; i++) {std::cout<<v[i]<<std::endl;}
    for (int i=L*dim; i<dim*(L+1); i++) // we are only interested in the last section of v (see definition of output of pm_polynomial)
    {
      pl(K,counter) = v[i];			// fill in the values into our output matrix
      counter++;
    }
    delete []v;		// free the memory
  }
  delete []x;x = NULL;
}

void plkt(int L_max, MatrixXd &P, MatrixXd &P_T) // calculate the derivative of P w.r.t. theta and store in P_T
{
  Eigen::MatrixXd tmp, t1, t2;
  double fac1, fac2;
  int ix;
  P_T.resize(P.rows(), P.cols());
  P_T.fill(0.0);
  int ia = 1;
  Eigen::MatrixXd dpnm2 = P.col(1);
  for(int L = 0;L<=L_max;L++)
  {
    ia +=L;
    tmp = P.col(ia-1);
    P_T.col(ia-1) = (-1)*sqrt(double(ia-1))*P.col(ia);
    //if(true){std::cout<<"(-1)*sqrt(double(ia-1))*P.col(ia)"<<std::endl<<(-1)*sqrt(double(ia-1))*P.col(ia)<<std::endl;}
    fac1 = sqrt(double(2.0*L*(L+1)));
    for(int K = 1;K<=(L-1);K++)
    {
      ix = ia+K;
      fac2 = sqrt(double((L-K)*(L+K+1)));
      //if(true){std::cout<<"t1"<<std::endl<<fac1*tmp.array()<<std::endl;}
      //if(true){std::cout<<"t2"<<std::endl<<(fac2*P.col(ix).array())<<std::endl;}
      t1 = (fac1*tmp.array());
      t2 = (fac2*P.col(ix).array());
      P_T.col(ix-1) = 0.5*(t1-t2);
      //if(true){std::cout<<"0.5*(t1-t2)"<<std::endl<<0.5*(t1-t2)<<std::endl;}
      //if(true){std::cout<<"P_T --- up to Lmax = 1"<<std::endl<<"length: "<<P_T.rows()<<std::endl<<P_T<<std::endl;}
      tmp = P.col(ix-1);
      fac1 = fac2;
    }
    P_T.col(ia+L-1) = sqrt(double(L)/2)*tmp;
    //if(true){std::cout<<"P_T --- up to Lmax = 1"<<std::endl<<"length: "<<P_T.rows()<<std::endl<<P_T<<std::endl;}
  }
  P_T.col(2) = dpnm2;
  //if(true){std::cout<<"P_T --- up to Lmax = 1"<<std::endl<<"length: "<<P_T.rows()<<std::endl<<P_T<<std::endl;}
}
double N_LK_bosh(int L, int K)
{
  K = std::abs(K);
  if(K>L){ return 0;}
  else{ return std::sqrt((2-int(K==0))*(2*L+1)*factorial(L-K)/factorial(L+K));}
}

//[[Rcpp::export(invisible = true)]]
List L_Ylm(int L_max, Eigen::Map<Eigen::MatrixXd> t, Eigen::Map<Eigen::MatrixXd> p)
{
  MatrixXd ct, plk_mat, tmp, tmp2, tmp3,theta,phi;
  int counter = 0;
  int pcounter = 0;
  double NLK;
  int dim = t.size();
  MatrixXd YLK ;
  MatrixXd PLK ;
  theta=t;
  phi=p;
  double pdim = ((double(L_max)*double(L_max))/2) + (3*double(L_max)/2) + 1;
  PLK.resize(dim,int(pdim));  // will hold the associated Legendre function values
  ct = t.array().cos();				// precalculate cosine theta for passing to legendre
  YLK.resize(dim,(L_max+1)*(L_max+1));
  for(int L = 0;L<=L_max;L++)			// loop over the L values up to L_max
  {
    legendre(L, ct, plk_mat);				// works just like Matlabs legendre
    for(int K = -L;K<=L;K++)						// loop over K
    {
      NLK = N_LK_bosh(L,K);
      if(K >=0)
      {
        tmp = K*phi;
        tmp2 = tmp.array().cos();
        tmp3 = plk_mat.row(K);
        PLK.col(pcounter) = NLK*(tmp3.transpose());
        YLK.col(counter) = PLK.col(pcounter).array()*tmp2.array();
        pcounter++;
      }
      else if(K<0)
      {
        tmp = std::abs(K)*phi;
        tmp2 = tmp.array().sin();
        tmp3 = plk_mat.row(std::abs(K));
        YLK.col(counter) = NLK*tmp3.transpose().array() * tmp2.array();
      }
      counter++;
    }
  }
  return(List::create(Named("Ylm")=YLK,Named("PLK")=PLK));
}

// P = PLK from call of ylk_cos_sin_bosh
// generates a spherical harmonics basis according to the Bosh 2000 definition and normalization (as in custom Matlab code)
// Input: L_max, p and t are matrices of equal size
// Output:	YLK matrix dimensions [p.cols()xp.rows()] x (L_max + 1)^2 and holds the basis vectors
//			in the sequence 0,0   1,-1   1,0   1,1   2,-2   2,-1 ... etc.
//			PLK matrix of same dimensions as YLK, only holds the associated Legendre function values for use in derivative calculations
//[[Rcpp::export(invisible = true)]]
MatrixXd Ylm_v(int L_max,Eigen::Map<Eigen::MatrixXd> t,Eigen::Map<Eigen::MatrixXd> p,Eigen::Map<Eigen::MatrixXd> PP)
{
  MatrixXd tmp, tmp2, tmp3, phi, theta;
  MatrixXd Y_P;
  int counter = 0;
  int pindx;
  int dim = p.size();
  phi = p;
  reshape(phi,dim,1);
  theta = t;
  reshape(theta,dim,1);
  Y_P.resize(dim, (L_max + 1)*(L_max + 1));		// make sure Y_P has the right size
  double pdim = ((double(L_max)*double(L_max))/2) + (3*double(L_max)/2) + 1;
  PP.resize(dim,int(pdim));  // will hold the associated Legendre function values
  for(int L = 0;L<=L_max;L++)						// loop over the L values up to L_max
  {
    for(int K = -L;K<=L;K++)						// loop over K
    {
      if(K==0)
      {
        Y_P.col(counter) = Y_P.col(counter).array()*0.0;
      }
      else if(K>0)
      {
        pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
        tmp = (K*phi).array().sin();
        Y_P.col(counter) = PP.col(pindx).array()*(-1)*K*tmp.array();
      }
      else if(K<0)
      {
        pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
        tmp = std::abs(K)*phi;
        tmp2 = std::abs(K) * tmp.array().cos();
        tmp3 = PP.col(pindx);
        Y_P.col(counter) = tmp3.array() * tmp2.array();
      }
      counter++;
    }
  }
  return Y_P;
}


//[[Rcpp::export(invisible = true)]]
MatrixXd Ylm_vv(int L_max, Eigen::Map<Eigen::MatrixXd> t,Eigen::Map<Eigen::MatrixXd> p,Eigen::Map<Eigen::MatrixXd> PP)
{
  MatrixXd tmp, tmp2, tmp3, phi, theta, Y_PP;
  int counter = 0;
  int pindx;
  int dim = p.rows()*p.cols();
  phi = p;
  reshape(phi,dim,1);
  theta = t;
  reshape(theta,dim,1);
  Y_PP.resize(dim, (L_max + 1)*(L_max + 1));		// make sure Y_P has the right size
  for(int L = 0;L<=L_max;L++)						// loop over the L values up to L_max
  {
    for(int K = -L;K<=L;K++)						// loop over K
    {
      if(K==0)
      {
        Y_PP.col(counter) = Y_PP.col(counter).array()*0.0;
      }
      else if(K>0)
      {
        pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
        tmp = (K*phi).array().cos();
        Y_PP.col(counter) = PP.col(pindx).array()*(-1)*K*K*tmp.array();
      }
      else if(K<0)
      {
        pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
        tmp = std::abs(K)*phi;
        tmp2 = -K*K* tmp.array().sin();
        tmp3 = PP.col(pindx);
        Y_PP.col(counter) = tmp3.array() * tmp2.array();
      }
      counter++;
    }
  }
  return Y_PP;
}


//[[Rcpp::export(invisible = true)]]
List L_Ylm_u(int L_max, Eigen::Map<Eigen::MatrixXd> t,Eigen::Map<Eigen::MatrixXd> p,Eigen::Map<Eigen::MatrixXd> PP )
{
  MatrixXd tmp, tmp2, tmp3, phi, theta,Y_T, P_T, P;
  int counter = 0;
  int pindx;
  int dim = p.rows()*p.cols();
  P = PP;
  phi = p;
  reshape(phi,dim,1);
  theta = t;
  reshape(theta,dim,1);
  Y_T.resize(dim, (L_max + 1)*(L_max + 1) );		// make sure Y_T has the right size
  plkt(L_max, P, P_T);							// calculate derivative P_T of associated Legendre polynomials
  for(int L = 0;L<=L_max;L++)						// loop over the L values up to L_max
  {
    for(int K = -L;K<=L;K++)						// loop over K
    {
      if(K>=0)
      {
        pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
        tmp = (K*phi).array().cos();
        Y_T.col(counter) = P_T.col(pindx).array()*tmp.array();
      }
      else if(K<0)
      {
        pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
        tmp = std::abs(K)*phi;
        tmp2 = tmp.array().sin();
        tmp3 = P_T.col(pindx);
        Y_T.col(counter) = tmp3.array() * tmp2.array();
      }
      counter++;
    }
  }
  return(List::create(Named("Ylm_u")=Y_T,Named("P_T")=P_T));
}



//[[Rcpp::export(invisible = true)]]
MatrixXd Ylm_uv(int L_max, Eigen::Map<Eigen::MatrixXd> t,Eigen::Map<Eigen::MatrixXd> p,Eigen::Map<Eigen::MatrixXd> PP,Eigen::Map<MatrixXd> PP_T  )
{
  MatrixXd tmp, tmp2, tmp3, phi, theta, Y_TP, P, P_T;
  int counter = 0;
  int pindx;
  int dim = p.rows()*p.cols();
  phi = p;
  reshape(phi,dim,1);
  theta = t;
  reshape(theta,dim,1);
  P=PP;
  P_T=PP_T;
  Y_TP.resize(dim, (L_max + 1)*(L_max + 1));		// make sure Y_T has the right size
  for(int L = 0;L<=L_max;L++)						// loop over the L values up to L_max
  {
    for(int K = -L;K<=L;K++)						// loop over K
    {
      if(K>=0)
      {
        pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
        tmp = (K*phi).array().sin();
        Y_TP.col(counter) = P_T.col(pindx).array()*-K*tmp.array();
      }
      else if(K<0)
      {
        pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
        tmp = std::abs(K)*phi;
        tmp2 = std::abs(K) * tmp.array().cos();
        tmp3 = P_T.col(pindx);
        Y_TP.col(counter) = tmp3.array() * tmp2.array();
      }
      counter++;
    }
  }
  return Y_TP;
}


//[[Rcpp::export(invisible = true)]]
MatrixXd Ylm_uu(int L_max,Eigen::Map<Eigen::MatrixXd> t,Eigen::Map<Eigen::MatrixXd> p, Eigen::Map<Eigen::MatrixXd> PP_T  )
{
  MatrixXd tmp, tmp2, tmp3, phi, theta, P_T, P_TT, Y_TT;
  int counter = 0;
  int pindx;
  int dim = p.rows()*p.cols();
  phi = p;
  reshape(phi,dim,1);
  theta = t;
  reshape(theta,dim,1);
  Y_TT.resize(dim, (L_max + 1)*(L_max + 1));		// make sure Y_T has the right size // crashes when compiled as 32 bit with L_max 48 and 10242 points --- the matrix then is
  P_T=PP_T;
  plkt(L_max, P_T, P_TT);							// calculate derivative P_T of associated Legendre polynomials
  for(int L = 0;L<=L_max;L++)						// loop over the L values up to L_max
  {
    for(int K = -L;K<=L;K++)						// loop over K
    {
      if(K>=0)
      {
        pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
        tmp = (K*phi).array().cos();
        Y_TT.col(counter) = P_TT.col(pindx).array()*tmp.array();
      }
      else if(K<0)
      {
        pindx = int((double(L)-1)*(double(L)-1)/2+3*(double(L)-1)/2 + 1 + std::abs(double(K)));
        tmp = std::abs(K)*phi;
        tmp2 = tmp.array().sin();
        tmp3 = P_TT.col(pindx);
        Y_TT.col(counter) = tmp3.array() * tmp2.array();
      }
      counter++;
    }
  }
  return Y_TT;
}
//[[Rcpp::export]]
double dot2(NumericVector x, NumericVector y) {
  return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}

//[[Rcpp::export]]
NumericVector IntegM(NumericVector q, List grd,List bas)
{
  int Aimax=bas["Ai_max"];
  int Nuv=grd["ndof"];
  NumericVector wx=grd["wx"];
  NumericVector wy=grd["wy"];
  NumericVector I(Aimax*3); // result gradient to return, then CONV2()
  int i,j,k,l,m,n;
  int nx=wx.length();
  int ny=wy.length();
  double v,w;
  { // easy to include for-loop parallelism
    for (j=0;j<3;j++) {
      n=j*Aimax;// n: start in I(); l: start in q; wx x wy or wy x wx? # alternative l + m*ny + k // k*nx+m
      l=j*Aimax*Nuv;
#define ADR l+k*nx+m
      for(i=0;i<Aimax;i++,l+=Nuv) {
        for (v=0,m=0; m<nx; m++ ) {
          for (k=0,w=0; k<ny; k++) {w+= q(ADR)*wy(k);}
          v+=w*wx(m);
        }
        I(n)=v; n++;
      } // spectral mode i~lm
    } // j ~ X,Y,Z
  } //
  return I;
}

// [[Rcpp::export]]
double IntegS(NumericVector q, List grd)
{
  NumericVector wx = grd["wx"];
  NumericVector wy = grd["wy"];
  double I = 0; // result gradient to return, then CONV2()
  int k, m, l=0 ;
  int nx = wx.length();
  int ny = wy.length();
  double w;
  for (m=0; m<ny; m++ ) {
    for (k=0,w=0; k<nx; k++,l++) { w+= q(l)*wx(k); }
    I += w*wy(m);
  }
  return I;
}

// [[Rcpp::export]]
List E_SCM_cxx( NumericMatrix A, List grd, List bas, List C,
                double C0, double K_b, double K_ADE ) //# dbg=TRUE means no stop on NA
{ int i,k;
  int Nuv=grd["ndof"];// # number of gridpoints for integrands
  NumericMatrix Xum=C["X_u"];
  NumericVector Xu=Xum(_,0);
  NumericVector Yu=Xum(_,1);
  NumericVector Zu=Xum(_,2);
  NumericMatrix Xvm=C["X_v"];
  NumericVector Xv=Xvm(_,0);
  NumericVector Yv=Xvm(_,1);
  NumericVector Zv=Xvm(_,2);
  NumericVector E(Nuv) ;
  NumericVector FF(Nuv) ;
  NumericVector G(Nuv) ;
  NumericVector L(Nuv) ;
  NumericVector M(Nuv) ;
  NumericVector NN(Nuv) ;
  NumericVector dA(Nuv) ;
  NumericVector dV(Nuv) ;
  NumericVector curv(Nuv) ;
  NumericVector curv_sq(Nuv) ;
  NumericVector E_SCM_dens(Nuv);
  NumericMatrix normal(Nuv,3) ;
  NumericVector Nn(Nuv);
  // the normal does not match normal from R-code; unknown reason
  normal(_,0) = Yu *Zv;   normal(_,0) =normal(_,0) - Zu *Yv ;
  normal(_,1) = Zu *Xv;   normal(_,1) =normal(_,1) - Xu *Zv ;
  normal(_,2) = Xu *Yv;   normal(_,2) =normal(_,2) - Yu *Xv ;
  NumericMatrix n = normal ;// # tb normalized
  for (i=0;i<Nuv;i++) Nn(i) = sqrt(dot2(normal(i,_),normal(i,_))) ;
  //    if (sum(Nn==0)>0) stop("ERROR: zeros in ||normal||")
  NumericVector iNn=1/Nn ;
  E  = Xu*Xu + Yu*Yu + Zu*Zu ;
  FF = Xu*Xv + Yu*Yv + Zu*Zv ;
  G  = Xv*Xv + Yv*Yv + Zv*Zv ;
  dA = sqrt(E*G-FF*FF) ;
  double Area = IntegS(dA,grd) ;
  NumericVector inn=1/dA ;
  NumericVector inn2=inn*inn ;
  NumericMatrix X=C["X"] ;
  for (k=0;k<3;k++) n(_,k) = normal(_,k) * iNn ;
  dV =  1/3.0*(X(_,0)*n(_,0)+X(_,1)*n(_,1)+X(_,2)*n(_,2)) * dA ; //# dA for integration
  double Volume = IntegS(dV,grd);
  NumericMatrix Xuu=C["X_uu"];
  NumericMatrix Xuv=C["X_uv"];
  NumericMatrix Xvv=C["X_vv"];
  L=  Xuu(_,0)*n(_,0) + Xuu(_,1)*n(_,1) + Xuu(_,2)*n(_,2);
  M=  Xuv(_,0)*n(_,0) + Xuv(_,1)*n(_,1) + Xuv(_,2)*n(_,2);
  NN= Xvv(_,0)*n(_,0) + Xvv(_,1)*n(_,1) + Xvv(_,2)*n(_,2);
  curv= -(E*NN+G*L-2.0*FF*M) * inn ; //# work with positive curvature; H=(EN+GL-2FM)^2/dA^2; one dA canceled with dA=sqrt(...) du dv
  //    curv[is.na(curv)]<-0# warning("one or more Inf in curv!")
  double Curv =  IntegS(curv,grd);
  curv_sq = curv*curv*inn;
  //    curv_sq[is.na(curv_sq)]<-0# warning("one or more Inf in 2curv_sq!")
  double H2_BC = IntegS(curv_sq,grd);
  double H2 = K_b/2.0 * (H2_BC - 2.0*C0*Curv) + //# + K_b/2*C0^2*Area + # constant terms out
    + K_ADE*(Curv*Curv ) /Area;
    E_SCM_dens = K_b/2.0 * (curv_sq - 2.0*C0*curv) + // # + K_b/2*C0^2*Area + //# constant terms out
    + K_ADE*curv*Curv  / Area;
    List  retL=List::create(Named("Wb")=H2, _["H2_BC"]=H2_BC,
                            _["Area"]=Area, _["Volume"]=Volume, _["Curv"]=Curv,
                            _["E"]=E, _["FF"]=FF, _["G"]=G,
                            _["L"]=L, _["M"]=M, _["NN"]=NN, _["Nn"]=Nn,
                            _["dA"]=dA, _["dV"]=dV, _["curv"]=curv, _["curv_sq"]=curv_sq,
                            _["normal"]=normal, _["inn"]=inn, _["inn2"]=inn2, _["n"]=n,
                            _["E_SCM_dens"]=E_SCM_dens);
    return(retL);
} // E_SCM_cxx


//[[Rcpp::export]]
List Grad_SCM_cxx(
    List h2, List grd, List B, List CC, double C0,
    int ncores, double K, double KADE)
{
  NumericMatrix X_u=CC["X_u"];
  NumericMatrix X_v=CC["X_v"];
  NumericMatrix X_vv=CC["X_vv"];
  NumericMatrix X_uv=CC["X_uv"];
  NumericMatrix X_uu=CC["X_uu"];
  NumericMatrix Ylm=B["Ylm"];
  NumericMatrix Ylm_u=B["Ylm_u"];
  NumericMatrix Ylm_v=B["Ylm_v"];
  NumericMatrix Ylm_uv=B["Ylm_uv"];
  NumericMatrix Ylm_uu=B["Ylm_uu"];
  NumericMatrix Ylm_vv=B["Ylm_vv"];
  NumericMatrix normal=h2["normal"];
  NumericVector inn=h2["inn"];
  NumericVector inn2=h2["inn2"];
  NumericVector E=h2["E"];
  NumericVector FF=h2["FF"];
  NumericVector G=h2["G"];
  NumericVector L=h2["L"];
  NumericVector M=h2["M"];
  NumericVector NN=h2["NN"];
  NumericVector curv=h2["curv"];
  NumericVector curv_sq=h2["curv_sq"];
  double Area=h2["Area"];
  double Curv=h2["Curv"];
  int Nuv=grd["ndof"];
  int Ai_max=B["Ai_max"];
  NumericVector wx=grd["wx"];
  NumericVector wy=grd["wy"];
  NumericVector  dNaz(Nuv*Ai_max*3); // Nuv to be integrated out
  NumericVector  dNay(Nuv*Ai_max*3);
  NumericVector  dNax(Nuv*Ai_max*3);
  NumericVector    dE(Nuv*Ai_max*3);
  NumericVector    dF(Nuv*Ai_max*3);
  NumericVector    dG(Nuv*Ai_max*3);
  NumericVector    dL(Nuv*Ai_max*3);
  NumericVector   dNN(Nuv*Ai_max*3);
  NumericVector    dM(Nuv*Ai_max*3);
  NumericVector      ddV (Nuv*Ai_max*3);
  NumericVector      ddA (Nuv*Ai_max*3);
  NumericVector    dcurv (Nuv*Ai_max*3) ;
  NumericVector dcurv_sq (Nuv*Ai_max*3);
  NumericVector gradH2 (Ai_max*3); // Nuv integrated out in these:
  NumericVector gradH2BC (Ai_max*3);
  NumericVector gradC (Ai_max*3);
  NumericVector gradV (Ai_max*3);
  NumericVector gradA (Ai_max*3);

  int offs=Nuv*Ai_max;
  int offs2=2*offs;
  int i,j,k,l;
#pragma omp parallel num_threads(ncores) private(l,j,i,k)
{
#pragma omp for
  for (j=0;j<Ai_max;j++) {l=j*Nuv;for (int i=0;i<Nuv;i++,l++) {
    dNay(l) = Ylm_u(i,j)*X_v(i,2) - Ylm_v(i,j)*X_u(i,2);
    dNaz(l) = Ylm_v(i,j)*X_u(i,1) - Ylm_u(i,j)*X_v(i,1);}
  }
#pragma omp for
  for (j=0;j<Ai_max;j++) {l=j*Nuv;for (int i=0;i<Nuv;i++,l++) {
    dNax(l+offs) = Ylm_v(i,j)*X_u(i,2) - Ylm_u(i,j)*X_v(i,2);
    dNaz(l+offs) = Ylm_u(i,j)*X_v(i,0) - Ylm_v(i,j)*X_u(i,0);}
  }
#pragma omp for
  for (j=0;j<Ai_max;j++) {l=j*Nuv;for (int i=0;i<Nuv;i++,l++) {
    dNax(l+offs2) = Ylm_u(i,j)*X_v(i,1) - Ylm_v(i,j)*X_u(i,1);
    dNay(l+offs2) = Ylm_v(i,j)*X_u(i,0) - Ylm_u(i,j)*X_v(i,0);}
  }
#pragma omp for
  for(k=0;k<3;k++) {l=k*offs;for (j=0;j<Ai_max;j++) for (i=0;i<Nuv;i++,l++)
    dE(l) = 2 * X_u(i,k) * Ylm_u(i,j) ;
  }
#pragma omp for
  for(k=0;k<3;k++) {l=k*offs; for (j=0;j<Ai_max;j++) for (i=0;i<Nuv;i++,l++)
    dF(l) =     X_u(i,k) * Ylm_v(i,j) + X_v(i,k) * Ylm_u(i,j) ;
  }
#pragma omp for
  for(k=0;k<3;k++) {l=k*offs;for (j=0;j<Ai_max;j++) for (i=0;i<Nuv;i++,l++)
    dG(l) = 2 * X_v(i,k) * Ylm_v(i,j) ;
  }
#pragma omp for
  for (k=0;k<3;k++) {l=k*offs;for(j=0;j<Ai_max;j++) for (i=0;i<Nuv;i++,l++)
    ddV(l) = normal(i,k)*Ylm(i,j);
  }
  {
#pragma omp barrier
  }
#pragma omp for
  for(k=0;k<3;k++) {l=k*offs;for (j=0;j<Ai_max;j++) for (i=0;i<Nuv;i++,l++)
    ddA(l) = 0.5 * ( dE(l) *G(i)  + E(i)*dG(l) - 2*FF(i)*dF(l) ) * inn(i) ;
  }
#pragma omp for
  for (j=0;j<Ai_max;j++) {l=j*Nuv;for (i=0;i<Nuv;i++,l++)
    dL(l) =  X_uu(i,1)*dNax(l+offs)+X_uu(i,2)*dNax(l+offs2) + Ylm_uu(i,j)*normal(i,0);
  }
#pragma omp for
  for (j=0;j<Ai_max;j++) {l=j*Nuv;for (i=0;i<Nuv;i++,l++)
    dL(l+offs) =  X_uu(i,0)*dNay(l)+X_uu(i,2)*dNay(l+offs2) + Ylm_uu(i,j)*normal(i,1);
  }
#pragma omp for
  for (j=0;j<Ai_max;j++) {l=j*Nuv;for (i=0;i<Nuv;i++,l++)
    dL(l+offs2) =  X_uu(i,0)*dNaz(l)+X_uu(i,1)*dNaz(l+offs)  + Ylm_uu(i,j)*normal(i,2);
  }
#pragma omp for
  for (j=0;j<Ai_max;j++) {l=j*Nuv;for (i=0;i<Nuv;i++,l++)
    dM(l) = X_uv(i,1)*dNax(l+offs)+X_uv(i,2)*dNax(l+offs2) + Ylm_uv(i,j)*normal(i,0);
  }
#pragma omp for
  for (j=0;j<Ai_max;j++) {l=j*Nuv;for (i=0;i<Nuv;i++,l++)
    dM(l+offs) = X_uv(i,0)*dNay(l)+ X_uv(i,2)*dNay(l+offs2) + Ylm_uv(i,j)*normal(i,1);
  }
#pragma omp for
  for (j=0;j<Ai_max;j++) {l=j*Nuv;for (i=0;i<Nuv;i++,l++)
    dM(l+offs2) = X_uv(i,0)*dNaz(l)+X_uv(i,1)*dNaz(l+offs) + Ylm_uv(i,j)*normal(i,2);
  }
#pragma omp for
  for (j=0;j<Ai_max;j++) {l=j*Nuv;for (i=0;i<Nuv;i++,l++)
    dNN(l) = X_vv(i,1)*dNax(l+offs)+X_vv(i,2)*dNax(l+offs2)+ Ylm_vv(i,j)*normal(i,0);
  }
#pragma omp for
  for (j=0;j<Ai_max;j++) {l=j*Nuv; for (i=0;i<Nuv;i++,l++)
    dNN(l+offs) = X_vv(i,0)*dNay(l)     +X_vv(i,2)*dNay(l+offs2)+ Ylm_vv(i,j)*normal(i,1);
  }
#pragma omp for
  for (j=0;j<Ai_max;j++) {l=j*Nuv;for (i=0;i<Nuv;i++,l++)
    dNN(l+offs2) = X_vv(i,0)*dNaz(l)     +X_vv(i,1)*dNaz(l+offs)+ Ylm_vv(i,j)*normal(i,2);
  }
  {
#pragma omp barrier
  }
#pragma omp for schedule (static)
  for(k=0;k<3;k++) {l=k*offs; for(j=0;j<Ai_max;j++) for (i=0;i<Nuv;i++,l++)
    dL(l) = dL(l)* inn(i) -  L[i]*inn[i]*ddA(l);
  }
#pragma omp for schedule (static)
  for(k=0;k<3;k++) {l=k*offs; for(j=0;j<Ai_max;j++) for (i=0;i<Nuv;i++,l++)
    dM(l) = dM(l)* inn(i) -  M[i]*inn[i]*ddA(l);
  }
#pragma omp for schedule (static)
  for(k=0;k<3;k++) {l=k*offs; for(j=0;j<Ai_max;j++) for (i=0;i<Nuv;i++,l++)
    dNN(l)=dNN(l)* inn(i) -  NN[i]*inn[i]*ddA(l);
  }
#pragma omp for schedule (static)
  for (k=0;k<3;k++) {l=k*offs;for (j=0;j<Ai_max;j++)  for (i=0;i<Nuv;i++,l++)
    dcurv(l) =  -(dE(l)*NN(i) + E(i)*dNN(l)  + dG(l)*L(i) + G(i)*dL(l)
                    - 2*dF(l)*M(i) - 2*FF[i]*dM(l)) * inn(i) +
                      - curv(i) * ddA(l) * inn(i)  ;
  }  // dcurv now ready for use in dcurv_sq
  {
#pragma omp barrier
  }
#pragma omp for schedule (static)
  for (k=0;k<3;k++) {l=k*offs; for (j=0;j<Ai_max;j++)  for (i=0;i<Nuv;i++,l++)
    dcurv_sq(l) =  ( 2*curv(i)*dcurv(l) * inn(i)
                       -  curv_sq(i)*ddA(l)*inn(i) );
  }
  {
#pragma omp barrier
  }
}//parallel
// the following 4 lines failed as parallel sections
//  due to stack imbalances in R
// #pragma omp parallel num_threads(4)
{
// #pragma omp sections
{
// #pragma omp section
  gradH2BC = IntegM(dcurv_sq, grd, B);
// #pragma omp section
  gradC = IntegM(dcurv, grd,B);
// #pragma omp section
  gradV = IntegM(ddV, grd,B);
// #pragma omp section
  gradA = IntegM(ddA, grd,B);
}
}

gradH2 = K/2 * (gradH2BC  - 2*C0*gradC )  + KADE * (2 * Curv * gradC / Area - gradA * Curv*Curv / (Area*Area) );
Function CONV2("vectomat_cxx");
Function CONV3("vectoarr_cxx");
List Ret= List::create(
  Named("grad_H2BC")=CONV2(gradH2BC,Ai_max),
  _["grad_SCM"] =CONV2(gradH2,Ai_max),
  _["gradV"] =CONV2(gradV,Ai_max),
  _["gradA"] =CONV2(gradA,Ai_max),
  _["gradC"] =CONV2(gradC,Ai_max),
  _["ddV"] =CONV3(ddV,Nuv,Ai_max),
  _["ddA"]=CONV3(ddA,Nuv,Ai_max),
  _["dE"]=CONV3(dE,Nuv,Ai_max),
  _["dF"]=CONV3(dF,Nuv,Ai_max),
  _["dG"]=CONV3(dG,Nuv,Ai_max)
);
return( Ret );
}; // Grad_SCM_cxx

NumericVector MatVec_cxx(NumericMatrix M, NumericVector V)
{
  Environment base("package:base");
  Function matvec = base["%*%"];
  return(matvec(M,V));
}

NumericMatrix MatMat_cxx(NumericMatrix M, NumericMatrix N)
{
  Environment base("package:base");
  Function matmat = base["%*%"];
  return(matmat(M,N));
}

// same low speed like with etest inlining
//[[Rcpp::export]]
MatrixXd matmatE(Eigen::Map<Eigen::MatrixXd> tm, Eigen::Map<Eigen::MatrixXd> tm2)
{
  Eigen::MatrixXd prod = tm*tm2;
  return(prod);
}

// not needed like this - use Function updateX_only("updateX_only") in cpp-codes
//   if needed
// however, faster code by faster mat-mat could be interesting

//[[Rcpp::export]]
List updateX_only_cxx(NumericMatrix A, List grd, List bas)
{
 Function etest("etest"); //for etest see: https://stackoverflow.com/questions/37191673/matrix-multiplication-in-rcpp
 List C;
 NumericMatrix Y=bas["Ylm"];
 NumericMatrix X=etest(Y,A);
 // this:
// NumericMatrix X=MatMat_cxx(bas["Ylm"],A);
// C["X"]=X;
// or this:
//  NumericMatrix X(grd["ndof"],3);
//  int i,j;
//  int imax=grd["ndof"];
//  int Aimax=bas["Ai_max"];
//  NumericMatrix Y=bas["Ylm"];
//  Y=transpose(Y);
//  for (i=0;i<imax;i++){
//      X(i,0) = sum(Y(_,i)*A(_,0));
//      X(i,1) = sum(Y(_,i)*A(_,1));
//      X(i,2) = sum(Y(_,i)*A(_,2));    }
// loop code needs 6 sec;
//  etest takes 3 sec. vs 2.7 with R "%*%"
//
  C["X"]=X;
  return(C);
}

//[[Rcpp::export]]
List Hessian_SCM_SEN_cxx(NumericMatrix A,List grd, List bas, List Ref, double del, int ncores, double C0, double K_b, double K_ADE)
{ int Aimax=bas["Ai_max"];
  NumericMatrix H(Aimax*3,Aimax*3);
  // NumericMatrix H1(Aimax*3,Aimax*3);
  int i,j,k,l;
  NumericMatrix A0(Aimax,3);
  NumericMatrix A1(Aimax,3);
  Function updateX("updateX");
  Function SEN("SEN");
  Function E_SCM("E_SCM");
  Function E_SEN("E_SEN");
  Function Grad_SEN("Grad_SEN");
  Function Grad_SCM("Grad_SCM");
  Function synth12("synth12"); //    function(A,C,i,j,k) # spatial k
  Function mat2vec("mat2vec");
  Function matdiff2vec("matdiff2vec");
  Function matadd2vec("matadd2vec");
  Function symmetrize("symmetrize");
  List C=updateX(A,grd,bas);
  List h20=E_SCM(A,grd,bas,C);
  List  S0=SEN(A,grd,bas,Ref,h20);
  // double  ES=E_SEN(A,grd,bas,S,Ref);
  List Gh20=Grad_SCM(h20,grd,bas,C);
  List GS0=Grad_SEN(A,grd,bas,Gh20,S0,Ref);
//  NumericVector G0 = matadd2vec(Gh20["grad_SCM"] , GS0["grad_SEN"]);
  NumericVector G0 = NumericVector(Gh20["grad_SCM"])+NumericVector(GS0["grad_SEN"]);
  A0=A;
  l=0;
  List G;
  List h2;
  List Gh2;
  List GS;
  List S;
  NumericVector D;
  NumericVector G1;
//
// #define PARHESS
// calling R-functions seems to be not thread safe, crashes in the following
//
#ifdef PARHESS
#pragma omp parallel num_threads(ncores) private(k,l,i,A,C,h2,S,Gh2,GS,G,D,G1)
#endif
{
#ifdef PARHESS
#pragma omp for
#endif
  for (k=0;k<3;k++){l=k*Aimax;
    for (i=0;i<Aimax;i++,l++)
    {
      A(i,k)=A(i,k)+del;
      C=updateX(A,grd,bas);
      A(i,k)=A(i,k)-del;
      // sweeping through C$X_uv instead of full update
      //  C=synth12(A,C,l,l+1, i, i+1,del); // better inline this; use R-numbering
      //  if (i>0) C$X[,k] = C$X[,k] - bas$Ylm[,i]*del # remove del term
      //  if(j<dim(A)[1]) C$X[,k] = C$X[,k] + bas$Ylm[,j]*del
      h2=E_SCM_cxx(A,grd,bas,C,C0,K_b,K_ADE);
      S=SEN(A,grd,bas,Ref,h2);
      // double  ES=E_SEN(A,grd,bas,S,Ref);
      Gh2=Grad_SCM_cxx(h2,grd,bas,C,C0,ncores,K_b,K_ADE);
      GS=Grad_SEN(A,grd,bas,Gh2,S,Ref);
      G1 = NumericVector(Gh2["grad_SCM"]) + NumericVector(GS["grad_SEN"]);
      D = NumericVector(G1-G0);
      H(_,l)=D;
    }
  }
} // parallel

return(List::create(Named("H")=H,_("G")=G0,
                    _("gradA")=Gh20["gradA"],_("gradV")=Gh20["gradV"],_("gradC")=Gh20["gradC"],
                    _["g2"]=Gh20,_["h2"]=h20));
}

//[[Rcpp::export]]
List Hessian_SCM_cxx(NumericMatrix A,List grd, List bas, List Ref, double del, int ncores, double C0, double K_b, double K_ADE)
{ int Aimax=bas["Ai_max"];
  NumericMatrix H(Aimax*3,Aimax*3);
  // NumericMatrix H1(Aimax*3,Aimax*3);
  int i,j,k,l;
  NumericMatrix A0(Aimax,3);
  NumericMatrix A1(Aimax,3);
  Function updateX("updateX");
  Function SEN("SEN");
  Function E_SCM("E_SCM");
  Function E_SEN("E_SEN");
  Function Grad_SEN("Grad_SEN");
  Function Grad_SCM("Grad_SCM");
  Function synth12("synth12"); //    function(A,C,i,j,k) # spatial k
  Function mat2vec("mat2vec");
  Function matdiff2vec("matdiff2vec");
  Function matadd2vec("matadd2vec");
  Function symmetrize("symmetrize");
  List C=updateX(A,grd,bas);
  List h20=E_SCM(A,grd,bas,C);
//  List  S0=SEN(A,grd,bas,Ref,h20);
  // double  ES=E_SEN(A,grd,bas,S,Ref);
  List Gh20=Grad_SCM(h20,grd,bas,C);
//  List GS0=Grad_SEN(A,grd,bas,Gh20,S0,Ref);
  NumericVector G0 = NumericVector(Gh20["grad_SCM"]); // +NumericVector(GS0["grad_SEN"]);
  A0=A;
  l=0;
  List G;
  List h2;
  List Gh2;
  List GS;
  List S;
  NumericVector D;
  NumericVector G1;
//
// #define PARHESS 1
// calling R-functions is not thread safe in the following
//
#ifdef PARHESS
#pragma omp parallel num_threads(ncores) private(k,l,i,A,C,h2,S,Gh2,GS,G,D,G1)
#endif
{
#ifdef PARHESS
#pragma omp for
#endif
  for (k=0;k<3;k++){l=k*Aimax;
    for (i=0;i<Aimax;i++,l++)
    {
      A(i,k)=A(i,k)+del;
      C=updateX(A,grd,bas);
      A(i,k)=A(i,k)-del;
      // sweeping through C$X_uv instead of full update
      //  C=synth12(A,C,l,l+1, i, i+1,del); // better inline this; use R-numbering
      //  if (i>0) C$X[,k] = C$X[,k] - bas$Ylm[,i]*del # remove del term
      //  if(j<dim(A)[1]) C$X[,k] = C$X[,k] + bas$Ylm[,j]*del
      h2=E_SCM_cxx(A,grd,bas,C,C0,K_b,K_ADE);
 //     S=SEN(A,grd,bas,Ref,h2);
      // double  ES=E_SEN(A,grd,bas,S,Ref);
      Gh2=Grad_SCM_cxx(h2,grd,bas,C,C0,ncores,K_b,K_ADE);
  //    GS=Grad_SEN(A,grd,bas,Gh2,S,Ref);
      G1 = NumericVector(Gh2["grad_SCM"]) ;//+ NumericVector(GS["grad_SEN"]);
      D = NumericVector(G1-G0);
      H(_,l)=D;
    }
  }
} // parallel
return(List::create(Named("H")=H,_("G")=G0,
                  _("gradA")=Gh20["gradA"],
                  _("gradV")=Gh20["gradV"],
                  _("gradC")=Gh20["gradC"],
                  _["g2"]=Gh20,_["h2"]=h20));
}
')
} # sourceCPP

citation.MemRBC();
utils::data(M.mu,M.C0,M.mu,M.Ka,M.K_b,M.K_ADE,M.Es,M.rho,M.a2,M.a3,M.a4,M.b0,M.b1,M.b2,M.rho,M.Rcpp,M.Rcpp_ncores,package = "MemRBC",envir = .GlobalEnv)
M.scr1=M.scr2=-1
}

HAVE_DEPRECATED=FALSE

M.TEST=FALSE
# main code must define M.TEST before sourcing

# for tests: show severe deviation from zero as error
#' @export
severe<-function(q, what="some test", tol=1e-12)
{ testthat::test_that(what,{testthat::expect_equal(q,0)})}
#{cat(q,":");if (all(abs(q)>tol)) {stop("severe imprecision, STOP")} else
#  cat(crayon::green("OK:",what,"\n"))
#}

#
# Functionals

# Spontaneous Curvature Model SCM
#  Bending Energy

#   functionals to compute from precomputed coordinates X and derivatives X.. stored in object C,
#   by spectral approach
#
#   + utilities for curvature and more
#
#
# compute energy and other quantities
# with pre-Updated coords and derivatives in the object often denoted C
#

#
#   no curvature-constraint, instead:
#  -- Helfrich spontaneous curvature with parameter H_0 in global variable C0
#     (note:   C0 = 2*H0_Surface_Evolver )
#
# and ADE, but delta_A_0=0 (delta_A_0 not implemented)
# H := c1 + c2
# W = K_b/2*int[(2H-C0)^2] + K_ADE / A * (int[2H])^2

#
# for Helfrich spontaneous curvature model
#  where W = int (H-C0)^2 da
#  and w = W reduced reported relative to unit sphere
#  and H = c1 + c2 = 2 * H_mean ; H_mean in the usual sense like in Rvcg
#

# Bending energy , returns list with "Wb" as bending energy
#' @export
E_SCM <- function(A,grd,bas,C,plt=FALSE,dbg=FALSE,clp=FALSE) # dbg=TRUE means no stop on NA
{ if(!M.Rcpp) return(E_SCM_R(A,grd,bas,C,plt=FALSE,dbg=FALSE,clp=FALSE))
  return(E_SCM_cxx(A,grd,bas,C,M.C0,M.K_b,M.K_ADE)) # no plotting in cxx
}

# Bending energy , returns list with "Wb" as bending energy
#' @export
E_SCM_R <- function(A,grd,bas,C,plt=FALSE,dbg=FALSE,clp=FALSE) # dbg=TRUE means no stop on NA
{ # returns some quantities and fundamentals in a list
  # reference values for unit sphere
  int2d_s=int2d_scalar_GLS

  Nuv<-grd$ndof # number of gridpoints for integrands
  # copy for easier coding - updateX()->C must have been called before
  Xu<-C$X_u[,1]; Yu<-C$X_u[,2]; Zu<-C$X_u[,3];
  Xv<-C$X_v[,1]; Yv<-C$X_v[,2]; Zv<-C$X_v[,3];
  E<-FF<-G<-L<-M<-NN<-dA<-dV<-curv<-curv_sq<-rep(1.0,Nuv) #fundamentals/ integrands
  normal<-matrix(1.0,grd$ndof,3) # non-normalized surface normals!
  # not faster: for (i in 1:grd$ndof) normal[i,]<-cross(C$X_u[i,],C$X_v[i,])
  normal[,1] <- Yu *Zv  - Zu *Yv
  normal[,2] <- Zu *Xv  - Xu *Zv
  normal[,3] <- Xu *Yv  - Yu *Xv
  n<-normal # tb normalized in n hereafter
  #  if(!dbg & any(is.na(n))) stop("one or more NA in normal!")
  # normalization of normals / some may be 0 (at poles, where dv=0)
  Nn <- apply(normal,1,pracma::Norm)
  if (sum(Nn==0)>0) stop("ERROR: zeros in ||normal||")
  iNn<-1/Nn
  #  iNn[is.infinite(iNn)]<-0 # a hack for problematic points; suppress by * inn[i]
  # normalize to unit normals
  for (k in 1:3) n[,k]<-n[,k]*iNn # second fundamental needs unit n
  #n[is.na(n)]<-0 # replace NAs
  #n[is.infinite(n)]<-0 # replace Infs
  # fundamental form I
  E  <- Xu*Xu + Yu*Yu + Zu*Zu
  FF <- Xu*Xv + Yu*Yv + Zu*Zv
  G  <- Xv*Xv + Yv*Yv + Zv*Zv
  #  if (!dbg & any(is.na(E))) stop("one or more NA in E!")
  #  if (!dbg & any(is.na(FF))) stop("one or more NA in F!")
  #  if (!dbg & any(is.na(G))) stop("one or more NA in G!")
  dA <- sqrt(E*G-FF*FF) # eqals Nn, this is sqrt(det(g)) with metric tensor g
  #  area element dA = sqrt(det(g)) * du * dv
  #cat("max diff dA-Nn:",max(abs(dA-Nn)),"\n")
  #cat("min dA:",min(dA),"\n")
  # if(any(dA==0)) {cat(which(dA==0)); warning("one or more zeros in dA!") }
  (Area=int2d_s(dA,grd))
  # inverse area element = 1/(sqrt(det(g))
  ##w=which(dAsqrt(grd$dU*grd$dV))
  ##matrix(dA,grd$nv+1,grd$nu+1)
  inn<-1/dA
  # if (sum(dA==0)>0) warning(yellow("one or more dA == 0; set inverse to zero."))
  # inn[dA<bas$del]<-0 # otherwise problems with infinity/NAs and outliers in gradient H2
  # inn[is.na(inn)]<-0
  #inn[is.infinite(inn)]<-0
  inn2<-inn^2
  dV <-   1/3*(C$X[,1]*n[,1]+C$X[,2]*n[,2]+C$X[,3]*n[,3]) * dA # dA for integration
  (Volume<-int2d_s(dV,grd))
  if(dbg & any(is.na(dV))) warning("one or more NA in dV!")
  # second fundamental  coeffs. = vectorized x*y
  L<-  C$X_uu[,1]*n[,1] + C$X_uu[,2]*n[,2] + C$X_uu[,3]*n[,3]
  M<-  C$X_uv[,1]*n[,1] + C$X_uv[,2]*n[,2] + C$X_uv[,3]*n[,3]
  NN<- C$X_vv[,1]*n[,1] + C$X_vv[,2]*n[,2] + C$X_vv[,3]*n[,3]
  if ( any(is.na(L))) warning(crayon::red("one or more NA in L!"))
  if ( any(is.na(M))) warning("one or more NA in M!")
  if ( any(is.na(NN))) warning("one or more NA in N!")
  if ( any(is.infinite(L))) warning(crayon::red("one or more Inf in L!"))
  if ( any(is.infinite(M))) warning("one or more Inf in M!")
  if ( any(is.infinite(NN))) warning("one or more Inf in N!")
  curv<- -(E*NN+G*L-2*FF*M) * inn # work with positive curvature; H=(EN+GL-2FM)^2/dA^2; one dA canceled with dA=sqrt(...) du dv
  #plot(curv[!is.infinite(curv)])
  curv[is.na(curv)]<-0# warning("one or more Inf in curv!")
  (Curv<-int2d_s(curv,grd))
  curv_sq<- curv^2*inn
  curv_sq[is.na(curv_sq)]<-0# warning("one or more Inf in 2curv_sq!")
  (H2_BC <- int2d_s(curv_sq,grd))
  # dA0f = deltaA0/A set zero, since not yet tested
  H2 <- M.K_b/2 * (H2_BC - 2*M.C0*Curv) + # + K_b/2*C0^2*Area + # constant terms out
    + M.K_ADE*(Curv^2 ) /Area

  E_SCM_dens<-M.K_b/2 * (curv_sq - 2*M.C0*curv) + # + K_b/2*C0^2 + # constant terms out
    + M.K_ADE*(curv*Curv ) /Area

  retL<-list(Wb=H2,H2_BC=H2_BC,Area=Area,Volume=Volume,Curv=Curv,
             E=E,FF=FF,G=G,L=L,M=M,NN=NN,Nn=Nn,dA=dA,dV=dV,
             curv=curv,curv_sq=curv_sq,normal=normal,
             inn=inn,inn2=inn2,n=n,p.na=which(is.na(n)),
             E_SCM_dens=E_SCM_dens )

  if(plt) plot3a(C$X,grd,clip=clp)
  return(retL)
}

# for Hessian_par only use serial/R code for Cpp-code does not load in PSocket
#
# gradients wrt amplitides A, to be integrated with int2d_...() over ndof points of the grd
#   -> dA = sqrt(EG-F^2) du dv to be differentiated as well
#  Amplitudes A are two-dimensional, A[lm,k], lm for l-m-parameters 1..Ai_max to Ylm(u,v)
#       and k for spatial dimension 1...3     (lm=(l+^)^2-l+m
#

#' @export
Grad_SCM <- function(h2,grd,bas,C,int2d_m=IntegM)
{ if (!M.Rcpp) return(Grad_SCM_R(h2,grd,bas,C))
  G2<-Grad_SCM_cxx(h2,grd, bas, C, M.C0, M.Rcpp_ncores, M.K_b, M.K_ADE)
  return(list(ddA=array(G2$ddA,c(grd$ndof,bas$Ai_max,3)),
              ddV=array(G2$ddV,c(grd$ndof,bas$Ai_max,3)),
              grad_SCM=G2$grad_SCM,gradV=G2$gradV,gradA=G2$gradA,gradC=G2$gradC,
              dE=array(G2$dE,c(grd$ndof,bas$Ai_max,3)),
              dF=array(G2$dF,c(grd$ndof,bas$Ai_max,3)),
              dG=array(G2$dG,c(grd$ndof,bas$Ai_max,3))))
}

# vectorizes first dimension (i = spatial)
#' @export
Grad_SCM_R <- function(Wb, grd, bas,C,int2d_m=int2d_matrix)
{ # flag=1 # was for ddA factors from chain rule with metric dA, relevant for H2 and Curv
  Nuv<-grd$ndof
  Ai_max<-bas$Ai_max
  H2<-Wb
  #  copies for easier programming
  Xu<-C$X_u[,1]; Yu<-C$X_u[,2]; Zu<-C$X_u[,3];
  Xv<-C$X_v[,1]; Yv<-C$X_v[,2]; Zv<-C$X_v[,3];
  Xuu<-C$X_uu[,1]; Yuu<-C$X_uu[,2]; Zuu<-C$X_uu[,3];
  Xuv<-C$X_uv[,1]; Yuv<-C$X_uv[,2]; Zuv<-C$X_uv[,3];
  Xvv<-C$X_vv[,1]; Yvv<-C$X_vv[,2]; Zvv<-C$X_vv[,3];
  # diff of normals
  #  normal[,1] <- Yu *Zv  - Zu *Yv (e1)
  #  normal[,2] <- Zu *Xv  - Xu *Zv  (e2)
  #  normal[,3] <- Xu *Yv  - Yu *Xv   (e3)
  dNaz<-dNay<-dNax<-array(0.0,c(Nuv,Ai_max,3))
  for (j in 1:Ai_max) {
    #  dNax[i,,1] <- 0 # normal [,1] is indep of X[,1]-derivs ; terms like this eliminated in code below
    dNay[,j,1] <- bas$Ylm_u[,j]*Zv - bas$Ylm_v[,j]*Zu #(e1)
    dNaz[,j,1] <- bas$Ylm_v[,j]*Yu - bas$Ylm_u[,j]*Yv #(e1)

    dNax[,j,2] <- bas$Ylm_v[,j]*Zu - bas$Ylm_u[,j]*Zv #(e2)
    #  dNay[i,,2] <- 0 no y-dependence in normal[,2]
    dNaz[,j,2] <- bas$Ylm_u[,j]*Xv - bas$Ylm_v[,j]*Xu #(e2)

    dNax[,j,3] <- bas$Ylm_u[,j]*Yv - bas$Ylm_v[,j]*Yu # (e3)
    dNay[,j,3] <- bas$Ylm_v[,j]*Xu - bas$Ylm_u[,j]*Xv # (e3)
    #  dNaz[i,,3] <- 0 no z-dep in normal[,3]
  }
  ddV<-array(0.0,c(Nuv,Ai_max,3)) # dV=1/3.0*(X[,1]*n[,1]+X[,2]*n[,2]+X[,3]*n[,3])*dA
  #   for(k in 1:3) for (i in 1:Nuv) ddV[i,,k] <- (H2$normal[i,k]*Ylm[i,]) # worked in aij-simpson
  #   correct by nloptr::check.derivatives
  for(k in 1:3) for (j in 1:Ai_max) ddV[,j,k] <- (H2$normal[,k]*bas$Ylm[,j]) # normal = n*dA
  # second term ~ dNormal * X checked to be zero

  # gradients of first fundamentals
  dE<-dF<-dG<-array(0.0,c(Nuv,Ai_max,3))
  for (k in 1:3)
    for (j in 1:Ai_max){
      dE[,j,k] <- 2 * C$X_u[,k] * bas$Ylm_u[,j]
      dF[,j,k] <-     C$X_u[,k] * bas$Ylm_v[,j] + C$X_v[,k] * bas$Ylm_u[,j]
      dG[,j,k] <- 2 * C$X_v[,k] * bas$Ylm_v[,j]
    }
  ddA<-array(0.0,c(Nuv,Ai_max,3))
  for(k in 1:3)  # dA = sqrt( E G - F^2 ) == 1/inn => dA= 1/sqrt()*0.5*(dE...dF) = 0.5*inn*(dE...dF)
    for (j in 1:Ai_max)  ddA[,j,k] <- 0.5*( dE[,j,k]*H2$G + H2$E*dG[,j,k] - 2*H2$FF*dF[,j,k] ) * H2$inn

  dL<-dM<-dNN<-array(0.0,c(Nuv,Ai_max,3)) # L=Xuu*nx+Yuu*ny+Zuu*nz ; n unit normal

  nn<-H2$normal # orig: nn<-normal: like dNax is not dn but d(n*dA), normal is needed, not the unit vectors
  for (j in 1:Ai_max)  {
    dL[,j,1] <- #Xuu[i]*dNax[,j,1] # is zero
      +Yuu*dNax[,j,2] + Zuu*dNax[,j,3] + bas$Ylm_uu[,j]*nn[,1];
    dL[,j,2] <-  Xuu*dNay[,j,1] +  #Yuu[i]*dNay[,j,2] # is zero
      +Zuu*dNay[,j,3] + bas$Ylm_uu[,j]*nn[,2];
    dL[,j,3] <-  Xuu*dNaz[,j,1] + Yuu*dNaz[,j,2]+ #Zuu[i]*dNaz[,j,3] # is zero
      + bas$Ylm_uu[,j]*nn[,3];
  }
  for (j in 1:Ai_max)  {
    dM[,j,1] <- # Xuv[i]*dNax[,j,1] # is zero
      + Yuv*dNax[,j,2] + Zuv*dNax[,j,3] + bas$Ylm_uv[,j]*nn[,1];
    dM[,j,2] <-   Xuv*dNay[,j,1] + # Yuv[i]*dNay[,j,2] # is zero
      + Zuv*dNay[,j,3] + bas$Ylm_uv[,j]*nn[,2];
    dM[,j,3] <-   Xuv*dNaz[,j,1] + Yuv*dNaz[,j,2]+ # Zuv[i]*dNaz[,j,3]  # is zero
      + bas$Ylm_uv[,j]*nn[,3]; }

  for (j in 1:Ai_max)  {
    dNN[,j,1] <- # Xvv[i]*dNax[,j,1] # is zero
      +Yvv*dNax[,j,2] + Zvv*dNax[,j,3] + bas$Ylm_vv[,j]*nn[,1];
    dNN[,j,2] <- Xvv*dNay[,j,1] + # Yvv[i]*dNay[,j,2] # is zero
      +Zvv*dNay[,j,3] + bas$Ylm_vv[,j]*nn[,2];
    dNN[,j,3] <- Xvv*dNaz[,j,1] + Yvv*dNaz[,j,2]+ # Zvv[i]*dNaz[,j,3] # is zero
      + bas$Ylm_vv[,j]*nn[,3]; }

  for (j in 1:Ai_max) { # extend by /dA (for the used dN and nn=normal) and last term
    dL[,j,] <-dL[,j,]* H2$inn -  H2$L *H2$inn*ddA[,j,]
    dM[,j,] <-dM[,j,]* H2$inn -  H2$M *H2$inn*ddA[,j,]
    dNN[,j,]<-dNN[,j,]*H2$inn -  H2$NN*H2$inn*ddA[,j,]
  }

  dcurv<-dcurv_sq<-array(0.0,c(Nuv,Ai_max,3))# checked OK for dcurv

  for(k in 1:3)   # curv = - (E*NN+G*L-2*F*M) * inn
    for (j in 1:Ai_max) dcurv[,j,k] <- ( -(dE[,j,k]*H2$NN + H2$E*dNN[,j,k]  + dG[,j,k]*H2$L + H2$G*dL[,j,k]
                                           - 2*dF[,j,k]*H2$M - 2*H2$FF*dM[,j,k]) * H2$inn +
                                           - H2$curv * ddA[,j,k] * H2$inn )   # ;last from chain rule del_inn *
  # concluded that dE,dF,dG,dL,dM,dN are correct by checking dcurv

  for(k in 1:3)
    for (j in 1:Ai_max) # curv_sq = curv^2*inn ;

      dcurv_sq[,j,k] <-  ( 2*H2$curv*dcurv[,j,k] *H2$inn
                           -  H2$curv_sq*ddA[,j,k]*H2$inn )

  gradH2BC<-int2d_m(dcurv_sq,grd) # no C0, no ADE
  gradC   <-int2d_m(dcurv,grd)
  gradA   <- int2d_m(ddA,grd)
  gradH2  <-  M.K_b/2 * (gradH2BC  - 2*M.C0*gradC ) +
    + M.K_ADE * (2 * H2$Curv * gradC / H2$Area - gradA * H2$Curv^2 / H2$Area^2 )

  # pure Spontaneous Curvature Model,
  # keep deltaA_0 zero everywhere
  #
  return(list(grad_SCM=gradH2,
              grad_W_BC=gradH2BC, # int2d_matrix needs 1:ndof spatial points in leading dimension
              gradV =int2d_m(ddV,grd,bas),
              gradA =gradA,
              gradC =gradC,
              ddV=ddV,ddA=ddA,dcurv=dcurv,
              dW_BC=dcurv_sq,dE=dE,dF=dF,dG=dG)) # these e.g. for debugging
} # adding gradients of shear terms outside, by giving dE,dF,dG,ddA etc back

#
# only area and volume gradients
#
#' @export
Grad_SCM_av <- function(Wb,grd,bas,C,int2d_m=IntegM)
{ Nuv<-grd$ndof
  Ai_max<-bas$Ai_max
  H2<-Wb
  #  copies for easier programming
  Xu<-C$X_u[,1]; Yu<-C$X_u[,2]; Zu<-C$X_u[,3];
  Xv<-C$X_v[,1]; Yv<-C$X_v[,2]; Zv<-C$X_v[,3];
  Xuu<-C$X_uu[,1]; Yuu<-C$X_uu[,2]; Zuu<-C$X_uu[,3];
  Xuv<-C$X_uv[,1]; Yuv<-C$X_uv[,2]; Zuv<-C$X_uv[,3];
  Xvv<-C$X_vv[,1]; Yvv<-C$X_vv[,2]; Zvv<-C$X_vv[,3];
  dNaz<-dNay<-dNax<-array(0.0,c(Nuv,Ai_max,3))
  for (i in 1:Nuv) {
    dNay[i,,1] <- bas$Ylm_u[i,]*Zv[i] - bas$Ylm_v[i,]*Zu[i] #(e1)
    dNaz[i,,1] <- bas$Ylm_v[i,]*Yu[i] - bas$Ylm_u[i,]*Yv[i] #(e1)

    dNax[i,,2] <- bas$Ylm_v[i,]*Zu[i] - bas$Ylm_u[i,]*Zv[i] #(e2)
    dNaz[i,,2] <- bas$Ylm_u[i,]*Xv[i] - bas$Ylm_v[i,]*Xu[i] #(e2)

    dNax[i,,3] <- bas$Ylm_u[i,]*Yv[i] - bas$Ylm_v[i,]*Yu[i] # (e3)
    dNay[i,,3] <- bas$Ylm_v[i,]*Xu[i] - bas$Ylm_u[i,]*Xv[i] # (e3)
  }
  ddV<-array(0.0,c(Nuv,Ai_max,3)) # dV=1/3.0*(X[,1]*n[,1]+X[,2]*n[,2]+X[,3]*n[,3])*dA
  for(k in 1:3) for (i in 1:Nuv) ddV[i,,k] <- (H2$normal[i,k]*bas$Ylm[i,]) # normal = n*dA
  dE<-dF<-dG<-array(0.0,c(Nuv,Ai_max,3))
  for (k in 1:3)
    for (i in 1:Nuv){
      dE[i,,k] <- 2 * C$X_u[i,k] * bas$Ylm_u[i,]
      dF[i,,k] <-     C$X_u[i,k] * bas$Ylm_v[i,] + C$X_v[i,k] * bas$Ylm_u[i,]
      dG[i,,k] <- 2 * C$X_v[i,k] * bas$Ylm_v[i,]
    }
  ddA<-array(0.0,c(Nuv,Ai_max,3))
  for(k in 1:3)  # dA = sqrt( E G - F^2 ) == 1/inn => dA= 1/sqrt()*0.5*(dE...dF) = 0.5*inn*(dE...dF)
    for (i in 1:Nuv)  ddA[i,,k] <- 0.5*( dE[i,,k]*H2$G[i] + H2$E[i]*dG[i,,k] - 2*H2$FF[i]*dF[i,,k] ) * H2$inn[i]
  return(list(gradV =int2d_m(ddV,grd.bas),
              gradA =int2d_m(ddA,grd,bas),
              ddV=ddV,ddA=ddA)) # these e.g. for debugging
} # support adding gradients of shear terms outside, by giving dE,dF,dG back


#
# this gives back vertex-wise curvatures, Wb precomputed by E_SCM()
#
#' @export
Membrane_Curvatures<-function(Wb, grd, plt.K=FALSE)
{
  int2d_s=IntegS
  H2=Wb
  K=( H2$L*H2$NN - H2$M^2) * H2$inn2
  chk=int2d_s(K*H2$dA,grd)/4/pi # 1 for genus 0; need a dA in integrand
  # need to remove det(g) from H2$curv, which is full integrand
  k1=H2$curv/2*H2$inn + sqrt(H2$curv^2/4*H2$inn2 - K ) # allow complex squareroot
  k2=k1 - 2*sqrt(H2$curv^2/4*H2$inn2 - K )
  if (plt.K) plot(K,pch=19,col=densCols(K),cex=0.45)
  return(list(k1=k1,k2=k2,K=K,chk=chk,comment="curvatures: total K + principal k1,k2" ))
}


# Shear-Elastic

#
# MembraneR3 Shear-Elastic Network model components
#

#
# sum of squares of principal stretches from trace of  right Cauchy-Green-tensor
#   trace can to be differentiated analytically
#

#
# computes reference H2 and derives inv(metric_tensor) to be
#   ready for computing not individual principal stretches lambda_1,2,
#   but (lambda1^2+lambda2^2), which is needed for shear parameter \beta.
#    Area expansion parameter \alpha is computed from local area elements dA,
#    already computed in H2_etc() and reference values in returned by Ref4Cauchy... in $h2ref.
#
# a reference object is needed with precomputed metrics etc., to be stored in $h2ref
#   Lim et al. (2002) use an oblate shape of v0=0.95 for this with reference area
#    and volume c(140,100) (no other information in the paper).
#

#
# vectorization helpers, not all used
#
mat4<-function(f,g,h,j)
{ matrix(c(f,g,h,j),2,2) }
vmat4<- Vectorize(mat4,SIMPLIFY=TRUE)
#vmat4mat(rep(1,10000),rep(2,10000),rep(3,10000),rep(4,10000))->d
vmat4_2lmat<- Vectorize(mat4,SIMPLIFY=FALSE) # keeps matrices in list

matmat <- function(A,B) A%*%B # works for A and B being lists of matricess
vmatmat<-Vectorize(matmat,SIMPLIFY=TRUE)
vmatmat_2lmat<-Vectorize(matmat,SIMPLIFY=FALSE)# keeps matrices in list

Vec2Spat<-function(V,grd)
{ return(matrix(V,grd$nu,grd$nv))
}


#' @export
Ref4CauchyGreen <- function(A,grd,bas,loop=FALSE)
{ ndof=grd$ndof;Ai_max=bas$Ai_max
# compute reference fields and gradients
updateX(A,grd,bas)->C
h2ref<-E_SCM(A,grd,bas,C)
h2ref_grad<-Grad_SCM(h2ref,grd,bas,C)

gPrep<-vmat4_2lmat(h2ref$E,h2ref$FF,h2ref$FF,h2ref$G)
Dg<-sapply(gPrep,det) # det g for inversion derivative

# prepare global reference inverse metric (per spatial point [[i]])
giPrep<-lapply(gPrep,pracma::inv)

dgiPrep<-array(0,c(2,2,ndof,Ai_max,3))
cat("I dgiPrep 1..3\n")
for (k in 1:3) for (j in 1:Ai_max)  { cat(k," ",round(j/Ai_max*100,2),"\r");for (i in 1:grd$ndof)
  dgiPrep[,,i,j,k]<- # d inv G = 1/det dG~ -1/det^2 G~ * ddet G
  1/Dg[i]*matrix(c(h2ref_grad$dG[i,j,k],-h2ref_grad$dF[i,j,k],-h2ref_grad$dF[i,j,k],h2ref_grad$dE[i,j,k]),2,2) -
  1/Dg[i]^2 * matrix(c(h2ref$G[i],-h2ref$FF[i],-h2ref$FF[i],h2ref$E[i]),2,2) *
  (h2ref_grad$dE[i,j,k] * h2ref$G[i] + h2ref$E[i] * h2ref_grad$dG[i,j,k] - 2 * h2ref_grad$dF[i,j,k] * h2ref$FF[i]) #  d detG
}
dgiPrepL<-list()
l=1
cat("II dgiPrep-List 1...3 \n")
for (k in 1:3) {cat(k,"\r");for (j in 1:Ai_max) for (i in 1:ndof) {dgiPrepL[[l]]<-dgiPrep[,,i,j,k];l=l+1}}
gi<-array(unlist(giPrep),c(2,2,ndof)) # must have matrix in front!
tgi<-array(0,c(ndof,2,2))
cat("III reformat gi 2 x 2 \n")
for (k in 1:2) for (j in 1:2) {cat(i,j,"\r");tgi[,j,k]<-gi[j,k,]}
return(list( tgi=tgi, giPrep=giPrep, h2ref=h2ref, v=h2ref$Volume ,a=h2ref$Area ,c=h2ref$Curv,ARef=A ))
}

# compute stretch parameters (alpha, beta) from SEN
#' @export
SEN<-function(A,grd,bas,Ref,Wb_cur)
{
  h2cur<-Wb_cur # current E_SCM needed
  G <- vmat4_2lmat(h2cur$E, h2cur$FF, h2cur$FF, h2cur$G) # needed for gradients
  RightCauchyGreen <- vmatmat_2lmat(Ref$giPrep,G)
  m <- 0.5*sapply(RightCauchyGreen, function(x) c(x[1,1]+x[2,2])) # half trace = m
  alpha <- h2cur$dA/Ref$h2ref$dA - 1
  beta <- m/(alpha+1) - 1
#  H2_glob$alpha<<-alpha # avoid globals
#  H2_glob$beta<<-beta
  return(list(alpha=alpha, beta=beta, m=m, h2cur=h2cur,A=A, G=G))
}

# energy of Shear-Elastiv Network (SEN)
#' @export
E_SEN<-function(A,grd,bas,S,Ref)
{
  WS<- M.Ka/2 * IntegS((M.a2*S$alpha^2+M.a3*S$alpha^3+M.a4*S$alpha^4)*Ref$h2ref$dA,grd) +
    + M.mu*IntegS( ( (M.b0+M.b1*S$alpha)*S$beta + M.b2*S$beta^2)*Ref$h2ref$dA, grd)
  return(WS)
}


# needs Ref$tgi and Ref$hsref, S$m , S$alpha, S$beta
#' @export
Grad_SEN<-function(A, grd, bas, h2cur_grad, S, Ref, int2d_m = IntegM){
  #          d_beta = dm/(alpha+1) - d_alpha*m/(alpha+1)^2
  #            beta = 2m / (2(alpha+1)) - 1 -> S$beta
  #
  Ai_max=bas$Ai_max;  ndof=grd$ndof
  gradAlpha=array(0,c(ndof,Ai_max,3)) #  local grad alpha
  gradBeta=array(0,c(ndof,Ai_max,3))  #  local grad beta
  gradS=array(0,c(Ai_max,3))     # integrated grad Stretch_Energy

  dm<-array(0,c(ndof,Ai_max,3))
  for (k in 1:3)
    for (j in 1:Ai_max)
    { dm[,j,k] <- 0.5*(Ref$tgi[,1,1]*h2cur_grad$dE[,j,k]+ Ref$tgi[,1,2]*h2cur_grad$dF[,j,k] + # =  1/2 tr gi %*% dG
                         Ref$tgi[,2,1]*h2cur_grad$dF[,j,k]+ Ref$tgi[,2,2]*h2cur_grad$dG[,j,k])
    }
  gradAlpha = h2cur_grad$ddA / Ref$h2ref$dA
  gradBeta = dm / (S$alpha+1) - gradAlpha*S$m / (S$alpha+1)^2
  gradS <- (M.Ka*0.5 * int2d_m( Ref$h2ref$dA *
                              gradAlpha*(M.a2*2*S$alpha + 3*M.a3*S$alpha^2 +
                                                   + 4*M.a4*S$alpha^3),grd,bas)  +
    +   M.mu * int2d_m( Ref$h2ref$dA *
                          ( gradBeta*( M.b0 + M.b1*S$alpha + 2*M.b2*S$beta) +
                                          + M.b1*gradAlpha*S$beta), grd, bas)
  )
  return(list(grad_SEN=gradS, gradAlpha=gradAlpha,
              gradBeta=gradBeta, dm=dm ))
}


#
# utilities
#


# annotates A with L-M-strings as rownames, + some attributes
#' @export
LM2A<-function(A,bas)
{
nm=apply(bas$LM,1,paste,sep=";",collapse=";")
if (!is.matrix(A)) A=matrix(A,ncol=3)
rownames(A)<-nm;
colnames(A)<-LETTERS[24:26];
attr(A,"C0")<-M.C0
attr(A,"V0")<-bas$Target["Volume"]
attr(A,"A0")<-bas$Target["Area"]
if (bas$Nc>2) attr(A,"Ct")<-bas$Target["Curv"]

attr(A,"Target")<-bas$Target
return(A)
}

# save A as object named Alm with additional attributes
#' @export
saveAlm<-function(A,bas,file)
{
  Alm=LM2A(A,bas)
  save(Alm,file=file) # name is Alm
}

# return A from Alm-file
#' @export
loadAlm<-function(file,bas)
{
  load(file)
  return(Alm)
}


#
# Total energy from SCM and SEN
#  SCM stored in h2cur in stretches S
#

#' @export
TotalEnergyDensity<-function(S)
{
  return( S$h2cur$E_SCM_dens +
            M.Ka/2 * (S$alpha^2+M.a3*S$alpha^3+M.a4*S$alpha^4) + #*Ref$h2ref$dA +
            + M.mu* ( (1+M.b1*S$alpha)*S$beta + M.b2*S$beta^2) # *Ref$h2ref$dA
  )
}

int2d_matrix_3<-function (field_m,grd,bas=NULL)
{ dms=dim(field_m); res=array(0.0,dms[2:4]);
  for (i in 1:dms[2]) for (j in 1:dms[3]) for (k in 1:dms[4]) res[i,j,k]=int2d_scalar_GLS(field_m[,i,j,k],grd)
  return(res)
}


# plot coeffs
#' @export
plotA<-function(A,...){matplot(A,type="l",lty=1,xlab=expression((l+1)^2-l+m),...)}

#
# plot coeffs
#   and l-values as color bar on x-axis
#     only for full basis from MakeBasis...(...exclude=c(), include_m=c())
#
#' @export
plotA_l<-function(A,bas,bar=FALSE,xlab=ifelse(all(diff(bas$Lset)==1),
                                          expression((l+1)^2-l+m),"k"),
                                    scale_up=FALSE,...){
  A1=A
  if(scale_up) for (k in 1:3) A1[,k]=A[,k]*sqrt(bas$G.tk)

  matplot(A1,type="h",lty=1,xlab=xlab,...)
  if(bar)
  {
    x1=cumsum(table(bas$LM[,1]))
    x0=c(0,x1[-length(x1)])
    l=length(x0)
    segments(x0=x0,x1=x1,y0=rep(0,length(x0)),y1=rep(0,length(x0)),col=1:length(bas$Lset),lwd=3)
    text((x0+x1)/2,rep(0,length(x0)),labels=bas$Lset,cex=2,pos=1,col=1:length(bas$Lset))

  }
  matplot(A1,type="l",lty=1,add = TRUE)
}


#
# evaluate a series of coefficients in the list LA, give back quantities and coordinates in a list
#   also give two sets of contours back, fast plot by segments3d
#   todo: could be changed to give back a list of objects.
#

# helper routines for 3d graphics objects
#' @export
Obj2X<-function(O) # extract vertex coordinates
{ return(t(O$vb[1:3,]))
}

#' @export
Obj2X_centre<-function(O)
{ X=t(O$vb[1:3,]) # extract vertex coordinates and centre
for (k in 1:3 ) X[,k]=X[,k]-mean(X[,k])
return(X)
}

#' @export
X2Obj<-function(O,X)
{ O$vb=rbind(t(X),1) # ingest coordinates in 3d-graphics object
return(Rvcg::vcgClean(Rvcg::vcgUpdateNormals(O),silent=TRUE))
# clean helps smooth color at zero meridian
}

#  build vertex areas from triangle areas using Rvcg
#' @export
VertexAreasOBJ<-function(O)
{
  VF=Rvcg::vcgVFadj(O)
  FA=Rvcg::vcgArea(O,perface = TRUE)
  a=FA$pertriangle
  a_v=sapply(VF,function(x) 2*mean(a[x]))
  print(sum(a_v)-FA$area/2)
  return(a_v)
}

# synthesize X coordinates
#' @export
synthX<-function(Y,A) # Y are precomputed spherical harmonics
{ return(Y%*%A)
}


# the following imag.obj... are showing data as color code on the 3d object
#


#' @export
imag.obj.colorbar.simple<-function(obj,f,clr=TRUE,...) {
  #  if(is.matrix(f)) f<-t(f)
  if(clr) rgl::clear3d()
  cols=rainbow(100);
  # f=c(f)
  rgl::shade3d(obj,meshcolor="vertices",color=cols[(f-min(f))/diff(range(f))*99+1],...)
  rgl::bgplot3d(fields::imagePlot(legend.only = TRUE, zlim = range(f), col = cols) )
}


# allow for limits to suppress outliers (color black)
#' @export
imag.obj.colorbar<-function(obj,f,limits=range(f),clr=FALSE,pal=heat.colors,width=550,height=480,par=TRUE,...) {
  # if(is.matrix(f)) f<-t(f)
  if (min(f)==max(f)) limits=c(f-f[1]/100,f+f[1]/100)

  f[limits[1]>f  ]<-limits[1]
  f[  limits[2]<f]<-limits[2]
  cols=pal(100);
  if(par) rgl::par3d(windowRect=c(1,30,width+1,height+30));
  if (clr) rgl::clear3d();
  col=cols[(f-limits[1])/diff(limits)*99+1];
  col[limits[1]>f]="#000000"
  col[limits[2]<f]="#000000"
  rgl::shade3d(obj,meshcolor="vertices",col=col,...)
  rgl::bgplot3d(fields::imagePlot(legend.only = TRUE,add=TRUE, new=FALSE,zlim = limits, col = cols) )
}


# this has more dense points next to 0 meridian
#  alternative with regular v-spacing is GaussLegendreSimpson

#' @export
MakeGrid_GaussLegendre<-function(n=25,uv_fac=1,comment="spherical coordinates Gauss-Legendre grid, type GL",check_plt=FALSE) # assume spherical coordinates
{
  grd=list(ua=0,ub=pi,va=0,vb=2*pi,n=n,
           nu=n,nv=round(n*uv_fac,0)) # to be filled further before return
  nu=n
  nv=round(n*uv_fac,0)
  cx <- pracma::gaussLegendre(nu, 0, pi)
  x <- cx$x
  wx <- cx$w
  cy <- pracma::gaussLegendre(nv, 0, 2*pi)
  y <- cy$x
  wy <- cy$w

  grd$xg <- x
  grd$yg <- y
  mesh=pracma::meshgrid(x,y) # results are [nv,nu] !!!! important also for setting colors correctly
  grd$u=t(mesh$X) # 2D
  grd$v=t(mesh$Y) # 2D

  #dim(grd$U); range(grd$U)
  (dm=dim(grd$u))
  grd$ndof=prod(dm)

  # for vectorization we also haveu and v as vectors U,V
  grd$U=c(grd$u)
  grd$V=c(grd$v)

  range(grd$u)
  range(grd$v)

  grd$wx=wx
  grd$wy=wy
  grd$UV=cbind(grd$U,grd$V)
  # object creation
  nx=grd$nu;ny=grd$nv;
  q=matrix(NA,3,nx*ny*2);k=0
  for (i in 1:(nx-1))  for (j in 1:(ny-1)){
    k=k+1;l=(j-1)*nx+i
    q[1,k]=l
    q[2,k]=l+1
    q[3,k]=l+1+nx
    k=k+1
    q[1,k]=l
    q[2,k]=l+nx+1
    q[3,k]=l+nx
  }
  q=q[,1:k]
  x=sin(grd$u)*cos(grd$v);y=sin(grd$u)*sin(grd$v);z=cos(grd$u)
  rgl::mesh3d(x=x,y=y,z=z,triangles=q) -> M

  if(check_plt){
    # correct:
    rgl::clear3d()
    imag.obj.colorbar.simple(M,grd$v)
    rgl::contourLines3d(M,grd$v)
    rgl::title3d("looks correct for v")

    rgl::open3d()
    imag.obj.colorbar.simple(M,grd$u)
    rgl::contourLines3d(M,grd$u)
    rgl::title3d("colors in imag.obj for u")

    Rvcg::vcgBorder(M)->b
    rgl::spheres3d(x[b$bordervb],y[b$bordervb],z[b$bordervb],radius=0.15,col="black")
  }

  #  vcgUpdateNormals(M)->M
  grd$Obj<-M
  grd$comment<-comment
  grd$type="GL"

  Obj2ObjQ(grd$Obj,grd)->grd$ObjQ
  return(grd)
}


# display f data on object obj
#' @export
imag.obj<-function(obj,f) {
  cols=rainbow(100); rgl::shade3d(obj,meshcolor="vertices",col=cols[as.integer(1+99*(f-min(f))/diff(range(f)))]) }

#' imag
#'
#' plot a scalar in 2D
#' @param field scalar field, size of grd$nv x grd$nu
#' @param grd the grid of (u,v) on which the scalar is defined
#' @export
imag<-function(field,grd,nx=grd$nv,ny=grd$nu/2,...)
  {fields::quilt.plot(grd$V,grd$U,field[],nx,ny,xlab="v",ylab="u",...);invisible()}

#' MakeBasis_UV
#'
#' compute the basis functions on given (u,v) values
#' (u,v) are usually from grd$U, grd$V, but you may also use irregular (u,v)
#' for example for fitting coefficients to a 3d-object, in which (u,v)-values
#' may be given as texture coordinates texccords.
#'
#'  @param L_max (=4) spectral order of basis
#'  @param u,v vectors of spehrical angles u (in 0..pi) and v (in 0..2*pi).
#'  @return basis, to be stored as $bas in membrane object or for other use.
#'  The basis also contains standard constraints for area and volume, which may be modified with SetConstraints()
#' @export
MakeBasis_UV<-function(L_max=4,u,v)
{ # not valid for Functionals
  # give back Boschs Ylm with normalization 1/sqrt(4*pi)
  Ai_max=(L_max+1)^2-1
  LM = data.frame(l=rep(0L,Ai_max),m=rep(0L,Ai_max))
  AI=1
  for (l in 1:L_max) for (m in (-l):l) {LM[AI,]=c(l,m);AI<-AI+1}
  n.v=length(u)
  if(length(v)!=n.v) stop("u not same length like v")
  L_Ylm=L_Ylm(L_max, u,  v)
  Ylm=L_Ylm$Ylm[,-1] / sqrt(4*pi)
  Ylm_v=Ylm_v(L_max, u,  v, L_Ylm$PLK)[,-1] / sqrt(4*pi)
  Ylm_vv=Ylm_vv(L_max, u,  v, L_Ylm$PLK)[,-1] / sqrt(4*pi)
  L_Y_u=L_Ylm_u(L_max,u,v,L_Ylm$PLK)
  Ylm_u=L_Y_u$Ylm_u[,-1] / sqrt(4*pi)
  Ylm_uu=Ylm_uu(L_max,u,v,L_Y_u$P_T)[,-1] / sqrt(4*pi)
  Ylm_uv=Ylm_uv(L_max,u,v,L_Ylm$PLK,L_Y_u$P_T)[,-1] / sqrt(4*pi)
  l=LM[,1];m=LM[,2]
  bas=list(n.v=n.v,uv=cbind(u,v),
           Ylm=Ylm,Ylm_u=Ylm_u,Ylm_v=Ylm_v,
           Ylm_uu=Ylm_uu,Ylm_uv=Ylm_uv,Ylm_vv=Ylm_vv,
           LM=LM,Ai_max=Ai_max,
           l=l, m=m,
           A=matrix(0,Ai_max,3), # zero amplitudes to start with
           L_max=L_max,G.tk=l^2*(l+1)^2,
           comment="irregular or Gauss-Legendre/Simpson basis from (u,v), computed with W. Bosch/ K. Khairy codes excluding l=0, A/V constraints set",
           Nupd=0,
           Lset=unique(l),
           Mset=unique(m),
           Nc=2, Cons=c("gradA","gradV"), # default constraints; unit sphere, spontaneous curvature model
           QCons=c("Area","Volume"),# could be enhanced by Curv in SetConstraints
           Target=c(140,100), # sphere values
           TNorm =c(140,100) # ,8*pi)) # keep Norms for A, V, set C separately in SetConstraints
  )

  names(bas$Cons)=names(bas$QCons)=names(bas$TNorm)=names(bas$Target)=c("Area","Volume")
  bas$A=LM2A(bas$A,bas)
  return(bas)
} # G.tk is diag Tikhonov, Gamma^T Gamma

#' updateX
#'
#' updates coordinates C$X and their partial derivatives wrt. u,v, like C$Xu.
#' @param A coefficients of shape
#' @param grd grid from on which the basis is computed
#' @param bas basis function values $Ylm and their derivatives, like $Ylm_u.
#' @return basis object, type is list.
#' @export
updateX<-function(A,grd,bas )
{
  bas$Ylm %*% A -> X
  bas$Ylm_u %*% A -> X_u
  bas$Ylm_v %*% A -> X_v
  bas$Ylm_uu %*% A -> X_uu
  bas$Ylm_uv %*% A -> X_uv
  bas$Ylm_vv %*% A -> X_vv
  # returned value usually stored as "C" for "Coordiinates"
  return(list(X=X,X_u=X_u,X_v=X_v,X_uu=X_uu,X_vv=X_vv,X_uv=X_uv,Amp=A))
}

# really needed? maybe for plotlseies?
updateX_subset<-function(A,grd,bas,Ai_max=dim(A)[1])
{
  s=1:Ai_max
  bas$Ylm[,s] %*% A[s,] -> X
  bas$Ylm_u[,s] %*% A[s,] -> X_u
  bas$Ylm_v[,s] %*% A[s,] -> X_v
  bas$Ylm_uu[,s] %*% A[s,] -> X_uu
  bas$Ylm_uv[,s] %*% A[s,] -> X_uv
  bas$Ylm_vv[,s] %*% A[s,]-> X_vv
  return(list(X=X,X_u=X_u,X_v=X_v,X_uu=X_uu,X_vv=X_vv,X_uv=X_uv,Amp=A))
}

# compute coordinates and derivatives
# could be used instead of updateX in several places

#' updateX_only
#'
#' updates coordinates C$X, NOT the derivatives.
#' @param A coefficients of shape
#' @param grd grid from on which the basis is computed
#' @param bas basis function values $Ylm
#' @return basis object, type is list.
#' @export
updateX_only<-function(A,grd,bas)
{
  bas$Ylm %*% A -> X
  return(list(X=X,Amp=A))
}

#
#  make A coeffs for a unit sphere
#

#' MakeSphere
#'
#' compute unit sphere for given grid and basis
#' @param grd given grid
#' @param bas given basis
#' @param r (=1) for radius of output sphere
#' @return coefficient matrix, derived from its prototype bas$A
#' @export
MakeSphere<-function(grd,bas,r=1)
{
  A=bas$A; A[,]=0
  A[1,"X"]=1/0.48860251190292*r
  A[2,"Z"]=1/0.48860251190291*r
  A[3,"Y"]=1/0.48860251190292*r
  return(A)
}

# the central driver to 2D integration
# use IntegS from cpp instead
int2d_scalar_GLS<-function(F2,grd)
{ Z <- matrix(F2,grd$nu,grd$nv)
  Q <- grd$wx %*% Z %*% as.matrix(grd$wy)
return(Q[,])
}

# plot X/object from segments in 3d
# with Cartesian Wireframe as (X,Y,Z)-contour levels for cont=TRUE
#' plot3a
#'
#' plot coordinates C$X (after updateX)
#' @param X Coordiates, e.g. in C$X from updateX()
#' @param grd grid with a basic rgl-object grd$Obj
#' @examples
#' data("M1")
#'  # take required data from M1
#' plot3a(updateX(M1$A,M1$grd,M1$bas,M1$)$X,M1$grd)
#'
#' @export
plot3a<-function(X,grd,pnts=FALSE,clip=FALSE,col="black",alpha=1,cont=TRUE,cont.grid=FALSE,fill=TRUE,fn="z",fn_data="z",...)
{
  X2Obj(grd$Obj,X)->O
  if (cont.grid)   {id<-rgl::contourLines3d(O,grd$U,col=col,levels=pracma::linspace(0,pi,12));
  rgl::contourLines3d(O,grd$V,col=col,levels=pracma::linspace(0,2*pi,13))}
  if (cont)   {id<-rgl::contourLines3d(O,X[,1],col=1);rgl::contourLines3d(O,X[,2],col=2);rgl::contourLines3d(O,X[,3],col=3)}
  if(fill) {#Rvcg::vcgUpdateNormals(O)->O;
    if(fn=="z")id<-rgl::filledContour3d(O,fn=X[,3],alpha=alpha)
  else id<-rgl::filledContour3d(O,fn=fn_data,alpha=alpha,...)}
  if(pnts) id<-rgl::points3d(X[,1],X[,2],X[,3],col="red",cex=2)
  if(clip) rgl::clipplanes3d(c(0.5,0.5,0))
  invisible()
}

#' plot3q
#'
#' plot coordinates C$X (after updateX) with wireframe of quadrilaterals
#' @param X Coordiates, e.g. in C$X from updateX()
#' @param grd grid with a basic rgl-object grd$Obj
#' @examples
#' data("M1")
#'  # take required data from M1
#' plot3q(updateX(M1$A,M1$grd,M1$bas,M1$)$X,M1$grd)
#' @export
plot3q<-function(X,grd,col="black",alpha=1,...)
{
  plot3b(X,grd,wire=FALSE,...)
  if (is.null(grd$ObjQ)) Obj2ObjQ(grd$Obj,grd)->Q else Q=grd$ObjQ
  X2ObjQ(Q,X)->Q
  rgl::wire3d(Q,col="black",specular="black")
}
#' plot3qs
#'
#' plot with coordinates C$X (after updateX) and scalar as color code
#' @param X Coordiates, e.g. in C$X from updateX()
#' @param grd grid with a basic rgl-object grd$Obj
#' @param s scalar to plot as color code on shape
#' @examples
#' data("M1"); SetParams(M1)
#'  # take required data from M1
#'  update(M1,"dA")->M1
#'  #plot area sizes as color code
#' plot3qs(updateX(M1$A,M1$grd,M1$bas,M1$)$X,M1$grd,M1$dA)
#'
#' @export
plot3qs<-function(X,grd,s,alpha=1,specular="black",...)
{ X2Obj(grd$Obj,X)->O;
  Rvcg::vcgUpdateNormals(O)->O
  col=heat.colors(100)[1+99*(s-min(s))/diff(range(s))]
  rgl::shade3d(O,col=col,specular=specular,...)
  if (is.null(grd$ObjQ)) Obj2ObjQ(grd$Obj,grd)->Q else Q=grd$ObjQ
  X2ObjQ(Q,X)->Q
    rgl::wire3d(Q,col="gray",specular="black")
}


#' @export
plot3b<-function(X,grd,col="white",specular="black",wire=TRUE,...)
{X2Obj(grd$Obj,X)->O;
  Rvcg::vcgUpdateNormals(O)->O
  rgl::shade3d(O,col=col,specular=specular,...)
if(wire) rgl::wire3d(O,col="black",lwd=2,specular="black")
}


#
# this should be the standard!
#
#' @export
MakeGrid_GaussLegendreSimpson<-function(n=20,ua=0,ub=pi,va=0,vb=2*pi,del_Ylm=1e-6,comment="spherical coordinates Gauss-Legendre-Simpson grid, type GLS",check_plt=FALSE) # assume spherical coordinates
{
  if (n%%2==1) n=n+1
  grd=list(ua=ua,ub=ub,va=va,vb=vb) # to be filled further and returned
  nu=n; nv=n+1 # we double the first v data point (v=0) at v=2pi
  cx <- pracma::gaussLegendre(nu, ua, ub)
  x <- cx$x
  wx <- cx$w
  y <- pracma::linspace(va,vb,n=nv)# keep last [-(nv+1)] # last is too much; equals first (periodics only)
  wy <- rep(0,nv)
  n1=n/2
  wy[1]=wy[nv]=1 #2 # dont have last point extra, therefore double first points weight (periodics only)
  wy[(1:n1)*2]=4 # "half x"
  wy[(1:(n1-1))*2+1]=2 # "full x"
  # wy
  length(wy)
  h=(vb-va)/n1
  wy=wy*h/6
  sum(wy)-(vb-va) # must be zero (integrand function == 1)
  grd$xg <- x
  grd$yg <- y
  mesh=pracma::meshgrid(x, y)
  grd$u=t(mesh$X) # 2D # for [nu,nv] adressing
  grd$v=t(mesh$Y) # 2D
  (dm=dim(grd$u))
  grd$ndof=prod(dm)
  #  vectorization we keep  u and v as vectors U,V
  grd$U=as.vector(grd$u)
  grd$V=as.vector(grd$v)
  #  range(grd$u)
  #  range(grd$v) # should have last < 2 pi
  grd$wx=wx;  grd$wy=wy
  grd$UV=cbind(grd$U,grd$V)
  # object creation not yet ready
  nx=nu;ny=nv;
  q=matrix(NA,3,nx*ny*2);k=0 # 2 trinagles per quad (i,j)
  for (i in 1:(nx-1))  for (j in 1:(ny-1)){
    k=k+1;l=(j-1)*nx+i
    q[,k]=c(l,l+1,l+1+nx)
    k=k+1
    q[ ,k]=c(l,l+nx+1,l+nx)
  }
  q=q[,1:(k)]
  print(k)
  x=sin(grd$u)*cos(grd$v);y=sin(grd$u)*sin(grd$v);z=cos(grd$u)
  rgl::mesh3d(x,y,z,triangles=q, normals = list(x=x,y=y,z=z) ) -> M
 # str(M)
  #  clear3d();
  # vcgUpdateNormals(M)->M # creates an artificial contrast at 0/2pi v-boundary
  grd$Obj<-M
  grd$comment<-comment
  grd$type="GLS"
  grd$n=n;grd$nu=nu;grd$nv=nv
  #  M=grd$Obj
  if(check_plt){
    rgl::clear3d()
    imag.obj.colorbar.simple(M,grd$v)
    rgl::contourLines3d(M,grd$v)
    rgl::title3d("looks correct for v")
    rgl::open3d()
    imag.obj.colorbar.simple(M,grd$u)
    rgl::contourLines3d(M,grd$u)
    rgl::title3d("colors in imag.obj for u")
  }
  Obj2ObjQ(grd$Obj,grd)->grd$ObjQ
  return(grd)
}


# return last n elements from v; v may be list
#' @export
last<-function(v,n=1)
{ n=min(n,length(v))
  r=v[(length(v)-n+1):length(v)]
  if (is.list(r) & n==1) r=r[[1]] # if you want last element, you'd have to add [[1]]
  return(r)}

#
# this is for non-scalar fields, e.g. gradients,
#   where gradients in all grid points 1...Nuv with respect to upto (max(l)+1)^2 coefficients are evaluated
#

# you may use IntegM from c++ instead
int2d_matrix <- function (field_m,grd,bas=NULL)
{ #str(field_m)
  dms=dim(field_m);
#  print(dms)
  res=array(NA,dms[2:3]);
for (i in 1:dms[2]) for (j in 1:dms[3]) res[i,j]=int2d_scalar_GLS(field_m[,i,j],grd)
return(res)
}

# not  used
#' @export
int2d_matrix_cxx <- function (field_m,grd,bas) # needed for C++-code of gradients
{ F2<-array(field_m, c( grd$ndof,bas$Ai_max,3))
  dms=dim(F2);res=array(NA,dms[2:3]);
 for (i in 1:dms[2]) for (j in 1:dms[3]) res[i,j]=int2d_scalar_GLS(F2[,i,j],grd)
 return(res)
}

#' @export
vectomat_cxx<-function(x,Aimax)  return(matrix(x,Aimax,3))

#' @export
vectoarr_cxx<-function(x,ndof,Aimax)  return(array(x,c(ndof,Aimax,3)))


#
# create a basis on (u,v)-data, also useful for irregular arranged (u,v)
# if you dont want l=0, put 0 in exclude=c(0,...)
# works also for Gauss-Legendre-Simpson mixed grid
#    checked smaller del; produces artifacts and zeros problems in dA
#


# old finite difference derivatives of Ylm
# always exclude constant terms, i.e. shifts in mean(X)
if (HAVE_DEPRECATED){
MakeBasis_GaussLegendre_old<-function(L_max=4,u,v,exclude=c(0),include_m=c(),del=1e-5,altern=TRUE,old=FALSE)
{ source("sphericalharmonics.R")
  .sh=sphericalharmonics
  # not valid for Functionals
  inclm=(length(include_m)>0)
  Ai_max<-0; for (l in 0:L_max) if(!l %in% exclude) for (m in (-l):l) if(inclm &&  m%in% include_m) Ai_max<-Ai_max+1 else {if (!inclm) Ai_max<-Ai_max+1}
  LM = data.frame(l=rep(0L,Ai_max),m=rep(0L,Ai_max))

  AI<-1; for (l in 0:L_max) if(!l %in% exclude) for (m in (-l):l) if(inclm &&  m%in% include_m) {LM[AI,]=c(l,m);AI<-AI+1} else {if (!inclm) {LM[AI,]=c(l,m);AI<-AI+1}}
  n.v=length(u)
  if(length(v)!=n.v) stop("u not same length like v")
  uv=cbind(u,v)
  Ylm<-Ylm_u<-Ylm_v<-Ylm_uu<-Ylm_uv<-Ylm_vv<-array(0.0,c(n.v,Ai_max)) # store all Ylm_u derivative values as central differences

  for (i in 1:Ai_max) {    l=LM[i,1];m=LM[i,2]
  Ylm[,i]=.sh(l,m,uv) }

  uminus=u-del;  uplus=u+del
  vminus=v-del;  vplus=v+del
  uminus2=u-2*del;  uplus2=u+2*del
  vminus2=v-2*del;  vplus2=v+2*del
  if (old)   for (i in 1:Ai_max) {
    l=LM[i,1];m=LM[i,2]
    Ylm_u[,i]<-(.sh(l,m,cbind(uplus,v),"real")-.sh(l,m,cbind(uminus,v),"real"))/2/del
    Ylm_v[,i]<-(.sh(l,m,cbind(u,vplus),"real")-.sh(l,m,cbind(u,vminus),"real"))/2/del
    Ylm_uu[,i]<-(.sh(l,m,cbind(uplus,v),"real")-2*.sh(l,m,uv,"real")+.sh(l,m,cbind(uminus,v),"real"))/del^2
    Ylm_vv[,i]<-(.sh(l,m,cbind(u,vplus),"real")-2*.sh(l,m,uv,"real")+.sh(l,m,cbind(u,vminus),"real"))/del^2
    Ylm_uv[,i]<-(.sh(l,m,cbind(uplus,vplus),"real")-.sh(l,m,cbind(uplus,vminus),"real") -
                   .sh(l,m,cbind(uminus,vplus),"real")+.sh(l,m,cbind(uminus,vminus),"real"))/4/del^2
  } else
    for (i in 1:Ai_max) {
      l=LM[i,1];m=LM[i,2]

      if(altern) Ylm_u[,i]<-(.sh(l,m,cbind(uminus2,v))-
                               8*.sh(l,m,cbind(uminus,v))+
                               8*.sh(l,m,cbind(uplus,v))-
                               .sh(l,m,cbind(uplus2,v)))/12/del else
                                 Ylm_u[,i]<-(.sh(l,m,cbind(uplus,v))
                                             -.sh(l,m,cbind(uminus,v)))/2/del

                               if(altern) Ylm_v[,i] <-  -m* .sh(l,-m,uv) else
                                 Ylm_v[,i]<-(.sh(l,m,cbind(u,vplus))-
                                               .sh(l,m,cbind(u,vminus)))/2/del
                               if(altern) Ylm_uu[,i]<-  ( -.sh(l,m,cbind(uminus2,v))
                                                          +16*.sh(l,m,cbind(uminus,v)) +
                                                            -30*.sh(l,m,uv)
                                                          +16*.sh(l,m,cbind(uplus,v)) +
                                                            -.sh(l,m,cbind(uplus2,v)) )/12/del^2  else
                                                              Ylm_uu[,i]<-(.sh(l,m,cbind(uplus,v))-2*.sh(l,m,uv)+.sh(l,m,cbind(uminus,v)))/del^2
                                                            if(altern) Ylm_vv[,i]<- - m*m*.sh(l,m,uv) else
                                                              Ylm_vv[,i]<-( -.sh(l,m,cbind(u,vminus2))
                                                                            +16*.sh(l,m,cbind(u,vminus)) +
                                                                              -30*.sh(l,m,uv)
                                                                            +16*.sh(l,m,cbind(u,vplus)) +
                                                                              -.sh(l,m,cbind(u,vplus2)) )/12/del^2

                                                            if(altern) Ylm_uv[,i]<- -m*(+.sh(l,-m,cbind(uminus2,v)) +
                                                                                          -8*.sh(l,-m,cbind(uminus ,v))+
                                                                                          +8*.sh(l,-m,cbind(uplus,v))  +
                                                                                          -.sh(l,-m,cbind(uplus2,v)) )/12/del else
                                                                                            Ylm_uv[,i]<-(  .sh(l,m,cbind(uplus,vplus)) +
                                                                                                             .sh(l,m,cbind(uminus,vminus)) -
                                                                                                             .sh(l,m,cbind(uplus,vminus)) -
                                                                                                             .sh(l,m,cbind(uminus,vplus)) )/4/del^2
    }
  # basis return:
  l=LM[,1];m=LM[,2]
  bas=list(n.v=n.v,uv=uv,
           Ylm=Ylm,Ylm_u=Ylm_u,Ylm_v=Ylm_v,
           Ylm_uu=Ylm_uu,Ylm_uv=Ylm_uv,Ylm_vv=Ylm_vv,
           LM=LM,Ai_max=Ai_max,del=del,
           l=l, m=m,
           A=matrix(0,Ai_max,3), # zero amplitudes to start with
           L_max=L_max,G.tk=l^2*(l+1)^2,
           comment=paste(c("irregular or Gauss-Legendre/Simpson basis from (u,v), excluding l %in% ",exclude,"and including m %in% ",include_m),collapse=" " ),uv=cbind(uv),Nupd=0,
           Lset=unique(l),
           Mset=unique(m),
           Nc=2, Cons=c("gradA","gradV"), # default constraints; unit sphere, spontaneous curvature model
           QCons=c("Area","Volume"), # could be enhanced by Curv in SetConstraints
           Target=c(4*pi,4/3*pi),    # sphere values
           TNorm =c(4*pi,4/3*pi)     # ,8*pi)) # keep Norms for A, V, set C separately in SetConstraints
  )
  names(bas$Cons)=names(bas$QCons)=names(bas$TNorm)=names(bas$Target)=c("Area","Volume")
  bas$A=LM2A(bas$A,bas)
  return(bas)
} # G.tk is diag Tikhonov, Gamma^T Gamma
}



#
# if you have an initial form in X2fit (been computed with u,v according to basis for Ylm!!)
#   Fit without regularization
#
#' @export
FitAlm <- function(X2fit,bas)
{ A<-array(0,c(bas$Ai_max,3)) # amplitudes for x,y,z
for (k in 1:3) { m<-lm(X2fit[,k]~bas$Ylm[,]);m$coefficients[-1]->A[,k]}
return(A)
}

# fit with regularization:
# filtering high frequencies in least squares
#   lambda should be tested systematically by L curve discussion
#' @export
FitAlm_Tikhonov<-function(X,bas,lambda=0.0385) # this lambda is random
{
  A<-array(0,c(dim(bas$Ylm)[2],3))
  B=t(bas$Ylm)%*%bas$Ylm
  # l=bas$LM[,1]
  InvB1=pracma::inv(B + lambda * diag(bas$G.tk))
  for (k in 1:3) A[,k]=InvB1%*%t(bas$Ylm)%*%X[,k]
  return(A)
}


# give angles from spherical 3D coordinates
#' @export
inv_sph<-function(X)
{ r=sqrt(sum(X^2)) # see https://mathworld.wolfram.com/SphericalCoordinates.html
return(c(acos(X[3]/r),atan2(X[2],X[1]))) # atan2 takes care of octants
}

# give angles u,v from a starlike 3d-obnoject, to be centred
#' @export
radial_uv<-function(starlike_obj)
{
  for (k in 1:3) starlike_obj$vb[k,]=starlike_obj$vb[k,]-mean(starlike_obj$vb[k,])
 # str(starlike_obj)
  rgl::wire3d(starlike_obj)
  uv=apply(starlike_obj$vb[1:3,],2,inv_sph)
  uv[,2]=uv[,2]+pi
  print(dim(uv))
  plot(t(uv),pch=".")
  return(t(uv))
}

#
# synthesize, but only for one spatial component in coefficients  A==A[,k] !!!
#   for dim in 1..3: see synthX
#synth<-function(Y,A) { res<-A[1]*Y[,1];for (i in 2:length(A)) res <- res + A[i] * Y[,i];return(res)}
#' @export
synth<-function(Y,A) { return(Y%*%A)}

#synth_s<-function(Y,A,mx) { res<-A[1]*Y[,1];for (i in 2:mx) res <- res + A[i] * Y[,i];return(res)}
synth_s<-function(Y,A,mx) { return(Y[,1:mx]%*%A[1:mx])}


#' @export
deltaX_norm<-function(C1,C2)
{ return(apply(abs(C2$X-C1$X),2,pracma::Norm))
}


# save a set of coefficients A to file
#' @export
saveA<-function(A,file)
{
  save(A,file=file)
}

# load a set of coeffs, make it compatible for the current basis bas.
#' @export
loadA<-function(file,bas) # loads amplitudes and stores according to basis bas
{   # if bas is larger, empty amplitudes are kept zero
  load(file) # was saved from A
  A2<-A
  A<-array(0,c(bas$Ai_max,3));
  L=min(dim(A)[1],dim(A2)[1])
  for (k in 1:3) A[1:L,k]<-A2[1:L,k] # copy only first coefficients
  return(A)
}

#' @export
loadAlm<-function(file,bas) # loads amplitudes and stores according to basis bas
{   # if bas is larger, empty amplitudes are kept zero
  load(file) # was saved from A
  A2<-Alm
  A<-array(0,c(bas$Ai_max,3));
  L=min(dim(A)[1],dim(A2)[1])
  for (k in 1:3) A[1:L,k]<-A2[1:L,k] # copy only first coefficients
  return(A)
}

# rotation first around x by px, then y then z-axis
#' @export
rotateX<-function(X,px=pi,py=0,pz=0)
{
  R.mz=matrix(c(cos(pz),-sin(pz),0,sin(pz),cos(pz),0,0,0,1),3,3)
  R.mx=matrix(c(1,0,0, 0,cos(px),-sin(px), 0,sin(px),cos(px)),3,3)
  R.my=matrix(c(cos(py),0,-sin(py), 0,1,0, sin(py),0,cos(py)),3,3)
  R<-R.mz%*%R.my%*%R.mx
  return(R%*%X)
}

rotateXbyM<-function(X,M)
{
  return(M%*%X)
}
#
# rotate coefficients by rotating coordinates
#   gives back a rotation error (new from fit vs. rotated coord)
#  rotation order is by X,Y,Z-axis
#
#' @export
rotateA <- function(A,bas,grd,px=pi/2,py=-pi/2,pz=pi/2,plt=FALSE)
{
  updateX(A,grd,bas)->C
  if(plt) id1<-plot3b(C$X,grd)
  R.mz=matrix(c(cos(pz),-sin(pz),0,sin(pz),cos(pz),0,0,0,1),3,3)
  R.mx=matrix(c(1,0,0, 0,cos(px),-sin(px), 0,sin(px),cos(px)),3,3)
  R.my=matrix(c(cos(py),0,-sin(py), 0,1,0, sin(py),0,cos(py)),3,3)
  R<-R.mz%*%R.my%*%R.mx
  X1=C$X;for (i in 1:dim(X1)[1]) X1[i,]=R%*%C$X[i,]
  A1<-A
  # fit new amplitudes A1 from rotated coords X1 against basis for each sp. dim.
  for (k in 1:3)
    lm(X1[,k]~bas$Ylm)$coefficients[-1]->A1[,k]
  updateX(A1,grd,bas)->C1
  if(plt){open3d();id2<-plot3b(C1$X,grd ) }
  return(list(A=A1,C=C1,rot_err=pracma::Norm(C1$X-X1)))
}

#
# perturb coeffs a bit, but damped by Tikhonov diagonal
#


# A must not be matrix
#' @export
pertA_Gauss<-function(A,bas,sd,flt=FALSE){
 N=bas$Ai_max;
 sd1=sd/sqrt(bas$G.tk);
 dA<-rnorm(3*N,sd=rep(sd1,3))
 dA[c(1+N,1+2*N)]=0
 dA[c(2,2+N)]=0
 dA[c(3,3+2*N)]=0
 A[]=A[]+dA
 return(A)
}

#' @export
pertA_Unif<-function(A,bas,sd,flt=FALSE){
  N=bas$Ai_max;
  sd1=sd/sqrt(bas$G.tk);
dA=(runif(3*N)-0.5)/2*rep(sd1,3) # > version 14: changed to sd1

# >version 14: filter outside, e.g. ApplyFilter_L1(A,bas), with bas$flt set
#dA[c(1+N,1+2*N)]=0
#dA[c(2,2+N)]=0
#dA[c(3,3+2*N)]=0

A[]=A[]+dA
return(A)
}

#
# shows a series of up to nr x nc images and "rep" values for the l>0 present in bas
#
#' @export
plotLseries<-function(nr=4,nc=5,A,C,grd,bas,reduced=FALSE,rec=FALSE,fill=TRUE, rep="H2",stretch=TRUE,S=S,Ref=Ref )
{
 # rgl::open3d()
  rgl::mfrow3d(nr,nc,sharedMouse = TRUE)
  if(rec) M=matrix(0,length(unique(bas$LM[,1])),7) # to return summed up quantities and Norm of coeffs
  k=1
  nrm=c(1,1,1,1,1,1,sum(abs(A[bas$LM[,1]==1,])))
  names(nrm)=c("Area","Volume","Curv","H2","H2_BC","ES","Norm")
  for (l in bas$Lset) # works with "exclude" in basis, but 0 may not be excluded
  {w=which(bas$LM[,1]==l)
  last=w[length(w)] #; print(last)
  Y=synthX(bas$Ylm[,1:last],A[1:last,])
  #plot3d(0,0,0,deco=FALSE,col="white")
  A1=A;A1[]=0;A1[1:last,]=A[1:last,]
  updateX(A1,grd,bas)->C2
  h2<-E_SCM(A1,grd,bas,C2)
  if (stretch) E=E_SEN(A,grd,bas,S,Ref)
  if(rec){
   M[k,1]=h2$Area
   M[k,2]=h2$Volume
   M[k,3]=h2$Curv
   M[k,4]=h2$Wb
   M[k,5]=h2$H2_BC
   if(stretch) M[k,6]=E
   M[k,7]=sum(abs(A[bas$LM[,1]==l,]))
  }
  plot3qs(Y,grd,S$alpha)
  rgl::title3d(paste(l,": E",round((h2$Wb+E)/M.Es,3)," C",round(h2$Curv,3)))
  k=k+1
  if(k>nr*nc) break() # exit loop if screen is full
  if(l<last(bas$Lset)) rgl::next3d()
  }
  if(rec){
    M1=apply(M,2,diff)
    M2=apply(M1,2,function(x) (x)/diff(range(x)))
    matplot(M2,type="l",lty=1,lwd=2,ylab=expression(Delta*Q),pch=20)
    matpoints(M2,pch=20)
    legend("topright",col=1:6,lty=1,lwd=2,leg=c("A","V","C","H2","H2_BC","Nrm"),pch=20)
    rownames(M)=sort(unique(bas$LM[,1]))[-1]
  if (reduced) for (i in 1:5) M[,i]=M[,i]/nrm[i]
  if(!reduced) colnames(M)=c("A","V","C","H2","H2_BC","Norm") else colnames(M)=c("a","v","c","h2","h2_BC","norm")
  M=as.data.frame(M)
  return(M)}  else return("Plot L Series done")
}

#
# fill e.g. grd$Obj  with new coords
#  and return 3D object, plot if wanted
#
#' @export
MakeOBJ <- function(A,grd,bas,C,col_wire="white",col="black",plot=FALSE)
{
  O=grd$Obj;
  O=Rvcg::vcgUpdateNormals(X2Obj(O,C$X))
  if(plot) {  rgl::clear3d();
    rgl::shade3d(O, col=col); rgl::wire3d(O,col=col_wire)
    rgl::title3d("made object")}
  return(O)
}

#
#  Laplacian matrix from graph
#
#' @export
Membrane_LaplacianOBJ <- function(X)
{
  # this procedure doubles all entries in L
  ia=c(X$it[1,],X$it[2,],X$it[3,]);ja=c(X$it[2,],X$it[3,],X$it[1,])
  N<-max(max(ia),ja); TN <- Matrix::sparseMatrix(dims=c(N,N),i=ia,j=ja,x=rep(1,length(ia)),
                                         use.last.ij=FALSE)
  ig <- igraph::graph.adjacency(TN,mode="undirected",diag = FALSE)
  # remove doubles by /2
  return(igraph::laplacian_matrix(ig,normalized = FALSE,sparse = TRUE)/2)
}

# from 3D-object X compute mesh Laplacian L, normalized Ln, Diagonal D
#' @export
Membrane_LaplaciansOBJ <- function(X)
{
  ia=c(X$it[1,],X$it[2,],X$it[3,]);ja=c(X$it[2,],X$it[3,],X$it[1,])
  N<-max(max(ia),ja); TN <- sparseMatrix(dims=c(N,N),i=ia,j=ja,x=rep(1,length(ia)),
                                         use.last.ij=FALSE)
  ig <- igraph::graph.adjacency(TN,mode="undirected",diag = FALSE)
  # recompute  u directly from Laplacian
  L=igraph::laplacian_matrix(ig,normalized = FALSE,sparse = TRUE)/2
  D=L;diag(D)<-0
  return(list(L=L,Ln=igraph::laplacian_matrix(ig,normalized = TRUE,sparse = TRUE),D=D))
}

#' @export
Membrane_Laplacian_cotan <- function(x,M) # input object x, M input matrix to fill
{ # todo: better directly go through matrix M in sparse format,  then point back to mesh for cotan
  (N.v=dim(x$vb)[2])
  (N.t=dim(x$it)[2])
  M[,]=0
  I=matrix(c(1,2,3,3,1,2,2,3,1),3,3)
  for (tr in 1:N.t) {
    # go through 3 triangles sides
    #  I<-1:3
    for (i1 in 1:3)
    {
      i<-x$it[I[i1,1],tr]
      j<-x$it[I[i1,2],tr]
      k<-x$it[I[i1,3],tr]
      u1<-x$vb[1:3,i]-x$vb[1:3,k]
      v1<-x$vb[1:3,j]-x$vb[1:3,k]
      cot<-sum(u1*v1)/sqrt(sum(cross(u1,v1)^2))
      M[i,j] <- M[i,j] + cot;
      M[j,i] <- M[j,i] + cot;
      M[i,i] <- M[i,i] - cot;
      M[j,j] <- M[j,j] - cot
    }
    if (tr%% 100 ==0 )cat(tr/N.t,"  progress cotan-Laplacian  \r")
  };
  return(M)
}



#
# get Eigensystem for mesh Laplacian (created from an rgl object -> Membrane_LaplacianOBJ())
#
#' @export
Membrane_Eig<-function(L,which=1:4,kind="SM")
{
  k=max(which)
  ee=RSpectra::eigs(L,which = kind,k=k,opts = list() )
  print("order of Eigenvalues/ Vectors:")
  print(ee$values)
  R=ee$vectors[,which]
  attr(R,"values")<-ee$values
  return(R)
}

#
# plot arrows of change between two shape states
#
#' @export
imag.delta.arrowws.pca <- function(A1,A2,grd2,bas2,O1)
{
  if(!is.matrix(A1)) A1=matrix(A1,bas2$Ai_max,3);
  print(dim(A1))
  if(!is.matrix(A2)) A2=matrix(A2,bas2$Ai_max,3);
  C1<-updateX(A1,grd2,bas2)
  C2<-updateX(A2,grd2,bas2)
  C1$X=princomp(C1$X)$scores
  C2$X=princomp(C2$X)$scores

  S=matrix(0,dim(C1$X)[1]*2,3)
  for (i in (1:dim(C1$X)[1])) {S[i*2-1,]=C1$X[i,]; S[i*2,]=C2$X[i,] }
  rgl::clear3d();rgl::segments3d(S)
  rgl::shade3d(O1,alpha=0.56,col="white")
}


#
# svd alignment, see
# https://www.cse.wustl.edu/~taoju/cse554/lectures/lect07_Alignment.pdf
#

#' @export
imag.delta.aligned <- function(A1,A2,grd2,bas2,shade2=FALSE)
{

  if(!is.matrix(A1)) A1=matrix(A1,bas2$Ai_max,3);
  #print(dim(A1))
  if(!is.matrix(A2)) A2=matrix(A2,bas2$Ai_max,3);
  C1<-updateX(A1,grd2,bas2)
  centr1=apply(C1$X,2,mean)
  C2<-updateX(A2,grd2,bas2)
  centr2=apply(C2$X,2,mean)
  X1=t(apply(C1$X,1,function(x) x-centr1))
  X2=t(apply(C2$X,1,function(x) x-centr2))

  h=t(X1)%*%X2
  s=svd(h)
  d=sign(det(s$v %*% t(s$u)))
  r=diag(c(1,1,d))
  R=s$v %*% r %*% t(s$u)

  X1p=t(apply(X1,1,function(x) R%*%x))

  C2$X=X2
  C1$X=X1p
  MakeOBJ(A1,grd2,bas2,C1)->O1
  MakeOBJ(A2,grd2,bas2,C2)->O2

  S=matrix(0,dim(X1)[1]*2,3)
  for (i in (1:dim(X1)[1])) {S[i*2-1,]=X1p[i,]; S[i*2,]=X2[i,] }
  rgl::clear3d();rgl::segments3d(S)

  rgl::points3d(X1p,col=1)
  rgl::points3d(X2,col=2)
  rgl::shade3d(O1,alpha=0.6,col="white")
  d=apply((X2-X1p),1,pracma::Norm)
  if (shade2) {rgl::open3d();rgl::shade3d(O1,alpha=0.5,col="red");  rgl::shade3d(O2,alpha=0.5,col="white")}

  return(list(obj1=O1,obj2=O2,C1=C1,C2=C2,r=r,seg3d=S,disp=d))

}

#
# read data from surface evolver dump file
#  quite slow for renumbering
#
#' @export
SE_2_OBJ <- function(f.in,f.out,comment="o")
{

  r<-readLines(f.in)
  r
  f=grep("^vertices",r)
  f1=r[f[1]]

  (nv=as.numeric(substr( f1,20,nchar(f1))))
  #read.table(f.in,header=FALSE,skip=f[2],nrows=nv-1,comment=comment)[1:10,]

  X=read.table(f.in,header=FALSE,skip=f[2],nrows=nv-1,comment=comment)
  tail(X)
  dim(X)
  summary(X)

  rgl::plot3d(X[,2:4],aspect=FALSE)
  range(X[,2:4])

  (f=grep("^facets",r))
  (nf=as.numeric(substr( r[f],20,nchar(r[f]))))
  (f=grep("^faces",r))
  (b=grep("^bodies",r)[2])
  b-f
  nf

  F=read.table(f.in,skip=f,nrows = nf-1,header=FALSE,comment=comment)
  tail(F)
  head(F)
  rownames(F)=F[,1]
  head(F)
  F=F[,-1]

  (f=grep("^edges",r))
  (ne=as.numeric(substr( r[f[1]],20,nchar(r[f[1]]))))
  E=read.table(f.in,skip=f[2],nrows =ne-1,fill=TRUE,comment=c("o"))
  tail(E)
  rownames(E)=E$V1
  E=E[,2:3]

  rownames(X)=X[,1]
  X=X[,-1]
  dim(X)

  rgl::mesh3d(NA)->obj

  #str(obj)
  obj$vb=rbind(X[,1],X[,2],X[,3],rep(1,dim(X)[1]))
  rgl::plot3d(t(obj$vb),aspect=FALSE)

  F[1:2,] # indices into E rownames, if not negative
  F=F[,1:3]
  tri=c(0,0,0) # 3 indices (must start from 1) into vertices, correctly ordered
  which(is.na(F))

  range(abs(F[,1:3])) # should equal::
  range(as.numeric(rownames(E)))

  range(abs(E[,1:2])) # should equal::
  range(as.numeric(rownames(X)))

  # lookup table for vertex names
  v.ind=1:dim(X)[1]
  names(v.ind)=rownames(X)
  v.ind[1:10]

  obj$it<-matrix(NA,3,dim(F)[1])
  for (i in 1:dim(F)[1]){
    e=F[i,]
    es=as.matrix(E[as.character(abs(e)),])
    #  if(sum(e<0)>0) stop("e<0 found")
    #if(e[1]<0) es[1,]=es[1,2:1] # swap vertices; maybe repairable in vcgClean
    #if(e[2]<0) es[2,]=es[2,2:1] # swap vertices; maybe repairable in vcgClean
    #if(e[3]<0) es[3,]=es[3,2:1] # swap vertices; maybe repairable in vcgClean
    if (sum(e<0)) es[which(e<0),]=es[which(e<0),2:1]
    #    tri=unlist(unique(c(es[1,],es[2,],es[3,])))
    tri=unique(unlist(c(es))) # c(es[1],es[3],es[5])#
    if(length(unique(tri))<3) {     print(tri);stop();}
    obj$it[1:3,i]<-v.ind[as.character(tri)]
    if (i %% 20 == 0 ) cat(round(100*i/dim(F)[1],1),"\r")

  }


  obj<-Rvcg::vcgClean(obj,sel=1:7)
  rgl::par3d(windowRect=c(1,20,80,82))
  rgl::clear3d()
  #  wire3d(obj,col="white",aspect=FALSE)
  rgl::shade3d(obj,col="red",alpha=0.34)
  rgl::writeOBJ(f.out)
  return(obj)
}

#
# read data from surface evolver dump file
#  special case: faster if SE renumbered before dump
#
#' @export
SErenumbered_2_OBJ <- function(f.in,f.out,comment="o")
{
  r<-readLines(f.in)
  r
  f=grep("^vertices",r)
  f1=r[f[1]]
  (nv=as.numeric(substr( f1,20,nchar(f1))))
  #read.table(f.in,header=FALSE,skip=f[2],nrows=nv-1,comment=comment)[1:10,]
  X=read.table(f.in,header=FALSE,skip=f[2],nrows=nv-1,comment=comment)
  tail(X)
  dim(X)
  summary(X)

  rgl::plot3d(X[,2:4],aspect=FALSE)
  range(X[,2:4])

  (f=grep("^facets",r))
  (nf=as.numeric(substr( r[f],20,nchar(r[f]))))
  (f=grep("^faces",r))
  (b=grep("^bodies",r)[2])
  b-f
  nf

  F=read.table(f.in,skip=f,nrows = nf-1,header=FALSE,comment=comment)
  tail(F)
  head(F)
  rownames(F)=F[,1]
  head(F)
  F=F[,-1]

  (f=grep("^edges",r))
  (ne=as.numeric(substr( r[f[1]],20,nchar(r[f[1]]))))
  E=read.table(f.in,skip=f[2],nrows =ne-1,fill=TRUE,comment=c("o"))
  tail(E)
  rownames(E)=E$V1
  E=E[,2:3]

  rownames(X)=X[,1]
  X=X[,-1]
  dim(X)

  rgl::mesh3d(NA)->obj

  #str(obj)
  obj$vb=rbind(X[,1],X[,2],X[,3],rep(1,dim(X)[1]))
  rgl::plot3d(t(obj$vb),aspect=FALSE)

  F[1:2,] # indices into E rownames, if not negative
  F=F[,1:3]
  tri=c(0,0,0) # 3 indices (must start from 1) into vertices, correctly ordered
  which(is.na(F))

  range(abs(F[,1:3])) # should equal::
  range(as.numeric(rownames(E)))

  range(abs(E[,1:2])) # should equal::
  range(as.numeric(rownames(X)))
  #print(str(E))
  obj$it<-matrix(NA,3,dim(F)[1]) # facets = 3 vertices
  f=as.matrix(abs(F))
  for (i in 1:dim(F)[1]){
    e=f[i,] # 3 edges
    v=E[e,1:2] # fetch vertex ids
    #  print(v)
    #  print(unique(unlist(v)))

    obj$it[1:3,i]<-unique(unlist(v))
    if (i %% 20 == 0 ) cat(round(100*i/dim(F)[1],1),"\r")
  }

  obj<-Rvcg::vcgClean(obj,sel=1:7)
  rgl::par3d(windowRect=c(1,20,80,82))
  rgl::clear3d()
  #  wire3d(obj,col="white",aspect=FALSE)
  rgl::shade3d(obj,col="red",alpha=0.34)
  rgl::writeOBJ(f.out)
  return(obj)
}

# not used
remVertOBJ<-function(O,w) # remove vertices from triangulation
{
  N=dim(O$vb)[2]
  if (max(w)>N ) stop("higher vertex index to remove than present in object")
  length(w)
  O2=O3=O
  O2$vb<-O2$vb[,-(w)]

  # make a lookup table for renumbering
  wn=1:dim(O2$vb)[2] # for surviving vertices to renumber to
  names(wn)=as.character(1:dim(O$vb)[2]) [-w] # filtered original indices

  wa=apply(O2$it,2,function(x) any(x %in% (w)))
  O2$it<-O2$it[,-which(wa)]

  O2$it[]=wn[as.character(c(O2$it[]))]

  Rvcg::vcgClean(O2,sel=1:7)->O2
  return(O2)
}

#
# implements Brechbuehler solver to find initial (u,v)
#
#
# Brechbuehler

#' @export
Brechbuehler.Init.uv.2<-function(X1, Fit_order=12, InitFit=FALSE,poles.axis=2, mat.mode=c("cotan.lapl","cotan.chi","euclid"),file.out="Brechbuehler-init-uv-2.obj",rxy=0.01)
{
  X=X1
  {
    for (k in 1:3) X$vb[k,]=X$vb[k,]-mean(X$vb[k,])
    x=t(X$vb[1:3,])
    #
    # original Brechbuehler 1995
    #
    tri=X$it
    str(tri)
    deg.v=max(table(X$it)) # maximum number of connections
    n.v=dim(X$vb)[2]
    nn=array(NA,c(n.v,deg.v)) #  upto nn neighbours ( == first shell addresses)
    dim(nn)
    for (i in 1:n.v){
      nn.i=c()
      w.i=c(
        tri[,which(tri[1,]==i)],
        tri[,which(tri[2,]==i)],
        tri[,which(tri[3,]==i)])
      nn.i=unique(w.i)
      #cat(i," ",i %in% nn.i,"\n")
      nn.i=setdiff(nn.i,i)
      #nn[i,1:(length(nn.i)-1)] <- nn.i[-which(nn.i==i)]
      nn[i,1:(length(nn.i))] <- nn.i
    }

    X.tri=x

    k.which=poles.axis

    if (k.which<4)
    { x1<-princomp(t(X1$vb[1:3,]))$scores[,1:3]
    plot(x1[,1:2])
    if (k.which==-3) {
      w=(sqrt(x1[,1]^2+x1[,2]^2)<rxy)
      f1p=which.min(x1[w,3])
      f2p=which.max(x1[w,3])
      f1=which(x1[,3]==x1[w,3][f1p])
      f2=which(x1[,3]==x1[w,3][f2p])

      points(x1[w,1:2],col=3,pch=19)
      points(x1[c(f1,f2),1:2],col=2,pch=19,cex=1.3)
    } else {
      f1=which.min(x1[,k.which])
      f2=which.max(x1[,k.which])}

    } else
    {
      cat("\a")
      print("PICK NORTH AND SOUTH POLE centrally in plot")
      princomp((t(X1$vb[1:3,])))$scores->mds
      w.upper=which(mds[,1]>0)
      w.lower=which(mds[,1]<0)
      dev.new();plot(mds[w.upper,2:3],main="PICK NORTH AND SOUTH POLE from PCA (x,y) upper z here")
      p=locator(1)
      f1=which.min( abs(mds[w.upper,2]-p$x[1]) + abs(mds[w.upper,3]-p$y[1]))
      plot(mds[w.lower,2:3],main="PICK NORTH AND SOUTH POLE from PCA (x,y) lower z here")
      points(mds[w.upper,2:3],col="green",pch=".",cex=2.5)
      points(mds[w.upper[f1],2:3],pch=20,cex=1.5,col="red")
      p=locator(1)
      f2=which.min( abs(mds[w.lower,2]-p$x[1]) + abs(mds[w.lower,3]-p$y[1]))
      cat("PICKED (upper/lower indices)",f1," ",f2,"\n")
      g1=c(which(mds[,1]==mds[w.upper[f1],1]),which(mds[,2]==mds[w.upper[f1],2]),which(mds[,3]==mds[w.upper[f1],3]))
      g2=c(which(mds[,1]==mds[w.lower[f2],1]),which(mds[,2]==mds[w.lower[f2],2]),which(mds[,3]==mds[w.lower[f2],3]))
      if (table(g1)==3) g1=g1[1] else stop("Error picking g1")
      if (table(g2)==3) g2=g2[1] else stop("Error picking g2")
      f1=g1;f2=g2;
      cat("PICKED (vertex indices) ",f1," ",f2,"\n")
    }
    (p.fix=c(f1,f2))
    n.m=setdiff(1:n.v,p.fix) # unconstraint points for u

    n=n.v-2
    A=matrix(0,n,n)

    renum=1:n # new index
    names(renum)=as.character(n.m) # from old index

    NN=NN1=list()
    for (i in 1:n.v) {NN[[i]]<-na.omit(nn[i,])
    NN1[[i]]<-setdiff(NN[[i]],p.fix)} # exclude poles from numbering of neighbours
    k=1;
    for (i in n.m)
    { A[k,k]=length( NN[[i]] )
    k=k+1
    }

    k=1;
    for (i in n.m) # original indices n.m
    { A[k,renum[as.character(NN1[[i]])]] <- -1
    k=k+1}
    #image(A)
    range(A-t(A)) # symmetry : 0 0

    b=rep(0,length(n.m)) # rhs for non-poles
    nn.sp=NN[[p.fix[2]]]  # neighb of southpole
    b[renum[as.character(nn.sp)]] <- pi # last line on page 157 (Brechbuehler et al 1995)

    u.s<-solve(A,b)
    u<-rep(0,n.v)
    u[f1]<-0
    u[f2]<-pi
    u[n.m]<-u.s # non-pole solution
    plot(u)

    u.bb.numb=c(0,u.s,pi) # for second phase to compute v

    u.bb=u
  }

  # save u from above
  u.0=u
  # compare against direct solve
  # plot(u.0,u.bb) # perfect

  ia=c(X$it[1,],X$it[2,],X$it[3,]);ja=c(X$it[2,],X$it[3,],X$it[1,])
  N<-max(max(ia),ja); TN <- Matrix::sparseMatrix(dims=c(N,N),i=ia,j=ja,x=rep(1,length(ia)),
                                         use.last.ij=FALSE)
  ig <-igraph::graph_from_adjacency_matrix(TN,mode="undirected",diag = FALSE)
  # recompute  u directly from Laplacian
  B=rep(0,n.v)
  L=igraph::laplacian_matrix(ig,normalized = FALSE,sparse = TRUE)
  # if (mat.mode=="cotan.lapl")   L<-cotan.Laplacian.Matrix(X,L) # overwrite sparse matrix with cotan-Laplacian
  #  if (mat.mode=="cotan.chi")   L<-cotan.Chi.Matrix(X,L) # or with Chi-Laplacian
  L.ret=L# give full Laplacian matrix back for later use
  L[f1,]=0 ; L[f1,f1]=1 ; B[f1]=0  # to fix pole to zero
  L[f2,]=0 ; L[f2,f2]=1 ; B[f2]=pi # to fix pole to pi ;
  # later this is 2pi for East reference in v computation

  solve(L,B)->u ##    re-solve for u
  plot(u,u.0,main="compare Brechbuehler and Laplace-solver") # compare - diffenernce due to alternative Laplacian
  #
  # solve v with NS on v=0 and West on v=2*pi-epsilon

  # create set of vertices on data-line N-pole to S-pole
  (north.south <- igraph::shortest_paths(ig,from=f1,to=f2,output="vpath")$vpath[[1]])
  # better take largest increases in u as north-south
  #   check halo for better path
  #
  # no improvements by the following observed:
  if(TRUE)
    for (i in 1:3){
      here=f1
      north.south.new=as.numeric(north.south)
      cnt=2 # start to set 2nd
      while(here != f2)
      {
        n.here=as.numeric(neighbors(ig,here))
        take=n.here[which.max(u[n.here])]
        north.south.new[cnt]=take
        here=take
        cnt=cnt+1 # next
      }
      north.south=north.south.new
    }

  # check visually
  halo=unique(unlist(ego(ig,order=1,north.south)))
  pure.halo=setdiff(halo,north.south)
  col=rep("black",n.v) ; col[north.south]="green"
  col[pure.halo]="red"
  rgl::wire3d(X,col=col,lwd=2)

  # try to do it similar to Brechbuehler
  EW=rep("",n.v)
  here=f1
  cnt=2 # start to set 2nd
  done=f1
  while(here != f2)
  {
    n.here=as.numeric(neighbors(ig,here))
    take=n.here[which.max(u[n.here])] # next
    n.take=neighbors(ig,take)
    done=c(done,take)

    (check=intersect(n.here,n.take))
    #check=setdiff(check,EW!="") # already marked not to check
    dV=X$vb[1:3,take]-X$vb[1:3,here]
    dN=X$normals[1:3,here]
    for (chk in check)
    {dS=X$vb[1:3,chk]-X$vb[1:3,here]
    D=det(matrix(c(dV,dS,dN),3,3)) # this is my solution to classify East vs. West
    EW[chk]=ifelse(D<0,"E","W")
    }
    #  EW[setdiff(n.here,check) ] = " " # mark as not to check; maybe only for f1 (N)?
    here=take
    cnt=cnt+1 # next
  }
  table(EW)
  sum(EW!="")

  rgl::clear3d()
  rgl::wire3d(X,col=col,lwd=2)
  rgl::spheres3d(t(X$vb[1:3,EW=="W"]),rad=0.05,col="green")
  rgl::spheres3d(t(X$vb[1:3,EW=="E"]),rad=0.05,col="blue")

  all(which(EW!="") %in% halo) # some more in halo to check for their neighbour on the date line!
  rest=setdiff(pure.halo,which(EW!=""))
  length(rest)
  length(pure.halo)

  # just fill up with "E" between two "E"
  East=which(EW=="E")
  West=which(EW=="W")

  # remove north/southpoles neighbours that are not to check
  n.n=neighbors(ig,f1)
  mark=setdiff(n.n,which(EW!=""))
  EW[mark] =" "
  n.s=neighbors(ig,f2)
  mark=setdiff(n.s,which(EW!=""))
  EW[mark] =" "

  #  open3d()
  #  wire3d(X,lwd=1)

  #  spheres3d(t(X$vb[1:3,West]),rad=0.05,col="green")
  #  spheres3d(t(X$vb[1:3,East]),rad=0.05,col="blue")
  #  spheres3d(t(X$vb[1:3,rest]),rad=0.06,col="black")

  # for vertices with only a single linkage to north wets:
  #    take EW from the neighbours EW
  #

  rest=setdiff(rest,c(f1,f2))

  for (i in rest) # rest has no poles
  { print(i)
    n.i=neighbors(ig,i)
    ew=EW[n.i];
    tb=table(ew)
    print(tb)
    if (any(names(tb) %in% c("W","E")))
    {
      if (("E" %in% names(tb)) & ("W" %in% names(tb)))
      {cat("cannot decide EW in rest (chose W): ",i,"\n");EW[i]="W"} else
        EW[i]=names(tb[names(tb) %in% c("W","E")])
    }
  }

  # open3d()
  #  wire3d(X,col=col,lwd=2)
  #  spheres3d(t(X$vb[1:3,EW=="W"]),rad=0.05,col="green")
  #  spheres3d(t(X$vb[1:3,EW=="E"]),rad=0.05,col="blue")

  #update East and West indices
  East=which(EW=="E")
  West=which(EW=="W")


  #
  # now we have marked points from pure halo (that are connected with north-south-line) for being left or right from line
  # time to use EW
  #

  # assemble linear problem from scratch for v
  L=igraph::laplacian_matrix(ig,sparse=TRUE) # or re-use L.ret from above

  B=rep(0,n.v)
  for (i in c(f1,f2)) {L[i,which(L[i,]<0)]<-0;L[i,i]<-1} # fix poles v
  B[f2]=pi # not actually needed; set all West to 2pi-epsilon
  B[f1]=pi
  epsilon=pi/length(north.south)


  inner=setdiff(north.south,c(f1,f2))
  for (i in inner) {L[i,which(L[i,]<0)] <-0 ; L[i,i]<-1; B[i]=0} # fix north-south to v=0 (B[inner]=0)

  for (i in West){ L[i,which(L[i,]<0)] <-0 ; L[i,i]<-1; B[i]=2*pi-epsilon}

  # TO CHECK : THIS CODE REMOVAL
  #    w=which(L[i,]<0) # what is i connected with?
  #    w1=intersect(w,north.south) # restrict to north-south
  #    cat(w1,"\n")

  #  if (length(w1)>0) L[i,w1]<-0; L[i,f2]=-length(w1) # instead of ns point to south,
  #  where v=2 pi implies that some connections get lost because
  #    so use |w1| as entry to compensate
  # this is would not work with cotan matrices?
  #
  #alternative:
  #   move weights into L_j,f2 one by one

  B[f2]=pi

  solve(L,B)->v
  v=as.numeric(v)
  cat("RANGE RESID:")
  print(abs(range(L%*%v-B)))
  plot(u,v,main="Brechbuehler initial")
  points(u[East],v[East],pch=19,col=2)
  points(u[West],v[West],pch=19,col=3)
  points(u[inner],v[inner],pch=19,col=4)

  X$texcoords=rbind(as.numeric(u),as.numeric(v))

  #shade3d(X, col="white",meshColor = "vertices")
  #  contourLines3d(X,fn=X$texcoords[1,],100,lwd=1)
  #open3d()
  #shade3d(X, col="white",meshColor = "vertices",texcoords=cbind(u,v))
  #  contourLines3d(X,fn=X$texcoords[2,],100,col="red")

  # jump at date line from inner to East generates quite dense countours
  #   which is correct !!!

  #  spheres3d(t(X$vb[1:3,EW=="W"]),rad=0.05,col="green")
  #  spheres3d(t(X$vb[1:3,EW=="E"]),rad=0.05,col="blue")
  #  open3d()
  #  clear3d()
  #  shade3d(X, col="white",meshColor = "vertices",textu="checkers_fine.png")

  # try a fit then look for problems/overhanging triangles
  if (InitFit)
  {
    uv=t(X$texcoords)
    dim(uv)

    bas.i<-MakeBasis_0_irreg(Fit_order,uv[,1],uv[,2])
    dim(bas.i$Ylm)

    a_v<-VertexAreasOBJ(X)

    A.init<-FitAlm_Tikhonov(x,bas.i,lambda = 0.1)#,weights=a_v^3)
    grd=MakeGrid_GaussLegendreSimpson(50)
    bas=MakeBasis_UV(Fit_order,grd$U,grd$V)
    updateX(A.init,grd,bas)->C
    plot3d(C$X,asp=F)
    range(A.init)
    open3d()
    wire3d(X)
    writeOBJ(file.out)
  } else A.init<-NULL

  return(list(uv=cbind(u,v),OBJ=X,L.full=L.ret,A.init=A, East=East, West=West,Inner=inner,f1=f1,f2=f2,igraph=ig,A.init=A.init,poles=p.fix))
}


#' @export
Obj2ObjQ<-function(O,grd)
{
  nx=grd$nu;ny=grd$nv;
  q=matrix(NA,4,nx*ny);k=0
  for (i in 1:(nx-1))  for (j in 1:(ny-1)){
    k=k+1;
    l=(j-1)*nx+i
    q[1,k]=l
    q[2,k]=l+1
    q[3,k]=l+1+nx
    q[4,k]=l+nx
  }
  # could make a closed form removing double vertices, but this woould complicate updating coords
  rgl::qmesh3d(O$vb,indices = q[,1:k], normals = matrix(0,3,k) )-> M
  return(M)
}

#' @export
X2ObjQ<-function(O,X)
{ O$vb=rbind(t(X),1) # ingest coordinates in 3d-graphics objec
return(O)  # quads dont allow for vcgUpdatenormals for unknown reason
}



#
# create an Unduloid 3d object for a fraction or multiple of periods
#
#' @export
TriMesh_Unduloid<-function(a=1,c=2,periods=1.0,nx=40,ny=40,shade=TRUE,wire=FALSE,clean=TRUE)
{ m=(c^2-a^2)/2;n=(c^2+a^2)/2
  mu=2/(a+c);k2=(c^2-a^2)/c^2
  ulimup=(pi/2+pi/4)*2/mu; ulimdown=pi/4*2/mu
  p=ulimup+ulimdown
  u=pracma::linspace(-ulimdown,-ulimdown+p*periods,nx); phi=u*mu/2-pi/4
  x=Re(a*Carlson::elliptic_F(phi,k2)+c*Carlson::elliptic_E(phi,k2))
  z=Re(sqrt(m*sin(mu*u)+n))

  d=2*pi/ny
  r_mat=matrix(c(cos(d),sin(d),-sin(d),cos(d)),2,2)
  X1=X=cbind(0,z)
  Y=array(NA,c(length(x),ny+1,3))
  for (i in 1:(ny+1))
  { Y[,i,1]=x
  Y[,i,2]=X1[,1]
  Y[,i,3]=X1[,2]
  X1=X1%*%r_mat
  }
  q=matrix(1,3,nx*ny*2);k=0
  for (i in 1:(nx-1))  for (j in 1:ny){
    k=k+1;l=(j-1)*nx+i
    q[1,k]=l
    q[2,k]=l+1
    q[3,k]=l+1+nx
    k=k+1
    q[1,k]=l
    q[2,k]=l+nx+1
    q[3,k]=l+nx
  }
  # shift to x=0 for leftmost slice
  rgl::mesh3d(x=Y[,,1]-min(Y[,,1]),y=Y[,,2],z=Y[,,3],triangles=q)->M
  Rvcg::vcgUpdateNormals(M)->M

  if (shade) rgl::shade3d(M,alpha=0.6,col="white")
  if (wire) rgl::wire3d(M)
  if (clean) Rvcg::vcgClean(M)->M
  attr(M,"H")=1/(a+c) # store theoretical value of constant mean curvature
  # the attribute is probably not inherited in derived objects.
  return(M)
}

#
#    use polynomials of degree 12 to fit to lowess fits with other data
#
#' @export
Lowess_vcg_meanvbOBJ<-function(x,O) # return corrected
{
  LW=list()
  LM=list()
  crv=Rvcg::vcgCurve(O)
  tb=table(O$it)
  kr=unique(tb)
  print(kr)
  plot(0,0,col=0,xlim=range(x),ylim=range(-crv$meanvb))
  for (k in kr){
    LW[[k]] = lowess(x[tb==k],-crv$meanvb[tb==k],f=0.1)
    points(LW[[k]],type="l",lwd=3,col=k-3)
    xx=LW[[k]]$x;yy=LW[[k]]$y
    LM[[k]] = polyfit(xx,yy,12)
    points(xx,polyval(LM[[k]],xx),col=k-3,cex=1.5,lwd=2)
  }
  X=Y=Z=rep(0,length(x)) # unordered return
  LX=LY=LK=list()
  for (k in kr){
    if (k!=6) LY[[k]]=-crv$meanvb[tb==k] - polyval(LM[[k]],x[tb==k]) + polyval(LM[[6]],x[tb==k]) else LY[[6]]=-crv$meanvb[tb==6]
    points(x[tb==k],LY[[k]],col=k-3)
    Y[tb==k]=LY[[k]]
    LX[[k]]=x[tb==k]
    X[tb==k]=LX[[k]]
    LK[[k]]=rep(k,sum(tb==k))
    Z[tb==k]=k
  }

  return(list(x=unlist(LX),y=unlist(LY),k=unlist(LK),X=X,Y=Y,K=Z)) # k-sorted output, not good for coloring 3d object
}





# Filters:

# must be:
# Z is zero mode
# X is +1 mode
# Y is -1 mode

# keep only lowest cos(i*v) and sin(i*v) terms, i.e. i=1
# and make circular crosssections
#' @export
Filter_Z_AxiSymm<-function(A,bas)
{if (!is.matrix(A)) A=LM2A(A,bas)
  A1=A; M=bas$LM[,2]; L=bas$LM[,1]
  for (l in bas$Lset){
    #take the mean of cos/sin coeffs, make structure
    #  really axisymmetric
    a=(A[L==l & M==-1,2] + A[L==l & M==+1,1])/2
    A1[  L==l & M==-1,2] = A1[L==l & M==+1,1] = a

    A1[L==l & M!=-1, 2]=0 # sparsity pattern for axisymm.
    A1[L==l & M!=+1, 1]=0 # sparsity "
    A1[L==l & M!= 0, 3]=0 # pole on Z axis
  }
  return(A1)
}


# needs more thought with signs of a ...
#' @export
Filter_Reduced_AxiSymm<-function(A,bas)
{ A=matrix(A,ncol=3)
  for (i in 1:bas$L_max)
  {
   a=(A[(i-1)*3+1,2] - A[(i-1)*3+3,1])/2
   A[(i-1)*3+1,2] = a
   A[(i-1)*3+3,1] = a
  }
  A[c(1,2),3]   =0
  A[c(1,3),2]   =0
  A[c(2,3),1]   =0
return(A)
}





#
# keep only abs-maximal entries per (l,m)
#   - removes partially correlated coordinates like X=1*sin, Y=0.1*sin
#
#' @export
Filter_A_m<-function(A,bas,max_per_l=1)
{
  A1=A; A1[]<-0
  L=bas$LM[,1];M=bas$LM[,2]
  Ls=which(diff(L)>0)
  Ls=c(0,Ls)
  ll=0;
  for (l in bas$Lset){
    ll=ll+1
    for (k in 1:3)
    {w=rev(order(abs(A[L==l,k])))[1:max_per_l]
    #     print(A[w+Ls[ll],k]);
    A1[w+Ls[ll],k]=A[w+Ls[ll],k]
    }
    # cat("\n")
  }
  return(A1)
}


# reduce basis according to zero-rows in filtered A1
#' @export
Filter_Basis<-function(bas,A,A1)
{
# modify bas according to zero entries in A1
# also give back a reduced version of A to continue work with
del=c() # deletion candidates
for(i in 1:bas$Ai_max) if (all(A1[i,]==0)) del=c(del,i)
#cat("Filtering for ",del,"\n")
rownames(bas$A)[del]->exclude
# now reduce
bas$Ylm <-bas$Ylm[,-del]
bas$Ylm_u <-bas$Ylm_u[,-del]
bas$Ylm_v <-bas$Ylm_v[,-del]
bas$Ylm_vv<-bas$Ylm_vv[,-del]
bas$Ylm_uu<-bas$Ylm_uu[,-del]
bas$Ylm_uv<-bas$Ylm_uv[,-del]
bas$LM<-bas$LM[-del,]
bas$l=bas$LM[,1]
bas$m=bas$LM[,2]
bas$Lset=unique(bas$l)
bas$Mset=unique(bas$m)
bas$Ai_max=dim(bas$LM)[1]
bas$comment=paste("m-filtered basis, excluded ",exclude)
A<-A1[-del,]
bas$A<-A
bas$G.tk<-bas$G.tk[-del]
return(list(bas=bas,A=A))
}

#
# dense regions in u,v can be stretched by this
#  - coserves triangulation quality
#   needed for postprocessing Brechbühler
#
#' @export
spreadout.uv<-function(uv)
{
  N=dim(uv)[1]

  u.order=order(uv[,1])
  v.order=order(uv[,2])

  uv.sort=apply(uv,2,sort)
  matplot(uv.sort,type="l")

  uv1=uv.sort
  min(diff(uv1[,1]))->small
  diff(range(diff(uv1[,1]))) -> scale
  k=0
  while(min(diff(uv1[,1]))<scale)
  {  i=which.min(diff(uv1[,1]))
  uv1[1:N > i,1] = uv1[1:N > i,1]+scale# *scale.fac
  cat("cycle U ",k,"\r");k=k+1
  }
  plot(uv1[,1])
  uv1[,1]=uv1[,1]/diff(range(uv1[,1]))*pi
  i.order=1:N
  i.order[u.order]=1:N
  # uv=t(X1$texcoords)
  uv=uv1
  uv[1:N,1]=uv1[i.order,1]
  plot(uv)

  min(diff(uv1[,2]))->small
  diff(range(diff(uv1[,2]))) -> scale
  k=0
  while(min(diff(uv1[,2]))<scale)
  {  i=which.min(diff(uv1[,2]))
  uv1[1:N > i,2] = uv1[1:N > i,2]+scale#*scale.fac
  cat("cycle V ",k,"\r");k=k+1
  }
  plot(uv1[,2])
  uv1[,2]=uv1[,2]/diff(range(uv1[,2]))*2*pi
  i.order=1:N
  i.order[v.order]=1:N
  uv[1:N,2]=uv1[i.order,2]
  plot(uv)
  return(uv)
}


#' @export
ConsIter<-function(A,grd,bas,C,g2, Ctol=1e-3, nsteps=20,
                   prn=FALSE,del_cons=0.3,nm=FALSE,do_one=TRUE)
{
  l=0; Nc=bas$Nc
  updateX(A,grd,bas)->C
  h2=E_SCM(A,grd,bas,C)
  Cons_RHS <-ConsRHS(h2,bas)
  NCons<-max(abs(Cons_RHS)/bas$Target)
  if(prn) cat(crayon::yellow(l,": Cons_%:"), crayon::cyan(round(100*Cons_RHS/bas$Target,4)),crayon::yellow(" |Cons|:"),ifelse(NCons>Ctol,crayon::red(NCons),crayon::green(NCons)),"\n")

  if (NCons<=Ctol) {
    E_SCM(A,grd,bas,C) -> h2
    Grad_SCM(h2,grd,bas,C) -> g2

    Cons_RHS <-ConsRHS(h2,bas)
    NCons<-max(abs(Cons_RHS)/bas$TNorm)
  };
  sol=rep(0,Nc) # default to return
  while ((NCons>Ctol & l<nsteps) | ( l==0 )){
    #if (l==0 & !do_one) break; # not impolemented;
    NCons1=NCons
    dFm<-c(g2[[ bas$Cons[1] ]])
    for (ii in 2:Nc) dFm=cbind(dFm,c(g2[[bas$Cons[ii]]])) # additional constraints
    M_c<-matrix(0.0,Nc,Nc); for (i in 1:Nc) for (j in 1:Nc) M_c[i,j]<-dot2(dFm[,i],dFm[,j])
    Pm<-pracma::pinv(M_c)
    sol<- (- Pm %*%Cons_RHS)[,1]
    names(sol)=names(bas$Cons)
    delta <- sol[1]*dFm[,1] ; for (i in 2:Nc) delta<-delta + sol[i]*dFm[,i]
  #  delta <- Filter_1_delta(delta) # fix the orientation of L=1 (ellipsoid)
    A <- A + del_cons * delta
    # A<-Filter_1_A(A) # not needed if delta was filtered and A is filtered before

    updateX(A,grd,bas) -> C
    #  plotA_l(delta,bas,bar=TRUE,ylab=expression(~delta*A[cons]),main="CONS")
    E_SCM(A,grd,bas,C) -> h2
    Grad_SCM(h2,grd,bas,C) -> g2

    Cons_RHS <-ConsRHS(h2,bas)
    NCons<-max(abs(Cons_RHS)/bas$TNorm)
    if (prn) cat(crayon::green(l," ",round(NCons,9)),"\r");
    l=l+1
    if(NCons1==NCons) break
  }
  if(prn) cat(crayon::green(l,": Cons_%:"), crayon::blue(round(100*Cons_RHS/bas$Target,4)),crayon::green(" |Cons|:"),ifelse(NCons>Ctol,crayon::red(NCons),crayon::green(NCons)),"\n")

  return(list(A=A,h2=h2,g2=g2,sol=sol,NCons=NCons,cons_iter=l,
               Lambda=sol,Cons_RHS=Cons_RHS,C=C))
}

#
# Hessian comp. from finite differences of gradients (forward)
#   constraints are treated outside with Langrangian, involving Jacobians of constraint functions, registered in basis as strings of gradient names (bas$QCons)
#

# have seen in a L=12 stomatocyte CNM that axis is tilting
# so retry with this flag TRUE
glob_fltZ=TRUE

#' @export
FullModelHessian<-function(A,grd,bas,Ref,del=1e-5,Ctol=1e-3,nm=FALSE)
{ tictoc::tic()
  Ai_max=bas$Ai_max
  C=updateX(A,grd,bas)
  h20=E_SCM(A,grd,bas,C)
  S=SEN(A,grd,bas,Ref,h20)
  ES=E_SEN(A,grd,bas,S,Ref)
  Gh20=Grad_SCM(h20,grd,bas,C)
  GS0=Grad_SEN(A,grd,bas,Gh20,S,Ref)
  if (glob_fltZ) G0=c(Filter_1_delta(Gh20$grad_SCM) + Filter_1_delta(GS0$grad_SEN)) else G0=c(Gh20$grad_SCM + GS0$grad_SEN) # given back as $G
  if(nm) G0=make_delta_normal_to_surface(G0,grd,bas,h20$n)
  H=matrix(0, Ai_max*3, Ai_max*3)
  for (j in 1:(3*Ai_max))
  { if (j %% 10==0) cat(" ",round(j/(3*Ai_max)*100,1),"\r")
    A1=A; A1[j]=A1[j]+del
    C=updateX(A1,grd,bas)
    h2=E_SCM(A1,grd,bas,C)
    S=SEN(A1,grd,bas,Ref,h2)
    Gh2=Grad_SCM(h2,grd,bas,C)
    GS=Grad_SEN(A1,grd,bas,Gh2,S,Ref)
    if (glob_fltZ) Gj=c(Filter_1_delta(Gh2$grad_SCM) + Filter_1_delta(GS$grad_SEN)) else
          Gj=c(Gh2$grad_SCM + GS$grad_SEN)
    if(nm) Gj=make_delta_normal_to_surface(Gj,grd,bas,h20$n)
    H[,j]= (Gj-G0)/del
  }
  tictoc::toc()
  return(list(H=(H+t(H))/2,H_fd=H,G=G0,GS=GS,g2=Gh20,
              E=h2$Wb+ES,ES=ES,Wb=h2$Wb,h2=h20,
              C=C,gradC=Gh20$gradC,gradA=Gh20$gradA,
              gradV=Gh20$gradV,A=A))
}

#' @export
ID<-function(A,bas)
{ return(A)
}


#' @export
FullModelHessian_Par<-function(A, grd, bas, Ref, del=5e-6, Ctol=1e-3, nm=FALSE, Mem_mc.cores = 4, filter_grad=ID, timing=TRUE )
{
#  cat("M.C0(Driver) ",M.C0,"P",eval(M.C0,envir = .GlobalEnv),"\n")

  Ai_max=bas$Ai_max
  C=updateX(A,grd,bas)
  .M.C0<<-C
  .M.Ref<<-Ref
  .M.bas<<-bas
  .M.grd<<-grd # these will be shared by cluster export; hence must reside in upper environment
  h20=E_SCM(A,grd,bas,C)
  S=SEN(A,grd,bas,Ref,h20)
  ES=E_SEN(A,grd,bas,S,Ref)

  Gh20=Grad_SCM(h20,grd,bas,C)
  GS0=Grad_SEN(A,grd,bas,Gh20,S,Ref)
  G0=c( Gh20$grad_SCM + GS0$grad_SEN ) # given back as $G
  if(nm) G0=make_delta_normal_to_surface(G0,grd,bas,h20$n)
  H=matrix(0, Ai_max*3, Ai_max*3) # to be assembled from parallel vectors
  Lpar=list()
  if(timing)tictoc::tic()
  for (j in 1:(3*Ai_max))
  { A1=A; A1[j]=A1[j]+del
#  C=updateX(A1,grd,bas) # update now runs in client
  Lpar[[j]]=list(A=A1,#grd=grd,bas=bas,Ref=Ref,
                 M.C0=M.C0,M.K_ADE=M.K_ADE,M.K_b=M.K_b,
                 M.mu=M.mu,M.Ka=M.Ka,M.a3=M.a3,M.a4=M.a4,
                 M.b1=M.b1,M.b2=M.b2,
                 M.Rcpp=TRUE,M.Rcpp_ncores=M.Rcpp_ncores,index=j)
  }

  if(timing){cat("paralleliz. preperation  ");  tictoc::toc();}
  #cat(names(Lpar[[1]]),"\n")
  #LINUX:    mclapply(Lpar,FullHessian_Client,mc.cores=Mem_mc.cores)->LH # not done
  #WINDOWS:  make a socket cluster
  if(timing)tictoc::tic()
  {
  cat("setup Cluster\n")
  cl<-parallel::makeCluster(Mem_mc.cores)
  # need to communicate int2d_...cxx for H2_grads (c++) to work on cluster
  parallel::clusterExport(cl,varlist=c("M.C0","M.K_ADE","M.K_b","M.mu","M.Ka","M.a3","M.a4","M.b1","M.b2",
                                       "M.Rcpp","M.Rcpp_ncores",".M.grd",".M.bas",".M.Ref","int2d_matrix_cxx","vectomat_cxx","vectoarr_cxx"))

    if(timing){cat("cluster startup ");  tictoc::toc();
    tictoc::tic()}
    LH<-parallel::parLapply(cl, Lpar, FullHessian_Client , FALSE)
    parallel::stopCluster(cl)
  }
  if(timing){cat("parallel lapply ");  tictoc::toc()}
  for (i in 1:(3*Ai_max))
  {
    H[,i]= c( filter_grad( LH[[i]] - G0, bas) )/del
  }

  return(list(H=(H+t(H))/2,H_fd=H,G=G0,g2=Gh20,
              E= h20$Wb + ES,ES=ES,Wb=h20$Wb,h2=h20,
              C=C,gradC=Gh20$gradC,gradA=Gh20$gradA,
              gradV=Gh20$gradV,A=A))
}

# no export - internal to FullHessian_Par
FullHessian_Client<-function(L,DBG=FALSE) # Lmax=13 takes 4 seconds per call
{
  M.Rcpp=L$M.Rcpp;M.Rcpp_ncores=L$M.Rcpp_ncores
  if(DBG)cat(tictoc::toc()[[4]],":",paste(err),"\n",file="setup.txt",append = TRUE);

  Rcpp::sourceCpp("data/MembraneRBC.cpp")

  C=updateX(L$A,.M.grd,.M.bas) # no longer in L
  h2=E_SCM(L$A,.M.grd,.M.bas,C)
  S=SEN(L$A,.M.grd,.M.bas,.M.Ref,h2)

  # decide by M.Rcpp for openmp-parallel code
  #       cant directly use Grad_SCM on cluster due to scattered objects names
  if(!M.Rcpp) {Gh2=Grad_SCM_R(h2,.M.grd,.M.bas,C)} else {
    G2=Grad_SCM_cxx(h2,.M.grd, .M.bas, C, M.C0, L$M.Rcpp_ncores, L$M.K_b, L$M.K_ADE)

    # needed for SEN() below;
    Gh2=list(ddA=array(G2$ddA,c(.M.grd$ndof,.M.bas$Ai_max,3)),
             ddV=array(G2$ddV,c(.M.grd$ndof,.M.bas$Ai_max,3)),
             grad_SCM=G2$grad_SCM,gradV=G2$gradV,gradA=G2$gradA,gradC=G2$gradC,
             dE=array(G2$dE,c(.M.grd$ndof,.M.bas$Ai_max,3)),
             dF=array(G2$dF,c(.M.grd$ndof,.M.bas$Ai_max,3)),
             dG=array(G2$dG,c(.M.grd$ndof,.M.bas$Ai_max,3)))
    }

  GS=Grad_SEN(L$A,.M.grd,.M.bas,Gh2, S,.M.Ref)
  Gj=c(Gh2$grad_SCM + GS$grad_SEN)

  return(Gj)
}


# could become interesting to updateX only the terms required for Hesssian computation
# rather than re-updating for every row of H.
#' @export
synth12<-function(A,C,i,j,k,del) # spatial k ~ X,Y,Z
{
 cat(sum(C$X_u),"-> \t")
 n=dim(A)[1]
 if (i>0)
   {
    C$X[,k] = C$X[,k] - bas$Ylm[,i]*del # remove del term
    C$X_u[,k] = C$X_u[,k] - bas$Ylm_u[,i]*del # remove del term
    C$X_v[,k] = C$X_v[,k] - bas$Ylm_v[,i]*del # remove del term
    C$X_uu[,k] = C$X_uu[,k] - bas$Ylm_uu[,i]*del # remove del term
    C$X_uv[,k] = C$X_uv[,k] - bas$Ylm_uv[,i]*del # remove del term
    C$X_vv[,k] = C$X_vv[,k] - bas$Ylm_vv[,i]*del # remove del term
   }
  if(j<=n)
    {
    C$X[,k] = C$X[,k] + bas$Ylm[,j]*del # add del term
    C$X_u[,k] = C$X_u[,k] + bas$Ylm_u[,j]*del # add del term
    C$X_v[,k] = C$X_v[,k] + bas$Ylm_v[,j]*del # add del term
    C$X_uu[,k] = C$X_uu[,k] + bas$Ylm_uu[,j]*del # add del term
    C$X_uv[,k] = C$X_uv[,k] + bas$Ylm_uv[,j]*del # add del term
    C$X_vv[,k] = C$X_vv[,k] + bas$Ylm_vv[,j]*del # add del term
  }
#else {
#    C$X[,k+1] = C$X[,k+1] + bas$Ylm[,j]*del # add del term
#    C$X_u[,k+1] = C$X_u[,k+1] + bas$Ylm_u[,j]*del # add del term
#    C$X_v[,k+1] = C$X_v[,k+1] + bas$Ylm_v[,j]*del # add del term
#    C$X_uu[,k+1] = C$X_uu[,k+1] + bas$Ylm_uu[,j]*del # add del term
#    C$X_uv[,k+1] = C$X_uv[,k+1] + bas$Ylm_uv[,j]*del # add del term
#    C$X_vv[,k+1] = C$X_vv[,k+1] + bas$Ylm_vv[,j]*del # add del term
#  }

 cat(sum(C$X_u),"\n")
 return(C)
}
#' @export
mat2vec<-function(m) return(c(m))
#' @export
matdiff2vec<-function(m1,m2)
{return(c(m1-m2))}

#' @export
symmetrize<-function(m)
{return((m+t(m))/2.0)}
#' @export
matadd2vec<-function(m1,m2)
{return(c(m1+m2))}

#' SetConstraints
#'
#' set the constraint target values and store information in basis.
#' Without parameters Cons, ..., the standard area and volume constraints are set with values (140, 100).
#' If you set a third constraint, for CNM you have to give M$Lambda a third component.
#' @param bas basis to modify constraints, e.g. M$bas (M MemRBC object)
#' @param Cons vector of character of constraint gradients, from "gradA","gradV" and "gradC"
#' @param QCons vector of constraint names, from "Area","Volume","Curv", same order as Cons
#' @param Target vector of constraint values, like c(140,100,88) for Area, Volume and Curv
#' @param TNorm vector of normalization constants for calculating the norm, i.e. degree of total constraint violation
#' @examples
#' data(D5); SetParams(D5)
#'
#' # add curvature constraint:
#' SetConstraints(D5$bas,Cons=c("gradA","gradV","gradC"),
#'   QCons=c("Area","Volume","Curv"),
#'   Target=c(140,100,121),
#'   TNorm=c(140,100,121)) -> D5$bas  # store modified basis back into membranes D5 basis
#'
#' # minimize with steepest descend under Rosen Constraint Projection
#'
#' SDRC(D5,1000)->D5sdrc
#'
#' # pair-plots of target quantities and energy E
#' plot(D5sdrc$SDRC_Sample[c("E","A","V","C")])
#'
#' @export
SetConstraints<-function(bas,Cons=c("gradA","gradV"),
                         QCons=c("Area","Volume"),
                         Target=c(140, 100),
                         TNorm=c(140, 100)) # any further (implemented) gradients allowed
{
  if(class(bas)=="MemRBC)") {M=bas;bas=M$bas;toMemRBC=TRUE} else toMemRBC=FALSE
  bas$Cons=Cons
  bas$Nc=length(Cons)
  bas$Target=Target
  bas$QCons=QCons
  bas$TNorm=TNorm
  names(bas$Cons)= QCons
  names(bas$QCons)= QCons
  names(bas$TNorm)= QCons
  names(bas$Target)= QCons
  message("SetConstraints: for CNM, remember to set M$Lambda with bas$Nc, eg M$Lambda=rep(0.1,M$bas$Nc) \n")
  if (!toMemRBC) return(bas) else {M$bas<-bas; return(M)}
}

#' @export
ConstraintHessian<-function(H,bas,Lambda,filter=ID)
{ dH=dim(H$H)[1] # should be equal Ai_max*3

Nc=bas$Nc
if(Nc==0) {message("No constraints for ConstraintHessian - return Hessian as is\n");return(H)}
if(length(Lambda)!=Nc) stop("wrong number of Lagrangian lambdas vs. registered constraints bas$Nc")
M=matrix(0,dH+Nc,dH+Nc)
M[1:dH,1:dH]=H$H # main block
for (i in 1:bas$Nc) {
  M[dH+i,1:dH]=M[1:dH,dH+i]=filter(H[[ bas$Cons[i] ]],bas)
  #else
  #  M[dH+i,1:dH]=M[1:dH,dH+i]=H[[ bas$Cons[i] ]]# off-diagonal blocks = Jacobian of constraint functions
}
return(M)
}

#' @export
ConsJacobian<-function(g2,bas)
{
  if (glob_fltZ) RosenA<-c(Filter_1_delta(g2[[bas$Cons[1]]])) else RosenA<-c(g2[[bas$Cons[1]]]) # first gradient, bas$Cons has names of gradients
  for (ii in 2:bas$Nc) # further gradients to add
    if (glob_fltZ) RosenA<-rbind(RosenA,c(Filter_1_delta(g2[[ bas$Cons[ii] ]]) )) else RosenA<-rbind(RosenA,c(g2[[ bas$Cons[ii] ]] ))
    return(RosenA)
}

#
# Rosen (1961) projection of gradient
#
#' @export
RosenProjection<-function(G,g2,bas) # output the projected energy gradient G; g2 contains needed gradients of constraint functions
{
  RosenA<-ConsJacobian(g2,bas)
  RosenAAt<-RosenA%*%t(RosenA)
  IRosenAAt<-pracma::pinv(RosenAAt) # pinv tolerates also gradients to be zero
  lambda<- t(RosenA)%*%IRosenAAt%*%RosenA
  Gprime <- c(G) - (lambda %*% c(G)) [,1]
  lambdaG=c()
  for (i in 1:bas$Nc) lambdaG[i]=dot2(c(G),RosenA[i,]) # how much from each constraint gradient  is projected out of Gradient?
  names(lambdaG)=bas$Qcons
  return(list(Gprime=G,Eigs=eigen(RosenAAt)$values,lambdaG=lambdaG))
}


#' @export
ConsRHS<-function(h2,bas)
{  Cons_RHS <- rep(0,bas$Nc);names(Cons_RHS)=bas$QCons
for (k in 1:bas$Nc) Cons_RHS[k] <- h2[[bas$QCons[k]]] - bas$Target[k]
return(Cons_RHS)
}


# if you decide to erase dot2 from Rcpp-code:
if(!exists("dot2")) dot2=function(x,y) sum(x*y)

  #
  # filters for a fixed orientation of L=1 ellipsoid shape
  #   keeps Z-rotation zero (phase=0)
  #
#' @export
Filter_1_A<-function(A)
  {
    A[1,c(2,3)]=0 # m=-1
    A[2,c(1,2)]=0 # m=0
    A[3,c(1,3)]=0    # m=+1
    return(A)
  }

#' @export
Filter_1_delta<-function(delta)
  { d=matrix(delta,ncol=3)
  d[1,c(2,3)]=0
  d[2,c(1,2)]=0
  d[3,c(1,3)]=0
  return(c(d))
  }

#' @export
Cons_filter_delta<-function(delta,grd,bas,H2,nm=TRUE)
  {
    if (nm)  delta=make_delta_normal_to_surface(delta,grd,bas,H2$n)
    delta=delta/(1+ll.a*sqrt(bas$G.tk))
    return(matrix(delta,bas$Ai_max,3))
  }

  Rosen_filter_delta<-function(delta,grd,bas,H2,nm=TRUE)
  { if(nm)  delta=make_delta_normal_to_surface(delta,grd,bas,H2$n)
  delta=delta/(1+ll.a*sqrt(bas$G.tk))
  return(matrix(delta,bas$Ai_max,3))
  }

  make_delta_normal_to_surface <- function(delta, grd, bas, n){
    result <- matrix(0,bas$Ai_max,3)
    d=matrix(delta,bas$Ai_max,3)
    delta_pointwise <- synthX(bas$Ylm,d)
    incomplete_integrand <- n * dot(t(n), t(delta_pointwise)) * sin(grd$U)
    for( i in 1:bas$Ai_max){
      Yi <- bas$Ylm[,i]
      integrand <- incomplete_integrand * Yi
      for (k in 1:3){
        result[i,k] <- int2d_scalar_GLS(integrand[,k], grd)
      }
    }
    return(result[])
  }


#' @export
Filter_Ortho<-function(A,bas)
{k=0
A1=bas$A;A1[]=0
for (l in bas$Lset)
  for (m in bas$LM[bas$LM[,1]==l,2])
  {k=k+1;  w.dim=which.max(abs(A[k,]))
  cat(l," ",m," ",w.dim,"\n")
  A1[k,w.dim]=A[k,w.dim]
  }
attr(A1,"orthified")=TRUE
attr(A1,"hash")=hash(A1)
return(A1)
}



#' @export
  E_FullModel_Penalty_AV<-function(A,grd,bas,Ref)
{
  updateX(A,grd,bas)->C
  h2<-E_SCM(A,grd,bas,C) # Wb contains full SCM energy, but not + K_b/2* M.C0^2 * Area
  S<-SEN(A,grd,bas,Ref,h2)
  e<-E_SEN(A,grd,bas,S,Ref)
  E<-h2$Wb + e + M.rho*((h2$Volume - bas$Target["Volume"])^2 + (h2$Area - bas$Target["Area"])^2) #+ K_b/2*C0^2*140
  names(E)=NULL # otherwise Volume is taken as name
  return(list(E=E, Wb=h2$Wb, Ws=e, E_uncons=h2$Wb + e, dA=h2$dA, S=S, Area=h2$Area, Volume=h2$Volume, Curv=h2$Curv))
}

#' @export
  Grad_FullModel_Penalty_AV<-function(A,grd,bas,Ref,S)
{
  updateX(A,grd,bas)->C
  h2<-E_SCM(A,grd,bas,C)
  Grad_SCM(h2,grd,bas,C)->G_SCM
  # gradH2=grad_SCM=  M.K_b/2 * (gradH2BC  - 2*M.C0*gradC ) +
  #  + M.K_ADE * (2 * H2$Curv * gradC / H2$Area - gradA * H2$Curv^2 / H2$Area^2 )
  Grad_SEN(A,grd,bas,G_SCM,S,Ref)->G_SEN
  G <- G_SCM$grad_SCM + G_SEN$grad_SEN + 2*M.rho*( G_SCM$gradV*(h2$Volume-bas$Target["Volume"]) +
                                                  + G_SCM$gradA*(h2$Area-bas$Target["Area"]))

#  W=h2$Wb + e + M.rho*((h2$Volume - bas$Target["Volume"])^2 + (h2$Area - bas$Target["Area"])^2) #+ K_b/2*C0^2*140
#  str(G)
  return(G)
}

  #' @export
  E_FullModel_Penalty_AVC<-function(A,grd,bas,Ref)
  {
    updateX(A,grd,bas)->C
    h2<-E_SCM(A,grd,bas,C) # Wb contains full SCM energy, but not + K_b/2* M.C0^2 * Area
    S<-SEN(A,grd,bas,Ref,h2)
    e<-E_SEN(A,grd,bas,S,Ref)
  #  print(names(bas$Target))
    E<-h2$Wb + e + M.rho*((h2$Volume - bas$Target["Volume"])^2 + (h2$Area - bas$Target["Area"])^2 + (h2$Curv - bas$Target["Curv"])^2) #+ K_b/2*C0^2*140
    names(E)=NULL # otherwise Volume is taken as name
    return(list(E=E, Wb=h2$Wb, Ws=e, E_uncons=h2$Wb + e, dA=h2$dA, S=S, Area=h2$Area, Volume=h2$Volume, Curv=h2$Curv))
  }

  #' @export
  Grad_FullModel_Penalty_AVC<-function(A,grd,bas,Ref,S)
  {
    updateX(A,grd,bas)->C
    h2<-E_SCM(A,grd,bas,C)
    Grad_SCM(h2,grd,bas,C)->G_SCM
    # gradH2=grad_SCM=  M.K_b/2 * (gradH2BC  - 2*M.C0*gradC ) +
    #  + M.K_ADE * (2 * H2$Curv * gradC / H2$Area - gradA * H2$Curv^2 / H2$Area^2 )
    Grad_SEN(A,grd,bas,G_SCM,S,Ref)->G_SEN
    G <- G_SCM$grad_SCM + G_SEN$grad_SEN + 2*M.rho*( G_SCM$gradV*(h2$Volume-bas$Target["Volume"]) +
                                                       + G_SCM$gradA*(h2$Area-bas$Target["Area"])+
                                                       + G_SCM$gradC*(h2$Curv-bas$Target["Curv"]))

    #  W=h2$Wb + e + M.rho*((h2$Volume - bas$Target["Volume"])^2 + (h2$Area - bas$Target["Area"])^2) #+ K_b/2*C0^2*140
    #  str(G)
    return(G)
  }


#' @export
Hessian_FullModel <- function(A,grd,bas,Ref,del,ncores)
{H=Hessian_SCM_SEN_cxx(A,grd,bas,Ref,del,ncores,M.C0,M.K_b,M.K_ADE)
 H$H=(H$H+t(H$H))/2
  return(list(H=H$H,gradA=H$gradA,gradV=H$gradV,gradC=H$gradC,g2=H$g2,h2=H$h2,G=H$G,comment="Hessian_SCM_SEN_cxx"))
}

#' @export
Spicules_as_Minima<-function(Q,grd,sgn=1)
{ iy=function(x){if (x<1) return(x+grd$nu) else return(x)}
jy=function(x){if (x>grd$nu) return(x-grd$nu) else return(x)}

ix=function(x){if (x<1) return(x+grd$nv) else return(x)}
jx=function(x){if (x>grd$nv) return(x-grd$nv) else return(x)}

idxx=function(p,grd){
  return(apply(p,1,function(x) x[1]+grd$nv*(x[2]-1)))
}
idx=function(x,grd){  return( x[1]+ grd$nv*(x[2]-1)) }

E=t(matrix(sgn*Q,grd$nu,grd$nv))
image(E)
p=matrix(0,grd$ndof,2);k=1
M=E;M[]=0# mark pixels in M by looking into 8 neigbours E-values
for (x in 1:(grd$nv))
  for (y in 1:(grd$nu))
  {        M[x,y]=
    as.integer(E[x,y]<=E[jx(x+1), y])+
    as.integer(E[x,y]<=E[x, jy(y+1)])+
    as.integer(E[x,y]<=E[x,               iy(y-1)]) +
    as.integer(E[x,y]<=E[ix(x-1),             y]) +
    as.integer(E[x,y]<=E[jx(x+1), iy(y-1)]) +
    as.integer(E[x,y]<=E[jx(x+1), jy(y+1)]) +
    as.integer(E[x,y]<=E[ix(x-1),             jy(y+1)]) +
    as.integer(E[x,y]<=E[ix(x-1),             iy(y-1)])
  if (M[x,y]==8) {p[k,]=c(x,y);k=k+1}
  }
image(M)
N=matrix((M>7),ncol=ncol(M))

p=p[1:(k-1),]
L=idxx(p,grd)
EE=E[L]

locate<-function(EE,E) {w=c();for (i in seq_along(EE)) w=c(w,which(c(E)==EE[i]));return(w)  }
EE=sort(EE)
I=locate(EE,t(E)) # back to original (i,j) by t()
I1=locate(EE,E)
E[I1]=NA
image(E)

plot(sort(EE))
image(t(E))
o=order(EE)
if(sgn==-1) o=rev(o)
plot(EE[o])
DE=c(0)
for (i in 2:(length(EE)-1)) DE[i]=EE[i-1]+EE[i+1]-2*EE[i]
plot(DE)
nspic=which.min(DE)-1
cat(crayon::red("Nspic ",nspic,"\n"))
return(list(M=M,E=E,p=p,L=L,o=o,Eo=EE[o],nspic=nspic))
}

#' @export
history_MemRBC<-function(M)
{
  print(paste(M$history))
}

#' @export
SetFlt_L1 <- function(M)
{ # build max-abs filter for L=1
  mx=c(which.max(abs(M$A[1,])), which.max(abs(M$A[2,])), which.max(abs(M$A[3,])))
  flt=c( (1:3)[-mx[1]], (1:3)[-mx[2]], (1:3)[-mx[3]])
  M$bas$flt=flt
  return(M)
}

#' @export
ApplyFlt_L1<-function(A,bas)
{ if (is.integer(bas$flt)) {
  A[1,bas$flt[1:2]]<-0
  A[2,bas$flt[3:4]]<-0
  A[3,bas$flt[5:6]]<-0
}
  else stop("Filter not set for ApplyFlt_L1: use SetFlt_L1(M) before")
  return(A)
}
