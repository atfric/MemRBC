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
load_param_RBC<-function()
{  data(M.mu,M.C0,M.mu,M.Ka,M.K_b,M.K_ADE,M.Es,M.rho,M.a2,M.a3,M.a4,M.b0,M.b1,M.b2,M.rho,M.Rcpp,M.Rcpp_ncores,M.scr1,M.scr2,M.muk,M.lam,M.cl,package = "MemRBC",envir = MemRBC_env) 
}

#' @export
print_param_RBC<-function()
{ sapply(names(MemRBC_env),function(x) {cat(x,MemRBC_env[[x]],"\n");return(NULL)}) ; return()}

#
# MemRBC_env$M.xxx parameters to be modified in main program
#

#' @export
citation.MemRBC<-function() {cat("when using this software for publications you must cite it as:\n Frickenhaus, S. (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340 \n")}


#' @export
MemRBC_env <- new.env(parent = emptyenv())
MemRBC_env$M.C0=0
MemRBC_env$M.mu=2.5
MemRBC_env$M.Ka=5
MemRBC_env$M.K_ADE=0.2
MemRBC_env$M.a3=-2
MemRBC_env$M.a2=1
MemRBC_env$M.b0=1
MemRBC_env$M.a4=8
MemRBC_env$M.b1=0.7
MemRBC_env$M.b2=0.75
MemRBC_env$M.cl=NULL
MemRBC_env$M.Es=112
MemRBC_env$M.K_b=0.2
MemRBC_env$M.rho=50

MemRBC_env$M.lam=c(1,1,1)
MemRBC_env$M.muk=10
MemRBC_env$M.scr1=-1
MemRBC_env$M.scr2=-1
MemRBC_env$M.Rcpp=TRUE
MemRBC_env$M.Rcpp_ncores=4
MemRBC_env$M.


.onLoad <- function(libname, pkgname) {
  # startup message
  msg="MemRBC, the red blood cell shape modeling R package... \n... with compiled SCM energy, gradient and Hessian by Rcpp.\n Needs time for compiling...\n use get_data_ZENODO() to load datasets, eg. for data(ss42denovo)\n"

  if(!interactive()) {
    msg[1] <- paste("Package 'MemRBC'")
  }

  packageStartupMessage(msg)
  # to have cached objects available on parallel cluster, 
  # place it in tempdir/../cachedir, not in tempdir, which may differ in seperate cluster R.
  { 
    if (base::.Platform$OS.type=="windows") {
      strsplit(tempdir(),"\\\\")[[1]]->s;
      paste(c(s[1:(length(s)-1)],"MemRBC_cache"),sep="",collapse="\\")->cdir;
      } else cdir="/tmp/MemRBC_cache"; # should work for MacOS and Linux
    if (!dir.exists(cdir)) {dir.create(cdir);}
    MemRBC_env$cache<-cdir;
    Rcpp::sourceCpp( cacheDir = cdir, rebuild = FALSE,
      code='\\   
#include <Rcpp.h>
#include<RcppEigen.h>
using namespace std;
using namespace Rcpp;
#include <omp.h>

// [[Rcpp::plugins(RcppEigen]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]

// integration with Rcpp, see RcppCore,Bates/Eddelbuettel, https://github.com/RcppCore/RcppEigen/blob/master/README.md
using Eigen::Map;         // maps rather than copies
using Eigen::MatrixXd;    // variable size matrix, double precision
using Eigen::VectorXd;    // variable size vector, double precision

// for a faster version of cotan Laplacian GEMINI-Pro was asked
// to derive a R to C++ translation for Rcpp

//GEMINI
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>

using namespace Eigen;

// Typedefs not possible in src/*.cpp 
// Typedef for the sparse matrix (double precision)
typedef SparseMatrix<double> SpMat;
typedef Triplet<double> T;

// [[Rcpp::export(.GEMINI_get_cotan_Laplacian_cxx)]]
SpMat GEMINI_get_cotan_Laplacian_cxx(List mesh) {

  // 1. Extract data from R rgl mesh object
  // mesh$vb is a 4xN matrix (homogeneous coords), we need the top 3 rows
  MatrixXd V_full = as<MatrixXd>(mesh["vb"]);
  int n_verts = V_full.cols();

  // mesh$it is a 3xM matrix of triangle indices (1-based)
  MatrixXi F = as<MatrixXi>(mesh["it"]);
  int n_faces = F.cols();

  // 2. Prepare triplets for Sparse Matrix
  // Each face has 3 edges. Each edge affects 4 entries in the matrix
  // (2 diagonal, 2 off-diagonal) per face.
  // Reserve memory to prevent re-allocation.
  std::vector<T> tripletList;
  tripletList.reserve(n_faces * 12);

  // 3. Loop over faces
  for(int i = 0; i < n_faces; ++i) {
    // Get indices (convert 1-based R to 0-based C++)
    int idx0 = F(0, i) - 1;
    int idx1 = F(1, i) - 1;
    int idx2 = F(2, i) - 1;

    // Get vertices
    Vector3d v0 = V_full.block<3,1>(0, idx0);
    Vector3d v1 = V_full.block<3,1>(0, idx1);
    Vector3d v2 = V_full.block<3,1>(0, idx2);

    // Compute edge vectors
    Vector3d e0 = v1 - v0;
    Vector3d e1 = v2 - v1;
    Vector3d e2 = v0 - v2;

    // Compute Cotangents of the angles opposite to edges
    // Cot(angle at v0) uses edges connected to v0: (v1-v0) and (v2-v0)
    // Careful with directions for dot product: u.v / |uxv|
    // We can use the squared edge lengths or cross products.
    // Explicit vector formulation is safest:

    // Angle at 0 (between v1-v0 and v2-v0) -> Opposite edge e1 (v1-v2)
    Vector3d u0 = v1 - v0;
    Vector3d v0_vec = v2 - v0;
    double cot0 = u0.dot(v0_vec) / u0.cross(v0_vec).norm();

    // Angle at 1 (between v2-v1 and v0-v1) -> Opposite edge e2 (v2-v0)
    Vector3d u1 = v2 - v1;
    Vector3d v1_vec = v0 - v1;
    double cot1 = u1.dot(v1_vec) / u1.cross(v1_vec).norm();

    // Angle at 2 (between v0-v2 and v1-v2) -> Opposite edge e0 (v0-v1)
    Vector3d u2 = v0 - v2;
    Vector3d v2_vec = v1 - v2;
    double cot2 = u2.dot(v2_vec) / u2.cross(v2_vec).norm();

    // Handle degenerate triangles (area ~ 0 -> inf cotan)
    if (!std::isfinite(cot0)) cot0 = 0.0;
    if (!std::isfinite(cot1)) cot1 = 0.0;
    if (!std::isfinite(cot2)) cot2 = 0.0;

    // 4. Push triplets
    // Weight w = 0.5 * cot(angle)
    // The Laplacian L = D - W.
    // For an edge connecting i,j with weight w:
    // L(i,j) -= w
    // L(j,i) -= w
    // L(i,i) += w
    // L(j,j) += w

    double w0 = 0.5 * cot0; // Associated with Edge opposite vertex 0 (between 1 and 2)
    double w1 = 0.5 * cot1; // Associated with Edge opposite vertex 1 (between 0 and 2)
    double w2 = 0.5 * cot2; // Associated with Edge opposite vertex 2 (between 0 and 1)

    // Edge (1,2) gets w0
    tripletList.push_back(T(idx1, idx2, -w0));
    tripletList.push_back(T(idx2, idx1, -w0));
    tripletList.push_back(T(idx1, idx1, w0));
    tripletList.push_back(T(idx2, idx2, w0));

    // Edge (0,2) gets w1
    tripletList.push_back(T(idx0, idx2, -w1));
    tripletList.push_back(T(idx2, idx0, -w1));
    tripletList.push_back(T(idx0, idx0, w1));
    tripletList.push_back(T(idx2, idx2, w1));

    // Edge (0,1) gets w2
    tripletList.push_back(T(idx0, idx1, -w2));
    tripletList.push_back(T(idx1, idx0, -w2));
    tripletList.push_back(T(idx0, idx0, w2));
    tripletList.push_back(T(idx1, idx1, w2));
  }

  // 5. Construct sparse matrix
  SpMat L(n_verts, n_verts);
  L.setFromTriplets(tripletList.begin(), tripletList.end());

  return L;
}
//GEMINI END


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

// [[Rcpp::export]]
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
// [[Rcpp::export]]
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


// [[Rcpp::export]]
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


// [[Rcpp::export]]
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


// [[Rcpp::export]]
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


// [[Rcpp::export]]
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
// [[Rcpp::export(.dot2)]]
double dot2(NumericVector x, NumericVector y) {
  return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}

// [[Rcpp::export(.IntegM)]]
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

// [[Rcpp::export(.IntegS)]]
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


// [[Rcpp::export]]
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

 ')
  } # sourceCPP
  
  
citation.MemRBC();

utils::data(M.mu,M.C0,M.mu,M.Ka,M.K_b,M.K_ADE,M.Es,M.rho,M.a2,M.a3,M.a4,M.b0,M.b1,M.b2,M.rho,M.Rcpp,M.Rcpp_ncores,M.cl,M.muk,M.lam,M.scr1,M.scr2,cache,package = "MemRBC",envir = MemRBC_env)
}



HAVE_DEPRECATED=FALSE

M.TEST=FALSE
# main code must define M.TEST before sourcing

# for tests: show severe deviation from zero as error
#' severe
#' @description
#' test two quantities for equality using test_that
#' 
#' @export
severe<-function(q1,q2, what="some test", tol=1e-12)
{ testthat::test_that(what,{testthat::expect_equal(q1,q2,tolerance=tol)})}
#{cat(q,":");if (all(abs(q)>tol)) {stop("severe imprecision, STOP")} else
#  cat(crayon::green("OK:",what,"\n"))
#}


# @useDynLib MemRBC, .registration = TRUE
# @importFrom Rcpp sourceCpp
# NULL

#
# not working when Rcpp-code is loaded
#cleanCache<-function(){cat("clean up ",MemRBC_env$cache,"\n");file.remove(MemRBC_env$cache)}


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

#' E_SCM
#' Bending energy , returns list with "Wb" as bending energy
#' @description
#' compute bending energy
#' @param A,grd,bas,C : coefficients, grid, basis and coordinates
#' @param plt (=FALSE) for plotting 3d
#' @param clp (=FALSE) for plotting with a central clipping-plane
#' @export
E_SCM <- function (A, grd, bas, C, plt = FALSE, dbg = FALSE, clp = FALSE)
{
  if (!MemRBC_env$M.Rcpp)
    return(
      structure(class="MemESCM",
      E_SCM_R(A, grd, bas, C, plt = FALSE, dbg = FALSE,
                   clp = FALSE)) 
       )
  return(structure(class="MemESCM",
    E_SCM_cxx(A, grd, bas, C, MemRBC_env$M.C0, MemRBC_env$M.K_b, MemRBC_env$M.K_ADE))
    )
}


#' Bending energy in R code, returns list with "Wb" as bending energy
#' @description
#' compute bending energy, no C-Code
#' @param A,grd,bas,C : coefficients, grid, basis and coordinates
#' @param plt (=FALSE) for plotting 3d
#' @param clp (=FALSE) for plotting with a central clipping-plane
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
  # normalization of normals / some may be 0 (at poles, where dv=0)
  Nn <- apply(normal,1,pracma::Norm)
  if (sum(Nn==0)>0) stop("ERROR: zeros in ||normal||")
  iNn<-1/Nn
  # normalize to unit normals
  for (k in 1:3) n[,k]<-n[,k]*iNn # second fundamental needs unit n
  # fundamental form I
  E  <- Xu*Xu + Yu*Yu + Zu*Zu
  FF <- Xu*Xv + Yu*Yv + Zu*Zv
  G  <- Xv*Xv + Yv*Yv + Zv*Zv
  dA <- sqrt(E*G-FF*FF) # eqals Nn, this is sqrt(det(g)) with metric tensor g
  (Area=int2d_s(dA,grd))
  inn<-1/dA

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
  H2 <- MemRBC_env$M.K_b/2 * (H2_BC - 2*MemRBC_env$M.C0*Curv) + # + K_b/2*C0^2*Area + # constant terms out
    + MemRBC_env$M.K_ADE*(Curv^2 ) /Area

  E_SCM_dens<-MemRBC_env$M.K_b/2 * (curv_sq - 2*MemRBC_env$M.C0*curv) + # + K_b/2*C0^2 + # constant terms out
    + MemRBC_env$M.K_ADE*(curv*Curv ) /Area

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

#' Grad_SCM
#' @description
#' compute bending energy gradient
#' @param h2 : result from E_SCM()
#' @param grd,bas,C : coefficients, grid, basis and coordinates
#' @param int2d_m : integration function, e.g. .IntegM (hidden c++ R-function)
#' @export
Grad_SCM <- function(h2,grd,bas,C,int2d_m=.IntegM)
{ if (!MemRBC_env$M.Rcpp) return(Grad_SCM_R(h2,grd,bas,C))
  G2<-Grad_SCM_cxx(h2,grd, bas, C, MemRBC_env$M.C0, MemRBC_env$M.Rcpp_ncores, MemRBC_env$M.K_b, MemRBC_env$M.K_ADE)
  return(list(ddA=array(G2$ddA,c(grd$ndof,bas$Ai_max,3)),
              ddV=array(G2$ddV,c(grd$ndof,bas$Ai_max,3)),
              grad_SCM=G2$grad_SCM,gradV=G2$gradV,gradA=G2$gradA,gradC=G2$gradC,
              dE=array(G2$dE,c(grd$ndof,bas$Ai_max,3)),
              dF=array(G2$dF,c(grd$ndof,bas$Ai_max,3)),
              dG=array(G2$dG,c(grd$ndof,bas$Ai_max,3))))
}

#' Grad_SCM_R 
#' R version of SCM energy gradient, vectorizes first dimension (i = spatial)
#' @description
#' compute bending energy gradient in R (slow)
#' @param Wb : result from E_SCM
#' @param grd,bas,C :  grid, basis and coordinates
#' @param plt (=FALSE) for plotting 3d
#' @param clp (=FALSE) for plotting with a central clipping-plabe
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
  for(k in 1:3) for (j in 1:Ai_max) ddV[,j,k] <- (H2$normal[,k]*bas$Ylm[,j]) # normal = n*dA

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
  gradH2  <-  MemRBC_env$M.K_b/2 * (gradH2BC  - 2*MemRBC_env$M.C0*gradC ) +
    + MemRBC_env$M.K_ADE * (2 * H2$Curv * gradC / H2$Area - gradA * H2$Curv^2 / H2$Area^2 )

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
#' gradient of area and volume terms for constraint Jacobian
#' @description
#' compute area and volume gradient in R (slow)
#' @param Wb : result from E_SCM
#' @param grd,bas,C :  grid, basis and coordinates
#' @export
Grad_SCM_av <- function(Wb,grd,bas,C,int2d_m=.IntegM)
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
#' from SCM energy compute principal curvatures
#' @description
#' compute principle curvatures k1,k2 from precomputed E_SCM
#' @param Wb : result from E_SCM
#' @param grd :  grid for integration
#' @param plt.K (=FALSE) for plotting 2d
#' @return list of curvature k1, k2, computed from K=L*N-M^2 (second fundamental)
#' @export
Membrane_Curvatures<-function(Wb, grd, plt.K=FALSE)
{
  int2d_s=.IntegS
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
# and none exported

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

#' MakeRef
#' @description
#' Build reference from input coefficients Aref and store in M.
#' If missing and demanded, M$Ref and M$SEN are computed here.
#' This helps plotting SEN colors like in PNEM and other Apps.
#' usually Aref comes from an slightly oblate shape of 90-95% of the target volume 100.
#' @param M MemRBC input object
#' @param Aref Reference shape coefficients for SEN
#' @return MemRBC object with Ref and SEN set.
#' @export
MakeRef<-function(M,Aref){
  M$Ref=Ref4CauchyGreen(Aref,M$grd,M$bas)
  h2=E_SCM(Aref,M$grd,M$bas,updateX(Aref,M$grd,M$bas))
  M$SEN=SEN(Aref,M$grd,M$bas,M$Ref,h2)
}

#' Ref4CauchyGreen
#' @description
#' compute lists of reference tensors for fast SEN computation
#' @param A, grd, bas: coefficients, grid and basis
#' @return Reference object containing pieaces of Cauchy-Green-tensors
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
return(structure(class="MemRef",list( tgi=tgi, giPrep=giPrep, h2ref=h2ref, v=h2ref$Volume ,a=h2ref$Area ,c=h2ref$Curv,ARef=A )))
}

#' SEN
#' compute stretch parameters (alpha, beta) from SEN
#' @description
#' Compute the stretches alpha and shear beta from SEN parameters.
#' Needed before E_SEN is called.
#' @param A,grd,bas,Ref Coefficients, grid, basis and reference objects
#' @param Wb_cur : current E_SCM() result
#' @return SEN object with alpha, beta and mean trace of right Cauchy-Green
#' @export 
SEN<-function (A, grd, bas, Ref, Wb_cur)
{
  h2cur <- Wb_cur
  G <- vmat4_2lmat(h2cur$E, h2cur$FF, h2cur$FF, h2cur$G)
  RightCauchyGreen <- vmatmat_2lmat(Ref$giPrep, G)
  m <- 0.5 * sapply(RightCauchyGreen, function(x) c(x[1, 1] +
                                                      x[2, 2]))
  alpha <- h2cur$dA/Ref$h2ref$dA - 1
  beta <- m/(alpha + 1) - 1
  return(structure(class="MemSEN",list(alpha = alpha, beta = beta, m = m, h2cur = h2cur,
              A = A, G = G)))
}

#' E_SEN
#' energy of Shear-Elastic Network (SEN)
#' @description
#' Compute the SEN energy
#' @param A,grd,bas,S, Ref Coefficients, grid, basis, SEN and reference objects
#' @param Wb_cur : current E_SCM() result
#' @return list with alpha and beta as well as trace m of right Cauchy-Green
#' @export
E_SEN<-function(A,grd,bas,S,Ref)
{ if (!is.null(Ref)){
  WS<- MemRBC_env$M.Ka/2 * .IntegS((MemRBC_env$M.a2*S$alpha^2+MemRBC_env$M.a3*S$alpha^3+MemRBC_env$M.a4*S$alpha^4)*Ref$h2ref$dA,grd) +
    + MemRBC_env$M.mu*.IntegS( ( (MemRBC_env$M.b0+MemRBC_env$M.b1*S$alpha)*S$beta + MemRBC_env$M.b2*S$beta^2)*Ref$h2ref$dA, grd)
  return(WS)} else return(0)
}

#' Gradient of SEN energy
#' @description
#' Compute gradients of SEN energy
#' @param A,grd,bas Coefficients, grid, basis 
#' @param h2cur_grad current SCM gradient (Grad_SCM())
#' @param S : current SEN() result
#' @param Ref : static reference of SEN
#' @return gradS for SEN energy gradient, gradAlpha, gradBeta, dm for gradients of stretches and shear.
#' @export
Grad_SEN<-function(A, grd, bas, h2cur_grad, S, Ref, int2d_m = .IntegM){
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
  gradS <- (MemRBC_env$M.Ka*0.5 * int2d_m( Ref$h2ref$dA *
                              gradAlpha*(MemRBC_env$M.a2*2*S$alpha + 3*MemRBC_env$M.a3*S$alpha^2 +
                                                   + 4*MemRBC_env$M.a4*S$alpha^3),grd,bas)  +
    +   MemRBC_env$M.mu * int2d_m( Ref$h2ref$dA *
                          ( gradBeta*( MemRBC_env$M.b0 + MemRBC_env$M.b1*S$alpha + 2*MemRBC_env$M.b2*S$beta) +
                                          + MemRBC_env$M.b1*gradAlpha*S$beta), grd, bas)
  )
  return(list(grad_SEN=gradS, gradAlpha=gradAlpha,
              gradBeta=gradBeta, dm=dm ))
}


#
# utilities
#


#' annotates A with L-M-strings as rownames, + some attributes
#' @export
LM2A<-function(A,bas)
{
nm=apply(bas$LM,1,paste,sep=";",collapse=";")
if (!is.matrix(A)) A=matrix(A,ncol=3)
rownames(A)<-nm;
colnames(A)<-LETTERS[24:26];
attr(A,"C0")<-MemRBC_env$M.C0
attr(A,"V0")<-bas$Target["Volume"]
attr(A,"A0")<-bas$Target["Area"]
if (bas$Nc>2) attr(A,"Ct")<-bas$Target["Curv"]

attr(A,"Target")<-bas$Target
return(A)
}

#' save A as object named Alm with additional attributes
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
  load(file)->n
  return(get(n))
}

#
# Total energy from SCM and SEN
#  SCM stored in h2cur in stretches S
#
#' TotalEnergyDensity
#' @description
#' compute total energy density, to visualize 
#' 
#' @param S SEN, containing also h2cur$E_SCM_dens
#' @return energy value
#' @export
TotalEnergyDensity<-function(S)
{
  return( S$h2cur$E_SCM_dens +
            MemRBC_env$M.Ka/2 * (S$alpha^2+MemRBC_env$M.a3*S$alpha^3+MemRBC_env$M.a4*S$alpha^4) + #*Ref$h2ref$dA +
            + MemRBC_env$M.mu* ( (1+MemRBC_env$M.b1*S$alpha)*S$beta + MemRBC_env$M.b2*S$beta^2) # *Ref$h2ref$dA
  )
}

# not used
int2d_matrix_3<-function (field_m,grd,bas=NULL)
{ dms=dim(field_m); res=array(0.0,dms[2:4]);
  for (i in 1:dms[2]) for (j in 1:dms[3]) for (k in 1:dms[4]) res[i,j,k]=int2d_scalar_GLS(field_m[,i,j,k],grd)
  return(res)
}


# plot coeffs
#' @export
plotA<-function(A,...){matplot(A,type="l",lty=1,xlab=expression((l+1)^2-l+m),...)}

#' plotA_l
#' plot coeffs
#'   and l-values as color bar on x-axis
#' @param A,bas : standard coefficient amd basis functions objects
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

#' helper routines for 3d graphics objects
#' @export
Obj2X<-function(O) # extract vertex coordinates
{ return(t(O$vb[1:3,]))
}

#' make coordinates from centered object
#' @export
Obj2X_centre<-function(O)
{ X=t(O$vb[1:3,]) # extract vertex coordinates and centre
for (k in 1:3 ) X[,k]=X[,k]-mean(X[,k])
return(X)
}

#' make coordinates from object
#' @export
X2Obj<-function(O,X)
{ O$vb=rbind(t(X),1) # ingest coordinates in 3d-graphics object
return(Rvcg::vcgClean(Rvcg::vcgUpdateNormals(O),silent=TRUE))
# clean helps smooth color at zero meridian
}

#'  build vertex areas from triangle areas using Rvcg
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

#' synthesize X coordinates by matrix-matrix-product
#' @export
synthX<-function(Y,A) # Y are precomputed spherical harmonics
{ return(Y%*%A)
}

# the following imag.obj... are showing data as color code on the 3d object
#

#' plot on 3D object surface a field f as color code
#' @export
imag.obj.colorbar.simple<-function(obj,f,clr=TRUE,...) {
  if(clr) rgl::clear3d()
  cols=rainbow(100);
  rgl::shade3d(obj,meshcolor="vertices",color=cols[(f-min(f))/diff(range(f))*99+1],...)
  rgl::bgplot3d(own.imagePlot(legend.only = TRUE, zlim = range(f), col = cols) )
}


#' imag.obj.colorbar
#' allow for limits to suppress outliers (color black)
#' @export
imag.obj.colorbar<-function(obj,f,limits=range(f),clr=FALSE,pal=heat.colors,width=550,height=480,par=FALSE,...) {
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
  rgl::bgplot3d(own.imagePlot(legend.only = TRUE,add=TRUE,zlim = limits, col = cols) )
}

#' MakeGrid_GaussLegendre
#' this has more dense points next to 0 meridian
#'  alternative with regular v-spacing is GaussLegendreSimpson
#' @param n : number of points along each dimension
#' @param uv_fac (=1) : for experiemnts with half periods in v set uv_fac=0.5
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

  # for vectorization we also have u and v as vectors U,V
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
  return(structure(class="MemGrd",grd))
}


#' display f data on object obj
#' @export
imag.obj<-function(obj,f,pal=rainbow) {
  cols=pal(100); rgl::shade3d(obj,meshcolor="vertices",col=cols[as.integer(1+99*(f-min(f))/diff(range(f)))]) }

#' imag
#' @description
#' plot a scalar in 2D
#' @param field scalar spatial field
#' @param grd the grid of (u,v) on which the scalar is defined
#' @param nx,ny dimension of output matrix plot
#' @export
imag<-function (field, grd, nx = grd$nv, ny = grd$nu/2, ...)
{
  fields::quilt.plot(grd$V, grd$U, field[], nx, ny, xlab = "v",
                     ylab = "u", ...)
}

#' MakeBasis_UV
#'
#' compute the basis functions on given (u,v) values
#' (u,v) are usually from grd$U, grd$V, but you may also use irregular (u,v)
#' for example for fitting coefficients to a 3d-object, in which (u,v)-values
#' may be given as texture coordinates texccords.
#'
#'  @param L_max (=4) spectral order of basis
#'  @param u,v vectors of spehrical angles u (in 0..pi) and v (in 0..2*pi).
#'  @param only_Ylm (=FALSE) for excluding derivative computation.
#'  @return basis, to be stored as $bas in membrane object or for other use.
#'  The basis also contains standard constraints for area and volume, which may be modified with SetConstraints()
#' @export
MakeBasis_UV<-function (L_max = 4, u, v, Pointsymmetry = FALSE,
                        only_Ylm=FALSE, kind="Ylm", KleinBottle=FALSE)
{ if (kind=="Ylm") {
  Ai_max = (L_max + 1)^2 - 1
  LM = data.frame(l = rep(0L, Ai_max), m = rep(0L, Ai_max))
  AI = 1
  for (l in 1:L_max) for (m in (-l):l) {
    LM[AI, ] = c(l, m)
    AI <- AI + 1
  }
  L = LM[, 1]
  M = LM[, 2]
  LM1 = LM
  if (Pointsymmetry) {
    w = which(L%%2 > 0)
    LM <- LM[w, ]
    L = LM[, 1]
    M = LM[, 2]
  }
  else w = (1:Ai_max)
  Ai_max = length(L)
  L_max = max(L)
  n.v = length(u)
  if (length(v) != n.v)
    stop("u not same length like v")
  L_Ylm_ = L_Ylm(L_max, u, v)
  Ylm_ = L_Ylm_$Ylm[, -1]/sqrt(4 * pi)
  if (!only_Ylm){
   Ylm_v_ = Ylm_v(L_max, u, v, L_Ylm_$PLK)[, -1]/sqrt(4 * pi)
   Ylm_vv_ = Ylm_vv(L_max, u, v, L_Ylm_$PLK)[, -1]/sqrt(4 * pi)
   L_Y_u_ = L_Ylm_u(L_max, u, v, L_Ylm_$PLK)
   Ylm_u_ = L_Y_u_$Ylm_u[, -1]/sqrt(4 * pi)
   Ylm_uu_ = Ylm_uu(L_max, u, v, L_Y_u_$P_T)[, -1]/sqrt(4 * pi)
   Ylm_uv_ = Ylm_uv(L_max, u, v, L_Ylm_$PLK, L_Y_u_$P_T)[, -1]/sqrt(4 * pi)
   }
  l = LM[, 1]
  m = LM[, 2]
  bas = list(n.v = n.v, uv = cbind(u, v), Ylm = Ylm_[, w], LM = LM, Ai_max = Ai_max, l = l,
             m = m, A = matrix(0, Ai_max, 3), L_max = L_max, G.tk = l^2 *
               (l + 1)^2, Wt = l * (l + 1), comment = "(for double entries masked) irregular or Gauss-Legendre-Simpson basis from (u,v), computed with W. Bosch/ K. Khairy codes excluding l=0, A/V constraints set",
             Nupd = 0, Lset = unique(l), Mset = unique(m), Nc = 2,
             Cons = c("gradA", "gradV"), QCons = c("Area", "Volume"),
             Target = c(140, 100), Pointsymmetry = Pointsymmetry)
  names(bas$Cons) = names(bas$QCons) = names(bas$Target) = c("Area","Volume")
  if (!only_Ylm){
    bas$Ylm_u = Ylm_u_[,w];
    bas$Ylm_v = Ylm_v_[,w];
    bas$Ylm_uu = Ylm_uu_[,w];
    bas$Ylm_uv = Ylm_uv_[,w];
    bas$Ylm_vv = Ylm_vv_[,w];
  }
  bas$A = LM2A(bas$A, bas)
  mask = double_uv_ind(u, v)
  bas$mask <- ifelse(is.numeric(mask), mask, 1) # minimum mask needed
  bas$kind=kind
  return(bas) }
  if (kind=="Fourier") {
    Ai_max=L_max^2*4+1
    n.v = length(u)
    if (n.v!=length(v)) stop("u and v have different length in MakeBasis_UV")
    Ylm=matrix(0.0,length(u),Ai_max)
    if (!only_Ylm) Ylm_u=Ylm_v=Ylm_uu=Ylm_vv=Ylm_uv=matrix(0.0,length(u),Ai_max)
    LM=matrix(NA,Ai_max,2)
    colnames(LM)=c("L","M")
    if (KleinBottle) uspace= (1:L)/2 else uspace= 1:L_max
    cat("U factors:",uspace,"\n")
    Ylm[,1]=1
    if (!only_Ylm) {Ylm_u_[,1]=Ylm_v_[,1]=Ylm_vv_[,1]=Ylm_uu_[,1]=Ylm_uv_[,1]=0}
    k=2
    for (i in uspace)
     for (j in 1:L_max){
      Ylm_[,k]  = sin(i*u)*sin(j*v)
    Ylm_[,k+1] = sin(i*u)*cos(j*v)
    Ylm_[,k+2] = cos(i*u)*cos(j*v)
    Ylm_[,k+3] = cos(i*u)*sin(j*v)

    if (!only_Ylm) {

    Ylm_u_[,k]  =  i*cos(i*u)*sin(j*v)
    Ylm_u_[,k+1]=  i*cos(i*u)*cos(j*v)
    Ylm_u_[,k+2]= -i*sin(i*u)*cos(j*v)
    Ylm_u_[,k+3]= -i*sin(i*u)*sin(j*v)

    Ylm_v_[,k]  =  j*sin(i*u)*cos(j*v)
    Ylm_v_[,k+1]= -j*sin(i*u)*sin(j*v)
    Ylm_v_[,k+2]= -j*cos(i*u)*sin(j*v)
    Ylm_v_[,k+3]=  j*cos(i*u)*cos(j*v)

    Ylm_uu_[,k]  =  i^2*Ylm_[,k]
    Ylm_uu_[,k+1]=  i^2*Ylm_[,k+1]
    Ylm_uu_[,k+2]=  i^2*Ylm_[,k+2]
    Ylm_uu_[,k+3]=  i^2*Ylm_[,k+3]

    Ylm_vv_[,k]  =  j^2*Ylm_[,k]
    Ylm_vv_[,k+1]=  j^2*Ylm_[,k+1]
    Ylm_vv_[,k+2]=  j^2*Ylm_[,k+2]
    Ylm_vv_[,k+3]=  j^2*Ylm_[,k+3]

    Ylm_uv_[,k]  = -j*i*cos(i*u)*cos(j*v)
    Ylm_uv_[,k+1]= -j*i*cos(i*u)*sin(j*v)
    Ylm_uv_[,k+2]= -i*j*sin(i*u)*sin(j*v)
    Ylm_uv_[,k+3]= -i*j*sin(i*u)*cos(j*v)
    }
     LM[k:(k+3),]=c(i,j)

     k=k+4

     }

    l=LM[,1]
    m=LM[,2]
  w=1:Ai_max # no symmetries
  mask = double_uv_ind(u, v)
  bas = list(n.v = n.v, uv = cbind(u, v), Ylm = Ylm_[, w], LM = LM, Ai_max = Ai_max, l = l,
             m = m,  L_max = L_max, G.tk = l^2*m^2, Wt = l*m, comment = "Fourier basis, no cos(0)",
             Nupd = 0, Lset = unique(l), Mset = unique(m), Nc = 2,
             Cons = c("gradA", "gradV"), QCons = c("Area", "Volume"),
             Target = c(140, 100),  
             Pointsymmetry = NA)
  if (!only_Ylm) {bas$Ylm_u=Ylm_u_[,w];bas$Ylm_v=Ylm_v_[,w];
                  bas$Ylm_uu=Ylm_uu_[,w];bas$Ylm_uv=Ylm_uv_[,w];
                  bas$Ylm_vv=Ylm_vv_[,w]}

  bas$LM=LM
  bas$uv=cbind(c(u),c(v))
  bas$mask <- ifelse(is.numeric(mask) == 0, mask, 1)
  bas$kind=kind
  bas$KleinB=KleinB
  bas$u=u;bas$v=v
  bas$A=matrix(0,Ai_max,3)
  bas$A=LM2A(bas$A,bas)
  return(structure(class="MemBas",bas))
  }
  warning ("Maybe wrong kind specified! return NULL.")
  return(NULL)
}


#' updateX
#'
#' updates coordinates C$X and their partial derivatives wrt. u,v, like C$X_u.
#' @param A coefficients of shape
#' @param grd grid from on which the basis is computed
#' @param bas basis function values $Ylm and their derivatives, like $Ylm_u.
#' @return Coord object of class MemC, with X and derivatives X_u, ... and input Coeff A
#' @export
updateX<-function (A, grd, bas)
{
  X <- bas$Ylm %*% A
  X_u <- bas$Ylm_u %*% A
  X_v <- bas$Ylm_v %*% A
  X_uu <- bas$Ylm_uu %*% A
  X_uv <- bas$Ylm_uv %*% A
  X_vv <- bas$Ylm_vv %*% A
  return(structure(class="MemC",list(X = X, X_u = X_u, X_v = X_v, X_uu = X_uu, X_vv = X_vv,
              X_uv = X_uv, Coeff = A)))
}


# compute coordinates only
# could be used instead of updateX in several places

#' updateX_only
#' updates coordinates C$X, NOT the derivatives.
#' @param A coefficients of shape
#' @param grd grid from on which the basis is computed
#' @param bas basis function values $Ylm
#' @return class MemC_X object, containing X and input Coeff A
#' @export
updateX_only<-function(A, grd, bas)
{
  bas$Ylm %*% A -> X
  return(structure(class="MemC_X",list(X=X, Coeff=A)))
}


#' MakeSphere
#' @description
#' Compute unit sphere for given grid and basis
#' 
#' @param grd given grid
#' @param bas given basis
#' @param r (=1) for radius of output sphere
#' @return coefficient matrix, derived from its prototype bas$A
#' @export
MakeSphere<-function(grd,bas,r=1)
{
  A=bas$A; A[,]=0
  A[1,"X"]=1/0.48860251190292*r
  A[2,"Z"]=1/0.48860251190292*r
  A[3,"Y"]=1/0.48860251190292*r
  return(structure(class="MemA",A))
}

# the central driver to 2D integration
# use faster .IntegS from cpp instead
int2d_scalar_GLS<-function(F2,grd)
{ Z <- matrix(F2,grd$nu,grd$nv)
  Q <- grd$wx %*% Z %*% as.matrix(grd$wy)
return(Q[,])
}

# plot X/object from segments in 3d
# with Cartesian Wireframe as (X,Y,Z)-contour levels for cont=TRUE
#' plot3a
#' @description
#' plot coordinates C$X (after updateX)
#' @param X Coordiates, e.g. in C$X from updateX()
#' @param grd grid with a basic rgl-object grd$Obj
#' @examples
#' data("M4",package = "MemRBC")  # take required data from M5
#' plot3a(updateX(M4$A,M4$grd,M4$bas)$X,M4$grd)
#' @export
plot3a<-function (X, grd, pnts = FALSE, clip = FALSE, col = "black",
                  alpha = 1, cont = TRUE, cont.grid = FALSE, fill = TRUE, fn = "z",
                  fn_data = "z", ...)
{
  O <- X2Obj(grd$Obj, X)
  if (cont.grid) {
    id <- rgl::contourLines3d(O, grd$U, col = col, levels = pracma::linspace(0,
                                                                             pi, 12))
    rgl::contourLines3d(O, grd$V, col = col, levels = pracma::linspace(0,
                                                                       2 * pi, 13))
  }
  if (cont) {
    id <- rgl::contourLines3d(O, X[, 1], col = 1)
    rgl::contourLines3d(O, X[, 2], col = 2)
    rgl::contourLines3d(O, X[, 3], col = 3)
  }
  if (fill) {
    if (fn == "z")
      id <- rgl::filledContour3d(O, fn = X[, 3], alpha = alpha)
    else id <- rgl::filledContour3d(O, fn = fn_data, alpha = alpha,
                                    ...)
  }
  if (pnts)
    id <- rgl::points3d(X[, 1], X[, 2], X[, 3], col = "red",
                        cex = 2)
  if (clip)
    rgl::clipplanes3d(c(0.5, 0.5, 0))
  invisible()
}


#' plot3q
#'
#' plot coordinates C$X (after updateX) with wireframe of quadrilaterals
#' @param X Coordiates, e.g. in C$X from updateX()
#' @param grd grid with a basic rgl-object grd$Obj
#' @examples
#' data("M4",package = "MemRBC")  # take required data from M4
#' plot3q(updateX(M4$A,M4$grd,M4$bas)$X,M4$grd)
#' @export
plot3q<-function (X, grd, col = "black", alpha = 1, ...)
{
  plot3b(X, grd, wire = FALSE, ...)
  if (is.null(grd$ObjQ))
    Q <- Obj2ObjQ(grd$Obj, grd)
  else Q = grd$ObjQ
  Q <- X2ObjQ(Q, X)
  rgl::wire3d(Q, col = "black", specular = "black")
}

#' plot3qs
#' @description
#' plot with coordinates C$X (after updateX) and scalar as color code
#' 
#' @param X Coordiates, e.g. in C$X from updateX()
#' @param grd grid with a basic rgl-object grd$Obj
#' @param s scalar to plot as color code on shape
#' @param pal (=heat.colors) color palette to use 
#' @examples
#' \dontrun{
#' data("M4",package = "MemRBC"); 
#' SetParams(M4)
#'  update(M4,"dA")->M4
#'  #plot area sizes as color code
#' plot3qs(updateX(M4$A,M4$grd,M4$bas)$X,M4$grd,M4$dA)
#' }
#' @export
plot3qs<-function (X, grd, s, alpha = 1, specular = "black", pal=heat.colors, ...)
{
  O <- X2Obj(grd$Obj, X)
  O <- Rvcg::vcgUpdateNormals(O)
  col = pal(100)[1 + 99 * (s - min(s))/diff(range(s))]
  rgl::shade3d(O, col = col, specular = specular, ...)
  if (is.null(grd$ObjQ))
    Q <- Obj2ObjQ(grd$Obj, grd)
  else Q = grd$ObjQ
  Q <- X2ObjQ(Q, X)
  rgl::wire3d(Q, col = "gray", specular = "black")
}


#' plot3b
#' @description
#' plot a 3d shape from coordinates C$X
#' @param X,grd : 3d coordinates and 2d grid object
#' @param ... : further plotting options, e.g. alpha=0.5 for semi-transparency
#' @export
plot3b<-function (X, grd, col = "white", specular = "black", wire = TRUE,                  ...)
{
  O <- X2Obj(grd$Obj, X)
  O <- Rvcg::vcgUpdateNormals(O)
  rgl::shade3d(O, col = col, specular = specular, ...)
  if (wire)
    rgl::wire3d(O, col = "black", lwd = 2, specular = "black")
}



#' MakeGrid_GaussLegendreSimpson
#' @description
#' compute an integration grid
#' @param n (=20) : number of points along one dimension
#' @param  ua,ub : usually 0,pi for u-interval
#' @param  va,vb : usually 0,2*pi for v-interval
#' @param comment : give your own comment (as.character)
#' @export
MakeGrid_GaussLegendreSimpson<-function(n=20,ua=0,ub=pi,va=0,vb=2*pi, comment="spherical coordinates Gauss-Legendre-Simpson grid, type GLS",check_plt=FALSE) # assume spherical coordinates
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
    k=k+1;l=(j-1)*nx+i; q[,k]=c(l,l+1,l+1+nx);
    k=k+1;  q[ ,k]=c(l,l+nx+1,l+nx)
  }
  q=q[,1:(k)]
  x=sin(grd$u)*cos(grd$v);y=sin(grd$u)*sin(grd$v);z=cos(grd$u)
  rgl::mesh3d(x,y,z,triangles=q, normals = list(x=x,y=y,z=z) ) -> M
  grd$Obj<-M
  grd$comment<-comment
  grd$type="GLS"
  grd$n=n;grd$nu=nu;grd$nv=nv
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

#'last -  return last n elements from v; v may be list
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
{ dms=dim(field_m);
  res=array(NA,dms[2:3]);
for (i in 1:dms[2]) for (j in 1:dms[3]) res[i,j]=int2d_scalar_GLS(field_m[,i,j],grd)
return(res)
}

#  used within c++
#' @export
int2d_matrix_cxx <- function (field_m,grd,bas) # needed for C++-code of gradients
{ F2<-array(field_m, c( grd$ndof,bas$Ai_max,3))
  dms=dim(F2);res=array(NA,dms[2:3]);
 for (i in 1:dms[2]) for (j in 1:dms[3]) res[i,j]=int2d_scalar_GLS(F2[,i,j],grd)
 return(res)
}

# needed in c++
#' @export
vectomat_cxx<-function(x,Aimax)  return(matrix(x,Aimax,3))

#' @export
vectoarr_cxx<-function(x,ndof,Aimax)  return(array(x,c(ndof,Aimax,3)))

#' FitAlm
#' @description
#' fit coefficients from 3d-coordinates
#' @param X,bas : input data and basis; bas$mask must be set to ecluded X points indices
#' @param WX (=1) : spatial weights, could better be sin(grd$U)
#' @export
FitAlm <- function (X, bas, WX = rep(1, nrow(X)))
{ A=FitAlm_Tikhonov(X,bas,lambda=0, WX = WX)
  return(A)
}


#' Weighted SPHARM fit
#' @export
Weighted_FitAlm <- function (X, bas, sigma = 0.001, CL = 0.95)
{
  mask = bas$mask
  if (is.null(mask)) {
    warning("no mask given in basis")
    X1 = X
    Y = bas$Ylm
  }
  else {
    X1 = X[-mask, ]
    Y = bas$Ylm[-mask, ]
  }
  M <- apply(X1, 2, mean)
  n <- dim(X1)[1]
  alpha <- 1 - CL
  p <- matrix(0, bas$L_max, 2)
  if (any(abs(M) > 1e-15))
    warning("masked X seems not centered")
  X1 <- apply(X1, 2, function(x) x - mean(x))
  A <- array(0, c(dim(bas$Ylm)[2], 3))
  Asmooth <- A
  Y <- Y[, 1:3]
  Ycommon <- pracma::inv(t(Y) %*% Y) %*% t(Y)
  for (k in 1:3) A[1:3, k] <- Ycommon %*% X1[, k]
  Asmooth[1:3, ] <- A[1:3, ]
  Xsmooth <- Y %*% A[1:3, ]
  Xestim <- Xsmooth
  L <- NA
  for (l in 2:bas$L_max) {
    cat("Fit order l=", l, "\r")
    Xresid <- X1 - Xestim
    s <- which(bas$LM[["l"]] == l)
    Y <- Y[, s]
    Ycommon <- pracma::inv(t(Y) %*% Y) %*% t(Y)
    betal <- matrix(0, 2 * l + 1, 3)
    for (k in 1:3) betal[, k] <- Ycommon %*% Xresid[, k]
    A[s, ] <- betal
    Asmooth[s, ] <- betal * exp(-l * (l + 1) * sigma)
    Xsmooth <- Y %*% betal
    Xestim <- Xestim + exp(-l * (l + 1) * sigma) * Xsmooth
    SSEk <- sum((X1 - Xestim)^2)
    SSEkm1 <- sum(Xresid^2)
    Fk <- (SSEkm1 - SSEk)/(2 * l + 1)/(SSEkm1/(n - (l + 1)^2))
    p_val <- pf(Fk, df1 = 2 * l + 1, df2 = (n - (l + 1)^2),
                lower.tail = FALSE)
    if (p_val > 1 - CL) {
      p[l, 1] = 1
      if (is.na(L)) {
        L = l
        Xresid_L = Xresid
      }
    }
    p[l, 2] = (SSEkm1 - SSEk)/(2 * l + 1)
  }
  cat("\n")
  colnames(p) = c("sig", "delta_SSE")
  A <- LM2A(A, bas)
  attr(A, "sigma_Weighted") <- sigma
  attr(A, "conf.level_Weighted") <- CL
  attr(Xestim, "sigma_Weighted") <- sigma
  attr(Xestim, "conf.level_Weighted") <- CL
  cat("sufficient L:", L, "\n")
  return(list(A = A, X = Xestim, Asmooth = Asmooth, p.value = p,
              sufficient_L = L, Xresid_L = ifelse(is.na(L), NA, Xresid_L),
              sigma = sigma))
}

#' FitAlm_Tikhonov
#' fit with regularization and weights:
#' filtering high frequencies in least squares
#'   lambda should be tested systematically by L curve discussion
#' @export
FitAlm_Tikhonov<-function (X, bas, lambda = 0,
                           WX = rep(1, nrow(X)), keepIM = FALSE,
                           newIM = FALSE)
{
  mask = bas$mask
  if (is.null(mask))
    stop("no mask in bas; create at least bas$mask=1")
  A <- matrix(0, bas$Ai_max + 1, 3)
  if (any(WX < 0))
    warning("negative spatial weights in FitAlm_Tikhonov - abs(WX) is taken")
  WX = abs(WX)
  WX = WX[-mask]
  X1 <- X[-mask, ]
  if (is.null(bas$IM))
    newIM = TRUE
  if (newIM) {
    Y = cbind(1, bas$Ylm[-mask, ])
    YtW = t(Y) ; # dim 1 is coeff., dim 2 is spatial
    for (i in 1:dim(YtW)[1]) YtW[i,]=YtW[i,]*WX[i]
    B = YtW %*% Y
    InvB1 = pracma::inv(B + lambda * diag(c(0, bas$G.tk)))
    IM = InvB1 %*% YtW
  }
  else IM = bas$IM
  for (k in 1:3) A[, k] = IM %*% X1[, k]
  cat("dropped A(l=0):", A[1, ], "\n")
  A <- LM2A(A[-1, ], bas)
  attr(A, "lambda_Tikhonov") = lambda
 # attr(A, "Fit spatial weights") = WX
  if (keepIM)
    attr(A, "IM") <- IM
  return(A)
}

#' inv_sph
#' @description
#' give angles  (u,v) from spherical or star-like shape 3D point X
#' @export
inv_sph<-function(X)
{ r=sqrt(sum(X^2)) # see https://mathworld.wolfram.com/SphericalCoordinates.html
return(c(acos(X[3]/r),atan2(X[2],X[1]))) # atan2 takes care of octants
}

#' radial_uv
#' give angles u,v from a starlike 3d-object, to be centred
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

#' synth
#' synthesize, but only for one spatial component in coefficients  A==A[,k] !!!
#'   for dim in 1..3: see synthX
#' @export
synth<-function(Y,A) { return(Y%*%A)}

# partial synthesis
synth_s<-function(Y,A,mx) { return(Y[,1:mx]%*%A[1:mx])}

#' synth12
#' a helper for faster finite differences, not yet used
#' @export
synth12<-function (A, C, i, j, k, del)
{
  cat(sum(C$X_u), "-> \t")
  n = dim(A)[1]
  if (i > 0) {
    C$X[, k] = C$X[, k] - bas$Ylm[, i] * del
    C$X_u[, k] = C$X_u[, k] - bas$Ylm_u[, i] * del
    C$X_v[, k] = C$X_v[, k] - bas$Ylm_v[, i] * del
    C$X_uu[, k] = C$X_uu[, k] - bas$Ylm_uu[, i] * del
    C$X_uv[, k] = C$X_uv[, k] - bas$Ylm_uv[, i] * del
    C$X_vv[, k] = C$X_vv[, k] - bas$Ylm_vv[, i] * del
  }
  if (j <= n) {
    C$X[, k] = C$X[, k] + bas$Ylm[, j] * del
    C$X_u[, k] = C$X_u[, k] + bas$Ylm_u[, j] * del
    C$X_v[, k] = C$X_v[, k] + bas$Ylm_v[, j] * del
    C$X_uu[, k] = C$X_uu[, k] + bas$Ylm_uu[, j] * del
    C$X_uv[, k] = C$X_uv[, k] + bas$Ylm_uv[, j] * del
    C$X_vv[, k] = C$X_vv[, k] + bas$Ylm_vv[, j] * del
  }
  cat(sum(C$X_u), "\n")
  return(C)
}

#' synth_update
#' @description
#' for forward differences, only add relevant del*Ylm
#' use this in Hessian computation by finite difference: 
#'     fwd: C=synth_update...(C,bas,i,k, +1e-6)
#'   reset: C=synth_update...(C,bas,i,k, -1e-6)
#' @param   k: spatial dim.
#' @param   i: spectral order as from (l,m)
#' @param del (=1e-6) : steplength of update
#' @export
synth_update<-function (C, bas, i, k,  del=1e-6)
{
    C$X[, k]    = C$X[, k] + bas$Ylm[, i] * del
    C$X_u[, k]  = C$X_u[, k] + bas$Ylm_u[, i] * del
    C$X_v[, k]  = C$X_v[, k] + bas$Ylm_v[, i] * del
    C$X_uu[, k] = C$X_uu[, k] + bas$Ylm_uu[, i] * del
    C$X_uv[, k] = C$X_uv[, k] + bas$Ylm_uv[, i] * del
    C$X_vv[, k] = C$X_vv[, k] + bas$Ylm_vv[, i] * del
  return(C)
}

#' synth_update_inplace
#' experimental, not faster
#' k: spatial dim.
#' i: spectral order (l,m)
#' use: fwd: C=synth_update...(A,C,bas,i,k, +1e-6)
#'      reset: C=synth_update...(A,C,bas,i,k, -1e-6)
#' not faster than returning full changed C
#' @export
synth_update_inplace<-function (C, bas, i, k, del=1e-6)
{  assign(deparse(substitute(C)), {

      C$X[, k]    = C$X[, k] + bas$Ylm[, i] * del
      C$X_u[, k]  = C$X_u[, k] + bas$Ylm_u[, i] * del
      C$X_v[, k]  = C$X_v[, k] + bas$Ylm_v[, i] * del
      C$X_uu[, k] = C$X_uu[, k] + bas$Ylm_uu[, i] * del
      C$X_uv[, k] = C$X_uv[, k] + bas$Ylm_uv[, i] * del
      C$X_vv[, k] = C$X_vv[, k] + bas$Ylm_vv[, i] * del
      C
  }, envir = rlang::env_parent())

}

# experimental, not faster
#' scale_inplace
#' @export
scale_inplace<-function (m, s)
{
  assign(deparse(substitute(m)), {
    ifelse(is.numeric(s), {
      m$A = m$A * s
      "scaled m in place"
    }, NULL)
  }, envir = rlang::env_parent())
}

#' deltaX_norm
#' @description
#' Norm of coordinate change between two coordinate objects 
#' @param C1,C2 : coordinate objects 
#' @export
deltaX_norm<-function(C1,C2)
{ return(apply(C2$X-C1$X,2,pracma::Norm))
}
#' saveA
#' save a set of coefficients A to file
#' @export
saveA<-function(A,file)
{  save(A,file=file) }

# load a set of coeffs, make it compatible for the current basis bas.
#' @export
loadA<-function(file, bas) # loads amplitudes and stores according to basis bas
{   # if bas is larger, empty amplitudes are kept zero
  load(file)->n # was saved from A
  A2<-get(n)
  A<-array(0,c(bas$Ai_max,3));
  L=min(dim(A)[1],dim(A2)[1])
  for (k in 1:3) A[1:L,k]<-A2[1:L,k] # copy only first coefficients
  class(A)<-"MemA"
  return(A)
}

#' loadAlm
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

#' rotateX
#' rotation first around x by px, then y then z-axis
#' @export
rotateX<-function (X, px = pi, py = 0, pz = 0, transpose = FALSE)
{
  R.mz = matrix(c(cos(pz), -sin(pz), 0, sin(pz), cos(pz), 0,
                  0, 0, 1), 3, 3)
  R.mx = matrix(c(1, 0, 0, 0, cos(px), -sin(px), 0, sin(px),
                  cos(px)), 3, 3)
  R.my = matrix(c(cos(py), 0, -sin(py), 0, 1, 0, sin(py), 0,
                  cos(py)), 3, 3)
  R <- R.mz %*% R.my %*% R.mx
  if (transpose)
    R <- t(R)
  if (is.vector(X))
    return(R %*% X)
  else if (is.matrix(X))
    return(t(apply(X, 1, function(x) R %*% x)))
  else stop("X is not a 3d vector, nor a matrix of 3d vectors")
}

# internal use
rotateXbyM<-function(X,M)
{ return(M%*%X)}

#' rotateA
#' rotate coefficients by rotating coordinates
#'  gives back a rotation error (new from fit vs. rotated coord)
#'  rotation order is by X,Y,Z-axis
#'
#' @export
rotateA <-function (A, bas, grd, px = pi/2, py = -pi/2, pz = pi/2, plt = FALSE)
{ C <- updateX(A, grd, bas)
  if (plt)
    id1 <- plot3b(C$X, grd)
  R.mz = matrix(c(cos(pz), -sin(pz), 0, sin(pz), cos(pz), 0,
                  0, 0, 1), 3, 3)
  R.mx = matrix(c(1, 0, 0, 0, cos(px), -sin(px), 0, sin(px),
                  cos(px)), 3, 3)
  R.my = matrix(c(cos(py), 0, -sin(py), 0, 1, 0, sin(py), 0,
                  cos(py)), 3, 3)
  R <- R.mz %*% R.my %*% R.mx
  X1 = C$X
  for (i in 1:dim(X1)[1]) X1[i, ] = R %*% C$X[i, ]
  A1 <- A
  if (!is.matrix(bas$IM))
  A1[, k] <- FitFast(bas,X1)
  C1 <- updateX(A1, grd, bas)
  if (plt) {
    rgl::open3d()
    id2 <- rgl::plot3b(C1$X, grd)
  }
  return(list(A = A1, C = C1,
              rot_err = pracma::Norm(C1$X-X1)))
}

# plotLseries
# shows a series of up to nr x nc images and "rep" values for the l>0 present in bas
# (truncation plot)
# not exported, but PlotLSeries 
plotLseries<-function (nr = 4, nc = 5, A, C, grd, bas,
                       Vals = TRUE,
                       fill = TRUE, rep = "H2", S = NULL,
                       Ref = NULL, stretch = !is.null(Ref))
{
  rgl::mfrow3d(nr, nc, sharedMouse = TRUE)
  if (Vals)
    M = matrix(0, length(unique(bas$LM[, 1])), 7)
  k = 1
  for (l in bas$Lset) {
    w = which(bas$LM[, 1] == l)
    last = w[length(w)]
    Y = synthX(bas$Ylm[, 1:last], A[1:last, ])
    A1 = A
    A1[] = 0
    A1[1:last, ] = A[1:last, ]
    C2 <- updateX(A1, grd, bas)
    h2 <- E_SCM(A1, grd, bas, C2)
    if (stretch)
      E = E_SEN(A, grd, bas, S, Ref)
    else E = 0
    if (Vals) {
      M[k, 1] = h2$Area
      M[k, 2] = h2$Volume
      M[k, 3] = h2$Curv
      M[k, 4] = h2$Wb
      M[k, 5] = h2$H2_BC
      if (stretch)
        M[k, 6] = E
      M[k, 7] = sum(abs(A[bas$LM[, 1] == l, ]))
    }
    if (is.null(S))
      S = list(alpha = 0)
    plot3qs(Y, grd, S$alpha)
    rgl::title3d(paste(l, ": E",
                       round((h2$Wb + E)/MemRBC_env$M.Es,
                              3), " C", round(h2$Curv, 3)))
    k = k + 1
    if (k > nr * nc)
      (break)()
    if (l < last(bas$Lset))
      rgl::next3d()
  }
  if (Vals) {
    M1 = apply(M, 2, diff)
    M2 = apply(M1, 2, function(x) (x)/diff(range(x)))
    matplot(M2, type = "l", lty = 1, lwd = 2, ylab = expression(Delta *
                                                                  Q), pch = 20)
    matpoints(M2, pch = 20)
    legend("topright", col = 1:6, lty = 1, lwd = 2, leg = c("A",
                                                            "V", "C", "H2", "H2_BC", "Nrm"), pch = 20)
    rownames(M) = paste("l=",as.character(bas$Lset))
    colnames(M) = c("A", "V", "C", "H2", "H2_BC", "E","||A_l||")
    M = as.data.frame(M)
    return(M)
  }
  else return("Plot L Series done")
}


#
# fill e.g. grd$Obj  with new coords
#  and return 3D object, plot if wanted
#' MakeOBJ
#' make a rgl-object from standard objects A,grd,bas,C
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

#' Membrane_LaplacianOBJ
#'  Laplacian matrix from graph
#' @export
Membrane_LaplacianOBJ <- function(X)
{
  # this procedure doubles all entries in L
  ia=c(X$it[1,],X$it[2,],X$it[3,]);ja=c(X$it[2,],X$it[3,],X$it[1,])
  N<-max(max(ia),ja); TN <- Matrix::sparseMatrix(dims=c(N,N),i=ia,j=ja,x=rep(1,length(ia)),
                                         use.last.ij=FALSE)
  ig <- igraph::graph_from_adjacency_matrix(TN,mode="undirected",diag = FALSE)
  # remove doubles by /2
  return(igraph::laplacian_matrix(ig,normalization="unnormalized",sparse = TRUE)/2)
}

#' Membrane_LaplaciansOBJ
#' from 3D-object X compute mesh Laplacian L, normalized Ln, Diagonal D
#' @export
Membrane_LaplaciansOBJ <- function(X)
{
  ia=c(X$it[1,],X$it[2,],X$it[3,]);ja=c(X$it[2,],X$it[3,],X$it[1,])
  N<-max(max(ia),ja); TN <- Matrix::sparseMatrix(dims=c(N,N),i=ia,j=ja,x=rep(1,length(ia)),
                                         use.last.ij=FALSE)
  ig <- igraph::graph.adjacency(TN,mode="undirected",diag = FALSE)
  # recompute  u directly from Laplacian
  L=igraph::laplacian_matrix(ig,normalization = "unnormalized",sparse = TRUE)/2
  D=L;diag(D)<-0
  return(list(L=L,Ln=igraph::laplacian_matrix(ig,normalized = TRUE,sparse = TRUE),D=D))
}

#' Membrane_Laplacian_cotan
#' slow implementation in R
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
      cot<-sum(u1*v1)/sqrt(sum(pracma::cross(u1,v1)^2))
      M[i,j] <- M[i,j] + cot;
      M[j,i] <- M[j,i] + cot;
      M[i,i] <- M[i,i] - cot;
      M[j,j] <- M[j,j] - cot
    }
    if (tr%% 100 ==0 )cat(tr/N.t,"  progress cotan-Laplacian  \r")
  };
  return(-M/2) # comparable to L from GEMINI-implementation-draft
}

#' GEMINI_cotan_Laplacian_II
#' GEMINI created cotan Laplacian from rgl-mesh object
#' @export
GEMINI_cotan_Laplacian_II <- function(mesh) {
  verts <- t(mesh$vb[1:3, ])
  faces <- t(mesh$it)

  n_verts <- nrow(verts)
  idx_A <- faces[, 1]
  idx_B <- faces[, 2]
  idx_C <- faces[, 3]
  A <- verts[idx_A, ]
  B <- verts[idx_B, ]
  C <- verts[idx_C, ]

  u_A <- B - A
  v_A <- C - A
  u_B <- C - B
  v_B <- A - B
  u_C <- A - C
  v_C <- B - C
  dot_A <- rowSums(u_A * v_A)
  dot_B <- rowSums(u_B * v_B)
  dot_C <- rowSums(u_C * v_C)

  cross_x <- u_A[,2] * v_A[,3] - u_A[,3] * v_A[,2]
  cross_y <- u_A[,3] * v_A[,1] - u_A[,1] * v_A[,3]
  cross_z <- u_A[,1] * v_A[,2] - u_A[,2] * v_A[,1]

  two_area <- sqrt(cross_x^2 + cross_y^2 + cross_z^2)
  two_area[two_area < 1e-12] <- 1e-12
  cot_A <- dot_A / two_area
  cot_B <- dot_B / two_area
  cot_C <- dot_C / two_area
  I <- c(idx_B, idx_C,  idx_C, idx_A,  idx_A, idx_B)
  J <- c(idx_C, idx_B,  idx_A, idx_C,  idx_B, idx_A)
  W_vals <- c(cot_A, cot_A, cot_B, cot_B, cot_C, cot_C) * 0.5
  W <- Matrix::sparseMatrix(i = I, j = J, x = W_vals, dims = c(n_verts, n_verts))
  diag_vals <- Matrix::rowSums(W)
  D <- Matrix::sparseMatrix(i = 1:n_verts, j = 1:n_verts, x = diag_vals, dims = c(n_verts, n_verts))
  L <- D - W
  return(L)
}

#' Gemini_cotan_Laplacian_I
#' alternative version
#' @export
Gemini_cotan_Laplacian_I <- function(mesh) {
  V <- t(mesh$vb[1:3, ])
  F <- t(mesh$it)
  n_verts <- nrow(V)
  n_faces <- nrow(F)
  compute_cotans <- function(i_A, i_B, i_C) {
    u <- V[i_B, ] - V[i_A, ]
    v <- V[i_C, ] - V[i_A, ]
    dot_uv <- rowSums(u * v)
    cross_prod <- cbind(
      u[,2]*v[,3] - u[,3]*v[,2],
      u[,3]*v[,1] - u[,1]*v[,3],
      u[,1]*v[,2] - u[,2]*v[,1]
    )
    norm_cross <- sqrt(rowSums(cross_prod^2))
    return(dot_uv / norm_cross)
  }
  i1 <- F[, 1]
  i2 <- F[, 2]
  i3 <- F[, 3]
  cot1 <- compute_cotans(i2, i3, i1) # Angle at 1 is between (1-2) and (1-3)?
  c1 <- compute_cotans(i1, i2, i3) # Angle at i1
  c2 <- compute_cotans(i2, i3, i1) # Angle at i2
  c3 <- compute_cotans(i3, i1, i2) # Angle at i3
  I <- c(i1, i2, i2, i3, i3, i1)
  J <- c(i2, i1, i3, i2, i1, i3)
  W <- c(c3, c3, c1, c1, c2, c2) * 0.5
  L_off <- Matrix::sparseMatrix(i = I, j = J, x = W, dims = c(n_verts, n_verts))
  diagonal_vals <- Matrix::rowSums(L_off)
  L_diag <- Matrix::Diagonal(x = diagonal_vals)
  L <- L_diag - L_off
  return(L)
}

#' Membrane_Eig
#' @description
#' get Eigensystem for mesh Laplacian (created from an rgl object -> Membrane_LaplacianOBJ())
#' @param L: mesh Laplacian
#' @param which (=1:4) : which values to return
#' @param kind (=SM) : e.g. for smallest magnitude order 
#' @return matrix of demanded eigenvectors (columns); eigenvalues given as attribute "values"
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

#' imag.delta.arrowws.pca
#' plot arrows of change between two shape states
#' @param A1,A2 coefficients for the shapes
#' @param grd2,bas2 grid and basis of second shape A2
#' @param O1 3d object for A1
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


#' imag.delta.aligned
#' @description
#' svd alignment of shape coefficients into 3d-objects, see
#' https://www.cse.wustl.edu/~taoju/cse554/lectures/lect07_Alignment.pdf
#' and plotting of spatial difference vectors
#'
#' @param A1,A2 input coefficients
#' @param grd2,bas2 grid and basis of A2
#' @returns obj1,obj2 aligned 3d-objects
#' @examples
#' data("M4",package = "MemRBC")
#' M4p=M4
#' M4p$A=M4$A+rnorm(M4$bas$Ai_max*3,sd=0.4)
#' imag.delta.aligned(M4p$A,M4$A,M4$grd,M4$bas)->L
#' rgl::open3d()
#' rgl::plot3d(L$obj1,alpha=0.6,col=1,aspect=FALSE)
#' rgl::plot3d(L$obj2,alpha=0.6,col=2,add=TRUE)
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
#' SE_2_OBJ
#' @description
#' Read data from surface evolver dump file.
#' The vertices are renumbered here, taking
#' more computing time.
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

#' SErenumbered_2_OBJ
#' @description
#' Read data from surface evolver dump file
#'  special case: faster if SE renumbered all before the dump
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

#' GEMINI_disk_conformal_map 
#' @description
#' Map a 3D Mesh to a 2D Unit Circle (Disk Conformal Map)
#'
#' @param mesh An Rvcg mesh3d object (must be an open mesh, not a closed sphere)
#' @param L The Cotangent Laplacian matrix. If NULL, it is computed.
#' @param spherical (=TRUE) for giving not uv on a disk but spherical coordinates (u,v)
#' @return A list containing the 2D coordinates and the flattened mesh.
#' @export
GEMINI_disk_conformal_map <- function(mesh, L = NULL,plt=FALSE,spherical=TRUE) {
  num_verts <- ncol(mesh$vb)

  if (is.null(L)) {
    L <- .GEMINI_get_cotan_Laplacian_cxx(mesh)
  }
  all_edges <- Rvcg::vcgGetEdge(mesh)
  border_edges <- all_edges[all_edges$border == 1, 1:2]

  if (nrow(border_edges) == 0) {
    stop("Mesh has no boundary! (Is it a closed sphere? You must cut it first.)")
  }
  g_border <- igraph::graph_from_edgelist(as.matrix(border_edges), directed = FALSE)

  start_node <- border_edges[1, 1]
  dfs_res <- igraph::dfs(g_border, root = start_node, unreachable = FALSE)
  b_indices <- as.numeric(dfs_res$order)
  b_indices <- b_indices[!is.na(b_indices)]

  n_b <- length(b_indices)
  cum_len <- 0
  thetas <- numeric(n_b)
  coords <- t(mesh$vb[1:3, ])

  for (i in 2:n_b) {
    v1 <- coords[b_indices[i-1], ]
    v2 <- coords[b_indices[i], ]
    dist <- sqrt(sum((v1 - v2)^2))
    cum_len <- cum_len + dist
    thetas[i] <- cum_len
  }
  total_len <- cum_len + sqrt(sum((coords[b_indices[n_b],] - coords[b_indices[1],])^2)) # Close loop
  thetas <- (thetas / total_len) * 2 * pi
  u_fixed <- cos(thetas)
  v_fixed <- sin(thetas)
  if(spherical)
  { v_fixed[v_fixed<0]=-1
    v_fixed[v_fixed>0]=1
    v_fixed[1]=0
    v_fixed[which.min(u_fixed)]=0
  }
  all_indices <- 1:num_verts
  interior_indices <- setdiff(all_indices, b_indices)
  L_ii <- L[interior_indices, interior_indices]
  L_ib <- L[interior_indices, b_indices]
  str(L_ib)
  str(u_fixed)
  rhs_u <- -L_ib %*% u_fixed
  rhs_v <- -L_ib %*% v_fixed
  u_in <- solve(L_ii, rhs_u[,1])
  v_in <- solve(L_ii, rhs_v[,1])
  u_full <- numeric(num_verts)
  v_full <- numeric(num_verts)

  u_full[b_indices] <- u_fixed
  u_full[interior_indices] <- as.vector(u_in)

  v_full[b_indices] <- v_fixed
  v_full[interior_indices] <- as.vector(v_in)
  flat_mesh <- mesh
  flat_mesh$vb[1,] <- u_full
  flat_mesh$vb[2,] <- v_full
  flat_mesh$vb[3,] <- 0
  if(plt)  plot(u_full,v_full,pch=".",cex=2)
  if(spherical) { u_full=u_full*pi/2+pi/2;v_full=v_full*pi+pi}
  return(list(uv = data.frame(u=u_full, v=v_full), mesh_flat = flat_mesh, bnd=b_indices))
}

#' GEMINI_cut_mesh_along_path
#' @description
#' Cut a Mesh Along a Path of Vertices
#'
#' Opens a closed mesh by "unzipping" it along a sequence of vertices.
#' Interior vertices of the path are duplicated. Endpoints remain as "hinges".
#'
#' @param mesh An Rvcg mesh3d object.
#' @param path A vector of integer vertex indices representing the cut path.
#'             Must be a connected sequence of edges.
#' @return A new mesh3d object with the cut applied (more vertices, updated faces).
#' @export
GEMINI_cut_mesh_along_path <- function(mesh, path) {
  path=as.integer(path)
  if(length(path) < 3) stop("Path must have at least 3 vertices to define a cut.")
  p_inner <- path[2:(length(path)-1)]
  n_dupes <- length(p_inner)
  new_idx_map <- rep(0, ncol(mesh$vb))
  new_coords <- mesh$vb[, p_inner]
  start_new_idx <- ncol(mesh$vb) + 1
  new_indices <- start_new_idx:(start_new_idx + n_dupes - 1)
  faces <- t(mesh$it)
  get_faces_on_edge <- function(u, v, face_mat) {
    w1 <- which(face_mat == u, arr.ind = TRUE)[,1]
    w2 <- which(face_mat == v, arr.ind = TRUE)[,1]
    intersect(w1, w2)
  }
  get_right_sector_faces <- function(v, start_face_idx, start_edge_v1, stop_edge_v2, all_faces) {
    sector_faces <- c()
    current_face_idx <- start_face_idx
    iter <- 0
    max_iter <- 50
    while(iter < max_iter) {
      sector_faces <- c(sector_faces, current_face_idx)
      f_verts <- all_faces[current_face_idx, ]
      if (stop_edge_v2 %in% f_verts) {
        return(sector_faces)
      }
      other_v <- f_verts[ !f_verts %in% c(v, ifelse(iter==0, start_edge_v1, prev_neighbor_v)) ]
      candidates <- get_faces_on_edge(v, other_v, all_faces)
      next_face <- setdiff(candidates, current_face_idx)
      if(length(next_face) == 0) {
        break
      }
      prev_neighbor_v <- other_v
      start_edge_v1 <- other_v
      current_face_idx <- next_face[1]
      iter <- iter + 1
    }
    return(sector_faces)
  }

  faces_to_update <- list()
  pivot <- path[2]
  prev_v <- path[1]
  next_v <- path[3]

  init_faces <- get_faces_on_edge(prev_v, pivot, faces)
  current_right_face <- init_faces[1]
  updates <- matrix(0, nrow=0, ncol=3)
  for (i in 1:length(p_inner)) {
    v_idx <- p_inner[i]
    v_new_idx <- new_indices[i]
    v_prev <- path[i]
    v_prev <- path[i]
    v_curr <- path[i+1] # This is v_idx
    v_next <- path[i+2]
    sector_faces <- get_right_sector_faces(v_curr, current_right_face, v_prev, v_next, faces)
    if(length(sector_faces) > 0) {
      new_rows <- cbind(sector_faces, rep(v_curr, length(sector_faces)), rep(v_new_idx, length(sector_faces)))
      updates <- rbind(updates, new_rows)
      current_right_face <- sector_faces[length(sector_faces)]
    }
  }
  mesh$vb <- cbind(mesh$vb, new_coords)
  new_it <- mesh$it
  if(nrow(updates) > 0) {
    for(r in 1:nrow(updates)) {
      f_idx <- updates[r, 1]
      old_v <- updates[r, 2]
      new_v <- updates[r, 3]
      col_vals <- new_it[, f_idx]
      match_pos <- which(col_vals == old_v)
      if(length(match_pos) > 0) {
        new_it[match_pos, f_idx] <- new_v
      }
    }
  }
  mesh$it <- new_it
  mesh <- Rvcg::vcgUpdateNormals(mesh)
  return(mesh)
}

#' NorthSouth
#' find two extremal vertices along PC1 and some 5 paths to select from
#' @param O 3d object
#' @param w (=1) : <15, integer, pre-selecting which proposed path to take 
#' @return NS=c(N,S), North and South vertex indices
#' @return Path: list of index vectors connecting N and S by k-shortest path
#' @examples
#' data("SF4lr",package = "MemRBC")
#' SF4lr->O
#' NorthSouth(O,13)->P # take 13th path
#' P$Paths # candidate paths of vertices for cutting
#' P$LVPath # 13th lowest var(z) path
#' PlotPaths(O,P)
#' O1 <- GEMINI_cut_mesh_along_path(O,P$LVPath)
#' @export
NorthSouth<-function(O,w=1){
 Y=princomp(Obj2X(O))$scores
 N=which.max(Y[,1])
 S=which.min(Y[,1])
 ia = c(O$it[1, ], O$it[2, ], O$it[3, ])
 ja = c(O$it[2, ], O$it[3, ], O$it[1, ])
  N <- max(max(ia), ja)
 TN <- Matrix::sparseMatrix(dims = c(N, N), i = ia, j = ja,
                           x = rep(1, length(ia)), use.last.ij = FALSE)
 ig <- igraph::graph_from_adjacency_matrix(TN, mode = "undirected",
                                          diag = FALSE)
 P <- igraph::k_shortest_paths(ig, from = N, to = S, k=15)
 v=sort(sapply(P$vpaths,function(x) {var(Y[x,3])}))
 return(list(NS=c(N,S),Paths=P,vars=v, LVPath=(P$vpaths[[ order(v)[w] ]])))
}

#' PlotPaths
#' @description
#'   draw some paths P on object O
#' @export
PlotPaths<-function(O,P,LVPath=TRUE)
{
  rgl::shade3d(O,alpha=0.6,col="white")
  if (LVPath) rgl::plot3d(Obj2X(O)[P$LVPath,],type="l",add=TRUE,col=1,lwd=6)
  for (i in 1:length(P$Paths$vpaths))
   rgl::plot3d(Obj2X(O)[unlist(P$Paths$vpaths[[i]]),],type="l",add=TRUE,col=i+1,lwd=1)
}

#
#' Brechbuehler.Init.uv.2
#' @description
#' computes spherical coordinates (u,v) for vertices of a 3D-object X
#' The original algorithm from Brechbühler (1995) is implemented,
#' but spherical areas are not iterated for refinement.
#' The resulting (u,v) may not be optimal for fitting.
#' A better result comes from algorithm GEMINI_disk_conformal_map()
#' @export
Brechbuehler.Init.uv.2<-function (X1, Fit_order = 12, InitFit = FALSE, poles.axis = 2,
                                  mat.mode = c("cotan.lapl", "cotan.chi", "euclid"), file.out = "Brechbuehler-init-uv-2.obj",
                                  rxy = 0.01)
{
  X = X1
  {
    for (k in 1:3) X$vb[k, ] = X$vb[k, ] - mean(X$vb[k, ])
    x = t(X$vb[1:3, ])
    tri = X$it
    str(tri)
    deg.v = max(table(X$it))
    n.v = dim(X$vb)[2]
    nn = array(NA, c(n.v, deg.v))
    dim(nn)
    for (i in 1:n.v) {
      nn.i = c()
      w.i = c(tri[, which(tri[1, ] == i)], tri[, which(tri[2,
      ] == i)], tri[, which(tri[3, ] == i)])
      nn.i = unique(w.i)
      nn.i = setdiff(nn.i, i)
      nn[i, 1:(length(nn.i))] <- nn.i
    }
    X.tri = x
    k.which = poles.axis
    if (k.which < 4) {
      x1 <- princomp(t(X1$vb[1:3, ]))$scores[, 1:3]
      plot(x1[, 1:2])
      if (k.which == -3) {
        w = (sqrt(x1[, 1]^2 + x1[, 2]^2) < rxy)
        f1p = which.min(x1[w, 3])
        f2p = which.max(x1[w, 3])
        f1 = which(x1[, 3] == x1[w, 3][f1p])
        f2 = which(x1[, 3] == x1[w, 3][f2p])
        points(x1[w, 1:2], col = 3, pch = 19)
        points(x1[c(f1, f2), 1:2], col = 2, pch = 19,
               cex = 1.3)
      }
      else {
        f1 = which.min(x1[, k.which])
        f2 = which.max(x1[, k.which])
      }
    }
    else {
      cat("\a")
      print("PICK NORTH AND SOUTH POLE centrally in plot")
      mds <- princomp((t(X1$vb[1:3, ])))$scores
      w.upper = which(mds[, 1] > 0)
      w.lower = which(mds[, 1] < 0)
      dev.new()
      plot(mds[w.upper, 2:3], main = "PICK NORTH AND SOUTH POLE from PCA (x,y) upper z here")
      p = locator(1)
      f1 = which.min(abs(mds[w.upper, 2] - p$x[1]) + abs(mds[w.upper,
                                                             3] - p$y[1]))
      plot(mds[w.lower, 2:3], main = "PICK NORTH AND SOUTH POLE from PCA (x,y) lower z here")
      points(mds[w.upper, 2:3], col = "green", pch = ".",
             cex = 2.5)
      points(mds[w.upper[f1], 2:3], pch = 20, cex = 1.5,
             col = "red")
      p = locator(1)
      f2 = which.min(abs(mds[w.lower, 2] - p$x[1]) + abs(mds[w.lower,
                                                             3] - p$y[1]))
      cat("PICKED (upper/lower indices)", f1, " ", f2,
          "\n")
      g1 = c(which(mds[, 1] == mds[w.upper[f1], 1]), which(mds[,
                                                               2] == mds[w.upper[f1], 2]), which(mds[, 3] ==
                                                                                                   mds[w.upper[f1], 3]))
      g2 = c(which(mds[, 1] == mds[w.lower[f2], 1]), which(mds[,
                                                               2] == mds[w.lower[f2], 2]), which(mds[, 3] ==
                                                                                                   mds[w.lower[f2], 3]))
      if (table(g1) == 3)
        g1 = g1[1]
      else stop("Error picking g1")
      if (table(g2) == 3)
        g2 = g2[1]
      else stop("Error picking g2")
      f1 = g1
      f2 = g2
      cat("PICKED (vertex indices) ", f1, " ", f2, "\n")
    }
    (p.fix = c(f1, f2))
    n.m = setdiff(1:n.v, p.fix)
    n = n.v - 2
    A = matrix(0, n, n)
    renum = 1:n
    names(renum) = as.character(n.m)
    NN = NN1 = list()
    for (i in 1:n.v) {
      NN[[i]] <- na.omit(nn[i, ])
      NN1[[i]] <- setdiff(NN[[i]], p.fix)
    }
    k = 1
    for (i in n.m) {
      A[k, k] = length(NN[[i]])
      k = k + 1
    }
    k = 1
    for (i in n.m) {
      A[k, renum[as.character(NN1[[i]])]] <- -1
      k = k + 1
    }
    range(A - t(A))
    b = rep(0, length(n.m))
    nn.sp = NN[[p.fix[2]]]
    b[renum[as.character(nn.sp)]] <- pi
    u.s <- solve(A, b)
    u <- rep(0, n.v)
    u[f1] <- 0
    u[f2] <- pi
    u[n.m] <- u.s
    plot(u)
    u.bb.numb = c(0, u.s, pi)
    u.bb = u
  }
  u.0 = u
  ia = c(X$it[1, ], X$it[2, ], X$it[3, ])
  ja = c(X$it[2, ], X$it[3, ], X$it[1, ])
  N <- max(max(ia), ja)
  TN <- Matrix::sparseMatrix(dims = c(N, N), i = ia, j = ja,
                             x = rep(1, length(ia)), use.last.ij = FALSE)
  ig <- igraph::graph_from_adjacency_matrix(TN, mode = "undirected",
                                            diag = FALSE)
  B = rep(0, n.v)
  L = igraph::laplacian_matrix(ig, normalized = FALSE, sparse = TRUE)
  L.ret = L
  L[f1, ] = 0
  L[f1, f1] = 1
  B[f1] = 0
  L[f2, ] = 0
  L[f2, f2] = 1
  B[f2] = pi
  u <- solve(L, B)
  plot(u, u.0, main = "compare Brechbuehler and Laplace-solver")
  (north.south <- igraph::shortest_paths(ig, from = f1, to = f2,
                                         output = "vpath")$vpath[[1]])
  if (TRUE)
    for (i in 1:3) {
      here = f1
      north.south.new = as.numeric(north.south)
      cnt = 2
      while (here != f2) {
        n.here = as.numeric(igraph::neighbors(ig, here))
        take = n.here[which.max(u[n.here])]
        north.south.new[cnt] = take
        here = take
        cnt = cnt + 1
      }
      north.south = north.south.new
    }
  halo = unique(unlist(igraph::ego(ig, order = 1, north.south)))
  pure.halo = setdiff(halo, north.south)
  col = rep("black", n.v)
  col[north.south] = "green"
  col[pure.halo] = "red"
  rgl::wire3d(X, col = col, lwd = 2)
  EW = rep("", n.v)
  here = f1
  cnt = 2
  done = f1
  while (here != f2) {
    n.here = as.numeric(igraph::neighbors(ig, here))
    take = n.here[which.max(u[n.here])]
    n.take = igraph::neighbors(ig, take)
    done = c(done, take)
    (check = intersect(n.here, n.take))
    dV = X$vb[1:3, take] - X$vb[1:3, here]
    dN = X$normals[1:3, here]
    for (chk in check) {
      dS = X$vb[1:3, chk] - X$vb[1:3, here]
      D = det(matrix(c(dV, dS, dN), 3, 3))
      EW[chk] = ifelse(D < 0, "E", "W")
    }
    here = take
    cnt = cnt + 1
  }
  table(EW)
  sum(EW != "")
  rgl::clear3d()
  rgl::wire3d(X, col = col, lwd = 2)
  rgl::spheres3d(t(X$vb[1:3, EW == "W"]), rad = 0.05, col = "green")
  rgl::spheres3d(t(X$vb[1:3, EW == "E"]), rad = 0.05, col = "blue")
  all(which(EW != "") %in% halo)
  rest = setdiff(pure.halo, which(EW != ""))
  length(rest)
  length(pure.halo)
  East = which(EW == "E")
  West = which(EW == "W")
  n.n = igraph::neighbors(ig, f1)
  mark = setdiff(n.n, which(EW != ""))
  EW[mark] = " "
  n.s = igraph::neighbors(ig, f2)
  mark = setdiff(n.s, which(EW != ""))
  EW[mark] = " "
  rest = setdiff(rest, c(f1, f2))
  for (i in rest) {
    print(i)
    n.i = igraph::neighbors(ig, i)
    ew = EW[n.i]
    tb = table(ew)
    print(tb)
    if (any(names(tb) %in% c("W", "E"))) {
      if (("E" %in% names(tb)) & ("W" %in% names(tb))) {
        cat("cannot decide EW in rest (chose W): ", i,
            "\n")
        EW[i] = "W"
      }
      else EW[i] = names(tb[names(tb) %in% c("W", "E")])
    }
  }
  East = which(EW == "E")
  West = which(EW == "W")
  L = igraph::laplacian_matrix(ig, sparse = TRUE)
  B = rep(0, n.v)
  for (i in c(f1, f2)) {
    L[i, which(L[i, ] < 0)] <- 0
    L[i, i] <- 1
  }
  B[f2] = pi
  B[f1] = pi
  epsilon = pi/length(north.south)
  inner = setdiff(north.south, c(f1, f2))
  for (i in inner) {
    L[i, which(L[i, ] < 0)] <- 0
    L[i, i] <- 1
    B[i] = 0
  }
  for (i in West) {
    L[i, which(L[i, ] < 0)] <- 0
    L[i, i] <- 1
    B[i] = 2 * pi - epsilon
  }
  B[f2] = pi
  v <- solve(L, B)
  v = as.numeric(v)
  cat("RANGE RESID:")
  print(abs(range(L %*% v - B)))
  plot(u, v, main = "Brechbuehler initial")
  points(u[East], v[East], pch = 19, col = 2)
  points(u[West], v[West], pch = 19, col = 3)
  points(u[inner], v[inner], pch = 19, col = 4)
  X$texcoords = rbind(as.numeric(u), as.numeric(v))
  if (InitFit) {
    uv = t(X$texcoords)
    dim(uv)
    bas.i <- MakeBasis_UV(Fit_order, uv[, 1], uv[, 2])
    dim(bas.i$Ylm)
    a_v <- VertexAreasOBJ(X)
    A.init <- FitAlm_Tikhonov(x, bas.i, lambda = 0.1)
    grd = MakeGrid_GaussLegendreSimpson(50)
    bas = MakeBasis_UV(Fit_order, grd$U, grd$V)
    C <- updateX(A.init, grd, bas)
    rgl::plot3d(C$X, asp = F)
    range(A.init)
    rgl::open3d()
    rgl::wire3d(X)
    rgl::writeOBJ(file.out)
  }
  else A.init <- NULL
  return(list(uv = cbind(u, v), OBJ = X, L.full = L.ret, A.init = A,
              East = East, West = West, Inner = inner, f1 = f1, f2 = f2,
              igraph = ig, A.init = A.init, poles = p.fix,
              north.south=north.south))
}


#' Obj2ObjQ
#' @description
#' makes a 3d-object of quads from a 3d-object of triangles (requires regular grid grd)
#' @export
Obj2ObjQ<-function (O, grd)
{
  nx = grd$nu
  ny = grd$nv
  q = matrix(NA, 4, nx * ny)
  k = 0
  for (i in 1:(nx - 1)) for (j in 1:(ny - 1)) {
    k = k + 1
    l = (j - 1) * nx + i
    q[1, k] = l
    q[2, k] = l + 1
    q[3, k] = l + 1 + nx
    q[4, k] = l + nx
  }
  M <- rgl::qmesh3d(O$vb, indices = q[, 1:k], normals = matrix(0,
                                                               3, k))
  return(M)
}

#' X2ObjQ
#' @description
#' puts coordinates X into object O, which may be of quads
#' @export
X2ObjQ<-function(O,X)
{ O$vb=rbind(t(X),1) # ingest coordinates in 3d-graphics object
return(O)  # quads dont allow for vcgUpdatenormals for unknown reason
}

#' MakeMemRBC
#' @description
#' Make a MemRBC object from coefficients, grid and basis
#' @param A,grd,bas coeffs, grid and dasis data for MemRBC object
#' @return MemRBC object
#' @export
MakeMemRBC <- function(A,G,B)
{ return(structure(class="MemRBC",list(A=LM2A(A,B),grd=G,bas=B)))}

#' CenterX
#' @description
#' gives back centered coordinates X
#'  but not center of mass, if X contains doubles
#'  
#' @export
CenterX<-function(X)
{  return(apply(X,2,function(x) x-mean(x))) }

#' TriMesh_Unduloid
#' @description
#' create an Unduloid 3d object for a fraction or multiple periods.
#' Unduloids may be interesting shapes to fit, see example.
#' @examples
#' # a special grid is made with hole at north and south pole,
#' # ie, u starts not at zero
#' open=0.4 # also higher work, but objects boundary is never fitted
#' g<-MakeGrid_GaussLegendreSimpson(180,ua=open,ub=pi-open)
#' g$ndof
#' range(g$U)
#' (U<-TriMesh_Unduloid(periods=4,nx=g$nu,ny=g$nv-1, a=1,c=0.1,clean=FALSE))
#' attr(U,"H_theor") # theoretical mean curvature
#' attr(U,"H_vcg_6") # mean from vertices with 6 neighbors
#' #  mean curvature for other vertices is problematic in vcg
#' b<-MakeBasis_UV(23,g$U,g$V)
#' # one should exclude double coordinates for
#' # the fit, so mask is needed
#' b$mask<-double_uv_ind(b$uv[,1],b$uv[,2])
#' CenterX(Obj2X(U)) -> X
#' cat(dim(X)[1],"?=", g$ndof,"\n")
#' if (dim(X)[1] == g$ndof)
#' { # if not matching repeat with alternative n in grid
#' rgl::plot3d(X,col=2,aspect=FALSE)
#' X2Obj(U,X) -> U1
#' rgl::contourLines3d(U1,b$uv[,1],levels=(0:100)*pi/100)
#' rgl::contourLines3d(U1,b$uv[,2],40)
#' rgl::shade3d(U1,col="grey",alpha=0.5)
#' A<-FitAlm_Tikhonov(X,b,lambda=0) # , WX=sin(g$U))
#' A[,3]<--A[,3] # wrong orientation correction
#' # you may try the fit with weights, WX=sin(g$U)
#' MakeMemRBC(A,g,b)->M
#' rgl::open3d()
#' rgl::plot3d(X,aspect=FALSE,alpha=0.45)
#' plot(M,alpha=0.45,col="cyan",wire=FALSE)
#' E=E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))
#' -E$Curv/E$Area/2 # mean curvature from integral over area
#' attr(U,"H_vcg_6") # comparison with original
#' imag.obj.colorbar(U1,E$curv)
#' rgl::title3d("curvature density")
#' X1<-updateX(M$A,M$grd,M$bas)$X
#' M$grd$Obj<-X2Obj(M$grd$Obj,X1)
#' rgl::shade3d(M$grd$Obj,alpha=0.2)
#' mean(E$curv/E$dA/2)
#' MemRBC_env$M.C0<-0;MemRBC_env$M.mu<-0;MemRBC_env$M.Ka<-0
#' M$bas$Target[1:2]=c(E$Area,E$Volume)
#' save_MemRBC(M,"Unduloid.rdat")
#' MMC(M,1000,plt=TRUE,pltfreq=2,C0=0)
#' }
#' @export
TriMesh_Unduloid<-function(a=1,c=2,periods=1.0,nx=40,ny=40,shade=TRUE,wire=FALSE,clean=TRUE)
{ m=(c^2-a^2)/2;n=(c^2+a^2)/2
  mu=2/(a+c);k2=(c^2-a^2)/c^2
  ulimup=(pi/2+pi/4)*2/mu; ulimdown=pi/4*2/mu
  p=ulimup+ulimdown
  u=pracma::linspace(-ulimdown,-ulimdown+p*periods,nx); phi=u*mu/2-pi/4
  x=Re(a*Carlson::elliptic_F(phi,k2)+c*Carlson::elliptic_E(phi,k2))
  z=Re(sqrt(m*sin(mu*u)+n))
  uv=matrix(NA,ny*nx,2)
  d=2*pi/ny
  r_mat=matrix(c(cos(d),sin(d),-sin(d),cos(d)),2,2)

  uv[,2]=rep(pracma::linspace(0,2*pi,ny),nx)
  uv[,1]=rep(pracma::linspace(0,pi,nx),each=ny)

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
  if (clean) Rvcg::vcgClean(M,1:7)->M
  attr(M,"H_theor")=1/(a+c) # store theoretical value of constant mean curvature
  tb=table(M$it) # trick to get degree 6 vertices:
  attr(M,"H_vcg_6")=mean(Rvcg::vcgCurve(M)$meanvb[tb==6])
  # the attribute is probably not inherited in derived objects.
  M$uv=uv
  return(M)
}

#
# use polynomials of degree 12 to fit to lowess fits with other data
#
#' @export
Lowess_vcg_meanvbOBJ<-function(x,O) # return mean curvature
{
  LW=list()
  LM=list()
  crv=Rvcg::vcgCurve(O)
  tb=table(O$it)# how often a vertex is addressed
  kr=unique(tb)
  print(kr)
  plot(0,0,col=0,xlim=range(x),ylim=range(-crv$meanvb))
  for (k in kr){
    length(x[tb==k])
    length(crv$meanvb[tb==k])

    LW[[k]] = lowess(x[tb==k],-crv$meanvb[tb==k],f=0.1)
    points(LW[[k]],type="l",lwd=3,col=k)
    xx=LW[[k]]$x;yy=LW[[k]]$y
    LM[[k]] = pracma::polyfit(xx,yy,12)
    points(xx,pracma::polyval(LM[[k]],xx),col=k,cex=1.5,lwd=2)
  }
  X=Y=Z=rep(0,length(x)) # unordered return
  LX=LY=LK=list()
  for (k in kr){
    if (k!=6) LY[[k]]=-crv$meanvb[tb==k] - pracma::polyval(LM[[k]],x[tb==k]) + pracma::polyval(LM[[6]],x[tb==k]) else LY[[6]]=-crv$meanvb[tb==6]
    points(x[tb==k],LY[[k]],col=k)
    Y[tb==k]=LY[[k]]
    LX[[k]]=x[tb==k]
    X[tb==k]=LX[[k]]
    LK[[k]]=rep(k,sum(tb==k))
    Z[tb==k]=k
  }
  return(list(x=unlist(LX),y=unlist(LY),k=unlist(LK),X=X,Y=Y,K=Z))
  # k-sorted output, not good for coloring 3d object
}


#
# dense regions in u,v can be stretched by this
#  - coserves triangulation quality
#   needed for postprocessing Brechbühler
#
#' spreadout.uv
#' @description Brechbuhler initial uv may be redistributed
#' Sparse regions are contracted.
#' @param uv : n x 2 matrix of (u,v), like $uv returned from Brechbuhler
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

#' ConsIter
#' @description
#' Iterates coefficients to fulfill constraints.
#' Needed in Rosen Projection methods.
#' only SCM gradient is computed for constraint Jacobian
#' @export
ConsIter<-function(A,grd,bas,C,g2, Ctol=1e-3, nsteps=20,
                   prn=FALSE,del_cons=0.3,nm=FALSE,do_one=TRUE)
{
  l=0; Nc=bas$Nc
  updateX(A,grd,bas)->C
  h2=E_SCM(A,grd,bas,C)
  Cons_RHS <-ConsRHS(h2,bas)
  NCons<-max(abs(Cons_RHS))
  if(prn) cat(crayon::yellow(l,": Cons_%:"), crayon::cyan(round(100*Cons_RHS/bas$Target,4)),crayon::yellow(" |Cons|:"),ifelse(NCons>Ctol,crayon::red(NCons),crayon::green(NCons)),"\n")

  if (NCons<=Ctol) {
    E_SCM(A,grd,bas,C) -> h2
    Grad_SCM(h2,grd,bas,C) -> g2

    Cons_RHS <-ConsRHS(h2,bas)
    NCons<-max(abs(Cons_RHS))
  };
  sol=rep(0,Nc) # default to return
  while ((NCons>Ctol & l<nsteps) | ( l==0 )){
    #if (l==0 & !do_one) break; # not implemented;
    NCons1=NCons
    dFm<-c(g2[[ bas$Cons[1] ]])
    for (ii in 2:Nc) dFm=cbind(dFm,c(g2[[bas$Cons[ii]]])) # additional constraints
    M_c<-matrix(0.0,Nc,Nc); for (i in 1:Nc) for (j in 1:Nc) M_c[i,j]<-.dot2(dFm[,i],dFm[,j])
    Pm<-pracma::pinv(M_c)
    sol<- (- Pm %*%Cons_RHS)[,1]
    names(sol)=names(bas$Cons)
    delta <- sol[1]*dFm[,1] ; for (i in 2:Nc) delta<-delta + sol[i]*dFm[,i]
    A <- A + del_cons * delta

    updateX(A,grd,bas) -> C
    E_SCM(A,grd,bas,C) -> h2
    Grad_SCM(h2,grd,bas,C) -> g2

    Cons_RHS <-ConsRHS(h2,bas)
    NCons<-max(abs(Cons_RHS))
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

#' FullModelHessian
#' @description
#' computes the full model Hessian of the energy, as used in CNM.
#' It is based on finite differences of gradients
#' with a symmetrization; includes full Lagrangian.
#' @param A, grd, bas, Ref : standard objects of a MemRBC
#' @param del (=1e-6) finite difference delta for coefficients
#' @returns Full Hessian matrix, ie also constraints, but not bordered
#' @export
FullModelHessian<-function (A, grd, bas, Ref, del = 1e-06, Ctol = 0.001)
{
  L0 = L = list()
  tictoc::tic()
  Ai_max = bas$Ai_max
  C = updateX(A, grd, bas)
  h20 = E_SCM(A, grd, bas, C)
  S = SEN(A, grd, bas, Ref, h20)
  ES = E_SEN(A, grd, bas, S, Ref)
  Gh20 = Grad_SCM(h20, grd, bas, C)
  GS0 = Grad_SEN(A, grd, bas, Gh20, S, Ref)
  G0 = c(Gh20$grad_SCM + GS0$grad_SEN)
  for (i in bas$Cons[1:bas$Nc]) L[[i]] = matrix(0, Ai_max *
                                                  3, Ai_max * 3)
  H = matrix(0, Ai_max * 3, Ai_max * 3)
  j = 0
  A1 = A
  for (k in 1:3) for (i in 1:Ai_max) {
    j = j + 1
    if (j%%10 == 0)
      cat(" ", round(j/(3 * Ai_max) * 100, 1), "\r")
    C = synth_update(C, bas, i, k, +del)
    h2 = E_SCM(A, grd, bas, C)
    Gh2 = Grad_SCM(h2, grd, bas, C)
    C = synth_update(C, bas, i, k, -del)
    A1[j] = A1[j] + del
    S = SEN(A1, grd, bas, Ref, h2)
    GS = Grad_SEN(A1, grd, bas, Gh2, S, Ref)
    A1[j] = A1[j] - del
    Gj = c(Gh2$grad_SCM + GS$grad_SEN) - G0
    for (m in bas$Cons[1:bas$Nc]) {
      L[[m]][, j] = (Gh2[[m]] - Gh20[[m]])/del
    }
    H[, j] = Gj/del
  }
  for (m in bas$Cons[1:bas$Nc]) L[[m]] = 0.5 * (L[[m]] + t(L[[m]]))
  tictoc::toc()
  return(list(H = (H + t(H))/2, L = L, H_fd = H, G = G0, GS = GS,
              g2 = Gh20, E = h2$Wb + ES, ES = ES, Wb = h2$Wb, h2 = h20,
              C = C, gradC = Gh20$gradC, gradA = Gh20$gradA, gradV = Gh20$gradV,
              A = A))
}


#' ID (not used anymore)
#' @description a template filter for modification of
#' coefficients.
#' Such a filter can be used in only a few apps, like MMC.
#' @export
ID<-function(A,bas)
{ return(A)
}

#' FullModelHessian_Par
#' @description parallel computed Hessian including constraint gradients; serial gradients but a list-parallel approach to assemble H
#' @param cl : parallel cluster-ID from previous result (H$cl in CNM)
#' @export
FullModelHessian_Par <- function(A, grd, bas, Ref, del = 5e-06, 
                                 Mem_mc.cores = 4, timing = TRUE, 
                                 stopdown = TRUE, cl=NULL)
{ pt0 = proc.time()
  startup=is.null(cl)
  Ai_max = bas$Ai_max
  C = updateX(A, grd, bas)
  h20 = E_SCM(A, grd, bas, C)
  S = SEN(A, grd, bas, Ref, h20)
  ES = E_SEN(A, grd, bas, S, Ref)
  Gh20 = Grad_SCM(h20, grd, bas, C)
  GS0 = Grad_SEN(A, grd, bas, Gh20, S, Ref)
  G0 = c(Gh20$grad_SCM + GS0$grad_SEN)
  L = L0 = list()
  for (m in bas$Cons[1:bas$Nc]) 
    L[[m]] = L0[[m]] = matrix(0, Ai_max * 3, Ai_max * 3)
  for (m in bas$Cons[1:bas$Nc]) for (j in 1:(3 * Ai_max)) 
    L0[[m]][,j] = Gh20[[m]]
  H = matrix(0, Ai_max * 3, Ai_max * 3)
  Lpar = list()
  L = list()
  for (i in bas$Cons[1:bas$Nc]) L[[i]] = matrix(0, Ai_max * 3, Ai_max * 3)
  if (timing)
    tictoc::tic()
  A1 = A
  cat("Par: ENV C0",MemRBC_env$M.C0,"\n")
  for (j in 1:(3 * Ai_max)) {
    A1[j] = A1[j] + del
    Lpar[[j]] = list(A = A1, 
                     M.C0 = MemRBC_env$M.C0, 
                     M.K_ADE = MemRBC_env$M.K_ADE,
                     M.K_b = MemRBC_env$M.K_b, 
                      M.mu = MemRBC_env$M.mu, 
                      M.Ka = MemRBC_env$M.Ka, 
                      M.a3 = MemRBC_env$M.a3,
                      M.a4 = MemRBC_env$M.a4, 
                      M.b1 = MemRBC_env$M.b1, 
                      M.b2 = MemRBC_env$M.b2, 
                    M.Rcpp = TRUE,
             M.Rcpp_ncores = MemRBC_env$M.Rcpp_ncores, 
                     index = j)
    A1[j] = A1[j] - del
  }
  if (timing) {
    cat("paralleliz. preperation  ")
    tictoc::toc()
    tictoc::tic()
  }
  { if (is.null(cl)) {
      cat("setup Cluster on ",Mem_mc.cores," cores\n")
      cl <- parallel::makeCluster(Mem_mc.cores,type="PSOCK", outfile = "tmp_cluster.txt")
  }  
    # variables to ship to cluster nodes:
    M.C0 = MemRBC_env$M.C0 
    M.K_ADE = MemRBC_env$M.K_ADE
    M.K_b = MemRBC_env$M.K_b
    M.mu = MemRBC_env$M.mu
    M.Ka = MemRBC_env$M.Ka 
    M.a3 = MemRBC_env$M.a3
    M.a4 = MemRBC_env$M.a4 
    M.b1 = MemRBC_env$M.b1 
    M.b2 = MemRBC_env$M.b2 
    M.Rcpp = TRUE
    M.Rcpp_ncores = MemRBC_env$M.Rcpp_ncores
    M.Ref=Ref
    M.bas=bas
    M.grd=grd
    a=pi
#    parallel::clusterEvalQ(cl,{print("Node lives")})
#    parallel:::checkCluster(cl)
#    print("Cluster checked")
    parallel::clusterExport(cl, varlist=
                          c("a","M.bas", "M.grd","M.Ref",
                            "M.C0", "M.K_ADE",
                            "M.K_b", "M.mu", "M.Ka", 
                            "M.a3", "M.a4", "M.b1",
                            "M.b2", "M.Rcpp", "M.Rcpp_ncores",
                            "int2d_matrix_cxx", 
                            "vectomat_cxx", 
                            "vectoarr_cxx"
                           ),envir = environment())
#    parallel::clusterEvalQ(cl,{library(MemRBC);print(MemRBC_env$M.C0)})
    if (timing) {
      cat("cluster startup ");  tictoc::toc()
      tictoc::tic()
    } #parLapply
    LH <- parallel::parLapply(cl, Lpar, FullHessian_Client)
    if (stopdown) {
      parallel::stopCluster(cl)
      cl=NULL # returns in H that cluster is dead
    }
    }
  if (timing) {
    cat("parallel lapply for ", length(Lpar), " calls:")
    tictoc::toc()
  }
  pt = rep(0, 5)
  for (i in 1:(3 * Ai_max)) {
    H[, LH[[i]]$index] = c(LH[[i]]$G - G0)/del
    pt = pt + LH[[i]]$time
  }
  for (m in bas$Cons[1:bas$Nc]) for (i in 1:(3 * Ai_max)) {
    L[[m]][, LH[[i]]$index] = LH[[i]]$Lj[[m]]
  }
  for (m in bas$Cons[1:bas$Nc]) L[[m]] = (L[[m]] - L0[[m]])/del
  for (m in bas$Cons[1:bas$Nc]) L[[m]] = 0.5 * (L[[m]] + t(L[[m]]))
  pt0 = proc.time() - pt0
  return(list(H = (H + t(H))/2, L = L, H_fd = H, G = G0, g2 = Gh20,
              E = h20$Wb + ES, ES = ES, Wb = h20$Wb, h2 = h20, C = C,
              gradC = Gh20$gradC, gradA = Gh20$gradA, gradV = Gh20$gradV,
              A = A, proc_time_clients = pt, proc_time_total = pt0,
              Nclients = Mem_mc.cores, Par = pt0/pt, cl=cl ))
}

# export - but internal to FullHessian_Par
#' FullHessian_Client
#' @description
#' code exported for your inspection/replacement
#' @export
FullHessian_Client<-function(L,DBG=FALSE) # Lmax=13 takes 4 seconds per call/core
{ pt=proc.time()
  cat("CLIENT C0 before ENV",MemRBC_env$M.C0,"\t",L$M.C0,"\n")
  
  MemRBC_env$M.Rcpp<-L$M.Rcpp;
  MemRBC_env$M.Rcpp_ncores<-L$M.Rcpp_ncores
  MemRBC_env$M.C0<-L$M.C0
  MemRBC_env$M.K_ADE<-L$M.K_ADE 
  MemRBC_env$M.K_b<-L$M.K_b 
  MemRBC_env$M.mu<-L$M.mu  
  MemRBC_env$M.Ka<-L$M.Ka  
  MemRBC_env$M.a3<-L$M.a3 
  MemRBC_env$M.a4<-L$M.a4 
  MemRBC_env$M.b1<-L$M.b1
  MemRBC_env$M.b2<-L$M.b2
  Ai_max=M.bas$Ai_max
  
  cat("CLIENT C0 from ENV",MemRBC_env$M.C0,"\n")
  C=updateX(L$A,M.grd,M.bas) # no longer in L
  h2=E_SCM(L$A,M.grd,M.bas,C)
  S=SEN(L$A,M.grd,M.bas,M.Ref,h2)

  # decide by M.Rcpp for openmp-parallel code
  #       cant directly use Grad_SCM on cluster due to scattered objects names
  if(!MemRBC_env$M.Rcpp) {Gh2=Grad_SCM_R(h2,M.grd,M.bas,C)} else {
    G2=Grad_SCM_cxx(h2, M.grd, M.bas, C, L$M.C0, L$M.Rcpp_ncores, L$M.K_b, L$M.K_ADE)
    Gh2=list(ddA=array(G2$ddA,c(M.grd$ndof,M.bas$Ai_max,3)),
             ddV=array(G2$ddV,c(M.grd$ndof,M.bas$Ai_max,3)),
             grad_SCM=G2$grad_SCM,gradV=G2$gradV,gradA=G2$gradA,gradC=G2$gradC,
             dE=array(G2$dE,c(M.grd$ndof,M.bas$Ai_max,3)),
             dF=array(G2$dF,c(M.grd$ndof,M.bas$Ai_max,3)),
             dG=array(G2$dG,c(M.grd$ndof,M.bas$Ai_max,3)))
    }
  GS=Grad_SEN(L$A,M.grd,M.bas,Gh2, S, M.Ref)
  Gj=c(Gh2$grad_SCM + GS$grad_SEN)
  Lj=list()
  for (k in M.bas$Cons[1:M.bas$Nc])  Lj[[k]] = G2[[k]]

  cat("client done with",L$index," at C0=",MemRBC_env$M.C0," in ",proc.time()-pt,"\n")
  return(list(G=Gj,Lj=Lj,time=proc.time()-pt,index=L$index))
}

# unlockBinding("M.C0", env = MemRBC_env)

#'mat2vec
#' @export
mat2vec<-function(m) return(c(m))
#' @export
matdiff2vec<-function(m1,m2)
{return(c(m1-m2))}

#'symmetrize
#' @export
symmetrize<-function(m)
{return((m+t(m))/2.0)}

#'matadd2vec
#' @export
matadd2vec<-function(m1,m2)
{return(c(m1+m2))}

#' SetConstraints
#' @description
#' set the constraint target values and store information in basis.
#' Without parameters Cons, ..., the standard area and volume constraints are set with values (140, 100).
#' If you set a third constraint, for CNM you have to give M$Lambda a third component.
#' Default constraints are 140 area, 100 volume
#' @param bas basis to modify constraints, e.g. M$bas (M MemRBC object)
#' @param Cons vector of character of constraint gradients, from "gradA","gradV" and "gradC"
#' @param QCons vector of constraint names, from "Area","Volume","Curv", same order as Cons
#' @param Target vector of constraint values, like c(140,100,88) for Area, Volume and Curv
#' @examples
#' data("M4",package = "MemRBC"); SetParams(M4)
#'
#' # add curvature constraint:
#' SetConstraints(M4$bas, Cons=c("gradA","gradV","gradC"),
#'   QCons=c("Area","Volume","Curv"),
#'   Target=c(140,100,121) ) -> M4$bas  # store modified basis back into membranes D5 basis
#'
#' # minimize with steepest descend under Rosen Constraint Projection
#' \donttest{
#' SDRC(M4,100)->M4sdrc
#' # pair-plots of target quantities and energy E
#' plot(M4sdrc$SDRC_Sample[c("E","A","V","C")])
#' }
#'@param bas: either basis like M$bas or MemRBC-object M
#'@return MemRBC, if bas is a MemRBC-object, or updated basis 
#' @export
SetConstraints<-function(bas,Cons=c("gradA","gradV"),
                         QCons=c("Area","Volume"),
                         Target=c(140, 100)) # any further (implemented) gradients allowed
{
  if(is(bas,"MemRBC")) {M=bas;bas=M$bas;toMemRBC=TRUE} else toMemRBC=FALSE
  bas$Cons=Cons
  bas$Nc=length(Cons)
  bas$Target=Target
  bas$QCons=QCons
 
  names(bas$Cons)= QCons
  names(bas$QCons)= QCons

  names(bas$Target)= QCons
  message("SetConstraints: for CNM, remember to set M$Lambda with bas$Nc, eg M$Lambda=rep(0.1,M$bas$Nc) \n")
  if (!toMemRBC) return(bas) else {M$bas<-bas; return(M)}
}

#' ConstraintHessian
#' @description 
#' constraint (=bordered) Hessian for CNM() Newton minimizer
#' @export
ConstraintHessian<-function(H,bas,Lambda,filter=ID)
{
 dH=dim(H$H)[1] # should be equal Ai_max*3
 Nc=bas$Nc
 if(Nc==0) {message("No constraints for ConstraintHessian - return Hessian as is\n");return(H)}
 if(length(Lambda)!=Nc) stop("wrong number of Lagrangian lambdas vs. registered constraints bas$Nc")
 M=matrix(0,dH+Nc,dH+Nc)
 M[1:dH,1:dH]=H$H # main block
 for (i in 1:bas$Nc)  M[dH+i,1:dH]=M[1:dH,dH+i]=filter(H[[ bas$Cons[i] ]],bas)
 return(M)
}

#' ConsJacobian
#' @description constraint Jacobian for Rosen projection
#' @export
ConsJacobian<-function(g2,bas)
{
  RosenA<-c(g2[[bas$Cons[1]]]) # first gradient, bas$Cons has names of gradients
  for (ii in 2:bas$Nc) # further gradients to add
  RosenA<-rbind(RosenA,c(g2[[ bas$Cons[ii] ]] ))
  return(RosenA)
}

#' RosenProjection
#' @description Rosen (1961) projection of gradient
#' @param G input gradient
#' @param g2 current energy
#' @param bas basis
#' @export
RosenProjection<-function(G,g2,bas) # output the projected energy gradient G; g2 contains needed gradients of constraint functions
{
  RosenA<-ConsJacobian(g2,bas)
  RosenAAt<-RosenA%*%t(RosenA)
  IRosenAAt<-pracma::pinv(RosenAAt) # pinv tolerates also gradients to be zero
  lambda<- t(RosenA)%*%IRosenAAt%*%RosenA
  Gprime <- c(G) - (lambda %*% c(G)) [,1]
  lambdaG=c()
  for (i in 1:bas$Nc) lambdaG[i]=.dot2(c(G),RosenA[i,]) # how much from each constraint gradient  is projected out of Gradient?
  names(lambdaG)=bas$Qcons
  return(list(Gprime=G,Eigs=eigen(RosenAAt)$values,lambdaG=lambdaG))
}

#' ConsRHS
#' @description
#' Returns constraint violation values, e.g. for rhs of Rosen solver
#' @param h2 current energy from E_SCM
#' @param bas basis
#' @returns a named vector of constrained violations
#' @export
ConsRHS<-function(h2,bas)
{  Cons_RHS <- rep(0,bas$Nc);names(Cons_RHS)=bas$QCons
for (k in 1:bas$Nc) Cons_RHS[k] <- h2[[bas$QCons[k]]] - bas$Target[k]
return(Cons_RHS)
}


# if you decide to erase dot2 from Rcpp-code:
if(!exists(".dot2")) .dot2=function(x,y) sum(x*y)


#' Cons_filter_delta
#' @description remove constraint-violating components from coefficient step delta
#' @param H2 current E_SCM
#' @returns projected delta as nx3 matrix
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

#' make_delta_normal_to_surface
#' @description
#' normal motion filter by integrals
#' (contributed by C Woelper, 2024, University Bremen)
#' @param delta step in coefficients to project
#' @param grd,bas grid and basis
#' @param n current spatial normals, e.g. returned as $n in E_SCM
#' @export
  make_delta_normal_to_surface <- function(delta, grd, bas, n){
    result <- matrix(0,bas$Ai_max,3)
    d=matrix(delta,bas$Ai_max,3)
    delta_pointwise <- synthX(bas$Ylm,d)
    incomplete_integrand <- n * pracma::dot(t(n), t(delta_pointwise)) * sin(grd$U)
    for( i in 1:bas$Ai_max){
      Yi <- bas$Ylm[,i]
      integrand <- incomplete_integrand * Yi
      for (k in 1:3){
        result[i,k] <- int2d_scalar_GLS(integrand[,k], grd)
      }
    }
    return(result[])
  }

#' E_FullModel_Penalty_AV
#' @description
#' Area-Volume-constraint energy computation
#'  works also without SEN (Ref=NULL)
#' @export
  E_FullModel_Penalty_AV<-function (A, grd, bas, Ref)
  {
    C <- updateX(A, grd, bas)
    h2 <- E_SCM(A, grd, bas, C)
    if (!is.null(Ref)) {
      S <- SEN(A, grd, bas, Ref, h2)
      e <- E_SEN(A, grd, bas, S, Ref)
    }
    else {
      e = 0
      S = NULL
    }
    E <- h2$Wb + e + MemRBC_env$M.rho * ((h2$Volume - bas$Target["Volume"])^2 +
                                (h2$Area - bas$Target["Area"])^2)
    names(E) = NULL
    return(list(E = E, Wb = h2$Wb, Ws = e, E_uncons = h2$Wb +
                  e, dA = h2$dA, S = S, Area = h2$Area, Volume = h2$Volume,
                Curv = h2$Curv, n = h2$n))
  }
#' Grad_FullModel_Penalty_AV
#' @export
Grad_FullModel_Penalty_AV<-
  function (A, grd, bas, Ref, S)
  {
    C <- updateX(A, grd, bas)
    h2 <- E_SCM(A, grd, bas, C)
    G_SCM <- Grad_SCM(h2, grd, bas, C)
    if (!is.null(Ref)) {
      G_SEN <- Grad_SEN(A, grd, bas, G_SCM, S, Ref)
      G <- G_SCM$grad_SCM + G_SEN$grad_SEN + 2 * MemRBC_env$M.rho * (G_SCM$gradV *
                                                            (h2$Volume - bas$Target["Volume"]) + +G_SCM$gradA *
                                                            (h2$Area - bas$Target["Area"]))
    }
    else G <- G_SCM$grad_SCM + 2 * MemRBC_env$M.rho * (G_SCM$gradV * (h2$Volume -
                                                             bas$Target["Volume"]) + +G_SCM$gradA * (h2$Area - bas$Target["Area"]))
    return(G)
  }

#' E_FullModel_Penalty_AVC
#' @description
#' Area-Volume-Curvature-constraint energy computation
#' @export
E_FullModel_Penalty_AVC<-function (A, grd, bas, Ref)
  {
    C <- updateX(A, grd, bas)
    h2 <- E_SCM(A, grd, bas, C)
    if (!is.null(Ref)) {
      S <- SEN(A, grd, bas, Ref, h2)
      e <- E_SEN(A, grd, bas, S, Ref)
    }
    else {
      e = 0
      S = NULL
    }
    E <- h2$Wb + e + MemRBC_env$M.rho * ((h2$Volume - bas$Target["Volume"])^2 +
                                (h2$Area - bas$Target["Area"])^2 + (h2$Curv - bas$Target["Curv"])^2)
    names(E) = NULL
    return(list(E = E, Wb = h2$Wb, Ws = e, E_uncons = h2$Wb +
                  e, dA = h2$dA, S = S, Area = h2$Area, Volume = h2$Volume,
                Curv = h2$Curv,n=h2$n))
  }


#' Grad_FullModel_Penalty_AVC
#' @export
Grad_FullModel_Penalty_AVC<-function (A, grd, bas, Ref, S)
{
  C <- updateX(A, grd, bas)
  h2 <- E_SCM(A, grd, bas, C)
  G_SCM <- Grad_SCM(h2, grd, bas, C)
  if (!is.null(Ref)) {
    G_SEN <- Grad_SEN(A, grd, bas, G_SCM, S, Ref)
    G <- G_SCM$grad_SCM + G_SEN$grad_SEN + 2 * MemRBC_env$M.rho * (G_SCM$gradV *
                                                          (h2$Volume - bas$Target["Volume"]) + +G_SCM$gradA *
                                                          (h2$Area - bas$Target["Area"]) + +G_SCM$gradC * (h2$Curv -
                                                                                                             bas$Target["Curv"]))
  }
  else G <- G_SCM$grad_SCM + 2 * MemRBC_env$M.rho * (G_SCM$gradV * (h2$Volume -
                                                           bas$Target["Volume"]) + +G_SCM$gradA * (h2$Area - bas$Target["Area"]) +
                                            +G_SCM$gradC * (h2$Curv - bas$Target["Curv"]))
  return(G)
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

#' history_MemRBC
#' @export
history_MemRBC<-function(M)
{
  print(paste(M$history))
}

#' MakeIM
#' @export
MakeIM<-function (bas, lambda = 0, WX = rep(1, dim(bas$Ylm)[1]))
{
  mask = bas$mask
  Y = cbind(1, bas$Ylm[-mask, ])
  YtW = t(Y) %*% diag(WX[-mask])
  B = YtW %*% Y
  InvB1 = pracma::inv(B + lambda * diag(c(0, bas$G.tk)))
  IM = InvB1 %*% YtW
  bas$IM <- IM
  return(bas)
}

#' double_uv_ind
#' @export
double_uv_ind<-function (u, v, digits = 7)
{
  l = length(u)
  p = round(cbind(sin(u), cos(u), sin(v), cos(v)), digits)
  pc = apply(p, 1, paste, sep = "", collapse = "")
  pt = table(pc)
  dbl = names(pt)[which(pt == 2)]
  d = rep(NA, length(dbl))
  for (i in 1:length(dbl)) d[i] = which(pc == dbl[i])[2]
  return(d)
}

#' FitFast
#' @description
#' fits coefficients to X from least squares solver matrix IM in bas
#' This is used for fitting e.g. displacements
#' like in pertA_complex.
#' @param bas,X basis with IM, 3d-coordinates to fit
#' @export
FitFast<-function (bas, X)
{
  if (is.null(bas$IM))
    stop("first create bas$IM with MakeIM(bas)->bas")
  mask = bas$mask
  return((bas$IM %*% X[-mask, ])[-1, ])
}

#' Movie
#' @description
#' show shapes along recorded coefficients in M$LA
#' @export
Movie<-function (M,x=400,y=400,from=1,to=length(M$LA),sleep=0.025,phi=pi/100,phi1=phi/2,...)
{
  ii1=i1=rgl::open3d(); rgl::par3d(dev=i1,windowRect=c(10,10,x+10,y+10))
  rgl::points3d(0,col="white")
  rgl::rgl.bringtotop()
  ii2=i2=rgl::open3d()
  rgl::points3d(0,col="white");rgl::par3d(dev=i2,windowRect=c(10,10,x+10,y+10))
  k=from
  Sys.sleep(sleep*3)
  R1=matrix(c(1,0,0, 0, cos(phi1),sin(phi1),  0,-sin(phi1),cos(phi1) ),3,3)
  R=matrix(c(cos(phi),0,sin(phi),  0,1,0, -sin(phi),0,cos(phi) ),3,3)
  for (A in M$LA[from:to]) {
    i3=i2;i2=i1;i1=i3;
    p=rgl::par3d();
    rgl::set3d(i1) # switch to background for updating
    p$userMatrix[1:3,1:3]<-R1%*%R%*%p$userMatrix[1:3,1:3]
    rgl::par3d( modelMatrix=p$modelMatrix,scale=p$scale,
                observer=p$observer,
                userMatrix=p$userMatrix,
                projMatrix=p$projMatrix,
                windowRect=p$windowRect,
                viewport=p$viewport,
                scale=p$scale,zoom=p$zoom ) # copy parms

    M$A <- A
    i = rgl::ids3d()
    rgl::pop3d(id = i$id[i$type %in% c("triangles", "quads","text")])
    plot(M, wire_col = "red",...)

    rgl::title3d(k);k=k+1
    rgl::rgl.bringtotop() # display background

    Sys.sleep(sleep)
  }
}

#' MemCols
#' @export
MemCols<-function(data,pal=rainbow,n=100)
{return(pal(n)[1+(n-1)*(data-min(data))/diff(range(data))])}

#' DistBasedWX
#' @description
#' computes spatial weights WX from 8 neares neighbors.
#' The weights may be used for subsampling or weighted fits.
#' Background: too dense gridpoints lead to underfitting in less denser regions.
#'
#' @export
DistBasedWX<-function(M)
{ X=updateX_only(M$A,M$grd,M$bas)$X
  NN=RANN::nn2(X,9)
  d=apply(NN$nn.dists[,-1],1,mean)
  WX=d/sum(d)
  return(WX)
}

#' @export
color2d<-function(col1=c(0,0,0),col2=c(1,0,0),col3=c(0,0,1), col4=c(0,1,0),n=25){
  RGB=matrix("",n,n)
  plot(0,xlim=c(1,n),ylim=c(1,n))
  f<-function(x,y) col1*(1-x)*(1-y)+col2*(1-x)*y+col3*x*(1-y)+col4*x*y
  x=pracma::linspace(0,1,n=n)
  for (i in 1:n) for (j in 1:n) {C=f(x[i],x[j]);
  RGB[i,j]=rgb(red=C[1],green=C[2],blue=C[3])
  points(i,j,col=RGB[i,j],pch=15,cex=4)}
  return(RGB)
}

# remove any App-cpmputational results from a memRBC object.
# The object is reduced to a minimum of data for further work.
# Only grd,bas,Ref,A, Params and flag SEN are returned.
#' @export
Empty<-function(M)
{
  for (i in names(M)[-which(names(M) %in% c("grd","bas","Ref","A","Params","SEN"))])
    M[[i]]=NULL
  return(M)
}

# rewind data and A in M by some frames from LA, and erase LA later than that
#' @export
Rewind<-function(M,last=2)
{ M$A=M$LA[[length(M$LA)-2]]
  M$LA[(length(M$LA)-last):length(M$LA)] <- NULL
  return(M)
}

own.imagePlot<-function (..., add = FALSE, breaks = NULL, nlevel = 64, col = NULL, 
                     horizontal = FALSE, legend.shrink = 0.9, legend.width = 1.2, 
                     legend.mar = ifelse(horizontal, 3.1, 5.1), legend.lab = NULL, 
                     legend.line = 2, graphics.reset = FALSE, bigplot = NULL, 
                     smallplot = NULL, legend.only = FALSE, lab.breaks = NULL, 
                     axis.args = NULL, legend.args = NULL, legend.cex = 1, midpoint = FALSE, 
                     border = NA, lwd = 1, verbose = FALSE) 
{ pn=(rgl::cur3d()>0)
  old.par <- par(no.readonly = TRUE)
  if (is.null(col)) {
    col <- tim.colors(nlevel)
  }
  else {
    nlevel <- length(col)
  }
  info <- fields::imagePlotInfo(..., breaks = breaks, nlevel = nlevel)
  breaks <- info$breaks
  if (verbose) {
    print(info)
  }
  if (add) {
    big.plot <- old.par$plt
  }
  if (legend.only) {
    graphics.reset <- TRUE
  }
  if (is.null(legend.mar)) {
    legend.mar <- ifelse(horizontal, 3.1, 5.1)
  }
  temp <- fields::imageplot.setup(add = add, legend.shrink = legend.shrink, 
                          legend.width = legend.width, legend.mar = legend.mar, 
                          horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
  smallplot <- temp$smallplot
  bigplot <- temp$bigplot
  if (!legend.only) {
    if (!add) {
      par(plt = bigplot)
    }
    if (!info$poly.grid) {
      image(..., breaks = breaks, add = add, col = col)
    }
    else {
      poly.image(..., add = add, breaks = breaks, col = col, 
                 midpoint = midpoint, border = border, lwd.poly = lwd)
    }
    big.par <- par(no.readonly = TRUE)
  }
  if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
    par(old.par)
    stop("plot region too small to add legend\n")
  }
  ix <- 1:2
  iy <- breaks
  nBreaks <- length(breaks)
  midpoints <- (breaks[1:(nBreaks - 1)] + breaks[2:nBreaks])/2
  iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints))
  if (verbose) {
    print(breaks)
    print(midpoints)
    print(ix)
    print(iy)
    print(iz)
    print(col)
  }
  par(new = FALSE, pty = "m", plt = smallplot, err = -1)
  if (!horizontal) {
    image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
          ylab = "", col = col, breaks = breaks)
  }
  else {
    image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
          ylab = "", col = col, breaks = breaks)
  }
  if (!is.null(lab.breaks)) {
    axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2), 
                        at = breaks, labels = lab.breaks), axis.args)
  }
  else {
    axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)), 
                   axis.args)
  }
  do.call("axis", axis.args)
  if (!is.null(legend.lab)) {
    legend.args <- list(text = legend.lab, side = ifelse(horizontal, 
                                                         1, 4), line = legend.line, cex = legend.cex)
  }
  if (!is.null(legend.args)) {
    do.call(mtext, legend.args)
  }
  mfg.save <- par()$mfg
  if (graphics.reset | add) {
    par(old.par)
    par(mfg = mfg.save, new = FALSE)
    invisible()
  }
  else {
    par(big.par)
    par(plt = big.par$plt, xpd = FALSE)
    par(mfg = mfg.save, new = FALSE)
    invisible()
  }
}

#' unified gradient computation with penalties
#' @export
Grad_FullModel_Penalty <- function (A, grd, bas, Ref, S)
{
  C <- updateX(A, grd, bas)
  h2 <- E_SCM(A, grd, bas, C)
  G_SCM <- Grad_SCM(h2, grd, bas, C)
  if (!is.null(Ref)) {
    G_SEN <- Grad_SEN(A, grd, bas, G_SCM, S, Ref)
    G <- G_SCM$grad_SCM + G_SEN$grad_SEN
  }
  else G <- G_SCM$grad_SCM
  for (i in 1:bas$Nc) {
    G <- G + 2 * MemRBC_env$M.rho * G_SCM[[bas$Cons[i]]] * (h2[[bas$QCons[i]]] -
                                                   bas$Target[[bas$QCons[i]]])
  }
  return(G)
}

#' standardized model gradient with penalties
#' @export
E_FullModel_Penalty <- function (A, grd, bas, Ref)
{
  C <- updateX(A, grd, bas)
  h2 <- E_SCM(A, grd, bas, C)
  if (!is.null(Ref)) {
    S <- SEN(A, grd, bas, Ref, h2)
    e <- E_SEN(A, grd, bas, S, Ref)
  }
  else {
    e = 0
    S = NULL
  }
  E <- h2$Wb + e
  for (i in bas$QCons[1:bas$Nc]) E <- E + MemRBC_env$M.rho * (h2[[i]] -
                                                     bas$Target[i])^2
  names(E) = NULL
  return(list(E = E, Wb = h2$Wb, Ws = e, Wuncons = h2$Wb +
                e, dA = h2$dA, S = S, Area = h2$Area, Volume = h2$Volume,
              Curv = h2$Curv, n = h2$n))
}

#' create data
#' @description
#' use this to create some data for experimenting further...
#' 
#' @export
create_data <- function()
{ print("this is a long-runner to create data for MemRBC")
  p=paste(.libPaths(),"MemRBC")
  w=which(file.exists(p))
  n=p[w]
  print("creating data in ",n,"\n")
  MakeStandardRBC(L=16)->SS
  SDRC(SS,20,del_cons = 1e-1, del_min=1e-4,max_iter = 20,prn_ci = TRUE)->SSsdrc
  SDRC(SSsdrc,20,del_cons = 1e-1, del_min=1e-5,max_iter = 20,prn_ci = TRUE)->SSsdrc   
  SDRC(SSsdrc,20,del_cons = 1e-1, del_min=3e-5,max_iter = 20,prn_ci = TRUE)->SSsdrc   
  SDRC(SSsdrc,20,del_cons = 1e-1, del_min=5e-5,max_iter = 20,prn_ci = TRUE)->SSsdrc   
  SDRC(SSsdrc,20,del_cons = 1e-1, del_min=1e-4,max_iter = 20,prn_ci = TRUE)->SSsdrc   
  plot(SSsdrc$SDRC_Sample$E,type="l")
  plot(SSsdrc$SDRC_Sample$C,type="l")
  # the following CNM can have massive speedup with cluster=TRUE, MemRBC_env$M.Rcpp_ncores=2, 
  #       ncores=60, e.g. on a fat node with 120 cores and enough memory 
  # 
  CNM(SSsdrc,2,cluster = FALSE)-> SSsdrc_cnm # takes >100minutes on Ryzen 7 Pro laptop
  save_MemRBC(SSsdrc_cnm,"SSsdrc_cnm",qs2=TRUE) # fast compressed save

}

#' get_data_ZENODO
#' @description
#' get the additional data by downloading from ZENODO
#' and moving to data-folder
#' @param local (=FALSE) to download in current dir, not packages data
#' @param folder (="data") name of local download folder (is kept on exit)
#' @param L (=NULL) specific files for download (vector of character)
#' @export
get_data_ZENODO <- function(local=FALSE,folder="data",L=NULL){
  if(!dir.exists(folder)) dir.create(folder)
  Pre="https://zenodo.org/records/18667916/files/"
  L0=c("D5.rda","D6.rda","D6_cnm_C00.rda","D6_eq.rda","KleinBottle6.rda","L5_stoma_4.rda","ss42denovo.rda",
       "L5_stoma_6.rda","L5_stomatocyte_equilib.rda","L9_Stomatocyce_6.rda",
       "M_stomatocyte_L12.rda","M_stomatocyte_L12_rotated.rda","S18-C0--5.qs2",
       "S18-C0-17.qs2","S18-C0-23.qs2","S18-C0-27.qs2","SF4obj.rda","SS.rda","SS20.rda","SS42.rda","ss42denovo_mmc.rda",
       "ss42denovo_pnem.rda","U17R.rda")
  if (is.null(L)) L=L0 else if(! all(L %in% L0)) stop("req. data not in repo.\n")
  for (x in L){
    destfile=paste(folder,x,sep="/")
    cat(destfile,"\r")
    if (!file.exists(destfile)) utils::download.file(paste(Pre,x,"?download=1",sep=""),
                  destfile=destfile,
                  method="curl") 
  }
  if (!local){
  p=.libPaths()
   if (dir.exists(paste(p[1],"MemRBC/data",sep="/")))
    for (x in L){
      cat(x,"\r")
      cpfile=paste(folder,x,sep="/")
      target=paste(p[1],"MemRBC/data",x,sep="/")
      if (!file.exists(target)) try(file.copy(cpfile,target))
    }
  }
  cat("data left in folder ",folder,"\n")
}