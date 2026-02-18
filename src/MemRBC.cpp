#include<Rcpp.h>
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

// for a faster version of cotan Laplacian GEMINI-Pro was asked
// to derive a R to C++ translation for Rcpp

//GEMINI
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>

using namespace Eigen;

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

//[[Rcpp::export(.L_Ylm)]]
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
//[[Rcpp::export(.Ylm_v)]]
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


//[[Rcpp::export(.Ylm_vv)]]
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


//[[Rcpp::export(.L_Ylm_u)]]
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


//[[Rcpp::export(.Ylm_uv)]]
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


//[[Rcpp::export(.Ylm_uu)]]
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
//[[Rcpp::export(.dot2)]]
double dot2(NumericVector x, NumericVector y) {
  return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}

//[[Rcpp::export(.IntegM)]]
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

// [[Rcpp::export(.E_SCM_cxx)]]
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


//[[Rcpp::export(.Grad_SCM_cxx)]]
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

