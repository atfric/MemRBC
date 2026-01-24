#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/


// [[Rcpp::depends(RcppCGAL)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]  
#include <Rcpp.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;

using namespace Rcpp;

// [[Rcpp::export]]
List cpp_delaunay2(NumericMatrix pts) {
  Delaunay dt;
  for(int i = 0; i < pts.nrow(); ++i) {
    dt.insert(K::Point_2(pts(i,0), pts(i,1)));
  }
  // extract triangles
  std::vector< std::vector<int> > triangles;
  int idx = 0;
  std::map<K::Point_2,int> index_of;
  for(auto it = dt.finite_vertices_begin(); it != dt.finite_vertices_end(); ++it) {
    index_of[it->point()] = idx++;
  }
  for(auto fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); ++fit) {
    triangles.push_back({
      index_of[fit->vertex(0)->point()],
              index_of[fit->vertex(1)->point()],
                      index_of[fit->vertex(2)->point()]
    });
  }
  return List::create(
    _["triangles"] = triangles
  );
}
