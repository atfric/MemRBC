library(Rcpp)
library(RcppCGAL)
# Source the C++ code (assuming it's in a file called 'triangle_intersection.cpp')
sourceCpp("../../../MemRBC/triangle_intersection.cpp")

# Define the points of the two triangles in R
# Triangle 1: (0, 0, 0), (1, 0, 0), (0, 1, 0)
t1 <- c(0, 0, 0)
t2 <- c(1, 0, 0)
t3 <- c(0, 1, 0)

# Triangle 2: (0, 0, 1), (1, 0, 1), (0, 1, 1)
t4 <- c(0, 0, 1)
t5 <- c(1, 0, 1)
t6 <- c(0, 1, 1)


Rcpp::sourceCpp("../../../MemRBC/del.cpp")

pts <- matrix(c(
  0,0,
  1,0,
  1,1,
  0,1,
  0.5, 0.5
), ncol=2, byrow=TRUE)

tri <- cpp_delaunay2(pts)
print(tri)


remotes::install_github(
  "stla/RCGAL", dependencies = TRUE, build_opts = "--no-multiarch"
)
Rcpp::sourceCpp("../../../MemRBC/hilb.cpp")


