# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2025 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, S. (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# Frickenhaus S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#


#' Make a standard RBC object
#'
#' run MakeStandardRBC() to create a membrane object
#' @param L spectral order l<=L
#' @param A0 to be constraint area
#' @param V0 to be constraint volume
#' @param V0_Ref volume of reference shape (148 ~ v_0=0.95)
#' @param n ~sqrt(spatial grid points)
#' @param prn boolean, print infos on constraint iteration
#' @param del delta>0 for iterations; decrease if convergence fails
#' @return membrane object with updates from MMC with data:
#' @return LA: list of recorded coefficients A
#' @return A: last coefficients
#' @return bas: basis functions
#' @return grd: Gauss-Legendre-Simpson-spatial grid
#' @return Ref: Cauchy-Riemann tensor data for reference shape
#' @return ARef:  coefficients of reference shape
#' @return history: history of App-calls that created the result
#' @examples
#' library(MemRBC)
#' load_param_MemRBC()
#' M.C0=0
#' # for fast experiments take a low-order L=5:
#' M <- MakeStandardRBC(L=5,prn=TRUE)
#' plot(M)
#' PSD(M,plt=TRUE)
#' M
#' @export
MakeStandardRBC <- function(A0=140, V0=100, V0_Ref=148, L=9,C0=-1,
                                n=(L+1) * 5 + 2, prn=FALSE,
                                del_cons=0.1,del=2e-7,dt=1e-5)
{ # symm=Axi for reduced basis
  # if(L<4) warning("MakeStandardRBC: probably no solution for L>5")
  t0=proc.time()
  if (M.C0!=C0) warning("M.C0 not equal demanded C0")
  cl <- match.call()

  grd <- MakeGrid_GaussLegendreSimpson(n)
  bas.axi <- MakeBasis_UV(1,grd$U,grd$V)
  A=LM2A(MakeSphere(grd,bas.axi,3.34),bas.axi)

  updateX(A,grd,bas.axi)->C
  Wb<-E_SCM(A,grd,bas.axi,C)
  g2<-Grad_SCM(Wb,grd,bas.axi,C)
  # re-iterate to new V0_Ref if needed
  bas.axi$Nc=2
  bas.axi$Target=bas.axi$TNorm=c(A0,V0_Ref)
  bas.axi$QCons=c("Area","Volume")
  bas.axi$Cons=c("gradA","gradV")
  ConsRHS(Wb,bas.axi)

  # compute Reference shape for A0, V0_Ref
  CI=ConsIter(A,grd,bas.axi,C,g2,Ctol=1e-4,del_cons = del_cons,
              nsteps = 1000,prn=prn)
  Aref=CI$A
  attr(Aref,"V0")=V0_Ref
  Filter_1_A<-function(A) {A[1,2:3]=A[2,1:2]=A[3,c(1,3)]=0;return(A)}
  Aref=Filter_1_A(Aref)

  updateX(Aref,grd,bas.axi)->C

  rgl::plot3d(C$X,aspect=FALSE)
  # now the target L shape by energy minimization
  bas <- MakeBasis_UV(L,grd$U,grd$V)

  A=bas$A
  A[]=0;A[1:3,]=Aref[]

  tictoc::tic()
  cat("Reference Tensors ...");
  Ref <- Ref4CauchyGreen( A, grd, bas )
  cat(" took ");tictoc::toc()

  bas$Target=c(A0,V0)
  names(bas$Target)=c("Area","Volume")
  bas$Nc=2
  A=A+rnorm(length(A),sd=0.01)
  M=structure(list(grd=grd, bas=bas, A=A, Ref=Ref, ARef=Aref,
                     hash_Ref=rlang::hash(Ref), Lambda=c(1,1),
                     LA=list(A), history=list(cl)), class="MemRBC")
  plot(M)
  Quantities(M)

  data("D5")
  SetParams(D5)
  transplant(D5,M)->M1
  rgl::open3d();plot(M1)
  t1=proc.time()
  M1$proc_time <-  t1-t0
  M1$comment="created from data(D5)"
  StoreParams(M1)->M1
  SetParams(M1)
  return(M1)
} # end of MakeStandardRBC


#' @export
MakeDiscocyteRBC<-function(L)
{
  data(D5, envir = environment())
  update(D5,what=c("Grid","Basis","Ref"),n=(L+1)*5+2,L=L)->M
  transplant(D5,M)->M1
  M1$proc_time=D5$proc_time
  StoreParams(M1)->M1
  return(M1)
}

#' make a stomatocyte from present data
#' @param L (=5) spectral order
#' @export
MakeStomatocyteRBC<-function(L=5 )
{
  data("L5_stomatocyte_equilib", envir = environment())
  SetParams(L5_stomatocyte_equilib)
  S5<-L5_stomatocyte_equilib
  update(S5,what=c("Grid","Basis","Ref"),n=(L+1)*5+2,L=L) -> M
  return(M)
}

#' @export
invag_N<-function (X, width = 0.1, depth = 0.1)
{
  u <- apply(X, 1, inv_sph)[1, ]
  w <- which(exp(-(u^2)/width) * abs(depth) > 0.001)
  X[w, 3] <- X[w, 3] - exp(-(u[w]^2)/width) * depth
  attr(X, "w") <- w
  return(X)
}

#
#' Create a spiculated cell from a sphere
#' via invag_N() with negative depth. The basic object started from is a L=16 standard sphere. The resulting shape is rescaled to obtain the target area specified in data(SS) (usually 140). You may also set a different volume target (see example code).
#'  WARNING: the reference is taken from unscaled sphere - may need updates!
#' NOTE: spatial points are weighted by sin(u) for the fit
#' @param N (=20) (6,8,12,20,32,42) number of regularly arranged spikes
#' @param w (=0.01) width of spikes; 0.01 is good for N=42
#' @param d (=-0.35) negative elevation of spikes; 0.35 is good for N=42
#' @param r (=1) radius of initial sphere. Helps to control initial volume.
#' @param rz (=1) helps to scale z, rz<1 : oblate, rz>1 : prolate shape
#' @examples
#' # example code for spiculated sphere ss
#' MakeSpiculated(N=12,w=0.06,d=-0.6,r=0.65) -> ss
#' ss
#' ss$bas$Target
#' #instead of tuning r,tz one may use ss$bas$Target["Volume"]=Volume(ss)
#'
#' rgl::open3d()
#' plot(ss)
#' M.C0 <- 20
#' MMC(ss,100000,plt=TRUE,pltfreq=100,LAfreq=1000,C0=20) -> ss_mmc
#' @export
MakeSpiculated<-function (N = 42, w = 0.01, d = -0.35, r = 1, rz = 1)
{
  data(SS)
  g = SS$grd
  b = SS$bas
  if (r <= 0)
    stop("non-positive radius; use e.g. r=1")
  a = MakeSphere(g, b) * r
  a[2, 3] = a[2, 3] * rz
  if (N == 42)
    S = Rvcg::vcgSphere(1)
  else if (N == 12)
    S = Rvcg::vcgIcosahedron()
  else if (N == 6)
    S = Rvcg::vcgOctahedron()
  else if (N == 32)
    S = Rvcg::vcgDodecahedron()
  if (N %in% c(12, 42, 6, 32))
    s = S$vb[1:3, ]
  else if (N == 20)
    s = t(Rvcg::vcgBary(Rvcg::vcgSphere(0)))
  else if (N == 8)
    s = t(Rvcg::vcgBary(Rvcg::vcgOctahedron()))
  else return(NULL)
  W = c()
  uv = (apply(s, 2, inv_sph))
  range(uv[2, ])
  X <- updateX_only(a, g, b)$X
  rgl::mfrow3d(3, 4, sharedMouse = TRUE)
  for (i in 1:dim(uv)[2]) {
    X <- rotateX(X, uv[1, i], 0, uv[2, i], transpose = TRUE)
    X <- invag_N(X, w, d)
    rgl::plot3d(X, aspect = FALSE)
    X <- rotateX(X, uv[1, i], 0, uv[2, i], transpose = FALSE)
  }
  col = rep(1, dim(X)[1])
  rgl::open3d()
  rgl::plot3d(X, aspect = FALSE, col = col)
  O = SS$grd$Obj
  O <- X2Obj(O, X)
  rgl::wire3d(O)
  b$mask
  A <- FitAlm_Tikhonov(X = X, bas = b, WX=sin(SS$grd$U), lambda = 0)
  SS$A = A
  a0 = Area(SS)
  (la = (SS$bas$Target["Area"]/a0)^(1/2))
  SS$A <- SS$A * la
  v0 = Volume(SS)
  a0 = Area(SS)
  c0 = Curv(SS)
  cat("Quantities A,V,C:", a0, v0, c0, "\n")
  SS$history <- match.call()
  SS$proc_time = NULL
  return(SS)
}

#' Create a invaginated cell from a sphere
#' via invag_N() positive depth.
#' NOTE: spatial points are weighted by sin(u) for the fit
#' @param w (=0.38) width of grove
#' @param d (=1.4) depth of grove
#' @param r (=1.045) radius to tune resulting volume
#' @examples
#' # example code for invaginated sphere ivs
#' MakeInvaginated() -> ivs
#' ivs
#' ivs$bas$Target
#'
#' rgl::open3d()
#' plot(ivs)
#' M.C0 <- -6.5
#' MMC(ivs,100000,plt=TRUE,pltfreq=100,LAfreq=1000,C0=-6.5) -> ivs_mmc
#' @export
MakeInvaginated <- function (w = 0.38, d = 1.4, r = 1.045, f = 0.75, uv = c(0, 0),
          plt = FALSE)
{
  data(SS)
  g = SS$grd
  b = SS$bas
  if (r <= 0)
    stop("non-positive radius; use e.g. r=1")
  a = MakeSphere(g, b) * r
  a[2, 3] = f * a[2, 3]
  X <- updateX_only(a, g, b)$X
  X <- rotateX(X, uv[1], 0, uv[2], transpose = TRUE)
  X <- invag_N(X, w, d)

  if (plt)
    rgl::plot3d(X, aspect = FALSE)
  X <- rotateX(X, uv[1], 0, uv[2], transpose = FALSE)
  col = rep(1, dim(X)[1])
#  col[as.numeric(names(table(W)))] = table(W)
  col="red"
  if (plt) {
    rgl::open3d()
    rgl::plot3d(X, aspect = FALSE, col = col)
  }
  O = SS$grd$Obj
  O <- X2Obj(O, X)
  rgl::wire3d(O)
  A <- FitAlm_Tikhonov(X = X, bas = b, WX=sin(g$U), lambda = 0)
  SS$A = A
  a0 = Area(SS)
  (la = (SS$bas$Target["Area"]/a0)^(1/2))
  SS$A <- SS$A * la
  v0 = Volume(SS)
  a0 = Area(SS)
  c0 = Curv(SS)
  cat("Quantities A,V,C:", a0, v0, c0, "\n")
  SS$history <- match.call()

  SS$proc_time = NULL
  return(SS)
}


#' @export
KleinB <- function(uv,b=2,h=6,plt=FALSE){
  u=2*uv[,1];v=uv[,2]
  r=2-cos(u)
  x = b*(1-sin(u))*cos(u)+r*cos(v)*(2*exp(-(u/2-pi)^2)-1)
  y = r*sin(v)
  z = h*sin(u) + 0.5*r*sin(u)*cos(v)*exp(-(u-1.5*pi)^2)
  X=cbind(x,y,z)
if (plt)  {rgl::open3d();rgl::plot3d(X,aspect=FALSE)}
  return(X)
}


#' @export
MakeGrid_Fourier<-function(n=30,r=1,R=2.5,check_plt=FALSE)
{
  grd=list(kind="Fourier",ua=0,ub=2*pi,va=0,vb=2*pi) # to be filled further and returned
  nu=n+1; nv=n+1 # we double the first u,v data point (u,v=0) at u,v=2pi
  x <- pracma::linspace(0,2*pi,n=nu)
  wx <- rep(2*pi/n,nu)
  y <- pracma::linspace(0,2*pi,n=nv)
  wy <- wx
  grd$xg <- x
  grd$yg <- y
  mesh=pracma::meshgrid(x, y)
  grd$u=t(mesh$X) # 2D # for [nu,nv] adressing
  grd$v=t(mesh$Y) # 2D
  (dm=dim(grd$u))
  grd$ndof=prod(dm)
  grd$U=as.vector(grd$u)
  grd$V=as.vector(grd$v)
  grd$wx=wx;  grd$wy=wy # for general integration scheme
  grd$UV=cbind(grd$U,grd$V)
  nx=nu;ny=nv;
  q=matrix(NA,3,nx*ny*2);k=0
  for (i in 1:(nx-1))  for (j in 1:(ny-1)){
    k=k+1;l=(j-1)*nx+i
    q[,k]=c(l,l+1,l+1+nx)
    k=k+1
    q[ ,k]=c(l,l+nx+1,l+nx)
  }
  q=q[,1:(k)]
  print(k)
  x=cos(grd$v)*(R+r*cos(grd$u));
  y=sin(grd$v)*(R+r*cos(grd$u));
  z=r*sin(grd$u)
  rgl::mesh3d(x,y,z,triangles=q, normals = list(x=x,y=y,z=z) ) -> M
  #  clear3d();
  grd$Obj<-M
  grd$comment<-comment
  grd$type="FOURIER"
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
  Rvcg::vcgUpdateNormals(grd$Obj)->grd$Obj
  return(grd)
}

#' Make a Klein bottle that is a bit open at the neck
#' the integratiön domain is not deltau...pi-deltau, i.e. caps of the spectral fit are suppressed
#' @param detlau (=0.15) cap-cutting parameter
#' @param L (=17) spectral order
#' @examples
#' M.C0=0;M.Ka=M.mu=0
#' M <- MakeKleinBottle(L=12,deltau=0.1)
#' #plot(M)
#' Quantities(M)
#' StoreParams(M)->M
#' M$bas$Target
#' M
#' update(M,"Obj") -> M
#' imag.obj.colorbar(M$grd$Obj,M$grd$v,pal=topo.colors,alpha=0.65)
#' @export
MakeKleinBottle<-function(L=17,n=L*4,b=2,h=6,deltau=0.15,plt=FALSE,SEN=FALSE)
{
  G1=MakeGrid_GaussLegendreSimpson(n=n,ua=0,ub=pi)
# rgl::plot3d(G$Obj,aspect=FALSE,col="white")
  B1=MakeBasis_UV(L,G1$U,G1$V,kind="Ylm")
  X=KleinB(B1$uv,h=h,b=b,plt=plt)
  X2Obj(G1$Obj,X)->O
if(plt)  rgl::shade3d(O,col="red",alpha=0.5)
  G=MakeGrid_GaussLegendreSimpson(n=n,ua=deltau,ub=pi-deltau)
  B=MakeBasis_UV(L,G$U,G$V,kind="Ylm")
  B$mask=mask=1
  lm(X[  ,] ~ B$Ylm[  ,])$coefficients[-1,] -> A
  MakeMemRBC(LM2A(A,B),G,B) -> res
  unlist(Quantities(M))->q
  res$bas$Nc=2
  res$bas=SetConstraints(B,Cons = c("gradA","gradV"), QCons = c("Area","Volume"), Target = q[1:2],TNorm = c(1,1,1))
if (SEN)  MakeRef(res,A) -> res$Ref
if(plt) plot(res,color="white")
 return(res)
}

# Make torus membrane, i.e. genus 1
# but something is wrong here...
## @export
MakeTorus<-function(L=5,plt=FALSE)
{
  G<-MakeGrid_Fourier(n=L*5,check_plt = FALSE)
  B<-MakeBasis_UV(L,u=G$U,v=G$V,kind = "Fourier")
  X<-Obj2X(G$Obj)
  if(plt)rgl::plot3d(X)
  lm( X ~ B$Ylm - 1 )$coefficients -> A # intercept in Ylm[,1]
  dim(A)
  dim(B$Ylm)
  A[,2]<- -A[,2]
  MakeMemRBC(A,G,B)->M
  A=MakeSphere(G,B)
  (q<-unlist(Quantities(M)))
  M$bas<-SetConstraints(M$bas,Cons = c("gradA","gradV"), QCons = c("Area","Volume"), Target = q[1:2],TNorm = c(1,1,1))
  plot(M)
  return(M)
}
