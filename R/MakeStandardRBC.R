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
#' @param C0 spontaneous curvature
#' @param del_cons (=0.1) constraint iterator step size
#' @param A0 to be constraint area
#' @param V0 to be constraint volume
#' @param V0_Ref volume of reference shape (148 ~ v_0=0.95)
#' @param n ~sqrt(spatial grid points)
#' @param prn boolean, print infos on constraint iteration
#' @param plt (=FALSE) for control plot 
#' @return membrane object with updates from MMC with:
#' @return LA: list of recorded coefficients A
#' @return A: last coefficients
#' @return bas: basis functions
#' @return grd: Gauss-Legendre-Simpson-spatial grid
#' @return Ref: Cauchy-Riemann tensor data for reference shape
#' @return ARef: coefficients of reference shape
#' @return history: history of App-calls that created the result
#' @examplesIf exists("L_Ylm")
#' MemRBC_env$M.C0 <- 0
#' # for fast experiments take a low-order L=5:
#' M <- MakeStandardRBC(L=5,prn=TRUE)
#' plot(M)
#' PSD(M,plt=TRUE)
#' M
#' 
#' @export
MakeStandardRBC <- function(A0=140, V0=100, V0_Ref=148, L=9,
                            C0=-1,n=(L+1) * 5 + 2, prn=FALSE,
                            del_cons=0.1, plt=FALSE)
{ t0=proc.time()
  data("Mempty",package="MemRBC",envir = environment())
  if (MemRBC_env$M.C0!=C0) warning("MemRBC_env$M.C0 not equal demanded C0")
  cl <- match.call()

  grd <- MakeGrid_GaussLegendreSimpson(n)
  bas.axi <- MakeBasis_UV(1,grd$U,grd$V)
  A=LM2A(MakeSphere(grd,bas.axi,3.34),bas.axi)

  updateX(A,grd,bas.axi)->C
  Wb<-E_SCM(A,grd,bas.axi,C)
  g2<-Grad_SCM(Wb,grd,bas.axi,C)
  # re-iterate to new V0_Ref if needed
  bas.axi$Nc=2
  bas.axi$Target=c(A0,V0_Ref)
  bas.axi$QCons=c("Area","Volume")
  bas.axi$Cons=c("gradA","gradV")
  ConsRHS(Wb,bas.axi)

  # compute Reference shape for A0, V0_Ref
  CI=ConsIter(A,grd,bas.axi,C,g2,Ctol=1e-6,del_cons = del_cons,
              nsteps = 1000,prn=prn)
  Aref=CI$A
  attr(Aref,"V0")=V0_Ref
  Filter_1_A<-function(A) {A[1,2:3]=A[2,1:2]=A[3,c(1,3)]=0;return(A)}
  Aref=Filter_1_A(Aref)
  updateX(Aref,grd,bas.axi)->C
  if(plt)rgl::plot3d(C$X,aspect=FALSE)
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
  if(plt)plot(M)
  Quantities(M)

  Mempty->D5
  SetParams(D5)
  D5$ARef
  D5<-update(D5,c("Grid","Basis","Class","Ref","Coor"))
  transplant(D5,M)->M1
  if(plt) {rgl::open3d();plot(M1)}
  t1=proc.time()
  M1$proc_time <-  t1-t0
  M1$comment="created from scratch"
  StoreParams(M1)->M1
  SetParams(M1)
  update(M1,"Class")->M1
  return(M1)
} # end of MakeStandardRBC

#' MakeDiscocyteRBC
#' @description
#' build from data M4 a higher spectral order discocyte object
#' @param L : spectral order of created object, L>5
#' @export
MakeDiscocyteRBC<-function(L)
{
  data("M4", package = "MemRBC",envir = environment())
  update(M4,what=c("Grid","Basis","Ref"),n=(L+1)*5+2,L=L)->M
  transplant(M4,M)->M
  M$proc_time=M4$proc_time
  StoreParams(M)->M
  return(M)
}

#' make invaginated shape MemRBC object
#' at the north pole by adding a gaussian displacement along Z
#' @param X spatial3d- coordinates
#' @param width (=0.1) width in micrometer
#' @param depth (=0.1) depth in micrometer
#' @return modified spatial 3d coordinates 
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
#' Create a spiculated cell from a L=16 sphere.
#' via invag_N() with negative depth. The basic object started from is a L=16 
#' standard sphere. The resulting shape is rescaled to obtain the target area specified in data(SS) (usually 140). You may also set a different volume target (see example code).
#'  WARNING: the reference is taken from unscaled sphere - may need updates!
#' data "SS.rda" is loaded from ZENODO to a local folder ./data/.
#' NOTE: spatial points are weighted by sin(u) for the fit
#' @param N (=20) from (6,8,12,20,32,42) number of regularly arranged spikes
#' @param w (=0.01) width of spikes; 0.01 is good for N=42
#' @param d (=-0.35) negative elevation of spikes; 0.35 is good for N=42
#' @param r (=1) radius of initial sphere. Helps to control initial volume.
#' @param rz (=1) helps to scale z, rz<1 : oblate, rz>1 : prolate shape
#' @return MemRBC object with spiculi, derived from published data "SS.rda"
#' @examplesIf exists("L_Ylm")
#' MakeSpiculated(N=12,w=0.06,d=-0.6,r=0.65) -> ss
#' ss
#' ss$bas$Target
#' #instead of tuning r,tz one may use ss$bas$Target["Volume"]=Volume(ss)
#'
#' rgl::open3d()
#' plot(ss)
#' MemRBC_env$M.C0 <- 20
#' MMC(ss,100000,plt=TRUE,pltfreq=100,LAfreq=1000,C0=20) -> ss_mmc
#' 
#' @export
MakeSpiculated<-function (N = 42, w = 0.01, d = -0.35, r = 1, rz = 1)
{ get_data_ZENODO(L="SS.rda",local=TRUE)
  load_MemRBC("data/SS.rda")->SS
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
#' @param f (=0.75) factor for coeff A[2,3]
#' @param uv (=c(0,0)) rotation before invagination
#' @param plt (=FALSE) for control plot
#' @return MemRBC object with invagination coded in the coefficients $A
#' @examplesIf exists("L_Ylm")
#' MakeInvaginated() -> ivs
#' ivs
#' ivs$bas$Target
#' rgl::open3d()
#' plot(ivs)
#' MemRBC_env$M.C0 <- -6
#' MMC(ivs,100000,plt=TRUE,pltfreq=100,LAfreq=1000,C0=-6) -> ivs_mmc
#' 
#' @export
MakeInvaginated <- function (w = 0.38, d = 1.4, r = 1.045, f = 0.75, uv = c(0, 0),
          plt = FALSE)
{ get_data_ZENODO(L="SS.rda",local=TRUE)
  load_MemRBC("data/SS.rda")->SS
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

#' Kleins bottle parameterization
#' @param uv n x 2 matrix of spherical angles (u,v) (u is internally scaled *2). uv can be created from `MakeBasis_UV()`
#' @param b (=2) width
#' @param h (=6) approx. height
#' @param plt (=FALSE) if 3d-plot is wanted
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

#' Make a Klein bottle that is a bit open at the neck
#' the integratiön domain is not deltau...pi-deltau, i.e. caps of the spectral fit are suppressed
#' @param deltau (=0.15) cap-cutting parameter
#' @param L (=17) spectral order
#' @param n (=L*4) grid points per 2d-dimension
#' @param b,h approx. width and height of bottle
#' @param plt (=FALSE) TRUE for 3d plot
#' @param SEN (=FALSE) TRUE to compute reference parameters for SEN in $Ref and SEN-parameters in $SEN
#' @return MemRBC object with Klein Bottle data
#' @examplesIf exists("L_Ylm")
#' MemRBC_env$M.C0=0
#' MemRBC_env$M.Ka=0
#' MemRBC_env$M.mu=0
#' M <- MakeKleinBottle(L=12,deltau=0.1)
#' plot(M)
#' Quantities(M)
#' StoreParams(M)->M
#' M$bas$Target
#' M
#' update(M,"Obj") -> M
#' imag.obj.colorbar(M$grd$Obj,M$grd$v,pal=topo.colors,alpha=0.65)
#' 
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
  unlist(Quantities(res))->q
  res$bas$Nc=2
  res$bas=SetConstraints(B,Cons = c("gradA","gradV"), QCons = c("Area","Volume"), Target = q[1:2])
if (SEN) {MakeRef(res,A) -> res$Ref; update(res,"SEN")}
if(plt) plot(res,color="white")
 return(res)
}

#' Make torus membrane, i.e. genus 1
#' set SEN reference to initial shape, i.e., E_SEN=0
#' @param L spectral order
#' @param n (=L*6) grid dimension
#' @param r (=1) smaller radius of torus
#' @param R (=2.5) greater radius of torus
#' @param plt (=FALSE) for control plots
#' @return MemRBC object for a torus
#' @export
MakeTorus<-function(L=5, n=L*6, R=2.5, r=1, plt=FALSE)
{
  G<-MakeGrid_Fourier(n=n,R=R,r=r,check_plt = plt)
  #rgl::shade3d(G$Obj)
  B<-MakeBasis_UV(L, u=G$U, v=G$V, kind = "Fourier")
  dim(B$A)
  dim(B$Ylm)
  X=Obj2X(G$Obj)
  if(plt) {rgl::open3d()
   rgl::shade3d(G$Obj,col="red",alpha=0.5)
  }
  B$mask=double_uv_ind(G$U,G$V)
  if(plt) rgl::plot3d(X,aspect=FALSE)
  lm( X[-B$mask,] ~ B$Ylm[-B$mask,] - 1 )$coefficients -> A # intercept in Ylm
  
  w=which(is.na(A[,1]))
  print(var(B$Ylm[,3]))
  w=c(3,w) # remove intercept from basis as well
  length(w)
  A=A[-w,]
  B$Ylm=B$Ylm[,-w]
  B$Ylm_u=B$Ylm_u[,-w]
  B$Ylm_v=B$Ylm_v[,-w]
  B$Ylm_uu=B$Ylm_uu[,-w]
  B$Ylm_uv=B$Ylm_uv[,-w]
  B$Ylm_vv=B$Ylm_vv[,-w]
  B$Ai_max=B$Ai_max-length(w)
  B$LM=B$LM[-w,]; B$l=B$l[-w];B$m=B$m[-w]
  B$G.tk=B$G.tk[-w];B$Wt=B$Wt[-w]
  LM2A(A,B)->A
  head(A,8)
  dim(A)
  MakeMemRBC(A,G,B)->M
  print(q<-unlist(Quantities(M)))
  M<-SetConstraints(M,Cons = c("gradA","gradV"), QCons = c("Area","Volume"), Target = q[1:2])
  update(M,c("Coor","Class","Obj"))->M
  if(plt) {plot(M,alpha=0.5,col="white");
   rgl::open3d()
   rgl::shade3d(M$grd$Obj,col="red",alpha=0.5)}
  M=MakeRef(M,M$A)
  M$ARef=M$A
  StoreParams(M)->M
  M$history=list(match.call())
  M$proc_time=0
  return(M)
}
