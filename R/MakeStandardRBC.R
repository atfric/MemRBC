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
  M.C0<<-C0
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
                     Ref_hash=rlang::hash(Ref), Lambda=c(1,1),
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

#' @export
MakeStomatocyteRBC<-function(L=5 )
{
  data("L5_stomatocyte_equilib", envir = environment())
  SetParams(L5_stomatocyte_equilib)
  S5<-L5_stomatocyte_equilib
  update(S5,what=c("Grid","Basis","Ref"),n=(L+1)*5+2,L=L) -> M
  return(M)
}
