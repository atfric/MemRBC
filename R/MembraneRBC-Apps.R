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
# high level driver "Apps" working with MemRBC object
#   are stored in seperate R-scripts
#
# remarks:
#   for statistics output per iteration step, create a text-output via sink()
#

# tests in this part to be run

Apps.TEST=FALSE

# some generics

#' image
#'
#' 2D plot of color coded variable on the shape
#' @param M Membrane to use for shape plot
#' @param q a 2D field to project onto shape, eg. M$SEN$alpha
#' @export
image.MemRBC<-function(M,q=M$SEN$alpha,main=expression(alpha),...)
  {
    imag(q,M$grd,main=main,...)
    invisible()
}

#' plot
#'
#' 3D plot of membrane shape
#' @param M Membrane to use for shape plot
#' @param wire (boolean) FALSE if  to suppress mesh quads plotted
#' @param ... further plotting options, e.g. col="green"
#' @export
plot.MemRBC <- function(R,wire=TRUE,wire_col="black",...){
    updateX_only(R$A,R$grd,R$bas)->C
  #  Rvcg::vcgUpdateNormals(C$Obj)->O
    plot3b(C$X,R$grd,wire=FALSE,...)
    if (is.null(R$grd$ObjQ)) Obj2ObjQ(R$grd$Obj,R$grd)->Q else Q=R$grd$ObjQ
    X2ObjQ(Q,C$X)->Q
   # Q$normals=O$normals
    if (wire) rgl::wire3d(Q,col=wire_col,specular="black")
}

#' PlotDiff
#' @description
#' plot the spatial difference vectors as wireframe.
#' A quite simmilar plot can be obtained by plot(M2-M1).
#' @param M1,M2 membrane objects; displayed delta is M2$A-M1$A
#' @param On3d (=FALSE) to plot on M1 surface the |dX| in color code
#' @param ... graphics parameters to wire3d() or imag.obj.colorbar()
#' @return named vector of changes of Quantities()
#' @examples
#' data(ss42denovo)
#' data(ss42denovo_pnem)
#' dQ <- PlotDiff(ss42denovo_pnem,ss42denovo,col="red")
#' rgl::open3d()
#' plot(ss42denovo_pnem-ss42denovo,col="green")
#' @export
PlotDiff<-function(M1,M2,On3d=FALSE,...)
{  updateX(M2$A-M1$A,M1$grd,M1$bas)->C

   X2Obj(M2$grd$Obj,C$X) -> O
   if(!On3d)rgl::wire3d(O,...) else {
     updateX(M1$A,M1$grd,M1$bas)->C1
     X2Obj(M2$grd$Obj,C1$X) -> O1
     imag.obj.colorbar(O1,apply(C$X,1,pracma::Norm),par=FALSE)
   }
   return( c(sapply(Quantities(M2),function(x) x)-sapply(Quantities(M1),function(x) x),apply(C$X,1,pracma::Norm)))
}

#' print
#'
#' print metadata from membrane object
#' @param M Membrane to print data from
#' @export
print.MemRBC<-function(R){
    cat("MembraneRBC Object:","\n")
    cat("ndof_grid=",R$grd$ndof," [",R$grd$nu,"x",R$grd$nv,"]\n")
    if(!is.null(R$bas$Ylm_u))
    {cat("L=",R$bas$L_max,", coeff. rows",R$bas$Ai_max,"\n")
     cat(R$bas$Nc," constraints target:\n");
     print(R$bas$Target)
     cat("Membranes recent coefficients A computed at C0=",attr(R$A,"C0"),"\n")
     if (!rlang::is_named(R$bas$Target)) message("Object constraint $Target is unnamed!")
     cat("quantities for current coeffs A: Area",Area(R)," Volume",Volume(R)," Curv",Curv(R),"\n")
    } else cat(crayon::cyan("The Basis in this Membrane is reduced to X!\nUse update(M,\"Basis\") to recreate.\n"))
    if (!is.null(R$LA)) cat("List of coeffs LA length:",length(R$LA),"\n")

    if (!is.null(R$kT)) cat("MMC kT:",R$kT,"\n")

    if (!is.null(R$visc)) if(length(c(R$visc))==1) cat("PNEM viscosity:",R$visc,"\n")
                          else cat("PNEM viscosity (head of diag of matrix)",diag(R$visc)[1:6],"\n")
    if (!is.null(R$rho)) cat("PNEM mass density:",R$rho,"\n")
    if (!is.null(R$Sample)) cat("MMC Samples:",dim(R$Sample)[1],"\n")
    if (!is.null(R$proc_time)) cat("Membrane has overall process time ",R$proc_time,"\n") else cat("Membrane has no process_time recorded.\n")
    if (!is.null(R$comment)) cat("Membrane comment ",R$comment,"\n") else cat("Membrane has no comment.\n")
    # cat(" Coordinates hash",attr(R$C,"hash"),"\n")
    #  cat("Coefficients hash",rlang::hash(R$A),"\n")
    #  cat("      Object hash",rlang::R$hash)
    cat("Membrane has a history of length ",length(R$history),"\n")
    if(is.null(R$Params)) cat(crayon::red("No Parameters - inject by StoreParams(M)->M\n")) else {
      ul=unlist(R$Params); cat("Params:\n"); print(ul)
    }
}



#' Plot a stability analysis visually on a 3x3 canvas of 3d-plots
#'
#' run PlotStabGallery on a membrane object with a Stab in it from MemStab()
#' @param M The input membrane with initial data and reference
#' @param which_n vector of Ids of eigenvalues, 1:9 is default for lowest nine
#' @param plt_scale scale of eigenvectors
#' @param col1 color of unperturbed shapes
#' @param col2 color of scaled eigenvector perturbed shapes
#' @param sharedMouse (boolean) if FALSE, each subplot can be controlled by mouse individually
#' @return none
#'@examples
#' M <- MakeStandardRBC(L=5)
#' plot(M)
#' M <- PSD(M,nsteps=100,del=1e-6,)
#' plot(M)
#' MemStab(M)->M
#' PlotStabGallery(M)
#' @export
PlotStabGallery <- function(M,which_n=1:9,sharedMouse=TRUE,wire=FALSE,...)
{ if (is.null(M$Stab)) stop("No stability data in MemRBC object")
  rgl::open3d();rgl::mfrow3d(3,3,sharedMouse = sharedMouse);

  for (i in which_n)
  {if (i>1) rgl::next3d();
    plotStab(M,plt_n=i,wire=wire,...)
    N=length(M$Stab$EigH$values)
    rgl::title3d(paste(i,"ev",round(M$Stab$EigH$values[N+1-i],6)))
  }
}

#' @export
plotStab<-function(M,plt_n=1,plt_scale=0.5,col1="blue",col2="red",alpha1=0.5,alpha2=0.5,wire=FALSE)
{
  plot(M,alpha=alpha1,col=col1,wire=wire);
  N=dim(M$Stab$EigH$vectors)[1]
  v=M$Stab$EigH$vectors[,N+1-plt_n]
  v=v/max(v)*plt_scale
  M$A=M$A+v
  plot(M,alpha=alpha2,col=col2,wire=wire)
}

#
# Membrane stability analysis: unconstrained Hessian computation + Eigenvalues/-vectors
#
#' MemStab
#'
#' analyse membrane shape stability in terms of eigenvectors of Hessian
#' @return M$Stab with eigensystem of Hessian of energy without constraints
#' @export
MemStab <- function(M, mc.cores = 4, plt=FALSE, plt_mode=1, plt_scale=0.5, serial = FALSE)
  { t0=proc.time()
    if(is.null(M$proc_time)) M$proc_time=0
    cl=match.call()
    if (serial) H=FullModelHessian(M$A,M$grd,M$bas,M$Ref) else
      H=FullModelHessian_Par(M$A,M$grd,M$bas,M$Ref, del = 5e-06,
                             Mem_mc.cores = mc.cores, timing = TRUE, startup = TRUE,
                             stopdown = TRUE)
    E=eigen(H$H)
    M$Stab=list(Hessian=H, EigH=E )
    if (plt) plotStab(M,plt_mode,plt_scale)
    M$last_App_called="MemStab"
    M$history=append(M$history,list(cl))
    t1=proc.time()
    M$proc_time <- M$proc_time + t1-t0

    return(M)
}

#' massmatrix
#'
#' compute the mass matrix for PNEM
#' @param M The input membrane with initial data and reference
#' @param rho area density, default 1
#' @return massmatrix in spectral space
#'@examples
#' data("M4")
#' mass <- massmatrix(m4,rho=2)
#' image(mass)
#' @export
massmatrix <- function(M,rho=1)
{ t0=proc.time()
  if (!is.null(M$dA)) q=M$dA else
{q=E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))$dA}
if (length(q)!=M$grd$ndof)   {q=E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))$dA}
if (length(q)!=M$grd$ndof) stop("wrong size in mass matrix dA vs. grid size")
cat("Integration of mass matrix density\n")
tictoc::tic()

mass=matrix(0,M$bas$Ai_max,M$bas$Ai_max)

for (i in 1:M$bas$Ai_max) {if(i%%15==0) cat(round(i/M$bas$Ai_max*100,1),"\r")
  for (j in i:M$bas$Ai_max) {YY=M$bas$Ylm[,j]*M$bas$Ylm[,i];
    mass[j,i] <- mass[i,j] <- rho*.IntegS( q * sin(M$grd$U) * YY, M$grd)
  }
  }
cat("\r      ");tictoc::toc()
attr(mass,"rho")=rho
attr(mass,"single_proc_time")=proc.time()-t0
return(mass) # could be stored as density M$rho_PNEM
}

#' Quantities
#'
#' report some quantities from a membrane object
#' @param M The input membrane with initial data and reference
#' @return quantities Area, Volume, curvature Curv and bending energy Wb as named vector
#'@examples
#' data("M1")
#' Quantities(M1)
#' @export
Quantities<-function(M)
{r=E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))[c("Area","Volume","Curv","Wb")]
 names(r)=c("Area","Volume","Curv","Wb")
 return(unlist(r))
}


#' Energy
#'
#' report  energy values for bending and SEN from a membrane object
#' @param M The input membrane with initial data and reference
#' @return vector of energies Wb (bending), Es (stress-shear), E (potetntial energy), Ekin (kinetic energy, optional, if a PNEM was run)
#'@examples
#' data("M1")
#' Energy(M1)
#' @export
Energy<-function (M)
{
  C <- updateX(M$A, M$grd, M$bas)
  h2 = E_SCM(M$A, M$grd, M$bas, C)
  if (!is.null(M$Ref)) {
    S = SEN(M$A, M$grd, M$bas, M$Ref, h2)
    ES = E_SEN(M$A, M$grd, M$bas, S, M$Ref)
  }
  else ES = 0
  r = c(h2$Wb, ES, h2$Wb + ES)
  names(r) = c("Wb", "Es", "E")
  if (!is.null(M$mass) & !is.null(M$Av)) {
    Ekin <- 0
    for (j in 1:3) Ekin <- Ekin + (0.5 * M$Av[, j] %*% M$mass %*%
                                     M$Av[, j])[1, 1]
    r = c(r, Ekin)
    names(r)[4] = "Ekin"
  }
  if(M$bas$Nc>0) {r[5]=sum((unlist(Quantities(M))[M$bas$QCons]-M$bas$Target)^2)
   names(r)[5]="Econstr"}
  return(r)
}



#' @export
Area<-function(M)
  {return(E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))$Area)}

#' @export
Volume<-function(M)
  {return(E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))$Volume)}

#' @export
Curv<-function(M)
  {return(E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))$Curv)}

# re-orient shape by prinicpal axis of initertia
#' MemPCA
#'
#' rotate Membrane coefficients to coordinates of principle axes of inertia.
#' Remark: the reference SEN in $Ref remains unchanged.
#' @param M the input membrane to be rotated
#' @return  membrane object with new coefficients of PCA-rotated spatial coordinates.
#'
#' @examples
#' data("M_stomatocyte_L12"); SetParams(M_stomatocyte_L12)
#' plot(M_stomatocyte_L12)
#' MemPCA(M_stomatocyte_L12)->M2
#' plot(M2)
#' @export
MemPCA<-function(M,WX=rep(1,M$grd$ndof))
  { cl=match.call()
    updateX_only(M$A,M$grd,M$bas)->M$C
    princomp(M$C$X)$scores->X
    A=FitAlm_Tikhonov(X = X,bas=M$bas,lambda=0,WX=WX)
    M$A=A
    M$history=append(M$history,list(cl))
    return(M)
  }

# open two empty screens, keep devs as global M.scr1 and M.scr2
#' two_screens3d
#'
#' open two screens for 3D-plots by two_draw3d()
#' @param x,y dimensions in pixels
#' @return global variables M.scr1 and M.scr2 are set
#' @examples
#' two_screens3d(x=650,y=300)
#' @export
two_screens3d<-function(x=400,y=400){
  if (M.scr2 %in% rgl::rgl.dev.list()) {rgl::open3d();assign("M.scr2",rgl::cur3d(),envir=.GlobalEnv);rgl::par3d(windowRect=c(x-30,30,2*x-30,y+30))}; # draw right first
  if (M.scr1 %in% rgl::rgl.dev.list()) {rgl::open3d();assign("M.scr1",rgl::cur3d(),envir=.GlobalEnv);rgl::par3d(windowRect=c(1,30,x+1,y+30))}; # then left (for 30 pixel overlap)
  }

# draw stress and shear; re-open screens if not opened (not in rgl.dev.list())
#' two_draw3d
#'
#' plot membrane object with area stress alpha and shear beta as color
#' @param A coefficients to use
#' @param M the input membrane to be plotted
#' @param cont (=FALSE) draw gridlines u=const, v=const
#' @return -
#' @examples
#' data("M_stomatocyte_L12"); M<-M_stomatocyte_L12
#'  perturb original coefficients
#' A <- pertA_Gauss(M$A, M$bas, sd=0.2)
#' two_screens(); two_draw3d(A, M)
#' @export
two_draw3d<-function(A,M,cont=FALSE,title="",x=400,y=400) # requires S (Stretches) as global variable
  { grd=M$grd;
    C<-updateX(A,grd,M$bas)
    Wb<-E_SCM(A,grd,M$bas,C)
    S<-SEN(A,grd,M$bas,M$Ref,Wb)
    X2Obj(grd$Obj,C$X)->O
    Rvcg::vcgClean(O,sel=1:7,silent = TRUE)->O
    if (M.scr2 %in% rgl::rgl.dev.list()) {rgl::set3d(M.scr2);rgl::clear3d();} else { assign("M.scr2",rgl::open3d(),envir=.GlobalEnv);rgl::par3d(windowRect=c(x-30,30,2*x-30,y+30)); }

    imag.obj.colorbar(O,f=S$beta,clr = FALSE,par=FALSE,specular="black"); rgl::title3d(paste("beta",title))
if(cont){    rgl::contourLines3d(O,grd$U,nlev=15,lwd=2)
  rgl::contourLines3d(O,grd$v,levels = pracma::linspace(0,2*pi,16)[-16], lwd=2)
}
    if (M.scr1 %in% rgl::rgl.dev.list()) {rgl::set3d(M.scr1);rgl::clear3d();} else { assign("M.scr1",rgl::open3d(),envir=.GlobalEnv);rgl::par3d(windowRect=c(x-30,30,2*x-30,y+30)); }
    imag.obj.colorbar(O,f=S$alpha,clr=FALSE,par=FALSE,specular="black");rgl::title3d(paste("alpha",title))
if(cont){    rgl::contourLines3d(O,grd$U,nlev=15,lwd=1)
  rgl::contourLines3d(O,grd$v,nlevels=15, levels = pracma::linspace(0,2*pi,16)[-16],lwd=1)
  }
}

#' PlotA
#'
#' plot amplitudes in matplot
#'
#' @export
#' @param M membrane to plot M$A from
#' @param bar (=TRUE)if L-levels should be color-plotted
#' @param scale_up (=TRUE) plot A scaled-up by sqrt(bas$G.tk)
#' @examples
#' data("M4");SetParams(M4)
#' plotA(M4,ylab=expression(A[scaled]))
#' plot(M4,scaled_up=FALSE,ylab=expression(A[unscaled]))
#'
PlotA<-function(M,bar=TRUE,scale_up=TRUE,...)
{ if (scale_up) yl=expression(A[upscaled]) else yl=expression(A[nonscaled])
  plotA_l(M$A,M$bas,bar=bar,scale_up=scale_up,ylab=yl,...)
}

#' update
#'
#' update data in membrane object
#' @param M the input membrane to be updated
#' @param what vector of character from "dA", "Quantities", "Basis", "Grid", "Ref", "curv", "SCM", "X", "SEN"
#' @param n for what="Grid": new number of grid points
#' @param L for what="Basis": spectral order L
#' @return updated membrane object
#' @examples
#' data("M_stomatocyte_L12")
#' update(M_stomatocyte_L12,"X")->M
#' plot3d(M$X, aspect=FALSE)
#' # make lower spectral order membrane from M
#' update(M, what=c("Grid","Basis","Ref"), n=30, L=8)->L8
#' two_screens3d(); two_draw3d(L8$A, L8)
#' rgl::open3d(); plot(L8)
#' Quantities(L8)
#' Energy(L8)
#' @export
update.MemRBC <- function(M, what=c("dA","Quantities","X"),n=0,L=5)
{
  if(is.null(M$grd)) message("Membrane has no Grid!\nUse what=\"Grid\" in update(...,n=ngrid)")
  if("Grid" %in% what){
    grd=MakeGrid_GaussLegendreSimpson(n);
    M$grd=grd;
    M$mass=NULL;M$Ref=NULL
    what=c(what,"Basis")
  }
  if(is.null(M$bas)) message("Membrane has no basis!\nUse what=\"Basis\" in update(,L=Lmax)")

  if("curv" %in% what)
  {q=E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))$curv
   M$curv=q;}
  if("curv_sq" %in% what)
  {q=E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))$curv_sq
   M$curv_sq=q;}
  
  if("SCM" %in% what)
  {q=E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))
   M$SCM=q;
  }
  if("Basis" %in% what){
    if(is.null(M$bas$Target)) stop("to update Basis you need $Target values in original $bas.")
    Target=M$bas$Target
    bas=MakeBasis_UV(L,M$grd$U,M$grd$V);
    bas$Target=Target
    Ain=LM2A(M$A,M$bas);
    M$bas=bas;
    A=LM2A(M$bas$A,M$bas)
    i=intersect(rownames(A),rownames(Ain))
    A[i,]=Ain[i,]
    M$A=A
    if(sum(A)==0) stop("A=0 in update Basis after taking over old coefficients M$A")
    M$Av=NULL # no velocities
    what=c(what,"Ref")
  }
  if("Ref" %in% what){
    A=M$bas$A
    Ain=M$ARef;
    i=intersect(rownames(A),rownames(Ain))
    A[i,]=Ain[i,]
    M$ARef=A
    M$Ref=Ref4CauchyGreen(M$ARef,M$grd,M$bas)
    M$mass=NULL
    M$Av=NULL
  }
  if ("X" %in% what) M$X=updateX_only(M$A,M$grd,M$bas)$X
  if ("Obj" %in% what)
  { #  the following update
  X2Obj(M$grd$Obj,updateX_only(M$A,M$grd,M$bas)$X)->M$grd$Obj
  # make Obj invalid for plot(M):
  #   C=updateX(M$A,M$grd,M$bas)
  #   M$grd$Obj=X2Obj(M$grd$Obj,C$X)
  #  M$grd$Obj=Rvcg::vcgUpdateNormals(M$grd$Obj)
  #   if (is.null(M$grd$ObjQ)) Obj2ObjQ(M$grd$Obj,M$grd)->Q else Q=M$grd$ObjQ
  #   X2ObjQ(Q,C$X)->Q
  #   Q$normals=M$grd$Obj$normals # no normals to Q; only for wire3d!
  #   M$grd$ObjQ=Q
  }
  if("dA" %in% what) {
   if (!is.null(M$dA)) q=M$dA else
   {q=E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))$dA}
   if (length(q)!=M$grd$ndof)  {q=E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))$dA}
   M$dA=q;
  }
  if("Quantities" %in% what) {
    M$Curv=Curv(M)
    M$Area=Area(M)
    M$Volume=Volume(M)
  }
  if("SEN" %in% what) {
    C=updateX(M$A,M$grd,M$bas)
    Wb=E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))
    M$SEN=SEN(M$A,M$grd,M$bas,M$Ref,Wb)
  }
  if("Class" %in% what)
  {class(M$bas) = "MemBas"
   class(M$A)   = "MemA"
   class(M$grd) = "MemGrd"}
  if("Time" %in% what){M$Timestamp=timestamp()}
  if("Mask" %in% what) M$bas$mask=double_uv_ind(M$grd$U,M$grd$V)
  return(M)
}


#' save_MemRBC
#'
#' save membrane object M to file
#' erases data from basis M$bas (recomputed when loaded back)
#' @param M The input membrane to be saved
#' @examples
#' data("M4")
#' MMC(M4,10)->M
#' save_MemRBC(M,"Mmmc.rdat")
#' @export
save_MemRBC<-function (M, file = "MemRBC.rdat",
                       reduce_basis = TRUE,
                       thinLA=TRUE, thinA=TRUE,
                       thinbas=TRUE, qs2=FALSE)
{
  if (reduce_basis)
    M$bas$Ylm_u = M$bas$Ylm_v = M$bas$Ylm_uu = M$bas$Ylm_vv = M$bas$Ylm_uv = NULL
  if(thinLA) M$LA=thin_LA(M$LA)
  if(thinA) {attr(M$A,"IM")<-NULL;attr(M$A,"Fit spatial weights")<-NULL}
  if(thinbas) M$bas$IM<-NULL
  M <- StoreParams(M)
  cat(crayon::cyan("stored Params in saved object are:\n"))
  print(unlist(M$Params))
  if (qs2) qs2::qs_save(M, file = paste(file,"qs2",sep=".") )
  else save(M, file = file)
}

#' thin_LA
thin_LA<-function(M){
  M$LA<-lapply(M$LA,function(x){n=names(attributes(x));attributes(x)[!(n %in% c("dim","dimnames"))]<-NULL;x})
  attr(M$A,"IM") <- NULL
  attr(M$A,"Fit spatial weights") <- NULL
  return(M)
}

#' load_MemRBC
#'
#' load membrane object from file and set global data from object
#' @return the membrane object loaded from file with $Params to set globally
#' @examples
#' M<-load_MemRBC("Mmmc.rdat",qs2=FALSE) # may also have saved with qs2=TRUE
#' M
#' @export
load_MemRBC<-function (file = "MemRBC.rdat", qs2 = FALSE)
{
  if (file.exists(file))
    if (qs2)
      M = qs2::qs_read(file)
    else {obj_name = load(file = file);M = get(obj_name)}
  else stop("Membrane file does not exists")
  
  if (is.null(M$bas$Ylm_u))
    M <- FillBasis_MemRBC(M)
  if (!is.null(M$Sample))
    if (!"Id" %in% names(M$Sample))
      M$Sample$Id <- 1
  if (is.null(M$Params)) {
    cat(crayon::red("No $Params found, now filled from environment\n"))
    M <- StoreParams(M)
  }
  else SetParams(M)
  if (is.null(M$bas$Pointymmetry))
    M$bas$Pointymmetry = FALSE
  return(M)
}

# recompute the basis (if deleted on save_MemRBC)
# not needed, see update(M,"Basis")

FillBasis_MemRBC<-function(M)
{  bas=M$bas; L_max=bas$L_max
    u=M$grd$u;v=M$grd$v
    Ai_max=bas$Ai_max
    LM = bas$LM
    L_Ylm=L_Ylm(L_max, u,  v)
    Ylm=L_Ylm$Ylm[,-1] / sqrt(4*pi)
    Ylm_v=Ylm_v(L_max, u,  v, L_Ylm$PLK)[,-1] / sqrt(4*pi)
    Ylm_vv=Ylm_vv(L_max, u,  v, L_Ylm$PLK)[,-1] / sqrt(4*pi)
    L_Y_u=L_Ylm_u(L_max,u,v,L_Ylm$PLK)
    Ylm_u=L_Y_u$Ylm_u[,-1] / sqrt(4*pi)
    Ylm_uu=Ylm_uu(L_max,u,v,L_Y_u$P_T)[,-1] / sqrt(4*pi)
    Ylm_uv=Ylm_uv(L_max,u,v,L_Ylm$PLK,L_Y_u$P_T)[,-1] / sqrt(4*pi)
    l=LM[,1];m=LM[,2]
    M$bas$Ylm=Ylm
    M$bas$Ylm_u=Ylm_u
    M$bas$Ylm_v=Ylm_v
    M$bas$Ylm_uu=Ylm_uu
    M$bas$Ylm_uv=Ylm_uv
    M$bas$Ylm_vv=Ylm_vv
    return(M)
}

#' set_A_to_Ref
#'
#' copy reference shape coefficients to current coeffs.
#' @param M membrane object to copy coefficients into
#' @return modified M2 with new coefficients from M$ARef; spectral order remains unchanged
#' @export
set_A_to_Ref<-function(M)
{ if (is.null(M$ARef)) stop("No reference coefficients in Membrane")
 M1=M;M1$A[]=0
 i=intersect(rownames(M$ARef),rownames(M$A))
 M1$A[i,]=M$ARef[i,];
 rgl::clear3d();plot(M1);rgl::title3d("Reference shape")
 M1$history=append(M1$history,match.call())
 return(M1)
}

#' transplant
#'
#' copy coefficients from one membrane to another
#' @param M1 membrane object to copy coefficients from
#' @param M2 membrane object to copy coefficients into
#' @return modified M2 with new coefficients
#' @examples
#' data("M4")
#' data("S5")
#' transplant(S5,M4),S4mod
#' S4mod
#' @export
transplant<-function(M1,M2,plt=FALSE)
{ i=intersect(rownames(M1$A),rownames(M2$A))
  M2$A[]=0;
  M2$A[i,]=M1$A[i,];
  if (plt) {rgl::clear3d();plot(M2);rgl::title3d("from transplated shape coeffs")}
  M2$history=append(M2$history,"--- transplanted history follows:")
  M2$history=append(M2$history,M1$history)
  M2$history=append(M2$history,"--- transplatned history end.")
  M2$history=append(M2$history,match.call())
  M2$Target=M1$Target
  return(M2)
}

#' StoreParams
#'
#' Store current global membrane parameters in $Param the object; useful for load and save
#' @param M Membrane to store actual parameters like M.C0, M.K_ADE etc. into
#' @examples
#' data("D5")
#' data("D6")
#' SetParams(D5) # set parameters stored in D5 to global variables
#' Storeparams(D6)->D6 # sets parameters from global variables to another object
#' # should be equivalent with
#' # D6$Params <- D5$Params
#' # but the SetParams(D5) sets global parameters for further use
#'
#' @export
StoreParams <- function(M)
{
  M$Params=list(M.C0=M.C0,M.K_b=M.K_b,M.K_ADE=M.K_ADE,M.mu=M.mu,M.Ka=M.Ka,M.a2=M.a2,M.a3=M.a3,M.a4=M.a4,M.b0=M.b0,M.b1=M.b1,M.b2=M.b2)
  return(M)
}

#' SetParams
#'
#' Set parameters from membrane object to global variables
#'  since data(...) cannot load multiple global data, after data() you should call SetParams on the loaded object.
#' @param M Object to take parameters from (stored in M$Params)
#' @export
SetParams<-function(M)
{
  for (i in names(M$Params)) assign(i,M$Params[[i]], envir = .GlobalEnv)
}

#' CheckC0
#'
#' check whether stored C0 equals global M.C0
#' @param M membrane object to check
#' @return -
#' @examples
#' data("M4")
#' CheckC0(M4)
#' @export
Check_C0<-function(M) {
  if(is.null(attr(M$A,"C0"))) stop("no C0 to check in M$A!")
  if (!exists("M.C0")) stop("no M.C0 globally found!")
  if (attr(M$A,"C0")!=M.C0) warning("C0 in M$A differs from global M.C0. Probably not indended!")
  }

# not yet ready, so no export
Replay<-function(M)
{ l0=ls() # actually one must analyse history[[i]] for object names used in apps and names of results.
  for (i in seq_along(M$history)){
    Exec(M$history[[i]]); l=ls()
    obj=setdiff(l,l0)
    if (length(obj)==1) M=get(obj)
    l0=l
  }
  if(exists("M"))  return(M) else stop("could not store last result;\n maybe final result already existed?\n")
}

#' PlotSample
#'
#' plot the curvature-energy-data sampled from MMC.
#' @param M membrane object with $Sample keeping the MMC recorded data
#' @param title title of the plot, placed in separate box
#' @return -
#' @examples
#' data("M4")
#' MMC(M4,1000)->M4mmc
#' PlotSample(M4mmc)
#' @export
PlotSample<-function(M,last=dim(M$Sample)[1],title="MMC sample plot",...)
{
if (!is.null(M$Sample$Id)) col=as.numeric(as.factor(M$Sample$Id)) else col=1
par(mfrow=c(2,2),mar=c(4,3.5,0.3,0.5),oma=c(0,0,1.5,0))
plot(last(M$Sample$Energy,last)/M.Es,pch=".",xlab="",ylab="",col=col,... )
title(ylab = "E", cex.lab = 1,
      line = 2)
title(xlab = "accepted steps", cex.lab = 1,
      line = 2)
plot(last(M$Sample$Curv,last),pch=".",xlab="",ylab="",col=col,...)
title(ylab = "C", cex.lab = 1,
      line = 2)
title(xlab = "accepted steps", cex.lab = 1,
      line = 2)

w = (dim(M$Sample)[1] - last+1) : dim(M$Sample)[1]

plot(Energy/M.Es~Curv,data=M$Sample[w,],pch=".",xlab="",ylab="",col=col,...)

title(ylab = "E", cex.lab = 1,
      line = 2)
title(xlab = "C", cex.lab = 1,
      line = 2)

plot(0,axes=FALSE,col=0,xlab="",ylab="")
title(title)
par(mfrow=c(1,1),mar=c(4,3.5,1,0))
invisible()
}

#' PlotLSeries
#'
#' plot series of shapes from a single object for truncated to lower spectral orders
#' @param M membrane object to plot truncation series from
#' @param nr,nc (=4,=3) number of rows and columns on screen
#' @examples
#' data("M4")
#' M4
#' PlotLSeries(M4,2,3)
#' @export
PlotLSeries<-function(M,nr=4,nc=3,...)
{
  rgl::open3d();
  updateX(M$A,M$grd,M$bas)->C
  SEN(M$A,M$grd,M$bas,M$Ref,E_SCM(M$A,M$grd,M$bas,C))->S
  return(plotLseries(nr,nc,M$A,C,M$grd,M$bas,S=S,Ref=M$Ref,...))
}

#
# fit mono-invaginated form
#   to spherical harmonics in X,Y,Z
#   with initial minimization (in total ~12000 steps) to fullfil constraints.

#' Fit coefficients to a stomatocyte (mono-invaginated) shape
#'
#' @param C0 (=-3) spontaneous curvature to create shape for (usually <0)
#' @param A0,V0 (=140,=100) target values for constraints on area and volume
#' @export
FitStomatocyte_L5<- function(C0=-3,A0=140,V0=100)
{ M.C0<<-C0
  data("M5_Ref", envir = environment())
  S=M5_Ref
  {
  n.grd=25
  grd=MakeGrid_GaussLegendreSimpson(n.grd)
  u=grd$u[,1]
  v=grd$v[1,]
x=u
plot(sin(x)-sin(2*x)/7+sin(4*x)/32+sin(3*x)/8,cos(x)-9*cos(2*x))
X=sin(x)-sin(2*x)/7+sin(4*x)/32+sin(3*x)/8
Z=(cos(x)-9*cos(2*x))/10
M=array(0,c(grd$nu,grd$nv,3))
M[,1,1]=X # initial curve for first v value
M[,1,3]=Z/2
M[,1,2]=0
p=2*pi/(grd$nv) # angle of rotation increment
m=matrix(c(cos(p),-sin(p),0,sin(p),cos(p),0,0,0,1),3,3)
for (j in 2:grd$nv)
  for (i in 1:grd$nu)  {
    M[i,j,]=m%*%M[i,j-1,]
  }
rgl::open3d();rgl::plot3d(x=M[,,1],y=M[,,2],z=M[,,3],aspect=FALSE)
Mflat=matrix(0,prod(dim(M)[1:2]),3)
for (k in 1:3) Mflat[,k]=t(M[,,k])
 bas=MakeBasis_UV(5,t(grd$u),t(grd$v))
 A=FitAlm(Mflat, bas )
 C=updateX(A,grd,bas)
 rgl::clear3d();plot3b(C$X,grd)
}
  E_SCM(A,grd,bas,C,plt=TRUE) -> H2
  (area0=H2$Area) # starting area; area is free, vol0, curv0 rescaled by R0
  (R0=sqrt(area0/4/pi))
  A<-A/R0*sqrt(140/4/pi)
  matplot(A,type="l")
  C=updateX(A,grd,bas  )
  E_SCM(A,grd,bas,C,plt=TRUE) -> H2
  H2$Wb
  H2$Curv
  H2$Area
  H2$Volume
  rgl::clear3d();plot3b(C$X,grd)
  S$A=A
  rgl::clear3d();plot(S)

  S$bas$Target[1:2]=c(A0,V0)

  PSD(S,600,del=1e-7,plt=TRUE)->S1
  S1
  PSD(S1,800,plt=TRUE)->S2

  PSD(S2,del=2e-6, 800,plt=TRUE)->S3
  PSD(S3,del=5e-6, 800,plt=TRUE)->S4

  plot(last(S4$E_PSD/M.Es,1900),type="l")
#  save_MemRBC(S4,"L5-stomatocyte-PSD.rdat")

  PNEM(S4,10000,dt=5e-3)->S4pnem
 # save_MemRBC(S4pnem,"L5-stomatocyte-PSD-PNEM.rdat")
  return(S4pnem)
}


#' Fit coefficients to a prolate ellipsoid
#'
#' Fit an ellipsoid of prolate shape e.g. with Undustick parameters.
#' Global parameters are set for Lipid modeling, i.e. SEN is switched off.
#' @param C0 (=2.562) Undustick spontaneous curvature
#' @param Vp (=85.66) Undustick volume, which is reduced volume 0.55
#' @return mebrane object of prolate shape to be optimized further, e.g. for Undustick
#' @examples
#' FitProlate_Ellipsoid_L5() -> M
#' plot(M)
#'
#' @export
FitProlate_Ellipsoid_L5<- function(C0=2.562,V0=85.66)
{ data("M5_Ref", envir = environment())
  SetParams(M5_Ref)
  M.C0<<-C0
  M.K_ADE<<-M.mu<<-M.Ka<<-0
  S<-M5_Ref
{
  n.grd=5*6+2
  grd=MakeGrid_GaussLegendreSimpson(n.grd)
  u=grd$u[,1]
  v=grd$v[1,]
  x=u
  X=sin(x)/2.2
  Z=cos(x)*8.1
  M=array(0,c(grd$nu,grd$nv,3))
  M[,1,1]=X # initial curve for first v value
  M[,1,3]=Z/2
  M[,1,2]=0
  p=2*pi/(grd$nv) # angle of rotation increment
  m=matrix(c(cos(p),-sin(p),0,sin(p),cos(p),0,0,0,1),3,3)
  for (j in 2:grd$nv)
    for (i in 1:grd$nu)  {
      M[i,j,]=m%*%M[i,j-1,]
    }
  rgl::open3d();rgl::plot3d(x=M[,,1],y=M[,,2],z=M[,,3],aspect=FALSE)
  Mflat=matrix(0,prod(dim(M)[1:2]),3)
  for (k in 1:3) Mflat[,k]=t(M[,,k])
  bas=MakeBasis_UV(5,t(grd$u),t(grd$v))
  bas$Target[2]=V0

  A=FitAlm(Mflat, bas )
  C=updateX(A,grd,bas)
  rgl::clear3d();plot3b(C$X,grd)
}

E_SCM(A,grd,bas,C,plt=TRUE) -> H2
(area0=H2$Area) # starting area; area is free, vol0, curv0 rescaled by R0
(R0=sqrt(area0/4/pi))
A<-A/R0*2.85

matplot(A,type="l")
C=updateX(A,grd,bas  )
E_SCM(A,grd,bas,C,plt=TRUE) -> H2
H2$Wb
H2$Curv
H2$Area
H2$Volume

r0=sqrt(140/4/pi) # here target area is 140
v0=4/3*pi*r0^3
v0_red=0.55*v0

rgl::clear3d();plot3b(C$X,grd)
S$A=A
rgl::clear3d();plot(S)
S$bas$Target[2] <- V0

save_MemRBC(S,"L5-prolate-Fit.rdat")
print(unlist(Quantities(S)))

PSD(S,1600,del=4e-8,plt=TRUE,pltfreq=100)->S1
PSD(S1,400,1e-7,plt=TRUE,pltfreq=100)->S2
PSD(S2,400,2e-7,plt=TRUE,pltfreq=100)->S2
PSD(S2,400,3e-7,plt=TRUE,pltfreq=100)->S3

plot(last(S3$E_PSD/M.Es,1900),type="l",ylab=expression(E[PSD]))
save_MemRBC(S3,"L5-prolate-PSD.rdat")

#PNEM(S3,1000,dt=1e-6,viscosity = 1000,rho=0.1)->S3
#save_MemRBC(S3pnem,"L5-prolate-PSD-PNEM.rdat")
return(S3)
}

#' FitProlate_L5
#'
#' Fit coefficients from a stick with semi-spherical caps to initiate Undustick.
#' @param C0 spontaneous curvature for the model (default 2.562 for Undustick)
#' @param V0 real volume (default 85.66 for Undustick)
#' @examples
#' # this sets also membrane parameters such that SEN is switched off
#' i.e. M.K_ADE=M.mu=M.Ka=0 # no SEN and ADE-term for Undustick
#' M.C0=2.562 # reference C0 value, see H0 in https://zenodo.org/records/13627757
#'
#' FitProlate_L5() -> P5
#' P5
#' save_MemRBC(P5,"Undustick-prototype.rdat")
#' PlotLSeries(P5,3,2)
#' MMC(P5,100000)->P5mmc
#' PlotLSeries(P5,3,2)
#' SDRC(P5,1000,plot=TRUE,pltfreq=1)->P5_sdrc
#' P5_sdrc
#' @export
FitProlate_L5<- function(C0=2.562,V0=85.66,filter_strength=10,no_minim=TRUE)
{
  data("M5_Ref", envir = environment())

  M.C0<<-C0
  M.K_ADE<<-M.mu<<-M.Ka<<-0
  S=M5_Ref
{
  n.grd=5*6+2 # usual grid paraneter for L=5
  grd=MakeGrid_GaussLegendreSimpson(n.grd)
  u=grd$u[,1]
  v=grd$v[1,]
  x=u
  X=Z=x # to initialize
  X[x<pi/4]=sin(2*x[x<pi/4])*0.1
  X[x>3*pi/4]=sin(pi/2+2*(x[x>3*pi/4]-3*pi/4))*0.1

  Z[x>pi/4 & x<3*pi/4] = 0.5 - (x[x>pi/4 & x<3*pi/4]-pi/4)/(pi/2)
  X[x>pi/4 & x<3*pi/4] = 0.1

  Z[x<pi/4]=0.1*cos(2*x[x<pi/4])+0.5
  Z[x>3*pi/4]=0.1*cos(pi/2+2*(x[x>3*pi/4]-3*pi/4))-0.5
 # plot(X,Z)

  M=array(0,c(grd$nu,grd$nv,3))
  M[,1,1]=X # initial curve for first v value
  M[,1,3]=Z
  M[,1,2]=0
  p=2*pi/(grd$nv) # angle of rotation increment
  m=matrix(c(cos(p),-sin(p),0,sin(p),cos(p),0,0,0,1),3,3)
  for (j in 2:grd$nv)
    for (i in 1:grd$nu)  {
      M[i,j,]=m%*%M[i,j-1,]
    }
  rgl::open3d();rgl::plot3d(x=M[,,1],y=M[,,2],z=M[,,3],aspect=FALSE)
  Mflat=matrix(0,prod(dim(M)[1:2]),3)
  for (k in 1:3) Mflat[,k]=t(M[,,k])
  bas=MakeBasis_UV(5,t(grd$u),t(grd$v))
  bas$Target[2]=V0

  A=FitAlm(Mflat, bas )
  C=updateX(A,grd,bas)
  rgl::clear3d();plot3b(C$X,grd)
}

E_SCM(A,grd,bas,C,plt=TRUE) -> H2
(area0=H2$Area) # starting area; area is free, vol0, curv0 rescaled by R0
(R0=sqrt(area0/4/pi))
A<-A/R0*2.85

matplot(A,type="l")
C=updateX(A,grd,bas  )
E_SCM(A,grd,bas,C,plt=TRUE) -> H2
H2$Wb
H2$Curv
H2$Area
H2$Volume

r0=sqrt(140/4/pi)
v0=4/3*pi*r0^3
v0_red=0.55*v0

rgl::clear3d();plot3b(C$X,grd)
S$A=A
rgl::clear3d();plot(S)
S$bas$Target[2] <- V0

save_MemRBC(S,"L5-prolate-Fit-b.rdat")
print(unlist(Quantities(S)))
rgl::open3d();plot(S);
if(!no_minim) PSD(S,100,del=4e-8,plt=TRUE,pltfreq=100,filter_strength = filter_strength)->S1 else S1<-S
S1$history=match.call()
StoreParams(S1)->S1
S1$comment="From FitProlate_L5, only 100 steps PSD minimization."
return(S1)
}


#
# fit invaginated form (Stomatocyte)
#   to spherical harmonics in X,Y,Z
#

#' FitStomatocyte_L
#'
#' Fit a stomatocyte shape from data, spectral order L
#' @param L spectral order of output membrane
#' @param C0 (=-3) C0-value to use for initial minimization
#' @param V0 (=100) target value of volume
#' @param dt (=1e-3) time step in PNEM minimizer
#' @return membrane object
#' @examples
#' FitStomatocyte_L(L=4,C0=0.5)->M
#' plot(M)
#' @export
FitStomatocyte_L<- function(L=7,C0=-3,V0=100,dt=1e-3)
{ M.C0<<-C0
  data("D5", envir = environment())
  D5$bas$Target["Volume"]=V0
  t0=proc.time()
  cl=match.call()
  {
  n.grd=(L+1)*5+2
  grd=MakeGrid_GaussLegendreSimpson(n.grd)
  u=grd$u[,1]
  v=grd$v[1,]
x=u
plot(sin(x)-sin(2*x)/7+sin(4*x)/32+sin(3*x)/8,cos(x)-9*cos(2*x))
X=sin(x)-sin(2*x)/7+sin(4*x)/32+sin(3*x)/8
Z=(cos(x)-9*cos(2*x))/10
M=array(0,c(grd$nu,grd$nv,3))
M[,1,1]=X # initial curve for first v value
M[,1,3]=Z/2
M[,1,2]=0
p=2*pi/(grd$nv) # angle of rotation increment
m=matrix(c(cos(p),-sin(p),0,sin(p),cos(p),0,0,0,1),3,3)
for (j in 2:grd$nv)
  for (i in 1:grd$nu)  {
    M[i,j,]=m%*%M[i,j-1,]
  }
rgl::open3d();rgl::plot3d(x=M[,,1],y=M[,,2],z=M[,,3],aspect=FALSE)
Mflat=matrix(0,prod(dim(M)[1:2]),3)
for (k in 1:3) Mflat[,k]=t(M[,,k])
 bas=MakeBasis_UV(L,t(grd$u),t(grd$v))

 A=FitAlm(Mflat, bas )
 C=updateX(A,grd,bas)
 rgl::clear3d();plot3b(C$X,grd)
}
  E_SCM(A,grd,bas,C,plt=TRUE) -> H2
  (area0=H2$Area) # starting area; area is free, vol0, curv0 rescaled by R0
  (R0=sqrt(area0/4/pi))
  A<-A/R0*sqrt(140/4/pi)
  matplot(A,type="l")
  C=updateX(A,grd,bas  )
  E_SCM(A,grd,bas,C,plt=TRUE) -> H2
  H2$Wb
  H2$Curv
  H2$Area
  H2$Volume
  rgl::clear3d();plot3b(C$X,grd)
  bas=MakeBasis_UV(L,grd$U,grd$V)
  bas$Target["Volume"]=V0
  A=LM2A(A,bas)
  M=structure(list(A=A,grd=grd,bas=bas),class="MemRBC")
  ARef=A;
  ARef[]=0;
  ARef[1:3,]=D5$Ref$ARef[1:3,]
  ARef<-LM2A(ARef,bas)
  M$ARef=ARef
  update(M,"Ref")->M
  rgl::clear3d();plot(M)
  PSD(M,300,del=1e-7,plt=TRUE)->M
  PSD(M,400,plt=TRUE)->M
  PSD(M,del=2e-6, 400,plt=TRUE)->M
  PSD(M,del=5e-6, 400,plt=TRUE)->M

  plot(last(M$E_PSD/M.Es,1000),type="l",ylab="E",xlab="step")
  M$comment="Stomatocyte shape from fit, weakly minimized"
  M$history=cl
  M$proc_time=proc.time()-t0
#  save_MemRBC(M,paste("L",L,"-stomatocyte-PSD.rdat",sep=""))
# no real dynamics wanted
  PNEM(M,5000,dt=dt,zero_Av=TRUE)->M
  M$comment="Stomatocyte shape from fit, stronger minimized"
  M$history=cl
  M$proc_time=proc.time()-t0
  StoreParams(M)->M
#  save_MemRBC(M,paste("L",L,"-stomatocyte-PSD-PNEM.rdat",sep=""))
  return(M)
}


#' FitDiscocyte_L5
#'
#' Fit a discocyte from data, spectral order L=5
#' @param C0 C0-value to use for initial minimization
#' @return membrane object
#' @examples
#' FitDiscocyte_L5(C0=0.5)->M
#' plot(M)
#' @export
FitDiscocyte_L5<- function(C0=0,V0=100)
{
  data("D5")
  SetParams(D5)
  M.C0=C0
  S=D5
  {
    n.grd=25
    grd=MakeGrid_GaussLegendreSimpson(n.grd)
    u=grd$u[,1]
    v=grd$v[1,]
    x=u
    plot(sin(x),cos(x)-cos(3*x)/1.7)
    X=sin(x)
    Z=cos(x)-cos(3*x)/1.7
    M=array(0,c(grd$nu,grd$nv,3))
    M[,1,1]=X*1.1 # initial curve for first v value
    M[,1,3]=Z/3.5
    M[,1,2]=0
    p=2*pi/(grd$nv) # angle of rotation increment
    m=matrix(c(cos(p),-sin(p),0,sin(p),cos(p),0,0,0,1),3,3)
    for (j in 2:grd$nv)
      for (i in 1:grd$nu)  {
        M[i,j,]=m%*%M[i,j-1,]
      }
    rgl::open3d();rgl::plot3d(x=M[,,1],y=M[,,2],z=M[,,3],aspect=FALSE)
    Mflat=matrix(0,prod(dim(M)[1:2]),3)
    for (k in 1:3) Mflat[,k]=t(M[,,k])
    bas=MakeBasis_UV(5,t(grd$u),t(grd$v))
    A=FitAlm(Mflat, bas )
    C=updateX(A,grd,bas)
    rgl::clear3d();plot3b(C$X,grd)
  }

  E_SCM(A,grd,bas,C,plt=TRUE) -> H2
  (area0=H2$Area) # starting area; area is free, vol0, curv0 rescaled by R0
  (R0=sqrt(area0/4/pi))
  A<-A/R0*sqrt(140/4/pi)
  matplot(A,type="l")
  C=updateX(A,grd,bas  )
  E_SCM(A,grd,bas,C,plt=TRUE) -> H2
  H2$Wb
  H2$Curv
  H2$Area
  H2$Volume
  rgl::clear3d();plot3b(C$X,grd)
  S$A=A
  rgl::clear3d();plot(S)
  S$comment="created from fit to discoid"
  S$bas$Target["Volume"]=V0
  S$Ref$v
  PSD(S,600,del=1e-7,plt=TRUE)->S1
  PSD(S1,200,del=1e-7,plt=TRUE)->S1
  S1$Ref$v
  PSD(S1,400,plt=TRUE)->S1
  PNEM(S1,2000,dt=1e-2)->S1pnem
  plot(S1pnem$E_total_PNEM,type="l")
  PNEM(S1pnem,1000,dt=1e-2,new_Av = FALSE)->S1pnem
  plot(last(S1pnem$E_total_PNEM,2000),type="l")
  save_MemRBC(S1pnem,"L5-discocyte-PSD-PNEM.rdat")
  return(S1pnem)
}

#' @export
Filter_1_per_lm<-function(A,bas)
{ # filter for largest magnitude entry per (X,Y,Z)
  A1=A;
  for (i in 1:bas$Ai_max)
  {w=which.max(abs(A[i,])); A1[i,setdiff(1:3,w)]<-0}
  return(A1)
}

non_sparsity_cost<-function(A,bas,thresh=1e-12)
{ # count non-zero entries of magnitude above thresh
  # normalized >0, <=1 ; all entries > thresh -> cost=1
  return(sum(abs(A)>thresh)/length(A[]))
}

#' Sparsify coefficient matrix
#'
#' zero out coefficients upto one (largest magnitude) entry per row (X,Y,Z).
#' For diagnosis, the energy values before and after Sparsify are printed and returned.
#' @param M membrane object to sparsify coefficients
#' @return membrane object, with diagnosis results in $Sparse; drop_norm reports the L2-norm of all zeroed out coefficients.
#' @examples
#' data("M_stomatocyte_L12")
#' Sparsify(M_stomatocyte_L12) -> M_sparse
#' plot(M_sparse)
#' M_sparse$SparseE
#' @export
Sparsify<- function(M)
{ M$Sparse=list("before" = Energy(M))
  Filter_1_per_lm(M$A,M$bas) -> A1
  drop_norm=pracma::Norm(M$A-A1)
  M$A <- A1
  M$Sparse[["after"]] <- Energy(M)
  M$Sparse[["drop_norm"]]<-drop_norm
  M$Sparse[["drop_norm_%"]]<-drop_norm/pracma::Norm(M$A)*100.0

  M$history=append(M$history,match.call())
  M$comment="Sparsified"
  M$last_App_called="Sparsify"
  return(M)
}

#' PlotRef
#'
#' plot SEN reference shape from membrane object
#' @param M membrane object to plot SEN reference shape from
#' @param Aref optional coefficients, eg if the objects ARef is not working
#' @examples
#' data("M4")
#' PlotRef(M4)
#' @export
PlotRef<-function(M,ARef=M$Ref$ARef)
{
  rgl::clear3d()
  plot3b(updateX(ARef,M$grd,M$bas)$X,M$grd)
  rgl::title3d("Reference shape")
}


#' PlotPNEM
#'
#' plot recorded energies from PNEM runs
#' @param M membrane object to plot energies from
#' @examples
#' data("M4")
#' PNEM(M4,1000)->M
#' PNEMVM(M,1000)->M1
#' PlotPNEM(M1)
#' @export
PlotPNEM<-function(M,from=1,to=length(M$C_PNEM))
{ p=par()$mfrow
  par(mfrow=c(2,1),mar=c(1,3,0,0),oma=c(3,1,0.5,0.5))
  if(is.null(M$E_kin_PNEM)) {message("No PNEM data to plot"); return();}
  if (!is.null(M$type_PNEM)) col=as.factor(M$type_PNEM) else col=1
  print(table(col))
  if (is.factor(col)) {
     plot(M$E_kin_PNEM[from:to]/M.Es,pch=".",col=col,xlab="",ylab="",axes=FALSE);

     legend("topright",legend = levels(col),pch=20,col=1:nlevels(col),cex=0.6)} else
    { plot(M$E_kin_PNEM[from:to]/M.Es,type="l",col=1,xlab="",ylab="",axes=FALSE);
      legend("topright",legend="PNEM",pch=20, col=1,cex=0.6 )}
    axis(3,labels = NA,tick = NA)
    axis(2,padj=0.2)
    title(ylab=expression(E[kin]),line=2 )
    box()
  if (is.factor(col)) {
    plot(M$E_total_PNEM[from:to]/M.Es,pch=".",col=col,xlab="",ylab="");
#    legend("topright",legend = levels(col),pch=20,col=1:nlevels(col),cex=0.6)
    } else
    {  plot(M$E_total_PNEM/M.Es,type="l",col=1,xlab="",ylab="");
#   legend("topright",legend="PNEM",pch=20, col=1,cex=0.6 )
    }
    title(xlab="steps",line=1,outer=TRUE,adj=0.57)
    title(ylab=expression(E[total]),line=2)

  par(mfrow=p)
  invisible()
}



#' PlotPSD
#'
#' plot recorded energies from PSD runs
#' @param M membrane object to plot PSD recorded energies from
#' @param lastE (length M$E_PSD) number of trailing energy points to plot
#' @param lastC (=lastE) number of trailing curvature points to plot
#' @examples
#' data("M4")
#' PSD(M4,1000) -> M
#' PSDC(M, curv=Curv(M), 1000) -> M1
#' PlotPSD(M1, lastE=750)
#' @export
PlotPSD<-function(M,lastE=length(M$E_PSD),lastC=lastE)
{
 if(lastE>length(M$E_PSD)) lastE=length(M$E_PSD)
 if(lastC>length(M$C_PSD)) lastC=length(M$C_PSD)
 wE=(length(M$E_PSD)-lastE+1):length(M$E_PSD)
 wC=(length(M$C_PSD)-lastC+1):length(M$C_PSD)
 p=par()$mfrow
 par(mfrow=c(2,1),mar=c(1,3,0,0),oma=c(3,1,0.5,0.5))
 if(is.null(M$E_PSD)) {message("No PSD data to plot"); return();}
 plot(M$E_PSD[wE]/M.Es,type="l",col=1,xlab="",ylab="",axes=FALSE);
 legend("topright",legend="PSD",pch=20, col=1,cex=0.6 )
 axis(3,labels = NA,tick = NA)
 axis(2,padj=0.2)
 title(ylab=expression(E[PSD]),line=2 )
 box()
 plot(M$C_PSD[wC],type="l",col=1,xlab="",ylab="");
 title(xlab="steps",line=1,outer=TRUE,adj=0.57)
 title(ylab=expression(C[PSD]),line=2)
 par(mfrow=p)
 invisible()
}

#' data_MemRBC
#'
#' load membrane data from package and set parameters
#'   - in contrast to data("name"), $Params from object are set to global params.
#'   - equivalent to : data("name");SetParams(name)
#'
#' @param name  name of data (in quotes) to load
#'
#' @export
data_MemRBC<-function(name)
{
  assign(name,load_MemRBC(paste(system.file(package = "MemBRC"),"data/",name,".rda",sep="")),envir = parent.env(environment()))
}

#' Rewind
#'
#' rewind coefficients in M$A by n records from the list of recorded coeffs M$LA
#'  - practically you do the same with
#'  M$A=M$LA[[length(M$LA)-n]]
#' @param M membrane object with list $LA of recorded coefficients
#' @param n (=1) how many records to rewind?
#' @return membrane object with M$A set to LA[[lengh(LA)-n]]
#' @export
Rewind<-function(M,n=1)
{ if (n<length(M$LA))
   M$A=M$LA[[length(M$LA)-n]] else message("no rewind; n too large\n")
  return(M)
}


#' ReduceM1
#'
#'
#' Reduce basis in M to leading cos v, sin v terms, i.e. |m|<2 in Y_lm.
#' Requires update(M,"Ref"), if SEN is used (i.e. M.Ka and/or M.mu not zero.
#' It is recommended to save a reduced model with save_MemRBC(...,reduce_basis=FALSE),
#' to keep the basis in reduced form in the file.
#' For loading, you must then use load_MemRBC(..., unreduce_basis=FALSE) to not overwrite the reduced basis.
#'
#' @param M membrane object to change
#' @return modified M with reduced basis, but $Ref unchanged.
#' @examples
#' data("US_L10_de_novo"); SetParams(US_L10_de_novo)
#' ReduceM1(US_L10_de_novo) -> L10fast
#' # should reduce from 120 to 30 rows in bas$A
#' #
#' # not needed, since paraneters switch off SEN::
#' update(L10,what="Ref")->L10fast
#' plot(L10fast)
#'
#' MMC(L10fast, 20000)->M10fast_mmc
#' PlotSample(L10fast_mmc,last=20000)
#'
#' @export
ReduceM1<-function(M)
{
  bas=M$bas
  w=which(abs(bas$LM[,2])>1)
  bas$A=bas$A[-w,];  bas$Ylm=bas$Ylm[,-w]
  bas$Ylm_u=bas$Ylm_u[,-w];  bas$Ylm_v=bas$Ylm_v[,-w]
  bas$Ylm_uu=bas$Ylm_uu[,-w];  bas$Ylm_uv=bas$Ylm_uv[,-w];  bas$Ylm_vv=bas$Ylm_vv[,-w]
  bas$LM=bas$LM[-w,];  bas$l=bas$LM[,1];  bas$m=bas$LM[,2]
  bas$Lset=unique(bas$l);  bas$Mset=unique(bas$m)
  bas$Ai_max=dim(bas$A)[1]
  bas$G.tk=bas$G.tk[-w]
  bas$comment=paste("reduced basis to |m|<2;  ", bas$comment)
  M$A=M$A[-w,]
  M$bas=bas
  return(M)
}

#' GEMINI_Obj_Obj_Intersect
#' @description
#' this code from GEMINI intends to locate spatial intersections of 3D-objects
#' @param S1, S2 MemRBC objects; S2=S1 is used for self-intersections
#' @return boolean, TRUE if intersection is probable
#' @export
GEMINI_Obj_Obj_Intersect <- function (S1, S2 = S1)
{
  GEMINI_tri_tri_intersect_3d <- function(T1, T2) {
    epsilon <- 1e-09
    cross_prod <- function(a, b) {
      c(a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] *
          b[3], a[1] * b[2] - a[2] * b[1])
    }
    v1 <- T2[, 2] - T2[, 1]
    v2 <- T2[, 3] - T2[, 1]
    N2 <- cross_prod(v1, v2)
    d2 <- -sum(N2 * T2[, 1])
    dists1 <- apply(T1, 2, function(p) sum(N2 * p) + d2)
    if (all(dists1 > epsilon) || all(dists1 < -epsilon))
      return(FALSE)
    v1 <- T1[, 2] - T1[, 1]
    v2 <- T1[, 3] - T1[, 1]
    N1 <- cross_prod(v1, v2)
    d1 <- -sum(N1 * T1[, 1])
    dists2 <- apply(T2, 2, function(p) sum(N1 * p) + d1)
    if (all(dists2 > epsilon) || all(dists2 < -epsilon))
      return(FALSE)
    D <- cross_prod(N1, N2)
    if (sum(D^2) < epsilon) {
      return(FALSE)
    }
    compute_interval <- function(Tri, dists, D) {
      projs <- apply(Tri, 2, function(p) sum(D * p))
      intervals <- c()
      for (i in 1:3) {
        j <- (i%%3) + 1
        s1 <- dists[i]
        s2 <- dists[j]
        if ((s1 > epsilon && s2 < -epsilon) || (s1 <
                                                -epsilon && s2 > epsilon)) {
          t <- s1/(s1 - s2)
          p <- projs[i] + t * (projs[j] - projs[i])
          intervals <- c(intervals, p)
        }
        else if (abs(s1) <= epsilon) {
          intervals <- c(intervals, projs[i])
        }
        else if (abs(s2) <= epsilon) {
          intervals <- c(intervals, projs[j])
        }
      }
      if (length(intervals) == 0)
        return(NULL)
      return(range(intervals))
    }
    r1 <- compute_interval(T1, dists1, D)
    r2 <- compute_interval(T2, dists2, D)
    if (is.null(r1) || is.null(r2))
      return(FALSE)
    if (r1[2] < r2[1] || r2[2] < r1[1]) {
      return(FALSE)
    }
    return(TRUE)
  }
  cand = GEMINI_rvcg_kdtree_candidates(S1, S2)
  for (i in 1:dim(cand)[1]) if (GEMINI_tri_tri_intersect_3d(S1$vb[1:3,
                                                           S1$it[, cand[i, 1]]], S2$vb[1:3, S2$it[, cand[i, 2]]]))
    return(TRUE)
  return(FALSE)
}

#' GEMINI_Intersect_Mem_Mem
#' @description
#' translate membranes to 3D objects and check for intersection
#'
#' @param M1,M2 MemRBC data
#' @export
GEMINI_Intersect_Mem_Mem<-function(M1,M2=M1)
{return(GEMINI_Obj_Obj_Intersect( update(M1,"Obj")$grd$Obj,update(M2,"Obj")$grd$Obj))
}

#' @export
"+.MemRBC"<-function(m1,m2)
{m1$A=m1$A+m2$A;return(m1)}
#' @export
"-.MemRBC"<-function(m1,m2)
{m1$A=m1$A-m2$A;return(m1)}

#' @export
"*.MemRBC"<-function(a,b)
{if (is(a,"MemRBC") & is.numeric(b)) {a$A=a$A*b;return(a)}
 if (is(b,"MemRBC") & is.numeric(a)) {b$A=b$A*a;return(b)}
 warning("not correct types - return NULL")
return(NULL)
}

GEMINI_rvcg_kdtree_candidates <- function(meshA, meshB) {

  # 1. Calculate a safe search radius
  # A safe r is the maximum distance from a triangle centroid to its furthest vertex
  # For simplicity, we can use the average edge length or a small user-defined epsilon
  max_edge <- max(Rvcg::vcgMeshres(meshB)$edgelength)
  # 2. Query the KD-tree
  # vcgKDtree finds indices of vertices in meshB closest to vertices in meshA
  # We use the 'radius' search to catch all potential overlaps
  kd_search <- Rvcg::vcgKDtree(meshB, meshA, k=1)

  # kd_search$index contains the indices of vertices in Mesh B
  # that are within 'max_edge' of vertices in Mesh A.

  # 3. Identify Candidate Triangles
  # Get indices of vertices in Mesh A that had at least one neighbor in B
  hits_in_A <- which(sapply(kd_search$index, length) > 0)

  if (length(hits_in_A) == 0) return(NULL)

  # Map these "close" vertices of A to their Triangles
  # A triangle is a candidate if ANY of its vertices are "hits"
  triA_candidates <- which(colSums(matrix(meshA$it %in% hits_in_A, nrow=3)) > 0)

  # For each candidate triangle in A, find the closest triangle in B
  # We use vcgClost here because it uses a fast AABB-tree/KD-tree internally
  # to find the EXACT closest face.
  final_pairs <- list()

  # Vectorized closest face search for the candidate centroids
  centroidsA <- (meshA$vb[1:3, meshA$it[1, triA_candidates]] +
                   meshA$vb[1:3, meshA$it[2, triA_candidates]] +
                   meshA$vb[1:3, meshA$it[3, triA_candidates]]) / 3

  closest_in_B <- Rvcg::vcgClost(t(centroidsA), meshB)

  candidates <- data.frame(
    triA = triA_candidates,
    triB = closest_in_B$faceptr
  )

  return(unique(candidates))
}
#' @export
print.MemA<-function(A)
{
   cat("MemA Coefficient object of size ",dim(A)[1],"x",dim(A)[2],"\n")
   cat("C0",attr(A,"C0"),"\t") 
   cat("A0",attr(A,"A0"),"\t") 
   cat("V0",attr(A,"V0"),"\n")
   cat("E",attr(A,"E"),"\n")
   cat("attributes:\n",names(attributes(A)),"\n")
}
#' @export
print.MemGrd<-function(G)
{
  cat("MemGrid object of ndof ",G$ndof,"\n")
  cat(G$comment,"\n")
  cat("Angles arrays dimension:",dim(G$u),"\n")
  cat("range u",range(G$u)," \nrange v",range(G$v),"\n")
}

#' @export
print.MemC<-function(C){
  cat("MemC coordinate object \n")
  cat("X    : ",dim(C$X),"\n")
  cat("Xu ...: ",dim(C$X_u),"\n")
}

#' @export
print.MemC_X<-function(C){
  cat("MemC coordinate only object \n")
  cat("X    : ",dim(C$X),"\n")
}

#' @export
print.MemBas<-function(B){
  cat("MemBas basis functions object \n")
  cat("Yml      : ",dim(B$Ylm),"\n")
  cat("Yml_u ...: ",dim(B$Ylm_u),"\n")
}


#' @export
UpdateM<-function(M)
{class(M$bas) = "MemBas"
 class(M$A)   = "MemA"
 class(M$grd) = "MemGrd"
 return(M)
}
