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

#' îmage from 2d data
#' @description
#' 2D plot of color coded variable defined on the shape.
#' x and y correspond to spherical angles
#' @param x Membrane to use for shape plot
#' @param q a 2D field to project onto shape, eg. M$SEN$alpha
#' @param main (=expression(alpha)) plot title
#' @param ... parameters to imag()
#' @return -
#' @export image.MemRBC
#' @export
image.MemRBC<-function(x,q=x$SEN$alpha,main=expression(alpha),...)
{
  imag(q,x$grd,main=main,...)
  return()
}

#' plot
#' @description
#' 3D plot of membrane shape
#' @param x Membrane to use for shape plot
#' @param wire (boolean) FALSE if  to suppress mesh quads plotted
#' @param wire_col (="black") 
#' @param ... further plotting options, e.g. col="green"
#' @export plot.MemRBC
#' @export
plot.MemRBC <- function(x,wire=TRUE,wire_col="black",...){
    R<-x
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
#' @examplesIf exists("L_Ylm")
#' #' # before using get:data_ZENODO, make sure the requested files really exist
#' get_data_ZENODO(L=c("ss42denovo_pnem.rda","ss42denovo_mmc.rda"),local=TRUE)
#' load_MemRBC("data/ss42denovo_pnem.rda")->ss42denovo_pnem
#' load_MemRBC("data/ss42denovo_mmc.rda")->ss42denovo_mmc
#' dQ <- PlotDiff(ss42denovo_pnem, ss42denovo_mmc, col="red",On3d=TRUE)
#' # PNEM removed two spicules to reduce volume, 
#' # while mmc reduced volume in favour of shrinking globally along Z
#' 
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

#' Plot a stability analysis visually on a 3x3 canvas of 3d-plots
#' @description
#' run PlotStabGallery on a membrane object with a Stab in it from MemStab()
#' @param M The input membrane with initial data and reference
#' @param which_n vector of Ids of eigenvalues, 1:9 is default for lowest nine
#' @param wire (=FALSE) to have wireframe plot on top of shade3d()
#' @param plt_scale (=0.3) scale of eigenvectors displacement
#' @param sharedMouse (boolean) if FALSE, each subplot can be controlled by mouse individually
#' @param ... for further plotting parameters, like setting semi-transparency alpha=0.5; col1 and col2 can be set here (blue,red are default)
#' @return none
#' @examplesIf exists("L_Ylm")
#' M <- MakeStandardRBC(L=5)
#' plot(M)
#' M <- PSD(M,nsteps=100,del=1e-6,)
#' plot(M)
#' MemStab(M)->M
#' PlotStabGallery(M)
#' @export
PlotStabGallery <- function(M, which_n=1:9, plt_scale=0.3, sharedMouse=TRUE,wire=FALSE, ...)
{ if (is.null(M$Stab)) stop("No stability data in MemRBC object")
  rgl::open3d();rgl::mfrow3d(3,3,sharedMouse = sharedMouse);
  cat("vectors for last eigenvalues used from =\n ")
  print(round(M$Stab$EigH$values,6))
  for (i in which_n)
  {if (i>1) rgl::next3d();
    plotStab(M,plt_n=i,wire=wire,plt_scale=plt_scale,...)
    N=length(M$Stab$EigH$values)
    rgl::title3d(paste(i,"ev",round(M$Stab$EigH$values[N+1-i],6)))
  }
}

plotStab<-function(M,plt_n=1,plt_scale=0.3,col1="blue",col2="red",alpha1=0.5,alpha2=0.5,wire=FALSE)
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
#' @description
#' analyse membrane shape stability in terms of eigenvectors of Hessian
#' @param M MemRBC object
#' @param mc.cores (=4) for parallel Hessian
#' @param plt (=FALSE) control plot, if wanted (TRUE)
#' @param plt_mode (=1) starting mode index, 1 for lowest eigenvalue
#' @param plt_scale (=0.5) for second shapes shift in coeffs along Hessians eigenvector
#' @param serial (=FALSE) TRUE, if serial Hessian computation is wanted.   
#' @return MemRBC object, M$Stab has eigensystem of Hessian of energy without constraints
#' @export
MemStab <- function(M, mc.cores = 4, plt=FALSE, plt_mode=1, plt_scale=0.5, serial = FALSE)
  { t0=proc.time()
    if(is.null(M$proc_time)) M$proc_time=0
    cl=match.call()
    if (serial) H=FullModelHessian(M$A,M$grd,M$bas,M$Ref) else
      H=FullModelHessian_Par(M$A,M$grd,M$bas,M$Ref, del = 5e-06,
                             Mem_mc.cores = mc.cores, timing = TRUE,
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
#'@description
#' compute the mass matrix for PNEM
#' @param M The input membrane with initial data and reference
#' @param rho area density, default 1
#' @return massmatrix in spectral space
#' @examplesIf exists("L_Ylm")
#' #' data("M4",package = "MemRBC")
#' mass <- massmatrix(M4,rho=2)
#' image(mass)
#' 
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
#' @description
#' report some quantities from a membrane object
#' @param M The input membrane with initial data and reference
#' @return quantities Area, Volume, curvature Curv and bending energy Wb as named vector
#' @examplesIf exists("L_Ylm")
#' data("M1",package = "MemRBC")
#' Quantities(M1)
#' @export
Quantities<-function(M)
{r=E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))[c("Area","Volume","Curv","Wb")]
 names(r)=c("Area","Volume","Curv","Wb")
 return(unlist(r))
}


#' Energy
#'@description
#' report  energy values for SCM and SEN from a membrane object
#' @param M The input membrane with initial data and reference
#' @return vector of energies Wb (bending), Es (stress-shear), E (potetntial energy), Ekin (kinetic energy, optional, if a PNEM was run)
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

#' compute area
#' @param M MemRBC object
#' @return area
#' @export
Area<-function(M)
  {return(E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))$Area)}

#' compute Volume
#' @param M MemRBC object
#' @return volume
#' @export
Volume<-function(M)
  {return(E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))$Volume)}

#' compute integrated curvature integ(k1+k2)
#' @param M MemRBC object
#' @return total curvature
#' @export
Curv<-function(M)
  {return(E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))$Curv)}

# re-orient shape by prinicpal axis of initertia
#' MemPCA
#'@description
#' rotate Membrane coefficients to coordinates of principle axes of inertia.
#' Remark: the reference SEN in $Ref remains unchanged.
#' @param M the input membrane to be rotated
#' @param WX (=1) spatial weights for the fit of coefficients, try sin(grd$U) to downweight pole influence
#' @return  membrane object with new coefficients of PCA-rotated spatial coordinates.
#' @examplesIf exists("L_Ylm")
#' data("M4",envir=environment())
#' MemPCA(M4)->M2
#' plot(M2)
#' 
#' @export
MemPCA<-function(M,WX=rep(1,M$grd$ndof))
  { cl=match.call()
    updateX_only(M$A,M$grd,M$bas)->M$C
    princomp(M$C$X)$scores->X
    A=FitAlm_Tikhonov(X = X,bas=M$bas,lambda=0,WX=WX)
    M$A<-A
    M$history<-append(M$history,cl)
    return(M)
  }

# open two empty screens, keep devs as  M.scr1 and M.scr2 in MemRBC_env
#' two_screens3d
#'@description
#' open two screens for 3D-plots by two_draw3d()
#' @param x,y dimensions in pixels
#' @return variables MemRBC_env$M.scr1 and MemRBC_env$M.scr2 are set to device ids
#' @examples
#' two_screens3d(x=650,y=300)
#' @export
two_screens3d<-function(x=400,y=400){
  if (MemRBC_env$M.scr2 %in% rgl::rgl.dev.list()) {rgl::open3d();MemRBC_env$scr2=rgl::cur3d();rgl::par3d(windowRect=c(x-30,30,2*x-30,y+30))}; # draw right first
  if (MemRBC_env$M.scr1 %in% rgl::rgl.dev.list()) {rgl::open3d();MemRBC_env$scr1=rgl::cur3d();rgl::par3d(windowRect=c(1,30,x+1,y+30))}; # then left (for 30 pixel overlap)
  }

# draw stress and shear; re-open screens if not opened (not in rgl.dev.list())
#' two_draw3d
#'@description
#' plot membrane object with area stress alpha and shear beta as color
#' @param A coefficients to use for plot
#' @param M the input membrane with $grid and $bas matching A
#' @param cont (=FALSE) draw gridlines u=const, v=const
#' @param x,y dimensions in pixels'
#' @param title title text for rgl::title3d()'
#' @return -
#' @examplesIf exists("L_Ylm")
#' data("M4",package = "MemRBC",envir=environment())
#'  # perturb original coefficients
#' A <- pertA_Gauss(M4$A, M4$bas, sd=0.2)
#' two_screens(); two_draw3d(A, M4)
#' @export
two_draw3d<-function(A,M,cont=FALSE,title="",x=400,y=400) # requires S (Stretches) as global variable
  { grd=M$grd;
    C<-updateX(A,grd,M$bas)
    Wb<-E_SCM(A,grd,M$bas,C)
    S<-SEN(A,grd,M$bas,M$Ref,Wb)
    X2Obj(grd$Obj,C$X)->O
    Rvcg::vcgClean(O,sel=1:7,silent = TRUE)->O
    if (MemRBC_env$M.scr2 %in% rgl::rgl.dev.list()) {rgl::set3d(MemRBC_env$M.scr2);rgl::clear3d();} else { assign("M.scr2",rgl::open3d(),envir=MemRBC_env);rgl::par3d(windowRect=c(x-30,30,2*x-30,y+30)); }

    imag.obj.colorbar(O,f=S$beta,clr = FALSE,par=FALSE,specular="black"); rgl::title3d(paste("beta",title))
if(cont){    rgl::contourLines3d(O,grd$U,nlev=15,lwd=2)
  rgl::contourLines3d(O,grd$v,levels = pracma::linspace(0,2*pi,16)[-16], lwd=2)
}
    if (MemRBC_env$M.scr1 %in% rgl::rgl.dev.list()) {rgl::set3d(MemRBC_env$M.scr1);rgl::clear3d();} else { assign("M.scr1",rgl::open3d(),envir=MemRBC_env);rgl::par3d(windowRect=c(x-30,30,2*x-30,y+30)); }
    imag.obj.colorbar(O,f=S$alpha,clr=FALSE,par=FALSE,specular="black");rgl::title3d(paste("alpha",title))
if(cont){    rgl::contourLines3d(O,grd$U,nlev=15,lwd=1)
  rgl::contourLines3d(O,grd$v,nlevels=15, levels = pracma::linspace(0,2*pi,16)[-16],lwd=1)
  }
}

#' update
#' @description
#' update data in membrane object
#' @param object the input membrane to be updated
#' @param what vector of character from "dA", "Quantities", "Basis", "Grid", "Ref", "curv", "SCM", "X", "SEN", "Mask", "Time", "Mask"
#' @param n for what="Grid": new number of grid points
#' @param L for what="Basis": spectral order L
#' @param ... not used
#' @return updated membrane object
#' @examplesIf exists("L_Ylm")
#' get_data_ZENODO(L="M_stomatocyte_L12.rda")
#' data("M_stomatocyte_L12",package = "MemRBC",envir=environment())
#' update(M_stomatocyte_L12,"X")->M
#' plot3d(M$X, aspect=FALSE)
#' # make lower spectral order membrane from M
#' update(M, what=c("Grid","Basis","Ref"), n=30, L=8)->L8
#' two_screens3d(); two_draw3d(L8$A, L8)
#' rgl::open3d(); plot(L8)
#' Quantities(L8)
#' Energy(L8)
#' @export update.MemRBC
#' @export
update.MemRBC <- function(object, what=c("dA","Quantities","X"),n=(L+1)*5+2,L=5,...)
{ M<-object
  if("Grid" %in% what){
    grd=MakeGrid_GaussLegendreSimpson(n*L);
    M$grd=grd;
    M$mass=NULL;M$Ref=NULL
    what=c(what,"Basis")
  }
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
    if(is.null(M$bas$Target)) {
      message("$Target created from scratch with 2 standard constraints")
      Target=c(140,100);names(Target)=c("Area","Volume")
    } else  Target=M$bas$Target
    bas=MakeBasis_UV(L,M$grd$U,M$grd$V);
    bas$Target=Target
    bas$Nc=length(Target)
    Ain=M$A;
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
    if (is.null(M$ARef)) {message("take as ref A[1:3,] from M$A, not M$ARef");ARef=M$A[1:3,]} else ARef=M$ARef
    M$ARef=ARef
    B=MakeBasis_UV(L_max=1,M$grd$U,M$grd$V)
    M$Ref=Ref4CauchyGreen(M$ARef,M$grd,B)
    M$mass=NULL
    M$Av=NULL
    class(M$Ref)="MemRef"
    class(M$A)="MemA"
    class(M$ARef)="MemA"
  }
  if ("X" %in% what) M$X=updateX_only(M$A,M$grd,M$bas)$X
  if ("Obj" %in% what)
  { #  the following update
  X2Obj(M$grd$Obj,updateX_only(M$A,M$grd,M$bas)$X)->M$grd$Obj
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
  if("Coor" %in% what) {
    C=updateX(M$A,M$grd,M$bas)
    M$C=C;class(M$C)="MemC"
  }
  if("Class" %in% what)
  {if (!is.null(M$Bas)) class(M$bas) = "MemBas"
   if (!is.null(M$A))   class(M$A)   = "MemA"
   if (!is.null(M$grd)) class(M$grd) = "MemGrd"
   if (!is.null(M$Ref)) class(M$Ref) = "MemRef"
   }
  
  if("Time" %in% what){M$Timestamp=timestamp()}
  if("Mask" %in% what) M$bas$mask=double_uv_ind(M$grd$U,M$grd$V)
  if(is.null(M$bas)) message("Membrane has no basis!\nUse what=\"Basis\" in update(,L=Lmax)")
  if(is.null(M$grd)) message("Membrane has no Grid!\nUse what=\"Grid\" in update(...,n=ngrid)")
  
  return(M)
}


#' save_MemRBC
#'@description
#' save membrane object M to file
#' erases data from basis M$bas (some recomputed when loaded back)
#' LA entries get class "MemA" (for older files)
#' @param M The input membrane to be saved
#' @param file  (="MemRBC.rdat") the filename for saving
#' @param reduce_basis (=TRUE) removes derivative data of basis functions
#' @param thinLA (=TRUE) removes additional attributes in LA coeff records
#' @param thinA (=TRUE) removes additional attributes in coeffs A
#' @param thinbas (=TRUE) removes IM matrix from basis
#' @param qs2 (=FALSE) for faster save in .qs2 file use this. In this case, ".qs2" is appended on filename.
#' @return name of saved file
#' @examplesIf exists("L_Ylm")
#' data("M4",package = "MemRBC",envir=environment())
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
  M$LA=lapply(M$LA,function(x) class(x)="MemA")
  if(thinA) {attr(M$A,"IM")<-NULL;attr(M$A,"Fit spatial weights")<-NULL}
  if(thinbas) M$bas$IM<-NULL
  M <- StoreParams(M)
  M <- update(M,"Class")
  cat(crayon::cyan("stored Params in saved object are:\n"))
  print(unlist(M$Params))
  if (qs2) qs2::qs_save(M, file = paste(file,"qs2",sep=".") )
  else save(M, file = file)
  return(file)
}


#' erase several attributes from coefficients
#' Only dim and dimnames are preserved.
#' If inversion matrix IM is attribute to M$A, it will be erased.
#' @param M MemRBC membrane object
#' @return MemRBC object, with fewer attributes
#' @examples
#' data(M4)
#' thin_LA(M4)->M4l
#' @export
thin_LA<-function(M){
  M$LA<-lapply(M$LA,function(x){n=names(attributes(x));attributes(x)[!(n %in% c("dim","dimnames"))]<-NULL;x})
  attr(M$A,"IM") <- NULL
  attr(M$A,"Fit spatial weights") <- NULL
  return(M)
}

#' load_MemRBC
#' @description
#' load membrane object from file and set global parameter data from object
#' @param file filename to load data from
#' @param qs2 TRUE if file has ending .qs2, FALSE otherwise, or set by user
#' @return the membrane object loaded from file with $Params to set globally
#' @examplesIf exists("L_Ylm")
#' M<-load_MemRBC("Mmmc.rdat",qs2=FALSE) # may also have saved with qs2=TRUE
#' M
#' @export
load_MemRBC<-function (file = "MemRBC.rdat", qs2 = length(grep("qs2$",file)==1))
{ 
  if (file.exists(file))
    if (qs2)
      M = qs2::qs_read(file)
    else {obj_name = load(file = file);M = get(obj_name)}
  else stop("Membrane file does not exists")
  
  if (is.null(M$bas$Ylm_u))
    M <- FillBasis_MemRBC(M)
  if (!is.null(M$Sample))
    if (!"Id" %in% names(M$Sample)) # correct missing Id data
      M$Sample$Id <- 1
  if (is.null(M$Params)) {
    cat(crayon::red("No $Params found, now filled from environment\n"))
    M <- StoreParams(M)
  } else SetParams(M)
  if (is.null(M$bas$Pointymmetry)) M$bas$Pointymmetry = FALSE # assume no point symmetry
  return(M)
}

# recompute the basis (if deleted on save_MemRBC)
# not needed, see update(M,"Basis")
FillBasis_MemRBC<-function(M)
{  bas=M$bas; L_max=bas$L_max
    u=M$grd$u;v=M$grd$v
    Ai_max=bas$Ai_max
    LM = bas$LM
    L_Ylm_=L_Ylm(L_max, u,  v)
    M$bas$Ylm=L_Ylm_$Ylm[,-1] / sqrt(4*pi)
    M$bas$Ylm_v=Ylm_v(L_max, u,  v, L_Ylm_$PLK)[,-1] / sqrt(4*pi)
    M$bas$Ylm_vv=Ylm_vv(L_max, u,  v, L_Ylm_$PLK)[,-1] / sqrt(4*pi)
    L_Y_u_=L_Ylm_u(L_max,u,v,L_Ylm_$PLK)
    M$bas$Ylm_u=L_Y_u_$Ylm_u[,-1] / sqrt(4*pi)
    M$bas$Ylm_uu=Ylm_uu(L_max,u,v,L_Y_u_$P_T)[,-1] / sqrt(4*pi)
    M$bas$Ylm_uv=Ylm_uv(L_max,u,v,L_Ylm_$PLK,L_Y_u_$P_T)[,-1] / sqrt(4*pi)
    l=LM[,1];m=LM[,2]
    return(M)
}

#' set_A_to_Ref
#'@description
#' copy reference shape coefficients $ARef to current coeffs $A.
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
#'@description
#' copy coefficients from one membrane to another
#' @param M1 membrane object to copy coefficients from
#' @param M2 membrane object to copy coefficients into
#' @param plt (=FALSE) TRUE for a 3d plot of resulting shape
#' @return modified M2 with new coefficients
#' @examplesIf exists("L_Ylm")
#' data("M4",package = "MemRBC",envir=environment())
#' get_data_ZENODO("SS.rda")
#' data("SS",package = "MemRBC",envir=environment())
#' transplant(SS,M4)->SSmod
#' SSmod
#' 
#' @export
transplant <- function(M1,M2,plt=FALSE)
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
#'@description
#' Store current MemRB_env membrane parameters in $Param of the object; useful for load and save
#' @param M Membrane to store actual parameters like MemRBC_env$M.0, M.K_ADE etc. into
#' @return MemRBC object, with Params from `MemRBC_env` 
#' @export
StoreParams <- function(M)
{
  M$Params=list(M.C0=MemRBC::MemRBC_env$M.C0,M.K_b=MemRBC::MemRBC_env$M.K_b,
                M.K_ADE=MemRBC::MemRBC_env$M.K_ADE,
                M.mu=MemRBC::MemRBC_env$M.mu,M.Ka=MemRBC::MemRBC_env$M.Ka,
                M.a2=MemRBC::MemRBC_env$M.a2,M.a3=MemRBC::MemRBC_env$M.a3,
                M.a4=MemRBC::MemRBC_env$M.a4,M.b0=MemRBC::MemRBC_env$M.b0,
                M.b1=MemRBC::MemRBC_env$M.b1,M.b2=MemRBC::MemRBC_env$M.b2)
  return(M)
}

#' SetParams
#'@description
#' Set parameters from membrane object to MemRBC_env variables
#'  since data(...) cannot load multiple global data, after data() you should call SetParams on the loaded object.
#' @param M Object to take parameters from (stored in M$Params)
#' @return -
#' @export
SetParams<-function(M)
{
  for (i in names(M$Params)) assign(i,M$Params[[i]], envir = MemRBC_env)
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
#'@description
#' plot the curvature-energy-data sampled from MMC.
#' @param M membrane object with $Sample keeping the MMC recorded data
#' @param last : integer, how many data from the tail should be plotted
#' @param title title of the plot, placed in separate box
#' @param ... further plot parameters
#' @return -
#' @examplesIf exists("L_Ylm")
#' data("M4",package = "MemRBC",envir=environment())
#' MMC(M4,1000)->M4mmc
#' PlotSample(M4mmc)
#' @export
PlotSample<-function(M,last=dim(M$Sample)[1],title="MMC sample plot",...)
{
if (!is.null(M$Sample$Id)) col=as.numeric(as.factor(M$Sample$Id)) else col=1
par(mfrow=c(2,2),mar=c(4,3.5,0.3,0.5),oma=c(0,0,1.5,0))
plot(last(M$Sample$Energy,last)/MemRBC_env$M.Es,pch=".",xlab="",ylab="",col=col,... )
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

plot(Energy/MemRBC_env$M.Es~Curv,data=M$Sample[w,],pch=".",xlab="",ylab="",col=col,...)

title(ylab = "E", cex.lab = 1,
      line = 2)
title(xlab = "C", cex.lab = 1,
      line = 2)

plot(0,axes=FALSE,col=0,xlab="",ylab="")
title(title)
par(mfrow=c(1,1),mar=c(4,3.5,1,0))
}

#' PlotLSeries
#' @description
#' plot series of shapes from a single object.
#' This visualizes truncation effects on shape and energy.
#' @param M membrane object to plot truncation series from
#' @param nr,nc (=4,=3) number of rows and columns of sub-plots on screen
#' @param ... parameters plot function 
#' @return -
#' @examplesIf exists("L_Ylm")
#' data("M4",package = "MemRBC",envir=environment())
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
#' FitStomatocyte_L5
#' @description
#' Fit coefficients to a stomatocyte (mono-invaginated) shape
#'
#' @param C0 (=-3) spontaneous curvature to create shape for (usually <0)
#' @param A0,V0 (=140,=100) target values for constraints on area and volume
#' @return MemRBC object
#' @export
FitStomatocyte_L5<- function(C0=-3,A0=140,V0=100)
{ MemRBC_env$M.C0<-C0
  data("M4",package = "MemRBC",envir=environment())
  S=M4
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

  plot(last(S4$E_PSD/MemRBC_env$MEs,1900),type="l")
#  save_MemRBC(S4,"L5-stomatocyte-PSD.rdat")

  PNEM(S4,10000,dt=5e-3)->S4pnem
 # save_MemRBC(S4pnem,"L5-stomatocyte-PSD-PNEM.rdat")
  return(S4pnem)
}


#' Fit coefficients to a prolate ellipsoid
#'@description
#' Fit an ellipsoid of prolate shape e.g. with Undustick parameters.
#' Parameters are set for Lipid modeling, i.e. SEN is switched 
#' off (mu=Ka=0).
#' This is a long-runner over 2800 PSD steps
#' @param C0 (=2.562) Undustick spontaneous curvature
#' @param V0 (=85.66) Undustick volume, which is reduced volume 0.55
#' @return MemRBC object of prolate shape to be optimized further, e.g. for Undustick
#' @export
FitProlate_Ellipsoid_L5<- function(C0=2.562,V0=85.66)
{ 
  data("M4", package = "MemRBC",envir=environment())
  SetParams(M4)
  MemRBC_env$M.C0<-C0
  MemRBC_env$M.K_ADE<-MemRBC_env$M.mu<-MemRBC_env$M.Ka<-0
  S<-M4
  StoreParams(M4)->M4
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

#matplot(A,type="l")
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
S$bas$Target["Volume"] <- V0

#save_MemRBC(S,"L5-prolate-Fit",qs2=TRUE)
print(Quantities(S))

PSD(S,300,del=4e-8,plt=TRUE,pltfreq=100)->S1
PSD(S1,400,2e-8,plt=TRUE,pltfreq=100)->S2
SDRC(S2,100,plt=TRUE,del_min=1e-8,del_cons = 0.1, cons_tol = 5,max_iter = 30)->S2r

SDRC(S2r,100,plt=TRUE)->S2r

PSD(S2,400,1e-8,plt=TRUE,pltfreq=100)->S2a
PSD(S2a,400,3e-7,plt=TRUE,pltfreq=100)->S3

#plot(last(S3$E_PSD/MemRBC_env$MEs,1900),type="l",ylab=expression(E[PSD]))
#save_MemRBC(S3,"L5-prolate-PSD",qs2=TRUE)

PNEM(S,1000,dt=1e-6,viscosity = 1000,rho=0.1)->S3
#save_MemRBC(S3pnem,"L5-prolate-PSD-PNEM.rdat")
return(S3)
}

#
# fit invaginated form (Stomatocyte)
#   to spherical harmonics in X,Y,Z
#

#' FitStomatocyte_L
#'@description
#' Fit a stomatocyte shape from data, spectral order L
#' @param L spectral order of output membrane
#' @param C0 (=-3) C0-value to use for initial minimization
#' @param V0 (=100) target value of volume
#' @param dt (=1e-3) time step in PNEM minimizer
#' @return MemRBC object
#' @examplesIf exists("L_Ylm")
#' FitStomatocyte_L(L=4,C0=0.5)->M
#' plot(M)
#' @export
FitStomatocyte_L<- function(L=7,C0=-3,V0=100,dt=1e-3)
{ MemRBC_env$M.C0<-C0
  data("M4" ,package = "MemRBC",envir=environment())
  M4$bas$Target["Volume"]=V0
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
  ARef[1:3,]=M4$Ref$ARef[1:3,]
  ARef<-LM2A(ARef,bas)
  M$ARef=ARef
  update(M,"Ref")->M
  rgl::clear3d();plot(M)
  PSD(M,300,del=1e-7,plt=TRUE)->M
  PSD(M,400,plt=TRUE)->M
  PSD(M,del=2e-6, 400,plt=TRUE)->M
  PSD(M,del=5e-6, 400,plt=TRUE)->M

  plot(last(M$E_PSD/MemRBC_env$M.Es,1000),type="l",ylab="E",xlab="step")
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
#'@description
#' Fit a discocyte from data, spectral order L=5
#' @param C0 (=0) C0-value (spontaneous curvature) to use for initial minimization
#' @param V0 (=100) target volume to use
#' @return MemRBC object
#' @export
FitDiscocyte_L5<- function(C0=0,V0=100)
{
  data("M4",package = "MemRBC",envir=environment())
  SetParams(M4)
  bas=M$bas
  MemRBC_env$M.0<-C0
  S=M4
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
  S$Params[["M.C0"]]=0
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

#' filter for largest coeff. per lm row
#' @param A,bas coefficient matrix ans basis
#' @return sparsified coefficient matrix
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
#'@description
#' zero out coefficients upto one (largest magnitude) entry per row (X,Y,Z).
#' For diagnosis, the energy values before and after Sparsify and other diagnostics are returned.
#' @param M membrane object to sparsify coefficients
#' @return membrane object, with diagnosis results in $Sparse; drop_norm reports the L2-norm of all zeroed out coefficients.
#' @examplesIf exists("L_Ylm")
#' data("M4",package = "MemRBC",envir=environment())
#' Sparsify(M4) -> M_sparse
#' 
#' M_sparse$Sparse
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
#' @description
#' plot SEN reference shape from membrane object
#' @param M membrane object to plot SEN reference shape from
#' @param ARef optional coefficients, eg if the objects ARef is not working
#' @return -
#' @examplesIf exists("L_Ylm")
#' data("M4",package = "MemRBC",envir=environment())
#' PlotRef(M4)
#' @export
PlotRef<-function(M,ARef=M$ARef)
{
  rgl::clear3d()
  plot3b(updateX(ARef,M$grd,M$bas)$X,M$grd)
  rgl::title3d("Reference shape")
}

#' PlotPNEM
#'@description
#' plot recorded energies from PNEM runs
#' @param M membrane object to plot energies from
#' @param from,to (=1), (=length) range of indices to plot data sets for 
#' @return -
#' @examplesIf exists("L_Ylm")
#' data("M4",package = "MemRBC",envir=environment())
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
     plot(M$E_kin_PNEM[from:to]/MemRBC_env$M.Es,pch=".",col=col,xlab="",ylab="",axes=FALSE);

     legend("topright",legend = levels(col),pch=20,col=1:nlevels(col),cex=0.6)} else
    { plot(M$E_kin_PNEM[from:to]/MemRBC_env$M.Es,type="l",col=1,xlab="",ylab="",axes=FALSE);
      legend("topright",legend="PNEM",pch=20, col=1,cex=0.6 )}
    axis(3,labels = NA,tick = NA)
    axis(2,padj=0.2)
    title(ylab=expression(E[kin]),line=2 )
    box()
  if (is.factor(col)) {
    plot(M$E_total_PNEM[from:to]/MemRBC_env$M.Es,pch=".",col=col,xlab="",ylab="");
#    legend("topright",legend = levels(col),pch=20,col=1:nlevels(col),cex=0.6)
    } else
    {  plot(M$E_total_PNEM/MemRBC_env$M.Es,type="l",col=1,xlab="",ylab="");
#   legend("topright",legend="PNEM",pch=20, col=1,cex=0.6 )
    }
    title(xlab="steps",line=1,outer=TRUE,adj=0.57)
    title(ylab=expression(E[total]),line=2)
  par(mfrow=p)
  invisible()
}



#' PlotPSD
#'@description
#' plot recorded energies from PSD runs
#' @param M membrane object to plot PSD recorded energies from
#' @param lastE (length M$E_PSD) number of trailing energy points to plot
#' @param lastC (=lastE) number of trailing curvature points to plot
#' @return -
#' @examplesIf exists("L_Ylm")
#' data("M4",package = "MemRBC",envir=environment())
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
 plot(M$E_PSD[wE]/MemRBC_env$M.Es,type="l",col=1,xlab="",ylab="",axes=FALSE);
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
#'@description
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

#' Rewind to set coefficient M$A by n records back from end of list of recorded coeffs M$LA
#' - practically you do the same with
#' M$A=M$LA[[length(M$LA)-n]]
#' ATTENTION: while M$A usually holds the lasrt coefficients from an App
#'  last entries in LA hold older coeffs, in particular if LAfreq>1 in the last App call 
#' @param M membrane object with list $LA of recorded coefficients
#' @param n (=1) how many records to rewind?
#' @return membrane object with M$A set to `LA[[lengh(LA)-n]]`
#' @export
Rewind<-function(M, n=1)
{ if (n<length(M$LA))
   M$A=M$LA[[length(M$LA)-n]] else message("no rewind; n too large\n")
  return(M)
}


#' ReduceM1
#'
#' @description
#' Reduce basis in M to leading cos v, sin v terms, i.e. |m|<2 in Y_lm.
#' Requires update(M,"Ref"), if SEN is used (i.e. M.Ka and/or M.mu not zero.
#' It is recommended to save a reduced model with save_MemRBC(...,reduce_basis=FALSE),
#' to keep the basis in reduced form in the file.
#' For loading, you must then use load_MemRBC(..., unreduce_basis=FALSE) to not overwrite the reduced basis.
#'  
#' The reduction on X and Y coeficients makes a Axisymmetric parameterization.
#'  Take care: The reduction is also on Z coeffs, 
#'  which drastically reduces the representable shape space.
#' @param M membrane object to change
#' @return modified M with reduced basis, but $Ref unchanged.
#' @examplesIf exists("L_Ylm")
#' #download data with 
#' if (interactive()) {
#' get_data_ZENODO(L="U17R.rda",local=TRUE)
#' load_MemRBC("data/U17R.rda")->U17R;
#' SetParams(U17R)
#' ReduceM1(U17R) -> L10fast
#' plot(L10fast)
#'  MMC(L10fast, 20000,C0=2.562)->M10fast_mmc
#'  PlotSample(L10fast_mmc,last=20000)
#' }
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
#' @param S1,S2 MemRBC objects; S2=S1 is used for self-intersections
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
#' compute from membrane data the 3D objects and check for intersection
#' (experimental, not checked intensively).
#' @param M1,M2 (M2=M1 for self-intersection) MemRBC membrane objects
#' @return TRUE, if M1 and M2 are intersecting; else: FALSE
#' @export
GEMINI_Intersect_Mem_Mem<-function(M1,M2=M1)
{return(GEMINI_Obj_Obj_Intersect( update(M1,"Obj")$grd$Obj,update(M2,"Obj")$grd$Obj))
}

#' add coefficients
#' @param m1,m2 MemRBC objects to add
#' @return result MemRBC object
#' @export
"+.MemRBC"<-function(m1,m2)
{ rn=intersect(rownames(m1$A),rownames(m2$A))
  m1$A[rn,]=m1$A[rn,]+m2$A[rn,];
  message("added in m1 ",rn,"\n");return(m1)}

#' subtract coefficients
#' @param m1,m2 MemRBC objects for M1$A-M2$A
#' @return result MemRBC object with coeffs difference
#' @export
"-.MemRBC"<-function(m1,m2)
{rn=intersect(rownames(m1$A),rownames(m2$A))
 m1$A[rn,]=m1$A[rn,]-m2$A[rn,];
 message("subtracted from m1 ",rn,"\n");return(m1)}

#' scale coefficients by a scalar
#' @param a,b scalar and MemRBC object for a*b$A (or the other way round)
#' @return result MemRBC object with scaled coeffs
#' @export
"*.MemRBC"<-function(a,b)
{if ( is.numeric(b)) {a$A=a$A*b;return(a)}
 if ( is.numeric(a)) {b$A=b$A*a;return(b)}
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



#' print
#' @description
#' print metadata from membrane object class MemRBC
#' @param x Membrane to print data from
#' @param ... not used
#' @export print.MemRBC
#' @export
print.MemRBC<-function(x,...){
  cat("MembraneRBC Object:","\n")
  cat("object size ",object.size(x))
  cat("ndof_grid=",x$grd$ndof," [",x$grd$nu,"x",x$grd$nv,"]\n")
  if(!is.null(x$bas$Ylm_u))
  {cat("L=",x$bas$L_max,", coeff. rows",x$bas$Ai_max,"\n")
    cat(x$bas$Nc," constraints target:\n");
    print(x$bas$Target)
    cat("Membranes recent coefficients A computed at C0=",attr(x$A,"C0"),"\n")
    if (!rlang::is_named(x$bas$Target)) message("Object constraint $Target is unnamed!")
    cat("quantities for current coeffs A: Area",Area(x)," Volume",Volume(x)," Curv",Curv(x),"\n")
  } else cat(crayon::cyan("The Basis in this Membrane is reduced to X!\nUse update(M,\"Basis\") to recreate.\n"))
  if (!is.null(x$LA)) cat("List of coeffs LA length:",length(x$LA),"\n")
  
  if (!is.null(x$kT)) cat("MMC kT:",x$kT,"\n")
  
  if (!is.null(x$visc)) if(length(c(x$visc))==1) cat("PNEM viscosity:",x$visc,"\n")
  else cat("PNEM viscosity (head of diag of matrix)",diag(x$visc)[1:6],"\n")
  if (!is.null(x$rho)) cat("PNEM mass density:",x$rho,"\n")
  if (!is.null(x$Sample)) cat("MMC Samples:",dim(x$Sample)[1],"\n")
  if (!is.null(x$proc_time)) cat("Membrane has overall process time ",x$proc_time,"\n") else cat("Membrane has no process_time recorded.\n")
  if (!is.null(x$comment)) cat("Membrane comment ",x$comment,"\n") else cat("Membrane has no comment.\n")
  cat("Membrane has a history of length ",length(x$history),"\n")
  if(is.null(x$Params)) cat(crayon::red("No Parameters - inject by StoreParams(M)->M\n")) else {
    ul=unlist(x$Params); cat("Params:\n"); print(ul)
  }
}


#' print some coefficient metadata class MemA
#' @param x coefficient matrix object class "MemA"
#' @param ... not needed
#' @return -
#' @export print.MemA
#' @export
print.MemA<-function(x,...)
{
   cat("MemA Coefficient object of size ",dim(x)[1],"x",dim(x)[2],"\n")
   cat("C0",attr(x,"C0"),"\t") 
   cat("A0",attr(x,"A0"),"\t") 
   cat("V0",attr(x,"V0"),"\n")
   cat("E",attr(x,"E"),"\n")
   cat("attributes:\n",names(attributes(x)),"\n")
   invisible(x)
}

#' print some grid metadata
#' @param x grid object
#' @param ... not needed
#' @return -
#' @export print.MemGrd
#' @export
print.MemGrd<-function(x,...)
{
  cat("MemGrid object of ndof ",x$ndof,"\n")
  cat(x$comment,"\n")
  cat("Angles arrays dimension:",dim(x$u),"\n")
  cat("range u",range(x$u)," \nrange v",range(x$v),"\n")
  invisible(x)
}

#' print some coordinate object metadata class MemC
#' @param x coordinate object of class "MemC", like from `updateX()`
#' @param ... not needed
#' @return -
#' @export print.MemC
#' @export
print.MemC<-function(x,...){
  cat("MemC coordinate object \n")
  cat("X    : ",dim(x$X),"\n")
  cat("Xu ...: ",dim(x$X_u),"\n")
  invisible(x)
}

#' print some coordinate object metadata, where only X is set (class MemC_X)
#' @param x grid object of class "MemC_X"
#' @param ... not needed
#' @return -
#' @export print.MemC_X
#' @export
print.MemC_X<-function(x,...){
  cat("MemC coordinate only object \n")
  cat("X    : ",dim(x$X),"\n")
  invisible(x)
}

#' print some basis metadata class MemBas
#' @param x basis object of class "MemBas"
#' @param ... not needed
#' @return -
#' @export print.MemBas
#' @export
print.MemBas<-function(x,...){
  cat("MemBas basis functions object \n")
  cat("Yml      : ",dim(x$Ylm),"\n")
  cat("Yml_u ...: ",dim(x$Ylm_u),"\n")
  invisible(x)
}

#' print some Reference metadata class MemRef
#' @param x grid object of class MemRef as create e.g. by `MakeRef()`
#' @param ... not needed
#' @return -
#' @export print.MemRef
#' @export
print.MemRef<-function(x,...)
{
  cat("MemRef memory size",object.size(x),"\n")
  cat("length tgi ",length(x$tgi),"\n")
  cat("length giPrep ",length(x$giPrep),"\n")
 invisible(x)  
}

#' print some SCM Energy metadata class MemESCM
#' @param x Energy object class "MemESCM"
#' @param ... not needed
#' @return -
#' @export print.MemESCM
#' @export
print.MemESCM<-function(x,...)
{
  cat("MemESCM energy Wb ",x$Wb,"\n")
  cat("length F ",length(x$FF),"\n")
  cat("length dA ",length(x$dA),"\n")
  invisible(x)
}



#' Plot coefficients by matplot(). 
#' Colors black, red, green for X, Y, Z coefficients.
#' @param M membrane to plot spectral coefficients M$A from
#' @param bar (=TRUE) for L-levels to be color-plotted (bar at y=0)
#' @param scale_up (=TRUE) plot A scaled-up by sqrt(bas$G.tk)
#' @param ... parameters handed to matplot
#' @return -
#' @examplesIf exists("L_Ylm")
#' data("M4",package = "MemRBC",envir=environment());
#' SetParams(M4)
#' PlotA(M4,scale_up=FALSE)
#' @export
PlotA<-function(M, bar=TRUE, scale_up=TRUE,...)
{  plotA_l(M$A, M$bas, bar=bar, scale_up=scale_up, ...) }


# no more