# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2025 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, S. (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# Frickenhaus S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#


#' Constraint Newton Minimizer
#'
#' run CNM on an existing membrane object to minimize energy under exact constraints on area and volume
#' You may include curvature constraint as well.
#' Temperature is controlled by parameter kT.
#' @param M The input membrane with initial data and reference
#' @param nsteps (=5) number of Newton steps to be performed
#' @param del (=0.3) update factor delte (<1)
#' @param thresh (=1e-4)threshold for pseudoinverse to suppress almost zero eigenvalues in matrix inversion
#' @param plt (=TRUE) for plotting two figures with alpha/beta stress-shear parameters
#' @param Mtol_Newton (=1e-4) stopping criterion for norm of gradient
#' @param LAfreq (=1) storage frequency into list of coefficients LA
#' @param cluster (=FALSE) for parallel Hessian using a cluster with ncores processes
#' @param ncores (=4) for parallel Hessian with ncores
#' @return membrane object with updated data from CNM:
#' @return LA: list of recorded coefficients A
#' @return A: last coefficients
#' @return history: history of App-calls that created the result
#' @return E_CNM: vector of nsteps energy values
#' @return CNMiter: number of iterations, incl. previous calls
#' @return proc_time: aggregated processing time of all Apps called before with this object
#' @examples
#' M <- MakeStandardRBC(L=5)
#' plot(M)
#' M <- MMC(M,nsteps=10000,kT=0.00411,kTfac=0.99,kTfreq=100)
#' plot(M)
#' M
#' @export
CNM<-function(M,nsteps=5,del=0.3, diag.reg=0.0,
              LAfreq = 1,
              Mtol_Newton=1e-4,
              thresh=1e-4,
              ncores=4, plt=TRUE,filter_delta=ID,cluster=FALSE) # returns a modified MemRBC
{ t0=proc.time()
  cl=match.call()
  if(is.null(M$proc_time)) M$proc_time<-0
  if(!exists("M.Rcpp")) stop("Cannot process - probably load_param_MemRBC has not been called.")
  M.Rcpp<<-TRUE
  M.Rcpp_ncores<<-ncores
  E0=1000
  if (is.null(M$LA)) LA=list(M$A) else LA=M$LA# to return the iterated solutions
  if (is.null(M$CNMiter)) M$CNMiter=0
  E_CNM=rep(0,nsteps)

  A=M$A;grd=M$grd;bas=M$bas; Ref1=M$Ref
  if( is.null(M$Lambda)) M$Lambda=rep(1,M$bas$Nc)

  Lambda=M$Lambda
  if(plt) two_screens3d()
  for (iter in (1:nsteps )) { # bas changed to bas everywher
    if(iter==1)tictoc::tic()
    A0=A
    Lambda0=Lambda
    if(iter==1)tictoc::tic()
    if (cluster)
       {H<-FullModelHessian_Par(A,grd,bas,Ref1,del=1e-6,Mem_mc.cores = ncores,timing=iter==1) } else
         {H<-Hessian_FullModel(A,grd,bas,Ref1,1e-6,ncores = ncores)} # del is forward finite difference parameter
    diag(H$H) = diag(H$H) * ( 1 + diag.reg)
    if(iter==1){cat("parallel Hessian ");tictoc::toc()}
    Cons_RHS <- ConsRHS(H$h2,bas)

    Hp=ConstraintHessian(H,bas,Lambda)
    eigen(H$H)$values -> Eig.H
    {cat("CNM:",iter,":Eig.H:",sort(Eig.H)[1:8],"\n")}
    eigen(Hp)$values -> Eig
    cat("CNM:",iter,":Eig.Hp:",sort(Eig)[1:8],"\n")
    # NEWTON STEP
    # indices of Lambda in full solution vector
    LambdaI=(dim(H$H)[1]+1):(dim(H$H)[1]+bas$Nc)
    ConsJ=ConsJacobian(H$g2,bas)
    #H$G has full energy gradient; add Lagr. terms here
    GG <- c(H$G) + ConsJ[1,]*Lambda[1]; for (ii in 2:bas$Nc) GG = GG + ConsJ[ii,]*Lambda[ii]
    RHS <- c(GG,Cons_RHS)
    RNorm<-pracma::Norm(Cons_RHS)
    cat("|Cons-f|:",crayon::red(RNorm),": Cons-f:",crayon::red(Cons_RHS),"\n")
    #  if (RNorm>RNorm1) {cat("exit by Cons Resid Norm increase \n");break} else RNorm1=RNorm

    delta=  (pracma::pinv(Hp,tol = thresh) %*% RHS)[,1]

    FullNorm=pracma::Norm(delta)
    NewtonNorm=pracma::Norm(delta[-LambdaI])
    dLambda=delta[LambdaI]

    delta[-LambdaI] <- filter_delta(delta[-LambdaI])

    A1<-c(c(A),c(Lambda)) - del * c(delta) # full update
    A[]<-A1[-LambdaI]   # now store the updated primal
    Lambda=A1[LambdaI]  #  and dual
    cat("|delta_Newton|:",NewtonNorm,":  |delta|:",crayon::red(FullNorm),"\n")
    C0 <- C # save coordinates for possible rejection / break
    C <- updateX(A,grd,bas)
    E_SCM(A,grd,bas,C)->h2
    SEN(A,grd,bas,Ref1,h2)->S
    E_SEN(A,grd,bas,S,Ref1)->w
    E <- h2$Wb + w; cat ("CNM:",iter,"E",E/M.Es,"Wb",h2$Wb/M.Es,"Ws",w/M.Es,"C0",M.C0,"C",h2$Curv,"NEWTON STEP",del,"\n",sep=":")
    E_CNM[iter]<-E
    if(plt){ two_draw3d(A,M,title = paste("CNM",iter,round(E/M.Es,3)))
#      X2Obj(grd$Obj,C$X)->O
#      rgl::set3d(M.scr2);rgl::clear3d();
#      imag.obj.colorbar(O,f=S$beta,clr = FALSE,par=FALSE,specular="black"); rgl::title3d(paste(iter,"beta"))

#      rgl::set3d(M.scr1);rgl::clear3d();
#      imag.obj.colorbar(O,f=S$alpha,clr=FALSE,par=FALSE,specular="black");rgl::title3d(paste(iter,"alpha"))
    }
    attr(A,"Lambda")<-Lambda
    attr(A,"C0")<-M.C0
    attr(A,"Target")<-bas$Target
    attr(A,"Eigs")<-Eig # only in CNM Eigenvalues of H' are part of A attributes
    attr(A,"|delta_Newton|")<-NewtonNorm
    attr(A,"FullNorm")<-FullNorm
    attr(A,"iter")<-iter
    attr(A,"CNM_Cons_Resid")<-RNorm
    attr(A,"E")<-E
    attr(A,"V")<-h2$Volume
    attr(A,"A")<-h2$Area
    attr(A,"C")<-h2$Curv
    attr(A,"method")<-"CNM"

    cat("dLambda:",dLambda,"\n")
    cat("Lambda:",Lambda," :mu: ",M.mu," :del:",del,"\n")
    FullNorm0<-FullNorm
    if(FullNorm<Mtol_Newton) {cat(crayon::green("exit by Mtol_Newton\n"));MEx=TRUE;break;} else MEx=FALSE
    if(E>E0*2.5) {cat(crayon::red("reject A and exit (by >1.5 energy increase)\n"));Lambda=Lambda0;A=A0;C=C0;break;} # indicate divergence by rapid energy increase
    E0 <- E
    if(iter==1){cat("Full Constrained Newton Step ");tictoc::toc()}
    if (iter%%LAfreq==0) LA[[length(LA)+1]] <- A
  } # end of constraint Newton iteration loop
  M$A <- A;
  M$Lambda <- Lambda
  M$SEN <- S;
  M$E <- E
  M$C <- h2$Curv
  M$Eig <- Eig
  if (is.null(M$E_CNM)) M$E_CNM<-E_CNM else M$E_CNM<-c(M$E_CNM,E_CNM)
  M$CNMiter <- M$CNMiter+iter
  M$LA <- LA # return iterated solutions list
  M$last_App_called <- "CNM"
  M$history <- append(M$history,cl)
  t1 <- proc.time()
  M$proc_time <- M$proc_time + t1 - t0

  return(M)
}# end of CNM
