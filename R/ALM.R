# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2025 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, S. (2025). MembraneRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#

#' Augmented Lagrangian Minimizer
#'
#' run ALM on an existing membrane object to adjust shape to target curvature, minimizing the residual norm RN.
#' @param M The input membrane with initial data and reference
#' @param curv target curvature
#' @param nsteps number of MMC steps to be run
#' @param plt (boolean) for plotting
#' @param method "SLSQP" for internal nloptr minimizer
#' @param kTfreq frequency for cooling factor applied on kT
#' @return membrane object with updates from MMC with data:
#' @return LA: list of recorded coefficients A
#' @return A: last coefficients
#' @return M.mu: last mu value
#' @return M.lam: last lambda vector
#' @return history: history of App-calls that created the result
#' @examples
#' M <- MakeStandardRBC(L=5)
#' plot(M)
#' M <- ALM(M,curv=86,20,plt=TRUE)
#' plot(M)
#' M
#' @export
ALM<-function(M, curv, nsteps=10,plt=TRUE,method="SLSQP",maxiter_solver=100,ncores=4,LAfreq=15)
{
 # energy and gradient for the solver:
 f_list_auglag=function(Ain,grd,bas,Ref1) #  needs bas, grd, Ref1 as global objects
 {
  A=matrix(Ain,ncol=3)
  C<-updateX(A,grd,bas)
  h2<-E_SCM(A,grd,bas,C)
  R<-ConsRHS(h2,bas)
  S<-SEN(A,grd,bas,Ref1,h2)
  SE<-E_SEN(A,grd,bas,S,Ref1)
  E_glob <<- h2$Wb + SE
  EAL_glob <<- h2$Wb + SE + M.muk/2 * sum(R^2) + sum(M.lam*R)
  g2<-Grad_SCM(h2,grd,bas,C)
  GS<-Grad_SEN(A,grd,bas,g2,S,Ref1)
  G <- c(g2$grad_SCM + GS$grad_SEN)
  cat(".")
  return( list("objective" = EAL_glob ,
               "gradient" = G + c( (M.lam + M.muk*R) %*% rbind(c(g2$gradA),
                                                               c(g2$gradV),
                                                               c(g2$gradC)))) )
}

# controls the solver and penalties:
AugLag_Step=function(A,tau=1.8,eta=1e-3,prec=1e-4,method="SLSQP",curv,maxiter_solver=300) # works with globals lam and muk for constraints and penalty
{ tictoc::tic()
  nlopt_opts <- list("algorithm"= paste("NLOPT_LD_",method,sep=""), # "NLOPT_LD_LBFGS", # "NLOPT_LD_SLSQP",
                     "xtol_rel"= prec,
                     "xtol_abs"= prec,
                     "maxeval" = maxiter_solver,

                     "print_level" = 0 )
  res_opt<-nloptr::nloptr( x0=c(A),
                   #      ub= rep(50/sqrt(bas$G.tk),3),
                   #      lb= -rep(50/sqrt(bas$G.tk),3),# more control needed to remain in reasonable A-ranges?
                   eval_f=f_list_auglag, # take objective and gradient from a list, saves computation of C and H2
                   opts= nlopt_opts , grd=M$grd, bas=bas,Ref1=M$Ref)
  cat("\n")
  A[]<-res_opt$solution
  C<-updateX(A,grd,bas)
  h2<-E_SCM(A,grd,bas,C)
  R<-ConsRHS(h2,bas)
  RN=pracma::Norm(R)
  cat("\nALM: Cons:",R,":RN:",RN,"\n")
  # here is the central augmented lagrangian update of lambda and mu
  if ( (RN)<eta ) M.lam<<-M.lam + M.muk*R else M.muk<<-M.muk*tau

  cat("\nChange Norm:",pracma::Norm(res_opt$x0-res_opt$solution),":|R|:",crayon::red(RN),":I:",crayon::red(res_opt$iterations),":S:", res_opt$status,":M:",res_opt$message,"\n")

  cat("AugLag_Step took ");tictoc::toc()
  return(list(A=A, R=R, RN=RN, iters=res_opt$iterations))
}

# actual ALM code starts here
if(!exists("M.Rcpp")) stop("Cannot process - probably load_param_MemRBC has not been called.")
M.Rcpp<<-TRUE # for faster SCM energy
M.Rcpp_ncores<<-ncores # and parallel SCM gradient
bas<-SetConstraints(M$bas,Cons=c("gradA","gradV","gradC"),
                    QCons = c("Area","Volume","Curv"),
                    Target=c(M$bas$Target[1:2],curv),
                    TNorm=c(M$bas$Target[1:2],curv)
)
t0=proc.time()
if(is.null(M$proc_time)) M$proc_time<-0
if (is.null(M$LA)) LA=list(M$A) else LA=M$LA # MMC records

cl <- match.call()
grd<-M$grd; A=M$A

C<-updateX(A,grd,bas)
h2<-E_SCM(A,grd,bas,C)
#if (plt) #two_screens3d()

S=SEN(A,grd,bas,M$Ref,h2)
if(plt) {rgl::clear3d();plot3b(C$X,grd)}# two_draw3d(A,M)

ll.a=0

if (is.null(M$ALM_RN)) RN0=1000 else  RN0<-M$ALM_RN; # allow first step
if (is.null(M$ALM_mu)) M.muk<<-50 #else  M.muk<<-M$ALM_mu;
if (is.null(M$ALM_lambda)) M.lam<<-c(1,1,1) #else  M.lam<<-M$ALM_lambda;

eta=0.6
prec=1e-4

iter_count=0

if(!is.null(M$M.muk)) M.muk<<-M$M.muk
if(!is.null(M$M.lam)) M.lam<<-M$M.lam
tictoc::tic()
for (iter in 1:nsteps) { # usually 10 cycles
  if (iter==1) tictoc::tic()
  if (prec>1e-8) prec=prec/2 # reduce by 2^10 in 10 steps
  AL = AugLag_Step(A, tau=1.9,eta=eta, prec=prec, method=method,maxiter_solver=maxiter_solver) # SLSQP
  iter_count=iter_count+AL$iters
  eta=eta*0.6
#  cat("H(A) ",rlang::hash(A),"H(AL$A)",rlang::hash(AL$A),"\n")
  A = LM2A(AL$A,bas) # update A from solver
  if(iter %% LAfreq==0) LA[[length(LA)+1]]<-A

  C<-updateX(A,grd,bas)
  h2<-E_SCM(A,grd,bas,C)
  R<-ConsRHS(h2,bas)
  S<-SEN(A,grd,bas,M$Ref,h2)
  SE<-E_SEN(A,grd,bas,S,M$Ref)

#  g2<-Grad_SCM(h2,grd,bas,C)
#  GS<-Grad_SEN(A,grd,bas,g2,S,M$Ref)

  E <- h2$Wb + SE
  cat("ALM:",iter,":E_AL:",EAL_glob/M.Es,":E:",E/M.Es,crayon::green(":mu:"),M.muk,crayon::green(":lam:"),M.lam,"\n")
  cat("ALM:",iter,":E:",E/M.Es,":Ct:",bas$Target[3],":C0:",M.C0,":C:",h2$Curv,"\n")

#  h2<-E_SCM(A,grd,bas,C) no update needed
#  S<-SEN(A,grd,bas,M$Ref,h2)

  if (plt) {rgl::clear3d();plot3b(C$X,grd);
    rgl::title3d(paste("C=",round(h2$Curv,3),"Ct=",
                       bas$Target["Curv"],"C0=",M.C0,"L=",bas$L_max))}# ,round(lambdaG[3],3)))
  #      attr(A,"Lambda")=lambdaG
  attr(A,"E")=E
  attr(A,"Target")=bas$Target
  attr(A,"C")=h2$Curv
  attr(A,"C0")=M.C0
  attr(A,"ALM_RN0")<-RN0
  attr(A,"iter")<-iter
  if (iter==1) {cat("ALM full cycle took ");tictoc::toc()}
  if (abs(AL$RN-RN0)<1e-19) { cat(crayon::yellow("Exit ALM by zero change of residual norm\n This reflects the experimental status of the current ALM implementation.\n"));break } else RN0=AL$RN

  #      attr(A,"lambdaRosen")<-lambdaG
  #      attr(A,"|G|_Rosen")<-GN

  LA[[length(LA)+1]]<-A;


}

cat("full ALM took "); tictoc::toc()
M$A=A
M$LA=LA;
if (is.null(M$ALMiter)) M$ALMiter=iter else M$ALMiter=M$ALMiter+iter
if (is.null(M$ALMiter_solve)) M$ALMiter_solve=iter_count else M$ALMiter_solve=M$ALMiter_solve+iter_count
M$ALM_lambda<-M.lam
M$ALM_mu<-M.muk
M$ALN_RN0=RN0
M$last_App_called="ALM"

M$history=append(M$history,list(cl))
t1=proc.time()
M$proc_time <- M$proc_time + t1-t0

return(M)
}# end of ALM

