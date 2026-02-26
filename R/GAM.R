# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2025 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, S. (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# Frickenhaus S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#

#' Genetic Algorithm Minizer (worked with GA package before 4.5.1)
#'
#' run GAM on a precomputed membrane object
#' @param M The input membrane with initial data and reference
#' @param nsteps number of GA-cycles to be run
#' @param maxiter number of iterations in each GA cycle
#' @return membrane with updates from GAM
#' @return LA: list of recorded coefficients A
#' @return GA: result of last GA cycle, needed for restart
#' @return A: last coefficients
#' @return history: history of App-calls that created the result
#' @examplesIf exists("L_Ylm")
#' data(M4,package = "MemRBC")
#' M <- GAM(M4,nsteps=3)
#' plot(M)
#' M 
#' 
#' @export
GAM<-function(M,nsteps=1,maxiter=50)
{
  # perturb only for the initial population
  cl=match.call()
  pertA<-function(A,bas,sd){ N=dim(A)[1];
  sd1=sd/sqrt(bas$G.tk);
<<<<<<< Updated upstream
  A[]=A[]+matrix(rnorm(3*N,sd=rep(sd1,3)),ncol=3)
=======
  A[]=A[] + matrix(rnorm(3*N,sd=rep(sd1,3)),ncol=3)
>>>>>>> Stashed changes
  return(A)
  }
  pop.restart<-function(O)
  {
    return(GA@population)
  }
  mon<-function(O)
  {
    print(O@fitness)
    print(O@bestSol)
  }
  plotA3d<-function(A)
  { updateX(matrix(A,ncol=3),grd,bas)->C
    h2<-E_SCM(A,grd,bas,C)
    plot3a(C$X,grd)
    print(h2$Curv)
    rgl::title3d(h2$Curv)
  }
  # perturbed restart population from a single coefficient matrix A
  pop<-function(O)
  { P=matrix(0, O@popSize, bas$Ai_max*3)
  for (i in 1:O@popSize) P[i,]=pertA(A,sd=0.02,bas)
  return(P)
  }
  mem_fit<-function(Ain) # fitness function is negative energy
  {A=matrix(Ain,ncol=3)
  updateX(A,grd,bas)->C
  h2<-E_SCM(A,grd,bas,C)
  S=SEN(A,grd,bas,Ref1,h2)
  e<-E_SEN(A,grd,bas,S,Ref1)
  W=h2$H2 + e + MemRBC_env$M.rho*((h2$Volume-bas$Target["Volume"])^2+(h2$Area-bas$Target["Area"])^2)
  return(-W)
  }
  # actual start of GA code
  t0=proc.time()
  if(is.null(M$proc_time)) M$proc_time<-0

  bas=M$bas;grd=M$grd;A=M$A;Ref1=M$Ref
  if (is.null(M$LA)) LA=list(M$A) else LA=M$LA
  if (is.null(M$GAMiter)) M$GAMiter=0
  sd <- 4*rep(1/sqrt(bas$G.tk),3) # to control upper/lower relative to current A
  for (i in 1:nsteps){
    if (i==1) pop_select=pop else pop_select=pop.restart
    GA::ga(type = "real-valued",mem_fit ,
       lower=c(A)-sd,upper=c(A)+sd, # fixed boundaries, require multiple runs to adapt
       maxiter=maxiter,
       mutation = GA::gareal_nraMutation,
       keepBest = TRUE,
       population = pop_select, # with pop.restart GA@population is used; pop needs A to vary
       monitor=TRUE, # no curvature output "C:" when parallel
       optim=FALSE,
       parallel=1 # speedup>2 for 4 worker-cores
    )->GA
    # read out results from GA
    M$GAMiter <- M$GAMiter + GA@iter
    A<-M$A;  A[]<-matrix(last(GA@bestSol)[[1]],ncol=3)
    M$A <- matrix(A,ncol=3)
    attr(M$A,"method")="GAM"
    LA[[ length(LA)+1 ]] <- M$A
  }
  M$GA <- GA # save whole result of last cycle population
  M$E <- (-GA@fitnessValue)
  M$E_GAM <- (-GA@fitnessValue)
  M$last_App_called<-"GAM"
  M$LA<-LA
  M$history=append(M$history,list(cl))
  t1=proc.time()
  M$proc_time <- M$proc_time + t1-t0

  return(M)
}# end of GAM
