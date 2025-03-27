# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2025 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, S. (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# Frickenhaus S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#

#'
#' run MMC on an existing membrane object under penalty constraints for volume and area.
#' Every 75 rejected steps, the interval size of random numbers added on coefficients is decreased by 5%.
#' If the interval size is below 1e-6, it is re-initialized to 0.01.
#' @param M The input membrane with initial data and reference
#' @param nsteps number of MMC steps to be run
#' @param plt (=FALSE) for plotting
#' @param plfreq (=5) for plotting shape every plotfreq number of accepted steps
#' @param LAfreq (=200) for storing shape coefficients every LAfreq accepted steps
#' @param sd (=0.004) standard deviation parameter to control steps
#' @param kT (=0.00411) Boltzmann-constant * temperature energy scale
#' @param nm (=FALSE) normal motion filter
#' @param kTfac (=1.0) cooling factor (0<kTfac<1) or heating factor (kTfac>1)
#' @param kTfreq (=100) number of accepted steps between two scalings of kT by kTfac
#' @return membrane object with updates from MMC with data:
#' @return LA: list of recorded coefficients A
#' @return A: final coefficients
#' @return kT:  value of kT on exit, important for restarting
#' @return MMCacceptanceRate: acceptance rate of the call
#' @return MMCiter: total MMC steps, incl. MMC from previous calls
#' @return history: history of App-calls that created the result
#' @examples
#' M <- MakeStandardRBC(L=5)
#' plot(M)
#' #  annealing simulation (decrease kT by kTfac every Ktfreq accepted steps)
#' M <- MMC(M,nsteps=100000, kT=0.00411, kTfac=0.99, kTfreq=100)
#' plot(M)
#' M
#' @export
MMC<-function(M, nsteps=1000, plt=FALSE, pltfreq=5, LAfreq=200, sd=0.004, kT=4.11e-3,
              nm = FALSE, kTfac=1.0, kTfreq=100, filter=ID)
{ # set kTfac~0.99 < 1 for cooling
  # kTfreq in terms of accepted steps aa
  # LAfreq in terms of accepted steps aa
  # kT at room temperature is 4.11e-21 J = 4.11 E-3 atto J
  t0=proc.time()
  if(is.null(M$proc_time)) M$proc_time<-0

  if(!exists("M.Rcpp")) stop("Cannot process - probably load_param_MemRBC has not been called.")
  M.Rcpp<<-TRUE; # for faster E_SCM
  cl=match.call()
  Cnt=rep(0L,nsteps)
  Ar=Cv=En=rep(0.0,nsteps); rec=1# records for sampling Area,Curv for Free Energy Perturbation
  bas=M$bas; grd=M$grd; Ref=M$Ref
  W=W0=1e8
  a=r=aa=rr=0 # acept/reject counters
  A=M$A # start coefficients
  if (is.null(M$LA)) LA=list(M$A) else LA=M$LA # MMC records
  if (is.null(M$MMCiter)) M$MMCiter=0 # MMC counter
  if (is.null(M$kT)) M$kT=kT else if(kT==0) kT=M$kT # MMC restart with kT
  run_id=rlang::hash(M)
  # idea: for a replica exchange one could save exchange candidates A regularly
  for (i in 1:nsteps){
    if (i %% 100==0) tictoc::tic()
    if (i %% 100==99) {cat("\n 100 steps ");tictoc::toc()}
    #  if (i==250000) kT=0.1
    A1=filter(pertA_Unif(A,bas,sd),bas)
    if (i>1 & nm) { # surface normal moves filter; only possible when E$n is known
      dA = A1-A; make_delta_normal_to_surface(dA,grd,bas,E$n)
      A1 = A + dA
    }

    FM<-E_FullModel_Penalty_AV(A1,grd,bas,Ref)
    W=FM$E
    if(i==1) {oldA=FM$Area;oldC=FM$Curv;oldE=W}
    attr(A1,"run_id")=run_id
    attr(A1,"E")=FM$E
    attr(A1,"kT")=kT
    attr(A1,"C")=FM$Curv
    attr(A1,"V")=FM$Volume
    attr(A1,"A")=FM$Area
    attr(A1,"Target")=bas$Target
    attr(A1,"C0")=M.C0
    attr(A1,"sd")=sd
    attr(A,"M.rho")<-M.rho
    attr(A1,"method")="MMC"

    if (min(1,exp(-(W-W0)/kT)) > runif(1)) # acceptance rule
    {  cat("a :EAVC:", W/M.Es,  FM$Area, FM$Volume, FM$Curv,":C0:", M.C0,":kT:", kT,"\n");
      W0=W;A=A1;
      Ar[rec]=oldA
      Cv[rec]=oldC
      En[rec]=oldE
      Cnt[rec] = r + 1 # counts also the acceptance
      oldA=FM$Area
      oldC=FM$Curv
      oldE=W
      rec= rec+1
      a=a+1;r=0;aa=aa+1;
      if (plt) if (aa%%pltfreq==0) {M$A=A;rgl::clear3d();plot(M,col="white");
      rgl::title3d(paste("MMC",round(W/M.Es,5),M$MMCiter+i,round(kT,7),M.C0,round(FM$Curv,4)))}
      if (aa%%250==0) save(A,file=paste("A_L",bas$L_max,"_C0_",M.C0,".rdat",sep=""))
      if (aa%%LAfreq==0) {  attr(A,"method")="MMC";LA[[length(LA)+1]] <- A;} # record coeffs
      if (aa%%kTfreq==0) kT <- kT*kTfac # cooling, if 0<kTfac<1
    } else {cat("r ");r=r+1;rr=rr+1;a=0} # rejection
    if (r>75) {sd=sd*0.95;r=0;a=0; cat("\nSD:",sd,"\n")}
    if (sd<1e-6) {sd=0.01;r=0;a=0; cat("\nSD-reset:",sd,"\n")}

  }
  M$kT=kT
  M$C0_MMC=M.C0
  M$E=FM$E
  M$A=A
  M$LA=LA
  df=data.frame(Area=Ar[1:(rec-1)],Curv=Cv[1:(rec-1)],Energy=En[1:(rec-1)],Cnt=Cnt[1:(rec-1)],Id=rep(run_id,(rec-1)))
  if(!is.null(M$Sample)) M$Sample=rbind(M$Sample,df) else M$Sample=df
  M$MMCiter=M$MMCiter+i
  M$MMCacceptanceRate=aa/(aa+rr)
  cat(crayon::green("Acceptance rate"),aa/(aa+rr),"\n")
  M$last_App_called="MMC"
  M$history=append(M$history,list(cl))
  t1=proc.time()
  M$proc_time <- M$proc_time + t1-t0

  return(M)
} # end of MMC


