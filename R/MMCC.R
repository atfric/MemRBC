# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2025 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, S. (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# Frickenhaus S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#


#' Metropolis Montecarlo with Curvature Constraint
#'
#' run MMC on an existing membrane object under penalty constraints for volume, area and curvature
#' @param M The input membrane with initial data and reference
#' @param curv the constraint curvature (default is current curvature of M)
#' @param nsteps number of MMC steps to be run
#' @param plt (boolean) for plotting
#' @param plofreq for plotting shape every plotfreq number of accepted steps
#' @param LAfreq for storing shape coefficients every LAfreq accepted steps
#' @param sd standard deviation parameter to control steps
#' @param kT Boltzmann-constant * temperature energy scale
#' @param nm (boolen) normal motion filter
#' @param kTfac cooling factor
#' @param kTfreq frequency for cooling factor applied on kT in units of accepted steps
#' @return membrane object with updates from MMC with data:
#' @return LA: list of recorded coefficients A
#' @return A: final coefficients
#' @return kT:  value of kT on exit, important for restarting
#' @return MMCCacceptanceRate: acceptance rate of the call
#' @return MMCCiter: total MMCC steps, incl. MMCC from previous calls
#' @return history: history of App-calls that created the result
#' @examples
#' M <- MakeStandardRBC(L=5)
#' plot(M)
#' #  annealing simulation (decrease kT by kTfac every Ktfreq accepted steps)
#' M <- MMCC(M, curv=Curv(M)+0.25, nsteps=100000, kT=0.00411, kTfac=0.99, kTfreq=100)
#' plot(M)
#' M
#' @export
MMCC<-function(M, curv=Curv(M), nsteps=1000, plt=FALSE, pltfreq=5, LAfreq=200, sd=0.004, kT=4.11e-3,
              nm = FALSE, kTfac=1.0, kTfreq=100)
{ # set kTfac~0.99 < 1 for cooling
  # kTfreq in terms of accepted steps aa
  # LAfreq in terms of accepted steps aa
  # kT at room temperature is 4.11e-21 J = 4.11 E-3 atto J
  t0=proc.time()
  if(is.null(M$proc_time)) M$proc_time<-0

  if(!exists("M.Rcpp")) stop("Cannot process - probably load_param_MemRBC has not been called.")
  M.Rcpp<<-TRUE; # for faster E_SCM
  cl=match.call()
  run_id=rlang::hash(M)
  bas=M$bas

  bas$Nc=3;bas$Target=c(bas$Target,curv);
  bas$Cons=c("gradA","gradV","gradC");
  bas$TNorm=c(bas$TNorm,curv);
  bas$Qcons=c("Area","Volume","Curv");
  names(bas$Target)=bas$Qcons

  Cnt=rep(0L,nsteps)
  Ar=Cv=En=rep(0.0,nsteps); rec=1# records for sampling Area,Curv for Free Energy Perturbation
  grd=M$grd; Ref=M$Ref
  W=W0=1e8
  a=r=aa=rr=0 # acept/reject counters
  A=M$A # start coefficients
  if (is.null(M$LA)) LA=list(M$A) else LA=M$LA # MMC records
  if (is.null(M$MMCCiter)) M$MMCCiter=0 # MMC counter
  if (is.null(M$kT)) M$kT=kT else if(kT==0) kT=M$kT # MMC restart with kT

  # idea: for a replica exchange one could save exchange candidates A regularly
  for (i in 1:nsteps){
    if (i %% 100==0) tictoc::tic()
    if (i %% 100==99) {cat("\n 100 steps ");tictoc::toc()}
    #  if (i==250000) kT=0.1
    A1=pertA_Unif(A,bas,sd)
    if (i>1 & nm) { # surface normal moves filter; only possible when E$n is known
      dA = A1-A; make_delta_normal_to_surface(dA,grd,bas,E$n)
      A1 = A + dA
    }

    FM<-E_FullModel_Penalty_AVC(A1,grd,bas,Ref)
    W=FM$E
    if(i==1) {oldA=FM$Area;oldC=FM$Curv;oldE=W}
    attr(A1,"E")=FM$E
    attr(A1,"kT")=kT
    attr(A1,"C")=FM$Curv
    attr(A1,"V")=FM$Volume
    attr(A1,"A")=FM$Area
    attr(A1,"C0")=M.C0
    attr(A1,"sd")=sd
    attr(A1,"Target")=bas$Target
    attr(A1,"method")="MMCC"
    attr(A,"M.rho")<-M.rho
    attr(A1,"run_id")=run_id


    if (min(1,exp(-(W-W0)/kT)) > runif(1)) # acceptance rule
    {  cat("a :EAVC:", W/M.Es,  FM$Area, FM$Volume, FM$Curv,":C0:", M.C0,":kT:", kT,"\n");
      W0=W;A=A1;
      Ar[rec]=oldA
      Cv[rec]=oldC
      En[rec]=oldE
      Cnt[rec]=r+1
      oldA=FM$Area
      oldC=FM$Curv
      oldE=W
      rec=rec+1
      a=a+1;r=0;aa=aa+1;
      if (plt) if (aa%%pltfreq==0) {M$A=A;rgl::clear3d();plot(M,col="white");
      rgl::title3d(paste("MMCC",round(W/M.Es,5),M$MMCCiter+i,"kT",round(kT,5),"C",round(FM$Curv,4)))}
      if (aa%%500==0) save(A,file=paste("A_MMCC_L",bas$L_max,"_C0_",M.C0,".rdat",sep=""))
      if (aa%%LAfreq==0) LA[[length(LA)+1]] <- A; # record coeffs
      if (aa%%kTfreq==0) kT <- kT*kTfac # cooling, if 0<kTfac<1
    } else {cat("r ");r=r+1;rr=rr+1;a=0} # rejection
    if (r>75) {sd=sd*0.95;r=0;a=0; cat("\nSD:",sd,"\n")}
    if (sd<1e-6) {sd=0.03;r=0;a=0; cat("\nSD-reset:",sd,"\n")}
    #  if (sd<0.007) sd=0.015
  }
  M$kT=kT
  M$C0_MMC=M.C0
  M$E=FM$E
  M$A=A
  M$LA=LA
  df=data.frame(Area=Ar[1:(rec-1)],Curv=Cv[1:(rec-1)],Energy=En[1:(rec-1)],Cnt=Cnt[1:(rec-1)],Id=rep(run_id,(rec-1)))
  if(!is.null(M$Sample)) M$Sample=rbind(M$Sample,df) else M$Sample=df

  M$MMCCiter=M$MMCCiter+i
  M$MMCCacceptanceRate=aa/(aa+rr)
  cat(crayon::green("MMCC Acceptance rate"),aa/(aa+rr),"\n")
  M$last_App_called="MMCC"
  M$history=append(M$history,list(cl))
  t1=proc.time()
  M$proc_time <- M$proc_time + t1-t0

  return(M)
} # end of MMC

