# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2026 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, Stephan (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# Frickenhaus S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#

<<<<<<< Updated upstream
#' MMC
#' @description
#'  Run Metropolis Montecarlo MMC on an existing membrane object under penalty constraints for volume and area.
#'  Constraints are set in M$bas$Target for penalization with strength M.rho.
#' @param M The input membrane with initial data and reference
#' @param nsteps number of MMC steps to be run
#' @param plt (=FALSE) for plotting
#' @param prn (=TRUE) for extensive printing
#' @param pltfreq (=5) for plotting shape every plotfreq number of accepted steps
=======
#'
#' run MMC on an existing membrane object under penalty constraints for volume and area.
#' Every 75 rejected steps, the interval size of random numbers added on coefficients is decreased by 5%.
#' If the interval size is below 1e-6, it is re-initialized to 0.01.
#' You can interrupt MMC without losing results by placing a file named "STOP_MMC.txt" in your working directory.
#'
#' @param M The input membrane with initial data and reference
#' @param nsteps number of MMC steps to be run
#' @param plt (=FALSE) for plotting
#' @param pertA (=pertA_unif) select move method (pertA_Unif, pertA_complex)
#' pertA_complex works with fitting some spatial normals, uses a global nearest-neighbor array in M.nn
#' defined in angle space. M.nn is computed on demand in MMC.
#' @param plfreq (=5) for plotting shape every plotfreq number of accepted steps
>>>>>>> Stashed changes
#' @param LAfreq (=200) for storing shape coefficients every LAfreq accepted steps
#' @param sd (=0.004) standard deviation parameter to control steps
#' @param kT (=0.00411) Boltzmann-constant * temperature energy scale
#' @param kTfac (=1.0) cooling factor (0<kTfac<1) or heating factor (kTfac>1)
<<<<<<< Updated upstream
#' @param kTfreq (=100) cooling/heating frequency
#' @param pertA (=pertA_Unif) perturbation scheme
#' @param record_dA (=FALSE) record perturbations
#' @param timing (=FALSE)
#' @param C0 (=M$C0) to verify wanted C0 against MemRBC_env$M.C0
#' @param ... further plot parameters
=======
#' @param kTfreq (=100) number of accepted steps between two scalings of kT by kTfac
>>>>>>> Stashed changes
#' @return membrane object with updates from MMC with data:
#' @return LA: list of recorded coefficients A
#' @return A: final coefficients
#' @return kT:  value of kT on exit, important for restarting
#' @return MMCacceptanceRate: acceptance rate of the call
#' @return MMCiter: total MMC steps, incl. MMC from previous calls
#' @return history: history of App-calls that created the result
#' @examplesIf exists("L_Ylm")
#' data(M4,package = "MemRBC"); M<-M4
#' plot(M)
#' #  annealing simulation (decrease kT by kTfac every Ktfreq accepted steps)
#' M1 <- MMC(M,nsteps=100000, kT=0.00411, kTfac=0.99, kTfreq=100, C0=-2)
#' plot(M1)
#' M1
#' 
#' @export
<<<<<<< Updated upstream
MMC<-function (M, nsteps = 1000, plt = TRUE, pltfreq = 10, prn = TRUE,
               LAfreq = 200, sd = 0.004, kT = 0.00411, kTfac = 1, kTfreq = 100,
               pertA = pertA_Unif, record_dA = FALSE, timing = FALSE, C0 = M$C0,
               ...)
{
  if (C0 != MemRBC_env$M.C0)
    stop("you should set global assign('M.C0',C0_value,envir = MemRBC::MemRBC_env) correctly by hand!")
  M$Params[["M.C0"]] <- C0
  t0 = proc.time()
  if (is.null(M$proc_time))
    M$proc_time <- 0
  if (is.null(MemRBC_env$M.Rcpp))
    stop("Cannot process - probably load_param_MemRBC has not been called.")
  MemRBC_env$M.Rcpp <- TRUE
  cl = match.call()
  if (identical(pertA, pertA_complex) & !exists("M.nn",envir = .GlobalEnv))
    M.nn <<- NNuv(M$grd$UV, 13)
  FM <- E_SCM(M$A, M$grd, M$bas, updateX(M$A, M$grd, M$bas))
  Cnt = rep(0L, nsteps)
  Ar = Cv = En = rep(0, nsteps)
  rec = 1
  bas = M$bas
  grd = M$grd
  Ref = M$Ref
  W = W0 = 1e+08
  a = r = aa = rr = 0
  A = M$A
  if (is.null(M$LA))
    LA = list(M$A)
  else LA = M$LA
  if (record_dA)
    Record = matrix(0, nsteps, length(c(A)) + 2)
  if (record_dA)
    attr(Record, "Remark") <- "A as vectors + last values are Energy and C0"
  if (is.null(M$MMCiter))
    M$MMCiter = 0
  if (is.null(M$kT))
    M$kT = kT
  else if (kT == 0)
    kT = M$kT
  run_id = rlang::hash(M)
  for (i in 1:nsteps) {
    if (i%%100 == 0)
      tictoc::tic()
    if (timing)
      if (i%%100 == 99) {
        cat("acceptance: ", round(aa/(aa + rr), 3), "\n 100 steps ")
        tictoc::toc()
      }
    A1 = pertA(A, bas, sd, FM$n)
    FM <- E_FullModel_Penalty(A1, grd, bas, Ref)
    W = FM$E
    if (i == 1) {
      oldA = FM$Area
      oldC = FM$Curv
      oldE = W
    }
    attr(A1, "run_id") = run_id
    attr(A1, "E") = FM$E
    attr(A1, "kT") = kT
    attr(A1, "C") = FM$Curv
    attr(A1, "V") = FM$Volume
    attr(A1, "A") = FM$Area
    attr(A1, "Target") = bas$Target
    attr(A1, "C0") = MemRBC_env$M.C0
    attr(A1, "sd") = sd
    attr(A1, "M.rho") <- MemRBC_env$M.rho
    attr(A1, "method") = "MMC"
    if (record_dA)
      Record[i, ] = c(c(A1), W, MemRBC_env$M.C0)
    if (min(1, exp(-(W - W0)/kT)) > runif(1)) {
      if (prn)
        cat("a :EAVC:", W/MemRBC_env$M.Es, FM$Area, FM$Volume, FM$Curv,
            ":C0:", MemRBC_env$M.C0, ":kT:", kT, "iter", i, "\n")
      W0 = W
      A = A1
      Ar[rec] = oldA
      Cv[rec] = oldC
      En[rec] = oldE
      Cnt[rec] = r + 1
      oldA = FM$Area
      oldC = FM$Curv
      oldE = W
      rec = rec + 1
      a = a + 1
      r = 0
      aa = aa + 1
      if (plt)
        if (aa%%pltfreq == 0) {
          M$A = A
          rgl::clear3d()
          plot(M, col = "white", ...)
          rgl::title3d(paste("MMC", round(W/MemRBC_env$M.Es, 5),
                             M$MMCiter + i, round(kT, 7), MemRBC_env$M.C0, round(FM$Curv,
                                                                      4)))
        }
      if (aa%%250 == 0)
        save(A, file = paste("A_L", bas$L_max, "_C0_",
                             MemRBC_env$M.C0, ".rdat", sep = ""))
      if (aa%%LAfreq == 0) {
        attr(A, "method") = "MMC"
        LA[[length(LA) + 1]] <- A
      }
      if (aa%%kTfreq == 0)
        kT <- kT * kTfac
    }
    else {
      if (prn)
        cat("r ")
      r = r + 1
      rr = rr + 1
      a = 0
    }
    if (r > 75) {
      sd = sd * 0.95
      r = 0
      a = 0
      cat("\nSD:", sd, "\n")
    }
    if (sd < 1e-06) {
      sd = 0.01
      r = 0
      a = 0
      cat("\nSD-reset:", sd, "\n")
    }
    if (file.exists("STOP_MMC.txt")) {
      message("break by STOP_MMC")
      file.remove("STOP_MMC.txt")
      break
    }
  }
  M$kT = kT
  M$C0_MMC = MemRBC_env$M.C0
  M$E = FM$E
  M$A = A
  M$LA = LA
  df = data.frame(Area = Ar[1:(rec - 1)], Curv = Cv[1:(rec -
                                                         1)], Energy = En[1:(rec - 1)], Cnt = Cnt[1:(rec - 1)],
                  Id = rep(run_id, (rec - 1)))
  if (!is.null(M$Sample))
    M$Sample = rbind(M$Sample, df)
  else M$Sample = df
  M$MMCiter = M$MMCiter + i
  M$MMCacceptanceRate = aa/(aa + rr)
  cat(crayon::green("Acceptance rate"), aa/(aa + rr), "\n")
  M$last_App_called = "MMC"
  M$history = append(M$history, list(cl))
  t1 = proc.time()
  M$proc_time <- M$proc_time + t1 - t0
  M$time_per_step_MMC = (t1 - t0)/i
  if (record_dA)
    M$Record = rbind(M$Record, Record[1:i, ])
=======
MMC<-function(M, nsteps=1000, plt=TRUE, pltfreq=10, LAfreq=200, sd=0.004, kT=4.11e-3,
              nm = FALSE, kTfac=1.0, kTfreq=100, filter=ID,
              pertA=pertA_Unif, record_dA=FALSE,timing=FALSE,C0=M$C0)
{ # set kTfac~0.99 < 1 for cooling
  # kTfreq in terms of accepted steps aa
  # LAfreq in terms of accepted steps aa
  # kT at room temperature is 4.11e-21 J = 4.11 E-3 atto J
  if (C0!=M.C0) stop("you should set M.C0 correctly by hand!")
  t0=proc.time()
  if(is.null(M$proc_time)) M$proc_time<-0
  if(!exists("M.Rcpp")) stop("Cannot process - probably load_param_MemRBC has not been called.")
  M.Rcpp<<-TRUE; # for faster E_SCM
  cl=match.call()
  if (identical(pertA,pertA_complex) & !exists("M.nn") ) M.nn<<- NNuv(M$grd$UV,13)
  FM<-E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas)) # for normals FM$n
  Cnt=rep(0L,nsteps)
  Ar=Cv=En=rep(0.0,nsteps); rec=1# records for sampling Area,Curv for Free Energy Perturbation
  bas=M$bas; grd=M$grd; Ref=M$Ref
  W=W0=1e8
  a=r=aa=rr=0 # acept/reject counters
  A=M$A # start coefficients
  if (is.null(M$LA)) LA=list(M$A) else LA=M$LA # MMC records
  if (record_dA ) Record=matrix(0,nsteps,length(c(A))+2) # MMC records mor ML
  if (record_dA) attr(Record,"Remark")<-"A as vectors + last values are Energy and C0"

  if (is.null(M$MMCiter)) M$MMCiter=0 # MMC counter
  if (is.null(M$kT)) M$kT=kT else if(kT==0) kT=M$kT # MMC restart with kT
  run_id=rlang::hash(M)
  # idea: for a replica exchange one could save exchange candidates A regularly
  for (i in 1:nsteps){
    if (i %% 100==0) tictoc::tic()
    if (timing) if (i %% 100==99) {cat("\n 100 steps ");tictoc::toc()}
    A1 = filter(pertA(A,bas,sd,FM$n),bas)
    if (i>1 & nm) { # surface normal moves filter; only possible when E$n is known
      dA = A1-A; make_delta_normal_to_surface(dA,grd,bas,E$n)->dA
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
    attr(A1,"M.rho")<-M.rho
    attr(A1,"method")="MMC"
    if(record_dA) Record[i,] = c(c(A1),W,M.C0)
  #  cat(W,W0,"ENER, logP=",log(min(0,(-(W-W0)/kT))),"\n")
    if (min(1,exp(-(W-W0)/kT)) > runif(1)) # acceptance rule
    {  cat("a :EAVC:", W/M.Es,  FM$Area, FM$Volume, FM$Curv,":C0:", M.C0,":kT:", kT,"\n");
      cat("i:",i,":");
      W0=W; A=A1;
      Ar[rec]=oldA
      Cv[rec]=oldC
      En[rec]=oldE
      Cnt[rec] = r + 1 # counts also the acceptance
      oldA=FM$Area
      oldC=FM$Curv
      oldE=W
      rec = rec+1
      a=a+1;r=0;aa=aa+1;
      if (plt) if (aa%%pltfreq==0) {M$A=A;rgl::clear3d();plot(M,col="white");
       rgl::title3d(paste("MMC",round(W/M.Es,5),M$MMCiter+i,round(kT,7),M.C0,round(FM$Curv,4)))}
      if (aa%%250==0) save(A,file=paste("A_L",bas$L_max,"_C0_",M.C0,".rdat",sep=""))
      if (aa%%LAfreq==0) {  attr(A,"method")="MMC";LA[[length(LA)+1]] <- A;} # record coeffs
      if (aa%%kTfreq==0) kT <- kT*kTfac # cooling, if 0<kTfac<1
    } else {cat("r ");r=r+1;rr=rr+1;a=0} # rejection
    if (r>75) {sd=sd*0.95;r=0;a=0; cat("\nSD:",sd,"\n")}
    if (sd<1e-6) {sd=0.01;r=0;a=0; cat("\nSD-reset:",sd,"\n")}
    if (file.exists("STOP_MMC.txt")) {message("break by STOP_MMC");file.remove("STOP_MMC.txt");break}
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
  M$time_per_step_MMC=(t1-t0)/i
  if (record_dA) M$Record=rbind(M$Record,Record[1:i,])
>>>>>>> Stashed changes
  return(M)
}


<<<<<<< Updated upstream
# A must not be matrix
#' perturb coeffs
#' @param A,bas coeffs and basis
#' @param sd standard dev. of gaussian
#' @param flt (=FALSE) to filter for only 3 entries in l=1 
#' @param n not used
#' @param ... not used
#' @export
pertA_Gauss<-function (A, bas, sd,  n, ...)
{
  N = bas$Ai_max
  sd1 = sd/sqrt(bas$G.tk)
  dA <- rnorm(3 * N, sd = rep(sd1, 3))
  if (flt) {
    dA[c(1 + N, 1 + 2 * N)] = 0
    dA[c(2, 2 + N)] = 0
    dA[c(3, 3 + 2 * N)] = 0
  }
  A[] = A[] + dA
  return(A)
}

#' perturb coeffs A by dA from uniform distribution, as needed in GAM and MMC
#' Perturbation dA is scaled down by `sqrt(bas$G.tk)`
#' @param A,bas coefficients and basis
#' @param sd width of perturbation
#' @param n normal vectors, e.g.from `E_SCM()`
#' @param ... not used
#' @return changed coeffs A+dA
#' @export
pertA_Unif<-function (A, bas, sd, n, ...)
{
  N = bas$Ai_max
  sd1 = sd/sqrt(bas$G.tk)
  dA = (runif(3 * N) - 0.5)/2 * rep(sd1, 3)
  A[] = A[] + dA
  return(A)
}


#' pertA_complex
#' @description
#'  Perturb with gaussian along spatial normals.
#'  BEFORE USE:
#' a) Normals of neighbours require to fill the global M.nn=Nuv(bas$uv)
#' b) the basis must contain a matrix IM for least squares inversion in FitFast();
#' Use bas<-MakeIM(bas,WX=sin(grd$U)) to create IM.
#' @param A,bas coefficients and basis
#' @param sd standard deviation, not needed here
#' @param n normal vectors
#' @param nn (=12) number of nearest neighbours to construct displacements from
#' @return modified coefficients A+dA
#' @export
pertA_complex<-function (A, bas, sd, n, nn = 12)
{ M.nn=get("M.nn",envir = .GlobalEnv)
  N = dim(bas$Ylm)[1] # ndof spatial
  s = sample(1:N, 3)
  dX = matrix(0, N, 3) # global effect to invert
  for (k in s) {
#    print(M.nn$nn.idx[k,])
    w = M.nn$nn.idx[k, 1:nn]
#    print(w)
    sc = max(M.nn$nn.dists[k, 1:nn]) * 4
    dX[w[1:nn], ] = dX[w[1:nn], ] + n[w[1:nn], ] * (runif(1) -
                                               0.5) * 0.3 * exp(-M.nn$nn.dists[k, 1:nn]^2/sc^2)
  }
  dA = FitFast(bas, dX)
  return(A + dA)
}

#' NNuv
#' @description
#'  compute nearest neighbors for pertA_complex
#'  store in global M.nn
#'  ATTENTION: returned distances $nn.dists are weighted sin(u), because in pole vicinity points get denser
#' @param uv (n x 2) matrix of (u,v) spherical angles
#' @param n (=13) number of nearest neighbours to compute
#' @return result from RANN::nn2: $nn.idx and $nn.dists
#' @export
NNuv <- function (uv, n = 13)
{
  nn <- RANN::nn2(uv, uv, n + 1)
#  nn$nn.dists = nn$nn.dists
  return(nn)
}

#' pertA_rndX
#' @description
#'  Perturb coeffs from spatial coordinates.
#'  BEFORE USE:
#' the basis mus contain a matrix IM for least squares inversion in FitFast();
#' Use bas<-MakeIM(bas,WX=sin(grd$U)) to create IM.
#' @param A,bas coeffs and basis
#' @param n not used
#' @param sd standard dev
#' @param sample_fraction (=0.65) fraction of points to be displaced
#' @return changed coeffs A+dA
#' @export
pertA_rndX<-function (A, bas, sd,n,sample_fraction=0.65)
{
  N = dim(bas$Ylm)[1]
  s = sample(1:N, floor(sample_fraction*N))
  dX = matrix(0, N, 3)
  dX[s, ] = rnorm(s*3,sd=sd)
  if (!is.null(bas$IM)) dA = FitFast(bas, dX) else try(stop("no IM in basis"))
  return(A + dA)
=======
#' @export
NNuv<-function(uv,n=8)
{ nn <- RANN::nn2(uv,uv,n+1) # include self
  nn$nn.dist=nn$nn.dist*sin(uv[,1])
  return(nn)
} # use: M.nn<-NNuv(M$grd$UV,12) to set global M.nn for pertA_complex

#' @export
pertA_complex<-function(A,bas,sd,n,nn=12)
{ N=dim(bas$Ylm)[1]
  s=sample(1:N,3)
  dX=matrix(0,N,3)
  for (k in s){
  w=M.nn$nn.idx[k,1:nn]
#  dX[w[1],]=n[w[1],]*sd[1] # central move
#  dX[w[2:5],]=n[w[2:5],]*sd[1]/2 # nn move
  sc=max(M.nn$nn.dists[k,1:nn])*4
  dX[w[1:nn],]=dX[w[1:nn],]+n[w[1:nn],] *(runif(1)-0.5)*0.3*exp(-M.nn$nn.dists[k,1:nn]^2/sc^2)
  }
  dA=FitFast(bas,dX) # intercept removal is in FitFast
  return(A+dA)
}


# A must not be matrix
#' @export
pertA_Gauss<-function(A,bas,sd,flt=FALSE,n){
  N=bas$Ai_max;
  sd1=sd/sqrt(bas$G.tk);
  dA<-rnorm(3*N,sd=rep(sd1,3))
  if (flt){
   dA[c(1+N,1+2*N)]=0
   dA[c(2,2+N)]=0
   dA[c(3,3+2*N)]=0}
  A[]=A[]+dA
  return(A)
}

#' @export
pertA_Unif<-function(A,bas,sd,flt=FALSE,n){
  N=bas$Ai_max;
  sd1=sd/sqrt(bas$G.tk);
  dA=(runif(3*N)-0.5)/2*rep(sd1,3) # > version 14: changed to sd1

  A[]=A[]+dA
  return(A)
}

#  erase unneccessary MMC data
#' @export
Streamline_MMC_data<-function(M){
M$LA=lapply(M$LA,function(x) {attr(x,"Fit spatial weights")=NULL;x})
M$LdA=lapply(M$LdA,function(x) {attr(x$dA,"Fit spatial weights")=NULL;x})
print(M)
return(M)
>>>>>>> Stashed changes
}
