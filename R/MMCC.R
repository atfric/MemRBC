# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2025 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, S. (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# Frickenhaus S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#

#' MMCC
#' @description
#'
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
#' if(exists("E_SCM_cxx")) { # catch problems in cran tests
#' data(M4)
#' M<-M4
#' plot(M)
#' #  annealing simulation (decrease kT by kTfac every Ktfreq accepted steps)
#' M <- MMCC(M, curv=Curv(M)+0.25, nsteps=10000, kT=0.00411, kTfac=0.99, kTfreq=100)
#' plot(M)
#' M }
#' @export
MMCC<-function (M, curv = Curv(M), nsteps = 1000, plt = FALSE, pltfreq = 5,
                LAfreq = 200, sd = 0.004, kT = 0.00411,  pertA = pertA_Unif,
                kTfac = 1, kTfreq = 100)
{
  t0 = proc.time()
  if (is.null(M$proc_time))
    M$proc_time <- 0
  if (is.null(MemRBC_env$M.Rcpp))
    stop("Cannot process - probably load_param_MemRBC has not been called.")
  MemRBC_env$M.Rcpp <- TRUE
  cl = match.call()
  run_id = rlang::hash(M)
  bas = M$bas
  bas$Nc = 3
  bas$Target = c(bas$Target, curv)
  bas$Cons = c("gradA", "gradV", "gradC")
  bas$QCons = c("Area", "Volume", "Curv")
  names(bas$Target) = bas$QCons
  Cnt = rep(0L, nsteps)
  Ar = Cv = En = rep(0, nsteps)
  rec = 1
  grd = M$grd
  Ref = M$Ref
  W = W0 = 1e+08
  a = r = aa = rr = 0
  A = M$A
  if (is.null(M$LA))
    LA = list(M$A)
  else LA = M$LA
  if (is.null(M$MMCCiter))
    M$MMCCiter = 0
  if (is.null(M$kT))
    M$kT = kT
  else if (kT == 0)
    kT = M$kT
  FM <- E_FullModel_Penalty(A, grd, bas, Ref)
  oldE = 1000
  for (i in 1:nsteps) {
    if (i%%100 == 0)
      tictoc::tic()
    if (i%%100 == 99) {
      cat("\n 100 steps ")
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
    attr(A1, "E") = FM$E
    attr(A1, "kT") = kT
    attr(A1, "C") = FM$Curv
    attr(A1, "V") = FM$Volume
    attr(A1, "A") = FM$Area
    attr(A1, "C0") = MemRBC_env$M.C0
    attr(A1, "sd") = sd
    attr(A1, "Target") = bas$Target
    attr(A1, "method") = "MMCC"
    attr(A, "M.rho") <- MemRBC_env$M.rho
    attr(A1, "run_id") = run_id
    if (min(1, exp(-(W - W0)/kT)) > runif(1)) {
      cat("a :EAVC:", W/MemRBC_env$M.Es, FM$Area, FM$Volume, FM$Curv,
          ":C0:", MemRBC_env$M.C0, ":kT:", kT, "\n")
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
          plot(M, col = "white")
          rgl::title3d(paste("MMCC", round(W/MemRBC_env$M.Es, 5),
                             M$MMCCiter + i, "kT", round(kT, 5), "C",
                             round(FM$Curv, 4)))
        }
      if (aa%%500 == 0)
        save(A, file = paste("A_MMCC_L", bas$L_max, "_C0_",
                             MemRBC_env$M.C0, ".rdat", sep = ""))
      if (aa%%LAfreq == 0)
        LA[[length(LA) + 1]] <- A
      if (aa%%kTfreq == 0)
        kT <- kT * kTfac
    }
    else {
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
      sd = 0.03
      r = 0
      a = 0
      cat("\nSD-reset:", sd, "\n")
    }
    if (file.exists("STOP_MMCC.txt")) {
      cat(crayon::red("exit by presence of file STOP_MMCC\n"))
      file.remove("STOP_MMCC.txt")
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
  M$MMCCiter = M$MMCCiter + i
  M$MMCCacceptanceRate = aa/(aa + rr)
  cat(crayon::green("MMCC Acceptance rate"), aa/(aa + rr),
      "\n")
  M$last_App_called = "MMCC"
  M$history = append(M$history, cl)
  t1 = proc.time()
  M$proc_time <- M$proc_time + t1 - t0
  return(M)
}

