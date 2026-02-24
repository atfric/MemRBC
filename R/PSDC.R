# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2025 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, S. (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# Frickenhaus S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#


#' Penalized Steepest Descend Curvature
#'
#' run PSDC on an existing membrane object to reduce energy by gradient steps.
#' Here, curvature is an additional constraint (by penalty).
#' @param M The input membrane with initial data and reference
#' @param curv (=Curv(M))target curvature
#' @param nsteps (=100) number of PSDC steps to be performed
#' @param del (=1e-6) step-size
#' @param plt (=FALSE) for plotting
#' @param pltfreq (=10) for plotting every pltfreq step
#' @param LAfreq (=100) storage frequency into list of coefficients $LA
#' @param ncores (=4) parallel threads in SCM energy and gradient
#' @return membrane object with updates from MMC with data:
#' @return LA: list of recorded coefficients A
#' @return A: last coefficients, with further data in it as attributes (see example)
#' @return C: last curvature
#' @return E: last total energy
#' @return E_PSD: all recorded total energies, including PSDC steps
#' @return PSDiter; total PSD steps, incl. previous calls
#' @return history: history of App-calls that created the result
#' @examplesIf exists("L_Ylm")
#' data("M4",package = "MemRBC")
#' M <- PSDC(M4,nsteps=10000,del=1e-6,)
#' plot(M)
#' M
#' attributes(last(M$LA))
#' attr(M$A,"method")
#' 
#' @export
PSDC <- function (M, curv = Curv(M), nsteps = 100, del = 1e-06, plt = FALSE,
                  pltfreq = 10, LAfreq = 100, ncores = 5)
{
  t0 = proc.time()
  if (is.null(M$proc_time))
    M$proc_time <- 0
  cl = match.call()
  run_id = rlang::hash(M)
  MemRBC_env$M.Rcpp <- TRUE
  MemRBC_env$M.Rcpp_ncores <- ncores
  E_PSD = C_PSD = rep(0, nsteps)
  A = M$A
  grd = M$grd
  bas = M$bas
  Ref = M$Ref
  bas$Nc = 3
  bas$Target = c(bas$Target, curv)
  bas$Cons = c("gradA", "gradV", "gradC")
  bas$QCons = c("Area", "Volume", "Curv")
  names(bas$Target) = bas$QCons
  if (is.null(M$LA))
    LA = M$LA
  else LA = list(M$A)
  tictoc::tic()
  for (i in 1:nsteps) {
    E <- E_FullModel_Penalty(A, grd, bas, Ref)
    C <- updateX(A, grd, bas)
    S <- SEN(A, grd, bas, Ref, E_SCM(A, grd, bas, C))
    G <- Grad_FullModel_Penalty(A, grd, bas, Ref, S)
    A = A - del * matrix(G, ncol = 3)
    cat("PSDC:", i, ":E:", E$E/MemRBC_env$M.Es, ":C:", E$Curv, ":del:",
        del, "\n")
    if (plt & (i%%pltfreq == 0)) {
      rgl::clear3d()
      plot3q(C$X, grd)
      rgl::title3d(paste("PSDC", i, "E", round(E$E/MemRBC_env$M.Es,
                                               4), "C", round(E$Curv, 4)))
    }
    if (i == 100) {
      cat("100 steps took ")
      tictoc::toc()
    }
    attr(A, "E") <- E$E
    attr(A, "C") <- E$Curv
    attr(A, "C0") <- MemRBC_env$M.C0
    attr(A, "Target") <- bas$Target
    attr(A, "run_id") <- run_id
    attr(A, "M.rho") <- MemRBC_env$M.rho
    attr(A, "method") = "PSDC"
    if (i%%LAfreq == 0)
      LA[[length(LA) + 1]] <- A
    E_PSD[i] = E$E
    C_PSD[i] = E$Curv
    if (file.exists("STOP_PSDC.txt")) {
      cat(crayon::red("exit by presence of file STOP_PSDC\n"))
      file.remove("STOP_PSDC.txt")
      break
    }
  }
  M$A = A
  E = E_FullModel_Penalty(A, grd, bas, Ref)
  M$E = E$E
  M$C = E$Curv
  if (is.null(M$E_PSD))
    M$E_PSD = E_PSD
  else M$E_PSD = c(M$E_PSD, E_PSD)
  if (is.null(M$C_PSD))
    M$C_PSD = C_PSD
  else M$C_PSD = c(M$C_PSD, C_PSD)
  if (is.null(M$PSDiter))
    M$PSDiter = i
  else M$PSDiter = M$PSDiter + i
  M$last_App_called = "PSDC"
  M$history = append(M$history, cl)
  M$LA = LA
  t1 = proc.time()
  M$proc_time <- M$proc_time + t1 - t0
  return(M)
}
# end of PSD
