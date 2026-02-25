# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2025 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, S. (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# Frickenhaus S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#


#' Penalized Steepest Descend
#'
#' run PSD on an existing membrane object to reduce energy by gradient of energy.
#' A more efficient minimizer may be the PNEM.
#' @param M The input membrane with initial data and reference
#' @param nsteps (=100) number of PSDC steps to be performed
#' @param del (=1e-6) step-size
#' @param filter_strength (=0) to scale down higher-l gradient components by `Filter = (1 + filter_strength * 1/sqrt(M$bas$G.tk))`
#' @param plt (=FALSE) for plotting
#' @param pltfreq (=10) for plotting every pltfreq step
#' @param LAfreq (=100) storage frequency into list of coefficients $LA
#' @param ncores (=4) parallel threads in SCM energy and gradient
#' @return membrane object with updates from MMC with data:
#' @return LA: list of recorded coefficients A
#' @return A: last coefficients
#' @return C: last curvature
#' @return E: last total energy
#' @return E_PSD: all recorded total energies
#' @return PSDiter; total PSD steps, incl. previous calls
#' @return history: history of App-calls that created the result
#' @examplesIf exists("L_Ylm")
#' data("M4",package = "MemRBC")
#' M<-M4
#' plot(M,alpha=0.65)
#' M <- PSD(M,nsteps=100,del=1e-6,LAfreq=10)
#' plot(M,add=TRUE,col="red")
#' M
#' attributes(last(M$LA))
#' 
#' @export
PSD <- function (M, nsteps = 100, del = 1e-06, plt = FALSE, pltfreq = 10,
                 LAfreq = 100, ncores = 4, filter_strength = 0)
{
  cl = match.call()
  run_id = rlang::hash(M)
  t0 = proc.time()
  if (is.null(M$proc_time))
    M$proc_time <- 0
  Filter = (1 + filter_strength * 1/sqrt(M$bas$G.tk))
  E_PSD = C_PSD = rep(0, nsteps)
  A = M$A
  grd = M$grd
  bas = M$bas
  Ref = M$Ref
  if (is.null(M$LA))
    LA = M$LA
  else LA = list(M$A)
  for (i in 1:nsteps) {
    E <- E_FullModel_Penalty(A, grd, bas, Ref)
    NCons <- sum((unlist(E[bas$QCons[1:bas$Nc]]) - bas$Target[1:bas$Nc])^2)
    C <- updateX(A, grd, bas)
    if (!is.null(M$Ref)) {
      S <- SEN(A, grd, bas, Ref, E_SCM(A, grd, bas, C))
    }
    else Ref = NULL
    G <- Grad_FullModel_Penalty(A, grd, bas, Ref, S)
    G <- G / Filter
    any(is.infinite(G))
    any(is.na(G))
    A = A - del * matrix(G, ncol = 3)
    any(is.na(A))
    any(is.infinite(A))

    cat("PSD:", i, ":E:", E$E/MemRBC_env$M.Es, ":C:", E$Curv, ":del:",
        del, ":C0", MemRBC_env$M.C0, ":F:", filter_strength, "|Cons|",
        NCons, "\n")
    if (plt & (i%%pltfreq == 0)) {
      rgl::clear3d()
      plot3b(C$X, grd)
      rgl::title3d(paste("PSD", i, "E", round(E$E/MemRBC_env$M.Es,
                                              4), "C", round(E$Curv, 4)))
    }
    attr(A, "method") = "PSD"
    attr(A, "E") <- E$E
    attr(A, "C") <- E$Curv
    attr(A, "C0") <- MemRBC_env$M.C0
    attr(A, "Target") <- bas$Target
    attr(A, "M.rho") <- MemRBC_env$M.rho
    attr(A, "run_id") <- run_id
    if (i%%LAfreq == 0)
      LA[[length(LA) + 1]] <- A
    E_PSD[i] = E$E
    C_PSD[i] = E$Curv
    if (file.exists("STOP_PSD.txt")) {
      message("break by STOP_PSD")
      file.remove("STOP_PSD.txt")
      break
    }
  } # step loop
  M$A = A
  E = E_FullModel_Penalty(A, grd, bas, Ref)
  M$E = E$E
  M$CurvPSD = E$Curv
  if (is.null(M$E_PSD))
    M$E_PSD = E_PSD
  else M$E_PSD = c(M$E_PSD, E_PSD)
  if (is.null(M$C_PSD))
    M$C_PSD = C_PSD
  else M$C_PSD = c(M$C_PSD, C_PSD)
  if (is.null(M$PSDiter))
    M$PSDiter = i
  else M$PSDiter = M$PSDiter + i
  M$last_App_called = "PSD"
  cl = paste(cl, "# steps done: ", i)
  M$history = append(M$history, list(cl))
  M$LA = LA
  M$PSD_filter_strength <- filter_strength
  t1 = proc.time()
  M$proc_time <- M$proc_time + t1 - t0
  return(M)
}
