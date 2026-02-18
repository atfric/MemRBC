# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2025 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, S. (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# Frickenhaus S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#

#' SDRC
#' @description
#' Steepest Descend with Rosen Constraint Solver
#'
#' Run SDRC on an existing membrane object to reduce energy by steepest descend.
#' Constraint violation is projected out by Rosens method (1961).
#' To stay near feasible points, constraint fix point iterations (Frickenhaus, Dissertation thesis, Humboldt-University Berlin, 1999) are performed every step.
#' @param M The input membrane with initial data and reference
#' @param nsteps number of SD-steps to be run (with varying number of constraint iteration steps CI in between)
#' @param del_min (=1e-7) step-size of SD
#' @param del_cons (=1e-3) step-size of fix point iteration
#' @param max_iter (=10) maximum number of constraint iterations
#' @param cons_tol (=1e-3) constraint iterations stop if max(abs(Cons_Values)/bas$Target)<Ctol, coded in ConsIter()
#' @param plt (boolean) for plotting
#' @param pltfreq number of iterations between plots
#' @param filter_strength weight of diagonal filter. i.e. 1+f/sqrt(M$bas$G.tk)
#' @param LAfreq storage frequency into list of coefficients LA
#' @return membrane object with updates from MMC with data:
#' @return LA: list of recorded coefficients A
#' @return A: last coefficients
#' @return C: last curvature
#' @return E: last total energy
#' @return SDRC_Sample: recorded quantities (E,A,V,C, CN, IC, ID), where CN ist constraint norm, IC is constraint iteration count, and ID is the run-id.
#' @return SDRCiter; total SDRC steps, incl. previous calls
#' @return history: history of App-calls that created the result
#' @examples
#' M <- MakeStandardRBC(L=5)
#' plot(M)
#' M.C0 <- 0
#' M <- SDRC(M, nsteps=10000, del=1e-6, LAfreq=100)
#' plot(M)
#' M
#' attributes(last(M$LA))
#' @export
SDRC <- function (M, nsteps = 100, del_min = 1e-07, del_cons = 0.001,
                  max_iter = 10, cons_tol = 0.001, plt = FALSE, pltfreq = 10,
                  LAfreq = 25, filter_strength = 0, prn_ci=FALSE)
{
  cl = match.call()
  run_id = rlang::hash(M)
  CIiter = 0
  t0 = proc.time()
  if (is.null(M$proc_time))
    M$proc_time <- 0
  filter = (1 + filter_strength * 1/sqrt(M$bas$G.tk))
  E_SD = C_SD = V_SD = A_SD = I_SD = CN = rep(0, nsteps)
  A = M$A
  grd = M$grd
  bas = M$bas
  Ref = M$Ref
  if (!is.null(M$LA))  LA = M$LA else LA = list(M$A)
  if (M.mu == 0 & M.Ka == 0)
    GS = list(grad_SEN = matrix(0, bas$Ai_max, 3))
  
  for (i in 1:nsteps) {
    C <- updateX(A, grd, bas)
    E <- E_SCM(A, grd, bas, C)
    if (M.mu != 0 | M.Ka != 0) {
      S <- SEN(A, grd, bas, Ref, E)
      ES <- E_SEN(A, grd, bas, S, Ref)
    }
    else ES = 0
    GE <- Grad_SCM(E, grd, bas, C)
    if (M.mu != 0 | M.Ka != 0)  GS <- Grad_SEN(A, grd, bas, GE, S, Ref)
      G <- RosenProjection((GE$grad_SCM + GS$grad_SEN) * filter,
                         GE, bas) # else: GS=0, set above
      
    A <- A - del_min * G$Gprime
    C <- updateX(A, grd, bas)
    E <- E_SCM(A, grd, bas, C)
    if (M.mu != 0 | M.Ka != 0) {
      S <- SEN(A, grd, bas, Ref, E)
      ES <- E_SEN(A, grd, bas, S, Ref)
    }
    else ES = 0
    GE <- Grad_SCM(E, grd, bas, C)
    CI <- ConsIter(A, grd, bas, C, GE, Ctol = cons_tol, del_cons = del_cons,
                   nsteps = max_iter,prn=prn_ci)
    A <- CI$A
    CIiter <- CIiter + CI$cons_iter
    Et <- E$Wb + ES
    E_SD[i] = Et
    C_SD[i] = E$Curv
    V_SD[i] = E$Volume
    A_SD[i] = E$Area
    I_SD[i] = CI$cons_iter
    CN[i] = pracma::Norm(CI$Cons_RHS)
    cat("SDRC:", i, "E", Et/M.Es, "C", E$Curv, "CN", CN[i],
        "del_min", del_min, "C0", M.C0, "F", filter_strength,
        "CI", CI$cons_iter, "\n")
    if (plt & (i%%pltfreq == 0)) {
      rgl::clear3d()
      plot3b(C$X, grd)
      rgl::title3d(paste("SDRC", i, "E", round(Et/M.Es,
                                               4), "C", round(E$Curv, 4)))
    }
    attr(A, "method") = "SDRC"
    attr(A, "E") <- Et
    attr(A, "C") <- E$Curv
    attr(A, "C0") <- M.C0
    attr(A, "Target") <- bas$Target
    attr(A, "run_id") <- run_id
    if (i%%LAfreq == 0)   {cat (crayon::red("REC LN",length(LA),"\n"));LA[[length(LA) + 1]] <- A }
    if (file.exists("STOP_SDRC.txt")) {
      file.remove("STOP_SDRC.txt")
      cat(crayon::red("STOP SDRC from extern\n"))
      break
    }
  }
  M$A = A
  E = E_FullModel_Penalty(A, grd, bas, Ref)
  M$E = E$E
  M$C = E$Curv
  if (is.null(M$SDRC_Sample))
    M$SDRC_Sample <- data.frame(E = E_SD[1:i], A = A_SD[1:i], V = V_SD[1:i],
                                C = C_SD[1:i], CN = CN[1:i], IC = I_SD[1:i], ID = run_id[1:i])
  else M$SDRC_Sample <- rbind(M$SDRC_Sample, data.frame(E = E_SD[1:i],
                                                        A = A_SD[1:i], V = V_SD[1:i], C = C_SD[1:i], CN = CN[1:i], IC = I_SD[1:i], ID = run_id[1:i]))
  if (is.null(M$SDRCiter))
    M$SDRCiter <- i
  else M$SDRCiter <- i + M$SDRCiter
  M$last_App_called = "SDRC"
  M$history = append(M$history, cl)
  M$LA = LA
  M$SDRC_filter_strength <- filter_strength
  M$ConsIter = CIiter
  t1 = proc.time()
  M$SDRC_Target = bas$Target
  M$proc_time <- M$proc_time + t1 - t0
  return(M)
}
# end of SDRC
