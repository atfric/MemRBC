# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2025 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, S. (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# Frickenhaus S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#


#' Penalized Newton Equation of Motion Viscosity Matrix
#'
#' run PNEMVM on an existing membrane object to simulate dynamics with spectral natrix viscosity
#' @param M The input membrane with initial data and reference
#' @param nsteps number of time steps to be run
#' @param dt time step
#' @param rho area mass density (scalar)
#' @param plt (boolean) for plotting
#' @param zero_Av set velocity zero for pure minimization
#' @param new_Av set velocity according to sd0 random numbers
#' @param sd0 standard deviation of Gaussian velocity initialization
#' @param LAfreq storage frequency into list of coefficients LA
#' @param visc_fac viscosity scaling factor for damping oscillations to zero
#' @param visc_matrix viscosity matrix of size $bas$Ai_max x $bas$Ai_max
#' @param diag_visc (boolean) viscosity tensor diagonal, and not ~ mass-matrix
#' @return membrane object with updates from MMC with data:
#' @return LA: list of recorded coefficients A
#' @return A: last coefficients
#' @return C: last curvature
#' @return E: last total energy
#' @return Av: last coefficient velocity, good for restart
#' @return E_kin_PNEM: kinetic energy of all previous PNEMs
#' @return E_total_PNEM: total energy of all previous PNEMs
#' @return PNEMiter; total PNEM step counter, incl. previous calls
#' @return history: history of App-calls that created the result
#' @examples
#' data(L5_stomatocyte_equilib)
#' M<-L5_stomatocyte_equilib
#' plot(M)
#' M.Rcpp=TRUE
#' M.Rcpp_ncores=4
#' M <- PNEMVM(M, nsteps=10000, del=1e-6)
#' plot(M)
#' M
#' # further data stored as attributes to coefficients A in LA
#' attributes(last(M$LA))
#' @export
PNEMVM <- function (M, nsteps = 100, dt = 5e-04, LAfreq = 50, plt = TRUE,
                    filter_delta = ID, rho = 1, sd0 = 0.001, pltfreq = 25, ncores = 4,
                    new_Av = FALSE, visc_fac = 1, visc_matrix = diag(sqrt(M$bas$G.tk[1]/M$bas$G.tk)),
                    zero_Av = FALSE, mass_update_freq = 100)
{
  t0 = proc.time()
  if (is.null(M$proc_time))
    M$proc_time <- 0
  run_id = rlang::hash(M)
  if (!exists("M.Rcpp"))
    stop("Cannot process - probably load_param_MemRBC has not been called.")
  M.Rcpp <<- TRUE
  M.Rcpp_ncores <<- ncores
  cl = match.call()
  E0 = 1000
  E_total = rep(0, nsteps)
  E_kin = rep(0, nsteps)
  type_PNEM = rep("PNEMVM", nsteps)
  C_PNEM = rep(0, nsteps)
  A = M$A
  grd = M$grd
  bas = M$bas
  Ref = M$Ref
  if (!is.null(M$LA))
    LA = M$LA
  else LA = list(A)
  C = updateX(A, grd, bas)
  if (is.null(M$LA))
    LA = list(M$A)
  else LA = M$LA
  if (is.null(M$PNEMiter))
    M$PNEMiter = 0
  if (is.null(M$Av) | new_Av) {
    Av = A
    Av[] = 0
    Av = pertA_Gauss(Av, bas, sd = sd0)
    print("Init Av")
  }
  else Av = M$Av
  if (zero_Av)
    Av[] = 0
  if (plt) {
    two_screens3d()
  }
  for (iter in (1:nsteps)) {
    A0 = A
    if (iter == 1)
      tictoc::tic()
    if (iter == 1) {
      E <- E_FullModel_Penalty(A, grd, bas, Ref)
      G <- Grad_FullModel_Penalty(A, grd, bas, Ref, E$S)
      tictoc::tic()
      M$dA = E$dA
      if (!is.null(M$mass)) {
        cat("re-use mass from given membrane\n")
        mass = M$mass
      }
      else mass = massmatrix(M, rho)
      inv_mass = pracma::inv(mass)
      cat("PNEM MassMatrix ")
      tictoc::toc()
      Aa <- inv_mass %*% (-G - visc_fac * visc_matrix %*%
                            Av)
    }
    else Aa <- Aa1
    Av <- Av + Aa * 0.5 * dt
    A <- A + Av * dt
    E <- E_FullModel_Penalty(A, grd, bas, Ref)
    G <- Grad_FullModel_Penalty(A, grd, bas, Ref, E$S)
    if (iter%%mass_update_freq == 0) {
      M$dA <- E$dA
      mass <- massmatrix(M, rho)
      inv_mass <- pracma::inv(mass)
    }
    Aa1 <- inv_mass %*% (-G - visc_fac * visc_matrix %*%
                           Av)
    Av <- Av + Aa1 * 0.5 * dt
    Xdot <- synthX(bas$Ylm, Av)
    Ekin_X <- 0
    for (j in 1:3) Ekin_X <- Ekin_X + 0.5 * rho * .IntegS(Xdot[,
                                                              j]^2 * E$dA, grd)
    Ekin <- 0
    for (j in 1:3) Ekin <- Ekin + (0.5 * Av[, j] %*% mass %*%
                                     Av[, j])[1, 1]
    Etot <- Ekin + E$E
    E_total[iter] <- Etot
    E_kin[iter] <- Ekin
    C_PNEM[iter] <- E$Curv
    cat("|dAv|:", pracma::Norm(Aa * dt), ":  |dA|:", pracma::Norm(Av *
                                                                    dt), ":  Ekin:", Ekin/M.Es, ": Ekin_X:", Ekin_X/M.Es,
        ":Ekin/Ekin_X:", Ekin/Ekin_X, ":Etot:", Etot/M.Es,
        "\n")
    cat("PNEM", iter, "E", E$E/M.Es, "Wb", E$Wb/M.Es, "Ws",
        E$Ws/M.Es, "C0", M.C0, "C", E$Curv, "A", E$Area,
        "V", E$Volume, "dt", dt, "v_f", visc_fac, "\n", sep = ":")
    if (plt & iter%%pltfreq == 0) {
      two_draw3d(A, M, title = paste("PNEM", iter))
    }
    attr(A, "method") <- "PNEM"
    attr(A, "C0") <- M.C0
    attr(A, "Target") <- bas$Target
    attr(A, "iter") <- iter
    attr(A, "E") <- E$E
    attr(A, "V") <- E$Volume
    attr(A, "A") <- E$Area
    attr(A, "C") <- E$Curv
    attr(A, "rho") <- rho
    attr(A, "visc_fac") <- visc_fac
    attr(A, "Etot") <- Etot
    attr(A, "rho_mass") <- rho
    attr(A, "M.rho") <- M.rho
    attr(A, "run_id") <- run_id
    E0 <- E$E
    if (iter == 1) {
      cat("Full Penalized Newton Equation Step ")
      tictoc::toc()
    }
    if (iter%%LAfreq == 0)
      LA[[length(LA) + 1]] <- A
  }
  M$A = A
  M$SEN = E$S
  if (!is.null(M$E_total_PNEM))
    M$E_total_PNEM = c(M$E_total_PNEM, E_total)
  else M$E_total_PNEM = E_total
  if (!is.null(M$E_kin_PNEM))
    M$E_kin_PNEM = c(M$E_kin_PNEM, E_kin)
  else M$E_kin_PNEM = E_kin
  if (!is.null(M$type_PNEM))
    M$type_PNEM = c(M$type_PNEM, type_PNEM)
  else M$type_PNEM = type_PNEM
  if (!is.null(M$C_PNEM))
    M$C_PNEM = c(M$C_PNEM, C_PNEM)
  else M$C_PNEM = C_PNEM
  M$mass = mass
  M$visc_matrix = visc_matrix
  M$E = E$E
  M$C = E$Curv
  M$PNEMiter = M$PNEMiter + iter
  M$LA = LA
  M$Av = Av
  M$rho_PNEM = rho
  M$last_App_called = "PNEMVM"
  M$history = append(M$history, list(cl))
  t1 = proc.time()
  M$proc_time <- M$proc_time + t1 - t0
  return(M)
}
# end of PNEMV

