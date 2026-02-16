# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2025 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, Stephan (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# Frickenhaus S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#

#' MMC
#' @description
#'  Run Metropolis Montecarlo MMC on an existing membrane object under penalty constraints for volume and area.
#'  Constraints are set in M$bas$Target for penalization with strength M.rho.
#' @param M The input membrane with initial data and reference
#' @param nsteps number of MMC steps to be run
#' @param plt (=FALSE) for plotting
#' @param prn (=TRUE) for extensive printing
#' @param plfreq (=5) for plotting shape every plotfreq number of accepted steps
#' @param LAfreq (=200) for storing shape coefficients every LAfreq accepted steps
#' @param sd (=0.004) standard deviation parameter to control steps
#' @param kT (=0.00411) Boltzmann-constant * temperature energy scale
#' @param nm (=FALSE) normal motion filter
#' @param kTfac (=1.0) cooling factor (0<kTfac<1) or heating factor (kTfac>1)
#' @param filter (=ID) no other filter implemented, do it yourself
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
#' M1 <- MMC(M,nsteps=100000, kT=0.00411, kTfac=0.99, kTfreq=100, C0=-2)
#' plot(M1)
#' M1
#' @export
MMC<-function (M, nsteps = 1000, plt = TRUE, pltfreq = 10, prn = TRUE,
               LAfreq = 200, sd = 0.004, kT = 0.00411, kTfac = 1, kTfreq = 100,
               pertA = pertA_Unif, record_dA = FALSE, timing = FALSE, C0 = M$C0,
               ...)
{
  if (C0 != M.C0)
    stop("you should set global M.C0 correctly by hand!")
  t0 = proc.time()
  if (is.null(M$proc_time))
    M$proc_time <- 0
  if (!exists("M.Rcpp"))
    stop("Cannot process - probably load_param_MemRBC has not been called.")
  M.Rcpp <<- TRUE
  cl = match.call()
  if (identical(pertA, pertA_complex) & !exists("M.nn"))
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
    attr(A1, "C0") = M.C0
    attr(A1, "sd") = sd
    attr(A1, "M.rho") <- M.rho
    attr(A1, "method") = "MMC"
    if (record_dA)
      Record[i, ] = c(c(A1), W, M.C0)
    if (min(1, exp(-(W - W0)/kT)) > runif(1)) {
      if (prn)
        cat("a :EAVC:", W/M.Es, FM$Area, FM$Volume, FM$Curv,
            ":C0:", M.C0, ":kT:", kT, "iter", i, "\n")
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
          rgl::title3d(paste("MMC", round(W/M.Es, 5),
                             M$MMCiter + i, round(kT, 7), M.C0, round(FM$Curv,
                                                                      4)))
        }
      if (aa%%250 == 0)
        save(A, file = paste("A_L", bas$L_max, "_C0_",
                             M.C0, ".rdat", sep = ""))
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
  M$C0_MMC = M.C0
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
  return(M)
}


# A must not be matrix
#' @export
pertA_Gauss<-function (A, bas, sd, flt = FALSE, n)
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


#' @export
pertA_Unif<-function (A, bas, sd, flt = FALSE, n)
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
#' a) Normals of neighbours require to fill global M.nn=Nuv(bas$uv)
#' b) the basis mus contain a matrix IM for least squares inversion in FitFast();
#' Use bas<-MakeIM(bas,WX=sin(grd$U)) to create IM.
#' @export
pertA_complex<-function (A, bas, sd, n, nn = 12)
{
  N = dim(bas$Ylm)[1]
  s = sample(1:N, 3)
  dX = matrix(0, N, 3)
  for (k in s) {
    w = M.nn$nn.idx[k, 1:nn]
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
#' @export
NNuv <- function (uv, n = 13)
{
  nn <- RANN::nn2(uv, uv, n + 1)
  nn$nn.dist = nn$nn.dist * sin(uv[, 1])
  return(nn)
}

#' pertA_rndX
#' @description
#'  Perturb coeffs from spatial coordinates.
#'  BEFORE USE:
#' the basis mus contain a matrix IM for least squares inversion in FitFast();
#' Use bas<-MakeIM(bas,WX=sin(grd$U)) to create IM.
#' @export
pertA_rndX<-function (A, bas, sd,n,sample_fraction=0.65)
{
  N = dim(bas$Ylm)[1]
  s = sample(1:N, floor(sample_fraction*N))
  dX = matrix(0, N, 3)
  dX[s, ] = rnorm(s*3,sd=sd)
  if (!is.null(bas$IM)) dA = FitFast(bas, dX) else try(stop("no IM in basis"))
  return(A + dA)
}
