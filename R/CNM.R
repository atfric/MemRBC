# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2025 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, S. (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# Frickenhaus S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#

# REWRITTEN for solving without bordered Hessian and corrected sign definition of multipliers Lambda

#' Constraint Newton Minimizer
#'
#' run CNM on an existing membrane object to minimize energy under exact constraints on area and volume
#' You may include curvature constraint as well.
#' Temperature is controlled by parameter kT.
#' @param M The input membrane with initial data and reference
#' @param nsteps (=5) number of Newton steps to be performed
#' @param del (=0.3) update factor delte (<1)
#' @param del_control (=FALSE) TRUE for use of statoíc step length del
#' @param plt (=TRUE) for plotting two figures with alpha/beta stress-shear parameters
#' @param Mtol_Newton (=1e-4) stopping criterion for KKT norm and constraints norm
#' @param LAfreq (=1) storage frequency into list of coefficients LA
#' @param cluster (=FALSE) for parallel Hessian using a cluster with ncores processes
#' @param ncores (=4) for parallel Hessian on ncores CPU-cores
#' @param keepEig2A (=FALSE) to keep eigensystem as attribute "EigsH" in $A (and thus LA[[]] -> high memory demands)
#' @param pinv (=TRUE) use pseudo-inverse in solver 
#' @param diag.reg (=0) diagonal regularization of to be inverted matrix
#' @return membrane object with updated data from CNM:
#' @return LA : list of recorded coefficients A
#' @return A : last coefficients
#' @return history : history of App-calls that created the result
#' @return E_CNM : vector of nsteps energy values
#' @return CNMiter : number of iterations, incl. previous calls
#' @return proc_time : aggregated processing time of all Apps called before with this object
#' @return Eig.H : eigensystem of last energy+constraints Hessian
#' @return Grad : final energy gradient
#' @return Jacobian : constraint Jacobian
#' @return Hessian : full energy Hessian matrix
#' @return ConsHessians : Lagrangian Hessians from active constraints
#' @examplesIf exists("L_Ylm")
#' data(M4)
#' MemRBC_env$M.C0 <- -4
#' plot(M4)
#' M4<-SDRC(M4,10,del_min=3e-6)
#' M4 <- CNM(M4,nsteps=20,cluster=TRUE,ncores=4,del_control=TRUE)
#' plot(M4)
#' M
#' 
#' @export
CNM <- function (M, nsteps = 5, del = 0.3, del_control = FALSE, diag.reg = 0,
                 pinv = TRUE, LAfreq = 1, Mtol_Newton = 1e-04, ncores = 4,
                 keepEig2A = FALSE, plt = TRUE, cluster = FALSE)
{ cl=NULL # no cluster yet
  t0 = proc.time()
  if (del_control)
    alpha = NA
  else del = NA
  if (is.null(M$Ref))
    stop("CNM needs SEN Reference; use MakeRef(M)->M;  If done so,  to switch off SEN, set MemRBC_env$M.mu=MemRBC_env$M.Ka=0.")
  if (is.null(M$proc_time))
    M$proc_time <- 0
  if (is.null(MemRBC_env$M.Rcpp))
    stop("Cannot process - probably load_param_MemRBC has not been called.")
  MemRBC_env$M.Rcpp <- TRUE
  MemRBC_env$M.Rcpp_ncores <- ncores
  E0 = 1000
  if (is.null(M$LA))
    LA = list(M$A)
  else LA = M$LA
  if (is.null(M$CNMiter))
    M$CNMiter = 0
  CNM_data = matrix(0, ncol = 10, nrow = nsteps)
  A = M$A
  grd = M$grd
  bas = M$bas
  Ref1 = M$Ref
  if (is.null(M$Lambda))
    Lambda = rep(1, bas$Nc)
  else Lambda = M$Lambda
  names(Lambda) = bas$Cons[1:bas$Nc]
  for (iter in (1:nsteps)) {
    if (iter == 1)   tictoc::tic()
    A0 = A
    Lambda0 = Lambda
    if (iter == 1)  {tictoc::tic(); cl=NULL}
    if (cluster) {
      H <- FullModelHessian_Par(A, grd, bas, Ref1, del = 1e-06,
                                Mem_mc.cores = ncores, 
                                timing = iter == 1, # startup if cl is NULL
                                stopdown = iter == nsteps, cl=cl)
      cl<-H$cl
    } else H <- FullModelHessian(A, grd, bas, Ref1, del = 1e-06)
    
#    else { # this is incomplete Hessian, not better performing than cluster version
#      if (use_serial_cpp) {
#        warning("no Hessian of constraints in serial C++ Hessian")
#        H <- Hessian_FullModel(A, grd, bas, Ref1, 1e-06,
#                               ncores = ncores)
#      }
    diag(H$H) <- diag(H$H) * (1 + diag.reg)
    if (iter == 1) {
      cat("parallel Hessian ")
      tictoc::toc()
    }
    ConsJ <- t(sapply(H$g2[bas$Cons[1:bas$Nc]], c))
    if (iter == 1) {
      Lambda <- (-pracma::pinv(ConsJ %*% t(ConsJ)) %*%
                   ConsJ %*% H$G)[, 1]
      cat(crayon::green("estimate Lambda0:", paste(Lambda,
                                                   collapse = " "), "\n"))
    }
    names(Lambda) <- names(Lambda0) <- bas$Cons[1:bas$Nc]
    Hfull = H$H
    for (i in bas$Cons[1:bas$Nc]) Hfull <- Hfull + Lambda[i] *
      H$L[[i]]
    Cons <- ConsRHS(H$h2, bas)
    Eig.H <- eigen(Hfull)
    if (all(Eig.H$values > 0))
      cat(crayon::white(crayon::bgGreen("all eigenvalues positive \n"))) else
        cat(sum(Eig.H$values > 0)," of ",length(Eig.H$values)," positive\n")
    if (diag.reg == 0)
      iH = mypinv(Hfull)
    else iH = pracma::inv(Hfull)
    if (attr(iH, "removed") > 0)
      cat(crayon::red("removed ", attr(iH, "removed"),
                      "\n"))
    v <- iH %*% H$G
    Lambda_QP <- (pracma::pinv(ConsJ %*% iH %*% t(ConsJ)) %*%
                    (Cons - ConsJ %*% v))[, 1]
    names(Lambda) <- names(Lambda_QP) <- bas$Cons[1:bas$Nc]
    GG <- c(H$G) + t(ConsJ) %*% Lambda_QP
    KKTnorm <- pracma::Norm(GG)
    cat(crayon::blue("KKT Norm:", KKTnorm))
    RNorm <- pracma::Norm(Cons)
    cat(" |Cons|:", crayon::red(RNorm), "\n")
    delta <- (-iH %*% GG)[, 1]
    if (del_control) {
      A <- A + del * delta
    }
    C0 <- C
    C <- updateX(A, grd, bas)
    h2 <- E_SCM(A, grd, bas, C)
    S <- SEN(A, grd, bas, Ref1, h2)
    w <- E_SEN(A, grd, bas, S, Ref1)
    E <- h2$Wb + w
    if (del_control)
      cat("CNM", iter, crayon::green("E"), crayon::green(E/MemRBC_env$M.Es),
          "Wb", h2$Wb/MemRBC_env$M.Es, "Ws", w/MemRBC_env$M.Es, "C0", MemRBC_env$M.C0, "C",
          h2$Curv, "del", del, "\n", sep = ":")
    if (!del_control) {
      mu = max(abs(Lambda)) + 1
      alpha = 1
      current_merit = E + mu * sum(abs(Cons))
      while (alpha > 1e-04) {
        A1 = A + alpha * delta
        C <- updateX(A1, grd, bas)
        h2 <- E_SCM(A1, grd, bas, C)
        S <- SEN(A1, grd, bas, Ref1, h2)
        w <- E_SEN(A1, grd, bas, S, Ref1)
        E <- h2$Wb + w
        Cons_next = ConsRHS(h2, bas)
        next_merit = E + mu * sum(abs(Cons_next))
        if (next_merit < current_merit)
          break
        alpha = alpha * 0.5
        cat("lowering alpha:", alpha, "\r")
      }
      A = A + alpha * delta
      Lambda = Lambda + alpha * (Lambda_QP - Lambda)
      cat("Lambda (updated):", Lambda, "Lambda_QP:", Lambda_QP,
          "\n")
      cat("\nCNM", iter, crayon::green("E"), crayon::green(E/MemRBC_env$M.Es),
          "Wb", h2$Wb/MemRBC_env$M.Es, "Ws", w/MemRBC_env$M.Es, "C0", MemRBC_env$M.C0, "C",
          h2$Curv, "alpha", alpha, "\n", sep = ":")
    }
    if (del_control)
      Lambda = Lambda + del * (Lambda_QP - Lambda)
    CNM_data[iter, ] <- c(E, h2$Wb, w, KKTnorm, RNorm, h2$Curv,
                          MemRBC_env$M.C0, iter, del, alpha)
    if (plt)
      two_draw3d(A, M, title = paste("CNM", iter, round(E/MemRBC_env$M.Es,
                                                        3)))
    attr(A, "Lambda") <- Lambda
    attr(A, "C0") <- MemRBC_env$M.C0
    attr(A, "Target") <- bas$Target
    if (keepEig2A) {
      attr(A, "EigsH") <- Eig.H
    }
    attr(A, "KKTNorm") <- KKTnorm
    attr(A, "iter") <- iter
    attr(A, "CNM_Cons_Resid") <- RNorm
    attr(A, "E") <- E
    attr(A, "V") <- h2$Volume
    attr(A, "A") <- h2$Area
    attr(A, "C") <- h2$Curv
    attr(A, "method") <- "CNM"
    MEx = FALSE
    if (max(KKTnorm, RNorm) < Mtol_Newton) {
      cat(crayon::green("exit by Mtol_Newton on max(KKTnorm, ConsNorm)\n"))
      MEx = TRUE
      break
    }
    if (file.exists("STOP_CNM.txt")) {
      cat(crayon::red("exit by presence of file STOP_CNM\n"))
      file.remove("STOP_CNM.txt")
      MEx = TRUE
      break
    }
    if (E > E0 * 1.5) {
      cat(crayon::red("reject A and exit (by > E*0.5 energy increase)\n"))
      Lambda = Lambda0
      A = A0
      C = C0
      break
    }
    E0 <- E
    if (iter == 1) {
      cat("Full Constrained Newton Step ")
      tictoc::toc()
    }
    if (iter%%LAfreq == 0)
      LA[[length(LA) + 1]] <- A
  }
  M$A <- A
  M$Lambda <- Lambda
  M$SEN <- S
  M$E <- E
  M$C <- h2$Curv
  M$Eig.H <- Eig.H
  M$Grad <- H$G
  M$Jacobian <- ConsJ
  M$Hessian <- H$H
  M$ConsHessians <- H$L
  M$CNMiter <- M$CNMiter + iter
  M$LA <- LA
  if (is.null(M$CNM_data))
    M$CNM_data = CNM_data
  else M$CNM_data = rbind(M$CNM_data, CNM_data)
  colnames(M$CNM_data) = c("E", "Wb", "Ws", "KKTnorm", "RNorm",
                           "C", "M.C0", "iter", "del", "alpha")
  M$last_App_called <- "CNM"
  M$history <- append(M$history, match.call() )
  t1 <- proc.time()
  M$proc_time <- M$proc_time + t1 - t0
  return(M)
}

# end of CNM


mypinv<-function (A, tol = .Machine$double.eps^(2/3))
{
  stopifnot(is.numeric(A) || is.complex(A), is.matrix(A))
  s <- svd(A)
  rem=0
  if (is.complex(A))
    s$u <- Conj(s$u)
  p <- (s$d > max(tol * s$d[1], 0))
  if (all(p)) {
    mp <- s$v %*% (1/s$d * t(s$u))
  }
  else if (any(p)) {
    mp <- s$v[, p, drop = FALSE] %*% (1/s$d[p] * t(s$u[,
                                                       p, drop = FALSE]))
    rem=sum(!p)
  }
  else {
    mp <- matrix(0, nrow = ncol(A), ncol = nrow(A))
  }
  attr(mp,"removed")=rem
  attr(mp,"diag_4")=s$d[order(abs(s$d))][1:4]
  attr(mp,"thresh")=max(tol * s$d[1], 0)
  return(mp)
}


