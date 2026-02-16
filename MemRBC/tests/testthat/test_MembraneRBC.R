data("ss42denovo")
ss42denovo->M
if (FALSE) { # synth_update_inplace not testable when in package; you could define it in your R code
  C=updateX(D5$A,D5$grd,D5$bas)
  synth_update_inplace(C,D5$bas, 4,3)
  C1=updateX(D5$A,D5$grd,D5$bas)
  C1=synth_update(C1,D5$bas, 4,3)

  severe(max(C1$X-C$X),0,"update vs update_inplace",1e-10)
  C=updateX(M$A,M$grd,M$bas)
#  microbenchmark::microbenchmark(synth_update_inplace(C,M$bas,sample(1:M$bas$Ai_max,1),sample(1:3,1)),times=1000)
#  microbenchmark::microbenchmark(C=synth_update(C,M$bas,sample(1:M$bas$Ai_max,1),sample(1:3,1)),times = 1000)
} # 1.5ms median

#microbenchmark::microbenchmark(C=updateX(M$A,M$grd,M$bas),times = 100) # 26 ms



 #a couple of checks, needs periodic integrands!!!
  grdGL=MakeGrid_GaussLegendreSimpson(120)
  # therefore, try with H2 evaluation of geometric integral quantities since integrands are periodic in v
  basGL=MakeBasis_UV(2,grdGL$U,grdGL$V)
  AGL=MakeSphere(grdGL,basGL)
  CGL=updateX(AGL,grdGL,basGL)
  M.C0=0
  H2=E_SCM(AGL,grdGL,basGL,CGL) 
  severe(H2$Volume,4/3*pi,"Volume",1e-13) # 8.6-14
  severe(H2$Area,4*pi,"Area",2e-13)    # 1.74e-13

  severe(H2$Curv,8*pi,"Curvature integral",1e-12)   # 1.522e-13
  severe(H2$H2_BC,16*pi,"Curvature_Square integral",1e-12)

    M.C0=0
    grd=MakeGrid_GaussLegendreSimpson(97)
    bas=MakeBasis_UV(18,grd$U,grd$V)
    A=MakeSphere(grd,bas)
    pertA_Gauss(A,bas,0.1)->A

    updateX(A,grd,bas)->C
    h2=E_SCM(A,grd,bas,C)

    tictoc::tic() # Ryzen7pro :: 8.2s
    g2<-Grad_SCM_R(h2,grd,bas,C )
    cat("serial Grad_SCM ")
    tictoc::toc()

    tictoc::tic() # Ryzen7pro :: 2.47 on 6 cores; 2.42 on 3 cores
    M.Rcpp=TRUE;M.Rcpp_ncores=3
    G2<-Grad_SCM(h2,grd,bas,C)

    cat("parallel Grad_SCM(",M.Rcpp_ncores,")")
    tictoc::toc()


    severe(pracma::Norm(g2$gradH2-G2$gradH2),0,"gradH compiled",3e-7)

    length(G2$ddA)==length(c(g2$ddA))

    severe(pracma::Norm(g2$ddA-G2$ddA),0,"ddA compiled",1e-12)
    severe(pracma::Norm(g2$gradC-G2$gradC),0,"gradC compiled",1e-12)
    severe(pracma::Norm(g2$gradV-G2$gradV),0,"gradV compiled",1e-12)
    severe(pracma::Norm(g2$gradA-G2$gradA),0,"gradA compiled",1e-12)

    severe(pracma::Norm(g2$dE-G2$dE),0,"dE compiled",1e-12)
    severe(pracma::Norm(g2$dF-G2$dF),0,"dF compiled",1e-12)
    severe(pracma::Norm(g2$dG-G2$dG),0,"dG compiled",1e-12)

  grd=MakeGrid_GaussLegendreSimpson(25)
  bas=MakeBasis_UV(5,grd$U,grd$V)
  M=MakeSphere(grd,bas)
  updateX(M,grd,bas)->C
  M.C0=1
  E_SCM(M,grd,bas,C)->h2
  h2$Volume/4/pi*3 # unit sphere reference
  h2$Curv/4/pi/2 # unit sphere reference curvature

  grd=MakeGrid_GaussLegendreSimpson(50)
  bas=MakeBasis_UV(1,grd$U,grd$V)
  A=MakeSphere(grd,bas)
  C=updateX(A,grd,bas)
  grd$Obj=X2Obj(grd$Obj,C$X)
  rgl::clear3d();imag.obj(grd$Obj,1)

  plot3b(C$X,grd)

  M.C0=0
  H2=E_SCM(A,grd,bas,C)

  cat(H2$Area/4/pi,H2$Curv/8/pi,H2$H2_BC/16/pi) # should be nearly 1 all three for sphere of radius 1

  dim(grd$u)
  rgl::wire3d(grd$Obj)
  rgl::shade3d(grd$Obj) # if you look carfully you see that it is not watertight around v=0 for GL; this is not so with  GLSimpson grid

  rgl::clear3d()
  #imag.obj.colorbar.simple(grd$Obj,(grd$v))  # look OK, poles identified
  rgl::clear3d()
  x11()
  imag(H2$dA,grd)
  x11()
  imag(H2$curv,grd)

  rgl::open3d()
  imag.obj(grd$Obj,(grd$u))    # OK results
  rgl::open3d()
  imag.obj(grd$Obj,(grd$V))

  grdGL=MakeGrid_GaussLegendreSimpson(15*5)
  dim(grdGL$u)
  rgl::wire3d(grdGL$Obj)
  severe((I1=int2d_scalar_GLS(grdGL$V,grdGL))-(2*pi)^2/2*pi,0,"I1-test",6e-14)
  severe((I2=int2d_scalar_GLS(grdGL$V^2,grdGL))-(2*pi)^3/3*pi,0,"I2-test",2e-13) # 1.7e-13
  severe((I3=int2d_scalar_GLS(grdGL$V*grdGL$U,grdGL))-pi^2/2*(2*pi)^2/2,0,"I3-test",9e-14) # 1.42e-14
  severe((I4=int2d_scalar_GLS(grdGL$V^3*grdGL$U^3,grdGL))-(2*pi)^4/4*pi^4/4,0,"I4-test",3e-11) # 2.182e-11

  basGL=MakeBasis_UV(15,grdGL$U,grdGL$V)
  severe(int2d_scalar_GLS(basGL$Ylm[,1],grdGL),0,"YLM1",8e-15) #-2.565e-15
  severe(int2d_scalar_GLS(basGL$Ylm[,2],grdGL),0,"YLM2",5e-15) #2.2426e-15
  severe(int2d_scalar_GLS(basGL$Ylm[,3],grdGL),0,"YLM3",2e-14) #-1.7484e-14
  severe(int2d_scalar_GLS(basGL$Ylm[,1]^2*sin(grdGL$U),grdGL)-1,0,"YLM_NORM",5e-15) #-1.7484e-14
  severe(int2d_scalar_GLS(basGL$Ylm[,basGL$Ai_max]^2*sin(grdGL$U),grdGL)-1,0,"YLM_NORM",2e-14) #-1.7484e-14
  severe(int2d_scalar_GLS(basGL$Ylm_uu[,basGL$Ai_max],grdGL),0,"Ylm_uu",3e-16)
  severe(int2d_scalar_GLS(basGL$Ylm_uv[,basGL$Ai_max],grdGL),0,"Ylm_vv",2e-15)
  severe(int2d_scalar_GLS(basGL$Ylm_vv[,basGL$Ai_max],grdGL),0,"Ylm_vv",23-12)
  severe(int2d_scalar_GLS(basGL$Ylm_u[,basGL$Ai_max],grdGL),0,"Ylm_u",1e-16)
  severe(int2d_scalar_GLS(basGL$Ylm_v[,basGL$Ai_max],grdGL),0,"Ylm_v",2e-14)


  grdGL=MakeGrid_GaussLegendreSimpson(120)
  basGL=MakeBasis_UV(2,grdGL$U,grdGL$V)
  AGL=MakeSphere(grdGL,basGL)
  CGL=updateX(AGL,grdGL,basGL)
  M.C0=0
  H2=E_SCM(AGL,grdGL,basGL,CGL)
  severe(H2$Volume-4/3*pi,0,"Vol",9e-14) # 8.08e-14
  severe(H2$Area-4*pi,0,"Area",1.6e-13)    # 1.51e-13
  severe(H2$Curv-8*pi,0,"Curv",1.5e-13)   # 1.492e-13
  severe(H2$H2_BC-16*pi,0,"BC",8.7e-14) # -8.527e-14

# TEST HESSIANS
  data(D5)
  M.C0=D5$Params[["M.C0"]] # or for all parans: SetParams(D5)
M.Rcpp_ncores=3
# serial gradients assembly version
  system.time({
    HR=FullModelHessian(D5$A,D5$grd,D5$bas,D5$Ref, del=1e-6) }) #  R, now used in CNM
# parallel grads assembly in cluster
  system.time({
    HP=FullModelHessian_Par(D5$A,D5$grd,D5$bas,D5$Ref,del=1e-6,Mem_mc.cores = 4) })#  R-cluster
  norm(HR$H-HP$H) # 2.8e-7 -> cluster-parallel problem? order of summations in C$X?

norm(HR$H-HP$H) # 2e-7 speaks for non-cluster dev due to fast synch_update of C

severe(norm(HP$H-HR$H),0,"HESSIAN",3e-7)

system.time(CNM(D5,2,cluster=FALSE,plt=FALSE)->cnmHR)
system.time(CNM(D5,2,cluster=TRUE,plt=FALSE)->cnmHP)

severe(norm(cnmHR$Hessian-cnmHP$Hessian),0,"CNM-Hessian",1e-6)


# speedup by synth_update: 26/1.5=17.3 ms
# i.e.
#(26-1.5) * M$bas$Ai_max*3/1000  # time waste in Hessian by full vs partial updates of C: 2.5 seconds for L=16

data("ss42denovo")
ss42denovo->M
system.time(FullModelHessian_Par(M$A,M$grd,M$bas,M$Ref,Mem_mc.cores = 5))

# with smaller load again:
system.time(CNM(D5,2,cluster=FALSE,plt=FALSE)->cnmHC)
system.time(CNM(D5,2,cluster=TRUE,plt=FALSE)->cnmHP)
severe(norm(cnmHC$H),norm(cnmHP$H),"CNM seq-par",1e-3) # 0.00187 -> cluster-parallel problem?

