#
# all tests
#
roxygen2::roxygenise()

#L=5 tests
data(D5)
SetParams(D5)
D5$mask

D5
D5$mask=NULL
D5$bas$mask=double_uv_ind(D5$grd$U,D5$grd$V)
usethis::use_data(D5,overwrite = TRUE)

PSD(D5,1)->D5.1

M.C0=-2
MMC(D5.1,1,C0=-2)->D5.2
#GAM(D5.2,1)->D5.2.gam
PNEM(D5.2,2)->D5.3
CNM(D5.3,1)->D5.4
MMCC(D5,curv = Curv(D5)+10,nsteps = 10000,plt=TRUE,pltfreq = 50) -> D5c

D5c$bas$mask=D5$bas$mask
M.nn=NNuv(D5$grd$UV,13)
D5c$bas<-MakeIM(D5c$bas)
PlotSample(D5c)


MMCC(D5c,pertA=pertA_complex,nsteps = 2000,
     curv=Curv(D5)+10,plt=TRUE,pltfreq=50)->D5c1
plot(D5c1)

PNEM(D5c1,220,viscosity = 50, dt=1e-4)->D5c1
PNEM(D5c1,200,viscosity = 50, dt=1e-4,zero_Av = TRUE)->D5c1
for (i in 1:3){
 PNEM(D5c1,200,viscosity = 150, dt=1e-3,zero_Av = TRUE)->D5c1
 PlotPNEM(D5c1)}

plot(D5c1)
image(D5c1)

M.C0=-6
MMC(D5c1,100,plt=TRUE,pltfreq = 50, C0=-6)->D5c2
MMC(D5c2,100,plt=TRUE,pltfreq = 50, C0=-6,
    pertA=pertA_complex,sd=0.2)->D5c3
image(D5c3)

PSDC(D5c3,curv=Curv(D5),100)->D5c4
plot(D5c4)
image(D5c4)

CNM(D5c4,10)->D5c5
image(D5c5)

plot(D5c5)
usethis::use_data(D5c5)

# starfish data
data(SF4)
SF4$bas$L_max
plot(SF4)
Energy(SF4)

M.C0=11
M.Ka=0
M.mu=0
StoreParams(SF4)->SF4
PSD(SF4,30,del=1e-7)->SF41
PNEM(SF41,1,dt=1e-4)->SF41
PNEM(SF41,10,dt=1e-4)->SF42
PlotPNEM(SF42)
PNEM(SF42,100,dt=2e-4)->SF42
plot(SF42)
PSD(SF42,100,plt=TRUE,pltfreq = 5,del=1e-7)->SF43

# make new mems
MakeStandardRBC(L=15)->R15
plot(R15)
Wb<-E_SCM(R15$A,R15$grd,R15$bas,updateX(R15$A,R15$grd,R15$bas))
S<-SEN(R15$A,R15$grd,R15$bas,R15$Ref,Wb)
R15$SEN<-S
imag(S$alpha,R15$grd)
image(R15,R15$SEN$beta)
R15
M.C0=-2
M.Ka=5
M.mu=2.5
save_MemRBC(R15,"R15.rdat")
usethis::use_data(R15,overwrite = TRUE)

# other routines
MakeGrid_GaussLegendreSimpson(40)->g
MakeBasis_UV(12,g$U,g$V)->b

updateX(SF4$A,SF4$grd,SF4$bas)->C
SF4$grd$Obj<-X2Obj(SF4$grd$Obj,C$X)
rgl::plot3d(SF4$grd$Obj,aspect=FALSE)

Energy(SF4)
E_SCM(SF4$A,SF4$grd,SF4$bas,updateX(SF4$A,SF4$grd,SF4$bas))->E
imag.obj.colorbar(SF4$grd$Obj,E$L)
imag.obj.colorbar(SF4$grd$Obj,E$dA)
imag.obj.colorbar(SF4$grd$Obj,E$curv_sq)
imag.obj.colorbar(SF4$grd$Obj,E$curv)
imag.obj.colorbar(SF4$grd$Obj,E$E_SCM_dens)
rgl::contourLines3d(SF4$grd$Obj,SF4$grd$U,30,lwd=2)

#CNM(SF4,1)
updateX(SF4$A,SF4$grd,SF4$bas)->C

H2=E_SCM(SF4$A,SF4$grd,SF4$bas,C)
Energy(SF4)
Membrane_Curvatures(H2,SF4$grd,plt.K = TRUE)->K
imag(K$k1,SF4$grd)
length(K$k1)

# other stuff
load_MemRBC("SS.rdat")->SS
SS
usethis::use_data(SS)
ss20=MakeSpiculated(20, w = 0.06, d = -0.6, r = 0.65)
data(SS20)
max(abs(ss20$A-SS20$A))
plot(ss20)
plot(SS20)
data(SS42)
plot(SS42)
SS20$history



data(D5); SetParams(D5)
SetConstraints(D5$bas,Cons=c("gradA","gradV","gradC"),
QCons=c("Area","Volume","Curv"),
Target=c(140,100,121),
TNorm=c(140,100,121)) -> D5$bas  # store modified basis back into membranes D5 basis

 # minimize with steepest descend under Rosen Constraint Projection

 SDRC(D5,100)->D5sdrc
 SDRC(D5sdrc,10000)->D5sdrc

 # pair-plots of target quantities and energy E
 plot(D5sdrc$SDRC_Sample[c("E","A","V","C")])

#### Unduloid case, lowres, L15
 open=0.4 # also higher work, but objects boundary is never fitted
  g<-MakeGrid_GaussLegendreSimpson(60,ua=open,ub=pi-open)
  g$ndof
  range(g$U)
  (U<-TriMesh_Unduloid(periods=2,nx=g$nu,ny=g$nv-1, a=1,c=0.1,clean=FALSE))
  attr(U,"H_theor") # theoretical mean curvature
  attr(U,"H_vcg_6") # mean from vertices with 6 neighbors
  #  mean curvature for other vertices is problematic in vcg
  b<-MakeBasis_UV(15,g$U,g$V)
  # one should exclude double coordinates for
  # the fit, so mask is needed
  b$mask<-double_uv_ind(b$uv[,1],b$uv[,2])
  CenterX(Obj2X(U)) -> X
  cat(dim(X)[1],"?=", g$ndof,"\n")
  if (dim(X)[1] == g$ndof)
  { # if not matching repeat with alternative n in grid
  rgl::plot3d(X,col=2,aspect=FALSE)
  X2Obj(U,X) -> U1
  rgl::contourLines3d(U1,b$uv[,1],levels=(0:100)*pi/100)
  rgl::contourLines3d(U1,b$uv[,2],40)
  rgl::shade3d(U1,col="grey",alpha=0.5)
  A<-FitAlm_Tikhonov(X,b,lambda=0) # , WX=sin(g$U))
  A[,3]<--A[,3] # wrong orientation correction
  # you may try the fit with weights, WX=sin(g$U)
  MakeMemRBC(A,g,b)->M
  rgl::open3d()
  rgl::plot3d(X,aspect=FALSE,alpha=0.45)
  plot(M,alpha=0.45,col="cyan",wire=FALSE)
  E=E_SCM(M$A,M$grd,M$bas,updateX(M$A,M$grd,M$bas))
  -E$Curv/E$Area/2 # mean curvature from integral over area
  attr(U,"H_vcg_6") # comparison with original
  imag.obj.colorbar(U1,E$curv)
  rgl::title3d("curvature density")
  X1<-updateX(M$A,M$grd,M$bas)$X
  M$grd$Obj<-X2Obj(M$grd$Obj,X1)
  rgl::shade3d(M$grd$Obj,alpha=0.2)
  mean(E$curv/E$dA/2)
  M.C0=0;M.mu=0;M.Ka=0
  M$bas$Target[1:2]=c(E$Area,E$Volume)
  save_MemRBC(M,"Unduloid_mmc.rdat")
  MMC(M,10000,plt=TRUE,pltfreq=10,C0=0) ->M
  }
