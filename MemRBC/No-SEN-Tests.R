
roxygen2::roxygenise()

M <- load_MemRBC("../Collection-MemRBC/M-sf4-L16.rdat")
M$A[,3]<- -M$A[,3]
Quantities(M)
M$bas$Target=c(4*pi, 2)
names(M$bas$Target)=c("Area", "Volume")
# no SEN possible
M$SEN=FALSE
M$Ref=NULL
M$S=NULL

Energy(M)

A=M$A;grd=M$grd;bas=M$bas
Ref=M$Ref

FM<-E_FullModel_Penalty_AV(A,grd,bas,Ref)
any(is.na(FM))
FM$n->NN
FG<-Grad_FullModel_Penalty_AV(A,grd,bas,S=NULL,Ref)
E_SEN(A,grd,bas,Ref)
any(is.na(FG))
M.C0=0;
M.mu=0;
M.K_b=0;
M.es=Energy(M)
PSD(M,20,del = 1e-9)->M_psd
PSD(M_psd,200,del = 1e-9,plt = TRUE,pltfreq = 10)->M_psd
PNEM(M,1)->M
PNEM(M,1)->M

MMC(M,10000,plt=TRUE,pltfreq = 100)->M_mmc
MMC(M_mmc,10000,plt=TRUE,pltfreq = 100,kT=0.1)->M_mmc_hot
MMC(M_mmc_hot,10000,plt=TRUE,pltfreq = 100,kT=0.1)->M_mmc_hot

PlotSample(M_mmc)
PlotPSD(M_psd)
plot(M_psd)

S=MakeSphere(grd,bas)
O=M$grd$Obj
X2Obj(O,updateX_only(S,M$grd,M$bas)$X)->O
rgl::wire3d(O)



data(D5)
D5$Ref=NULL

FullModelHessian(D5$A,D5$grd,D5$bas,Ref=D5$Ref,del=1e-5,Ctol=1e-3,nm=FALSE,SEN=FALSE)

CNM(D5,2,plt=FALSE) # calls SCM_SEN_cxx; no Ref==NULL treatment; see c++ code todo in MembraneRBC.R

# so set M.mu=M. =0
data(D5)
M.mu=0;M.Ka=0;StoreParams(D5)->D5
CNM(D5,2,plt=FALSE) # calls SCM_SEN_cxx; no Ref==NULL treatment

#
SelfIntersect <- function(M){
O=M$grd$Obj
X2Obj(O,updateX_only(M$A,M$grd,M$bas)$X)->O
rgl::wire3d(O)
return(Obj_Obj_Intersect(O,O))
}
SelfIntersect(M_psd)


