# there are copies of M1 and M4 in the testhat folder
load("M4.rda")
update(M4,c("Quantities","dA"))->M4
severe(M4$Curv,Curv(M4),"update(M)/ Curv")
severe(length(M4$dA),M4$grd$ndof,"length of dA from update(M)")
M4$SEN=NULL
update(M4,"SEN")->M4

image(M4,M4$SEN$alpha)
image(M4,M4$SEN$beta)

#rgl::clear3d();imag.obj.colorbar(M$grd$Obj,f = M$SEN$alpha)

print(M4)
image(M4)

M.Rcpp
M.Rcpp_ncores=3
MemStab(M4,serial=TRUE)->M4a

PlotStabGallery(M4a)

Quantities(M4)

Volume(M4)
Area(M4)
Curv(M4)
Energy(M4)

two_screens3d()
two_draw3d(M4$A,M4)
data("M4",package = "MemRBC")
M4$bas$mask = double_uv_ind(M4$bas$uv[,1],M4$bas$uv[,2])
M4$history
M4pca<-MemPCA(M4)
two_draw3d(M4pca$A,M4pca)

plot(M4pca)

PNEM(M4,30)    -> M4a
PNEMVM(M4a,30) -> M4b
PlotPNEM(M4b)
PlotDiff(M4,M4b,On3d = TRUE)

rotUV(M4b,du=pi/2,dv=pi/4)->M4b_rot
plot(M4b_rot)

CNM(M4b_rot,2,del = 1e-8)->M4cnm

PSD(M4b,5,del=5e-6)->M4c
SDRC(M4c,20)->M4d
MemRBC_env$M.C0<<-0
MMC(M4d,100,C0=MemRBC_env$M.C0)->M4e

system("touch STOP_MMC.txt")
MMC(M4d,100,C0=MemRBC_env$M.C0)->M4e

MMCC(M4e,curv=Curv(M4e),100)->M4f
plot(M4f)

PlotLSeries(M = M4f,nr = 2,nc=3)

PSDC(M4f,curv=Curv(M4f)+0.1,nsteps=20)->M4g
plot(M4g)

MemRBC_env$M.muk<<-MemRBC_env$M.lam<<-0.1
MemRBC_env$M.C0<<- -2
# no test for CRAN here:
#ALM(M4g,Curv(M4)+0.2,5,ncores=2)->M4h

data("SF4lr",package = "MemRBC")
SF4mr<-rgl::subdivision3d(SF4lr,1)
rgl::plot3d(SF4mr,aspect=FALSE,col="white")

NorthSouth(SF4mr)->P
GEMINI_cut_mesh_along_path(SF4mr,P$Paths[[1]][[1]]) -> SF4_open
GEMINI_disk_conformal_map(SF4_open,spherical = TRUE,plt=TRUE)->uv 
plot(uv$uv,pch=".",cex=2)

L=Membrane_LaplacianOBJ(SF4lr)
dim(L)
diag(as.matrix(L))
EV=Membrane_Eig(L,which = 1:10)
rgl::shade3d(SF4lr,col=MemCols(EV[,10],pal=rainbow))

L1=Membrane_Laplacian_cotan(SF4lr,L)
L2=GEMINI_cotan_Laplacian_II(SF4lr)
class(L2)
dL1=as.matrix(L1)
dL2=as.matrix(L2)
severe(norm(dL1-dL2),0,"Laplacians",1e-14)

L1=L2=M4
# rotate coeffs by 90deg
L1$A[,1:3]=L1$A[,c(1,3,2)]
rgl::open3d()
plot(L1,col=2)
plot(L2,col=3)
GEMINI_Intersect_Mem_Mem(L1,L2)

rgl::open3d()
plot(L1-0.5*L2,col=3)
GEMINI_Intersect_Mem_Mem(L1+0.5*L2)

