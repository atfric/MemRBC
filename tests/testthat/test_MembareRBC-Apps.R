

M.C0=-2 # for discoid shape

data(M4)
rgl::clear3d();plot(M4)

update(M4)->M4

severe(M4$Curv-Curv(M4),"update(M)/ Curv")
severe(length(M4$dA)-M4$grd$ndof,"length of dA from update(M)")
M4$SEN=NULL
update(M4,"SEN")->M4

image(M4,M4$SEN$alpha)
image(M4,M4$SEN$beta)

rgl::clear3d();imag.obj.colorbar(M$grd$Obj,f = M$SEN$alpha)

plot(M4)

print(M4)

image(M4)

M.Rcpp
M.Rcpp_ncores=3
MemStab(M4,serial=FALSE)->M4a

PlotStabGallery(M4a)

Quantities(M4)

Volume(M4)
Area(M4)
Curv(M4)

two_screens3d()
two_draw3d(M4$A,M4)

M4pca<-update(MemPCA(M4))
two_draw3d(M4pca$A,M4pca)

plot(M4)
plot(M4pca)

PNEM(M4,100)->M4a
PNEMVM(M4a,100)->M4b
PlotPNEM(M4b)

rotUV(M4b,du=pi/2,dv=pi/4)->M4b_rot
plot(M4b_rot)
