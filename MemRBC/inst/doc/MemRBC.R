## ----message=FALSE, warning=FALSE, include=FALSE, paged.print=FALSE-----------
knitr::opts_chunk$set(  collapse = TRUE,  comment = "#>"
                      , out.width = "100%", out.height="50%")


## ----start, echo =TRUE--------------------------------------------------------
library(MemRBC)

load_param_RBC() # usually already done on start-up
data("M1")
M1
SetParams(M1)

# modify target volume in basis
M1$bas$Target["Volume"]=95
save_MemRBC(M1,"M1_v95.rdat") # saves object with modified $Target
save_MemRBC(M1,"M1_v95_full_basis.rdat",reduce_basis = FALSE) # saves object with modified $Target and full basis, i.e. basis need not be recomputed on load_MemRBC.


## ----plotfirst,webgl=TRUE-----------------------------------------------------
plot(M1)
rgl::rglwidget() # only needed inside html-vignette

## ----prepareMMC, message=TRUE, warning=TRUE, webgl=TRUE-----------------------
MakeStandardRBC(L=7,prn=TRUE)->M
rgl::rglwidget(width=300,height=250)

## ----MMCCcompute, webgl=TRUE--------------------------------------------------
sink("MMCC_L7.txt") # redirect console output to file
M.C0=-4
for (i in 1:3){ # loop may go upto 40 
    MMCC(M,curv=Curv(M)+1,nsteps = 100 ,plt=FALSE,LAfreq=10)->M
  }
save_MemRBC(M,"M_L7_to_Echinocyte.rdat")
sink() # terminate sink()
plot(M) # only plot last
rgl::rglwidget(width=300,height=250)

## ----plot,fig.dim = c(6, 6)---------------------------------------------------
PlotSample(M) # can have differente colors for multiple runs

## ----stomatocyte--------------------------------------------------------------
data("M_stomatocyte_L12")
SetParams(M_stomatocyte_L12)
Curv(M_stomatocyte_L12)


## ----MMCstomatocyte-----------------------------------------------------------
sink("MMC_Stomatocyte.txt")
M.C0=-6.5
MMCC(M_stomatocyte_L12,curv=70,nsteps=50,sd=0.001,plt=FALSE)->M12
sink()
M12

## ----plotM12, webgl=TRUE------------------------------------------------------
rgl::clear3d();
plot(M12)
rgl::rglwidget(width=300,height=300)  

## ----LSeries, webgl=TRUE, echo=FALSE------------------------------------------
PlotLSeries(M12,nr=3,nc=4)->null
rgl::rglwidget(width=500,height=500)

## ----stomatocyte reduced model, fig.dim = c(6, 4)-----------------------------
data("L5_stomatocyte_equilib")
SetParams(L5_stomatocyte_equilib)
L5_stomatocyte_equilib -> M # sets C0 for next computations
M.Rcpp=TRUE
M.Rcpp_ncores=3
sink(file="PSD_L5s_250.txt")
PSD(M,250,del=1.5e-5,plt=FALSE)->Mpsd #reduced energy
sink()

sink(file="PNEM_L5s_250.txt")
PNEM(M,250,dt=0.01,plt=FALSE, mass_update_freq = 1500, zero_Av=TRUE)->Mpnem 
sink()

E0<-Energy(M)["E"]

Energy(Mpsd)["E"]/M.Es-E0/M.Es
Energy(Mpnem)["E"]/M.Es-E0/M.Es

Mpsd$proc_time-M$proc_time
Mpnem$proc_time-M$proc_time

plot(last(Mpnem$E_total_PNEM-Mpnem$E_kin_PNEM,250)/M.Es,type="l",col=2,xlab="steps",ylab="E")
points(last(Mpsd$E_PSD/M.Es,250),type="l")
legend("topright",lwd=2,col=1:2,c("PSD","PNEM"),cex=0.45)


## ----plots, webgl=TRUE--------------------------------------------------------
rgl::open3d();

imag.delta.aligned(M$A,Mpnem$A,M$grd,M$bas,shade2=TRUE )->O
rgl::title3d("red: PNEM-minimized; grey: initial shape")

rgl::rglwidget(width=300,height=300)

## ----imagequnatity2d,fig.dim = c(6, 3)----------------------------------------
data(L9_Stomatocyce_6); SetParams(L9_Stomatocyce_6)
update(L9_Stomatocyce_6,"SEN")->M
imag( TotalEnergyDensity(M$SEN),M$grd)
title("Energy density of L9 stomatocyte")

## ----parMMC-------------------------------------------------------------------
if (FALSE){
parallel::makeCluster(4) -> cl
startup<-parallel::clusterEvalQ(cl,library("MemRBC"))
# give each process an identifier tid:
startup<-parallel::clusterApply(cl, 1:4, function(x){ tid <<- x})
startup<-parallel::clusterEvalQ(cl,data("M4"))
startup<-parallel::clusterEvalQ(cl,SetParams(M4))
# set different C0, using tid and a vector of 4 values
startup<-parallel::clusterEvalQ(cl,M.C0<-c(-6,-4,-2,0)[tid])
#report that C0 is set
unlist(parallel::clusterEvalQ(cl,M.C0))
#perform only 30 mmc steps; 
M4$C0=M.C0
L_MMC<-parallel::clusterEvalQ(cl,M<-MMC(M4,30,C0=M.C0))
# save results
endup<-parallel::clusterEvalQ(cl,save_MemRBC(M,paste("M4-MMC-par-",tid,".rdat",sep="")))

#L_MMC contains all M objects from the cluster processes
sapply(L_MMC,Energy) 
}

## ----CNM, fig.dim = c(6, 4)---------------------------------------------------
data("S4");
SetParams(S4) # remember to set C0 and others correctly, e.g. to M.C0
S4
CNM(S4,5,plt=FALSE)->S4_cnm
plot(S4_cnm$CNM_data[,1]/M.Es,ylab=expression(E[CNM]),xlab="CNM iteration")

## ----rotate,webgl=TRUE--------------------------------------------------------
data("M_stomatocyte_L12")
SetParams(M_stomatocyte_L12)
rotUV(M_stomatocyte_L12,pi/2,3*pi/4)->M1
rgl::clear3d()
plot(M1)
rgl::rglwidget()
rgl::title3d("pole shifted shape")
# one may be interested to see the effects -
#  since the SEN referene shape is not rotated.
# if you have time and energy, try:
# PSD(M1,5000,plt=TRUE)->M2
# M2 # Energy 1.39
# CNM(M2,20)


## ----StomatoRotDiffm,webgl=TRUE-----------------------------------------------
data("M_stomatocyte_L12_rotated")
SetParams(M_stomatocyte_L12_rotated) # needed after data() to take over parameters
Energy(M_stomatocyte_L12_rotated)

data("M_stomatocyte_L12")
SetParams(M_stomatocyte_L12) # needed after data() to take over parameters
#plot(M_stomatocyte_L12)
Energy(M_stomatocyte_L12)

# norm of coefficient differences
pracma::Norm(M_stomatocyte_L12$A-M_stomatocyte_L12_rotated$A)

# plot SEN data alpha and beta
two_screens3d()
two_draw3d(M_stomatocyte_L12_rotated$A, M_stomatocyte_L12_rotated)
#rgl::rglwidget()



