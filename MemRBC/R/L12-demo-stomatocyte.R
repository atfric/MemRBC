#
# Demo L12
#  not run when building the package
if(FALSE)
{
library(MemRBC)
data("M_stomatocyte_L12")
SetParams(M_stomatocyte_L12)
M_stomatocyte_L12
data("M_stomatocyte_L12_rotated")
M_stomatocyte_L12_rotated$Params=M_stomatocyte_L12$Params

M=M_stomatocyte_L12

SetParams(M_stomatocyte_L12)
rotUV(M,pi/2,-pi/4)->M_stomatocyte_L12_rot
M_stomatocyte_L12_rot
save_MemRBC(M_stomatocyte_L12_rot,"data/M_stomatocyte_L12_rot.rda")
SetParams(M_stomatocyte_L12_rot)
M=M_stomatocyte_L12_rot

MMC(M,100000,plt=TRUE)->Mmmc

save(Mmmc,file=paste("M_stoma_L12_rot_MMC",0,".rda"))

load("M_stoma_L12_rot_MMC 0 .rda")

for (i in 1:10){
 MMC(Mmmc,10000,plt=TRUE)->Mmmc
 save(Mmmc,file=paste("M_stoma_L12_rot_MMC",i,".rda"))
}
PlotSample(Mmmc)
}
