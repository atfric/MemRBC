#
# Torus demo
#

SetParam("M.C0",1)
SetParam("M.mu",0)
SetParam("M.Ka",0)
print_param_RBC()


M<-MakeTorus(L = 5, r = 1.75)
T5_mmc<-MMC(M, nsteps = 10000, C0 = 1)
T5_mmc<-MMC(T5_mmc, nsteps = 90000, C0 = 1)
save_MemRBC(T5_mmc,"../T5_mmc",qs2=TRUE)

T5_mmc_sdrc<-SDRC(T5_mmc, nsteps = 500, cons_tol = 0.25)
save_MemRBC(T5_mmc_sdrc,"../T5_mmc",qs2=TRUE)

ct=0.25
for (i in 1:50)
{ct=ct*0.7
 T5_mmc_sdrc<-SDRC(T5_mmc_sdrc, 
                   nsteps = 1000,
                   Gtol = 1e-3
                   cons_tol = ct,
                   maxiter=150)
 save_MemRBC(T5_mmc_sdrc,"../T5_mmc",qs2=TRUE)
}