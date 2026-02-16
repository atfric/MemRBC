#
# speed talk demonstrator
#

data("SS20")

SS20$Params

SetParams(SS20)

M.C0=24

plot(SS20)
Quantities(SS20)

MMC(SS20,250,plt=TRUE,pltfreq = 10) -> SS20_1

PlotSample(SS20_1)

SS20_1

SS20_1$history

rotUV(SS20_1,-0.45,0,plt = TRUE)->SS20_2
plot(SS20_2)

History(SS20)

S42 <- load_MemRBC("../Collection-MemRBC/SS16_mmc_C035_rho5.rdat")
PlotSample(S42)

plot(S42+SS20)




data("L5_stomatocyte_equilib")
SetParams(L5_stomatocyte_equilib)
L5_stomatocyte_equilib -> M # sets C0 for next computations
M.Rcpp=TRUE
M.Rcpp_ncores=3
sink(file="PSD_L5s_250.txt")
PSD(M,250,del=1.5e-5,plt=FALSE)->Mpsd
sink()

sink(file="PNEM_L5s_250.txt")
PNEM(M,250,dt=0.01,plt=FALSE, mass_update_freq = 1500, zero_Av=TRUE)->Mpnem
sink()

E0<-Energy(M)["E"]

Energy(Mpsd)["E"]/M.Es-E0/M.Es
Energy(Mpnem)["E"]/M.Es-E0/M.Es

Mpsd$proc_time-M$proc_time # 25sec
Mpnem$proc_time-M$proc_time # 24sec

plot(last(Mpnem$E_total_PNEM-Mpnem$E_kin_PNEM,250)/M.Es,type="l",col=2,xlab="steps",ylab="E")
points(last(Mpsd$E_PSD/M.Es,250),type="l")
legend("topright",lwd=2,col=1:2,c("PSD","PNEM"),cex=0.45)

CNM(M,150,del=1.6e-7,plt=TRUE)->Mcnm
Mcnm$proc_time-M$proc_time # 2227 s
S5cnm=Mcnm1
usethis::use_data(S5cnm)

plot(last(Mpnem$E_total_PNEM-Mpnem$E_kin_PNEM,250)/M.Es,type="l",col=2,xlab="steps",ylab="E",ylim=c(0.52,0.75))
points(last(Mpsd$E_PSD/M.Es,250),type="l")
points(last(Mcnm$E_CNM/M.Es,250),type="l",col=3)

Mcnm1$proc_time-M$proc_time # 711 s

legend("topright",lwd=2,col=1:3,c("PSD","PNEM","CNM"),cex=1)

CNM(M,50,del=1e-7,plt=TRUE)->Mcnm1
plot(last(Mcnm1$E_CNM/M.Es,50),type="l",col=3)

plot(last(Mpnem$E_total_PNEM-Mpnem$E_kin_PNEM,250)/M.Es,type="l",col=2,xlab="steps",ylab="E",ylim=c(0.52,0.75))
points(last(Mpsd$E_PSD/M.Es,250),type="l")
points(last(Mcnm1$E_CNM/M.Es,50),type="l",col=3)
legend("topright",lwd=2,col=1:3,c("PSD","PNEM","CNM"),cex=1)
Mcnm1$proc_time-M$proc_time # 711 s

PNEM(Mpnem,250,dt=0.01,plt=FALSE, mass_update_freq = 50, zero_Av=TRUE)->Mpnem1

plot(Mpnem3$E_total_PNEM-Mpnem3$E_kin_PNEM/M.Es,type="l",col=2,xlab="steps",ylab="E")
plot(Mpnem1)

2227/35 # (Hessian gradients) = 63 sec.

data(D5)
bas=D5$bas
mask=1
X=Obj2X(D5$grd$Obj)
WX=rep(1,nrow(X))[-mask]
X1<-X[-mask,]
Y=cbind(1,bas$Ylm[-mask,])
YtW=t(Y)%*%diag(WX)
B=YtW %*% Y
InvB1=pracma::inv(B + lambda * diag(c(1,bas$G.tk)))
IM = InvB1 %*% YtW
A=bas$A
for (k in 1:3) A[,k] = (IM %*% X1[,k]) [-1,]

FitAlm_Tikhonov(X,bas,0)->A1


library(torch)


library(torch)

# Define the model architecture
multivariate_regressor <- nn_module(
  "MultivariateRegressor",

  initialize = function(input_size) {
    self$layer1 <- nn_linear(input_size, 64)
    self$layer2 <- nn_linear(64, 32)
    self$output_layer <- nn_linear(32, 1)
    self$relu <- nn_relu()
  },

  forward = function(x) {
    x %>%
      self$layer1() %>%
      self$relu() %>%
      self$layer2() %>%
      self$relu() %>%
      self$output_layer()
  }
)

input_dim <- dim(dA)[1]
model <- multivariate_regressor(input_dim)
criterion <- nn_mse_loss()
optimizer <- optim_adam(model$parameters, lr = 0.01)

# 2. Dummy Data (Creating Tensors)
# R matrices need to be converted to torch_tensor
x_train <- torch_tensor(t(dA))
y_train <- torch_tensor(dE)

# 3. Training Loop
for (epoch in 1:100) {
  # Forward pass
  preds <- model(x_train)
  loss <- criterion(preds, y_train)

  # Backward pass and optimization
  optimizer$zero_grad()
  loss$backward()
  optimizer$step()

  if (epoch %% 10 == 0) {
    cat(sprintf("Epoch %d, Loss: %.4f\n", epoch, loss$item()))
  }
}

preds-dE

