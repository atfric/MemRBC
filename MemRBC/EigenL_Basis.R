#
#
# Eigenvectors of cotan laplacian as basis
#
#


data(D5)
D5$grd$Obj->O
L<-GEMINI_cotan_laplacian_II(O)
eigen(L)$vectors->b
eigen(L)$values->v
plot(v)
bas=D5$bas
updateX_only(D5$A,D5$grd,D5$bas)->C
b=cbind(1,b)
m=lm(C$X ~ b - 1)
plot(sort(log(abs(c(m$coefficients)))),main="sorted log abs coefficients")
plot((log(abs(c(m$coefficients)))),main="unsorted log abs coefficients")

n=dim(b)[2]

X1=Obj2X(O)
thresh=0.1
co=m$coefficients
f=rep(FALSE,dim(co)[1])
for (i in 1:dim(co)[1])
  f[i] = any( abs(m$coefficients[i,]) > thresh )
f=f & !is.na(f)
which(f)

cat("leftover: [", sum(f),",3 ]\n")
cat("compare A:",dim(D5$A),"\n")

for (k in 1:3){
  X1[,k]=m$coefficients[f,k][1] %*% b[,-f][,1]
   for (l in 2:sum(f)) X1[,k]=X1[,k] + m$coefficients[f,k][l] %*% (b[,f][,l])
}
plot(D5,alpha=0.5)
rgl::plot3d(X1,col=3,aspect=FALSE,add=TRUE)

min(which(f))
n-min(which(f)) # to be kept basis vectors, compared to 35 Ylm

     S=data.frame(k=0,S=0)
     k=0
     for (i in seq(1,n-50,by=10)){
       m=lm(C$X ~ b[,i:n])
       k=k+1; S[k,]=c(n-i,sum((m$fitted.values-C$X)^2))
       cat(i/n,"\r")
     }
     plot(S$k,log10(S$S))
     dim(D5$A)
     m$fitted.values

     plot(D5,alpha=0.35)
     rgl::plot3d(m$fitted.values,col=2,add=TRUE)
print(sum(log10(S$S)< -3.66)) # can drop last 16*10 from SSE curve
n-sum(log10(S$S)< -3.66)*10
# keep 897 of n = 84%

m=lm(C$X ~ b[,(n-160):n])
plot(D5,alpha=0.5)
rgl::plot3d(m$fitted.values,col=3,aspect=FALSE,add=TRUE)

m=lm(C$X[-bas$mask,] ~ b[-bas$mask,(n-35):n])
plot(D5,alpha=0.5)
rgl::plot3d(m$fitted.values,col=3,aspect=FALSE,add=TRUE)

