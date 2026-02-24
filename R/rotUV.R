# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2025 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, S. (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# Frickenhaus S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#


#
# rotUV - uv rotation to re-orient poles
#    possible to orient poles to Z-Axis?
#    just find the right rotation angles du,dv
#  example:
# rotUV(M_stomatocyte_L12,pi/2,3*pi/4)->M
#

#' rotUV
#'
#' rotate membrane coordinates and re-compute coefficients
#' without modifying the SEN reference.
#' rotUV may help to move poles to Z-axis or for more homogeneous distribution of spatial points.
#'
#' @param Min membrane object to rotate
#' @param du angles of rotation relative to Z-axis
#' @param dv angle of rotation around Z-axis
#' @param plt (=FALSE) TRUE for 3d plotting
#' @param transpose (=FALSE) for backward rotation, not verified
#' @return membrane MemRBC object with coefficients after rotation, but original SEN reference
#' @examplesIf exists("L_Ylm")
#' data("M4",package = "MemRBC")
#' plot(M4)
#' rotUV(M4,pi/2,3*pi/4)->M
#' plot(M)
#' 
#' @export
rotUV<-function(Min,du,dv,plt=FALSE,transpose=FALSE)
{
 mask=Min$bas$mask
 if (is.null(mask))mask=double_uv_ind(Min$bas$uv[,1],Min$bas$uv[,2])
 A=Min$A
 a=du
 M=matrix(c(cos(a), 0,  sin(a),
            0, 1, 0,
         -sin(a), 0, cos(a)),3,3)
 b=dv
 M1=matrix(c(cos(b), sin(b),   0,
           -sin(b), cos(b), 0,
           0, 0, 1), 3, 3)
# rotate angles (u,v) on sphere
 MakeSphere(Min$grd,Min$bas)->AS
 Cs=updateX(AS,Min$grd,Min$bas)
 Xuv=Cs$X
 Xuvp=Xuv;
 if (!transpose) for (i in 1:Min$grd$ndof) Xuvp[i,]=M1%*%M%*%Xuv[i,] else for (i in 1:Min$grd$ndof) Xuvp[i,]=t(M1%*%M)%*%Xuv[i,]
 uvp=t(apply(Xuvp,1,inv_sph)) # back from X on sphere to angles uvp
#plot(uvp,pch=".")
# make rotated basis
 MakeBasis_UV(Min$bas$L_max,uvp[,1],uvp[,2]) -> bas1
 C=updateX(Min$A,Min$grd,Min$bas) # original
 Y=bas1$Ylm
 Cp=updateX(Min$A,Min$grd,bas1) # new coords
 Xn=Cp$X
 X=C$X
#rgl::plot3d(X,aspect=FALSE);
#rgl::plot3d(Xn,col=2,aspect=FALSE)
 lm(Xn[-mask,] ~ Y[-mask,] , weights = sin(Min$grd$U)[-mask])$coefficients[-1,]->Ar
 Cn=updateX(Ar,Min$grd,bas1)
 if (plt) rgl::plot3d(Cn$X,col=1,aspect=FALSE)
 Xn=Cn$X
 # backrotation :
 Xnp=Xn;if (!transpose) for (i in 1:Min$grd$ndof) Xnp[i,]=t(M)%*%t(M1)%*%Xn[i,] else for (i in 1:Min$grd$ndof) Xnp[i,]=t(t(M)%*%t(M1))%*%Xn[i,]
# rotated basis not useful for integration!
# ->   needs fitting rotated X back with original basis:
 lm( Xnp[-mask,] ~  Min$bas$Ylm[-mask,] , weights = sin(Min$grd$U)[-mask])$coefficients[-1,]->Arr

  Arr[,2]=-Arr[,2]
 M=Min;
 M$A=LM2A(Arr,Min$bas);
 print(unlist(Quantities(M)))
 M$comment="rotated"
 M$history=append(M$history,match.call())
 if (plt) {rgl::open3d();plot(Min,alpha=0.5,col="red");plot(M,alpha=0.6)}
 cat("relative error of rotation on Quantities:\n")
 print((unlist(Quantities(M))-unlist(Quantities(Min))) / (unlist(Quantities(M))+unlist(Quantities(Min))))

 return(M)
} # rotUV(M_stomatocyte_L12,pi/2,3*pi/4)->M
