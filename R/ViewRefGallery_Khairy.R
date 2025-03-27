# MemRBC the R package for spectral shape modeling of red blood cells
# (C) 2025 Stephan Frickenhaus

# CITATION
# when using this software for publications you must cite it as:
#  Frickenhaus, S. (2025). MemRBC - a numerical modeling laboratory for the stomatocyte-discocyte-echinocyte-transformation of Red Blood Cell shape, ZENODO, DOI: https://doi.org/10.5281/zenodo.13908340

#
# ORIGINAL CITATION (->first version)
# Frickenhaus S. (2024). MembraneR3 - A spectral model of membrane shape based on Helfrich spontaneous curvature in R. Zenodo. https://doi.org/10.5281/zenodo.13627757 ")}
#


if(FALSE)
{
library(MemRBC)
library(memoise)

MB<-memoise(MakeBasis_UV)

print(load("LA.Ref-shape_library_cs_old-full.rdat"))
print(load("LA.Ref-Shapes_Mat.rdat"))
print(load("LA.Ref-shape_library_cs_old.rdat"))

show3d=function(A)
{ L=sqrt(dim(A)[1])-1
  print(L)
  MakeGrid_GaussLegendreSimpson(L*(L+1)+2)->grd
  MB(L,grd$U,grd$V)->bas
  updateX(A[-1,],grd,bas)->C
  rgl::open3d();plot3q(C$X,grd);rgl::title3d(i); i<<-i+1
  A=A[,c(2,1,3)]
  M=structure(list(grd=grd,bas=bas,A=A[-1,]),class="MemRBC")
  return(M)
}
i=1;
LM=lapply(LA.Ref[1:10],show3d)
LC=sapply(LM,Quantities)
LC
rgl::plot3d(t(LC))
show3d(LA.Ref[[8]])
}
