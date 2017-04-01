## This script tests 
## - isotropic smoothing 
## - 1st order FEs 
## - C++ code
#setwd("/home/pacs_student/progetto_pacs/tests")
source('../R/datasets.R',chdir=T)
source('../R/fdaPDE.objects.R')
source("../R/fdaPDE.smoothing.R")
source("../R/fdaPDE.checkParameters.R")
source("../R/fdaPDE.smoothing_CPP.R")
source("../R/fdaPDE.smoothing_R.R")
source("../R/fdaPDE.locator.R")
source("../R/mesh.2D.R")
source("../R/zzz.R")
source("../R/fdaPDE.smoothing.manifold_CPP.R")
source("../R/fdaPDE.plot.mesh_R.R")
dyn.load("../src/fdaPDE.so") 



order = 1

#FEMbasis = create.FEM.basis(mesh)

lambda = c(1)
filename = 'Caramella.csv'

read.mesh<-function(filename){
  nnodes = read.table(filename, nrows = 2, header = FALSE, sep =' ', stringsAsFactors = FALSE)[1,2]
  ntriangles = read.table(filename, nrows = 2, header = FALSE, sep =' ', stringsAsFactors = FALSE)[2,2]
  geometry = read.csv(filename, skip=2, header= FALSE)
  nodes = geometry[1:nnodes,]
  triangles = geometry[(nnodes+1):(nnodes+ntriangles),]
  
  retlist = list(nnodes = nnodes,ntriangles = ntriangles,nodes = nodes,triangles = triangles)
  return(retlist)
}

mymesh=read.mesh(filename)


#locations=NULL
#locations=mymesh$nodes
locations = as.matrix(read.csv("locations.csv",sep=",",header = F))
###### Create our observations ########
nodes = read.csv("locations.csv",sep=",",header=F)
nnodes = dim(nodes)[1]
set.seed(6012017)
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)
func_evaluation = numeric(nnodes)
for (i in 1:nnodes){
  func_evaluation[i] = a1* sin(2*pi*nodes[i,1]) +  a2* sin(2*pi*nodes[i,2]) +  a3*sin(2*pi*nodes[i,3]) +1
}
data=func_evaluation
rm(nnodes,nodes)
#data = read.csv('observation_caramella_by_index.csv',header=T)[,2]
#head(data) 
covariates = NULL
BC = read.csv("bc_caramella.csv",header=T)
colnames(BC)=c("BC_indices","BC_values")
head(BC)

mesh <- create.surface.mesh(mymesh$nodes, mymesh$triangles, order=1)

FEMbasis <- create.FEM.basis(mesh)

#locations=as.matrix(locations[1008:1009,])
#data=data[1008:1009]

output_CPP = smooth.FEM.basis(locations = locations,
                              observations = data, 
                              FEMbasis = FEMbasis, lambda = lambda, 
                              BC = BC,
                              GCV = TRUE,
                              CPP_CODE = TRUE)




print(output_CPP$fit.FEM$coeff)
plot.surface.mesh(mesh,output_CPP$fit.FEM$coeff)

##############

triangles = mymesh$triangles
nodes = mymesh$nodes
punti_strani = read.csv(file="../src/punti_strani.csv",header=F,sep=",")
for (i in 1:dim(punti_strani[1])){
  points3d(punti_strani[i,1],punti_strani[i,2],punti_strani[i,3],col="red",pch=100)
}

triangles[1007,]
tail(punti_strani)
punto_1 = nodes[triangles[1008,1],]
punto_2 = nodes[triangles[1008,2],]
punto_3 = nodes[triangles[1008,3],]

A=matrix(rep(0,6),nrow=3)
A[,1]=t(punto_2-punto_1)
A[,2]=t(punto_3-punto_1)
b=t(locations[1,]-punto_1)
