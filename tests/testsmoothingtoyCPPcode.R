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
dyn.load("../src/fdaPDE.so") 



order = 1

#FEMbasis = create.FEM.basis(mesh)

lambda = c(0.1)

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


locations = NULL
data = read.csv('observation_caramella_by_index.csv',header=T)[,2]
head(data) 
covariates = NULL
BC = read.csv("bc_caramella.csv",header=T)
colnames(BC)=c("BC_indices","BC_values")
head(BC)

mesh <- create.surface.mesh(mymesh$nodes, mymesh$triangles, order=1)

FEMbasis <- create.FEM.basis(mesh)

output_CPP = smooth.FEM.basis(observations = data, 
                              FEMbasis = FEMbasis, lambda = lambda, 
                              BC = BC,
                              GCV = TRUE,
                              CPP_CODE = TRUE)

print(output_CPP$fit.FEM$coeff)


