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
dyn.load("../src/fdaPDE.so") 



order = 1

#FEMbasis = create.FEM.basis(mesh)

lambda = c(0.01,0.1,0.5)

filename = 'Caramella.csv'

read.mesh<-function(filename){
  nnodes = read.table(filename, nrows = 2, header = FALSE, sep =' ', stringsAsFactors = FALSE)[1,2]
  ntriangles = read.table(filename, nrows = 2, header = FALSE, sep =' ', stringsAsFactors = FALSE)[2,2]
  geometry = read.csv(filename, skip=2, header= FALSE)
  nodes = geometry[1:nnodes,]
  triangles = geometry[(nnodes+1):(nnodes+ntriangles),]
  
  retlist = list(nnodes=nnodes,ntriangles=ntriangles,nodes=c(t(nodes)),triangles=c(t(triangles)))
  return(retlist)
}

mesh=read.mesh(filename)


locations = NULL
data =read.csv("observation_caramella_by_index.csv", header=T)[,2]
head(data) 
covariates = NULL
BC = read.csv("bc_caramella.csv",header=T)
colnames(BC)=c("BC_indices","BC_values")
head(BC)
ndim=3
mydim=2
GCV=1


#output_CPP = smooth.FEM.basis(locations  = as.matrix(locations), 
#                              observations = data, 
#                             FEMbasis = FEMbasis, lambda = lambda, 
#                             covariates = covariates,ndim=ndim,mydim=mydim, 
#                             GCV = TRUE,
#                              CPP_CODE = TRUE)


  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = ndim)
  }
  
  if(is.null(BC$BC_indices))
  {
    BC$BC_indices<-vector(length=0)
  }else
  {
    BC$BC_indices<-as.vector(BC$BC_indices)-1
  }
  
  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }
  

  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  data <- as.vector(data)
  storage.mode(data) <- "double"
  storage.mode(order) <- "integer"
  storage.mode(mesh$nnodes) <- "integer"
  storage.mode(mesh$ntriangles) <- "integer"
  storage.mode(mesh$nodes) <- "double"
  storage.mode(mesh$triangles) <- "integer"
  covariates = as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(lambda) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values)  <- "double"
  
  GCV = as.integer(GCV)
  storage.mode(GCV)<-"integer"
  
  ## Call C++ function
  bigsol <- .Call("regression_Laplace", locations, data, mesh, 
                  order,mydim,ndim, lambda, covariates,
                  BC$BC_indices, BC$BC_values, GCV,
                  package = "fdaPDE")
  
write.table(bigsol[[1]],file="res_caramella_new.txt")


