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

read.mesh<-function(filename){
  nnodes = read.table(filename, nrows = 2, header = FALSE, sep =' ', stringsAsFactors = FALSE)[1,2]
  ntriangles = read.table(filename, nrows = 2, header = FALSE, sep =' ', stringsAsFactors = FALSE)[2,2]
  geometry = read.csv(filename, skip=2, header= FALSE)
  nodes = geometry[1:nnodes,]
  triangles = geometry[(nnodes+1):(nnodes+ntriangles),]
  
  retlist = list(nnodes = nnodes,ntriangles = ntriangles,nodes = nodes,triangles = triangles)
  return(retlist)
}


filename = '108100_noAne_noHole.csv'

mymesh=read.mesh(filename)

### Only for curve cyl ####
V = read.table(file="V_curveCyl.csv",header=F,sep=",")
T = read.table(file="T_curveCyl.csv",header=F,sep=",")
mymesh=list(nodes=V,triangles=T,nnodes=nrow(V),ntriangles=nrow(T))



### GENERATING DATA #####

set.seed(6042017)

a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)

nnodes=mymesh$nnodes
func_evaluation = numeric(nnodes)

for (i in 1:nnodes){
  func_evaluation[i] = a1* sin(2*pi*mymesh$nodes[i,1]) +  a2* sin(2*pi*mymesh$nodes[i,2]) +  a3*sin(2*pi*mymesh$nodes[i,3]) +1
}

observations=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)
lambda=c(1)
locations = NULL
covariates = NULL
BC = NULL
ndim = 3
mydim = 2
GCV=FALSE

mesh=create.surface.mesh(mymesh$nodes,mymesh$triangles)

if(!is.null(locations))
  locations = as.matrix(locations)
observations = as.matrix(observations)
lambda = as.matrix(lambda)
if(!is.null(covariates))
  covariates = as.matrix(covariates)
if(!is.null(BC))
{
  BC$BC_indices = as.matrix(BC$BC_indices)
  BC$BC_values = as.matrix(BC$BC_values)
}

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
  BC$BC_indices<-as.vector(BC$BC_indices)
  
}

if(is.null(BC$BC_values))
{
  BC$BC_values<-vector(length=0)
}else
{
  BC$BC_values<-as.vector(BC$BC_values)
}

## Set propr type for correct C++ reading
locations <- as.matrix(locations)
storage.mode(locations) <- "double"
data <- as.vector(observations)
storage.mode(observations) <- "double"
storage.mode(mesh$order) <- "integer"
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
library(microbenchmark)
res_cyl= microbenchmark(out = .Call("regression_Laplace", locations, data, mesh, 
                mesh$order, mydim, ndim, lambda, covariates,
                BC$BC_indices, BC$BC_values, GCV,
                package = "fdaPDE"),times=15)



                