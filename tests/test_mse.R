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

#### HUB ######

filename = 'hub.csv'

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

times=50
MSE_hub=numeric(times)

for(j in 1:times){
  ########## Generating Test data ############
  
  a1 = rnorm(1,mean = 1, sd = 1)
  a2 = rnorm(1,mean = 1, sd = 1)
  a3 = rnorm(1,mean = 1, sd = 1)
  
  nnodes=mymesh$nnodes
  func_evaluation = numeric(nnodes)
  
  for (i in 1:nnodes){
    func_evaluation[i] = a1* sin(2*pi*mymesh$nodes[i,1]) +  a2* sin(2*pi*mymesh$nodes[i,2]) +  a3*sin(2*pi*mymesh$nodes[i,3]) +1
  }
  
  data=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)
  
  
  mesh <- create.surface.mesh(mymesh$nodes, mymesh$triangles, order=1)
  
  FEMbasis <- create.FEM.basis(mesh)
  
  lambda=c(0.00375)
  output_CPP =smooth.FEM.basis(observations = data, 
                               FEMbasis = FEMbasis, lambda = lambda,
                               CPP_CODE = TRUE) 
  
  MSE_hub[j] = sum((as.vector(output_CPP$fit.FEM$coeff)-func_evaluation)^2)/nnodes
}

boxplot(MSE_hub)

####################### Caramella ################

filename = 'Caramella.csv'


mymesh=read.mesh(filename)

times=50
MSE_caramella=numeric(times)

for(j in 1:times){
  ########## Generating Test data ############
  
  a1 = rnorm(1,mean = 1, sd = 1)
  a2 = rnorm(1,mean = 1, sd = 1)
  a3 = rnorm(1,mean = 1, sd = 1)
  
  nnodes=mymesh$nnodes
  func_evaluation = numeric(nnodes)
  
  for (i in 1:nnodes){
    func_evaluation[i] = a1* sin(2*pi*mymesh$nodes[i,1]) +  a2* sin(2*pi*mymesh$nodes[i,2]) +  a3*sin(2*pi*mymesh$nodes[i,3]) +1
  }
  
  data=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)
  
  
  mesh <- create.surface.mesh(mymesh$nodes, mymesh$triangles, order=1)
  
  FEMbasis <- create.FEM.basis(mesh)
  
  lambda=c(0.005)
  output_CPP =smooth.FEM.basis(observations = data, 
                               FEMbasis = FEMbasis, lambda = lambda,
                               CPP_CODE = TRUE) 
  
  MSE_caramella[j] = sum((as.vector(output_CPP$fit.FEM$coeff)-func_evaluation)^2)/nnodes
}

boxplot(MSE_caramella)

MSE=cbind(MSE_hub,MSE_caramella)
MSE=as.data.frame(MSE)
colnames(MSE) = c("hub","caramella")
write.table(MSE,file="MSE_hub_caramella.csv",sep=",")



########### Carotide ##############

filename = '../paziente12348/geom_211.csv'


mymesh=read.mesh(filename)

times=25
MSE_carotide=numeric(times)

for(j in 1:times){  
  a1 = rnorm(1,mean = 1, sd = 1)
  a2 = rnorm(1,mean = 1, sd = 1)
  a3 = rnorm(1,mean = 1, sd = 1)
  
  nnodes=mymesh$nnodes
  func_evaluation = numeric(nnodes)
  
  for (i in 1:nnodes){
    func_evaluation[i] = a1* sin(2*pi*mymesh$nodes[i,1]) +  a2* sin(2*pi*mymesh$nodes[i,2]) +  a3*sin(2*pi*mymesh$nodes[i,3]) +1
  }
  
  data=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)
  
  
  mesh <- create.surface.mesh(mymesh$nodes, mymesh$triangles, order=1)
  
  FEMbasis <- create.FEM.basis(mesh)
  
  lambda=c(0.001,0.00375,0.005,0.01)
  output_CPP =smooth.FEM.basis(observations = data, 
                               FEMbasis = FEMbasis, lambda = lambda,
                               CPP_CODE = TRUE,GCV=TRUE) 
  
  MSE_carotide[j] = sum((as.vector(output_CPP$fit.FEM$coeff)-func_evaluation)^2)/nnodes
}

write.table(MSE_carotide,file="MSE_carotide.csv",sep=",")

########### Carotide No Aneurism ##########

filename = '108100_noAne_noHole.csv'

mymesh=read.mesh(filename)

times=25
MSE_carotide=numeric(times)

for(j in 1:times){  
  a1 = rnorm(1,mean = 1, sd = 1)
  a2 = rnorm(1,mean = 1, sd = 1)
  a3 = rnorm(1,mean = 1, sd = 1)
  
  nnodes=mymesh$nnodes
  func_evaluation = numeric(nnodes)
  
  for (i in 1:nnodes){
    func_evaluation[i] = a1* sin(2*pi*mymesh$nodes[i,1]) +  a2* sin(2*pi*mymesh$nodes[i,2]) +  a3*sin(2*pi*mymesh$nodes[i,3]) +1
  }
  
  data=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)
  
  
  mesh <- create.surface.mesh(mymesh$nodes, mymesh$triangles, order=1)
  
  FEMbasis <- create.FEM.basis(mesh)
  
  lambda=c(0.001,0.00375,0.005,0.01)
  output_CPP =smooth.FEM.basis(observations = data, 
                               FEMbasis = FEMbasis, lambda = lambda,
                               CPP_CODE = TRUE,GCV=TRUE) 
  
  MSE_carotide[j] = sum((as.vector(output_CPP$fit.FEM$coeff)-func_evaluation)^2)/nnodes
}

write.table(MSE_carotide,file="MSE_carotide2.csv",sep=",")


