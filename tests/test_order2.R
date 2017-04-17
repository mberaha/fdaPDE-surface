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


set.seed(13042017)

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
print('ci sono')

V = read.table(file="hubV.csv",header=F,sep=",")
T = read.table(file="hubT.csv",header=F,sep=",")
#mymesh=list(nodes=V,triangles=T,nnodes=nrow(V),ntriangles=nrow(T))
mymesh2<-second.order.mesh(V,T)
plot.surface.mesh(mymesh2)
cat('size of T: ',dim(mymesh2$triangles))
print('second order done\n')
hub_true = NULL
hub_noise = NULL

#mymesh=read.mesh(filename)

times=50
MSE_hub=numeric(times)

for(j in 1:times){
  ########## Generating Test data ############
print('ciclo for\n')  
  a1 = rnorm(1,mean = 1, sd = 1)
  a2 = rnorm(1,mean = 1, sd = 1)
  a3 = rnorm(1,mean = 1, sd = 1)
  
  nnodes=mymesh2$nnodes
	cat('numero nodi mesh ',nnodes)
  func_evaluation = numeric(nnodes)
  
  for (i in 1:nnodes){
    func_evaluation[i] = a1* sin(2*pi*mymesh2$nodes[i,1]) +  a2* sin(2*pi*mymesh2$nodes[i,2]) +  a3*sin(2*pi*mymesh2$nodes[i,3]) +1
  }
  print('func evaluuation done')
  hub_true = rbind(hub_true,func_evaluation)
  
  data=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)
  
  hub_noise = rbind(hub_noise,data)
  
  
  mesh <- create.surface.mesh(mymesh2$nodes, mymesh2$triangles, order=2)
  
  FEMbasis <- create.FEM.basis(mesh)
  
  lambda=c(0.00375)
  output_CPP =smooth.FEM.basis(observations = data, 
                               FEMbasis = FEMbasis, lambda = lambda,
                               CPP_CODE = TRUE) 
  print('cpp ended')
  MSE_hub[j] = sum((as.vector(output_CPP$fit.FEM$coeff)-func_evaluation)^2)/nnodes
}

boxplot(MSE_hub)

write.table(hub_true,file="hub_true_order2.csv",sep=",",row.names = F)
write.table(hub_noise,file="hub_noise_order2.csv",sep=",",row.names = F)
write.table(MSE_hub, file = "MSE_hub_R_order2.csv",sep=",",row.names=F)



