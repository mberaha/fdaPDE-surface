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

V = read.table(file="hubV.csv",header=F,sep=",")
T = read.table(file="hubT.csv",header=F,sep=",")

#locations=list(nodes=V,triangles=T,nnodes=nrow(V),ntriangles=nrow(T))
mymesh2<-second.order.mesh(V,T)

times=50
MSE_hub=numeric(times)

data_nonoise = read.table(file="hub_true_Matlab.csv",sep=",",header=F)
noise = read.table(file="hub_onlyNoise_Matlab.csv",sep=",",header = F)
observed = data_nonoise + noise

mesh <- create.surface.mesh(mymesh2$nodes, mymesh2$triangles, order=2)
FEMbasis <- create.FEM.basis(mesh)
lambda=c(0.00375)

for(j in 1:times){
  ########## Generating Test data ############

  func_evaluation = data_nonoise[,j]
  
  data=observed[,j]
  
  output_CPP =smooth.FEM.basis(locations=V, observations = data, 
                               FEMbasis = FEMbasis, lambda = lambda,
                               CPP_CODE = TRUE) 
  
  MSE_hub[j] = sum((as.vector(output_CPP$fit.FEM$coeff)-func_evaluation)^2)/nnodes
}

boxplot(MSE_hub)

write.table(MSE_hub, file = "MSE_hub_R_v2_order2.csv",sep=",",row.names=F)



#### CAROTIDE ####


#V = read.table(file="ICA_V.csv",header=F,sep=",")
#T = read.table(file="ICA_T.csv",header=F,sep=",")
#mymesh=list(nodes=V,triangles=T,nnodes=nrow(V),ntriangles=nrow(T))


#times=25
#MSE_ICA=numeric(times)

#data_nonoise = read.table(file="ICA_true_Matlab.csv",sep=",",header=F)
#noise = read.table(file="ICA_onlyNoise_Matlab.csv",sep=",",header = F)
#observed = data_nonoise + noise

#mesh <- create.surface.mesh(mymesh$nodes, mymesh$triangles, order=1)
#FEMbasis <- create.FEM.basis(mesh)
#lambda=c(2.5*10^(-3))
#nnodes=mesh$nnodes

for(j in 1:times){

 # func_evaluation = data_nonoise[,j]
  
 # data=observed[,j]
  
 # output_CPP =smooth.FEM.basis(observations = data, 
 #                              FEMbasis = FEMbasis, lambda = lambda,
 #                              CPP_CODE = TRUE) 
  
 # MSE_ICA[j] = sum((as.vector(output_CPP$fit.FEM$coeff)-func_evaluation)^2)/nnodes
}

#boxplot(MSE_ICA)

#write.table(MSE_ICA, file = "MSE_ICA_R_v3.csv",sep=",",row.names=F)
