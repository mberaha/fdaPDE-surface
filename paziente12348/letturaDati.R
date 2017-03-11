
library(R.matlab)

write_mesh_file <- function(filename=NULL,nodes,triangles){
  if(is.null(filename))
    filename='mymesh.csv'
  nnodes = dim(nodes)[1]
  ntriangles = dim(triangles)[1]
  cat("num_nodes",nnodes,"\nnum_triangles",ntriangles,"\n\n", file=filename)
  write.table(nodes,file=filename,sep=",",row.names = F, col.names = F,append = TRUE)
  cat("\n",file=filename, append = T)
  write.table(triangles,file=filename,sep=",",row.names = F, col.names = F,append = TRUE)
}
#setwd("C:\Users\Utente\Desktop\AnuRisk\12348\wss_vtk")
for(i in 211:220){
file =paste0("wss_",toString(i),".mat")
filemat<-readMat(file)
filename = paste0("wss_",toString(i),".csv")
write_mesh_file(filename=filename,nodes=filemat$V,triangles=filemat$T)
}













read.mesh<-function(filename){
  nnodes = read.table(filename, nrows = 2, header = FALSE, sep =' ', stringsAsFactors = FALSE)[1,2]
  ntriangles = read.table(filename, nrows = 2, header = FALSE, sep =' ', stringsAsFactors = FALSE)[2,2]
  geometry = read.csv(filename, skip=2, header= FALSE)
  nodes = geometry[1:nnodes,]
  triangles = geometry[(nnodes+1):(nnodes+ntriangles),]
  
  retlist = list(nnodes=nnodes,ntriangles=ntriangles,nodes=nodes,triangles=triangles)
  return(retlist)
}

mesh=read.mesh(filename)


#### CREATING A MESH FILE FROM NODES AND CONNECTIVITY #####

