plot.surface.mesh<-function(mesh,node_values=NULL){
  
  if(!require(rgl)){
    stop("The plot surface_mesh_function(...) requires the R package rgl, please install it and try again!")
  }
  
  #p <- colorRampPalette(c("#0E1E44","#3E6DD8","#68D061","#ECAF53", "#EB5F5F","#E11F1C"))(128)
  #palette(p)
  
  order=mesh$order
  nnodes=mesh$nnodes
  ntriangles=mesh$ntriangles
  
  if(is.null(node_values)){
    rgl.open()
    triangle = c(mesh$triangles[1]-1,mesh$triangles[2]-1,mesh$triangles[3]-1)
    vertices = as.numeric(c(
      mesh$nodes[3*order*triangle[1]+1],mesh$nodes[3*order*triangle[1]+2],mesh$nodes[3*order*triangle[1]+3],1,
      mesh$nodes[3*order*triangle[2]+1],mesh$nodes[3*order*triangle[2]+2],mesh$nodes[3*order*triangle[2]+3],1,
      mesh$nodes[3*order*triangle[3]+1],mesh$nodes[3*order*triangle[3]+2],mesh$nodes[3*order*triangle[3]+3],1))
    
    bg3d(color = "white")
    indices=c(1,2,3)
    wire3d(tmesh3d(vertices,indices) , col="black")
    
    for(i in 2:ntriangles){
      triangle = c(mesh$triangles[3*order*(i-1)+1]-1,mesh$triangles[3*order*(i-1)+2]-1,mesh$triangles[3*order*(i-1)+3]-1)
      vertices = as.numeric(c(
        mesh$nodes[3*order*triangle[1]+1],mesh$nodes[3*order*triangle[1]+2],mesh$nodes[3*order*triangle[1]+3],1,
        mesh$nodes[3*order*triangle[2]+1],mesh$nodes[3*order*triangle[2]+2],mesh$nodes[3*order*triangle[2]+3],1,
        mesh$nodes[3*order*triangle[3]+1],mesh$nodes[3*order*triangle[3]+2],mesh$nodes[3*order*triangle[3]+3],1))
      
      indices=c(1,2,3)
      wire3d(tmesh3d(vertices,indices) , col="black")
    }
  }else{
    diffrange = max(node_values)-min(node_values)
    rgl.open()
    triangle = c(mesh$triangles[1]-1,mesh$triangles[2]-1,mesh$triangles[3]-1)
    vertices = as.numeric(c(
      mesh$nodes[3*order*triangle[1]+1],mesh$nodes[3*order*triangle[1]+2],mesh$nodes[3*order*triangle[1]+3],1,
      mesh$nodes[3*order*triangle[2]+1],mesh$nodes[3*order*triangle[2]+2],mesh$nodes[3*order*triangle[2]+3],1,
      mesh$nodes[3*order*triangle[3]+1],mesh$nodes[3*order*triangle[3]+2],mesh$nodes[3*order*triangle[3]+3],1))
    indices=c(1,2,3)
    col = mean(node_values[triangle[1]+1],node_values[triangle[2]+1],node_values[triangle[3]+1])
    col= (col - min(node_values))/diffrange*127+1
    shade3d( tmesh3d(vertices,indices) , col=col)
    bg3d(color = "white")
    
    for(i in 2:ntriangles){
      triangle = c(mesh$triangles[3*order*(i-1)+1]-1,mesh$triangles[3*order*(i-1)+2]-1,mesh$triangles[3*order*(i-1)+3]-1)
      vertices = as.numeric(c(
        mesh$nodes[3*order*triangle[1]+1],mesh$nodes[3*order*triangle[1]+2],mesh$nodes[3*order*triangle[1]+3],1,
        mesh$nodes[3*order*triangle[2]+1],mesh$nodes[3*order*triangle[2]+2],mesh$nodes[3*order*triangle[2]+3],1,
        mesh$nodes[3*order*triangle[3]+1],mesh$nodes[3*order*triangle[3]+2],mesh$nodes[3*order*triangle[3]+3],1))
      
      indices=c(1,2,3)
      col = mean(node_values[triangle[1]+1],node_values[triangle[2]+1],node_values[triangle[3]+1])
      col= (col - min(node_values))/diffrange*127+1
      shade3d( tmesh3d(vertices,indices) , col= col)
    }
  }
  
}
