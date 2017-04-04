CPP_smooth.manifold.FEM.basis<-function(locations, observations, mesh, lambda, covariates = NULL, ndim, mydim, BC = NULL, GCV)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  # This is done in C++ now to optimize speed
    
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
  { cat("here BC \n")
    BC$BC_indices<-as.vector(BC$BC_indices)
    cat("here BC 2\n")
  }
  
  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  { cat("here BC 3\n")
    BC$BC_values<-as.vector(BC$BC_values)
    cat("here BC 4\n")
  }
  
  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  cat("here locations \n")
  storage.mode(locations) <- "double"
  data <- as.vector(observations)
  cat("here data \n")
  storage.mode(observations) <- "double"
  cat("here order \n")
  storage.mode(mesh$order) <- "integer"
  cat("here nnodes \n")
  storage.mode(mesh$nnodes) <- "integer"
  cat("here ntriangles \n")
  storage.mode(mesh$ntriangles) <- "integer"
  cat("here nodes \n")
  storage.mode(mesh$nodes) <- "double"
  cat("here triangles \n")
  storage.mode(mesh$triangles) <- "integer"
  covariates = as.matrix(covariates)
  storage.mode(covariates) <- "double"
  storage.mode(lambda) <- "double"
  cat("here ndim \n")
  storage.mode(ndim) <- "integer"
  cat("here mydim \n")
  storage.mode(mydim) <- "integer"
  cat("here bc_indices \n")
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values)  <- "double"
  cat("here GCV1 \n")
  GCV = as.integer(GCV)
  cat("here GCV2 \n")
  storage.mode(GCV)<-"integer"
  cat("here GCV3 \n")
  
  ## Call C++ function
  bigsol <- .Call("regression_Laplace", locations, data, mesh, 
                  mesh$order, mydim, ndim, lambda, covariates,
                  BC$BC_indices, BC$BC_values, GCV,
                  package = "fdaPDE")
  cat("smoothing CPP done \n")
  return(bigsol)
}

CPP_eval.manifold.FEM = function(FEM, locations, redundancy, ndim, mydim)
{
  FEMbasis = FEM$FEMbasis

  # Imposing types, this is necessary for correct reading from C++
  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$points) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  coeff = as.matrix(FEM$coeff)
  storage.mode(coeff) <- "double"
  storage.mode(ndim)<- "integer"
  storage.mode(mydim)<- "integer"
  storage.mode(locations) <- "double"
  storage.mode(redundancy) <- "integer"
  
  #Calling the C++ function "eval_FEM_fd" in RPDE_interface.cpp
  evalmat = matrix(0,nrow(locations),ncol(coeff))
  cat("entering for loop \n")
    for (i in 1:ncol(coeff)){
      cat(".Call \n")
      evalmat[,i] <- .Call("eval_FEM_fd", FEMbasis$mesh, locations[,1], locations[,2], locations[,3], coeff[,i], FEMbasis$order, redundancy, mydim, ndim,
                         package = "fdaPDE")
    }
  
  
  #Returning the evaluation matrix
  evalmat
}

