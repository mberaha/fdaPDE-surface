


#include "fdaPDE.h"
#include "regressionData.h"
#include "mesh_objects.h"
#include "mesh.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mixedFERegression.h"


int main{
	
	std::string filename("Caramella.csv");
	MeshHandler<1,2,3> mesh(filename);
	
	VectorXr observations;
	std::vector<Point> locations;
	UInt order=1;				
	std::vector<Real> lambda;
	lambda.push_back(1);
	MatrixXr covariates;
	std::vector<UInt> bc_indices;
	std::vector<Real> bc_values;
	
	// Mario legge tanti file :)

	RegressionData data(locations, observations, order, lambda, covariates, bc_indices, bc_values, 0);

	MixedFERegression< RegressionData, IntegratorTriangleP2, 1, 2, 3> regression(mesh, data);

	regression.smoothLaplace();
	
	const std::vector<VectorXr>& solution = regression.getSolution();

	


}
