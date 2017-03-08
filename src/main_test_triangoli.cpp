#include"mesh_objects.h"
#include"mesh.h"
#include <iostream>
#include <vector>

int main()
{ //! Ciao Anna :P
	/*Point a(1,1,0,0,0);
	Point b(2,2,0,1,0);
	Point c(3,3,0,0,1);

	a.print(std::cout);
	b.print(std::cout);
	c.print(std::cout);

	std::cout<<b.getId()<<std::endl;
	b.print(std::cout);

	std::vector<Point> points;
	points.push_back(c);
	points.push_back(b);
	points.push_back(a);

	Triangle<3,2,3> t(1,points);

	std::cout<<"I punti sono"<<std::endl;
	for (int i=0; i<3; ++i)
		points[i].print(std::cout);
	std::cout<<std::endl;

	//t.print(std::cout);

	Point P(4,4,0,0.25,0.25);
	P.print(std::cout);

	Eigen::Matrix<Real,3,1> BaryHope=t.getBaryCoordinates(P);

	for (int i=0; i<3; ++i)
		std::cout<<BaryHope[i]<<std::endl;

	bool verita = t.isPointInside(P);
	std::cout<< "il punto è dentro iff "<<verita<<std::endl;	*/

	/////
	//! Test sulla mesh

	std::string filename("Caramella.csv");
	std::cout<<"filename = "<<filename<<"\n";
	std::cout<<"calling the mesh"<<"\n";
	MeshHandler<1,2,3> mesh(filename);
	std::cout<<"finished importing the mesh"<<"\n";
	std::cout<<"num_nodes="<<mesh.num_nodes()<<"\n";
	std::cout<<"num_triangles="<<mesh.num_triangles()<<"\n";

	mesh.printPoints(std::cout);
	mesh.printTriangles(std::cout);

return 0;
}
