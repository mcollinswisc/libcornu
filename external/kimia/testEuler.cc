// Include standard C++ input/output library
#include <iostream>
#include "euler.hh"

// main is the first function to be called
int main() {
        std::cout << "Using euler spiral comletion from kimia2004..." << std::endl;

	//Point2D<double> start_pt = new Point2D<double> (1,2);
	Point2D<double> start_pt(1, 2);

	double start_angle = 0;
	//Point2D<double> end_pt = new Point2D<double> (30,40);
	Point2D<double> end_pt(40, 2);
	double end_angle = 0;
	//EulerSpiral spiral = new EulerSpiral ( start_pt, start_angle, end_pt, end_angle);
	EulerSpiral eulerspiral(start_pt, start_angle, end_pt, end_angle);

	eulerspiral.compute_es_params();

	std::vector<Point2D<double> > spiral;
	int N = 100;
	eulerspiral.computeSpiral(spiral, 0, 0);

	double eulerParams[3];
	eulerspiral.getParams(eulerParams);

	std::cout << "kappa0: " << eulerParams[0] << std::endl;
	std::cout << "gamma: " << eulerParams[1] << std::endl;
	std::cout << "L: " << eulerParams[2] << std::endl;

	for (int i = 0; i < N; i++) {
		Point2D<double> p;
		p = spiral.at(i);
		std::cout << "Point " << i << " x: " << p.getX() << "; y: " << p.getY()
				<< std::endl;
	}

	return 0;
}
