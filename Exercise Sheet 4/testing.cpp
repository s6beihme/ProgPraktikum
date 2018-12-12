#include "VectorClass.h"
//#include "CsrMatrixClass.h"
#include "SplineClass.h"
#include <complex>
int main() {
	Spline<double, 3> a;


	//double knots[4] = { 0,0.3,0.6,1 };
	std::unique_ptr<std::unique_ptr<double[]>[]> ctr = std::make_unique<std::unique_ptr<double[]>[]>(5);

	for (int i = 0; i < 5; i++) ctr[i] = std::make_unique<double[]>(3);
	ctr[0][0] = 1.11;
	ctr[0][1] = 0;
	ctr[0][2] = 3;
	ctr[1][0] = 6.5;
	ctr[1][1] = 3;
	ctr[1][2] = 5;
	ctr[2][0] = 2;
	ctr[2][1] = 4.9;
	ctr[2][2] = 2;
	ctr[3][0] = 1.4;
	ctr[3][1] = 5.4;
	ctr[3][2] = 3.4;
	ctr[4][0] = 0;
	ctr[4][1] = 0;
	ctr[4][2] = 0;
	
	std::unique_ptr<double[]> param = std::make_unique<double[]>(5);
	a.create_interpol_parameters(ctr, param, 5);
	a.create_knots(param,  5, 2);
	std::cout << "\n";
	a.print_knots();
	std::cout << a.find_intervall(0.6) << "\n";
	std::cout << a.eval_bspline(0.6, 2, 2) << "\n";
	
	
	/*Spline<double, 2> a(2, 5);
	double knots[4] = { 0,0.3,0.6,1 };
	double** ctr = new double*[5];
	for (int i = 0; i < 5; i++) ctr[i] = new double[2];
	ctr[0][0] = 0;
	ctr[0][1] = 0;
	ctr[1][0] = 6;
	ctr[1][1] = 3;
	ctr[2][0] = 2;
	ctr[2][1] = 2;
	ctr[3][0] = 1.4;
	ctr[3][1] = 5.4;
	ctr[4][0] = 0;
	ctr[4][1] = 0;
	a.assemble(4, 5, 2, knots, ctr);
	a.write_to_file("C:/Users/ASUS/Documents/Studium/Prog Praktikum/Exercise Sheet 4/splinefile1.txt", 50);*/
	std::cout << "\nEnd";
	return 0;
}