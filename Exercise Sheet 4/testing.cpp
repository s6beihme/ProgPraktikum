#include <iostream>
#include <vector>
#include "SplineClass.h"
int main() {
	std::vector<std::vector<double>> intpts(5);
	std::vector<double> parameters(5);
	for (int i = 0; i < 5; i++) intpts[i] = std::vector<double>(2);
	intpts[0][0] = 1;
	intpts[0][1] = 1;
	intpts[1][0] = 6;
	intpts[1][1] = 3;
	intpts[2][0] = 2;
	intpts[2][1] = 2;
	intpts[3][0] = 1.4;
	intpts[3][1] = 5.4;
	intpts[4][0] = 0;
	intpts[4][1] = 0;
	
	Spline<double, 2> s;
	s.create_ctrpts_from_inter_pts(intpts, 5, 1);
	s.write_to_file("splinefile1.txt", 50);
	
	/*std::cout << "\nSolving a small linear system of Equations using gs_solve\n";
	double valsm[19] = { -100,4,2,1,-50,4,2,-110,3,6,-20,2,1,2,3,-90,3,7,-90 };
	int row[19] = { 0,0,0,1,1,1,2,2,2,2,3,3,4,4,4,4,4,5,5 };
	int c[19] = { 0,1,3,0,1,3,1,2,3,4,3,5,1,2,3,4,5,2,5 };
	CsrMatrix<double> A(19, 6, 6);
	A.assemble(valsm, row, c, 19);
	std::cout << "A is:\n";
	std::cout << A;

	Vector<double> a(6), b(6);
	double valsa[6] = { 0,0,0,0,0,0 };
	double valsb[6] = { -96,-49,2,0,1,0 };
	a.assemble(valsa);
	b.assemble(valsb);
	A.gs_solve(b, a);
	std::cout << a;*/
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