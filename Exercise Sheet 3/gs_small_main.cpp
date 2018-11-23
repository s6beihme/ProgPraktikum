#include <iostream>
#include "VectorClass.h"
#include "CsrMatrixClass.h"




int main() {
	/*std::cout << "\nSolving a small linear system of Equations using gs_solve\n";
	double valsm[9] = { 5,2,2,2,5,2,2,2,5 };
	int row[9] =      { 0,0,0,1,1,1,2,2,2 };
	int c[9] =        { 0,1,2,0,1,2,0,1,2 };
	CsrMatrix A(3,3,9);
	A.csr_assemble(valsm, row, c, 9);
	std::cout << "A is:\n";
	A.print_csr_matrix();

	Vector a(3), b(3);
	double valsa[6] = { 0,0,0 };
	double valsb[6] = { 9,9,9 };
	a.vec_assemble(valsa);
	b.vec_assemble(valsb);
	A.gs_solve(a, b);
	a.print_vector();*/

	std::cout << "\nSolving a small linear system of Equations using gs_solve\n";
	double valsm[19] = { -100,4,2,1,-50,4,2,-110,3,6,-20,2,1,2,3,-90,3,7,-90 };
	int row[19] = { 0,0,0,1,1,1,2,2,2,2,3,3,4,4,4,4,4,5,5 };
	int c[19] = { 0,1,3,0,1,3,1,2,3,4,3,5,1,2,3,4,5,2,5 };
	CsrMatrix A(6, 6, 19);
	A.csr_assemble(valsm, row, c, 19);
	std::cout << "A is:\n";
	A.print_csr_matrix();
	
	Vector a(6), b(6);
	double valsa[6] = { 0,0,0,0,0,0 };
	double valsb[6] = { -96,-49,2,0,1,0 };
	a.vec_assemble(valsa);
	b.vec_assemble(valsb);
	A.gs_solve(b, a);
	a.print_vector();


	return 0;
}

