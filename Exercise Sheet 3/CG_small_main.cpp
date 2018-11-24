#include <iostream>
#include "VectorClass.h"
#include "CsrMatrixClass.h"

int main() {
	std::cout << "\nSolving a small linear system of Equations using CG_prec\n";
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
	A.CG_prec(b, a, 0);
	a.print_vector();
}