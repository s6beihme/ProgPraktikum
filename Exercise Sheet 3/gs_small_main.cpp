#include <iostream>
#include "VectorClass.h"
#include "CsrMatrixClass.h"




int main() {
	std::cout << "\nSolving a small linear system of Equations using gs_solve\n";
	double valsm[9] = { 5,2,2,2,5,2,2,2,5 };
	int row[9] =      { 0,0,0,1,1,1,2,2,2 };
	int c[9] =        { 0,1,2,0,1,2,0,1,2 };
	CsrMatrix A(3,3,9);
	A.csr_assemble(valsm, row, c, 9);

	Vector a(3), b(3);
	double valsa[6] = { 0,0,0 };
	double valsb[6] = { 9,9,9 };
	a.vec_assemble(valsa);
	b.vec_assemble(valsb);
	A.gs_solve(a, b);
	a.print_vector();




	return 0;
}

