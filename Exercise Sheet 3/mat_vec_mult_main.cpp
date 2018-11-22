#include <iostream>
#include "VectorClass.h"
#include "CsrMatrixClass.h"

int main() {
	std::cout << "MAtrix Vector multiplication\n";
	double valsm[19] = { -100,4,2,1,-50,4,2,-110,3,6,-20,2,1,2,3,-90,3,7,-90 };
	int row[19] = { 0,0,0,1,1,1,2,2,2,2,3,3,4,4,4,4,4,5,5 };
	int c[19] = { 0,1,3,0,1,3,1,2,3,4,3,5,1,2,3,4,5,2,5 };
	CsrMatrix A(6,6,19);
	A.csr_assemble(valsm, row, c, 19);
	A.print_csr_matrix();

	Vector v(6),r(9);
	double valsv[6] = { 1,1,0,0,0,0 };
	v.vec_assemble(valsv);
	v.print_vector();
	r = A * v;
	std::cout << "A*v:\n";
	r.print_vector();
	std::cout << "\nFinished\n";
	return 0;
}