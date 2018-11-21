#include <iostream>
#include "VectorClass.h"
#include "CsrMatrixClass.h"

int main() {
	std::cout << "\nmultiplying two vectors\n";
	Vector a(3), b(3);
	double valsa[6] = { 1,2,2.5 };
	double valsb[6] = { 2,2,4 };
	
	a.vec_assemble(valsa);
	b.vec_assemble(valsb);
	a.print_vector();
	b.print_vector();

	std::cout << "a*b = " << a * b << std::endl;
	return 0;
}