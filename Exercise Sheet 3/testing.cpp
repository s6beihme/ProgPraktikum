#include <iostream>
#include "VectorClass.h"
#include "CsrMatrixClass.h"

//teste ob a=b, dann b verändern a verädert
int main() {
	//TESTING ASSIGNING A MATRIX
	/*double valsm[19] = { -100,4,2,1,-50,4,2,-110,3,6,-20,2,1,2,3,-90,3,7,-90 };
	int row[19] = { 0,0,0,1,1,1,2,2,2,2,3,3,4,4,4,4,4,5,5 };
	int c[19] = { 0,1,3,0,1,3,1,2,3,4,3,5,1,2,3,4,5,2,5 };
	CsrMatrix A(6, 6, 19);
	A.csr_assemble(valsm, row, c, 19);*/
	/*CsrMatrix Minv= A.inv_diagonal();
	//Minv = (A.inv_diagonal()); //this doesnt work!!!!!! (ASSIGNING MINV TO A.INV_DIAGONAL DIRECTLY DOES)
	//Minv = A;
	std::cout << "\nafter assigning Minv";
	Minv.print_csr_matrix();*/

	//TESTING ASSIGNING A VECTOR
	/*double valsm[4] = { 1,1,1,1 };
	int row[4] = { 0,0,1,1 };
	int col[4] = { 0,1,0,1 };
	CsrMatrix A(2, 2, 4);
	A.csr_assemble(valsm, row, col, 4);
	A.print_csr_matrix();*/

	/*Vector v(2),r(2);
	double valsr[6] = { 1,1 };
	double valsv[6] = { -96,-49 };
	v.vec_assemble(valsv);
	r.vec_assemble(valsr);
	//v.print_vector();
	//r.print_vector();
	Vector v2(0);
	v2 = A*r; //DOESNT WORK!!!
	(A*r).print_vector(); // prints 2,2
	v2.print_vector(); //prints 2,0
	//v.print_vector();
	//r.print_vector();
	return 0;*/


	//TESTING Richardson
	/*std::cout << "\nSolving a small linear system of Equations using gs_solve\n";
	double valsm[19] = { -100,4,2,1,-50,4,2,-110,3,6,-20,2,1,2,3,-90,3,7,-90 };
	int row[19] = { 0,0,0,1,1,1,2,2,2,2,3,3,4,4,4,4,4,5,5 };
	int c[19] = { 0,1,3,0,1,3,1,2,3,4,3,5,1,2,3,4,5,2,5 };
	CsrMatrix A(6, 6, 19);
	A.csr_assemble(valsm, row, c, 19);
	std::cout << "A is:\n";
	A.print_csr_matrix();

	Vector x(6), b(6);
	double valsx[6] = { 0,0,0,0,0,0 };
	double valsb[6] = { -96,-49,2,0,1,0 };
	x.vec_assemble(valsx);
	b.vec_assemble(valsb);
	A.Richardson_prec(b, x, 0);
	std::cout << "\nbefore printing Result\n";
	x.print_vector();*/

	//TEST IF (A-B):NORM_SQUARED CHANGES A OR B
	/*Vector v(2), r(2);
	double valsr[6] = { 1,1 };
	double valsv[6] = { 5,5 };
	v.vec_assemble(valsv);
	r.vec_assemble(valsr);
	v.print_vector();
	r.print_vector();
	std::cout << (v - r).norm_squared() << std::endl; 
	v.print_vector();
	r.print_vector();*/
	
	//TESTING SOLVE ONLY DIAG
	/*double valsm[4] = { 2,1,1,2 };
	int row[4] = { 0,0,1,1 };
	int col[4] = { 0,1,0,1 };
	CsrMatrix A(2, 2, 4);
	A.csr_assemble(valsm, row, col, 4);
	A.print_csr_matrix();

	Vector b(2), x(2);
	double valsb[6] = { 1,1 };
	b.vec_assemble(valsb);
	std::cout << "\nResult:\n";
	A.solve_only_diag(b,x);
	x.print_vector();*/

	//TESTING ONLY UPPER TRIANG
	/*std::cout << "\nSolving a small linear system of Equations using gs_solve\n";
	double valsm[19] = { -7,4,2,1,-5,4,2,-11,3,6,-2,2,1,2,3,-9,3,7,-9 };
	int row[19] = { 0,0,0,1,1,1,2,2,2,2,3,3,4,4,4,4,4,5,5 };
	int c[19] = { 0,1,3,0,1,3,1,2,3,4,3,5,1,2,3,4,5,2,5 };
	CsrMatrix A(6, 6, 19);
	A.csr_assemble(valsm, row, c, 19);
	std::cout << "A is:\n";
	A.print_csr_matrix();

	Vector x(6), b(6);
	double valsx[6] = { 0,0,0,0,0,0 };
	double valsb[6] = { -1,-1,-2,0,-6,-9 };
	x.vec_assemble(valsx);
	b.vec_assemble(valsb);
	A.solve_only_upper_triang(b,x);
	x.print_vector();*/


	CsrMatrix A("C:/Users/ASUS/Documents/Studium/Prog Praktikum/Exercise Sheet 3/matfile.txt");
	Vector v("C:/Users/ASUS/Documents/Studium/Prog Praktikum/Exercise Sheet 3/vecfile.txt");
	
	Vector x4(v.get_size());
	A.CG_prec(v, x4, 0);
	
//std::fstream file;
//file.open("C:/Users/ASUS/Documents/Studium/Prog Praktikum/Exercise Sheet 3/matfile.txt");
//if (!file.is_open()) std::cout << "failed\n";
}