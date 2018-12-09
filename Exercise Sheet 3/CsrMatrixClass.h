#pragma once
#include "VectorClass.h"
#include <string>


//about the Solvers: gaus seidel yields reasonable results most constantly (compared to the otherss)


class CsrMatrix {
public:
	//allocates memory for CsrMatrix object
	//parameters:
		//rows: number of rows in Matrix
		//cols: -"- cols -"-
		//nnz: number of non-zero entries in matrix
	CsrMatrix(int rows, int cols, int _nnz);

	//reades Matrix from a file and creates matrix object accordingly
	//file has to be in following format:
		//line 1:	 row_count col_count nnz
		//line 2:	 row column entry
		//line 3:	 row column entry
		//...
		//line nnz+1:row column entry 
	//parameters:
		//filename: name of file in which matrix is written
	CsrMatrix(std::string filename);

	~CsrMatrix();

	//prints the entr, col_ind and new_row_ind array
	void print_csr_matrix() const;
	int get_row_count() const;
	int get_col_count() const;
	int get_nnz() const;

	//returns entr[i]
	double get_entr(int i) const;

	//fills entr, col_ind and new_row_ind according to input
	//parameters:
		//vals: array containing the values that should be stored in the matrix
		//row_inds: array containing the row indices of the values -"-
		//col_inds: array containing the column indices of the values-"-
		//nnz: number of given non-zero entries
	void csr_assemble(double* vals, int* row_inds, int*col_inds, int nnz);

	//matrix vector produkt
	//parameters:
		//v: vector with which the matrix is multiplied (has to have correct dimension)
		//result: functions fills result vector with (*this)*v
	void mat_vec_mult(const Vector &v, Vector& result) const; 
	

	//solves linear system of equations Mx=b for x algebraically where M is diagonal part of (*this)
	//parameters:
		//b: right hand side vector of equation
		//result: is filled with result such that M*result=b
	void solve_only_diag(const Vector& b, Vector& result) const;

	//solves linear system of equations Mx=b for x algebraically where M is upper triagonal part of (*this)
	//parameters:
		//b: right hand side vector of equation
		//result: is filled with result such that M*result=b
	void solve_only_upper_triang(const Vector& b, Vector& result);
	//solves linear system of equations Mx=b for x algebraically where M is the upper triangular part of (*this)
	//Vector solve_only_upper_triang(const Vector& b);

	//solves Linear system of Equations Au=b for x using gauss seidel iteration
	//parameters:
		//result: initial guess vector of result. it will be overwritten with the result
		//b: right hand side vector of Equation
	void gs_solve(const Vector& b, Vector& result) const;

	//solves linear system of equations Ax=b for x using preconditioned Conjugate Gradient Method
	//parameters:
		//b: right hand side vector of equation
		//x0: initial guess. result will be stored in here
		//preconditioner: if =0 then uses Jacobi preconditioner, if =1 uses gauss seidel preconditioner
	//returns: Vector type object x that satisfies Ax=b
	void CG_prec(const Vector& b, Vector& x0, int preconditioner);

	//solves linear system of equations Ax=b for x using preconditioned Richardson Method (gamma_n=1 for all n)
	//parameters:
		//b: right hand side vector of equation
		//x0: initial guess. result will be stored in here
		//preconditioner: if =0 then uses Jacobi preconditioner, if =1 uses gauss seidel preconditioner
	//returns: Vector type object x that satisfies Ax=b
	void Richardson_prec(const Vector& b, Vector& x0, int preconditioner);
	 
private:
	int nnz;
	int row_count;
	int col_count;
	//double* entr;
	std::unique_ptr<double[]> entr;
	//int* col_ind;
	std::unique_ptr<int[]> col_ind;
	//int* new_row_ind;
	std::unique_ptr<int[]> new_row_ind;
};






