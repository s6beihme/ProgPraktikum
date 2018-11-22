#pragma once
#include "VectorClass.h"
class CsrMatrix {
public:
	CsrMatrix(int rows, int cols, int _nnz);
	~CsrMatrix();
	void print_csr_matrix() const;
	int get_row_count() const;
	int get_col_count() const;
	int get_nnz() const;
	double get_entr(int i) const;
	void csr_assemble(double* vals, int* row_inds, int*col_inds, int nnz);
	Vector operator *(const Vector &v) const; //FRAGEN: ALS OPERATOR ÜBERLADEN UND VECTOR RETURNEN?
	void gs_solve(Vector& u, const Vector& b) const;
	CsrMatrix inv_diagonal() const;
	Vector CG_Jac_prec(const Vector b, Vector x0);
private:
	int nnz;
	int row_count;
	int col_count;
	double* entr;
	int* col_ind;
	int* new_row_ind;
};






