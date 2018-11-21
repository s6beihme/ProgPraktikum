#pragma once
#include "VectorClass.h"
class CsrMatrix {
public:
	CsrMatrix(int rows, int cols, int _nnz);
	~CsrMatrix();
	void print_csr_matrix();
	int get_row_count();
	int get_col_count();
	int get_nnz();
	double get_entr(int i);
	void csr_assemble(double* vals, int* row_inds, int*col_inds, int nnz);
	void mat_vec_multiply(Vector &v, Vector& res);
	void gs_solve(Vector& b, Vector& u);
private:
	int nnz;
	int row_count;
	int col_count;
	double* entr;
	int* col_ind;
	int* new_row_ind;
};






