#include <iostream>
#include "CsrMatrixClass.h"




CsrMatrix::CsrMatrix(int rows, int cols, int _nnz) {
	nnz = _nnz;
	row_count = rows;
	col_count = cols;

	entr = new double[nnz];
	if (entr == NULL) {
		std::cout << "failed to create CsrMAtrix object: failed to allocate memory for entr\n";
		return;
	}

	col_ind = new int[nnz];
	if (col_ind == NULL) {
		std::cout << "failed to create CsrMAtrix object: failed to allocate memory for entr\n";
		delete[] entr;
		return;
	}

	new_row_ind = new int[row_count + 1];
	col_ind = new int[nnz];
	if (new_row_ind == NULL) {
		std::cout << "failed to create CsrMAtrix object: failed to allocate memory for entr\n";
		delete[] entr;
		delete[] col_ind;
		return;
	}
}

CsrMatrix::~CsrMatrix() {
	delete[] entr;
	delete[] col_ind;
	delete[] new_row_ind;
}


int CsrMatrix::get_row_count() const {
	return row_count;
}
int CsrMatrix::get_col_count() const {
	return col_count;
}
int CsrMatrix::get_nnz() const {
	return nnz;
}
double CsrMatrix::get_entr(int i) const {
	return entr[i];
}


void CsrMatrix::print_csr_matrix() const {
	std::cout << "\nThe entries of the Matrix are:\n";
	for (int i = 0; i < nnz; i++) {
		std::cout << "\n" << entr[i] << ", in column:" << col_ind[i];
	}
	std::cout << "\n\nThe new_row_ind array is:\n";
	for (int i = 0; i < row_count + 1; i++) {
		std::cout << new_row_ind[i] << ", ";
	}
	std::cout << std::endl;
}

void CsrMatrix::csr_assemble(double* vals, int* row_inds, int*col_inds, int _nnz) {
	if (_nnz != nnz) {
		std::cout << "csr_assemble failed because given nnz didnt match nnz of matrix\n";
		return;
	}

	//copy vals and col_inds
	for (int i = 0; i < nnz; i++) {
		entr[i] = vals[i];
		col_ind[i] = col_inds[i];
	}

	new_row_ind[0] = 0;

	int cur_pos = 0;

	//fill up new_row_inds untill you reach first no-zero entry
	for (int i = 0; i < row_inds[0]; i++) {
		cur_pos++;
		new_row_ind[cur_pos] = 0;
	}

	for (int i = 1; i < nnz; i++) {
		if (row_inds[i] > row_inds[i - 1]) {
			cur_pos++;
			new_row_ind[cur_pos] = i;
		}
	}
	cur_pos++;
	new_row_ind[cur_pos] = nnz;
}


Vector CsrMatrix::operator *(const Vector &v) const {
	if (col_count != v.get_size()) {
		std::cout << "\nMatrix Vector multiplication failed bacause dimension of first argument didnt correspond to matrix" << std::endl;
		exit(0);
	}

	double* res_vals = new double[row_count];
	double sum = 0;
	for (int i = 0; i < row_count; i++) {
		sum = 0;
		for (int j = new_row_ind[i]; j < new_row_ind[i + 1]; j++) {
			sum += entr[j] * v.get_data(col_ind[j]);
		}
		res_vals[i] = sum;
	}
	Vector res(row_count);
	res.vec_assemble(res_vals);
	delete[] res_vals;
	return res;
}

void CsrMatrix::gs_solve(Vector& u,const Vector& b) const {
	Vector Au(b.get_size());

	double norm_res = 0;

	//variable to store diagonal entry
	double aii = 0;
	int iter = 0;

	//variable to sum up non diagonal*entries of u or u2
	double sum = 0;

	while (iter < 100) {
		iter++;

		for (int i = 0; i < row_count; i++) {
			sum = 0;
			aii = 0;
			for (int j = new_row_ind[i]; j < new_row_ind[i + 1]; j++) {
				if (col_ind[j] < i) {
					sum += entr[j] * u.get_data(col_ind[j]);
				}
				if (col_ind[j] == i) {
					aii = entr[j];
				}
				if (col_ind[j] > i) {
					sum += entr[j] * u.get_data(col_ind[j]);
				}
			}
			if (aii == 0) {
				std::cout << "\ngs_solve failed due to no diagonal entry in line" << i << std::endl;
				return;
			}
			u.set_data(i, (1 / aii)*(b.get_data(i) - sum));
		}
		norm_res = 0;
		Au = (*this)*u;
		//mat_vec_multiply(u, Au); //is the function callable like this?
		for (int i = 0; i < b.get_size(); i++) {
			norm_res += (Au.get_data(i) - b.get_data(i))*(Au.get_data(i) - b.get_data(i));
		}
		std::cout << "iteration " << iter << ", residual " << norm_res << std::endl;
		if (norm_res < 1e-20) return;
	}
}

// 
CsrMatrix CsrMatrix::inv_diagonal() const {
	if (row_count != col_count) {
		std::cout << "Matrix has to be square!\n";
		exit(0); //WAS GIBT ES BESSERES ALS EXIT 0?
	}
	double* diag_entr = new double[row_count];
	if (diag_entr == NULL) {
		std::cout << "inv_diagonal failed to allocate memory\n";
		exit(0);
	}
	int* cols = new int[row_count];
	if (cols == NULL) {
		std::cout << "inv_diagonal failed to allocate memory\n";
		delete[] diag_entr;
		exit(0);
	}

	//set rows and cols
	for (int i = 0; i < row_count; i++) {
		cols[i] = i;
	}

	double aii = 0;
	//go throgh rows and set diag_entr to inverse of diagonal entries of Matrix
	for (int i = 0; i < row_count; i++) {
		aii = 0;
		for (int j = new_row_ind[i]; j < new_row_ind[i + 1]; j++) {
			if (col_ind[j] == i) {
				aii = entr[j];
				break;
			}
		}
		if (aii == 0) {
			std::cout << "\nNo diagonal entry in row " << i << std::endl;
			exit(0);
		}
		diag_entr[i] = 1 / aii;
	}
	CsrMatrix inv_diag_mat(row_count, col_count, row_count);

	//row indices and column indices are the same
	inv_diag_mat.csr_assemble(diag_entr, cols, cols, row_count);
	delete[] diag_entr;
	delete[] cols;
	//inv_diag_mat.print_csr_matrix(); //printing here yields correct result
	return inv_diag_mat;
}

Vector CsrMatrix::CG_Jac_prec(const Vector b, Vector x0) {
	Vector r0(0), p0(0);
	CsrMatrix Minv(0, 0, 0);
	Minv = (*this).inv_diagonal();
	r0 = (b - (*this)*x0);

}