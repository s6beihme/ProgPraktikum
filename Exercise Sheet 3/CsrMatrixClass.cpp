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


int CsrMatrix::get_row_count() {
	return row_count;
}
int CsrMatrix::get_col_count() {
	return col_count;
}
int CsrMatrix::get_nnz() {
	return nnz;
}
double CsrMatrix::get_entr(int i) {
	return entr[i];
}


void CsrMatrix::print_csr_matrix() {
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


void CsrMatrix::mat_vec_multiply(Vector &v, Vector& res) {
	if (col_count != v.get_size()) {
		std::cout << "\nMatrix Vector multiplication failed bacause dimension of first argument didnt correspond to matrix" << std::endl;
		exit(0);
	}
	if (row_count != res.get_size()) {
		std::cout << "\nMatrix Vector multiplication failed, because dimension of second argument didnt correspond to matrix" << std::endl;
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
	res.vec_assemble(res_vals);
	delete[] res_vals;
}

void CsrMatrix::gs_solve(Vector& u, Vector& b) {
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
		mat_vec_multiply(u, Au); //is the function callable like this?
		for (int i = 0; i < b.get_size(); i++) {
			norm_res += (Au.get_data(i) - b.get_data(i))*(Au.get_data(i) - b.get_data(i));
		}
		std::cout << "iteration " << iter << ", residual " << norm_res << std::endl;
		if (norm_res < 1e-20) return;
	}
}

