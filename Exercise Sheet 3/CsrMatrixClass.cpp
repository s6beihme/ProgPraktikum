#include <iostream>
#include <memory>
#include <string>
#include <fstream>
#include "CsrMatrixClass.h"




CsrMatrix::CsrMatrix(int rows, int cols, int _nnz) {
	nnz = _nnz;
	row_count = rows;
	col_count = cols;

	entr = std::make_unique<double[]>(nnz);
	
	col_ind = std::make_unique<int[]>(nnz);

	new_row_ind = std::make_unique<int[]>(row_count + 1);
}

CsrMatrix::CsrMatrix(std::string filename) {
	std::ifstream myFile;
	myFile.open(filename);
	if (myFile.is_open()==false) {
		std::cout << "\nfailed to open file!\n";
		exit(0);
	}
	myFile >> row_count >> col_count >> nnz;
	if (row_count < 0 || col_count < 0 || nnz < 0) {
		std::cout << "dimensions and number of non zero entries in file have to be >=0\n";
		exit(0);
	}

	entr = std::make_unique<double[]>(nnz);
	col_ind = std::make_unique<int[]>(nnz);
	new_row_ind = std::make_unique<int[]>(row_count + 1);

	if (nnz == 0) return;
	new_row_ind[0] = 0;
	int cur_pos = 0;

	int row_last_entr;
	int row_cur_entr;
	myFile >> row_last_entr >> col_ind[0] >> entr[0];

	//fill up new_row_inds untill you reach first no-zero entry
	for (int i = 0; i < row_last_entr; i++) {
		cur_pos++;
		new_row_ind[cur_pos] = 0;
	}

	for (int i = 1; i < nnz; i++) {
		myFile >> row_cur_entr >> col_ind[i] >> entr[i];
		if (row_cur_entr > row_last_entr) {
			cur_pos++;
			new_row_ind[cur_pos] = i;
		}
		row_last_entr = row_cur_entr;
	}
	myFile.close();
	cur_pos++;
	new_row_ind[cur_pos] = nnz;
}

CsrMatrix::~CsrMatrix() {
	entr.release();
	col_ind.release();
	new_row_ind.release();
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


void CsrMatrix::mat_vec_mult(const Vector &v, Vector& result) const {
	if (col_count != v.get_size() || row_count!=result.get_size()) {
		std::cout << "\nMatrix Vector multiplication failed bacause dimension of argument didnt correspond to matrix" << std::endl;
		return;
	}

	double sum = 0;
	for (int i = 0; i < row_count; i++) {
		sum = 0;
		for (int j = new_row_ind[i]; j < new_row_ind[i + 1]; j++) {
			sum += entr[j] * v.get_data(col_ind[j]);
		}
		//res_vals[i] = sum;
		result.set_data(i, sum);
	}
}

void CsrMatrix::gs_solve(const Vector& b, Vector& result) const {
	Vector Au(b.get_size());
	Vector residual(b.get_size());
	long double norm_res = 0;

	//variable to store diagonal entry
	double aii = 0;
	int iter = 0;

	//variable to sum up non diagonal*entries of u or u2
	double sum = 0;

	while (iter < 1e+4) {
		iter++;

		for (int i = 0; i < row_count; i++) {
			sum = 0;
			aii = 0;
			for (int j = new_row_ind[i]; j < new_row_ind[i + 1]; j++) {
				if (col_ind[j] < i) {
					sum += entr[j] * result.get_data(col_ind[j]);
				}
				if (col_ind[j] == i) {
					aii = entr[j];
				}
				if (col_ind[j] > i) {
					sum += entr[j] * result.get_data(col_ind[j]);
				}
			}
			if (aii == 0) {
				std::cout << "\ngs_solve failed due to no diagonal entry in line" << i << std::endl;
				return;
			}
			result.set_data(i, (1 / aii)*(b.get_data(i) - sum));
		}
		//Au = (*this)*u;
		mat_vec_mult(result, Au);
		Au.subtr_vect(b, residual);
		norm_res = residual.norm_squared();
		//std::cout << "iteration " << iter << ", norm residual: " << norm_res << std::endl;
		
		if (norm_res < 1e-10) {
			std::cout << "\ngs_solve finished in iteration " << iter << std::endl;
			return;
		}
	}
	std::cout << "GS Iteration didnt reach close enough result, current (norm of residual)^2 is: " << norm_res << std::endl;
}


void CsrMatrix::solve_only_diag(const Vector& b, Vector& result) const{
	if (result.get_size() != col_count || b.get_size() != row_count || row_count != col_count) {
		std::cout << "\nsolve only diag failed because dimensions didnt fit\n";
		return;
	}
	for (int i = 0; i < row_count; i++) {
		for (int j = new_row_ind[i]; j < new_row_ind[i + 1]; j++) {
			if (col_ind[j] == i) {
				result.set_data(i, (b.get_data(i) / entr[j]));
				break;
			}
		}
	}
}

void CsrMatrix::solve_only_upper_triang(const Vector& b, Vector& result) {
	if (result.get_size() != col_count || b.get_size() != row_count || row_count != col_count) {
		std::cout << "\nsolve only upper triag failed because dimensions didnt fit\n";
		return;
	}

	double aii = 0;
	double sum = 0;
	for (int i = row_count - 1; i >= 0; i--) {
		aii = 0;
		sum = 0;
		for (int j = new_row_ind[i]; j < new_row_ind[i + 1]; j++) {
			if (col_ind[j] < i) continue;
			if (col_ind[j] == i) {
				aii = entr[j];
			}
			if (col_ind[j] > i) {
				sum +=entr[j]*result.get_data(col_ind[j]);
			}
		}
		if (aii == 0) {
			std::cout << "\nsolve only upper triag failed due to no diagonal entry in line " << i << std::endl;
			return;
		}
		result.set_data(i, (b.get_data(i)-sum)/aii);
	}
}

void CsrMatrix::CG_prec(const Vector& b, Vector& x0, int preconditioner) {
	if (preconditioner != 0 && preconditioner != 1) {
		std::cout << "CG_prec failed: integer giving preconditioner has to be 0 or one\n";
		return;
	}
	Vector rk(row_count), pk(row_count), zk(row_count), zk_1(row_count), rk_1(row_count);
	
	Vector temp1(row_count), temp2(row_count);

	long double res_norm_squared = 0;
	double alpha = 0;
	double betta = 0;

	//for explanation see https://en.wikipedia.org/wiki/Conjugate_gradient_method
	
	//temp1=(*this)*x0
	mat_vec_mult(x0, temp1);
	//rk=b-temp1
	b.subtr_vect(temp1, rk);
	
	//zk=M^-1*rk
	if (preconditioner == 0) {
		solve_only_diag(rk, zk);
	}
	if (preconditioner == 1) {
		solve_only_upper_triang(rk, zk);
	}

	//pk = zk
	zk.make_copy(pk);

	int k = 0;
	while (k < 1e+4) {

		//temp1=(*this*pk)
		mat_vec_mult(pk, temp1);

		alpha = (rk*zk)/(pk*temp1);

		//temp1=alpha*pk
		scalar_mult(alpha, pk, temp1);

		//temp2=xk+temp=xk+alpha*pk
		x0.add_vect(temp1, temp2); //this may couse problems

		//xk = temp2 = xk+alpha*pk
		temp2.make_copy(x0);

		//temp1=(*this)*pk
		mat_vec_mult(pk, temp1);

		//temp2=alpha*temp1=alpha*(*this)*pk
		scalar_mult(alpha, temp1, temp2);

		//rk_1=rk-alpha*(*this)*pk
		rk.subtr_vect(temp2, rk_1);

		//check for convergence
		res_norm_squared = rk.norm_squared();
		std::cout << "iteration " << k << ": norm_res " << res_norm_squared << std::endl;
		if (res_norm_squared < 1e-10) {
			std::cout << "CG_prec finished in iteration " << k << std::endl;
			return;
		}

		//zk_1=M^-1*rk_1
		if (preconditioner == 0) {
			solve_only_diag(rk_1, zk_1);
		}
		if (preconditioner == 1) {
			solve_only_upper_triang(rk_1, zk_1);
		}

		betta = (zk_1*rk_1) / (zk*rk);

		//temp1=betta*pk
		scalar_mult(betta, pk, temp1);

		//pk = zk_1 + (betta*pk)
		zk_1.add_vect(temp1, pk);
		

		//update zk and rk
		zk_1.make_copy(zk);
		rk_1.make_copy(rk);

		k++;
	}
	std::cout << "CG Iteration didnt reach close enough result, current (norm of residual)^2 is: " << res_norm_squared << std::endl;
}

void CsrMatrix::Richardson_prec(const Vector& b, Vector& x0, int preconditioner) {
	if (preconditioner != 0 && preconditioner != 1) {
		std::cout << "CG_prec failed: integer giving preconditioner has to be 0 or one\n";
		return;
	}

	Vector temp1(row_count), temp2(row_count);
	long double res_norm_squared;
	int k = 0;
	while (k < 1e+4) {
		//temp1=(*this)*x0
		mat_vec_mult(x0, temp1);

		//temp2=temp1-b=(*this)*x0-b
		temp1.subtr_vect(b, temp2);

		//temp1=P^-1*temp2=P^-1*((*this)*x0-b)
		if (preconditioner == 0) {
			solve_only_diag(temp2, temp1);
		}
		if (preconditioner == 1) {
			solve_only_upper_triang(temp2, temp1);
		}

		//temp2=x0-temp1=x0-P^-1*((*this)*x0-b)
		x0.subtr_vect(temp1, temp2);

		temp2.make_copy(x0);

		//test for convergence
		mat_vec_mult(x0, temp1);
		b.subtr_vect(temp1, temp2);
		
		res_norm_squared = temp2.norm_squared();
		std::cout << "iteration " << k << ", norm residual: " << res_norm_squared << std::endl;
		if (res_norm_squared < 1e-10) {
			std::cout << "Richardson_prec finished in iteratio " << k;
			return;
		}

		k++;
	}
	std::cout << "Richardson Iteration didnt reach close enough result, current (norm of residual)^2 is: " << res_norm_squared << std::endl;
}