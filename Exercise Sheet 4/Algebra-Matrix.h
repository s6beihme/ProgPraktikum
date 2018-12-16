#pragma once
#include <iostream>
#include <memory>
#include <string>
#include <fstream>
#include "VectorClass.h"

template <typename T>
class CsrMatrix {
public:
	CsrMatrix();
	//Copy Constructor
	CsrMatrix(const CsrMatrix& other);

	//allocates memory for CsrMatrix object
	//parameters:
		//_number_non_zero: number of non-zero entries in matrix
		//_row_count: number of rows in Matrix
		//_column_count: -"- cols -"-
	CsrMatrix(int _number_non_zero, int _row_count, int _column_count);

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

	//assignment operator
	CsrMatrix& operator=(const CsrMatrix& other);

	//test for equality
	bool operator ==(const CsrMatrix& other);
	//test for inequality
	bool operator !=(const CsrMatrix& other);

	int get_row_count() const { return row_count; }
	int get_column_count() const { return column_count; }
	int get_number_non_zero() const { return number_non_zero; }

	//fills entr, column_ind and new_row_ind according to input
	//parameters:
		//values: array containing the values that should be stored in the matrix
		//row_inds: array containing the row indices of the values -"-
		//column_inds: array containing the column indices of the values-"-
		//number_non_zero: number of given non-zero entries
	void assemble(T* values, int* row_inds, int*column_inds, int _number_non_zero);

	//matrix vector produkt
	//parameters:
		//v: vector with which the matrix is multiplied (has to have correct dimension)
		//result: functions fills result vector with (*this)*v
	void matrix_vector_multiply(Vector<T> &v, Vector<T>& result) const;

	//solves Linear system of Equations Au=b for x using gauss seidel iteration
	//parameters:
		//result: initial guess vector of result. it will be overwritten with the result
		//b: right hand side vector of Equation
	void gauss_seidel_solve(Vector<T>& b, Vector<T>& result) const;

	//solves linear system of equations Mx=b for x algebraically where M is diagonal part of (*this)
	//parameters:
		//b: right hand side vector of equation
		//result: is filled with result such that M*result=b
	void solve_only_diag(const Vector<T>& b, Vector<T>& result) const;

	//solves linear system of equations Mx=b for x algebraically where M is upper triagonal part of (*this)
	//parameters:
		//b: right hand side vector of equation
		//result: is filled with result such that M*result=b
	void solve_only_upper_triang(const Vector<T>& b, Vector<T>& result);

	//solves linear system of equations Ax=b for x using preconditioned Conjugate Gradient Method
	//parameters:
		//b: right hand side vector of equation
		//x0: initial guess. result will be stored in here
		//preconditioner: if =0 then uses Jacobi preconditioner, if =1 uses gauss seidel preconditioner
	//returns: Vector type object x that satisfies Ax=b
	void CG_preconditioner(const Vector<T>& b, Vector<T>& x0, int preconditioner);

	//solves linear system of equations Ax=b for x using preconditioned Richardson Method (gamma_n=1 for all n)
	//parameters:
		//b: right hand side vector of equation
		//x0: initial guess. result will be stored in here
		//preconditioner: if =0 then uses Jacobi preconditioner, if =1 uses gauss seidel preconditioner
	//returns: Vector type object x that satisfies Ax=b
	void Richardson_preconditioner(const Vector<T>& b, Vector<T>& x0, int preconditioner);

	template <typename T2>
	friend std::ostream& operator << (std::ostream& out, const CsrMatrix<T2>& A);
private:
	int number_non_zero;
	int row_count;
	int column_count;
	//double* entr;
	std::unique_ptr<T[]> entr;
	//int* column_ind;
	std::unique_ptr<int[]> column_ind;
	//int* new_row_ind;
	std::unique_ptr<int[]> new_row_ind;
};

template <typename T>
CsrMatrix<T>::CsrMatrix() :
	number_non_zero(0), row_count(0), column_count(0), entr(nullptr), column_ind(nullptr), new_row_ind(nullptr)
{}

template <typename T>
CsrMatrix<T>::CsrMatrix(const CsrMatrix<T>& other) :
	number_non_zero(other.number_non_zero), row_count(other.row_count), column_count(other.column_count)
{
	entr = std::make_unique<T[]>(number_non_zero);
	column_ind = std::make_unique<int[]>(number_non_zero);
	new_row_ind = std::make_unique<int[]>(row_count + 1);

	for (int i = 0; i < number_non_zero; i++) {
		entr[i] = other.entr[i];
		column_ind[i] = other.column_ind[i];
	}

	for (int i = 0; i < row_count + 1; i++) new_row_ind[i] = other.new_row_ind[i];
}

template <typename T>
CsrMatrix<T>::CsrMatrix(int _number_non_zero, int _row_count, int _column_count) :
	number_non_zero(_number_non_zero), row_count(_row_count), column_count(_column_count)
{
	entr = std::make_unique<T[]>(number_non_zero);
	column_ind = std::make_unique<int[]>(number_non_zero);
	new_row_ind = std::make_unique<int[]>(row_count + 1);
}

template <typename T>
CsrMatrix<T>::CsrMatrix(std::string filename) {
	std::ifstream myFile;
	myFile.open(filename);
	if (myFile.is_open() == false) {
		std::cout << "\nfailed to open file!\n";
		exit(0);
	}
	myFile >> row_count >> column_count >> number_non_zero;
	if (row_count < 0 || column_count < 0 || number_non_zero < 0) {
		std::cout << "dimensions and number of non zero entries in file have to be >=0\n";
		exit(0);
	}

	entr = std::make_unique<T[]>(number_non_zero);
	column_ind = std::make_unique<int[]>(number_non_zero);
	new_row_ind = std::make_unique<int[]>(row_count + 1);

	if (number_non_zero == 0) return;
	new_row_ind[0] = 0;
	int cursor_position = 0;

	int row_last_entr;
	int row_cursor_entr;
	myFile >> row_last_entr >> column_ind[0] >> entr[0];

	//fill up new_row_inds untill you reach first no-zero entry
	for (int i = 0; i < row_last_entr; i++) {
		cursor_position++;
		new_row_ind[cursor_position] = 0;
	}

	for (int i = 1; i < number_non_zero; i++) {
		myFile >> row_cursor_entr >> column_ind[i] >> entr[i];
		if (row_cursor_entr > row_last_entr) {
			cursor_position++;
			new_row_ind[cur_pos] = i;
		}
		row_last_entr = row_cursor_entr;
	}
	myFile.close();
	cursor_position++;
	new_row_ind[row_count] = number_non_zero;
}

template <typename T>
CsrMatrix<T>& CsrMatrix<T>::operator=(const CsrMatrix<T>& other) {
	number_non_zero = other.number_non_zero;
	column_count = other.column_count;
	row_count = other.row_count;
	entr = std::make_unique<T[]>(number_non_zero);
	column_ind = std::make_unique<int[]>(number_non_zero);
	new_row_ind = std::make_unique<int[]>(row_count + 1);

	for (int i = 0; i < number_non_zero; i++) {
		entr[i] = other.entr[i];
		column_ind[i] = other.column_ind[i];
	}

	for (int i = 0; i < row_count + 1; i++) new_row_ind[i] = other.new_row_ind[i];
	return *this;
}

template <typename T>
bool CsrMatrix<T>::operator ==(const CsrMatrix<T>& other) {
	if (number_non_zero != other.number_non_zero || row_count != other.row_count || column_count != other.column_count) return false;
	for (int i = 0; i < number_non_zero; i++) if (entr[i] != other.entr[i] || column_ind[i] != other.column_ind[i]) return false;
	for (int i = 0; i < row_count + 1; i++) if (new_row_ind[i] != other.new_row_ind[i]) return false;
	return true;
}

template <typename T>
bool CsrMatrix<T>::operator !=(const CsrMatrix<T>& other) {
	return !((*this) == other);
}

template <typename T>
std::ostream& operator << (std::ostream& out, const CsrMatrix<T>& A) { //binary operator, why this syntax?
	out << "\nMatrix with entries:\n";
	for (int i = 0; i < A.number_non_zero; i++) {
		out << "\n" << A.entr[i] << ", in column:" << A.column_ind[i];
	}
	out << "\n\nThe new_row_ind array is:\n";
	for (int i = 0; i < A.row_count + 1; i++) {
		out << A.new_row_ind[i] << ", ";
	}
	out << "\n";
	return out;
}

template<typename T>
void CsrMatrix<T>::assemble(T* values, int* row_inds, int*column_inds, int _number_non_zero) {
	if (_number_non_zero != number_non_zero) {
		std::cout << "csr_assemble failed because given number_non_zero didnt match number_non_zero of matrix\n";
		return;
	}

	//copy values and column_inds
	for (int i = 0; i < number_non_zero; i++) {
		entr[i] = values[i];
		column_ind[i] = column_inds[i];
	}

	new_row_ind[0] = 0;

	int cursor_position = 0;

	//fill up new_row_inds untill you reach first no-zero entry
	for (int i = 0; i < row_inds[0]; i++) {
		cursor_position++;
		new_row_ind[cursor_position] = 0;
	}

	for (int i = 1; i < number_non_zero; i++) {
		if (row_inds[i] > row_inds[i - 1]) {
			cursor_position++;
			new_row_ind[cursor_position] = i;
		}
	}
	cursor_position++;
	new_row_ind[cursor_position] = number_non_zero;
}

template<typename T>
void CsrMatrix<T>::matrix_vector_multiply(Vector<T> &v, Vector<T>& result) const {
	if (column_count != v.get_size() || row_count != result.get_size()) {
		std::cout << "\nMatrix Vector multiplication failed bacause dimension of argument didnt correspond to matrix" << std::endl;
		return;
	}

	double sum = 0;
	for (int i = 0; i < row_count; i++) {
		sum = 0;
		for (int j = new_row_ind[i]; j < new_row_ind[i + 1]; j++) {
			sum += entr[j] * v[column_ind[j]];
		}
		//res_values[i] = sum;
		result[i] = sum;
	}
}



template<typename T>
void CsrMatrix<T>::gauss_seidel_solve(Vector<T>& b, Vector<T>& result) const {
	Vector<T> Au(b.get_size());
	Vector<T> residual(b.get_size());
	T norm_residual = 0;

	//variable to store diagonal entry
	T aii = 0;
	int iteration = 0;

	//variable to sum up non diagonal*entries of u or u2
	T sum = 0;

	while (iteration < 1000) {
		iteration++;

		for (int i = 0; i < row_count; i++) {
			sum = 0;
			aii = 0;
			for (int j = new_row_ind[i]; j < new_row_ind[i + 1]; j++) {
				if (column_ind[j] < i) {
					sum += entr[j] * result[column_ind[j]];
				}
				if (column_ind[j] == i) {
					aii = entr[j];
				}
				if (column_ind[j] > i) {
					sum += entr[j] * result[column_ind[j]];
				}
			}
			if (aii == 0) {
				std::cout << "\ngauss_seidel_solve failed due to no diagonal entry in line" << i << std::endl;
				return;
			}
			result[i] = (1 / aii)*(b[i] - sum);
		}
		//Au = (*this)*u;
		matrix_vector_multiply(result, Au);
		subtr_vect(Au, b, residual);
		norm_residual = residual.norm_squared();
		//std::cout << "iteration " << iteration << ", norm residual: " << norm_residual << std::endl;

		if (norm_residual < 1e-17) {
			std::cout << "\ngauss_seidel_solve finished in iteration " << iteration << std::endl;
			return;
		}
	}
	std::cout << "gauss_seidel Iteration didnt reach close enough result, current (norm of residual)^2 is: " << norm_res << std::endl;
}

template<typename T>
void CsrMatrix<T>::solve_only_diag(const Vector<T>& b, Vector<T>& result) const {
	if (result.get_size() != column_count || b.get_size() != row_count || row_count != column_count) {
		std::cout << "\nsolve only diag failed because dimensions didnt fit\n";
		return;
	} // Assert
	for (int i = 0; i < row_count; i++) {
		for (int j = new_row_ind[i]; j < new_row_ind[i + 1]; j++) {
			if (column_ind[j] == i) {
				result.set_data(i, (b.get_data(i) / entr[j]));
				break;
			}
		}
	}
}

template<typename T>
void CsrMatrix<T>::solve_only_upper_triang(const Vector<T>& b, Vector<T>& result) {
	if (result.get_size() != column_count || b.get_size() != row_count || row_count != column_count) {
		std::cout << "\nsolve only upper triag failed because dimensions didnt fit\n";
		return;
	}

	T aii = 0;
	T sum = 0;
	for (int i = row_count - 1; i >= 0; i--) {
		aii = 0;
		sum = 0;
		for (int j = new_row_ind[i]; j < new_row_ind[i + 1]; j++) {
			if (column_ind[j] < i) continue;
			if (column_ind[j] == i) {
				aii = entr[j];
			}
			if (column_ind[j] > i) {
				sum += entr[j] * result.get_data(column_ind[j]);
			}
		}
		if (aii == 0) {
			std::cout << "\nsolve only upper triag failed due to no diagonal entry in line " << i << std::endl;
			return;
		}
		result.set_data(i, (b.get_data(i) - sum) / aii);
	}
}

template<typename T>
void CsrMatrix<T>::CG_preconditioner(const Vector<T>& b, Vector<T>& x0, int preconditioner) {
	if (preconditioner != 0 && preconditioner != 1) {
		std::cout << "CG_preconditioner failed: integer giving preconditioner has to be 0 or one\n";
		return;
	}
	Vector<T> rk(row_count), pk(row_count), zk(row_count), zk_1(row_count), rk_1(row_count);

	Vector<T> temp1(row_count), temp2(row_count);

	long double residual_norm_squared = 0;
	T alpha = 0;
	T beta = 0;

	//for explanation see https://en.wikipedia.org/wiki/Conjugate_gradient_method

	//temp1=(*this)*x0
	matrix_vector_multiply(x0, temp1);
	//rk=b-temp1
	b.subtr_vect(temp1, rk);

	//zk=M^-1*rk
	if (preconditioner == 0) {
		solve_only_diag(rk, zk);
	}
	if (preconditioner == 1) {
		solve_only_upper_triang(rk, zk);
	}


	pk = zk;

	int k = 0;
	while (k < 1e+4) {

		//temp1=(*this*pk)
		matrix_vector_multiply(pk, temp1);

		alpha = (rk*zk) / (pk*temp1);

		//temp1=alpha*pk
		scalar_mult(alpha, pk, temp1);

		//temp2=xk+temp=xk+alpha*pk
		x0.add_vect(temp1, temp2); //this may couse problems

		//xk = temp2 = xk+alpha*pk
		temp2.make_copy(x0);

		//temp1=(*this)*pk
		matrix_vector_multiply(pk, temp1);

		//temp2=alpha*temp1=alpha*(*this)*pk
		scalar_mult(alpha, temp1, temp2);

		//rk_1=rk-alpha*(*this)*pk
		rk.subtr_vect(temp2, rk_1);

		//check for convergence
		residual_norm_squared = rk.norm_squared();
		std::cout << "iteration " << k << ": norm_residual " << residual_norm_squared << std::endl;
		if (residual_norm_squared < 1e-10) {
			std::cout << "CG_preconditioner finished in iteration " << k << std::endl;
			return;
		}

		//zk_1=M^-1*rk_1
		if (preconditioner == 0) {
			solve_only_diag(rk_1, zk_1);
		}
		if (preconditioner == 1) {
			solve_only_upper_triang(rk_1, zk_1);
		}

		beta = (zk_1*rk_1) / (zk*rk);

		//temp1=betta*pk
		scalar_mult(beta, pk, temp1);

		//pk = zk_1 + (betta*pk)
		zk_1.add_vect(temp1, pk);


		//update zk and rk
		zk = zk_1;        //Copyconstricter?
		rk = rk_1;

		k++;
	}
	std::cout << "CG Iteration didnt reach close enough result, current (norm of residual)^2 is: " << residual_norm_squared << std::endl;
}

template<typename T>
void CsrMatrix<T>::Richardson_preconditioner(const Vector<T>& b, Vector<T>& x0, int preconditioner) {
	if (preconditioner != 0 && preconditioner != 1) {
		std::cout << "CG_prec failed: integer giving preconditioner has to be 0 or one\n";
		return;
	}

	Vector<T> temp1(row_count), temp2(row_count);
	long double residual_norm_squared;
	int k = 0;
	while (k < 1e+4) {
		//temp1=(*this)*x0
		matrix_vector_multiply(x0, temp1);

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

		x0 = temp2;

		//test for convergence
		matrix_vector_multiply(x0, temp1);
		b.subtr_vect(temp1, temp2);

		residual_norm_squared = temp2.norm_squared();
		std::cout << "iteration " << k << ", norm residual: " << residual_norm_squared << std::endl;
		if (residual_norm_squared < 1e-10) {
			std::cout << "Richardson_preconditioner finished in iteration " << k;
			return;
		}

		k++;
	}
	std::cout << "Richardson Iteration didnt reach close enough result, current (norm of residual)^2 is: " << residual_norm_squared << std::endl;
}