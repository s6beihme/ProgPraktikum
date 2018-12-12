#pragma once
#include <iostream>
#include <memory>
#include <string>
#include <fstream>

template <typename T>
class CsrMatrix {
public:
	CsrMatrix();
	CsrMatrix(const CsrMatrix& other);
	CsrMatrix(int _nnz, int _row_count, int _col_count);
	CsrMatrix(std::string filename);

	CsrMatrix& operator=(const CsrMatrix& other);

	bool operator ==(const CsrMatrix& other);
	bool operator !=(const CsrMatrix& other);

	int get_row_count() const { return row_count; }
	int get_col_count() const { return col_count; }
	int get_nnz() const { return nnz; }

	template <typename T2>
	friend std::ostream& operator << (std::ostream& out, const CsrMatrix<T2>& A);
private:
	int nnz;
	int row_count;
	int col_count;
	std::unique_ptr<T[]> entr;
	std::unique_ptr<int[]> col_ind;
	std::unique_ptr<int[]> new_row_ind;
};

template <typename T>
CsrMatrix<T>::CsrMatrix() :
	nnz(0), row_count(0), col_count(0), entr(nullptr), col_ind(nullptr), new_row_ind(nullptr)
{}

template <typename T>
CsrMatrix<T>::CsrMatrix(const CsrMatrix<T>& other) :
	nnz(other.nnz), row_count(other.row_count), col_count(other.col_count)
{
	entr = std::make_unique<T[]>(nnz);
	col_ind = std::make_unique<int[]>(nnz);
	new_row_ind = std::make_unique<int[]>(row_count + 1);

	for (int i = 0; i < nnz; i++) {
		entr[i] = other.entr[i];
		col_ind[i] = other.col_ind[i];
	}

	for (int i = 0; i < row_count + 1; i++) new_row_ind[i] = other.new_row_ind[i];
}

template <typename T>
CsrMatrix<T>::CsrMatrix(int _nnz, int _row_count, int _col_count) :
	nnz(_nnz), row_count(_row_count), col_count(_col_count)
{
	entr = std::make_unique<T[]>(nnz);
	col_ind = std::make_unique<int[]>(nnz);
	new_row_ind = std::make_unique<int[]>(row_count+1);
}

template <typename T>
CsrMatrix<T>::CsrMatrix(std::string filename) {
	std::ifstream myFile;
	myFile.open(filename);
	if (myFile.is_open() == false) {
		std::cout << "\nfailed to open file!\n";
		exit(0);
	}
	myFile >> row_count >> col_count >> nnz;
	if (row_count < 0 || col_count < 0 || nnz < 0) {
		std::cout << "dimensions and number of non zero entries in file have to be >=0\n";
		exit(0);
	}

	entr = std::make_unique<T[]>(nnz);
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
	new_row_ind[row_count] = nnz;
}

template <typename T>
CsrMatrix<T>& CsrMatrix<T>::operator=(const CsrMatrix<T>& other) {
	nnz = other.nnz;
	col_count = other.col_count;
	row_count = other.row_count;
	entr = std::make_unique<T[]>(nnz);
	col_ind = std::make_unique<int[]>(nnz);
	new_row_ind = std::make_unique<int[]>(row_count + 1);

	for (int i = 0; i < nnz; i++) {
		entr[i] = other.entr[i];
		col_ind[i] = other.col_ind[i];
	}

	for (int i = 0; i < row_count + 1; i++) new_row_ind[i] = other.new_row_ind[i];
	return *this;
}

template <typename T>
bool CsrMatrix<T>::operator ==(const CsrMatrix<T>& other) {
	if (nnz != other.nnz || row_count != other.row_count || col_count != other.col_count) return false;
	for (int i = 0; i < nnz; i++) if (entr[i] != other.entr[i] || col_ind[i] != other.col_ind[i]) return false;
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
	for (int i = 0; i < A.nnz; i++) {
		out << "\n" << A.entr[i] << ", in column:" << A.col_ind[i];
	}
	out << "\n\nThe new_row_ind array is:\n";
	for (int i = 0; i < A.row_count + 1; i++) {
		out << A.new_row_ind[i] << ", ";
	}
	out << "\n";
	return out;
}