#include <iostream>
#include <mpi.h>
#include <memory>

template <typename T>
void print_matrix(int row_count, int column_count, std::unique_ptr<std::unique_ptr<T[]>[]>& U) {
	for (int i = 0; i < row_count; i++) {
		for (int j = 0; j < column_count; j++) {
			std::cout << U[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

MPI_Datatype get_MPI_Datatype(double) {
	return MPI_DOUBLE;
}

MPI_Datatype get_MPI_Datatype(float) {
	return MPI_FLOAT;
}

//sends the highest row, that is not the boundary row, to Process target
//row is a type T array of length jmax
template <typename T>
void send_upper_row(std::unique_ptr<std::unique_ptr<T[]>[]>& M, T* row, int, int jmax, int target) {
	//copy the entries of matrix into array
	for (int j = 1; j < jmax + 1; j++) {
		row[j - 1] = M[1][j];
	}
	//send row to target process
	T a = 0;
	MPI_Send(row, jmax, get_MPI_Datatype(a), target, 110, MPI_COMM_WORLD);
}

//sends the lowest row, that is not the boundary row, to Process target
//row is a type T array of length jmax
template <typename T>
void send_lower_row(std::unique_ptr<std::unique_ptr<T[]>[]>& M, T* row, int imax, int jmax, int target) {
	//copy the entries of matrix into array
	for (int j = 1; j < jmax + 1; j++) {
		row[j - 1] = M[imax][j];
	}
	//send row to target process
	T a = 0;
	MPI_Send(row, jmax, get_MPI_Datatype(a), target, 110, MPI_COMM_WORLD);
}

//sends the most left column, that is not the boundary column to process target
//left is a type T array of length imax
template <typename T>
void send_left_column(std::unique_ptr<std::unique_ptr<T[]>[]>& M, T* left, int imax, int, int target) {
	//copy the entries of matrix into array
	for (int i = 1; i < imax + 1; i++) {
		left[i - 1] = M[i][1];
	}
	//send column to target process
	T a=0;
	MPI_Send(left, imax, get_MPI_Datatype(a), target, 110, MPI_COMM_WORLD);
}

//sends the most right column, that is not the boundary column to process target
//right is a type T array of length imax
template <typename T>
void send_right_column(std::unique_ptr<std::unique_ptr<T[]>[]>& M, T* right, int imax, int jmax , int target) {
	//copy the entries of matrix into array
	for (int i = 1; i < imax + 1; i++) {
		right[i - 1] = M[i][jmax];
	}
	//send column to target process
	T a = 0;
	MPI_Send(right, imax, get_MPI_Datatype(a), target, 110, MPI_COMM_WORLD);
}

//receives the lower boundary values of matrix M and sets them
//row is an array of type T of length jmax
template <typename T>
void receive_lower_boundary_and_set_matrix(std::unique_ptr<std::unique_ptr<T[]>[]>& M, T* row, int imax, int jmax, int source) {
	MPI_Status status;
	T a = 0;
	MPI_Recv(row, jmax, get_MPI_Datatype(a), source, 110, MPI_COMM_WORLD, &status);

	for (int j = 1; j < jmax + 1; j++) {
		M[imax+1][j] = row[j - 1];
	}
}

//receives the upper boundary values of matrix M and sets them
//row is an array of type T of length jmax
template <typename T>
void receive_upper_boundary_and_set_matrix(std::unique_ptr<std::unique_ptr<T[]>[]>& M, T* row, int , int jmax, int source) {
	MPI_Status status;
	T a = 0;
	MPI_Recv(row, jmax, get_MPI_Datatype(a), source, 110, MPI_COMM_WORLD, &status);

	for (int j = 1; j < jmax + 1; j++) {
		M[0][j] = row[j - 1];
	}
}

//receives the right boundary values of matrix M and sets them
//right is an array of type T of length imax
template <typename T>
void receive_right_boundary_and_set_matrix(std::unique_ptr<std::unique_ptr<T[]>[]>& M, T* right, int imax, int jmax, int source) {
	MPI_Status status;
	T a=0;
	MPI_Recv(right, imax, get_MPI_Datatype(a), source, 110, MPI_COMM_WORLD, &status);

	for (int i = 1; i < imax + 1; i++) {
		M[i][jmax + 1] = right[i - 1];
	}
}

//receives the left boundary values of matrix M and sets them
//left is an array of type T of length imax
template <typename T>
void receive_left_boundary_and_set_matrix(std::unique_ptr<std::unique_ptr<T[]>[]>& M, T* left, int imax, int , int source) {
	MPI_Status status;
	T a = 0;
	MPI_Recv(left, imax, get_MPI_Datatype(a), source, 110, MPI_COMM_WORLD, &status);

	for (int i = 1; i < imax + 1; i++) {
		M[i][0] = left[i - 1];
	}
}

int main(int argc, char **argv) {
	int  myrank;
	int imax = 3;
	int jmax = 3;
	double* row = new double[jmax];

	std::unique_ptr<std::unique_ptr<double[]>[]> M = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		M[i] = std::make_unique<double[]>(jmax + 2);
	}
	for (int i = 1; i < imax + 1; i++) M[3][i] = 3;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if (myrank == 0) {
		send_lower_row<double>(M, row, imax, jmax, 1);
	}
	else if (myrank == 1) {
		receive_upper_boundary_and_set_matrix<double>(M, row, imax, jmax, 0);
		print_matrix(imax + 2, jmax + 2, M);
	}
	MPI_Finalize();
	delete[] row;
}