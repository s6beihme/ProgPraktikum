#pragma once
#include <iostream>
#include <mpi.h>
#include "driven_cavity.h"


//A few things regarding this programm:

//Whenever I pass a T* row or T* column to a function, those arrays dont have to contain any
//relevant data. They are only used a a container to use in MPI_Send or MPI_Recv and by passing
//them, I avoid allocating memory for them in each call of one of the functions.
//row has to have length jmax (the number of INNER columns of the matrix) and 
//column has to have length imax (the number of INNER rows of the matrix)

//When I pass int* coords or int* dims to a function, I assume that they have already been 
//calculated in the main programm.

//Same as in my last programm, imax and jmax are the number of the INNER rows and INNER columns respectively of the 
//matrices, NOT the dimensions of the whole matrices.





MPI_Datatype get_MPI_Datatype(double) {
	return MPI_DOUBLE;
}

MPI_Datatype get_MPI_Datatype(float) {
	return MPI_FLOAT;
}

//sends the highest row, that is not the boundary row, to Process target
//parameters: 
	//M: matrix, whose upper row is to be sent
	//row: type T array of length jmax to store the upper row of matrix in to be sent
	//jmax: M has j+2 columns (jmax "real" ones)
	//target: the rank of the process, to which the data is sent
template <typename T>
void send_upper_boundary(std::unique_ptr<std::unique_ptr<T[]>[]>& M, T* row, int, int jmax, int target) {
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
void send_lower_boundary(std::unique_ptr<std::unique_ptr<T[]>[]>& M, T* row, int imax, int jmax, int target) {
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
void send_left_boundary(std::unique_ptr<std::unique_ptr<T[]>[]>& M, T* left, int imax, int, int target) {
	//copy the entries of matrix left boundary into array
	for (int i = 1; i < imax + 1; i++) {
		left[i - 1] = M[i][1];
	}
	//send column to target process
	T a = 0;
	MPI_Send(left, imax, get_MPI_Datatype(a), target, 110, MPI_COMM_WORLD);
}

//sends the most right column, that is not the boundary column to process target
//right is a type T array of length imax
template <typename T>
void send_right_boundary(std::unique_ptr<std::unique_ptr<T[]>[]>& M, T* right, int imax, int jmax, int target) {
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
		M[imax + 1][j] = row[j - 1];
	}
}

//receives the upper boundary values of matrix M and sets them
//row is an array of type T of length jmax
template <typename T>
void receive_upper_boundary_and_set_matrix(std::unique_ptr<std::unique_ptr<T[]>[]>& M, T* row, int, int jmax, int source) {
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
	T a = 0;
	MPI_Recv(right, imax, get_MPI_Datatype(a), source, 110, MPI_COMM_WORLD, &status);

	for (int i = 1; i < imax + 1; i++) {
		M[i][jmax + 1] = right[i - 1];
	}
}

//receives the left boundary values of matrix M and sets them
//left is an array of type T of length imax
template <typename T>
void receive_left_boundary_and_set_matrix(std::unique_ptr<std::unique_ptr<T[]>[]>& M, T* left, int imax, int, int source) {
	MPI_Status status;
	T a = 0;
	MPI_Recv(left, imax, get_MPI_Datatype(a), source, 110, MPI_COMM_WORLD, &status);

	for (int i = 1; i < imax + 1; i++) {
		M[i][0] = left[i - 1];
	}
}


template <typename T>
void exchange_boundary_values_UV_or_FG_between_processes(MPI_Comm cartesian_handle ,std::unique_ptr<std::unique_ptr<T[]>[]>& U, std::unique_ptr<std::unique_ptr<T[]>[]>& V, T* column, T* row, int imax, int jmax, int*coords, int* dims) {
	int target = 0;
	int target_coords[2] = { 0,0 };
	int source = 0;
	int source_coords[2] = { 0,0 };

	//send left boundary U
	if (coords[1] >= 1) {
		target_coords[0] = coords[0];
		target_coords[1] = coords[1] - 1;
		MPI_Cart_rank(cartesian_handle, target_coords, &target);
		send_left_boundary(U, column, imax, jmax, target);
	}
	//receive right boundary U
	if (coords[1] <= dims[1] - 2) {
		source_coords[0] = coords[0];
		source_coords[1] = coords[1] + 1;
		MPI_Cart_rank(cartesian_handle, source_coords, &source);
		receive_right_boundary_and_set_matrix(U, column, imax, jmax, source);
	}
	//send left boundary V
	if (coords[1] >= 1) {
		send_left_boundary(V, column, imax, jmax, target);
	}
	//receive right boundary V
	if (coords[1] <= dims[1] - 2) {
		receive_right_boundary_and_set_matrix(V, column, imax, jmax, source);
	}

	//send right boundary U 
	if (coords[1] <= dims[1] - 2) {
		target_coords[0] = coords[0];
		target_coords[1] = coords[1] + 1;
		MPI_Cart_rank(cartesian_handle, target_coords, &target);
		send_right_boundary(U, column, imax, jmax, target);
	}
	//receive left boundary U
	if (coords[1] >=1) {
		source_coords[0] = coords[0];
		source_coords[1] = coords[1] - 1;
		MPI_Cart_rank(cartesian_handle, source_coords, &source);
		receive_left_boundary_and_set_matrix(U, column, imax, jmax, source);
	}
	//send right boundary V
	if (coords[1] <= dims[1] - 2) {
		send_right_boundary(V, column, imax, jmax, target);
	}
	//receive left boundary V
	if (coords[1] >= 1) {
		receive_left_boundary_and_set_matrix(V, column, imax, jmax, source);
	}

	//send upper boundary U
	if (coords[0] >= 1) {
		target_coords[0] = coords[0] - 1;
		target_coords[1] = coords[1];
		MPI_Cart_rank(cartesian_handle, target_coords, &target);
		send_upper_boundary(U, row, imax, jmax, target);
	}
	//receive lower boundary U
	if (coords[0] <= dims[0] - 2) {
		source_coords[0] = coords[0] + 1;
		source_coords[1] = coords[1];
		MPI_Cart_rank(cartesian_handle, source_coords, &source);
		receive_lower_boundary_and_set_matrix(U, row, imax, jmax, source);
	}
	//send upper boundary V
	if (coords[0] >= 1) {
		send_upper_boundary(V, row, imax, jmax, target);
	}
	//receive lower boundary V
	if (coords[0] <= dims[0] - 2) {
		receive_lower_boundary_and_set_matrix(V, row, imax, jmax, source);
	}

	//send lower boundary U
	if (coords[0] <= dims[0] - 2) {
		target_coords[0] = coords[0] + 1;
		target_coords[1] = coords[1];
		MPI_Cart_rank(cartesian_handle, target_coords, &target);
		send_lower_boundary(U, row, imax, jmax, target);
	}
	//receive upper boundary U
	if (coords[0] >=1) {
		source_coords[0] = coords[0] - 1;
		source_coords[1] = coords[1];
		MPI_Cart_rank(cartesian_handle, source_coords, &source);
		receive_upper_boundary_and_set_matrix(U, row, imax, jmax, source);
	}
	//send lower boundary V
	if (coords[0] <= dims[0] - 2) {
		send_lower_boundary(V, row, imax, jmax, target);
	}
	//receive upper boundary V
	if (coords[0] >= 1) {
		receive_upper_boundary_and_set_matrix(V, row, imax, jmax, source);
	}
}


//applies boundary conditions to U and V on the boundaries of sub_matrices that correspond to the boundaries
//of the whole matrices
template <typename T>
void apply_boundary_conditions_UV_on_outer_processes(std::unique_ptr<std::unique_ptr<T[]>[]>& U, std::unique_ptr<std::unique_ptr<T[]>[]>& V, int imax, int jmax, int* coords, int* dims) {
	//change upper boundaries
	if (coords[0] == 0 ) {
		for (int j = 1; j < jmax + 1; j++) {
			U[0][j] = 0;
			V[0][j] = -V[1][j];
		}
	}
	//change lower boundaries
	if (coords[0] == dims[0] - 1) {
		for (int j = 1; j < jmax + 1; j++) {
			U[imax][j] = 0;
			V[imax + 1][j] = -V[imax][j];
		}
	}
	//change left boundaries
	if (coords[1] == 0 ) {
		for (int i = 1; i < imax + 1; i++) {
			U[i][0] = -U[i][1];
			V[i][0] = 0;
		}
	}
	//change right boundaries
	if (coords[1] == dims[1]-1) {
		for (int i = 1; i < imax + 1; i++) {
			U[i][jmax + 1] = 2-U[i][jmax];
			V[i][jmax] = 0;
		}
	}
}

//applies boundary conditions to the Pressure matrix on the boundaries of sub_matrices that correspond to the boundaries
//of the whole matrices
template <typename T>
void apply_boundary_conditions_Pressure_on_outer_processes(std::unique_ptr<std::unique_ptr<T[]>[]>& P, int imax, int jmax, int* coords, int* dims) {
	//change upper boundaries
	if (coords[0] == 0 && coords[1] < dims[1]) {
		for (int j = 1; j < jmax + 1; j++) {
			P[0][j] = P[1][j];
		}
	}
	//change lower boundaries
	if (coords[0] == dims[0] - 1 && coords[1] < dims[1]) {
		for (int j = 1; j < jmax + 1; j++) {
			P[imax + 1][j] = P[imax][j];
		}
	}
	//change left boundaries
	if (coords[1] == 0 && coords[0] < dims[0]) {
		for (int i = 1; i < imax + 1; i++) {
			P[i][0] = P[i][1];
		}
	}
	//change right boundaries
	if (coords[1] == dims[1] - 1 && coords[0] < dims[0]) {
		for (int i = 1; i < imax + 1; i++) {
			P[i][jmax + 1] = P[i][jmax];
		}
	}
}

//calculates delt and sets delt to that value in every process
template <typename T>
void calculate_delt_parallel(MPI_Comm cartesian_handle, std::unique_ptr<std::unique_ptr<T[]>[]>& U, std::unique_ptr<std::unique_ptr<T[]>[]>& V, int imax, int jmax, T tau, int Re, T delx, T dely, int myrank, T& delt) {
	if (tau < 0) return;
	T umax_local = find_maximal_absolute(imax + 2, jmax + 2, U);
	T vmax_local = find_maximal_absolute(imax + 2, jmax + 2, V);
	T umax = 0;
	T vmax = 0;
	MPI_Reduce(&umax_local, &umax, 1, get_MPI_Datatype(umax), MPI_MAX, 0, cartesian_handle);
	MPI_Reduce(&vmax_local, &vmax, 1, get_MPI_Datatype(vmax), MPI_MAX, 0, cartesian_handle);
	if (myrank == 0) {
		if (umax == 0) {
			if (vmax == 0) delt = tau * (Re / 2.0)*(1.0 / ((1.0 / (delx*delx)) + (1.0 / (dely*dely))));
			else delt = tau * min_2<T>((Re / 2.0)*(1.0 / ((1.0 / (delx*delx)) + (1.0 / (dely*dely)))), dely / vmax);
		}
		else {
			if (vmax == 0) delt = tau * min_2<T>((Re / 2.0)*(1.0 / ((1.0 / (delx*delx)) + (1.0 / (dely*dely)))), delx / umax);
			else delt = tau * min_3<T>((Re / 2.0)*(1.0 / ((1.0 / (delx*delx)) + (1.0 / (dely*dely)))), delx / umax, dely / vmax);
		}
	}
	MPI_Bcast(&delt, 1, get_MPI_Datatype(umax), 0, cartesian_handle);
}

//calculate F and G according to formula. 
template <typename T>
void calculate_F_G_parallel(int* coords, int* dims, int imax, int jmax, T delt, T delx, T dely, int Re, T alpha, T GX, T GY, std::unique_ptr<std::unique_ptr<T[]>[]>& U, std::unique_ptr<std::unique_ptr<T[]>[]>& V, std::unique_ptr<std::unique_ptr<T[]>[]>& F, std::unique_ptr<std::unique_ptr<T[]>[]>& G) {
	//apply boundary conditions first
	//change upper boundaries
	if (coords[0] == 0 && coords[1] < dims[1]) {
		for (int j = 1; j < jmax + 1; j++) {
			F[0][j] = U[0][j];
		}
	}
	//change lower boundaries
	if (coords[0] == dims[0] - 1 && coords[1] < dims[1]) {
		for (int j = 1; j < jmax + 1; j++) {
			F[imax][j] = U[imax][j];
		}
	}
	//change left boundaries
	if (coords[1] == 0 && coords[0] < dims[0]) {
		for (int i = 1; i < imax + 1; i++) {
			G[i][0] = V[i][0];
		}
	}
	//change right boundaries
	if (coords[1] == dims[1] - 1 && coords[0] < dims[0]) {
		for (int i = 1; i < imax + 1; i++) {
			G[i][jmax] = V[i][jmax];
		}
	}
	
	//now the rest (be aware, that the last row or column should only not be calculated, if the process 
	//contains the last row or column of global matrix)
	if (coords[0] == dims[0] - 1) {
		for (int i = 1; i < imax; i++) {
			for (int j = 1; j < jmax + 1; j++) {
				T ddudxx = (U[i + 1][j] - 2 * U[i][j] + U[i - 1][j]) / (delx*delx);
				T ddudyy = (U[i][j + 1] - 2 * U[i][j] + U[i][j - 1]) / (dely*dely);
				T duudx = (1 / delx)*((((U[i][j] + U[i + 1][j]) / 2)*((U[i][j] + U[i + 1][j]) / 2)) - (((U[i - 1][j] + U[i][j]) / 2)*((U[i - 1][j] + U[i][j]) / 2))) + alpha * (1 / (4 * delx))*((abs(U[i][j] + U[i + 1][j])*(U[i][j] - U[i + 1][j])) - (abs(U[i - 1][j] + U[i][j])*(U[i - 1][j] - U[i][j])));
				T duvdy = (1 / (dely * 4))*((V[i][j] + V[i + 1][j])*(U[i][j] + U[i][j + 1]) - (V[i][j - 1] + V[i + 1][j - 1])*(U[i][j - 1] + U[i][j])) + alpha * (1 / (dely * 4))*((abs(V[i][j] + V[i + 1][j])*(U[i][j] - U[i][j + 1])) - (abs(V[i][j - 1] + V[i + 1][j - 1])*(U[i][j - 1] - U[i][j])));
				F[i][j] = U[i][j] + delt * (((1.0 / Re)*(ddudxx + ddudyy)) - duudx - duvdy + GX);
			}
		}
	}
	else {
		for (int i = 1; i < imax+1; i++) {
			for (int j = 1; j < jmax + 1; j++) {
				T ddudxx = (U[i + 1][j] - 2 * U[i][j] + U[i - 1][j]) / (delx*delx);
				T ddudyy = (U[i][j + 1] - 2 * U[i][j] + U[i][j - 1]) / (dely*dely);
				T duudx = (1 / delx)*((((U[i][j] + U[i + 1][j]) / 2)*((U[i][j] + U[i + 1][j]) / 2)) - (((U[i - 1][j] + U[i][j]) / 2)*((U[i - 1][j] + U[i][j]) / 2))) + alpha * (1 / (4 * delx))*((abs(U[i][j] + U[i + 1][j])*(U[i][j] - U[i + 1][j])) - (abs(U[i - 1][j] + U[i][j])*(U[i - 1][j] - U[i][j])));
				T duvdy = (1 / (dely * 4))*((V[i][j] + V[i + 1][j])*(U[i][j] + U[i][j + 1]) - (V[i][j - 1] + V[i + 1][j - 1])*(U[i][j - 1] + U[i][j])) + alpha * (1 / (dely * 4))*((abs(V[i][j] + V[i + 1][j])*(U[i][j] - U[i][j + 1])) - (abs(V[i][j - 1] + V[i + 1][j - 1])*(U[i][j - 1] - U[i][j])));
				F[i][j] = U[i][j] + delt * (((1.0 / Re)*(ddudxx + ddudyy)) - duudx - duvdy + GX);
			}
		}
	}
	if (coords[1] == dims[1] - 1) {
		for (int i = 1; i < imax + 1; i++) {
			for (int j = 1; j < jmax; j++) {
				T duvdx = (1 / (delx * 4))*((U[i][j] + U[i][j + 1])*(V[i][j] + V[i + 1][j]) - (U[i - 1][j] + U[i - 1][j + 1])*(V[i - 1][j] + V[i][j])) + alpha * (1 / (delx * 4))*((abs(U[i][j] + U[i][j + 1])*(V[i][j] - V[i + 1][j])) - (abs(U[i - 1][j] + U[i - 1][j + 1])*(V[i - 1][j] - V[i][j])));
				T dvvdy = (1 / dely)*((((V[i][j] + V[i][j + 1]) / 2)*((V[i][j] + V[i][j + 1]) / 2)) - (((V[i][j - 1] + V[i][j]) / 2)*((V[i][j - 1] + V[i][j]) / 2))) + alpha * (1 / (4 * dely))*((abs(V[i][j] + V[i][j + 1])*(V[i][j] - V[i][j + 1])) - (abs(V[i][j - 1] + V[i][j])*(V[i][j - 1] - V[i][j])));
				T ddvdxx = (V[i + 1][j] - 2 * V[i][j] + V[i - 1][j]) / (delx*delx);
				T ddvdyy = (V[i][j + 1] - 2 * V[i][j] + V[i][j - 1]) / (dely*dely);
				G[i][j] = V[i][j] + delt * (((1.0 / Re)*(ddvdxx + ddvdyy)) - duvdx - dvvdy + GY);
			}
		}
	}
	else {
		for (int i = 1; i < imax + 1; i++) {
			for (int j = 1; j < jmax+1; j++) {
				T duvdx = (1 / (delx * 4))*((U[i][j] + U[i][j + 1])*(V[i][j] + V[i + 1][j]) - (U[i - 1][j] + U[i - 1][j + 1])*(V[i - 1][j] + V[i][j])) + alpha * (1 / (delx * 4))*((abs(U[i][j] + U[i][j + 1])*(V[i][j] - V[i + 1][j])) - (abs(U[i - 1][j] + U[i - 1][j + 1])*(V[i - 1][j] - V[i][j])));
				T dvvdy = (1 / dely)*((((V[i][j] + V[i][j + 1]) / 2)*((V[i][j] + V[i][j + 1]) / 2)) - (((V[i][j - 1] + V[i][j]) / 2)*((V[i][j - 1] + V[i][j]) / 2))) + alpha * (1 / (4 * dely))*((abs(V[i][j] + V[i][j + 1])*(V[i][j] - V[i][j + 1])) - (abs(V[i][j - 1] + V[i][j])*(V[i][j - 1] - V[i][j])));
				T ddvdxx = (V[i + 1][j] - 2 * V[i][j] + V[i - 1][j]) / (delx*delx);
				T ddvdyy = (V[i][j + 1] - 2 * V[i][j] + V[i][j - 1]) / (dely*dely);
				G[i][j] = V[i][j] + delt * (((1.0 / Re)*(ddvdxx + ddvdyy)) - duvdx - dvvdy + GY);
			}
		}
	}
}

//calculates the square of the global residual andsets it to that value in every process
template <typename T>
void calculate_residual_squared_parallel(MPI_Comm cartesian_handle, std::unique_ptr<std::unique_ptr<T[]>[]>& P, std::unique_ptr<std::unique_ptr<T[]>[]>& RHS, T& res_squared, int imax, int jmax, T delx, T dely, int myrank) {
	T sum = 0;
	for (int i = 1; i < imax + 1; i++) {
		for (int j = 1; j < jmax + 1; j++) {
			sum += ((P[i + 1][j] - 2 * P[i][j] + P[i - 1][j]) / (delx*delx) + (P[i][j + 1] - 2 * P[i][j] + P[i][j - 1]) / (dely*dely) - RHS[i][j])*((P[i + 1][j] - 2 * P[i][j] + P[i - 1][j]) / (delx*delx) + (P[i][j + 1] - 2 * P[i][j] + P[i][j - 1]) / (dely*dely) - RHS[i][j]);
		}
	}
	MPI_Reduce(&sum, &res_squared, 1, get_MPI_Datatype(sum), MPI_SUM, 0, cartesian_handle);
	if (myrank == 0) res_squared = res_squared / (imax*jmax);
	MPI_Bcast(&res_squared, 1, get_MPI_Datatype(sum), 0, cartesian_handle);
}

template <typename T>
void exchange_boundary_values_Pressure_between_processes(MPI_Comm cartesian_handle, std::unique_ptr<std::unique_ptr<T[]>[]>& P, T* column, T* row, int imax, int jmax, int*coords, int* dims) {
	int target = 0;
	int target_coords[2] = { 0,0 };
	int source = 0;
	int source_coords[2] = { 0,0 };

	//send left boundary P
	if (coords[1] >= 1) {
		target_coords[0] = coords[0];
		target_coords[1] = coords[1] - 1;
		MPI_Cart_rank(cartesian_handle, target_coords, &target);
		send_left_boundary(P, column, imax, jmax, target);
	}
	//receive right boundary P
	if (coords[1] <= dims[1] - 2) {
		source_coords[0] = coords[0];
		source_coords[1] = coords[1] + 1;
		MPI_Cart_rank(cartesian_handle, source_coords, &source);
		receive_right_boundary_and_set_matrix(P, column, imax, jmax, source);
	}

	//send right boundary P
	if (coords[1] <= dims[1] - 2) {
		target_coords[0] = coords[0];
		target_coords[1] = coords[1] + 1;
		MPI_Cart_rank(cartesian_handle, target_coords, &target);
		send_right_boundary(P, column, imax, jmax, target);
	}
	//receive left boundary P
	if (coords[1] >= 1) {
		source_coords[0] = coords[0];
		source_coords[1] = coords[1] - 1;
		MPI_Cart_rank(cartesian_handle, source_coords, &source);
		receive_left_boundary_and_set_matrix(P, column, imax, jmax, source);
	}

	//send upper boundary P
	if (coords[0] >= 1) {
		target_coords[0] = coords[0] - 1;
		target_coords[1] = coords[1];
		MPI_Cart_rank(cartesian_handle, target_coords, &target);
		send_upper_boundary(P, row, imax, jmax, target);
	}
	//receive lower boundary P
	if (coords[0] <= dims[0] - 2) {
		source_coords[0] = coords[0] + 1;
		source_coords[1] = coords[1];
		MPI_Cart_rank(cartesian_handle, source_coords, &source);
		receive_lower_boundary_and_set_matrix(P, row, imax, jmax, source);
	}

	//send lower boundary P
	if (coords[0] <= dims[0] - 2) {
		target_coords[0] = coords[0] + 1;
		target_coords[1] = coords[1];
		MPI_Cart_rank(cartesian_handle, target_coords, &target);
		send_lower_boundary(P, row, imax, jmax, target);
	}
	//receive upper boundary P
	if (coords[0] >= 1) {
		source_coords[0] = coords[0] - 1;
		source_coords[1] = coords[1];
		MPI_Cart_rank(cartesian_handle, source_coords, &source);
		receive_upper_boundary_and_set_matrix(P, row, imax, jmax, source);
	}
}


//pass P2, so it doesnt have to be created each time function is called
template<typename T>
int calculate_Pressure_with_SOR_parallel(MPI_Comm cartesian_handle, std::unique_ptr<std::unique_ptr<T[]>[]>& P, std::unique_ptr<std::unique_ptr<T[]>[]>& P2, std::unique_ptr<std::unique_ptr<T[]>[]>& RHS, int imax, int jmax, int itermax, T delx, T dely, T eps, T omg, T* column, T* row, int* coords, int* dims, int myrank, T&res) {
	apply_boundary_conditions_Pressure_on_outer_processes(P, imax, jmax, coords, dims);
	make_B_to_a_copy_of_A(imax, jmax, P2, P);
	T eps_squared = (eps * eps);
	T res_squared = eps_squared + 1;

	//use P and P2 in turn for the new and the old version of P
	//this way copying the matrix in each iteration is avoided
	int it = 0;
	while (it<itermax && res_squared>eps_squared) {
		//if it%2==0, P will become the new version and P2 the old one
		if (it % 2 == 0) {
			//apply boundary conditions from P2 to P

			//change upper boundaries
			if (coords[0] == 0 && coords[1] < dims[1]) {
				for (int j = 1; j < jmax + 1; j++) {
					P[0][j] = P2[1][j];
				}
			}
			//change lower boundaries
			if (coords[0] == dims[0] - 1 && coords[1] < dims[1]) {
				for (int j = 1; j < jmax + 1; j++) {
					P[imax + 1][j] = P2[imax][j];
				}
			}
			//change left boundaries
			if (coords[1] == 0 && coords[0] < dims[0]) {
				for (int i = 1; i < imax + 1; i++) {
					P[i][0] = P2[i][1];
				}
			}
			//change right boundaries
			if (coords[1] == dims[1] - 1 && coords[0] < dims[0]) {
				for (int i = 1; i < imax + 1; i++) {
					P[i][jmax + 1] = P2[i][jmax];
				}
			}
			//do 1 SOR cycle
			for (int i = 1; i < imax + 1; i++) {
				for (int j = 1; j < jmax + 1; j++) {
					P[i][j] = (1 - omg)*P2[i][j] + (omg / (2 * ((1 / (delx*delx)) + (1 / (dely*dely)))))*(((P2[i + 1][j] + P[i - 1][j]) / (delx * delx)) + ((P2[i][j + 1] + P[i][j - 1]) / (dely * dely)) - RHS[i][j]);
				}
			}

			exchange_boundary_values_Pressure_between_processes(cartesian_handle, P, column, row, imax, jmax, coords, dims);

			//calculate the squared norm of the residual of the now new P 
			calculate_residual_squared_parallel(cartesian_handle, P, RHS, res_squared, imax, jmax, delx, dely, myrank);
			it++;
		}
		//if it%2!=0, P2 woll become the new version and P is the old one
		else {
			//apply boundary conditions from P to P2
			//change upper boundaries
			if (coords[0] == 0 && coords[1] < dims[1]) {
				for (int j = 1; j < jmax + 1; j++) {
					P2[0][j] = P[1][j];
				}
			}
			//change lower boundaries
			if (coords[0] == dims[0] - 1 && coords[1] < dims[1]) {
				for (int j = 1; j < jmax + 1; j++) {
					P2[imax + 1][j] = P[imax][j];
				}
			}
			//change left boundaries
			if (coords[1] == 0 && coords[0] < dims[0]) {
				for (int i = 1; i < imax + 1; i++) {
					P2[i][0] = P[i][1];
				}
			}
			//change right boundaries
			if (coords[1] == dims[1] - 1 && coords[0] < dims[0]) {
				for (int i = 1; i < imax + 1; i++) {
					P2[i][jmax + 1] = P[i][jmax];
				}
			}

			//do 1 SOR cycle
			for (int i = 1; i < imax + 1; i++) {
				for (int j = 1; j < jmax + 1; j++) {
					P2[i][j] = (1 - omg)*P[i][j] + (omg / (2 * ((1 / (delx*delx)) + (1 / (dely*dely)))))*(((P[i + 1][j] + P2[i - 1][j]) / (delx * delx)) + ((P[i][j + 1] + P2[i][j - 1]) / (dely * dely)) - RHS[i][j]);
				}
			}

			exchange_boundary_values_Pressure_between_processes(cartesian_handle, P2, column, row, imax, jmax, coords, dims);

			//calculate the squared norm of the residual of the now new P2
			calculate_residual_squared_parallel(cartesian_handle, P2, RHS, res_squared, imax, jmax, delx, dely, myrank);
			it++;
		}
	}
	//if it%2==0, then P is the old version and has to be overwritten with the values from P2
	if (it % 2 == 0) make_B_to_a_copy_of_A(imax, jmax, P, P2);
	res = std::sqrt(res_squared);
	return it;
}


//calculate velocities U and V according to(10) and (11) 
//(be aware, that the last row or column should only not be calculated, if the process 
//contains the last row or column of global matrix)
template <typename T>
void calculate_U_and_V_parallel(int* coords, int* dims, int imax, int jmax, T delt, T delx, T dely, std::unique_ptr<std::unique_ptr<T[]>[]>& F, std::unique_ptr<std::unique_ptr<T[]>[]>& G, std::unique_ptr<std::unique_ptr<T[]>[]>& P, std::unique_ptr<std::unique_ptr<T[]>[]>& U, std::unique_ptr<std::unique_ptr<T[]>[]>& V) {
	if (coords[0] == dims[0] - 1) {
		for (int i = 1; i < imax; i++) {
			for (int j = 1; j < jmax + 1; j++) {
				U[i][j] = F[i][j] - (delt / delx)*(P[i + 1][j] - P[i][j]);
			}
		}
	}
	else {
		for (int i = 1; i < imax+1; i++) {
			for (int j = 1; j < jmax + 1; j++) {
				U[i][j] = F[i][j] - (delt / delx)*(P[i + 1][j] - P[i][j]);
			}
		}
	}
	if (coords[1] == dims[1] - 1) {
		for (int i = 1; i < imax + 1; i++) {
			for (int j = 1; j < jmax; j++) {
				V[i][j] = G[i][j] - (delt / dely)*(P[i][j + 1] - P[i][j]);
			}
		}
	}
	else {
		for (int i = 1; i < imax + 1; i++) {
			for (int j = 1; j < jmax + 1; j++) {
				V[i][j] = G[i][j] - (delt / dely)*(P[i][j + 1] - P[i][j]);
			}
		}
	}
}


template<typename T>
void write_output_into_file_parallel(std::string filename, MPI_Comm cartesian_handle, int myrank, int* coords, int* dims, T* row, std::unique_ptr<std::unique_ptr<T[]>[]>& U, std::unique_ptr<std::unique_ptr<T[]>[]>& V, std::unique_ptr<std::unique_ptr<T[]>[]>& P, std::unique_ptr<std::unique_ptr<T[]>[]>& U2, std::unique_ptr<std::unique_ptr<T[]>[]>& V2, int imax, int jmax, int imax_global, int jmax_global, T xlength, T ylength) {
	std::ofstream myfile;
	T a=0;
	int source = 0;
	int coords_source[2] = { 0,0 };
	MPI_Status status;

	//colculate U2 and V2 as the values of U and V in the middle of cells
	for (int i = 1; i < imax + 1; i++) {
		for (int j = 1; j < jmax + 1; j++) {
			U2[i][j] = (U[i][j - 1] + U[i][j]) / 2;
			V2[i][j] = (V[i][j] + V[i + 1][j]) / 2;
		}
	}
	int rank_top_left;
	
	int coords_top_left[2] = { 0,0 };
	MPI_Cart_rank(cartesian_handle, coords_top_left, &rank_top_left);
	if (myrank == rank_top_left) {
		myfile.open(filename);
		if (!myfile.is_open()) {
			std::cout << "ERROR: write_output_into_file failed to open file";
			return;
		}
		myfile << xlength << "\n" << ylength << "\n" << imax_global << "\n" << jmax_global << "\n";
	}
	//write U2 into file (first the rows that process top left contains a part of)
	if (coords[0] == 0) {
		for (int i = 1; i < imax + 1; i++) {
			if (myrank == rank_top_left) {
				for (int j = 1; j < jmax + 1; j++) myfile << U2[i][j] << " ";
			}
			//send row i to top left process 
			for (int column_process_matrix = 1; column_process_matrix < dims[1]; column_process_matrix++) {
				if (coords[1]==column_process_matrix) {
					for (int j = 0; j < jmax; j++) row[j] = U2[i][j + 1];
					MPI_Send(row, jmax, get_MPI_Datatype(a), rank_top_left, 110, cartesian_handle);
				}
				if (myrank == rank_top_left) {
					coords_source[1] = column_process_matrix;
					MPI_Cart_rank(cartesian_handle, coords_source, &source);
					MPI_Recv(row, jmax, get_MPI_Datatype(a), source, 110, cartesian_handle, &status);
					for (int j = 0; j < jmax; j++) myfile << row[j] << " ";
				}
			}
			if (myrank == rank_top_left) myfile << "\n";
		}
	}
	//now write remaining rows into file
	for (int row_process_matrix = 1; row_process_matrix < dims[0]; row_process_matrix++) {
		
		for (int i = 1; i < imax + 1; i++) {
			//send row i to top left process 
			for (int column_process_matrix = 0; column_process_matrix < dims[1]; column_process_matrix++) {
				if (coords[1] == column_process_matrix && coords[0]==row_process_matrix) {
					for (int j = 0; j < jmax; j++) row[j] = U2[i][j + 1];
					MPI_Send(row, jmax, get_MPI_Datatype(a), rank_top_left, 110, cartesian_handle);
				}
				if (myrank == rank_top_left) {
					coords_source[1] = column_process_matrix;
					coords_source[0] = row_process_matrix;
					MPI_Cart_rank(cartesian_handle, coords_source, &source);
					MPI_Recv(row, jmax, get_MPI_Datatype(a), source, 110, cartesian_handle, &status);
					for (int j = 0; j < jmax; j++) myfile << row[j] << " ";
				}
			}
			if (myrank == rank_top_left) myfile << "\n";
		}
	}

	//write V2 into file (first the rows that process top left contains a part of)
	if (coords[0] == 0) {
		coords_source[0] = 0;
		for (int i = 1; i < imax + 1; i++) {
			if (myrank == rank_top_left) {
				for (int j = 1; j < jmax + 1; j++) myfile << V2[i][j] << " ";
			}
			//send row i to top left process 
			for (int column_process_matrix = 1; column_process_matrix < dims[1]; column_process_matrix++) {
				if (coords[1] == column_process_matrix) {
					for (int j = 0; j < jmax; j++) row[j] = V2[i][j + 1];
					MPI_Send(row, jmax, get_MPI_Datatype(a), rank_top_left, 110, cartesian_handle);
				}
				if (myrank == rank_top_left) {
					coords_source[1] = column_process_matrix;
					MPI_Cart_rank(cartesian_handle, coords_source, &source);
					MPI_Recv(row, jmax, get_MPI_Datatype(a), source, 110, cartesian_handle, &status);
					for (int j = 0; j < jmax; j++) myfile << row[j] << " ";
				}
			}
			if (myrank == rank_top_left) myfile << "\n";
		}
	}
	//now write remaining rows into file
	for (int row_process_matrix = 1; row_process_matrix < dims[0]; row_process_matrix++) {

		for (int i = 1; i < imax + 1; i++) {
			//send row i to top left process 
			for (int column_process_matrix = 0; column_process_matrix < dims[1]; column_process_matrix++) {
				if (coords[1] == column_process_matrix && coords[0] == row_process_matrix) {
					for (int j = 0; j < jmax; j++) row[j] = V2[i][j + 1];
					MPI_Send(row, jmax, get_MPI_Datatype(a), rank_top_left, 110, cartesian_handle);
				}
				if (myrank == rank_top_left) {
					coords_source[1] = column_process_matrix;
					coords_source[0] = row_process_matrix;
					MPI_Cart_rank(cartesian_handle, coords_source, &source);
					MPI_Recv(row, jmax, get_MPI_Datatype(a), source, 110, cartesian_handle, &status);
					for (int j = 0; j < jmax; j++) myfile << row[j] << " ";
				}
			}
			if (myrank == rank_top_left) myfile << "\n";
		}
	}

	//write P into file (first the rows that process top left contains a part of)
	if (coords[0] == 0) {
		coords_source[0] = 0;
		for (int i = 1; i < imax + 1; i++) {
			if (myrank == rank_top_left) {
				for (int j = 1; j < jmax + 1; j++) myfile << P[i][j] << " ";
			}
			//send row i to top left process 
			for (int column_process_matrix = 1; column_process_matrix < dims[1]; column_process_matrix++) {
				if (coords[1] == column_process_matrix) {
					for (int j = 0; j < jmax; j++) row[j] = P[i][j + 1];
					MPI_Send(row, jmax, get_MPI_Datatype(a), rank_top_left, 110, cartesian_handle);
				}
				if (myrank == rank_top_left) {
					coords_source[1] = column_process_matrix;
					MPI_Cart_rank(cartesian_handle, coords_source, &source);
					MPI_Recv(row, jmax, get_MPI_Datatype(a), source, 110, cartesian_handle, &status);
					for (int j = 0; j < jmax; j++) myfile << row[j] << " ";
				}
			}
			if (myrank == rank_top_left) myfile << "\n";
		}
	}
	//now write remaining rows into file
	for (int row_process_matrix = 1; row_process_matrix < dims[0]; row_process_matrix++) {

		for (int i = 1; i < imax + 1; i++) {
			//send row i to top left process 
			for (int column_process_matrix = 0; column_process_matrix < dims[1]; column_process_matrix++) {
				if (coords[1] == column_process_matrix && coords[0] == row_process_matrix) {
					for (int j = 0; j < jmax; j++) row[j] = P[i][j + 1];
					MPI_Send(row, jmax, get_MPI_Datatype(a), rank_top_left, 110, cartesian_handle);
				}
				if (myrank == rank_top_left) {
					coords_source[1] = column_process_matrix;
					coords_source[0] = row_process_matrix;
					MPI_Cart_rank(cartesian_handle, coords_source, &source);
					MPI_Recv(row, jmax, get_MPI_Datatype(a), source, 110, cartesian_handle, &status);
					for (int j = 0; j < jmax; j++) myfile << row[j] << " ";
				}
			}
			if (myrank == rank_top_left) myfile << "\n";
		}
	}
}
