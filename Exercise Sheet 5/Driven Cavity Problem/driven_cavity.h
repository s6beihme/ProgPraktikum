#pragma once
#include <iostream>
#include <memory>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cmath>
//In here I write all the functions that will be used by the final algorithm

template <typename T> 
T convert_to(const std::string &str)
{
	std::istringstream ss(str);
	T num;
	ss >> num;
	return num;
}

template<typename T>
T use_partition_to_convert_to_T(std::string s) {
	std::size_t found = s.find("=");
	s.erase(0, found + 1);
	return convert_to<T>(s);
}

template <typename T>
void print_matrix(int row_count, int column_count, std::unique_ptr<std::unique_ptr<T[]>[]>& U) {
	for (int i = 0; i < row_count; i++) {
		for (int j = 0; j < column_count; j++) {
			std::cout << U[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

template<typename T>
void read_parameters_from_file(std::string filename, int& imax, int& jmax, T& xlength, T& ylength, T& delt, T& t_end, T& tau, T& del_vec, T& eps, T& omg, T& alpha, int& itermax, T& GX, T& GY, T& Re, T& UI, T& VI, T& PI){
	std::ifstream myfile;
	myfile.open(filename);

	if (!myfile.is_open()) { //open vs is_open?
		std::cout << "FILE DIDNT OPEN!\n";
		exit(0);
	}
	std::string temp;

	myfile >> temp;
	imax = use_partition_to_convert_to_T<int>(temp);
	myfile >> temp;
	jmax = use_partition_to_convert_to_T<int>(temp);
	myfile >> temp;
	xlength = use_partition_to_convert_to_T<T>(temp);
	myfile >> temp;
	ylength = use_partition_to_convert_to_T<T>(temp);
	myfile >> temp;
	delt = use_partition_to_convert_to_T<T>(temp);
	myfile >> temp;
	t_end = use_partition_to_convert_to_T<T>(temp);
	myfile >> temp;
	tau = use_partition_to_convert_to_T<T>(temp);
	myfile >> temp;
	del_vec = use_partition_to_convert_to_T<T>(temp);
	myfile >> temp;
	eps = use_partition_to_convert_to_T<T>(temp);
	myfile >> temp;
	omg = use_partition_to_convert_to_T<T>(temp);
	myfile >> temp;
	alpha = use_partition_to_convert_to_T<T>(temp);
	myfile >> temp;
	itermax = use_partition_to_convert_to_T<int>(temp);
	myfile >> temp;
	GX = use_partition_to_convert_to_T<T>(temp);
	myfile >> temp;
	GY = use_partition_to_convert_to_T<T>(temp);
	myfile >> temp;
	Re = use_partition_to_convert_to_T<int>(temp);
	myfile >> temp;
	UI = use_partition_to_convert_to_T<T>(temp);
	myfile >> temp;
	VI = use_partition_to_convert_to_T<T>(temp);
	myfile >> temp;
	PI = use_partition_to_convert_to_T<T>(temp);
}

//initializes U, V and P according to UI, VI, PI
//those matrices have to exist already and have right dimensions
template<typename T>
void initialize_fields_UVP(int imax, int jmax, T UI, T VI, T PI, std::unique_ptr<std::unique_ptr<T[]>[]>& U, std::unique_ptr<std::unique_ptr<T[]>[]>& V, std::unique_ptr<std::unique_ptr<T[]>[]>& P) {
	for (int i = 1; i < imax + 1; i++) {
		for (int j = 1; j < jmax + 1; j++) {
			U[i][j] = UI;
			V[i][j] = VI;
			P[i][j] = PI;
		}
	}
}

//applies boundary conditions(17, 18) to U and V
template<typename T>
void apply_boundary_conditions_UV_2(int imax, int jmax, std::unique_ptr<std::unique_ptr<T[]>[]>& U, std::unique_ptr<std::unique_ptr<T[]>[]>& V) {
	for (int j = 1; j < jmax + 1; j++) {
		U[0][j] = 0;
		U[imax][j] = 0;
		V[0][j] = -V[1][j];
		V[imax + 1][j] = -V[imax][j];
	}
	for (int i = 1; i < imax + 1; i++){
		U[i][0] = -U[i][1];
		U[i][jmax + 1] = -U[i][jmax];
		V[i][0] = 0;
		V[i][jmax] = 0;
	}
}

//applies boundary conditions(17, 18) to U and V
template <typename T>
void apply_boundary_conditions_UV(int imax, int jmax, std::unique_ptr<std::unique_ptr<T[]>[]>& U, std::unique_ptr<std::unique_ptr<T[]>[]>& V) {
	for (int j = 1; j < jmax + 1; j++) {
		U[0][j] = 0;
		U[imax][j] = 0;
		V[0][j] = -V[1][j];
		V[imax + 1][j] = -V[imax][j];
	}
	for (int i = 1; i < imax + 1; i++) {
		U[i][0] = -U[i][1];
		U[i][jmax + 1] = 2 - U[i][jmax];
		V[i][0] = 0;
		V[i][jmax] = 0;
	}
}


template<typename T>
T find_maximal_absolute(int row_count, int column_count, std::unique_ptr<std::unique_ptr<T[]>[]>& M) {
	T max = 0;
	for (int i =0 ; i < row_count; i++) {
		for (int j = 0; j < column_count; j++) {
			if (abs(M[i][j]) > max) max = abs(M[i][j]);
		}
	}
	return max;
}

//need to define min, max

template<typename T>
T min_2(T a, T b) {
	if (a <= b) return a;
	else return b;
}

template<typename T>
T min_3(T a, T b, T c) {
	if (a <= b && a <= c) return a;
	else {
		if (b <= c) return b;
		else return c;
	}
}

template<typename T>
void calculate_delt(int imax,int jmax, T tau,int Re, T delx,T dely, std::unique_ptr<std::unique_ptr<T[]>[]>& U, std::unique_ptr<std::unique_ptr<T[]>[]>& V, T& delt) {
	if (tau < 0) return;
	T umax = find_maximal_absolute(imax + 2, jmax + 2, U);
	T vmax = find_maximal_absolute(imax + 2, jmax + 2, V);
	if (umax == 0) {
		if (vmax == 0) delt = tau * (Re / 2.0)*(1.0 / ((1.0 / (delx*delx)) + (1.0 / (dely*dely)))); 
		else delt= tau * min_2<T>((Re / 2.0)*(1.0 / ((1.0 / (delx*delx)) + (1.0 / (dely*dely)))), dely / vmax);
	}
	else {
		if (vmax == 0) delt= tau * min_2<T>((Re / 2.0)*(1.0 / ((1.0 / (delx*delx)) + (1.0 / (dely*dely)))), delx / umax);
		else delt= tau * min_3<T>((Re / 2.0)*(1.0 / ((1.0 / (delx*delx)) + (1.0 / (dely*dely)))), delx / umax, dely / vmax);
	}
}

template <typename T>
void calculate_F_G(int imax,int jmax,T delt,T delx,T dely, T Re, T alpha,T GX,T GY, std::unique_ptr<std::unique_ptr<T[]>[]>& U, std::unique_ptr<std::unique_ptr<T[]>[]>& V, std::unique_ptr<std::unique_ptr<T[]>[]>& F, std::unique_ptr<std::unique_ptr<T[]>[]>& G) {
	//apply boundary conditions first
	for(int j=1; j<jmax+1; j++) {
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}
	for(int i=1; i<imax+1; i++) {
		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
	}
	//now the rest
	for(int i=1; i<imax; i++) {
		for(int j=1; j<jmax+1; j++) {
			T ddudxx = (U[i + 1][j] - 2 * U[i][j] + U[i - 1][j]) / (delx*delx);
			T ddudyy = (U[i][j + 1] - 2 * U[i][j] + U[i][j - 1]) / (dely*dely);
			T duudx = (1 / delx)*((((U[i][j] + U[i + 1][j]) / 2)*((U[i][j] + U[i + 1][j]) / 2)) - (((U[i - 1][j] + U[i][j]) / 2)*((U[i - 1][j] + U[i][j]) / 2))) + alpha * (1 / (4 * delx))*((abs(U[i][j] + U[i + 1][j])*(U[i][j] - U[i + 1][j])) - (abs(U[i - 1][j] + U[i][j])*(U[i - 1][j] - U[i][j])));
			T duvdy = (1 / (dely * 4))*((V[i][j] + V[i + 1][j])*(U[i][j] + U[i][j + 1]) - (V[i][j - 1] + V[i + 1][j - 1])*(U[i][j - 1] + U[i][j])) + alpha * (1 / (dely * 4))*((abs(V[i][j] + V[i + 1][j])*(U[i][j] - U[i][j + 1])) - (abs(V[i][j - 1] + V[i + 1][j - 1])*(U[i][j - 1] - U[i][j])));
			F[i][j] = U[i][j] + delt * (((1 / Re)*(ddudxx + ddudyy)) - duudx - duvdy + GX);
		}
	}
	for(int i=1; i<imax+1; i++) {
		for(int j=1; j<jmax; j++) {
			T duvdx = (1 / (delx * 4))*((U[i][j] + U[i][j + 1])*(V[i][j] + V[i + 1][j]) - (U[i - 1][j] + U[i - 1][j + 1])*(V[i - 1][j] + V[i][j])) + alpha * (1 / (delx * 4))*((abs(U[i][j] + U[i][j + 1])*(V[i][j] - V[i + 1][j])) - (abs(U[i - 1][j] + U[i - 1][j + 1])*(V[i - 1][j] - V[i][j])));
			T dvvdy = (1 / dely)*((((V[i][j] + V[i][j + 1]) / 2)*((V[i][j] + V[i][j + 1]) / 2)) - (((V[i][j - 1] + V[i][j]) / 2)*((V[i][j - 1] + V[i][j]) / 2))) + alpha * (1 / (4 * dely))*((abs(V[i][j] + V[i][j + 1])*(V[i][j] - V[i][j + 1])) - (abs(V[i][j - 1] + V[i][j])*(V[i][j - 1] - V[i][j])));
			T ddvdxx = (V[i + 1][j] - 2 * V[i][j] + V[i - 1][j]) / (delx*delx);
			T ddvdyy = (V[i][j + 1] - 2 * V[i][j] + V[i][j - 1]) / (dely*dely);
			G[i][j] = V[i][j] + delt * (((1 / Re)*(ddvdxx + ddvdyy)) - duvdx - dvvdy + GY);
		}
	}
}

template <typename T>
void calculate_RHS(int imax, int jmax, T delt, T delx, T dely, std::unique_ptr<std::unique_ptr<T[]>[]>& F, std::unique_ptr<std::unique_ptr<T[]>[]>& G, std::unique_ptr<std::unique_ptr<T[]>[]>& RHS) {
	for (int i = 1; i < imax + 1; i++) {
		for (int j = 1; j < jmax + 1; j++) {
			RHS[i][j] = (1 / delt)*(((F[i][j] - F[i - 1][j]) / delx) + ((G[i][j] - G[i][j - 1]) / dely));
		}
	}
}

template <typename T>
void apply_boundary_conditions_Pressure(int imax, int jmax, std::unique_ptr<std::unique_ptr<T[]>[]>& P) {
	for (int j = 1; j < jmax + 1; j++) {
		P[0][j] = P[1][j];
		P[imax + 1][j] = P[imax][j];
	}
	for (int i = 1; i < imax + 1; i++) {
		P[i][0] = P[i][1];
		P[i][jmax + 1] = P[i][jmax];
	}
}

template<typename T>
void make_B_to_a_copy_of_A(int imax, int jmax, std::unique_ptr<std::unique_ptr<T[]>[]>& B, std::unique_ptr<std::unique_ptr<T[]>[]>& A) {
	for (int i = 0; i < imax + 2; i++) {
		for (int j = 0; j < jmax + 2; j++) {
			B[i][j] = A[i][j];
		}
	}
}

template<typename T>
void calculate_residual_squared(int imax, int jmax, T delx, T dely, std::unique_ptr<std::unique_ptr<T[]>[]>& P, std::unique_ptr<std::unique_ptr<T[]>[]>& RHS, T& res_squared) {
	T sum = 0;
	for (int i = 1; i < imax + 1; i++) {
		for (int j = 1; j < jmax + 1; j++) {
			sum += ((P[i + 1][j] - 2 * P[i][j] + P[i - 1][j]) / (delx*delx) + (P[i][j + 1] - 2 * P[i][j] + P[i][j - 1]) / (dely*dely) - RHS[i][j])*((P[i + 1][j] - 2 * P[i][j] + P[i - 1][j]) / (delx*delx) + (P[i][j + 1] - 2 * P[i][j] + P[i][j - 1]) / (dely*dely) - RHS[i][j]);
		}
	}
	res_squared = sum / (imax*jmax);
}

//pass P_old, so it doesnt have to be created each time function is used
template<typename T>
int calculate_Pressure_with_SOR(int imax, int jmax, int itermax, T delx, T dely, T eps, T omg, std::unique_ptr<std::unique_ptr<T[]>[]>& RHS, std::unique_ptr<std::unique_ptr<T[]>[]>& P, std::unique_ptr<std::unique_ptr<T[]>[]>& P2, T& res) {
	apply_boundary_conditions_Pressure(imax, jmax, P);
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
			for (int j = 1; j < jmax + 1; j++) {
				P[0][j] = P2[1][j];
				P[imax + 1][j] = P2[imax][j];
			}
			for (int i = 1; i < imax + 1; i++) {
				P[i][0] = P2[i][1];
				P[i][jmax + 1] = P2[i][jmax];
			}

			//do 1 SOR cycle
			for (int i = 1; i < imax + 1; i++) {
				for (int j = 1; j < jmax + 1; j++) {
					P[i][j] = (1 - omg)*P2[i][j] + (omg / (2 * ((1 / (delx*delx)) + (1 / (dely*dely)))))*(((P2[i + 1][j] + P[i - 1][j]) / (delx * delx)) + ((P2[i][j + 1] + P[i][j - 1]) / (dely * dely)) - RHS[i][j]);
				}
			}
			//calculate the squared norm of the residual of the now new P 
			calculate_residual_squared(imax, jmax, delx, dely, P, RHS, res_squared);
			it++;
		}
		//if it%2!=0, P2 woll become the new version and P is the old one
		else {
			//apply boundary conditions from P to P2
			for (int j = 1; j < jmax + 1; j++) {
				P2[0][j] = P[1][j];
				P2[imax + 1][j] = P[imax][j];
			}
			for (int i = 1; i < imax + 1; i++) {
				P2[i][0] = P[i][1];
				P2[i][jmax + 1] = P[i][jmax];
			}

			//do 1 SOR cycle
			for (int i = 1; i < imax + 1; i++) {
				for (int j = 1; j < jmax + 1; j++) {
					P2[i][j] = (1 - omg)*P[i][j] + (omg / (2 * ((1 / (delx*delx)) + (1 / (dely*dely)))))*(((P[i + 1][j] + P2[i - 1][j]) / (delx * delx)) + ((P[i][j + 1] + P2[i][j - 1]) / (dely * dely)) - RHS[i][j]);
				}
			}
			//calculate the squared norm of the residual of the now new P2
			calculate_residual_squared(imax, jmax, delx, dely, P2, RHS, res_squared);
			it++;
		}
	}
	//if it%2==0, then P is the old version and has to be overwritten with the values from P2
	if (it % 2 == 0) make_B_to_a_copy_of_A(imax, jmax, P, P2);
	res = std::sqrt(res_squared);
	return it;
}



//calculate velocities U and V according to(10) and (11)
template <typename T>
void calculate_U_and_V(int imax, int jmax, T delt, T delx, T dely, std::unique_ptr<std::unique_ptr<T[]>[]>& F, std::unique_ptr<std::unique_ptr<T[]>[]>& G, std::unique_ptr<std::unique_ptr<T[]>[]>& P, std::unique_ptr<std::unique_ptr<T[]>[]>& U, std::unique_ptr<std::unique_ptr<T[]>[]>& V) {
	for (int i = 1; i < imax; i++) {
		for (int j = 1; j < jmax + 1; j++) {
			U[i][j] = F[i][j] - (delt / delx)*(P[i + 1][j] - P[i][j]);
		}
	}
	for (int i = 1; i < imax + 1; i++) {
		for (int j = 1; j < jmax; j++) {
			V[i][j] = G[i][j] - (delt / dely)*(P[i][j + 1] - P[i][j]);
		}
	}
}

//maybe write the python version for heuristic testing
/*template <typename T>
int calculate_Pressure_with_SOR_2(int imax, int jmax, int itermax, T delx, T dely, T eps, T omg, std::unique_ptr<std::unique_ptr<T[]>[]>& RHS, std::unique_ptr<std::unique_ptr<T[]>[]>& P, std::unique_ptr<std::unique_ptr<T[]>[]>& P_old, T& res) {
	apply_boundary_conditions_Pressure(imax, jmax, P);
	make_B_to_a_copy_of_A(imax, jmax, P_old, P);

	it = 0;
	res_squared = eps * eps + 1;
	while (it<itermax && res_squared>eps*eps) {
		apply_boundary_conditions_Pressure(imax, jmax, P);
		apply_boundary_conditions_Pressure(imax, jmax, P_old);
			//do 1 SOR cycle
		for (int i = 1; i < imax + 1; i++) {
			for (int j = 1; j < jmax + 1; j++) {
				P[i][j] = (1 - omg)*P_old[i][j] + (omg / (2 * ((1 / (delx*delx)) + (1 / (dely*dely)))))*(((P_old[i + 1][j] + P[i - 1][j]) / delx * *2) + ((P_old[i][j + 1] + P[i][j - 1]) / dely * *2) - RHS[i][j]);
			}
		}
				make_B_to_a_copy_of_A(imax, jmax, P_old, P)
				#calculate residuum(16)
				res_squared = sum([sum([(((P[i + 1][j] - 2 * P[i][j] + P[i - 1][j]) / delx * *2) + ((P[i][j + 1] - 2 * P[i][j] + P[i][j - 1]) / dely * *2) - RHS[i][j])**2 for j in range(1, jmax + 1)]) for i in range(1, imax + 1)]) / (imax*jmax)
				it += 1
# =============================================================================
#     if(res_squared>eps*eps):
				#         print("SOR didnt yield a result that was close enough")
# =============================================================================
				return it, sqrt(res_squared)
	}
}*/

//wites data into filename according to given format
//note that it doesnt write the boundary values into the file
template <typename T>
void write_output_into_file(std::string filename, T xlength, T ylength, int imax, int jmax, std::unique_ptr<std::unique_ptr<T[]>[]>& U, std::unique_ptr<std::unique_ptr<T[]>[]>& V, std::unique_ptr<std::unique_ptr<T[]>[]>& P, std::unique_ptr<std::unique_ptr<T[]>[]>& U2, std::unique_ptr<std::unique_ptr<T[]>[]>& V2) {
	//colculate U2 and V2 as the values of U and V in the middle of cells
	for (int i = 1; i < imax + 1; i++) {
		for (int j = 1; j < jmax + 1; j++) {
			U2[i][j] = (U[i][j - 1] + U[i][j]) / 2;
			V2[i][j] = (V[i][j] + V[i + 1][j]) / 2;
		}
	}

	//write data into file
	std::ofstream myfile;
	myfile.open(filename);
	if (!myfile.is_open()) {
		std::cout << "ERROR: write_output_into_file failed to open file";
		return;
	}
	myfile << xlength << "\n" << ylength << "\n" << imax << "\n" << jmax << "\n";
	for (int i = 1; i < imax + 1; i++) {
		for (int j = 1; j < jmax + 1; j++) {
			myfile << U2[i][j] << " ";
		}
		myfile << "\n";
	}
	for (int i = 1; i < imax + 1; i++) {
		for (int j = 1; j < jmax + 1; j++) {
			myfile << V2[i][j] << " ";
		}
		myfile << "\n";
	}
	for (int i = 1; i < imax + 1; i++) {
		for (int j = 1; j < jmax + 1; j++) {
			myfile << P[i][j] << " ";
		}
		myfile << "\n";
	}
}
