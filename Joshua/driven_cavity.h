#pragma once
#include <iostream>
#include <memory>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>

/*
T convert_to(const std::string &str)
T use_partition_to_convert_to_T(std::string s)
void input(std::string filename, int& imax, int& jmax, T& xlength, T& ylength, T& delt, T& t_end, T& tau, T& del_vec, T& eps, T& omg, T& alpha, int& itermax, T& GX, T& GY, T& Re, T& UI, T& VI, T& PI)
void output(std::string filename, T xlength, T ylength, int imax, int jmax, std::unique_ptr<std::unique_ptr<T[]>[]>& U, std::unique_ptr<std::unique_ptr<T[]>[]>& V, std::unique_ptr<std::unique_ptr<T[]>[]>& P, std::unique_ptr<std::unique_ptr<T[]>[]>& U2, std::unique_ptr<std::unique_ptr<T[]>[]>& V2)
 Mmax(std::unique_ptr<std::unique_ptr<T[]>[]> &Mat, int imax, int jmax)
 T _u2_x(std::unique_ptr<std::unique_ptr<T[]>[]> &U, int i, int j, T delx, T alpha)
 T _uv_y(std::unique_ptr<std::unique_ptr<T[]>[]> &U, std::unique_ptr<std::unique_ptr<T[]>[]> &V, int i, int j, T dely, T alpha)
 T _2u_x2(std::unique_ptr<std::unique_ptr<T[]>[]> &U, int i, int j, T delx)
 T _2u_y2(std::unique_ptr<std::unique_ptr<T[]>[]> &U, int i, int j, T dely)
 T _p_x(std::unique_ptr<std::unique_ptr<T[]>[]> &P, int i, int j, T delx)
 T _v2_y(std::unique_ptr<std::unique_ptr<T[]>[]> &V, int i, int j, T dely, T alpha)
 T _uv_x(std::unique_ptr<std::unique_ptr<T[]>[]> &U, std::unique_ptr<std::unique_ptr<T[]>[]> &V, int i, int j, T delx, T alpha)
 T _2v_x2(std::unique_ptr<std::unique_ptr<T[]>[]> &V, int i, int j, T delx)
 T _2v_y2(std::unique_ptr<std::unique_ptr<T[]>[]> &V, int i, int j, T dely)
 T _p_y(std::unique_ptr<std::unique_ptr<T[]>[]> &P, int i, int j, T dely)
 void init_fields(std::unique_ptr<std::unique_ptr<T[]>[]> &Mat, int imax, int jmax, T MI)
 T _delt(std::unique_ptr<std::unique_ptr<T[]>[]> &U, std::unique_ptr<std::unique_ptr<T[]>[]> &V, int imax, int jmax, T tau, T Re, T delt, T delx, T dely)
 void boundary_conditions_UV(std::unique_ptr<std::unique_ptr<T[]>[]> &U, std::unique_ptr<std::unique_ptr<T[]>[]> &V, int imax, int jmax)
 void compute_FG(std::unique_ptr<std::unique_ptr<T[]>[]> &F, std::unique_ptr<std::unique_ptr<T[]>[]> &G, std::unique_ptr<std::unique_ptr<T[]>[]> &U, std::unique_ptr<std::unique_ptr<T[]>[]> &V, int imax, int jmax, T delt, T delx, T dely, T Re, T GX, T GY, T alpha)
 void compute_RHS(std::unique_ptr<std::unique_ptr<T[]>[]> &F, std::unique_ptr<std::unique_ptr<T[]>[]> &G, std::unique_ptr<std::unique_ptr<T[]>[]> &RHS, int imax, int jmax, T delt, T delx, T dely)
 void SOR(std::unique_ptr<std::unique_ptr<T[]>[]> &RHS, std::unique_ptr<std::unique_ptr<T[]>[]> &P, int itermax, int imax, int jmax, T omg, T delx, T dely, T eps)
 void compute_UV(std::unique_ptr<std::unique_ptr<T[]>[]> &U, std::unique_ptr<std::unique_ptr<T[]>[]> &V, std::unique_ptr<std::unique_ptr<T[]>[]> &F, std::unique_ptr<std::unique_ptr<T[]>[]> &G, std::unique_ptr<std::unique_ptr<T[]>[]> &P, int imax, int jmax, T delt, T delx, T dely)
 */

template <typename T>
T convert_to(const std::string &str)
{
	std::istringstream ss(str);
	T num;
	ss >> num;
	return num;
}

template<typename T>
T use_partition_to_convert_to_T(std::string s) 
{
	std::size_t found = s.find("=");
	s.erase(0, found + 1);
	return convert_to<T>(s);
}

template<typename T>
void input(std::string filename, int& imax, int& jmax, T& xlength, T& ylength, T& delt, T& t_end, T& tau, T& del_vec, T& eps, T& omg, T& alpha, int& itermax, T& GX, T& GY, T& Re, T& UI, T& VI, T& PI) 
{
	std::ifstream myfile;
	myfile.open(filename);

	if (!myfile.is_open()) 
	{ 
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

template <typename T>
void output(std::string filename, T xlength, T ylength, int imax, int jmax, std::unique_ptr<std::unique_ptr<T[]>[]>& U, std::unique_ptr<std::unique_ptr<T[]>[]>& V, std::unique_ptr<std::unique_ptr<T[]>[]>& P, std::unique_ptr<std::unique_ptr<T[]>[]>& U2, std::unique_ptr<std::unique_ptr<T[]>[]>& V2) 
{
	//colculate U2 and V2 as the values of U and V in the middle of cells
	for (int i = 1; i < imax + 1; i++) 
	{
		for (int j = 1; j < jmax + 1; j++) 
		{
			U2[i][j] = (U[i][j - 1] + U[i][j]) / 2;
			V2[i][j] = (V[i][j] + V[i + 1][j]) / 2;
		}
	}

	std::ofstream myfile;
	myfile.open(filename);
	if (!myfile.is_open()) 
	{
		std::cout << "ERROR: write_output_into_file failed to open file";
		return;
	}
	myfile << xlength << "\n" << ylength << "\n" << imax << "\n" << jmax << "\n";
	for (int i = 1; i < imax + 1; i++) 
	{
		for (int j = 1; j < jmax + 1; j++) 
		{
			myfile << U2[i][j] << "\n";
		}
		
	}
	for (int i = 1; i < imax + 1; i++) 
	{
		for (int j = 1; j < jmax + 1; j++) 
		{
			myfile << V2[i][j] << "\n";
		}
	}
	for (int i = 1; i < imax + 1; i++) 
	{
		for (int j = 1; j < jmax + 1; j++) 
		{
			myfile << P[i][j] << "\n";
		}
	}
}

//max value in matrix
template <typename T>
T Mmax(std::unique_ptr<std::unique_ptr<T[]>[]> &Mat, int imax, int jmax)
{
	T Mmax = 0;
		for (int i = 0; i < imax + 2; i++)
		{
			for (int j = 0; j < imax + 2; j++)
			{
				if (abs(Mat[i][j]) > Mmax)
				{
					Mmax = abs(Mat[i][j]);
				}
			}
		}
	return Mmax;
}


//derivatives
template <typename T>
T _u2_x(std::unique_ptr<std::unique_ptr<T[]>[]> &U, int i, int j, T delx, T alpha)
{
	return 1 / delx * (((U[i][j] + U[i + 1][j]) / 2)*((U[i][j] + U[i + 1][j]) / 2) - ((U[i - 1][j] + U[i][j]) / 2)*((U[i - 1][j] + U[i][j]) / 2)) + alpha * 1 / delx * (((abs(U[i][j] + U[i + 1][j])) / 2)*((U[i][j] - U[i + 1][j]) / 2) - ((abs(U[i - 1][j] + U[i][j])) / 2)*((U[i - 1][j] - U[i][j]) / 2));
}

template <typename T>
T _uv_y(std::unique_ptr<std::unique_ptr<T[]>[]> &U, std::unique_ptr<std::unique_ptr<T[]>[]> &V, int i, int j, T dely, T alpha)
{
	return 1 / dely * ((V[i][j] + V[i + 1][j]) / 2 * (U[i][j] + U[i][j + 1]) / 2 - (V[i][j - 1] + V[i + 1][j - 1]) / 2 * (U[i][j - 1] + U[i][j]) / 2) + alpha * 1 / dely * ((abs(V[i][j] + V[i + 1][j])) / 2 * (U[i][j] - U[i][j + 1]) / 2 - (abs(V[i][j - 1] + V[i + 1][j - 1])) / 2 * (U[i][j - 1] - U[i][j]) / 2);
}

template <typename T>
T _2u_x2(std::unique_ptr<std::unique_ptr<T[]>[]> &U, int i, int j, T delx)
{
	return (U[i + 1][j] - 2 * U[i][j] + U[i - 1][j]) / (delx*delx);
}

template <typename T>
T _2u_y2(std::unique_ptr<std::unique_ptr<T[]>[]> &U, int i, int j, T dely)
{
	return (U[i][j + 1] - 2 * U[i][j] + U[i][j - 1]) / (dely*dely);
}

template <typename T>
T _p_x(std::unique_ptr<std::unique_ptr<T[]>[]> &P, int i, int j, T delx)
{
	return (P[i + 1][j] - P[i][j]) / delx;
}

template <typename T>
T _v2_y(std::unique_ptr<std::unique_ptr<T[]>[]> &V, int i, int j, T dely, T alpha)
{
	return 1 / dely * (((V[i][j] + V[i][j + 1]) / 2)*((V[i][j] + V[i][j + 1]) / 2) - ((V[i][j - 1] + V[i][j]) / 2)*((V[i][j - 1] + V[i][j]) / 2)) + alpha * 1 / dely * (((abs(V[i][j] + V[i][j + 1])) / 2)*((V[i][j] - V[i][j + 1]) / 2) - ((abs(V[i][j - 1] + V[i][j])) / 2)*((V[i][j - 1] - V[i][j]) / 2));
}

template <typename T>
T _uv_x(std::unique_ptr<std::unique_ptr<T[]>[]> &U, std::unique_ptr<std::unique_ptr<T[]>[]> &V, int i, int j, T delx, T alpha)
{
	return 1 / delx * ((U[i][j] + U[i][j + 1]) / 2 * (V[i][j] + V[i + 1][j]) / 2 - (U[i - 1][j] + U[i - 1][j + 1]) / 2 * (V[i - 1][j] + V[i][j]) / 2) + alpha * 1 / delx * ((abs(U[i][j] + U[i][j + 1])) / 2 * (V[i][j] - V[i + 1][j]) / 2 - (abs(U[i - 1][j] + U[i - 1][j + 1])) / 2 * (V[i - 1][j] - V[i][j]) / 2);
}

template <typename T>
T _2v_x2(std::unique_ptr<std::unique_ptr<T[]>[]> &V, int i, int j, T delx)
{
	return (V[i + 1][j] - 2 * V[i][j] + V[i - 1][j]) / (delx*delx);
}

template <typename T>
T _2v_y2(std::unique_ptr<std::unique_ptr<T[]>[]> &V, int i, int j, T dely)
{
	return (V[i][j + 1] - 2 * V[i][j] + V[i][j - 1]) / (dely*dely);
}

template <typename T>
T _p_y(std::unique_ptr<std::unique_ptr<T[]>[]> &P, int i, int j, T dely)
{
	return (P[i][j + 1] - P[i][j]) / dely;
}


//Init Fields
template <typename T>
void init_fields(std::unique_ptr<std::unique_ptr<T[]>[]> &Mat, int imax, int jmax, T MI) 
{
	for (int i = 0; i < imax + 2; i++) 
	{
		for (int j = 0; j < jmax + 2; j++) 
		{
			Mat[i][j] = MI;
		}
	}
}

//compute delt
template <typename T>
T _delt(std::unique_ptr<std::unique_ptr<T[]>[]> &U, std::unique_ptr<std::unique_ptr<T[]>[]> &V, int imax, int jmax, T tau, T Re, T delt, T delx, T dely)
{
	if (tau < 0)
	{
		return delt;
	}
	T u_max = Mmax(U, imax, jmax);
	if (u_max == 0)
	{
		u_max = delx / (Re / (2 * (1 / (delx*delx) + 1 / (dely*dely))));
	}
	T v_max = Mmax(V, imax, jmax);
	if (v_max == 0)
	{
		v_max = delx / (Re / (2 * (1 / (delx*delx) + 1 / (dely*dely))));
	}
	return std::min(Re / (2 * (1 / (delx*delx) + 1 / (dely*dely))), std::min(delx / u_max, dely / v_max));
}


template <typename T>
void boundary_conditions_UV(std::unique_ptr<std::unique_ptr<T[]>[]> &U, std::unique_ptr<std::unique_ptr<T[]>[]> &V, int imax, int jmax)
{
	for (int j = 1; j < jmax + 1; j++) 
	{
		U[0][j] = 0;
		U[imax][j] = 0;
		V[0][j] = -V[1][j];
		V[imax + 1][j] = -V[imax][j];
	}
	for (int i = 1; i < imax + 1; i++) 
	{
		U[i][0] = -U[i][1];
		U[i][jmax + 1] = -U[i][jmax];
		V[i][0] = 0.0;
		V[i][jmax] = 0.0;
	}
}

//copmute F, G
template <typename T>
void compute_FG(std::unique_ptr<std::unique_ptr<T[]>[]> &F, std::unique_ptr<std::unique_ptr<T[]>[]> &G, std::unique_ptr<std::unique_ptr<T[]>[]> &U, std::unique_ptr<std::unique_ptr<T[]>[]> &V, int imax, int jmax, T delt, T delx, T dely, T Re, T GX, T GY, T alpha)
{
	//F


	for (int i = 1; i < imax; i++)
	{
		for (int j = 1; j < jmax + 1; j++)
		{
			F[i][j] = U[i][j] + delt * (1 / Re * (_2u_x2(U, i, j, delx) + _2u_y2(U, i, j, dely)) - _u2_x(U, i, j, delx, alpha) - _uv_y(U, V, i, j, dely, alpha) + GX);
		}
	}

	//boundary condition
	for (int j = 1; j < jmax + 1; j++)
	{
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}

	//G


	for (int i = 1; i < imax + 1; i++)
	{
		for (int j = 1; j < jmax; j++)
		{
			G[i][j] = V[i][j] + delt * (1 / Re * (_2v_x2(V, i, j, delx) + _2v_y2(V, i, j, dely)) - _uv_x(U, V, i, j, delx, alpha) - _v2_y(V, i, j, dely, alpha) + GY);
		}
	}

	//boundary condition
	for (int i = 1; i < imax + 1; i++)
	{
		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
	}
}

//compute RHS
template <typename T>
void compute_RHS(std::unique_ptr<std::unique_ptr<T[]>[]> &F, std::unique_ptr<std::unique_ptr<T[]>[]> &G, std::unique_ptr<std::unique_ptr<T[]>[]> &RHS, int imax, int jmax, T delt, T delx, T dely)
{  
	for (int i = 1; i < imax + 1; i++) 
	{
		for (int j = 1; j < jmax + 1; j++) 
		{
			RHS[i][j] = (1 / delt) * ((F[i][j] - F[i - 1][j]) / delx + (G[i][j] - G[i][j - 1]) / dely);
		}
	}
}

//SOR Cycle
template <typename T>
void SOR(std::unique_ptr<std::unique_ptr<T[]>[]> &RHS, std::unique_ptr<std::unique_ptr<T[]>[]> &P, int itermax, int imax, int jmax, T omg, T delx, T dely, T eps)
{       
	int it = 0;
	T res =eps+1;
	T r;
	std::unique_ptr<std::unique_ptr<double[]>[]> P_old = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++)
	{
		P_old[i] = std::make_unique<double[]>(jmax + 2);
	}
	while (it < itermax && res > eps) 
	{
		
		T res = 0;
		
		for (int i = 0; i < imax + 2; i++)
		{
			for (int j = 0; j < imax + 2; j++)
			{
				P_old[i][j] = P[i][j];
			}
		}
		
		//boundary condition P
		for (int j = 1; j < jmax + 1; j++) 
		{            
			P[0][j] = P[1][j];
			P[imax + 1][j] = P[imax][j];
		}
		for (int i = 1; i < imax + 1; i++) 
		{
			P[i][0] = P[i][1];
			P[i][jmax + 1] = P[i][jmax];
		}
		
		
		
		for (int i = 1; i < imax + 1; i++) 
		{
			for (int j = 1; j < jmax + 1; j++) 
			{
				P[i][j] = (1 - omg)*P_old[i][j] + (omg / (2 * (1 / (delx*delx) + 1 / (dely*dely))))*((P_old[i + 1][j] + P[i - 1][j]) / (delx*delx) + (P_old[i][j + 1] + P[i][j - 1]) / (dely*dely) - RHS[i][j]);
			}
		}
		for (int i = 1; i < imax + 1; i++) 
		{
			for (int j = 1; j < jmax + 1; j++) 
			{
				
				r = ((P[i + 1][j] - 2 * P[i][j] + P[i][j - 1]) / (delx*delx) + (P[i][j + 1] - 2 * P[i][j] + P[i][j - 1]) / (dely*dely) - RHS[i][j]);
				res += (r*r) / (imax*jmax);
			}
		}
		res = sqrt(res);     
		it += 1;
	}
}

template <typename T>
void compute_UV(std::unique_ptr<std::unique_ptr<T[]>[]> &U, std::unique_ptr<std::unique_ptr<T[]>[]> &V, std::unique_ptr<std::unique_ptr<T[]>[]> &F, std::unique_ptr<std::unique_ptr<T[]>[]> &G, std::unique_ptr<std::unique_ptr<T[]>[]> &P, int imax, int jmax, T delt, T delx, T dely) 
{
	for (int i = 1; i < imax; i++) 
	{
		for (int j = 1; j < jmax + 1; j++) 
		{
			U[i][j] = F[i][j] - (delt / delx)*(P[i + 1][j] - P[i][j]);
		}
	}
	for (int i = 1; i < imax + 1; i++) 
	{
		for (int j = 1; j < jmax; j++) 
		{
			V[i][j] = G[i][j] - (delt / delx)*(P[i][j + 1] - P[i][j]);
		}
	}
}


						
																												

																											
