#include <iostream>
#include "driven_cavity.h"

int main() {
	/*int imax; int jmax; double xlength; double ylength; double delt; double t_end; double tau; double del_vec; double eps; double omg; double alpha; int itermax=0; double GX; double GY; double Re; double UI; double VI; double PI;
	read_parameters_from_file("param.txt", imax, jmax, xlength, ylength, delt, t_end, tau, del_vec, eps, omg, alpha, itermax, GX, GY, Re, UI, VI, PI);
	std::unique_ptr<std::unique_ptr<double[]>[]> U = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		U[i]= std::make_unique<double[]>(imax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> U2 = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		U2[i] = std::make_unique<double[]>(imax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> V = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		V[i] = std::make_unique<double[]>(imax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> V2 = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		V2[i] = std::make_unique<double[]>(imax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> P = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		P[i] = std::make_unique<double[]>(imax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> P2 = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		P2[i] = std::make_unique<double[]>(imax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> F = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		F[i] = std::make_unique<double[]>(imax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> G = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		G[i] = std::make_unique<double[]>(imax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> RHS = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		RHS[i] = std::make_unique<double[]>(imax + 2);
	}
	double res = 0;


	initialize_fields_UVP<double>(imax, jmax, UI, VI, PI, U, V, P);
	//print_matrix(imax + 2, jmax + 2, U);
	double delx = xlength / imax;
	double dely = ylength / jmax;
	std::cout << "tau, Re: " << tau << " " << Re << "\n";
	std::cout <<"delx, dely:\n"<< delx << " " << dely << std::endl;

	calculate_delt<double>(imax, jmax, tau, Re, delx, dely, U, V, delt);
	std::cout << "delt: "<<delt << std::endl;

	apply_boundary_conditions_UV(imax, jmax, U, V);
	std::cout << "\napplied boundary conditions\n";
	
	std::cout << "\nalpha: " << alpha << ", GX: " << GX << ", GY: " << GY << "\n";
	calculate_F_G(imax, jmax, delt, delx, dely,Re, alpha, GX, GY, U, V, F, G);
	std::cout << "\ncalculated F, G\n";
	//print_matrix(imax + 2, jmax + 2, F);
	calculate_RHS(imax, jmax, delt, delx, dely, F, G, RHS);
	std::cout << "\ncalculated RHS\n";
	//print_matrix(imax + 2, jmax + 2, RHS);
	//std::cout << "\neps = " << eps << "\n";
	calculate_Pressure_with_SOR(imax, jmax, itermax, delx, dely, eps, omg, RHS, P, P2, res);
	std::cout << "\ncalculated P with SOR\n";
	calculate_U_and_V(imax, jmax, delt, delx, dely, F, G, P, U, V);
	std::cout << "\ncalculated U and V\n";
	//print_matrix(imax + 2, jmax + 2, U);
	//std::cout << "\nit = " << it << ", res = " << res << std::endl;
	write_output_into_file("outputfile.txt", xlength, ylength, imax, jmax, U, V, P, U2, V2);
	std::cout << "\noutput was written into file\n";*/
	char a[] = "abs";
	std::string str = a;
	std::cout << str+"def" << std::endl;

	return 0;
}