#include "driven_cavity.h"

int main(int argc, char* argv[]) {
	
	//convert the input to std strings for easy handling
	if (argc < 3) {
		std::cout << "this program requires two arguments to be passed: the name of the input file and the base name of the output files";
		return 0;
	}

	std::string inputname = argv[1];
	std::string outputname = argv[2];
	
	//first set up all the neccessary variables
	int imax; int jmax; double xlength; double ylength; double delt; double t_end; double tau; double del_vec; double eps; double omg; double alpha; int itermax = 0; double GX; double GY; double Re; double UI; double VI; double PI;
	read_parameters_from_file(inputname, imax, jmax, xlength, ylength, delt, t_end, tau, del_vec, eps, omg, alpha, itermax, GX, GY, Re, UI, VI, PI);
	
	//now set up all the matrices
	std::unique_ptr<std::unique_ptr<double[]>[]> U = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		U[i] = std::make_unique<double[]>(jmax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> U2 = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		U2[i] = std::make_unique<double[]>(jmax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> V = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		V[i] = std::make_unique<double[]>(jmax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> V2 = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		V2[i] = std::make_unique<double[]>(jmax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> P = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		P[i] = std::make_unique<double[]>(jmax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> P2 = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		P2[i] = std::make_unique<double[]>(jmax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> F = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		F[i] = std::make_unique<double[]>(jmax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> G = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		G[i] = std::make_unique<double[]>(jmax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> RHS = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		RHS[i] = std::make_unique<double[]>(jmax + 2);
	}

	

	std::string appendix;
	int number_of_files = 0;
	double del_vec2 = del_vec;
	double res = 0;
	double delx = xlength / imax;
	double dely = ylength / jmax;
	
	initialize_fields_UVP<double>(imax, jmax, UI, VI, PI, U, V, P);
	std::cout << "\nHERE\n";
	double t = 0;
	
	while (t < t_end) {
		calculate_delt(imax, jmax, tau, Re, delx, dely, U, V, delt);

		apply_boundary_conditions_UV(imax, jmax, U, V);

		calculate_F_G(imax, jmax, delt, delx, dely, Re, alpha, GX, GY, U, V, F, G);

		calculate_RHS(imax, jmax, delt, delx, dely, F, G, RHS);

		calculate_Pressure_with_SOR(imax, jmax, itermax, delx, dely, eps, omg, RHS, P, P2, res);

		calculate_U_and_V(imax, jmax, delt, delx, dely, F, G, P, U, V);

		if (t > del_vec) {
			number_of_files += 1;
			//now find the right string to append to outputname
			if(number_of_files<10) appendix = "_00"+ std::to_string(number_of_files);
			else {
				if(number_of_files<100) appendix = "_0"+ std::to_string(number_of_files);
				else appendix = "_"+std::to_string(number_of_files);
			}
			write_output_into_file(outputname+appendix, xlength, ylength, imax, jmax, U, V, P, U2, V2);
			del_vec = t + del_vec2;
		}

		t += delt;
	}
	number_of_files += 1;
	//now find the right string to append to outputname
	if (number_of_files < 10) appendix = "_00" + std::to_string(number_of_files);
	else {
		if (number_of_files < 100) appendix = "_0" + std::to_string(number_of_files);
		else appendix = "_"+std::to_string(number_of_files);
	}
	write_output_into_file(outputname+appendix, xlength, ylength, imax, jmax, U, V, P, U2, V2);
}