#include "driven_cavity.h"

int main(int argc, char* argv[]) {

	std::cout << "\nSTART\n";

	//convert input to strings
	if (argc < 3) {
		std::cout << "this program requires two arguments to be passed: the name of the input file and the base name of the output files";
		return 0;
	}

	std::string inputname = argv[1];
	std::string outputname = argv[2];

	//variables
	int imax; 
	int jmax; 
	double xlength; 
	double ylength; 
	double delt; 
	double t_end; 
	double tau; 
	double del_vec; 
	double eps; 
	double omg; 
	double alpha; 
	int itermax = 0; 
	double GX; 
	double GY; 
	double Re; 
	double UI; 
	double VI; 
	double PI;
	input(inputname, imax, jmax, xlength, ylength, delt, t_end, tau, del_vec, eps, omg, alpha, itermax, GX, GY, Re, UI, VI, PI);

	//matries
	std::unique_ptr<std::unique_ptr<double[]>[]> U = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		U[i] = std::make_unique<double[]>(jmax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> U_copy = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		U_copy[i] = std::make_unique<double[]>(jmax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> V = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		V[i] = std::make_unique<double[]>(jmax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> V_copy = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		V_copy[i] = std::make_unique<double[]>(jmax + 2);
	}
	std::unique_ptr<std::unique_ptr<double[]>[]> P = std::make_unique<std::unique_ptr<double[]>[]>(imax + 2);
	for (int i = 0; i < imax + 2; i++) {
		P[i] = std::make_unique<double[]>(jmax + 2);
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
	double t_vis = del_vec;
	double delx = xlength / imax;
	double dely = ylength / jmax;

	init_fields(U, imax, jmax, UI);
	init_fields(V, imax, jmax, VI);
	init_fields(P, imax, jmax, PI);
	
	double t = 0;

	while (t < t_end) {
		delt = _delt(U, V, imax, jmax, tau, Re, delt, delx, dely);
		
		boundary_conditions_UV(U, V, imax, jmax);

		//movement on top of surface
		for (int i = 1; i < imax + 1; i++)
		{
			U[i][jmax + 1] =2.0-U[i][jmax];
		}
		
		compute_FG(F, G, U, V, imax, jmax, delt, delx, dely, Re, GX, GY, alpha);
		
		compute_RHS(F, G, RHS, imax, jmax, delt, delx, dely);
		
		SOR(RHS, P, itermax, imax, jmax, omg, delx, dely, eps);
		
		compute_UV(U, V, F, G, P, imax, jmax, delt, delx, dely);
		
		if (t > t_vis) 
		{
			number_of_files += 1;
			if (number_of_files < 10) appendix = "_00" + std::to_string(number_of_files);
			else {
				if (number_of_files < 100) appendix = "_0" + std::to_string(number_of_files);
				else appendix = "_" + std::to_string(number_of_files);
			}
			output(outputname + appendix, xlength, ylength, imax, jmax, U, V, P, U_copy, V_copy);
			t_vis = t + del_vec;
		}

		t += delt;
	}
	number_of_files += 1;
	if (number_of_files < 10) appendix = "_00" + std::to_string(number_of_files);
	else {
		if (number_of_files < 100) appendix = "_0" + std::to_string(number_of_files);
		else appendix = "_" + std::to_string(number_of_files);
	}
	output(outputname + appendix, xlength, ylength, imax, jmax, U, V, P, U_copy, V_copy);

	std::cout << "\nEND\n";

}