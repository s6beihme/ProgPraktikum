#include "driven_cavity_parallel.h"
#include "driven_cavity.h"


int main(int argc, char **argv) {
	int  myrank, numberoftasks;
	MPI_Comm cartesian_handle;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numberoftasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	//check, if inputfile and name of outputfiles were given
	if (argc < 3) {
		std::cout << "this program requires two arguments to be passed:" <<
			"the name of the input file and the base name of the output files";
		MPI_Finalize();
		return 0;
	}

	std::string inputname = argv[1];
	std::string outputname = argv[2];
	int imax_global; int jmax_global; double xlength; double ylength; double delt; double t_end; double tau; double del_vec; double eps; double omg; double alpha; int itermax = 0; double GX; double GY; int Re; double UI; double VI; double PI;
	read_parameters_from_file(inputname, imax_global, jmax_global, xlength, ylength, delt, t_end, tau, del_vec, eps, omg, alpha, itermax, GX, GY, Re, UI, VI, PI);
	
	//set up cartesian handle of processes
	int dims[2] = { 0,0 };
	MPI_Dims_create(numberoftasks, 2, dims);
	
	if (imax_global < jmax_global) {
		int temp = dims[0];
		dims[0] = dims[1];
		dims[1] = temp;
	}
	if (imax_global%dims[0] != 0 || jmax_global % dims[1] != 0) {
		if (myrank == 0) {
			std::cout << "for this program to run, it is required, that the data grid can be evenly devided " <<
				"into sub_grids of equal size, such that each process handles one of the sub_grids. " <<
				"That isnt the case here. Please use a more suitable number of processes or another sized grid";
		}
		MPI_Finalize();
		return 0;
	}

	int periods[2] = { 0,0 };
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cartesian_handle);
	int coords[2] = { 0,0 };
	MPI_Cart_coords(cartesian_handle, myrank, 2, coords);

	//set up necessary variables
	int imax, jmax;

	imax = imax_global / dims[0];
	jmax = jmax_global / dims[1];

	std::string appendix;
	int number_of_files = 0;
	double del_vec2 = del_vec;
	double res = 0;
	double delx = xlength / imax_global;
	double dely = ylength / jmax_global;

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
	double* row = new double[jmax];
	double* column = new double[imax];

	initialize_fields_UVP<double>(imax, jmax, UI, VI, PI, U, V, P);
	exchange_boundary_values_UV_or_FG_between_processes(cartesian_handle, U, V, column, row, imax, jmax, coords, dims);
	double t = 0;

	while (t < t_end) {
		
		calculate_delt_parallel<double>(cartesian_handle, U, V, imax, jmax, tau, Re, delx, dely, myrank, delt);

		apply_boundary_conditions_UV_on_outer_processes(U, V, imax, jmax, coords, dims);

		calculate_F_G_parallel(coords, dims, imax, jmax, delt, delx, dely, Re, alpha, GX, GY, U, V, F, G);
		//exchange boundary values of F and G so the RHS will be computed correctly
		exchange_boundary_values_UV_or_FG_between_processes(cartesian_handle, F, G, column, row, imax, jmax, coords, dims);
		calculate_RHS(imax, jmax, delt, delx, dely, F, G, RHS);

		calculate_Pressure_with_SOR_parallel(cartesian_handle, P, P2, RHS, imax, jmax, itermax, delx, dely, eps, omg, column, row, coords, dims, myrank, res);

		calculate_U_and_V_parallel(coords, dims, imax, jmax, delt, delx, dely, F, G, P, U, V);

		exchange_boundary_values_UV_or_FG_between_processes(cartesian_handle, U, V, column, row, imax, jmax, coords, dims);

		if (t > del_vec) {
			number_of_files += 1;
			//now find the right string to append to outputname
			if (number_of_files < 10) appendix = "_00" + std::to_string(number_of_files);
			else {
				if (number_of_files < 100) appendix = "_0" + std::to_string(number_of_files);
				else appendix = "_" + std::to_string(number_of_files);
			}
			write_output_into_file_parallel(outputname + appendix, cartesian_handle, myrank, coords, dims, row, U, V, P, U2, V2, imax, jmax, imax_global, jmax_global, xlength, ylength);
			del_vec = t + del_vec2;
		}

		t += delt;
	}
	number_of_files += 1;
	//now find the right string to append to outputname
	if (number_of_files < 10) appendix = "_00" + std::to_string(number_of_files);
	else {
		if (number_of_files < 100) appendix = "_0" + std::to_string(number_of_files);
		else appendix = "_" + std::to_string(number_of_files);
	}
	write_output_into_file_parallel(outputname + appendix, cartesian_handle, myrank, coords, dims, row, U, V, P, U2, V2, imax, jmax, imax_global, jmax_global, xlength, ylength);
	MPI_Finalize();
}