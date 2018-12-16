#pragma once
#include <iostream>
#include <memory>
#include <string>
#include <fstream>
#include <cmath>
#include <cassert>
#include <vector>
#include "VectorClass.h"
#include "CsrMatrixClass.h"

template <typename T, int Dim>
class Spline {
public:
	Spline();
	Spline(const Spline& other);
	Spline(int deg, int _num_ctrpts);
	void assemble(int _num_knots, int _num_ctrpts, int _Dim, T* _knots, T** _ctrpts);

	Spline& operator = (const Spline& other);
	
	bool operator ==(const Spline& other);
	bool operator !=(const Spline& other);


	//only with knots in [0,1]
	int find_intervall(T x);
	void de_boor(T x, std::vector<T>& result); 
	void write_to_file(std::string filename, int n);

	void create_interpol_parameters(std::vector<std::vector<T>>& inter_pts, std::vector<T>& res, int num_pts);
	void create_knots(std::vector<T>& param, int num_inter_pts, int degree);
	void print_knots();
	T eval_bspline(T x, int i);
	void create_ctrpts_from_inter_pts(std::vector<std::vector<T>>& inter_pts, int num_pts, int _degree);

private:
	//
	int degree;
	int num_knots;
	int num_ctrpts;
	
	std::vector<T> knots; //has length num_knots+2*degree
	std::vector<std::vector<T>> ctrpts;

};


template <typename T, int Dim>
Spline<T, Dim>::Spline() :
	degree(0), num_knots(0), num_ctrpts(0)
{}

template <typename T, int Dim>
Spline<T, Dim>::Spline(const Spline<T, Dim>& other) :
	degree(other.degree), num_knots(other.num_knots), num_ctrpts(other.num_ctrpts)
{
	knots = std::vector<T>(num_knots+2*degree);
	ctrpts = std::vector<std::vector<T>>(num_ctrpts);
	for (int i = 0; i < num_ctrpts; i++) ctrpts[i] = std::vector<T>(Dim);

	for (int i = 0; i < num_knots+2*degree; i++) knots[i] = other.knots[i];
	for (int i = 0; i < num_ctrpts; i++) {
		for (int j = 0; j < Dim; j++) ctrpts[i][j] = other.ctrpts[i][j];
	}
}

template <typename T, int Dim>
Spline<T, Dim>::Spline(int deg, int _num_ctrpts) :
	degree(deg)
{
	num_knots = _num_ctrpts + 1 - degree;
	num_ctrpts = _num_ctrpts;
	knots = std::vector<T>(num_knots+2*degree, 0);
	ctrpts = std::vector<std::vector<T>>(num_ctrpts);
	for (int i = 0; i < num_ctrpts; i++) ctrpts[i] = std::vector<T>(Dim, 0);

}

template <typename T, int Dim>
void Spline<T, Dim>::assemble(int _num_knots, int _num_ctrpts, int _Dim, T* _knots, T** _ctrpts) {
	assert(_num_knots == _num_ctrpts + 1 - degree);
	assert(_num_knots == num_knots && _num_ctrpts != num_ctrpts && Dim != _Dim);

	for (int i = 0; i < degree; i++) knots[i] = 0;
	for (int i = degree; i < degree + num_knots; i++) {
		assert(_knots[i - degree] >= 0 && _knots[i - degree] <=1);
		
		knots[i] = _knots[i - degree];
	}
	for (int i = degree + num_knots; i < num_knots + 2 * degree; i++) knots[i] = 1;
	for (int i = 0; i < num_ctrpts; i++) {
		for (int j = 0; j < Dim; j++) ctrpts[i][j] = _ctrpts[i][j];
	}
}

template <typename T, int Dim>
Spline<T, Dim>& Spline<T, Dim>::operator=(const Spline<T, Dim>& other) {
	degree = other.degree;
	num_ctrpts = other.num_ctrpts;
	num_knots = other.num_knots;

	knots = std::vector<T>(num_knots+2*degree);
	ctrpts = std::vector<std::vector<T>>(num_ctrpts);
	for (int i = 0; i < num_ctrpts; i++) ctrpts[i] = std::vector<T>(Dim);

	for (int i = 0; i < num_knots+2*degree; i++) knots[i] = other.knots[i];
	for (int i = 0; i < num_ctrpts; i++) {
		for (int j = 0; j < Dim; j++) ctrpts[i][j] = other.ctrpts[i][j];
	}

	return *this;
}

template <typename T, int Dim>
bool Spline<T, Dim>::operator==(const Spline<T, Dim>& other) {
	if (degree != other.degree || num_ctrpts != other.num_ctrpts || num_knots != other.num_knots) return false;
	for (int i = 0; i < num_knots+2*degree; i++) if (knots[i] != other.knots[i]) return false;
	for (int i = 0; i < num_ctrpts; i++) {
		for (int j = 0; j < Dim; j++) if (ctrpts[i][j] != other.ctrpts[i][j]) return false;
	}
	return true;
}

template <typename T, int Dim>
bool Spline<T, Dim>::operator!=(const Spline<T, Dim>& other) {
	return !((*this) == other);
}

template <typename T, int Dim>
int Spline<T, Dim>::find_intervall(T x) {
	if (num_knots <= 1) { std::cout << "less than two knots: cant find nonexistent intervall\n"; return -1; }
	
	//we look for k such that x in [t_k, t_k+1) (see De Boors algorithm)
	if (x < 0 || x > 1) { std::cout << "x has to be in [0,1]\n"; return -1; }

	//we let [t_(degree+num_nodes-2), 1] be a closed intervall, so we can evaluate the curve at 1
	if (x == 1) return degree + num_knots - 2;
	for (int k = degree; k < degree + num_knots; k++) {
		if (knots[k] > x) return (k - 1);
	}
	return -1;
}

template <typename T, int Dim>
void Spline<T,Dim>::de_boor(T x, std::vector<T>& result) {
	int k = find_intervall(x);
	if (k == -1) return;
	T alpha;

	std::vector<std::vector<T>> d(degree + 1);
	for (int i = 0; i < degree + 1; i++) d[i] = std::vector<T>(Dim);

	for (int j = 0; j < degree + 1; j++) {
		for (int i = 0; i < Dim; i++) d[j][i] = ctrpts[j + k - degree][i];
	}

	for (int r = 1; r <= degree; r++) {
		for (int j = degree; j >= r; j--) {
			alpha = (x - knots[j + k - degree]) / (knots[j + 1 + k - r] - knots[j + k - degree]);
			for (int i = 0; i < Dim; i++) d[j][i] = (1 - alpha)*d[j - 1][i] + alpha * d[j][i];
		}
	}

	for (int i = 0; i < Dim; i++) result[i] = d[degree][i];
}

template <typename T, int Dim>
void Spline<T, Dim>::write_to_file(std::string filename, int n) {
	std::ofstream myfile;
	myfile.open(filename);
	if (myfile.is_open() == false) {
		std::cout << "\nwrite to file failed: couldnt open file\n";
		return;
	}
	std::vector<T> res(Dim);
	myfile << Dim << "\n";
	for (int i = 0; i < n-1; i++) {
		de_boor(0 + (T)i / (T)n, res);
		for (int j = 0; j < Dim - 1; j++) myfile << res[j] << " ";
		myfile << res[Dim - 1] << "\n";
	}
	de_boor(1, res);
	for (int j = 0; j < Dim - 1; j++) myfile << res[j] << " ";
	myfile << res[Dim - 1];
	myfile.close();
}

//to be used in interpolation
template <typename T, int Dim>
void Spline<T,Dim>::create_interpol_parameters(std::vector<std::vector<T>>& inter_pts, std::vector<T>& res, int num_pts) {
	T d = 0;
	T r = 0;
	res[0] = 0;
	res[num_pts - 1] = 1;
	for (int k = 1; k < num_pts; k++) {
		r = 0;
		for (int i = 0; i < Dim; i++) r += (inter_pts[k][i] - inter_pts[k - 1][i])*(inter_pts[k][i] - inter_pts[k - 1][i]);
		d += sqrt(r);
	}

	for (int k = 1; k < num_pts - 1; k++) {
		r = 0;
		for (int i = 0; i < Dim; i++) r += (inter_pts[k][i] - inter_pts[k - 1][i])*(inter_pts[k][i] - inter_pts[k - 1][i]);

		res[k] = res[k - 1] + (sqrt(r)) / d;
	}
}

//sets num_knots, degree, num_ctrpts and creates knots
template <typename T, int Dim>
void Spline<T, Dim>::create_knots(std::vector<T>& param, int num_inter_pts, int _degree) {
	
	T sum = 0;
	num_ctrpts = num_inter_pts;
	degree = _degree;
	num_knots = num_inter_pts + 1 - _degree;
	knots = std::vector<T>(num_knots+2*degree);
	for (int i = 0; i <= degree; i++) knots[i] = 0;
	for (int i = num_inter_pts; i <= num_inter_pts+degree+1; i++) knots[i] = 1;
	for (int j = 1; j <= num_inter_pts - degree - 1; j++) {
		sum = 0;
		for (int i = j; i <= j + degree - 1; i++) sum += param[i];
		knots[j + degree] = sum / degree;
	}

	return;
}


template <typename T, int Dim>
void Spline<T, Dim>::print_knots() {
	std::cout << "\nknots of spline are:\n";
	for (unsigned int i = 0; i < knots.size(); i++) std::cout << knots[i] << ", ";
	std::cout << "\n";
}


//cox, de boor formula
//evaluates B spline i of degree k at position x
//needs knots and num_knots to be known
//NOTE: i refers to knots[i]. so to get the first "real" knot, call with i=degree
template <typename T, int Dim>
T Spline<T, Dim>::eval_bspline(T x, int i) {
	if ((i == 0 && x == knots[0]) || (i == num_knots + degree - 2 && x == knots[num_knots + 2 * degree - 1])) return 1;
	if (x < knots[i] || x >= knots[i + degree + 1]) return 0;
	if (i > num_knots + degree - 3) return 0;
	
	T saved = 0;
	T temp = 0;
	T Uleft = 0;
	T Uright = 0;

	std::vector<T> N(degree + 1);
	for (int j = 0; j <= degree; j++) {
		if (x >= knots[i + j] && x < knots[i + j + 1]) N[j] = 1;
		else N[j] = 0;
	}
	for (int k = 1; k <= degree; k++) {
		if (N[0] == 0) saved = 0;
		if (knots[i + k] == knots[i]) saved = 0;
		else saved = ((x - knots[i])*N[0]) / (knots[i + k] - knots[i]);
		for (int j = 0; j < degree - k + 1; j++) {
			Uleft = knots[i + j + 1];
			Uright = knots[i + j + k + 1];
			if (N[j + 1] == 0) {
				N[j] = saved;;
				saved = 0;
			}
			else {
				if (Uright != Uleft) temp = N[j + 1] / (Uright - Uleft);
				else temp = 0;
				N[j] = saved + (Uright - x)*temp;
				saved = (x - Uleft)*temp;
			}
			
		}
	}
	T result = N[0];
	return result;
}

template <typename T, int Dim>
void Spline<T, Dim>::create_ctrpts_from_inter_pts(std::vector<std::vector<T>>& inter_pts, int num_inter_pts, int _degree) {
	
	ctrpts = std::vector<std::vector<T>>(num_inter_pts);
	for (int i = 0; i < num_inter_pts; i++) ctrpts[i] = std::vector<T>(Dim);

	std::vector<T> parameters(num_inter_pts);
	create_interpol_parameters(inter_pts, parameters, num_inter_pts);
	std::cout << std::endl;
	std::cout << std::endl;
	create_knots(parameters, num_inter_pts, _degree);
	//now degree, num_knots, num_ctrpts and knots are set

	//check how many non zero entries matrix has
	int nnz=0;
	int num_entr_in_row = 0;
	for(int k=0; k<num_ctrpts; k++) {
		num_entr_in_row = 0;
		for (int i = 0; i < num_ctrpts; i++) {
			if (eval_bspline(parameters[k], i) != 0) {
				nnz++;
				num_entr_in_row++;
				if (num_entr_in_row == degree + 1) break;
			}
		}
	}

	T* vals = new T[nnz];
	int* row_inds = new int[nnz];
	int* col_inds = new int[nnz];

	//build matrix
	int pos = 0;
	for (int k = 0; k < num_ctrpts; k++) {
		for (int i = 0; i < num_ctrpts; i++) {
			if (eval_bspline(parameters[k], i) != 0) {
				vals[pos] = eval_bspline(parameters[k], i);
				col_inds[pos] = i;
				row_inds[pos] = k;
				pos++;
			}
		}
	}

	CsrMatrix<T> A(nnz, num_ctrpts, num_ctrpts);
	A.assemble(vals, row_inds, col_inds, nnz);

	Vector<T> temp_result(num_ctrpts), b(num_ctrpts);

	//solve for each dimension
	for (int i = 0; i < Dim; i++) {
		for (int j = 0; j < num_ctrpts; j++) {b[j] = inter_pts[j][i];  }
		A.gs_solve(b, temp_result);
		for (int j = 0; j < num_ctrpts; j++) ctrpts[j][i] = temp_result[j];
	}

}


