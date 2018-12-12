#pragma once
#include <iostream>
#include <memory>
#include <string>
#include <fstream>
#include <cmath>


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

	void f(int n);

	//only with knots in [0,1]
	int find_intervall(T x);
	void de_boor(T x, std::unique_ptr<T[]>& result); //warum muss result als referenz übergeben werden? ist doch ein pointer
	void write_to_file(std::string filename, int n);

	void create_interpol_parameters(std::unique_ptr<std::unique_ptr<T[]>[]>& inter_pts, std::unique_ptr<T[]>& res, int num_pts);
	void create_knots(std::unique_ptr<T[]>& param, int num_inter_pts, int degree);
	void print_knots();
	T eval_bspline(T x, int i, int _degree);

private:
	//
	int degree;
	int num_knots;
	int num_ctrpts;
	
	std::unique_ptr<T[]> knots; //has length num_knots+2*degree
	std::unique_ptr<std::unique_ptr<T[]>[]> ctrpts;

};

template <typename T, int Dim>
void Spline<T, Dim>::f(int n) {
	knots = std::make_unique<T[]>(n);
}


template <typename T, int Dim>
Spline<T, Dim>::Spline() :
	degree(0), num_knots(0), num_ctrpts(0), knots(nullptr), ctrpts(nullptr)
{}

template <typename T, int Dim>
Spline<T, Dim>::Spline(const Spline& other) :
	degree(other.degree), num_knots(other.num_knots), num_ctrpts(other.num_ctrpts)
{
	knots = std::make_unique<T[]>(num_knots+2*degree);
	ctrpts = std::make_unique<std::unique_ptr<T[]>[]>(num_ctrpts);
	for (int i = 0; i < num_ctrpts; i++) ctrpts[i] = std::make_unique<T[]>(Dim);

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
	knots = std::make_unique<T[]>(num_knots+2*degree);
	ctrpts = std::make_unique<std::unique_ptr<T[]>[]>(num_ctrpts);
	for (int i = 0; i < num_ctrpts; i++) ctrpts[i] = std::make_unique<T[]>(Dim);

	for (int i = 0; i < num_knots+2*degree; i++) knots[i] = 0;
	for (int i = 0; i < num_ctrpts; i++) {
		for (int j = 0; j < Dim; j++) ctrpts[i][j] = 0;
	}
}

template <typename T, int Dim>
void Spline<T, Dim>::assemble(int _num_knots, int _num_ctrpts, int _Dim, T* _knots, T** _ctrpts) {
	if (_num_knots != _num_ctrpts + 1 - degree) {
		std::cout << "number of ctr points, knots and the degree didnt correspond\n";
		exit(0);
	}
	if (_num_knots != num_knots || _num_ctrpts != num_ctrpts || Dim!=_Dim) {
		std::cout << "input sized didnt correspond with spline sizes\n";
		exit(0);
	}
	for (int i = 0; i < degree; i++) knots[i] = 0;
	for (int i = degree; i < degree + num_knots; i++) {
		if (_knots[i - degree] < 0 || _knots[i - degree] >1) {
			std::cout << "_knots entry is out of bounds\n";
			exit(0);
		}
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

	knots = std::make_unique<T[]>(num_knots+2*degree);
	ctrpts = std::make_unique<std::unique_ptr<T[]>[]>(num_ctrpts);
	for (int i = 0; i < num_ctrpts; i++) ctrpts[i] = std::make_unique<T[]>(Dim);

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
void Spline<T,Dim>::de_boor(T x, std::unique_ptr<T[]>& result) {
	int k = find_intervall(x);
	if (k == -1) return;
	T alpha;

	std::unique_ptr<std::unique_ptr<T[]>[]> d = std::make_unique<std::unique_ptr<T[]>[]>(degree + 1);
	for (int i = 0; i < degree + 1; i++) d[i] = std::make_unique<T[]>(Dim);

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
	std::unique_ptr<T[]> res = std::make_unique<T[]>(Dim);
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
void Spline<T,Dim>::create_interpol_parameters(std::unique_ptr<std::unique_ptr<T[]>[]>& inter_pts, std::unique_ptr<T[]>& res, int num_pts) {
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

//sets num_knots, degree and creates knots
template <typename T, int Dim>
void Spline<T, Dim>::create_knots(std::unique_ptr<T[]>& param, int num_inter_pts, int _degree) {
	std::cout << "\nbeginning of CREATE KNOTS\n";
	
	T sum = 0;
	degree = _degree;
	num_knots = num_inter_pts + 1 - _degree;
	knots = std::make_unique<T[]>(num_knots+2*degree);
	for (int i = 0; i <= degree; i++) knots[i] = 0;
	for (int i = num_inter_pts; i <= num_inter_pts+degree+1; i++) knots[i] = 1;
	for (int j = 1; j <= num_inter_pts - degree - 1; j++) {
		sum = 0;
		for (int i = j; i <= j + degree - 1; i++) sum += param[i];
		knots[j + degree] = sum / degree;
	}

	std::cout << "\nend of CREATE KNOTS\n";
	return;
}

template <typename T>
int find_intervall2(T x, std::unique_ptr<T[]>& knots, int num_knots, int degree) {
	if (num_knots <= 1) { std::cout << "less than two knots: cant find nonexistent intervall\n"; return -1; }

	//we look for k such that 
	//x in [t_k, t_k+1) (see De Boors algorithm)
	if (x < 0 || x > 1) { std::cout << "x has to be in [0,1]\n"; return -1; }

	//we let [t_(degree+num_nodes-2), 1] be a closed intervall, so we can evaluate the curve at 1
	if (x == 1) return degree + num_knots - 2;

	for (int k = degree; k < degree + num_knots; k++) {
		if (knots[k] > x) return (k - 1);
	}
	return -1;
}

template <typename T, int Dim>
void Spline<T, Dim>::print_knots() {
	std::cout << "\nknots of spline are:\n";
	for (int i = 0; i < num_knots + 2 * degree; i++) std::cout << knots[i] << ", ";
	std::cout << "\n";
}


//cox, de boor formula
//evaluates B spline i of degree k at position x
template <typename T, int Dim>
T Spline<T, Dim>::eval_bspline(T x, int i, int k) { //Was ist wenn i quasi nicht mehr relevant ist?
	std::cout << "\nbeginning EVAL BSPLINE\n";
	if (k == 0) {
		if (i == find_intervall(x)) return 1;
		else return 0;
	}
	T omega_ik;
	T omega_iplus1k;
	if (knots[i] != knots[i + k - 1]) omega_ik = (x - knots[i]) / (knots[i + k - 1] - knots[i]);
	else omega_ik = 0;
	if (knots[i + 1] != knots[i + k]) omega_iplus1k = (x - knots[i + 1]) / (knots[i + k] - knots[i + 1]);
	else omega_iplus1k = 0;
	
	else return omega_ik*eval_bspline(x, i, _degree-1)+omega_iplus1k*eval_bspline(x, i+1, _degree-1);
	
}



