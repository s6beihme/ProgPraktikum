#pragma once

class SplineCurves
{
public:
	SplineCurves(int _degree, int _knots, int _ctrlp);
	~SplineCurves();
	void sc_assemble(double* knots, double** conp, int _degree, int _knots, int _ctrlp);
	int knotintervall(double x);
	double * de_boor(double x);
private:
	int degree;
	int knots;
	int ctrlp;
	double* knotvec;
	double** ctrlp_vec;
};