#include <iostream>
#include "Splines.h"


SplineCurves::SplineCurves(int _degree, int _knots, int _ctrlp)
{
	degree = _degree;
	knots = _knots;
	ctrlp = _ctrlp;
	knotvec = new double[2 * (degree - 1) + knots];
	if (knotvec == NULL)
	{
		std::cout << "failed to create SplineCurve object: failed to allocate memory for knotvec\n";
		return;
	}
	ctrlp_vec = new double*[ctrlp];
	if (ctrlp_vec == NULL)
	{
		std::cout << "failed to create SplineCurve object: failed to allocate memory for ctrlpvec\n";
		delete[] knotvec;
		return;
	}
	for (int i = 0; i < ctrlp; i++)
	{
		ctrlp_vec[i] = new double[2];
		if (ctrlp_vec[i] == NULL)
		{
			std::cout << "failed to create SplineCurve object: failed to allocate memory for ctrlpvec\n";
			delete[] knotvec;
			for (int j = 0; j < i; j++)
			{
				delete[] ctrlp_vec[j];
			}
			delete[] ctrlp_vec;
			return;
		}
	}
}

SplineCurves::~SplineCurves()
{
	delete[] knotvec;
	for (int i = 0; i < ctrlp; i++)
	{
		delete[] ctrlp_vec[i];
	}
	delete[] ctrlp_vec;
}


void SplineCurves::sc_assemble(double* knotsv, double** conp, int _degree, int _knots, int _ctrlp)
{
	if (_degree != degree || _knots != knots || _ctrlp != ctrlp)
	{
		std::cout << "check your input sizes\n";
		return;
	}

	for (int i = 0; i < 2 * (degree - 1) + knots; i++)
	{
		knotvec[i] = knotsv[i];
	}
	for (int i = 0; i < _ctrlp; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			ctrlp_vec[i][j] = conp[i][j];
		}
	}
}

int SplineCurves::knotintervall(double x)
{
	//i ONLY WENT UP TO knots BUT THERE ARE THE ADDITIONAL BOUNDARY POINTS !!! (I fixed it in a not so optimal way)
	for (int i = 1; i < 2 * (degree - 1) + knots; i++)
	{
		if (knotvec[i] > x)
		{
			if (i == 0)
			{
				return (i - 1 + degree);
			}
			return (i - 1);
		}
	}
	std::cout << "x is out of bounds\n";
	return -1;
}

double* SplineCurves::de_boor(double x)
{
	double alpha;
	double ** d;
	int k = knotintervall(x);

	//I ADDED THIS
	if (k == -1) return NULL;

	d = new double*[degree + 1];
	if (d == NULL)
	{
		std::cout << "sorry, we are out of memory\n";
		delete[] d;
		return NULL;
	}
	for (int i = 0; i <= degree; i++)
	{
		d[i] = new double[2];
		if (d[i] == NULL)
		{
			std::cout << "sorry, we are out of memory\n";
			for (int j = 0; j < i; j++)
			{
				delete[] d[j];
			}
			delete[] d;
			return NULL;
		}
	}
	for (int i = 0; i <= degree; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			d[i][j] = ctrlp_vec[i + k - degree][j];
		}
	}
	for (int i = 1; i <= degree; i++)
	{
		for (int j = degree; j > i - 1; j--)
		{
			alpha = (x - knotvec[j + k - degree]) / (knotvec[j + 1 + k - i] - knotvec[j + k - degree]);
			for (int l = 0; l < 2; l++)
			{
				d[j][l] = (1 - alpha)*d[j - 1][l] + alpha * d[j][l];
			}
		}
	}
	return d[degree];
}

int main()
{
	int p = 2;
	double knots[8] = { 0,0,0,0.3,0.5,1,1,1 };

	//DOUBLE CONTR[5][2] ISNT A DOUBLE**,	IT IS BASICALLY THE SAME AS DOUBLE CONTR[5*2]
	//double contr[5][2] = { {0,1}, {1,-0.2}, {3,2}, {6,1.3}, {7,1.5} };
	double** contr = new double*[5];
	for (int i = 0; i < 5; i++) contr[i] = new double[2];
	contr[0][0]=0;
	contr[0][1]=1;
	contr[1][0]=1;
	contr[1][1]=-0.2;
	contr[2][0]=3;
	contr[2][1]=2;
	contr[3][0]=6;
	contr[3][1]=1.3;
	contr[4][0]=7;
	contr[4][1]=1.5;
	SplineCurves B(2, 4, 5);
	B.sc_assemble(knots, contr, p, 4, 5);
	std::cout << "B.de_boor(0.4)\n" << B.de_boor(0.4)[0] << "\n" << B.de_boor(0.4)[1];



	return 0;
}