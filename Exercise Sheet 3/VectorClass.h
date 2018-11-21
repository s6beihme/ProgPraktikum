#pragma once
#include <iostream>

class Vector {
public:
	Vector(int _size);
	~Vector();
	void vec_assemble(double* vals);
	double operator * (Vector& v2);
	void print_vector();

	int get_size();
	double get_data(int i);
	void set_data(int i, double val);

private:
	int size;
	double* data;
};


