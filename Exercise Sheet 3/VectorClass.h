#pragma once
#include <iostream>

class Vector {
public:
	Vector(int _size);
	~Vector();
	void vec_assemble(double* vals);
	double operator * (const Vector& v2) const;
	Vector operator - (const Vector&v2) const;
	void print_vector() const;

	int get_size() const;
	double get_data(int i) const;
	void set_data(int i, double val);

private:
	int size;
	double* data;
};


