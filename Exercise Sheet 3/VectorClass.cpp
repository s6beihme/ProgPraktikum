#include "VectorClass.h"
#include <iostream>


Vector::Vector(int _size) {
	size = _size;
	data = new double[size];
	if (data == NULL) {
		std::cout << "Vector Constructor failed. allocating memory for data failed" << std::endl;
		return;
	}
}

Vector::~Vector() {
	delete[] data;
}

int Vector::get_size() const {
	return size;
}
double Vector::get_data(int i) const {
	return data[i];
}
void Vector::set_data(int i, double val) {
	if (i < 0 || i >= size) {
		std::cout << "\nset_data failed, index out of bound\n";
		return;
	}
	data[i] = val;
}

void Vector::vec_assemble(double* vals) {
	for (int i = 0; i < size; i++) {
		data[i] = vals[i];
	}
}

double Vector::operator * (const Vector& v2) const {
	if (size != v2.get_size()) {
		std::cout << "vector multiplication failed, because vectors dont have same size" << std::endl;
		exit(0); //maybe return is better
	}

	double sum = 0;
	for (int i = 0; i < size; i++) {
		sum += data[i] * v2.get_data(i);
	}
	return sum;
}

Vector Vector::operator - (const Vector& v2) const {
	if (size != v2.get_size()) {
		std::cout << "\ntrying to subtract vectors of different size\n";
		exit(0);
	}
	Vector res(size);
	double* res_data = new double[size];
	if (res_data == NULL) {
		std::cout << "\nin operator : failed to allocate memory for entry array of result vector\n";
		exit(0);
	}
	for (int i = 0; i < size; i++) {
		res_data[i] = data[i] - v2.get_data(i);
	}
	res.vec_assemble(res_data);
	delete[] res_data;
	return res;
}


void Vector::print_vector() const {
	for (int i = 0; i < size; i++) {
		std::cout << "\n" << data[i];
	}
	std::cout << std::endl;
}
