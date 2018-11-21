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

int Vector::get_size() {
	return size;
}
double Vector::get_data(int i) {
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

double Vector::operator * (Vector& v2) {
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



void Vector::print_vector() {
	for (int i = 0; i < size; i++) {
		std::cout << "\n" << data[i];
	}
	std::cout << std::endl;
}
