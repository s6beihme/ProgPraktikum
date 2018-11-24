#include "VectorClass.h"
#include <iostream>
#include <fstream>
#include <string>

Vector::Vector(int _size) {
	size = _size;
	data = std::make_unique<double[]>(size);
	for (int i = 0; i < size; i++) {
		data[i] = 0;
	}
}

Vector::Vector(std::string filename) {
	std::ifstream myFile;
	myFile.open(filename);
	if (myFile.is_open() == false) {
		std::cout << "\nfailed to open file!\n";
		exit(0);
	}
	
	myFile >> size;
	if (size < 0) {
		std::cout << "first entry of file representing vector is the size. it has to be >= 0";
		exit(0);
	}
	data = std::make_unique<double[]>(size);
	for (int i = 0; i < size; i++) {
		myFile >> data[i];
	}
	myFile.close();
}

Vector::~Vector() {
	data.release();
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

void Vector::add_vect(const Vector& v2, Vector& result) const {
	if (size != v2.get_size()) {
		std::cout << "\ntrying to subtract vectors of different size\n";
		return;
	}
	for (int i = 0; i < size; i++) {
		result.set_data(i, data[i] + v2.get_data(i));
	}
}

void Vector::subtr_vect (const Vector& v2, Vector& result) const {
	if (size != v2.get_size()) {
		std::cout << "\ntrying to subtract vectors of different size\n";
		exit(0);
	}
	
	for (int i = 0; i < size; i++) {
		result.set_data(i, data[i] - v2.get_data(i));
	}
}


void Vector::print_vector() const {
	std::cout << "\nEntries of vector are:\n";
	for (int i = 0; i < size; i++) {
		std::cout << "\n" << data[i];
	}
	std::cout << std::endl;
}

long double Vector::norm_squared() const {
	long double res = 0;
	for (int i = 0; i < size; i++) {
		res += data[i] * data[i];
	}
	return res;
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

void Vector::make_copy(Vector& v2) const {
	if (size != v2.get_size()) {
		std::cout << "\nmake copy failed, because size of argument didnt fit\n";
		return;
	}
	for (int i = 0; i < size; i++) {
		v2.set_data(i, data[i]);
	}
}

void scalar_mult(double a, const Vector& v, Vector& result) {

	for (int i = 0; i < v.get_size(); i++) {
		result.set_data(i, a * v.get_data(i));
	}
}

