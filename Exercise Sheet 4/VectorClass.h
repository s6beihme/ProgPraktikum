#pragma once
#include <iostream>
#include <memory>
#include <string>
#include <fstream>

template <typename T>
class Vector {
public:
	Vector();
	Vector(int _size);
	Vector(std::string filename);
	Vector(const Vector& other);
	Vector& operator = (const Vector& other);

	void assemble(double* vals);

	bool operator ==(const Vector& other);
	bool operator !=(const Vector& other);
	T& operator [] (int i);
	T operator *(const Vector& other);
	T norm_squared();
	
	int get_size() const;
	
	template <typename T2>
	friend void add_vect (const Vector<T2>& v1, const Vector<T2>& v2, Vector<T2>& result);
	template <typename T2>
	friend void subtr_vect (const Vector<T2>& v1, const Vector<T2>& v2, Vector<T2>& result);
	template <typename T2>
	friend std::ostream& operator << (std::ostream& out, const Vector<T2>& v);

private:
	int size;
	std::unique_ptr<T[]> data;
};

// default construction
template <typename T>
Vector<T>::Vector() :
	size(0),
	data(nullptr)
{}

//allocates memory for Vector type object and fills data with zeros
	//parameters:
		//_size: size of vector;
template <typename T>
Vector<T>::Vector(int _size) :
	size(_size)
{
	data = std::make_unique<T[]>(size);
	for (int i = 0; i < size; i++) {
		data[i] = 0;
	}
}

//creates vector from given file
	//parameters:
		//filename: file representing vector in format:
			//line 1:	  size
			//line 2:	  data
			//...
			//line size+1:data
template <typename T>
Vector<T>::Vector(std::string filename) {
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
	data = std::make_unique<T[]>(size);
	for (int i = 0; i < size; i++) {
		myFile >> data[i];
	}
	myFile.close();
}

//copy construction
template <typename T>
Vector<T>::Vector(const Vector<T>& other) :
	size(other.size)
{
	data = std::make_unique<T[]>(size);
	for (int i = 0; i < size; i++) {
		data[i] = other.data[i];
	}
}

//copy assignment
template <typename T>
Vector<T>& Vector<T>::operator = (const Vector<T>& other) {
	size = other.size;
	data = std::make_unique<T[]>(size);
	for (int i = 0; i < size; i++) {
		data[i] = other.data[i];
	}
	return *this;
}

//fills this->data with given values (Note that this function doesnt check for correct size of argument)
	//parameters:
		//vals: array containing the values that should be stored in vector
template <typename T>
void Vector<T>::assemble(double* vals) {
	for (int i = 0; i < size; i++) data[i] = vals[i];
}

// equality
template <typename T>
bool Vector<T>::operator ==(const Vector<T>& other) {
	if (size != other.size) return false;

	for (int i = 0; i < size; i++) {
		if (data[i] != other.data[i]) return false; //this ok for doubles?
	}
	return true;
}

//inequality
template <typename T>
bool Vector<T>::operator !=(const Vector<T>& other) {
	return !((*this) == other);
}

template <typename T>
T& Vector<T>::operator [] (int i) {
	if (i<0 || i>size-1) { std::cout << "\nVector operator[] failed index out of bounds"; exit(0); }
	return data[i];
}

//dot product of two vectors:
	//parameters:
		//v2: Vector with whicht to multiply (*this)
	//returns: double value of dot product v2*(*this)
template <typename T>
T Vector<T>::operator *(const Vector<T>& other) {
	if (size != other.size) { std::cout << "\nvector multiplication failed: sizes didnt correspond"; exit(0); }
	T sum = 0;
	for (int i = 0; i < size; i++) sum += data[i] * other.data[i];
	return sum;
}

//squared norm of a vector
template <typename T>
T Vector<T>::norm_squared() {
	T result = 0;
	for (int i = 0; i < size; i++) {
		result += data[i] * data[i];
	}
	return result;
}

template <typename T>
int Vector<T>::get_size() const {
	return size;
}

//adds two vectors (doesnt change (*this))
	//parameters:
		//v2: vector to add with (*this)
		//result: data of result gets overwritten with added data of (*this) and v2
template <typename T>
void add_vect(const Vector<T>& v1, const Vector<T>& v2, Vector<T>& result) {
	if (v1.size != v2.size || v2.size != result.size) { std::cout << "\nadd_vects failed: vector sizes didnt correspond"; exit(0); }
	for (int i = 0; i < v1.size; i++) result.data[i] = v1.data[i] + v1.data[i];
}

//subtracts two vectors (doesnt change (*this))
	//parameters:
		//v2: vector to subtract (*this)
		//result: data of result gets overwritten with subtracted data of (*this) and v2
template <typename T>
void subtr_vect(const Vector<T>& v1, const Vector<T>& v2, Vector<T>& result) {
	if (v1.size != v2.size || v2.size != result.size) { std::cout << "\nadd_vects failed: vector sizes didnt correspond"; exit(0); }
	for (int i = 0; i < v1.size; i++) result.data[i] = v1.data[i] - v1.data[i];
}

template <typename T>
std::ostream& operator << (std::ostream& out, const Vector<T>& v) { 
	out << "Vector of size " << v.size << ":\n";
	for (int i = 0; i < v.size; i++) {
		out << v.data[i] << ", ";
	}
	out << "\n";
	return out;
}
