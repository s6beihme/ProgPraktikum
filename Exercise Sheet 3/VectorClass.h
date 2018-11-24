#pragma once
#include <iostream>
#include <memory>
#include <string>
class Vector {
public:
	//allocates memory for Vector type object and fills data with zeros
	//parameters:
		//_size: size of vector;
	Vector(int _size);

	//creates vector from given file
	//parameters:
		//filename: file representing vector in format:
			//line 1:	  size
			//line 2:	  data
			//...
			//line size+1:data
	Vector(std::string filename);
	~Vector();

	//fills this->data with given values (Note that this function doesnt check for correct size of argument)
	//parameters:
		//vals: array containing the values that should be stored in vector
	void vec_assemble(double* vals);

	//dot product of two vectors:
	//parameters:
		//v2: Vector with whicht to multiply (*this)
	//returns: double value of dot product v2*(*this)
	double operator * (const Vector& v2) const;

	//adds two vectors (doesnt change (*this))
	//parameters:
		//v2: vector to add with (*this)
		//result: data of result gets overwritten with added data of (*this) and v2
	void add_vect(const Vector& v2, Vector& result) const;

	//subtracts two vectors (doesnt change (*this))
	//parameters:
		//v2: vector to subtract (*this)
		//result: data of result gets overwritten with subtracted data of (*this) and v2
	void subtr_vect(const Vector& v2, Vector& result) const;

	//prints data of vector
	void print_vector() const;

	int get_size() const;
	double get_data(int i) const;
	void set_data(int i, double val);

	//copys data of (*this) into v2;
	void make_copy(Vector& v2) const;
	long double norm_squared() const;

private:
	int size;
	//double* data;
	std::unique_ptr<double[]> data;
};

void scalar_mult(double a, const Vector& v, Vector& result);


