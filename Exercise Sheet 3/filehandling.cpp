#include <iostream>
#include <string>
#include <fstream>

//NOTE: it is best for filename to be the complete path to the file, in path use '/' instead of '\'
//the file should be a txt file, if it doesnt exist yet, it will be created automatically

void write_value_pairs_to_file_array(std::string filename, int size, double* array1, double* array2) {
	if (size < 0) {
		std::cout << "\nwrite_value_pairs_to_file failed: size has to be >=0\n";
		return;
	}
	std::ofstream myFile;
	myFile.open(filename);
	if (myFile.is_open() == false) {
		std::cout << "\nwrite_value_pairs_to_file failed: couldnt open file\n";
		return;
	}
	myFile << size << "\n";
	for (int i = 0; i < size - 1; i++) {
		myFile << array1[i] << " " << array2[i] << "\n";
	}
	myFile << array1[size - 1] << " " << array2[size - 1];
	myFile.close();
}

void write_value_pairs_to_file_vec(std::string filename, const Vector& v1, const Vector& v2) {
	if (v1.get_size() != v2.get_size()) {
		std::cout << "\nwrite_value_pairs_to_file_vect failed: vectors have to have same size\n";
		return;
	}
	std::ofstream myFile;
	myFile.open(filename);
	if (myFile.is_open() == false) {
		std::cout << "\nwrite_value_pairs_to_file failed: couldnt open file\n";
		return;
	}
	myFile << v1.get_size() << "\n";
	for (int i = 0; i < v1.get_size() - 1; i++) {
		myFile << v1.get_data(i) << " " << v2.get_data(i) << "\n";
	}
	myFile << v1.get_data(v1.get_size() - 1) << " " << v2.get_data(v2.get_size() - 1);
	myFile.close();
}
