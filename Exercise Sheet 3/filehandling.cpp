#include <iostream>
#include <string>
#include <fstream>

//NOTE: it is best for filename to be the complete path to the file, in path use '/' instead of '\'. 
//(Example:"C:/Users/ASUS/Documents/Studium/Prog Praktikum/Exercise Sheet 3/testfile.txt")
//the file should be a .txt file, if it doesnt exist yet, it will be created automatically

//parameters:
	//filename: name of file (or path to file)
	//size: size of the arrays
	//array1: array of entries that will be written into the left hand column
	//array2: array of entries that will be written into the right hand column
//creates file of format:
	//line 1:	size
	//line 2:	  array1[0] array2[0]
	//line 3:	  array1[1] array2[1]
	//...
	//line size+1:array1[size-1] array2[size-1]
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


//parameters:
	//filename: name of file (or path to file)
	//v1: Vector of entries that will be written into the left hand column
	//v2: Vector of entries that will be written into the right hand column
	//vectors have to have same size
//creates file of format:
	//line 1:	     size
	//line 2:	     v1.data[o] v2.data[0]
	//line 3:	     v1.data[1] v2.data[1]
	//...
	//line v1.size+1:v1.data[v1.size-1] v2.data[v2.size-1]
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
