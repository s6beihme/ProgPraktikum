#pragma once
#include "VectorClass.h"
#include <string>
void write_value_pairs_to_file_array(std::string filename, int size, double* array1, double* array2);
void write_value_pairs_to_file_vec(std::string filename, const Vector& v1, const Vector& v2);
int number_of_lines_in_file(std::string filename);