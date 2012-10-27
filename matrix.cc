/**************************************************
 * Solving Linear systems using the Jordan method *
 *        Copyright: Dmitry Shachnev, 2012        *
 *      This file contains the main program.      *
 **************************************************/

#define NO_INCLUDE
#define EPS 1e-8
#define PRETTY(x) (((x) < EPS && (x) > -EPS) ? 0 : x)

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <sys/resource.h>
#include "block.hh"
#include "matrix.hh"
#include "matrix_solve.hh"
#include "function.hh"

double get_result(Matrix *matrix, double *rightcol, int i) {
	int j = 0;
	while (matrix->blockcolumnind[j] != i/matrix->blocksize) ++j;
	return rightcol[j + i % matrix->blocksize];
}

void print_result(Matrix *matrix, double *rightcol, int limit) {
	for (int i = 0; i < (matrix->size > limit ? limit : matrix->size); ++i)
		std::cout << "variable " << i+1 << " is: "
		<< PRETTY(get_result(matrix, rightcol, i)) << std::endl;
}

void matrix_apply_to_vector(Matrix *matrix, double *origmatrix,
double *vector, double *newvector) {
	int i, j;
	for (i = 0; i < matrix->size; ++i) {
		newvector[i] = 0;
		for (j = 0; j < matrix->size; ++j)
			newvector[i] += origmatrix[i*(matrix->size)+j] * get_result(matrix, vector, j);
	}
}

double get_diff_norm(Matrix *matrix, double *rightcol) {
	double sum = 0;
	for (int i = 0; i < matrix->size; ++i)
		if (i & 1)
			sum += pow(get_result(matrix, rightcol, i), 2);
		else
			sum += pow(get_result(matrix, rightcol, i) - 1, 2);
	return sqrt(sum);
}

double get_residual(int size, double *oldc, double *newc) {
	double result = 0;
	for (int i = 0; i < size; ++i)
		result += (oldc[i]-newc[i])*(oldc[i]-newc[i]);
	return sqrt(result);
}

int main(int argc, char **argv) {
	Matrix *matrix = new Matrix;
	int size, blocksize, i, j;
	
	if (argc >= 2) {
		size = atoi(argv[1]);
		blocksize = atoi(argv[2]);
	} else {
		std::cout << "Enter size and blocksize: ";
		std::cin >> size;
		std::cin >> blocksize;
	}
	matrix_new(matrix, size, blocksize);
	double *rightcol = new double[size];
	double *origmatrix = new double[size*size];
	double *origrightcol = new double[size];
	double *testrightcol = new double[size];
	if (argc < 4) {
		for (i = 0; i < size; ++i)
			for (j = 0; j < size; ++j)
				matrix_set_element(matrix, i, j, get_matrix_element(i, j));
		for (i = 0; i < size; ++i) {
			rightcol[i] = 0;
			for (j = 0; j < size; j += 2)
				rightcol[i] += matrix_get_element(matrix, i, j);
		}
	} else { 
		std::fstream matrixfile;
		if (argc >= 3)
			matrixfile.open(argv[3]);
		else
			matrixfile.open("matrix.txt");

		double buffer;
		matrixfile >> size;
		for (i = 0; i < size; ++i) {
			for (j = 0; j < size; ++j) {
				matrixfile >> buffer;
				matrix_set_element(matrix, i, j, buffer);
			}
			matrixfile >> rightcol[i];
		}
		matrixfile.close();
	}
	for (i = 0; i < size; ++i)
		origrightcol[i] = rightcol[i];
	for (i = 0; i < size; ++i)
		for (j = 0; j < size; ++j)
			origmatrix[i*size+j] = matrix_get_element(matrix, i, j);
	print_matrix(matrix, rightcol);
	matrix_solve(size, blocksize, matrix, rightcol);
	rusage resource_usage;
	getrusage(RUSAGE_SELF, &resource_usage);
	print_result(matrix, rightcol, 10);
	if (argc < 4)
		std::cout << "Norm of diff vector is "
		<< get_diff_norm(matrix, rightcol) << std::endl;
	matrix_apply_to_vector(matrix, origmatrix, rightcol, testrightcol);
	std::cout << "Residual is " <<
	get_residual(size, origrightcol, testrightcol) << std::endl;
	delete[] origmatrix;
	delete[] origrightcol;
	delete[] testrightcol;
	std::cout << "=============== TIME: "
	<< resource_usage.ru_utime.tv_sec << " seconds, "
	<< resource_usage.ru_utime.tv_usec / 1000
	<< " milliseconds ===============" << std::endl;
	matrix_free(matrix);
	delete[] rightcol;
	delete matrix;
	return 0;
}
