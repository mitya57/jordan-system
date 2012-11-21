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
#include <pthread.h>
#include <sys/resource.h>

#include "block.hh"
#include "matrix.hh"
#include "matrix_solve.hh"
#include "function.hh"

void print_result(Matrix *matrix, double *rightcol, int limit) {
	for (int i = 0; i < (matrix->size > limit ? limit : matrix->size); ++i)
		std::cout << "variable " << i+1 << " is: "
		<< PRETTY(rightcol[i]) << std::endl;
}

double get_diff_norm(Matrix *matrix, double *rightcol) {
	double sum = 0;
	for (int i = 0; i < matrix->size; ++i)
		if (i & 1)
			sum += pow(rightcol[i], 2);
		else
			sum += pow(rightcol[i] - 1, 2);
	return sqrt(sum);
}

double get_residual(int size, double *oldc, double *newc) {
	double result = 0;
	for (int i = 0; i < size; ++i)
		result += (oldc[i]-newc[i])*(oldc[i]-newc[i]);
	return sqrt(result);
}

int main(int argc, char **argv) {
	Matrix *matrix = new Matrix, *origmatrix = new Matrix;
	int size, blocksize, threads, i, j;
	
	if (argc >= 4) {
		size = atoi(argv[1]);
		blocksize = atoi(argv[2]);
		threads = atoi(argv[3]);
	} else {
		std::cerr << "Usage: matrix SIZE BLOCKSIZE THREADS [FILENAME]"
		<< std::endl;
		exit(0);
	}
	matrix_new(matrix, size, blocksize);
	matrix_new(origmatrix, size, blocksize);
	double *rightcol = new double[size];
	double *realrightcol = new double[size];
	double *origrightcol = new double[size];
	if (argc < 5) {
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
		if (argc >= 4)
			matrixfile.open(argv[4]);
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
			matrix_set_element(origmatrix, i, j, matrix_get_element(matrix, i, j));
	print_matrix(matrix, rightcol);
	matrix_solve(size, blocksize, threads, matrix, rightcol);
	rusage resource_usage;
	getrusage(RUSAGE_SELF, &resource_usage);
	for (i = 0; i < size; ++i) {
		j = 0;
		while (matrix->blockcolumnind[j] != i/matrix->blocksize) ++j;
		realrightcol[i] = rightcol[j + i % matrix->blocksize];
	}
	print_result(matrix, realrightcol, 10);
	if (argc < 5)
		std::cout << "Norm of diff vector is "
		<< get_diff_norm(matrix, realrightcol) << std::endl;
	matrix_apply_to_vector(origmatrix, realrightcol, rightcol);
	std::cout << "Residual is " <<
	get_residual(size, origrightcol, rightcol) << std::endl;
	std::cout << "=============== TIME: "
	<< resource_usage.ru_utime.tv_sec << " seconds, "
	<< resource_usage.ru_utime.tv_usec / 1000
	<< " milliseconds ===============" << std::endl;
	delete[] rightcol;
	delete[] realrightcol;
	delete[] origrightcol;
	matrix_free(matrix);
	matrix_free(origmatrix);
	delete matrix;
	delete origmatrix;
	return 0;
}
