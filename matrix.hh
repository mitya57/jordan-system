/****************************************************
 *  Solving Linear systems using the Jordan method  *
 *        Copyright: Dmitry Shachnev, 2012          *
 * This file contains the matrix-related functions. *
 ***************************************************/

#ifndef NO_INCLUDE
#include <algorithm>
#endif

// TODO: reduce size of numberOfBlock* arrays (?)

struct Matrix {
	double *elements;
	int *blockrowind, *blockcolumnind;
	int size, blocksize;
	int numberOfBlockRows, numberOfBlockColumns, numberOfBlocks;
};

void matrix_new(Matrix *matrix, int size, int blocksize) {
	matrix->size = size;
	matrix->elements = new double[size*size];
	matrix->blocksize = blocksize;
	matrix->numberOfBlockRows =
		(size % blocksize ? size / blocksize + 1 : size / blocksize);
	matrix->numberOfBlockColumns = matrix->numberOfBlockRows;
	matrix->numberOfBlocks = matrix->numberOfBlockRows *
		matrix->numberOfBlockColumns;
	matrix->blockrowind = new int[size];
	matrix->blockcolumnind = new int[size];
	for (int i=0; i<size; ++i) {
		matrix->blockrowind[i] = i / blocksize;
		matrix->blockcolumnind[i] = i / blocksize;
	}
}

void matrix_free(Matrix *matrix) {
	delete[] matrix->elements;
	delete[] matrix->blockrowind;
	delete[] matrix->blockcolumnind;
}

int matrix_get_block_width(Matrix *matrix, int index) {
	if (index % matrix->numberOfBlockColumns < matrix->size / matrix->blocksize)
		return matrix->blocksize;
	return matrix->size % matrix->blocksize;
}

int matrix_get_pos_width(Matrix *matrix, int col) {
	if (col < matrix->size / matrix->blocksize)
		return matrix->blocksize;
	return matrix->size % matrix->blocksize;
}

int matrix_get_block_height(Matrix *matrix, int index) {
	if (index / matrix->numberOfBlockColumns < matrix->size / matrix->blocksize)
		return matrix->blocksize;
	return matrix->size % matrix->blocksize;
}

int matrix_get_pos_height(Matrix *matrix, int row) {
	if (row < matrix->size / matrix->blocksize)
		return matrix->blocksize;
	return matrix->size % matrix->blocksize;
}

int matrix_get_pos_size(Matrix *matrix, int row, int col) {
	return matrix_get_pos_height(matrix, row) * matrix_get_pos_width(matrix, col);
}

int matrix_get_block_index(Matrix *matrix, int x, int y) {
	return (matrix->blockrowind[x])*(matrix->numberOfBlockRows) + matrix->blockcolumnind[y];
}

double *matrix_get_block(Matrix *matrix, int index) {
	int start = matrix->blocksize * matrix->size *
		(index / matrix->numberOfBlockColumns);
	start += matrix_get_block_height(matrix, index) *
		matrix->blocksize * (index % matrix->numberOfBlockColumns);
	return &((matrix->elements)[start]);
}

double *matrix_get_pos_block(Matrix *matrix, int row, int col) {
	int index = matrix_get_block_index(matrix, row*matrix->blocksize, col*matrix->blocksize);
	return matrix_get_block(matrix, index);
}

int matrix_get_element_index_in_block(Matrix *matrix, int x, int y) {
	int index = matrix_get_block_index(matrix, x, y);
	return (x % matrix->blocksize)*(matrix_get_block_width(matrix, index))
	+ (y % matrix->blocksize);
}

double matrix_get_element(Matrix *matrix, int x, int y) {
	int index = matrix_get_block_index(matrix, x, y);
	return matrix_get_block(matrix, index)[matrix_get_element_index_in_block(matrix, x, y)];
}

void matrix_set_element(Matrix *matrix, int x, int y, double value) {
	int index = matrix_get_block_index(matrix, x, y);
	matrix_get_block(matrix, index)[matrix_get_element_index_in_block(matrix, x, y)]
		= value;
}

void matrix_swap_block_columns(Matrix *matrix, int y1, int y2) {
	if (y1 == y2) return;
	if (y1 > y2)
		std::swap(y1, y2);
	int   start1 = y1 * matrix->blocksize,
	end1, start2 = y2 * matrix->blocksize;
	end1 = start1;
	while (matrix->blockcolumnind[end1] == matrix->blockcolumnind[start1]) ++end1;
	while (start1 < end1) {
		std::swap(matrix->blockcolumnind[start1], matrix->blockcolumnind[start2]);
		++start1;
		++start2;
	}
}

void matrix_swap_block_rows(Matrix *matrix, int x1, int x2) {
	if (x1 == x2) return;
	if (x1 > x2)
		std::swap(x1, x2);
	int   start1 = x1 * matrix->blocksize,
	end1, start2 = x2 * matrix->blocksize;
	end1 = start1;
	while (matrix->blockrowind[end1] == matrix->blockrowind[start1]) ++end1;
	while (start1 < end1) {
		std::swap(matrix->blockrowind[start1], matrix->blockrowind[start2]);
		++start1;
		++start2;
	}
}
