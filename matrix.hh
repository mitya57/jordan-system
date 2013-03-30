/****************************************************
 *  Solving linear systems using the Jordan method  *
 *          Author: Dmitry Shachnev, 2013           *
 * This file contains the matrix-related functions. *
 ***************************************************/

// TODO: reduce size of numberOfBlock* arrays (?)

#define numberOfBlockColumns numberOfBlockRows

struct MatrixInfo {
	int *blockrowind, *blockcolumnind;
	int size, blocksize;
	int numberOfBlockRows;
};

void matrix_clear(MatrixInfo *matrix) {
	for (int i = 0; i < matrix->size; ++i) {
		matrix->blockrowind[i] = i / matrix->blocksize;
		matrix->blockcolumnind[i] = i / matrix->blocksize;
	}
}

void matrix_new(MatrixInfo *matrix, int size, int blocksize) {
	matrix->size = size;
	matrix->blocksize = blocksize;
	matrix->numberOfBlockRows =
		(size % blocksize ? size / blocksize + 1 : size / blocksize);
	matrix->blockrowind = new int[size];
	matrix->blockcolumnind = new int[size];
	matrix_clear(matrix);
}

void matrix_free(MatrixInfo *matrix) {
	delete[] matrix->blockrowind;
	delete[] matrix->blockcolumnind;
}

void matrix_swap_block_columns(MatrixInfo *matrix, int y1, int y2) {
	if (y1 == y2)
		return;
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

void matrix_swap_block_rows(MatrixInfo *matrix, int x1, int x2) {
	if (x1 == x2)
		return;
	if (x1 > x2)
		std::swap(x1, x2);
	int start1 = x1 * matrix->blocksize,
	end1, start2 = x2 * matrix->blocksize;
	end1 = start1;
	while (matrix->blockrowind[end1] == matrix->blockrowind[start1]) ++end1;
	while (start1 < end1) {
		std::swap(matrix->blockrowind[start1], matrix->blockrowind[start2]);
		++start1;
		++start2;
	}
}
