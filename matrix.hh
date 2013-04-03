/****************************************************
 *  Solving linear systems using the Jordan method  *
 *          Author: Dmitry Shachnev, 2013           *
 * This file contains the matrix-related functions. *
 ***************************************************/

// TODO: reduce size of numberOfBlock* arrays (?)

#define numberOfBlockColumns numberOfBlockRows

struct MatrixInfo {
	int *blockRowInd, *blockColumnInd;
	int size, blocksize;
	int numberOfBlockRows;
};

void matrix_clear(MatrixInfo *matrix) {
	for (int x = 0; x < matrix->numberOfBlockColumns; ++x) {
		matrix->blockRowInd[x] = x;
		matrix->blockColumnInd[x] = x;
	}
}

void matrix_new(MatrixInfo *matrix, int size, int blocksize) {
	matrix->size = size;
	matrix->blocksize = blocksize;
	matrix->numberOfBlockRows =
		(size % blocksize ? size / blocksize + 1 : size / blocksize);
	matrix->blockRowInd = new int[matrix->numberOfBlockRows];
	matrix->blockColumnInd = new int[matrix->numberOfBlockColumns];
	matrix_clear(matrix);
}

void matrix_free(MatrixInfo *matrix) {
	delete[] matrix->blockRowInd;
	delete[] matrix->blockColumnInd;
}

void matrix_swap_block_columns(MatrixInfo *matrix, int y1, int y2) {
	if (y1 == y2)
		return;
	std::swap(matrix->blockColumnInd[y1], matrix->blockColumnInd[y2]);
}

void matrix_swap_block_rows(MatrixInfo *matrix, int x1, int x2) {
	if (x1 == x2)
		return;
	std::swap(matrix->blockRowInd[x1], matrix->blockRowInd[x2]);
}
