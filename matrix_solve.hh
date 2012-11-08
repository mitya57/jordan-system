/***************************************************
 * Solving Linear systems using the Jordan method  *
 *        Copyright: Dmitry Shachnev, 2012         *
 * This file contains the main solution procedure. *
 ***************************************************/

#ifndef NO_INCLUDE
#include <iostream>
#include <pthread.h>
#include "matrix.hh"
#include "block.hh"
#endif

static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

struct Arguments {
	Matrix *matrix;
	int i;
	int s;
	int threads;
	int *processedthreads;
	double *minnorm;
	double *mini;
	double *minj;
};

void *thread_function(void *argp) {
	Arguments *args = (Arguments *)argp;
	Matrix *matrix = args->matrix;
	
	double lminnorm = 1e20, norm;
	int lmini = args->s;
	int lminj = lmini;
	int l, t;
	
	double *lbuf = new double[matrix->blocksize * matrix->blocksize];
	double *rmatrix = new double[matrix->blocksize * matrix->blocksize];
	
	// On every line we want to process blocks:
	// t = s + i + k*threads    for each k
	
	for (l = args->s; l < matrix->size / matrix->blocksize; ++l)
		for (t = args->s + args->i; t < matrix->size / matrix->blocksize; t += args->threads) {
			// Process block (l, t)
			if (block_get_reverse(matrix->blocksize,
				matrix_get_pos_block(matrix, l, t),
				rmatrix,
				lbuf)
			) {
				norm = block_get_norm(matrix->blocksize, rmatrix);
				if (norm < lminnorm) {
					lminnorm = norm;
					lmini = l;
					lminj = t;
				}
			}
		}
	
	pthread_mutex_lock(&mutex);
	if (lminnorm < *(args->minnorm)) {
		*(args->minnorm) = lminnorm;
		*(args->mini) = lmini;
		*(args->minj) = lminj;
	}
	++(*(args->processedthreads));
	pthread_mutex_unlock(&mutex);
	return NULL;
}

void print_matrix(Matrix *matrix, double *rightcol) {
	int i, j;
	double el;
	for (i = 0; i < (matrix->size > 10 ? 10 : matrix->size); ++i) {
		for (j = 0; j < (matrix->size > 10 ? 10 : matrix->size); ++j) {
			el = matrix_get_element(matrix, i, j);
			el = double(int(el*1e4))/1e4;
			std::cout << el << "\t";
		}
		if (matrix->size > 10) std::cout << " ...";
		std::cout << "\t" << rightcol[i] << std::endl;
	}
	if (matrix->size > 10) std::cout << ". . . . ." << std::endl;
	std::cout << "-----------------------------" << std::endl;
}

void matrix_solve(int size, int blocksize, Matrix *matrix, double *rightcol) {
	int mini, minj, i, j, s, index;
	double *tempmatrix = new double[blocksize*blocksize];
	double *buf = new double[blocksize*blocksize];
	double *block, *tempvector = new double[blocksize];
	double minnorm, norm;
	
	for (s = 0; s < matrix->numberOfBlockColumns; ++s) {
		minnorm = 1e20;
		mini = s;
		minj = s;
		
		// TODO: calculate only norms of non-zero blocks
		for (i = s; i < size/blocksize; ++i)
			for (j = s; j < size/blocksize; ++j) {
				if (block_get_reverse(blocksize,
					matrix_get_pos_block(matrix, i, j),
					tempmatrix,
					buf)
				) {
				    norm = block_get_norm(blocksize, tempmatrix);
				    if (norm < minnorm) {
				    	minnorm = norm;
				    	mini = i;
				    	minj = j;
				    }
				}
			}
		matrix_swap_block_columns(matrix, minj, s);
		matrix_swap_block_rows   (matrix, mini, s);
		if (s != mini)
		for (j = 0; j < blocksize; ++j)
			std::swap(rightcol[blocksize*s+j],
				rightcol[blocksize*mini+j]);
		block_get_reverse(matrix_get_pos_height(matrix, s),
			matrix_get_pos_block(matrix, s, s),
			tempmatrix,
			buf
		);
		for (i = s+1; i < matrix->numberOfBlockColumns; ++i)
			block_left_multiply_alt(blocksize,
				matrix_get_pos_width(matrix, i),
				matrix_get_pos_block(matrix, s, i),
				tempmatrix,
				buf
			);
		block_apply_to_vector_alt(
			matrix_get_pos_height(matrix, s),
			tempmatrix,
			&(rightcol[blocksize*s]),
			buf
		);
		for (i = 0; i < matrix_get_pos_height(matrix, s); ++i)
			for (j = 0; j < matrix_get_pos_width(matrix, s); ++j)
				matrix_get_pos_block(matrix, s, s)
				[i*matrix_get_pos_height(matrix, s)+j] = (i == j);
		for (i = s+1; i < matrix->numberOfBlockColumns; ++i) {
			for (j = s+1; j < matrix->numberOfBlockColumns; ++j) {
				block_left_multiply(
					matrix_get_pos_height(matrix, i),
					matrix_get_pos_width(matrix, s), /* blocksize */
					matrix_get_pos_height(matrix, j),
					matrix_get_pos_block(matrix, s, j),
					matrix_get_pos_block(matrix, i, s),
					tempmatrix
				);
				block_subtract(
					matrix_get_pos_size(matrix, i, j),
					matrix_get_pos_block(matrix, i, j),
					tempmatrix
				);
			}
			block_apply_to_vector(
				blocksize,
				matrix_get_pos_width(matrix, i),
				matrix_get_pos_block(matrix, i, s),
				&(rightcol[s*blocksize]),
				tempvector
			);
			block_subtract(matrix_get_pos_height(matrix, i),
				&(rightcol[i*blocksize]),
				tempvector
			);
		}
		for (i = s+1; i < matrix->numberOfBlockColumns; ++i) {
			block = matrix_get_pos_block(matrix, i, s);
			index = blocksize * matrix_get_pos_height(matrix, i);
			for (j = 0; j < index; ++j)
				block[j] = 0;
		}
	}
	for (i = size-1; i >= 0; --i)
		for (j = i+1; j < size; ++j)
			rightcol[i] -= matrix_get_element(matrix, i, j)*rightcol[j];
	delete[] tempmatrix;
	delete[] buf;
	delete[] tempvector;
}
