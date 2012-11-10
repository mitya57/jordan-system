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

#include <cassert>

static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

struct Arguments {
	Matrix *matrix;
	int i;
	int ti;
	int s;
	int threads;
	int *processedthreads;
	double *minnorm;
	double *rightcol;
	int *mini;
	int *minj;
};

void *tf1(void *argp) {
	Arguments *args = (Arguments *)argp;
	Matrix *matrix = args->matrix;
	
	double lminnorm = -1, norm;
	int lmini = args->s;
	int lminj = lmini;
	int l, t;
	
	double *lbuf = new double[matrix->blocksize * matrix->blocksize];
	double *rmatrix = new double[matrix->blocksize * matrix->blocksize];
	
	// On every line we want to process blocks:
	// t = s + i + k*threads    for each k
	
	for (l = args->s; l < matrix->size / matrix->blocksize; ++l)
		for (t = args->s + args->ti; t < matrix->size / matrix->blocksize; t += args->threads) {
			// Process block (l, t)
			if (block_get_reverse(matrix->blocksize,
				matrix_get_pos_block(matrix, l, t),
				rmatrix,
				lbuf)
			) {
				norm = block_get_norm(matrix->blocksize, rmatrix);
				if (lminnorm < 0 || norm < lminnorm) {
					lminnorm = norm;
					lmini = l;
					lminj = t;
				}
			}
		}
	
	pthread_mutex_lock(&mutex);
	if (!*(args->processedthreads) || lminnorm < *(args->minnorm)) {
		*(args->minnorm) = lminnorm;
		*(args->mini) = lmini;
		*(args->minj) = lminj;
	}
	++(*(args->processedthreads));
	pthread_mutex_unlock(&mutex);
	delete[] lbuf;
	delete[] rmatrix;
	return NULL;
}

void *tf2(void *argp) {
	Arguments *args = (Arguments *)argp;
	int i = args->i, j, s = args->s;
	Matrix *matrix = args->matrix;
	double *tmp = new double [matrix->blocksize*matrix->blocksize];
	for (j = s + 1 + args->ti; j < matrix->numberOfBlockColumns; j += args->threads) {
#ifdef DEBUG
		pthread_mutex_lock(&mutex);
		std::cout << "thread " << args->ti << ": s = " << s << ", i := " << i << ", j := " << j << std::endl;
		pthread_mutex_unlock(&mutex);
#endif
		block_left_multiply(
			matrix_get_pos_height(matrix, i),
			matrix_get_pos_width(matrix, s),
			matrix_get_pos_height(matrix, j),
			matrix_get_pos_block(matrix, s, j),
			matrix_get_pos_block(matrix, i, s),
			tmp
		);
		block_subtract(
			matrix_get_pos_size(matrix, i, j),
			matrix_get_pos_block(matrix, i, j),
			tmp
		);
	}
	pthread_mutex_lock(&mutex);
	++(*(args->processedthreads));
	pthread_mutex_unlock(&mutex);
	delete[] tmp;
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

void matrix_solve(int size, int blocksize, int threads, Matrix *matrix, double *rightcol) {
	int mini, minj, i, j, s, t, index, processedthreads, retval;
	double *tempmatrix = new double[blocksize*blocksize];
	double *buf = new double[blocksize*blocksize];
	double *block;
	double numberofthreads;
	double minnorm;
	
	pthread_t *thr = new pthread_t[threads];
	Arguments *args = new Arguments[threads];
	
	for (s = 0; s < matrix->numberOfBlockColumns; ++s) {
#ifdef DEBUG
		std::cout << "s: " << s << std::endl;
		print_matrix(matrix, rightcol);
#endif
		if (s + 1 != matrix->numberOfBlockColumns) {
			processedthreads = 0;
			numberofthreads = threads;
			if (numberofthreads > size/blocksize - s)
				numberofthreads = size/blocksize - s;
			for (t = 0; t < numberofthreads; ++t) {
				args[t].matrix = matrix;
				args[t].ti = t;
				args[t].s = s;
				args[t].threads = numberofthreads;
				args[t].processedthreads = &processedthreads;
				args[t].minnorm = &minnorm;
				args[t].mini = &mini;
				args[t].minj = &minj;
				
				int res = pthread_create(
					&(thr[t]), // Thread identifier
					NULL,      // Thread attributes: using defaults
					&tf1,      // Thread start function
					&(args[t]) // Parameter to be passed to thread function
				);
				if (res) std::cerr << "Cannot create thread!" << std::endl;
			}
			
			for (t = 0; t < numberofthreads; ++t)
				pthread_join(thr[t], NULL);
			
			assert (processedthreads == numberofthreads);
#ifdef DEBUG
			std::cout << "minblock: (" << mini << ", " << minj << ")" << std::endl;
#endif
			matrix_swap_block_columns(matrix, minj, s);
			matrix_swap_block_rows   (matrix, mini, s);
			if (s != mini)
				for (j = 0; j < blocksize; ++j)
					std::swap(rightcol[blocksize*s+j],
						rightcol[blocksize*mini+j]);
		}
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
#ifdef DEBUG
		std::cout << "Before calling tf2:" << std::endl;
		print_matrix(matrix, rightcol);
#endif
		numberofthreads = threads;
		if (numberofthreads > matrix->numberOfBlockColumns - s)
			numberofthreads = matrix->numberOfBlockColumns - s;
		for (i = s+1; i < matrix->numberOfBlockColumns; ++i) {
			processedthreads = 0;
			for (t = 0; t < numberofthreads; ++t) {
				//args[t].matrix = matrix;
				//args[t].ti = t;
				//args[t].s = s;
				args[t].threads = numberofthreads;
				args[t].processedthreads = &processedthreads;
				args[t].ti = t;
				args[t].i = i;
				args[t].rightcol = rightcol;
				int res = pthread_create(
					&(thr[t]), // Thread identifier
					NULL,      // Thread attributes: using defaults
					&tf2,      // Thread start function
					&(args[t]) // Parameter to be passed to thread function
				);
				if (res) std::cerr << "Cannot create thread!" << std::endl;
			}
			for (t = 0; t < numberofthreads; ++t)
				pthread_join(thr[t], NULL);
			assert (processedthreads == numberofthreads);
			block_apply_to_vector(
				matrix->blocksize,
				matrix_get_pos_width(matrix, i),
				matrix_get_pos_block(matrix, i, s),
				&(rightcol[s * matrix->blocksize]),
				tempmatrix
			);
			block_subtract(matrix_get_pos_height(matrix, i),
				&(rightcol[i * matrix->blocksize]),
				tempmatrix
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
	pthread_mutex_destroy(&mutex);
	delete[] tempmatrix;
	delete[] buf;
	delete[] thr;
	delete[] args;
}
