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
static pthread_cond_t condvar_tomf = PTHREAD_COND_INITIALIZER;
static pthread_cond_t condvar_totf = PTHREAD_COND_INITIALIZER;

static bool threads = false;

struct Arguments {
	Matrix *matrix;
	int ti;
	int threads;
	int *processedthreads;
	double *minnorm;
	double *rightcol;
	int *mini;
	int *minj;
};

void *tf(void *argp) {
	Arguments *args = (Arguments *)argp;
	Matrix *matrix = args->matrix;
	
	double *tmp = new double[matrix->blocksize * matrix->blocksize];
	double *rmatrix = new double[matrix->blocksize * matrix->blocksize];
	
	double lminnorm, norm;
	int lmini;
	int lminj;
	int s, i, j;

	pthread_mutex_lock(&mutex);
	++(*(args->processedthreads));
	if (*(args->processedthreads) == args->threads)
		pthread_cond_signal(&condvar_tomf);
	pthread_cond_wait(&condvar_totf, &mutex);
	pthread_mutex_unlock(&mutex);

	for (s = 0; s < matrix->numberOfBlockColumns; ++s) {
		// Part 1: selecting main block
		// ----------------------------
		lminnorm = -1;
		
		for (i = s + args->ti; i < matrix->size / matrix->blocksize; i += args->threads) {
			for (j = s; j < matrix->size / matrix->blocksize; ++j)
				// Process block (i, j)
				if (block_get_reverse(matrix->blocksize,
					matrix_get_pos_block(matrix, i, j),
					rmatrix,
					tmp)
				) {
					norm = block_get_norm(matrix->blocksize, rmatrix);
					if (lminnorm < 0 || norm < lminnorm) {
						lminnorm = norm;
						lmini = i;
						lminj = j;
					}
				}
		}
		pthread_mutex_lock(&mutex);
		if (lminnorm >= 0 && (*(args->minnorm) < 0 || lminnorm < *(args->minnorm))) {
			*(args->minnorm) = lminnorm;
			*(args->mini) = lmini;
			*(args->minj) = lminj;
		}
		++(*(args->processedthreads));
		if (*(args->processedthreads) == args->threads)
			pthread_cond_signal(&condvar_tomf);

		// Part 2: subtracting block rows
		// ------------------------------
		pthread_cond_wait(&condvar_totf, &mutex);
		pthread_mutex_unlock(&mutex);
		for (i = s + args->ti + 1; i < matrix->numberOfBlockRows; i += args->threads) {
			for (j = s+1; j < matrix->numberOfBlockColumns; ++j) {
				block_left_multiply(
					matrix_get_pos_height(matrix, i),
					matrix_get_pos_width(matrix, s), /* blocksize */
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
			block_apply_to_vector(
				matrix->blocksize,
				matrix_get_pos_width(matrix, i),
				matrix_get_pos_block(matrix, i, s),
				&(args->rightcol[s * matrix->blocksize]),
				tmp
			);
			block_subtract(matrix_get_pos_height(matrix, i),
				&(args->rightcol[i * matrix->blocksize]),
				tmp
			);
		}
		pthread_mutex_lock(&mutex);
		++(*(args->processedthreads));
		if (*(args->processedthreads) == args->threads)
			pthread_cond_signal(&condvar_tomf);
		if (s + 1 != matrix->numberOfBlockColumns)
			pthread_cond_wait(&condvar_totf, &mutex);
		pthread_mutex_unlock(&mutex);
	}
	delete[] tmp;
	delete[] rmatrix;
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
	int mini, minj, i, j, s, t, index, processedthreads = 0;
	double *tempmatrix = new double[blocksize*blocksize];
	double *buf = new double[blocksize*blocksize];
	double *block;
	double minnorm;
	
	pthread_t *thr = new pthread_t[threads];
	Arguments *args = new Arguments[threads];
	
	pthread_mutex_lock(&mutex);
	for (t = 0; t < threads; ++t) {
		args[t].matrix = matrix;
		args[t].ti = t;
		args[t].threads = threads;
		args[t].rightcol = rightcol;
		args[t].processedthreads = &processedthreads;
		args[t].minnorm = &minnorm;
		args[t].mini = &mini;
		args[t].minj = &minj;
		int res = pthread_create(
			&(thr[t]), // Thread identifier
			NULL,      // Thread attributes: using defaults
			&tf,       // Thread start function
			&(args[t]) // Parameter to be passed to thread function
		);
		if (res) std::cerr << "Cannot create thread!" << std::endl;
	}
	
	// let the threads initialize
	pthread_cond_wait(&condvar_tomf, &mutex);
	pthread_mutex_unlock(&mutex);
	
	for (s = 0; s < matrix->numberOfBlockColumns; ++s) {
		processedthreads = 0;
		pthread_mutex_lock(&mutex);
		minnorm = -1;
		mini = s;
		minj = s;
		pthread_cond_broadcast(&condvar_totf);
		pthread_cond_wait(&condvar_tomf, &mutex);
		pthread_mutex_unlock(&mutex);
#ifdef DEBUG
		std::cout << "minblock: (" << mini << ", " << minj << ")" << std::endl;
#endif
		if (s != minj)
			matrix_swap_block_columns(matrix, minj, s);
		if (s != mini) {
			matrix_swap_block_rows(matrix, mini, s);
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
		processedthreads = 0;
		pthread_mutex_lock(&mutex);
		pthread_cond_broadcast(&condvar_totf);
		if (s + 1 != matrix->numberOfBlockColumns)
			pthread_cond_wait(&condvar_tomf, &mutex);
		pthread_mutex_unlock(&mutex);
		for (i = s+1; i < matrix->numberOfBlockColumns; ++i) {
			block = matrix_get_pos_block(matrix, i, s);
			index = blocksize * matrix_get_pos_height(matrix, i);
			for (j = 0; j < index; ++j)
				block[j] = 0;
		}
	}
	for (i = matrix->numberOfBlockRows; i >= 0; --i)
		for (j = i+1; j < matrix->numberOfBlockColumns; ++j) {
			block_apply_to_vector(
				matrix_get_pos_height(matrix, j),
				matrix_get_pos_height(matrix, i),
				matrix_get_pos_block(matrix, i, j),
				&(rightcol[j*blocksize]),
				tempmatrix
			);
			block_subtract(
				matrix_get_pos_height(matrix, i),
				&(rightcol[i*blocksize]),
				tempmatrix
			);
		}
	for (i = 0; i < threads; ++i)
		pthread_join(thr[i], NULL);
	pthread_mutex_destroy(&mutex);
	delete[] tempmatrix;
	delete[] buf;
	delete[] thr;
	delete[] args;
}
