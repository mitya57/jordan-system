/***************************************************
 *  Solving Linear systems using the Jordan method *
 *        Copyright: Dmitry Shachnev, 2012         *
 * This file contains the block-related functions. *
 ***************************************************/

#ifndef NO_INCLUDE
#include <cmath>
#endif

#ifndef EPS
#define EPS 1e-8
#endif

void block_clone(int size, double *block1, double *block2) {
	/* block2 = block1 */
	for (int i = 0; i < size; ++i)
		block2[i] = block1[i];
}

bool block_get_reverse(int size, double *origblock, double *reverse, double *block) {
	int i, j, s, maxi;
	double maxel, tempd;
	block_clone(size*size, origblock, block);
	for (i = 0; i < size; ++i)
		for (j = 0; j < size; ++j) reverse[i*size+j] = (i == j);
	for (s = 0; s < size; ++s) {
		maxi = s;
		maxel = 0;
		for (i = s; i < size; ++i)
			if (block[i*size+s] > maxel) {
				maxel = block[i*size+s];
				maxi = i;
			}
		for (j = 0; j < size; ++j) {
			std::swap(block[maxi*size+j], block[s*size+j]);
			std::swap(reverse[maxi*size+j], reverse[s*size+j]);
		}
		//for (i = s; i < size; ++i)
		//	std::swap(block[i*size+maxj], block[i*size+s]);
		if (block[s*size+s] > -EPS && block[s*size+s] < EPS)
		    return false;
		tempd = 1 / block[s*size+s];
		for (j = s+1; j < size; ++j)
			block[s*size+j] *= tempd;
		for (j = 0; j < size; ++j)
			reverse[s*size+j] *= tempd;
		//block[s*size+s] = 1;
		for (i = s+1; i < size; ++i) {
			for (j = s+1; j < size; ++j)
				block[i*size+j] -= block[i*size+s]*block[s*size+j];
			for (j = 0; j < size; ++j)
				reverse[i*size+j] -= block[i*size+s]*reverse[s*size+j];
			//block[i*size+s] = 0;
		}
	}
	for (i = size - 1; i >= 0; --i) {
		for (j = i+1; j < size; ++j) {
			for (s = 0; s < size; ++s)
				reverse[i*size+s] -= reverse[j*size+s]*block[i*size+j];
			//block[i*size+j] = 0;
		}
	}
	return true;
}

double block_get_norm(int size, double *block) {
	double maxblocknorm = 0, blocknorm;
	int i, j;
	for (i = 0; i < size; ++i) {
		blocknorm = 0;
		for (j = 0; j < size; ++j) blocknorm += block[j*size+i]*block[j*size+i];
		//blocknorm = sqrt(blocknorm);
		if (blocknorm > maxblocknorm) maxblocknorm = blocknorm;
	}
	return maxblocknorm;
}

void block_left_multiply(int n, int m, int k, double *block1, double *block2, double *result) {
	/* result = block2 * block1
	 * block2 : n rows, m cols
	 * block1 : m rows, k cols
	 * result : n rows, k cols */
	int i, j, l;
	for (i = 0; i < n; i += 2)
		for (j = 0; j < k; j += 2) {
			result[i*k+j] = 0;
			for (l = 0; l < m; ++l)
				result[i*k+j] += block2[i*m+l]*block1[l*k+j];
			if (i+1 < n) {
				result[i*k+k+j] = 0;
				for (l = 0; l < m; ++l)
					result[i*k+k+j] += block2[i*m+m+l]*block1[l*k+j];
			}
			if (j+1 < k) {
				result[i*k+j+1] = 0;
				for (l = 0; l < m; ++l)
					result[i*k+j+1] += block2[i*m+l]*block1[l*k+j+1];
			}
			if (i+1 < n && j+1 < k) {
				result[i*k+k+j+1] = 0;
				for (l = 0; l < m; ++l)
					result[i*k+k+j+1] += block2[i*m+m+l]*block1[l*k+j+1];
			}
		}
}

void block_left_multiply_alt(int n, int m, double *block1, double *block2, double *tempmatrix) {
	/* block1 = block2 * block1
	 * block2 : n rows, n cols
	 * block1 : n rows, m cols */
	block_clone(n*m, block1, tempmatrix);
	block_left_multiply(n, n, m, block1, block2, tempmatrix);
	block_clone(n*m, tempmatrix, block1);
}

void block_subtract(int size, double *block1, double *block2) {
	/* block1 -= block2 */
	for (int i = 0; i < size; ++i)
		block1[i] -= block2[i];
}

void block_apply_to_vector(int size, int newsize, double *block,
double *vector, double *newvector) {
	/* block : newsize rows * size cols */
	int i, j;
	for (i = 0; i < newsize; ++i) {
		newvector[i] = 0;
		for (j = 0; j < size; ++j)
			newvector[i] += block[i*size+j]*vector[j];
	}
}

void block_apply_to_vector_alt(int size, double *block, double *vector, double *tempvector) {
	int i, j;
	for (i = 0; i < size; ++i) {
		tempvector[i] = 0;
		for (j = 0; j < size; ++j)
			tempvector[i] += block[i*size+j]*vector[j];
	}
	block_clone(size, tempvector, vector);
}
