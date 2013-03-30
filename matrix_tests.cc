/**************************************************
 * Solving Linear systems using the Jordan method *
 *        Copyright: Dmitry Shachnev, 2012        *
 *         This file contains some tests.         *
 **************************************************/

#define EPS 1e-8
#define EQUAL(m, n) fabs(m-(n)) < EPS
#define NO_INCLUDE

#include <algorithm>
#include <cassert>
#include <cmath>
#include "block.hh"
#include "matrix.hh"

void test_block_functions() {
	double *block1, *block2;
	double *buf;
	int i;
	// test subtraction
	block1 = new double[9];
	block2 = new double[9];
	buf = new double[9];
	for (i=0; i<9; ++i) {
		block1[i] = 3*i;
		block2[i] = i;
	}
	block_subtract(9, block1, block2);
	for (i=0; i<9; ++i)
		assert(EQUAL(block1[i], 2*i));
	// test norm calculating
	assert(EQUAL(block_get_norm(3, block1), 372));
	// test multiplication
	block_left_multiply_alt(3, 2, block2, block1, buf);
	assert(EQUAL(block2[0], 20));
	assert(EQUAL(block2[1], 26));
	assert(EQUAL(block2[2], 56));
	assert(EQUAL(block2[3], 80));
	assert(EQUAL(block2[4], 92));
	assert(EQUAL(block2[5], 134));
	// test reversing
	block1[0] = 1;
	block1[1] = 1;
	block1[2] = 2;
	block1[3] = 3;
	block_get_reverse(2, block1, block2, buf);
	assert(EQUAL(block2[0], 3));
	assert(EQUAL(block2[1], -1));
	assert(EQUAL(block2[2], -2));
	assert(EQUAL(block2[3], 1));
	// another reversing test
	block1[0] = 2;
	block1[1] = -1;
	block1[2] = 0;
	block1[3] = -1;
	block1[4] = 2;
	block1[5] = -1;
	block1[6] = 0;
	block1[7] = -1;
	block1[8] = 2;
	block_get_reverse(3, block1, block2, buf);
	assert(EQUAL(block2[0], 0.75));
	assert(EQUAL(block2[1], 0.5));
	assert(EQUAL(block2[2], 0.25));
	assert(EQUAL(block2[3], 0.5));
	assert(EQUAL(block2[4], 1));
	assert(EQUAL(block2[5], 0.5));
	assert(EQUAL(block2[6], 0.25));
	assert(EQUAL(block2[7], 0.5));
	assert(EQUAL(block2[8], 0.75));
	delete[] buf;
	delete[] block1;
	delete[] block2;
}

int main(void) {
	test_block_functions();
	return 0;
}
