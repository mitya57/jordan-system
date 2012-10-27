/******************************************************
 *   Solving Linear systems using the Jordan method   *
 *         Copyright: Dmitry Shachnev, 2012           *
 * This file contains the matrix generation function. *
 ******************************************************/

double get_matrix_element(int i, int j) {
	return 1. / (1 + i + j);
}
