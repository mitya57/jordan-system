/******************************************************
 *   Solving Linear systems using the Jordan method   *
 *         Copyright: Dmitry Shachnev, 2012           *
 * This file contains the matrix generation function. *
 ******************************************************/

double get_matrix_element(int i, int j) {
	return fabs(i-j);
}
