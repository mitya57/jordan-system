/***************************************************************
 *       Solving linear systems using the Jordan method        *
 *                Author: Dmitry Shachnev, 2013                *
 * This file contains the main function of MPI implementation. *
 ***************************************************************/

#define PRETTY(x) (((x) < EPS && (x) > -EPS) ? 0 : x)

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "mpi.h"

#include "matrix.hh"
#include "block.hh"
#include "function.hh"
#include "mpi-util.hh"

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);

	MatrixInfo matrix_info;
	MPI_Data data;
	char *filename;
	int size, blocksize, rank, s;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	data.matrix_info = &matrix_info;

	if (argc >= 3) {
		size = atoi(argv[1]);
		blocksize = atoi(argv[2]);
	} else {
		if (!rank)
			std::cerr << "Usage: matrix SIZE BLOCKSIZE [FILENAME]"
			<< std::endl;
		MPI_Finalize();
		return 0;
	}

	matrix_new(&matrix_info, size, blocksize);
	mpi_initialize_data(&data);
	mpi_check_supported(&data);
	filename = (argc >= 4 ? argv[3] : 0);
	mpi_read_matrix(&data, filename);
	mpi_print_matrix(&data);

	double starttime = MPI_Wtime();
	for (s = 0; s < matrix_info.numberOfBlockColumns; ++s) {
		mpi_find_and_move_main_block(&data, s);
		mpi_process_main_block_row_and_subtract_rows(&data, s);
	}
	for (s = matrix_info.numberOfBlockColumns - 1; s >= 0; --s)
		mpi_reverse_subtract(&data, s);
	//mpi_print_matrix(&data);

	double timediff = MPI_Wtime() - starttime;
	double residual;

	mpi_get_rightcol(&data, data.mrightcol);
	matrix_clear(data.matrix_info);
	mpi_read_matrix(&data, filename);
	mpi_get_residual(&data, data.mrightcol, &residual);

	if (!rank) {
		std::cout << "RESULT:";
		for (s = 0; s < ((size < 10) ? size : 10); ++s)
			std::cout << " " << PRETTY(data.mrightcol[s]);
		if (size > 10)
			std::cout << " ...";
		std::cout << std::endl;
		std::cout << "TIME: " << timediff << " seconds" << std::endl;
		std::cout << "RESIDUAL: " << sqrt(residual) << std::endl << std::endl;
	}

	mpi_deinitialize_data(&data);
	matrix_free(&matrix_info);
	MPI_Finalize();
}
