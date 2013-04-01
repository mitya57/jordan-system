/**************************************************
 * Solving linear systems using the Jordan method *
 *         Author: Dmitry Shachnev, 2013          *
 *    This file contains MPI utility functions.   *
 **************************************************/

#define M_INFO data->matrix_info
#define BLOCKSIZE (M_INFO->blocksize)
#define REAL_ROW_HEIGHT (data->bpp * BLOCKSIZE)
#define ROW_HEIGHT ((data->processes > 1) ? REAL_ROW_HEIGHT : SIZE)
#define ROW_HEIGHT_LAST (((SIZE - 1) % REAL_ROW_HEIGHT) + 1)
#define ROW_HEIGHT_FOR_P(p) ((p < data->processes - 1) ? ROW_HEIGHT : ROW_HEIGHT_LAST)
#define ROW_SIZE (ROW_HEIGHT * SIZE)
#define REGULAR_PROCESS (data->rank + 1 < data->processes)
#define BLOCKS_PER_PROCESS_LAST ((ROW_HEIGHT_LAST - 1) / BLOCKSIZE + 1)
#define BLOCKS_PER_PROCESS_MY (REGULAR_PROCESS ? data->bpp : BLOCKS_PER_PROCESS_LAST)
#define SIZE (M_INFO->size)
#define PRINT_SIZE ((SIZE < 10) ? SIZE : 10)
#define PRINT_FORMAT(x) ((0 < x && x < 10) ? (double(int(x*1e4))/1e4) : (double(int(x*1e3))/1e3))

#define mpi_get_pos_width mpi_get_pos_height

/****************** NOTE *******************
 *  blocks  are denoted by x and y letters *
 * elements are denoted by i and j letters *
 *******************************************/

struct MPI_Data {
	int rank;
	int processes;
	int bpp; // blockrows per process
	double *matrix;
	double *buffer;
	double *tempbuf;
	double *blockbuf;
	double *rightcol;
	double *rightcolbuf;
	double *mrightcol;
	MatrixInfo *matrix_info;
};

struct MPI_Double_Int {
	double value;
	int pos;
};

void mpi_initialize_data(MPI_Data *data) {
	MPI_Comm_rank(MPI_COMM_WORLD, &(data->rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(data->processes));
	data->bpp = (M_INFO->numberOfBlockRows - 1) / data->processes + 1;
	data->matrix = new double[ROW_SIZE];
	data->buffer = new double[ROW_SIZE];
	data->tempbuf = new double[BLOCKSIZE * BLOCKSIZE];
	data->blockbuf = new double[BLOCKSIZE * BLOCKSIZE];
	data->rightcol = new double[ROW_HEIGHT];
	data->rightcolbuf = new double[ROW_HEIGHT];
	data->mrightcol = new double[SIZE];
}

void mpi_deinitialize_data(MPI_Data *data) {
	delete[] data->matrix;
	delete[] data->buffer;
	delete[] data->tempbuf;
	delete[] data->blockbuf;
	delete[] data->rightcol;
	delete[] data->rightcolbuf;
	delete[] data->mrightcol;
}

void mpi_abort_program(MPI_Data *data) {
	mpi_deinitialize_data(data);
	MPI_Abort(MPI_COMM_WORLD, 1);
}

void mpi_check_supported(MPI_Data *data) {
	if (ROW_HEIGHT * (data->processes - 1) + ROW_HEIGHT_LAST != SIZE) {
		if (!(data->rank))
			std::cerr << "ERROR: unsupported configuration; try"
			<< " running with less processes." << std::endl;
		mpi_abort_program(data);
	}
}

int mpi_get_pos_height(MPI_Data *data, int x) {
	// x is *global* here!
	// use mpi_localx_to_matrixx to convert.
	if (x < SIZE / M_INFO->blocksize)
		return M_INFO->blocksize;
	return SIZE % M_INFO->blocksize;
}

int mpi_get_pos_size(MPI_Data *data, int x, int y) {
	// x and y are *global* here!
	return mpi_get_pos_height(data, x) * mpi_get_pos_width(data, y);
}

double *mpi_get_pos_block(MPI_Data *data, int x, int y, bool _smode = false) {
	// x and y are *global* here!
	// use mpi_localx_to_matrixx and mpi_localy_to_matrixy to convert.
	int localx = M_INFO->blockrowind[x * BLOCKSIZE] % data->bpp;
	int localy = (_smode ? 0 : M_INFO->blockcolumnind[y * BLOCKSIZE]);
	int start = BLOCKSIZE * SIZE * localx;
	start += mpi_get_pos_height(data, x) * BLOCKSIZE * localy;
	return &((data->matrix)[start]);
}

double *mpi_get_pos_rightcol_block(MPI_Data *data, int x) {
	// x is *global* here!
	// use mpi_localx_to_matrixx to convert.
	int localx = M_INFO->blockrowind[x * BLOCKSIZE] % data->bpp;
	return &((data->rightcol)[localx * BLOCKSIZE]);
}

double *mpi_get_local_block(MPI_Data *data, int localx, int localy,
double *matrix = 0, int p = -1) {
	int blockheight;
	if (p == -1)
		p = data->rank;
	if (p + 1 < data->processes || (localx + 1) * BLOCKSIZE < ROW_HEIGHT_LAST)
		blockheight = BLOCKSIZE;
	else
		blockheight = ROW_HEIGHT_LAST - localx * BLOCKSIZE;
	int start = BLOCKSIZE * SIZE * localx;
	start += blockheight * BLOCKSIZE * localy;
	if (!matrix) matrix = data->matrix;
	return &(matrix[start]);
}

int mpi_get_local_width(MPI_Data *data, int localy) {
	if ((localy + 1) * BLOCKSIZE <= SIZE)
		return BLOCKSIZE;
	return SIZE % BLOCKSIZE;
}

int mpi_localx_to_matrixx(MPI_Data *data, int localx) {
	int pos = localx + data->rank * data->bpp;
	int i = 0;
	while (M_INFO->blockrowind[i] != pos)
		++i;
	return i / BLOCKSIZE;
}

int mpi_localy_to_matrixy(MPI_Data *data, int localy) {
	int j = 0;
	while (M_INFO->blockcolumnind[j] != localy)
		++j;
	return j / BLOCKSIZE;
}

double *mpi_local_element(MPI_Data *data, int locali, int localj,
double *matrix = 0, int p = -1) {
	double *block = mpi_get_local_block(data,
		locali / BLOCKSIZE, localj / BLOCKSIZE, matrix, p);
	int width = mpi_get_local_width(data, localj / BLOCKSIZE);
	return &(block[(locali % BLOCKSIZE) * width + (localj % BLOCKSIZE)]);
}

double *mpi_global_element(MPI_Data *data, int i, int j) {
	int localy = M_INFO->blockcolumnind[j];
	int localj = localy * BLOCKSIZE + (j % BLOCKSIZE);
	int localx = M_INFO->blockrowind[i] % data->bpp;
	int locali = localx * BLOCKSIZE + (i % BLOCKSIZE);
	return mpi_local_element(data, locali, localj);
}

int mpi_rank_for_i(MPI_Data *data, int i) {
	int x = M_INFO->blockrowind[i];
	return x / data->bpp;
}

void mpi_print_matrix(MPI_Data *data) {
	int i, j, p, vblock, locali;
	double el;
	MPI_Status status;
	//MPI_Barrier(MPI_COMM_WORLD);
	if (!data->rank)
		std::cout << "==========BEGIN=PRINT=MATRIX==========" << std::endl;
	for (i = 0; i < PRINT_SIZE; ++i) {
		p = mpi_rank_for_i(data, i);
		if (p) {
			if (!(data->rank)) {
				// receive and print
				for (j = 0; j < PRINT_SIZE; ++j) {
					MPI_Recv(&el, 1, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &status);
					std::cout << el << "\t";
				}
				if (SIZE > 10)
					std::cout << "...";
				MPI_Recv(&el, 1, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &status);
				std::cout << "\t" << el << std::endl;
			} else if (data->rank == p) {
				// send
				for (j = 0; j < PRINT_SIZE; ++j) {
					el = PRINT_FORMAT(*mpi_global_element(data, i, j));
					MPI_Send(&el, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
				}
				vblock = M_INFO->blockrowind[i];
				locali = (vblock % data->bpp) * BLOCKSIZE + (i % BLOCKSIZE);
				el = data->rightcol[locali % ROW_HEIGHT];
				MPI_Send(&el, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			}
		} else if (!(data->rank)) {
			// just print
			for (j = 0; j < PRINT_SIZE; ++j) {
				el = PRINT_FORMAT(*mpi_global_element(data, i, j));
				std::cout << el << "\t";
			}
			if (SIZE > 10)
				std::cout << "...";
			vblock = M_INFO->blockrowind[i];
			locali = (vblock % data->bpp) * BLOCKSIZE + (i % BLOCKSIZE);
			el = data->rightcol[locali % ROW_HEIGHT];
			std::cout << "\t" << el << std::endl;
		}
	}
	if (!data->rank) {
		if (SIZE > 10)
			std::cout << ". . . . . . . . . . . . . . . . . . ." << std::endl;
		std::cout << "===========END=PRINT=MATRIX===========" << std::endl;
	}
	//MPI_Barrier(MPI_COMM_WORLD);
}

void mpi_read_matrix(MPI_Data *data, char *filename) {
	int i, j, p;
	MPI_Status status;
	if (data->rank) {
		MPI_Recv(data->matrix, ROW_SIZE, MPI_DOUBLE, 0,
			data->rank, MPI_COMM_WORLD, &status);
		MPI_Recv(data->rightcol, ROW_HEIGHT, MPI_DOUBLE, 0,
			data->rank, MPI_COMM_WORLD, &status);
	} else {
		std::fstream matrixfile;
		if (filename) {
			matrixfile.open(filename);
			if (!matrixfile && !(data->rank))
				std::cerr << "ERROR: cannot open file " << filename
				<< std::endl;
			if (!matrixfile)
				mpi_abort_program(data);
		}
		// read our piece first
		for (i = 0; i < ROW_HEIGHT; ++i) {
			for (j = 0; j < SIZE; ++j) {
				if (filename)
					matrixfile >> *mpi_local_element(data, i, j);
				else *mpi_local_element(data, i, j) =
					get_matrix_element(i, j);
			}
			if (filename) matrixfile >> data->rightcol[i];
			else {
				data->rightcol[i] = 0;
				for (j = 0; j < SIZE; j += 2)
					data->rightcol[i] += get_matrix_element(i, j);
			}
		}
		// read the rest now
		for (p = 1; p < data->processes; ++p) {
			for (i = 0; i < ROW_HEIGHT_FOR_P(p); ++i) {
				for (j = 0; j < SIZE; ++j) {
					if (filename) matrixfile >>
						*mpi_local_element(data, i, j, data->buffer, p);
					else
						*mpi_local_element(data, i, j, data->buffer, p) =
						get_matrix_element(i + ROW_HEIGHT * p, j);
				}
				if (filename)
					matrixfile >> data->rightcolbuf[i];
				else {
					data->rightcolbuf[i] = 0;
					for (j = 0; j < SIZE; j += 2)
						data->rightcolbuf[i] +=
						get_matrix_element(i + ROW_HEIGHT * p, j);
				}
			}
			MPI_Send(data->buffer, ROW_SIZE, MPI_DOUBLE, p, p, MPI_COMM_WORLD);
			MPI_Send(data->rightcolbuf, ROW_HEIGHT, MPI_DOUBLE, p, p,
				MPI_COMM_WORLD);
		}
		if (filename)
			matrixfile.close();
	}
}

void mpi_find_and_move_main_block(MPI_Data *data, int s) {
	MPI_Double_Int sendbuf, recvbuf;
	int matrixx, matrixy;
	double norm;
	sendbuf.value = -1;
	sendbuf.pos = -1;
	// find the main block
	for (matrixx = s; (matrixx+1)*BLOCKSIZE <= ROW_HEIGHT_FOR_P(data->rank); ++matrixx)
		for (matrixy = s; (matrixy+1)*BLOCKSIZE <= SIZE; ++matrixy) {
			// calculate norm
			if (!block_get_reverse(BLOCKSIZE,
			mpi_get_pos_block(data, matrixx, matrixy),
			data->blockbuf, data->tempbuf))
				continue;
			norm = block_get_norm(BLOCKSIZE, data->blockbuf);
			// if norm is the smallest then calculate global block index
			if (norm < sendbuf.value || sendbuf.value < 0) {
				sendbuf.value = norm;
				sendbuf.pos = matrixx * M_INFO->numberOfBlockColumns + matrixy;
			}
		}
	MPI_Reduce(&sendbuf, &recvbuf, 1, MPI_DOUBLE_INT,
		MPI_MINLOC, 0, MPI_COMM_WORLD);
	// move it, move it
	if (!(data->rank) && recvbuf.value > 0) {
		int mini = recvbuf.pos / M_INFO->numberOfBlockColumns;
		int minj = recvbuf.pos % M_INFO->numberOfBlockColumns;
		matrix_swap_block_rows(M_INFO, s, mini);
		matrix_swap_block_columns(M_INFO, s, minj);
	}
	// let others know about the change
	MPI_Bcast(M_INFO->blockcolumnind, SIZE, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(M_INFO->blockrowind, SIZE, MPI_INT, 0, MPI_COMM_WORLD);
}

void mpi_process_main_block_row_and_subtract_rows(MPI_Data *data, int s) {
	int p = mpi_rank_for_i(data, s * BLOCKSIZE);
	int x, y, localy;
	double *sspointer, *bpointer, *spointer;
	// process main block row
	if (data->rank == p) {
		sspointer = mpi_get_pos_block(data, s, s);
		bpointer = mpi_get_pos_rightcol_block(data, s);
		block_get_reverse(
			mpi_get_pos_height(data, s),
			sspointer,
			data->tempbuf,
			data->blockbuf
		);
		for (y = s+1; y < M_INFO->numberOfBlockColumns; ++y)
			block_left_multiply_alt(
				BLOCKSIZE,
				mpi_get_pos_width(data, y),
				mpi_get_pos_block(data, s, y),
				data->tempbuf,
				data->blockbuf
			);
		block_apply_to_vector_alt(
			mpi_get_pos_height(data, s),
			data->tempbuf,
			bpointer,
			data->blockbuf
		);
		spointer = mpi_get_pos_block(data, s, 0, true);
		//width = mpi_get_pos_width(data, s);
		//for (i = 0; i < mpi_get_pos_height(data, s); ++i)
		//	for (j = 0; j < width; ++j)
		//		sspointer[i*width+j] = (i == j);
	} else {
		spointer = data->buffer;
		bpointer = data->rightcolbuf;
	}
	// broadcast it
	MPI_Bcast(spointer, BLOCKSIZE * SIZE, MPI_DOUBLE,
		p, MPI_COMM_WORLD);
	MPI_Bcast(bpointer, BLOCKSIZE, MPI_DOUBLE, p, MPI_COMM_WORLD);
	// subtract block rows
	for (x = s + 1; x < M_INFO->numberOfBlockRows; ++x) {
		p = mpi_rank_for_i(data, x * BLOCKSIZE);
		if (data->rank == p) {
			for (y = s+1; y < M_INFO->numberOfBlockColumns; ++y) {
				localy = M_INFO->blockcolumnind[y * BLOCKSIZE];
				block_left_multiply(
					mpi_get_pos_height(data, x),
					mpi_get_pos_width(data, s),
					mpi_get_pos_height(data, y),
					&(spointer[localy * BLOCKSIZE * BLOCKSIZE]),
					//mpi_get_pos_block(data, s, y)
					mpi_get_pos_block(data, x, s),
					data->tempbuf
				);
				block_subtract(
					mpi_get_pos_size(data, x, y),
					mpi_get_pos_block(data, x, y),
					data->tempbuf
				);
			}
			// When in DEBUG mode, empty the (x, s) block?
			block_apply_to_vector(
				BLOCKSIZE,
				mpi_get_pos_width(data, x),
				mpi_get_pos_block(data, x, s),
				bpointer, //mpi_get_pos_rightcol_block(data, s)
				data->tempbuf
			);
			block_subtract(
				mpi_get_pos_height(data, x),
				mpi_get_pos_rightcol_block(data, x),
				data->tempbuf
			);
		}
	}
}

void mpi_reverse_subtract(MPI_Data *data, int s) {
	// broadcast the block row to subtract
	int p = mpi_rank_for_i(data, s * BLOCKSIZE);
	double *bpointer;
	if (data->rank == p)
		bpointer = mpi_get_pos_rightcol_block(data, s);
	else
		bpointer = data->rightcolbuf;
	MPI_Bcast(bpointer, BLOCKSIZE, MPI_DOUBLE, p, MPI_COMM_WORLD);
	// now subtract
	for (int x = s - 1; x >= 0; --x) {
		p = mpi_rank_for_i(data, x * BLOCKSIZE);
		if (data->rank == p) {
			block_apply_to_vector(
				mpi_get_pos_height(data, s),
				mpi_get_pos_height(data, x),
				mpi_get_pos_block(data, x, s),
				bpointer,
				data->tempbuf
			);
			block_subtract(
				mpi_get_pos_height(data, x),
				mpi_get_pos_rightcol_block(data, x),
				data->tempbuf
			);
		}
	}
}

void mpi_get_residual(MPI_Data *data, double *vector, double *residual) {
	int i, x, y;
	double localresult = 0, diff;
	MPI_Bcast(vector, SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	int start = data->rank * data->bpp;
	for (x = start; x < start + BLOCKS_PER_PROCESS_MY; ++x) {
		for (i = 0; i < mpi_get_pos_height(data, x); ++i)
			data->rightcolbuf[i] = 0;
		for (y = 0; y < M_INFO->numberOfBlockColumns; ++y) {
			block_apply_to_vector(
				mpi_get_pos_width(data, y),
				mpi_get_pos_height(data, x),
				mpi_get_pos_block(data, x, y),
				&(vector[y * BLOCKSIZE]),
				data->tempbuf
			);
			for (i = 0; i < mpi_get_pos_height(data, x); ++i)
				data->rightcolbuf[i] += data->tempbuf[i];
		}
		for (i = 0; i < mpi_get_pos_height(data, x); ++i) {
			diff = data->rightcol[(x - start) * BLOCKSIZE + i] - data->rightcolbuf[i];
			localresult += diff * diff;
		}
	}
	MPI_Reduce(&localresult, residual, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

void mpi_get_rightcol(MPI_Data *data, double *realrightcol) {
	int i, j, k;
	MPI_Status status;
	for (i = 0; i < SIZE; ++i) {
		j = 0;
		while (M_INFO->blockcolumnind[j] != i/BLOCKSIZE) ++j;
		k = mpi_rank_for_i(data, j);
		j = M_INFO->blockrowind[j] * BLOCKSIZE;
		if (!k)
			realrightcol[i] = data->rightcol[j + i % BLOCKSIZE];
		else if (data->rank == k)
			MPI_Send(&(data->rightcol[(j + i % BLOCKSIZE) % ROW_HEIGHT]),
			1, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
		else if (!(data->rank))
			MPI_Recv(&realrightcol[i], 1, MPI_DOUBLE, k, i, MPI_COMM_WORLD, &status);
	}
}
