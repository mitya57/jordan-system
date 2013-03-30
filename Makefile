MPICXX := mpic++
CXXFLAGS := -Wall -o mpi-main

mpirelease: mpi-main.cc mpi-util.hh matrix.hh block.hh
	$(MPICXX) $(CXXFLAGS) -O3 mpi-main.cc && \
	touch mpirelease

mpidebug: mpi-main.cc mpi-util.hh matrix.hh block.hh
	$(MPICXX) $(CXXFLAGS) -g -DDEBUG mpi-main.cc && \
	touch mpidebug

mpidebug-alt: mpi-main.cc mpi-util.hh matrix.hh block.hh
	$(MPICXX) $(CXXFLAGS) -g mpi-main.cc && \
	touch mpidebug-alt

clean:
	rm -f *.o* *.gc* mpi-main mpidebug* mpirelease
