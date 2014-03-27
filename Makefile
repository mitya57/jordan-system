CXX := c++
CXXFLAGS := -lpthread

debug:
	$(CXX) $(CXXFLAGS) -g -o matrix matrix.cc

release:
	$(CXX) $(CXXFLAGS) -O3 -o matrix matrix.cc

profile-generate:
	$(CXX) $(CXXFLAGS) -O3 -fprofile-generate -o matrix matrix.cc

profile-use:
	$(CXX) $(CXXFLAGS) -O3 -fprofile-use -o matrix matrix.cc

clean:
	rm -f matrix matrix_tests *.o* *.gc*

pg: profile-generate
pu: profile-use
main: debug
