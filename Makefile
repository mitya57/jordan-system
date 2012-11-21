CXX := g++
CXXFLAGS := -lpthread

debug:
	$(CXX) $(CXXFLAGS) -g -DDEBUG -o matrix matrix.cc

release:
	$(CXX) $(CXXFLAGS) -O3 -o matrix matrix.cc

profile-generate:
	$(CXX) $(CXXFLAGS) -O3 -fprofile-generate -o matrix matrix.cc

profile-use:
	$(CXX) $(CXXFLAGS) -O3 -fprofile-use -o matrix matrix.cc

clean:
	rm matrix matrix_tests *.o* *.gc*

pg: profile-generate
pu: profile-use
main: debug
