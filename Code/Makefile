CXX = g++
CXXFLAGS = -Wextra -O3 -fopenmp -march=native -ftree-loop-im -ftree-loop-ivcanon -fvariable-expansion-in-unroller -g -MMD -MP

SRC = sbox.cpp poly.cpp precomputation.cpp decomp.cpp reconstitution.cpp test_main.cpp
OBJ = $(SRC:.cpp=.o)

PROGRAM = implem

SRC_QUAD = my_slp_heuristic_for_quad.cpp
SRC_LIN = my_slp_heuristic_for_lin.cpp
SLP_QUAD = my_slp_heuristic_quad
SLP_LIN = my_slp_heuristic_lin

all: $(PROGRAM) slp

$(PROGRAM): $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

slp:
	$(CXX) $(CXXFLAGS) $(SRC_QUAD) -o $(SLP_QUAD)
	$(CXX) $(CXXFLAGS) $(SRC_LIN) -o $(SLP_LIN)

clean:
	rm -f $(PROGRAM) $(SLP_QUAD) $(SLP_LIN) *.o *.d

launch: clean all
	./$(PROGRAM)

-include $(OBJ:.o=.d)
