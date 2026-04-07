CXX = g++
CXXFLAGS = -Wextra -O3 -fopenmp -march=native -ftree-loop-im -ftree-loop-ivcanon -fvariable-expansion-in-unroller -g -MMD -MP

SRC = sbox.cpp poly.cpp precomputation.cpp decomp.cpp reconstitution.cpp test_main.cpp
OBJ = $(SRC:.cpp=.o)

PROGRAM = implem

all: $(PROGRAM)

$(PROGRAM): $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(PROGRAM) *.o *.d

launch: clean all
	./$(PROGRAM)

-include $(OBJ:.o=.d)
