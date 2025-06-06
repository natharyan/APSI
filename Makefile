CXX = g++
CXXFLAGS = -std=c++17 -O2 -I/opt/homebrew/opt/gmp/include -I./src -I./src/hashing -I./src/utils
LDFLAGS = -L/opt/homebrew/opt/gmp/lib -lgmp

SRC_UTIL = src/utils/util.cpp
SRC_HASH = src/hashing/simple_hashing.cpp

all: alice bob

alice: src/psi/alice.cpp $(SRC_UTIL) $(SRC_HASH)
	$(CXX) $(CXXFLAGS) -o alice src/psi/alice.cpp $(SRC_UTIL) $(SRC_HASH) $(LDFLAGS)

bob: src/psi/bob.cpp $(SRC_UTIL) $(SRC_HASH)
	$(CXX) $(CXXFLAGS) -o bob src/psi/bob.cpp $(SRC_UTIL) $(SRC_HASH) $(LDFLAGS)

clean:
	rm -f alice bob