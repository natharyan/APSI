CXX = g++
GMP_PREFIX = /opt/homebrew/opt/gmp
OPENSSL_PREFIX = /opt/homebrew/opt/openssl@3
MCL_PREFIX = /usr/local

CXXFLAGS = -std=c++17 -O2 -Wall -Wno-shift-count-overflow -I$(GMP_PREFIX)/include -I$(OPENSSL_PREFIX)/include -I./src -I./src/hashing -I./src/utils
LDFLAGS = -L$(GMP_PREFIX)/lib -L$(OPENSSL_PREFIX)/lib -L$(MCL_PREFIX)/lib -lmcl -lgmp -lcrypto -Wl,-rpath,$(MCL_PREFIX)/lib

SRC_UTIL = src/utils/util.cpp
SRC_HASH = src/hashing/simple_hashing.cpp

all: client server

client: src/psi/client.cpp $(SRC_UTIL) $(SRC_HASH)
	$(CXX) $(CXXFLAGS) -o client src/psi/client.cpp $(SRC_UTIL) $(SRC_HASH) $(LDFLAGS)

server: src/psi/server.cpp $(SRC_UTIL) $(SRC_HASH)
	$(CXX) $(CXXFLAGS) -o server src/psi/server.cpp $(SRC_UTIL) $(SRC_HASH) $(LDFLAGS)

clean:
	rm -f client server

.PHONY: clean