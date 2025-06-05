CXX = g++
GMP_PREFIX = /opt/homebrew/opt/gmp
OPENSSL_PREFIX = /opt/homebrew/opt/openssl@3
MCL_PREFIX = /usr/local


CXXFLAGS = -std=c++17 -O2 -Wall -I$(GMP_PREFIX)/include -I$(OPENSSL_PREFIX)/include
LDFLAGS = -L$(GMP_PREFIX)/lib -L$(OPENSSL_PREFIX)/lib -L$(MCL_PREFIX)/lib -lmcl -lgmp -lcrypto -Wl,-rpath,$(MCL_PREFIX)/lib, /usr/local/lib/libmcl.dylib

SOURCES = main.cpp kzg.cpp polynomial.cpp ../utils/util.cpp
TARGET = kzg_test

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES) $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: clean