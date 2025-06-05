CXX = g++
GMP_PREFIX = /opt/homebrew/opt/gmp
CXXFLAGS = -std=c++17 -O2 -I$(GMP_PREFIX)/include
LDFLAGS = -L$(GMP_PREFIX)/lib -lgmp
SOURCES = src/psi/main.cpp src/utils/util.cpp
TARGET = main

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES) $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: clean