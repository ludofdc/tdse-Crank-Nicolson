# Makefile - TDSE Crank-Nicolson solver

CXX      = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

SRC_DIR  = src
OUT_DIR  = output

SRCS     = $(SRC_DIR)/tdse.cpp $(SRC_DIR)/params.cpp $(SRC_DIR)/physics.cpp
TARGET   = tdse

$(shell mkdir -p $(OUT_DIR))

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS)

clean:
	rm -f $(TARGET)