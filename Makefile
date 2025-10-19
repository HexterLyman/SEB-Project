CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -O2
DEBUGFLAGS = -g -DDEBUG

SRC_DIR = src
UTILS_DIR = src/Utils
HEADER_DIR = src/Header
DATA_DIR = src/Data

SRCS = $(SRC_DIR)/prototype_v1.cpp $(UTILS_DIR)/geometric_functions.cpp
OBJS = $(SRCS:.cpp=.o)
EXEC = src/program

.PHONY: all debug clean

all: $(EXEC)

debug: CXXFLAGS += $(DEBUGFLAGS)
debug: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(HEADER_DIR) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)
