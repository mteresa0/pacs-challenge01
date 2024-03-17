CXX = g++
CXXFLAGS += -std=c++20 
CPPFLAGS = -Wall -Wextra -Iinc -Iexternals/json/include -Isrc

INC_DIR = inc, submodules/json/include
SRC_DIR = src
OBJ_DIR = obj

SRCS_IN_DIR = $(wildcard $(SRC_DIR)/*.cpp) 
SRCS = $(wildcard *.cpp) $(SRCS_IN_DIR)
HDRS = $(wildcard $(INC_DIR)/*.hpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS_IN_DIR)) $(OBJ_DIR)/main.o
EXEC = main

all : $(EXEC) 

run : clean
	make $(EXEC)
	./main

$(EXEC) : $(SRCS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^

.PHONY : clean

clean : 
	rm -f main
