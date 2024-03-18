CXX = g++
CXXFLAGS += -std=c++20 
CPPFLAGS = -Wall -Wextra -Iinc -Iexternals/json/include

INC_DIR = inc
SRC_DIR = src
OBJ_DIR = obj

SRCS = $(wildcard $(SRC_DIR)/*.cpp) 
HDRS = $(wildcard $(INC_DIR)/*.hpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
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
