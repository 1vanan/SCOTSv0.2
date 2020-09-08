#
# compiler
#
#CC        = g++
CC       	  = clang++
CXXFLAGS 		= -g -Wall -Wextra -std=c++11
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG

# TODO: check this file
#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src -I$(SCOTSROOT)/utils


.PHONY: power_system simulate

TARGET = power_system simulate

all: $(TARGET)

power_system:
	$(CC) $(CXXFLAGS) $(SCOTSINC) power_system.cc -o power_system

simulate: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) simulate.cc -o simulate


clean:
	rm power_system simulate
