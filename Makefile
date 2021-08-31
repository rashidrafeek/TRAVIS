#*****************************************************************************
#    TRAVIS - Trajectory Analyzer and Visualizer
#    http://www.travis-analyzer.de/
#
#    Copyright (c) 2009-2021 Martin Brehm
#                  2012-2021 Martin Thomas
#                  2016-2021 Sascha Gehrke
#
#    This file written by Martin Brehm.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************

# TRAVIS Makefile

# Edit the follwing lines according to your needs

# Your C++ compiler
CXX         = g++

#**************************************************************************
# You should not need to edit the following part
#**************************************************************************

CFLAGS      = -g -Wall -Wextra -Wformat -Wformat-security -pedantic -std=c++14 -fsigned-char -pthread -O2
CFLAGS_PROF = -g -Wall -Wextra -Wformat -Wformat-security -pedantic -std=c++14 -fsigned-char -pthread -O2 -pg
CFLAGS_DEB  = -g -Wall -Wextra -Wformat -Wformat-security -pedantic -std=c++14 -fsigned-char -pthread -O0
CFLAGS_REL  = -g -Wall -Wextra -Wformat -Wformat-security -pedantic -std=c++14 -fsigned-char -pthread -O3 -march=athlon-fx
CFLAGS_OMP  = -g -Wall -Wextra -Wformat -Wformat-security -pedantic -std=c++14 -fsigned-char -pthread -O2 -fopenmp

SRC  = $(wildcard src/*.cpp)
SRC += $(wildcard src/libtravis/*.cpp)
SRC += $(wildcard src/libtravis/kiss_fft/*.cpp)
SRC += $(wildcard src/libtravis/voro++/*.cpp)
OBJ  = $(SRC:%.cpp=%.o)
BIN  = exe/travis

all: executable

debug: CFLAGS = $(CFLAGS_DEB)
debug: executable

release: CFLAGS = $(CFLAGS_REL)
release: executable

profile: CFLAGS = $(CFLAGS_PROF)
profile: executable

omp: CFLAGS = $(CFLAGS_OMP)
omp: executable

executable: $(OBJ)
	mkdir -p exe
	$(CXX) $(CFLAGS) -o $(BIN) $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJ)

.PHONY: distclean
distclean:
	rm -f $(OBJ)
	rm -rf exe


