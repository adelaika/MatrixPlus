CC = g++ 
SRC = s21_matrix_oop.cpp
OBJ = s21_matrix_oop.o
OBJ_TEST = s21_matrix_test.o
LIB = s21_matrix_oop.a
FLAGS= -Wall -Werror -Wextra -std=c++17
TESTFLAGS= -pthread -lgtest -lgmock
CHECK_FLAG = -lstdc++
LCOVFLAGS = --ignore-errors inconsistent --ignore-errors mismatch
OPEN= open 

OS := $(shell uname -s) 
ifeq ($(OS), Linux)
	OPEN = xdg-open
	TESTFLAGS= -pthread -lgtest -lgmock -lpthread
	FLAGS= -Wall -Werror -Wextra -std=c++17 -fPIC -fprofile-arcs -ftest-coverage -lstdc++fs
endif

all: clean gcov_report 
s21_matrix_oop.a:
	$(CC) -c $(FLAGS) $(SRC) -o $(OBJ)
	ar -rcs $(LIB) $(OBJ)
	ranlib $(LIB)
	rm -rf *.o
test: clean s21_matrix_oop.a
	$(CC) $(FLAGS) -c test.cpp -o test.o
	$(CC) $(FLAGS) $(CHECK_FLAG) test.o $(LIB) $(TESTFLAGS) -o test
	rm -rf *.o 
	./test
gcov_report: add_coverage_flag test
	./test
	lcov -t "test" -o test.info -c -d ./  --no-external $(LCOVFLAGS)
	lcov --remove test.info 'test.cpp' -o test.info
	genhtml -o report test.info $(LCOVFLAGS)
	$(OPEN) report/index.html
add_coverage_flag:
	$(eval FLAGS += --coverage)
clean: 
	rm -rf *.out *.o *.gcno *.info s21_matrix_oop.a *.gcda *.a *.o *.txt
	rm -rf report 
	rm -rf test
leaks: 
	leaks -atExit --  ./test
clang-format:
	clang-format -i *.cpp *.h
valgrind:
	valgrind --tool=memcheck --leak-check=yes ./test 2> output.txt

# CC = g++ 
# SRC = s21_matrix_oop.cpp
# OBJ = s21_matrix_oop.o
# OBJ_TEST = s21_matrix_test.o
# LIB = s21_matrix_oop.a
# FLAGS = -Wall -Werror -Wextra -std=c++17 -I$(brew --prefix googletest)/include
# TESTFLAGS = -pthread -lgtest -lgmock -L$(brew --prefix googletest)/lib
# CHECK_FLAG = -lstdc++
# LCOVFLAGS = --ignore-errors inconsistent --ignore-errors mismatch
# OPEN = open 

# OS := $(shell uname -s) 
# ifeq ($(OS), Linux)
#     OPEN = xdg-open
#     TESTFLAGS = -pthread -lgtest -lgmock -lpthread -L$(brew --prefix googletest)/lib
#     FLAGS = -Wall -Werror -Wextra -std=c++17 -fPIC -fprofile-arcs -ftest-coverage -lstdc++fs -I$(brew --prefix googletest)/include
# endif

# all: clean gcov_report 

# s21_matrix_oop.a:
# 	$(CC) -c $(FLAGS) $(SRC) -o $(OBJ)
# 	ar -rcs $(LIB) $(OBJ)
# 	ranlib $(LIB)
# 	rm -rf *.o

# test: clean s21_matrix_oop.a
# 	$(CC) $(FLAGS) -c test.cpp -o test.o
# 	$(CC) $(FLAGS) $(CHECK_FLAG) test.o $(LIB) $(TESTFLAGS) -o test
# 	rm -rf *.o 

# gcov_report: add_coverage_flag test
# 	./test
# 	lcov -t "test" -o test.info -c -d ./  --no-external $(LCOVFLAGS)
# 	lcov --remove test.info 'test.cpp' -o test.info
# 	genhtml -o report test.info $(LCOVFLAGS)
# 	$(OPEN) report/index.html

# add_coverage_flag:
# 	$(eval FLAGS += --coverage)

# clean: 
# 	rm -rf *.out *.o *.gcno *.info s21_matrix_oop.a *.gcda *.a *.o *.txt
# 	rm -rf report 
# 	rm -rf test

# leaks: 
# 	leaks -atExit --  ./test

# clang-format:
# 	clang-format -i *.cpp *.h

# valgrind:
# 	valgrind --tool=memcheck --leak-check=yes ./test 2> output.txt

# CC = g++ 
# SRC = s21_matrix_oop.cpp
# OBJ = s21_matrix_oop.o
# OBJ_TEST = s21_matrix_test.o
# LIB = s21_matrix_oop.a
# FLAGS = -Wall -Werror -Wextra -std=c++17 -I/opt/homebrew/Cellar/googletest/1.16.0/include
# TESTFLAGS = -pthread -lgtest -lgmock -L/opt/homebrew/Cellar/googletest/1.16.0/lib
# CHECK_FLAG = -lstdc++
# LCOVFLAGS = --ignore-errors inconsistent --ignore-errors mismatch
# OPEN = open 

# OS := $(shell uname -s) 
# ifeq ($(OS), Linux)
#     OPEN = xdg-open
#     TESTFLAGS = -pthread -lgtest -lgmock -lpthread -L/opt/homebrew/Cellar/googletest/1.16.0/lib
#     FLAGS = -Wall -Werror -Wextra -std=c++17 -fPIC -fprofile-arcs -ftest-coverage -lstdc++fs -I/opt/homebrew/Cellar/googletest/1.16.0/include
# endif

# all: clean gcov_report 

# s21_matrix_oop.a:
# 	$(CC) -c $(FLAGS) $(SRC) -o $(OBJ)
# 	ar -rcs $(LIB) $(OBJ)
# 	ranlib $(LIB)
# 	rm -rf *.o

# test: clean s21_matrix_oop.a
# 	$(CC) $(FLAGS) -c test.cpp -o test.o
# 	$(CC) $(FLAGS) $(CHECK_FLAG) test.o $(LIB) $(TESTFLAGS) -o test
# 	rm -rf *.o 

# gcov_report: add_coverage_flag test
# 	./test
# 	lcov -t "test" -o test.info -c -d ./  --no-external $(LCOVFLAGS)
# 	lcov --remove test.info 'test.cpp' -o test.info
# 	genhtml -o report test.info $(LCOVFLAGS)
# 	$(OPEN) report/index.html

# add_coverage_flag:
# 	$(eval FLAGS += --coverage)

# clean: 
# 	rm -rf *.out *.o *.gcno *.info s21_matrix_oop.a *.gcda *.a *.o *.txt
# 	rm -rf report 
# 	rm -rf test

# leaks: 
# 	leaks -atExit --  ./test

# clang-format:
# 	clang-format -i *.cpp *.h

# valgrind:
# 	valgrind --tool=memcheck --leak-check=yes ./test 2> output.txt

