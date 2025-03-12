CC=g++
FLAGS=-Wall -Werror -Wextra -std=c++17
FUNC=s21_matrix_oop.cpp
FUNC_O=s21_matrix_oop.o
TEST=test.cpp

all: gcov_report

test: s21_matrix_oop.a
	$(CC) $(FLAGS) -fprofile-arcs -ftest-coverage $(TEST) s21_matrix_oop.a -lgtest -o test

s21_matrix_oop.a:
	$(CC) -fprofile-arcs -ftest-coverage $(FLAGS) -c $(FUNC)
	ar rc s21_matrix_oop.a $(FUNC_O)
	ranlib s21_matrix_oop.a

gcov_report: test
	./test
	gcov --no-output $(TEST)
	lcov --no-external --capture --directory . --output-file coverage.info
	genhtml coverage.info -o report
	rm -rf *.o *.info *.gcda *.gcno *.gcov *.gch *.out *.a *.txt test

clean:
	rm -rf *.o report *.info *.gcda *.gcno *.gcov *.gch *.out *.a *.txt test
