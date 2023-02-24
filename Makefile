default:
	g++ -g main.cpp -larmadillo -o main

test: *.cpp *.h
	g++ -g test.cpp -larmadillo -lgtest -lgtest_main -pthread -o test
