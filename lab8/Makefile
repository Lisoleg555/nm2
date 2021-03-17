all: lib.o main.o scatter_matrix.o
	g++ -g main.o scatter_matrix.o lib.o -o lab8
	rm *.o
main.o:
	g++ -c -g main.cpp
lib.o:
	g++ -c -g lib.cpp
scatter_matrix.o:
	g++ ../common/scatter_matrix.cpp -c -g
clear:
	rm *.o
