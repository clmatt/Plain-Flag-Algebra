all: main.cpp
	g++ -g -Wall -o a.out main.cpp ./nauty27r2/nauty.a -I ./nauty27r2 -lm
	./a.out

