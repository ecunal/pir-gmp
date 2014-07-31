all:
	g++ -g Timing.cpp Server.cpp Client.cpp DamgardJurik.cpp -lrt -lgmp -lgmpxx -fopenmp -O3 -o main
