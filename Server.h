#include <string>
#include <gmp.h>
#include <iostream>
#include <unistd.h>
#include <cmath>
#include "DamgardJurik.h"

class Server {

private:
	DamgardJurik* dj;
	int max_s, bit_length /*dj parameters*/,
		f_size /*number of files*/;
	TreeType tree;
	gmp_randstate_t state; // for generating random files

	mpz_t *files; // file array

	void generate_files(bool debug);
	void init_random();

public:
	Server(int b_length, int s, int file_size, TreeType t, mpz_t n, mpz_t g);

	// if parallel = 0 then it is iterative.
	double get_file(mpz_t result, mpz_t s_bits[], int s_length, int parallel);
};
