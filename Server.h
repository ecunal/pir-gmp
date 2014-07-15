#include <string>
#include <gmp.h>
#include <iostream>
#include <unistd.h>
#include <cmath>
#include "DamgardJurik.h"

#define CORE_SIZE 4

class Server {

private:
	DamgardJurik* dj;
	int bit_length /*dj parameters*/,
		f_size /*number of files*/,
		subtree_size, scale /* for scalable method*/;
	TreeType tree;
	gmp_randstate_t state; // for generating random files

	mpz_t *files; // file array

	void generate_files(bool debug);
	void init_random();

public:



	/**
	 * s: maximum s to be used in encryptions
	 * b_length: bit length for DamgardJurik object
	 * file_size: number of files stored in the server
	 * t: BINARY, QUAD or OCTO
	 * n, g: public key for DamgardJurik, comes from client
	 */
	Server(int b_length, int file_size, TreeType t, mpz_t n, mpz_t g);
	Server(int b_length, int file_size, TreeType t, mpz_t n, mpz_t g, int sub);

	/* if parallel == 0 && extra_prl == 0 	-> 	iterative method
	 * if parallel == 1 && extra_prl == 0	-> 	first parallelization with only protocols in different threads
	 * if parallel == 1 && extra_prl == 1	-> 	second parallelization where
	 * 											encryption and exponentiation in each protocol are also parallelized.
	 *
	 *
	 * The return value is the computation time in seconds
	 * Ciphertext is returned in the variable result.
	 */

	double get_file(mpz_t result, mpz_t s_bits[], int parallel, int extra_prl);

	/**
	 * New method that distributes the tree into cores.
	 *
	 * The return value is the computation time in seconds
	 * Ciphertext is returned in the variable result.
	 */

	double get_file_new_p(mpz_t result, mpz_t s_bits[]);

	/**
	 * Scalable method
	 */
	double get_file_scalable(mpz_t result, mpz_t s_bits[], mpz_t subt_bits[], int subt_length);
};
