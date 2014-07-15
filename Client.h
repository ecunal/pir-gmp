#include <string>
#include <gmp.h>
#include <iostream>
#include <unistd.h>
#include "DamgardJurik.h"

class Client {

private:

	TreeType tree;
	int max_s, bit_length;
	DamgardJurik* dj;

public:

	Client(int m_s, int b_length, TreeType t);

	void get_pub_keys(mpz_t n, mpz_t g);

	/**
	 * Encrypts given selection bits.
	 *
	 * s_bits[]:		selection bit array, elements can be either 0 or 1
	 * s_bit_length:	size of s_bits array
	 * result[]:		encrypted selection bits array
	 * result_length:	size of result array.
	 * 					should be supplied according to tree type as follows:
	 * 				  	binary tree
	 * 				  		same as s_bit_length
	 * 				  	quadtree
	 * 				  		s_bit_length * 3/2
	 * 				  	octree
	 * 				  		s_bit_length * 7/3
	 *
	 */
	double encrypt_s_bits(mpz_t result[], int result_length, unsigned char s_bits[], int s_bit_length);


	double decr_file(mpz_t dec, mpz_t file);

};

