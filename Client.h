#include <string>
#include <gmp.h>
#include <iostream>
#include <unistd.h>
#include <cmath>
#include "DamgardJurik.h"

class Client {

private:

	TreeType tree;
	bool scalable;
	int max_s, bit_length;
	DamgardJurik* dj;

public:

	Client(int m_s, int b_length, TreeType t);
	Client(int m_s, int b_length, TreeType t, bool scalable);

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
	double encrypt_s_bits(mpz_t result[], int result_length,
			unsigned char s_bits[], int s_bit_length);

	/**
	 * Scalable method functions
	 */

	double get_scalable_s_bits(mpz_t encr_s_bits[], int encr_s_bit_length, unsigned char s_bits[], int s_bit_length,
			mpz_t encr_subtree[], int subtree_length, int selected_subtree);

	double decr_file(mpz_t dec, mpz_t file);

};

