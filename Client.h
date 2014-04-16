#include <string>
#include <gmp.h>
#include <iostream>
#include <unistd.h>
#include "DamgardJurik.h"

enum TreeType {
	BINARY, QUAD, OCTO
};

class Client {

private:

	TreeType tree;
	int max_s, bit_length;

public:

	DamgardJurik* dj;

	Client(int m_s, int b_length, TreeType t);

	void get_pub_keys(mpz_t n, mpz_t g);
	void encrypt_s_bits(mpz_t result[], int result_length, unsigned char s_bits[], int s_bit_length);

};

